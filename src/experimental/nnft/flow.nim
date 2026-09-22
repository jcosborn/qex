## Masked U(1) stages, with explicit neural precision and double gauge storage.
import qex, comms/halo
import ../../nn
import gauge/stoutsmear
import std/math
export halo, RealField, realField

const
  nnftFeatures* = 6
    ## Network inputs: masked sin/cos of the plaquette and both rectangles.
  nnftCoefs* = 12
    ## Network outputs: plaquette and rectangle coefficients of both directions.

type
  NnftNet*[T: SomeFloat] = object
    ## Convolutions built by convParams with GELU between them, then the
    ## output scales; the widths chain nnftFeatures -> ... -> nnftCoefs.
    layers*: seq[ConvParams[T]]
    scale*: seq[T]
  NnftParams*[T: SomeFloat] = seq[NnftNet[T]]
    ## One network per stage; stage s updates the class s mod 8.
  GaugeElem = DLatticeColorMatrixV.T
  RealElem = DLatticeRealV.T
  GaugeHalo = Halo[Layout[VLEN],DLatticeColorMatrixV,GaugeElem]
  LinkTerm = object
    mu, sign, tap: int
    off: array[2,int32]
  NnftGeom = ref object
    layout: HaloLayout[Layout[VLEN]]
    map: HaloMap[Layout[VLEN]]
    halo: array[2,GaugeHalo]
    back: array[2,GaugeHalo]
    offsets: seq[seq[int32]]
    index: seq[int32]
    loops: array[3,seq[LinkTerm]]
    staples: array[2,array[6,seq[LinkTerm]]]
  NnftStage*[T: SomeFloat] = ref object
    ## Tape of one stage. `m` = V_d ds^dagger is meaningful at the active sites
    ## of the stage direction only.
    input, loops, openStaples, ds: seq[DLatticeColorMatrixV]
    output*: seq[DLatticeColorMatrixV]
    features: seq[RealField[T]]
    convSum, conv, gelu: seq[seq[RealField[T]]] # per convolution; gelu skips the last
    scaled, atan, divPi: seq[RealField[T]]
    coef*: seq[RealField[T]]
    active*: SiteMask
    featureMasks*: array[3,SiteMask]
    stage: int
    geom: NnftGeom
    logdet*: float
    expa: DLatticeColorMatrixV
    m*: DLatticeColorMatrixV
    net: seq[ConvWorkspace[T]]
    dsGrad, loopGrad, stapleGrad: seq[DLatticeColorMatrixV]
    coefGrad, tmpGrad1, tmpGrad2: seq[RealField[T]]
    grad: seq[seq[RealField[T]]] # per convolution, its input cotangent

const
  stapleAnchors* = [
    [[0,0],[0,-1],[-1,0],[-1,-1],[0,0],[0,-1]],
    [[0,0],[-1,0],[0,-1],[-1,-1],[0,0],[-1,0]]]
    ## Loop origin relative to the active link of direction d for open staple k.
  coefSlots* = [[0,1,4,5,6,7],[2,3,8,9,10,11]]
    ## Network output channel weighting open staple k of direction d.
  nnftJacFloor* = 1e-8
    ## Floor of 1+Re M inside the stage log determinant. JAX max convention:
    ## slope 0 below, 1/2 at equality, 1 above.

proc stageClass*(s: int): tuple[dir, p0, p1, parity: int] =
  ## Stage s changes direction dir at sites of row parity p0 and column parity
  ## p1, inside the even/odd subset parity; the eight classes repeat.
  if s < 0:
    raise newException(ValueError, "NNFT stage index must be nonnegative")
  let c = s mod 8
  result.dir = c div 4
  result.p0 = (c mod 4) div 2
  result.p1 = c mod 2
  result.parity = (result.p0+result.p1) mod 2

proc gauges(lo: Layout[VLEN]; n: int): seq[DLatticeColorMatrixV] =
  result = newSeq[DLatticeColorMatrixV](n)
  for i in 0..<n: result[i] = lo.ColorMatrixD()

proc reals[T: SomeFloat](lo: Layout[VLEN]; n: int; precision: typedesc[T]): seq[RealField[T]] =
  result = newSeq[RealField[T]](n)
  for i in 0..<n: result[i] = realField(lo,T)

proc term(mu, r, c, sign: int): LinkTerm =
  LinkTerm(mu: mu, off: [int32(r),int32(c)], sign: sign)

proc newGeom(g: seq[DLatticeColorMatrixV]): NnftGeom =
  result.new
  result.loops[0] = @[term(0,0,0,1), term(1,1,0,1), term(0,0,1,-1), term(1,0,0,-1)]
  result.loops[1] = @[term(0,0,0,1), term(0,1,0,1), term(1,2,0,1), term(0,1,1,-1), term(0,0,1,-1), term(1,0,0,-1)]
  result.loops[2] = @[term(0,0,0,1), term(1,1,0,1), term(1,1,1,1), term(0,0,2,-1), term(1,0,1,-1), term(1,0,0,-1)]
  const signs = [[1,-1,1,-1,1,-1],[-1,1,-1,1,-1,1]]
  for d in 0..1:
    for k in 0..<6:
      let lp = if k < 2: 0 else: d+1
      var removed = 0
      for t in result.loops[lp]:
        var q = t
        q.off[0] += int32(stapleAnchors[d][k][0])
        q.off[1] += int32(stapleAnchors[d][k][1])
        if q.mu == d and q.off == [0'i32,0'i32]:
          doAssert q.sign == signs[d][k]
          inc removed
        else:
          # L = W^a Q, hence D = Q^(-a) and W D† = L^a.
          q.sign *= -signs[d][k]
          result.staples[d][k].add q
      doAssert removed == 1
  proc addPath(geo: NnftGeom; path: var seq[LinkTerm]) =
    for t in path.mitems:
      let off = @[t.off[0],t.off[1]]
      var i = geo.offsets.find(off)
      if i < 0:
        i = geo.offsets.len
        geo.offsets.add off
      t.tap = i
  for p in result.loops.mitems: result.addPath(p)
  for dir in result.staples.mitems:
    for p in dir.mitems: result.addPath(p)
  let lo = g[0].l
  result.layout = haloLayout(lo,result.offsets)
  result.map = haloMap(result.layout,lo.comm,result.offsets)
  for d in 0..1:
    result.halo[d] = result.layout.makeHalo(g[d])
    result.back[d] = result.layout.makeHalo(g[d])
  result.index = result.layout.haloIndex(result.offsets)

proc requireNnftGauge*(g: seq[DLatticeColorMatrixV]) =
  if g.len != 2 or g[0].l.nDim != 2 or g[0][0].ncols != 1:
    raise newException(ValueError, "NNFT learned flow requires two-dimensional U(1) fields")
  let lo = g[0].l
  for d in lo.physGeom:
    if d < 4 or d mod 2 != 0:
      raise newException(ValueError, "NNFT learned flow requires even extents of at least four")
  if g[1].l != lo:
    raise newException(ValueError, "NNFT gauge layouts differ")

proc requireNet*[T: SomeFloat](p: NnftNet[T]) =
  if p.layers.len == 0:
    raise newException(ValueError, "NNFT network requires convolutions")
  var cin = nnftFeatures
  for q in p.layers:
    if q.cin != cin:
      raise newException(ValueError, "NNFT convolution input channels must chain from " & $nnftFeatures)
    cin = q.cout
  if cin != nnftCoefs or p.scale.len != nnftCoefs:
    raise newException(ValueError, "NNFT network must end in " & $nnftCoefs & " channels and scales")

proc nnftMasks*(lo: Layout[VLEN]; s: int): tuple[active: SiteMask, featureMasks: array[3,SiteMask]] =
  ## Masks for the learned map's validated two-dimensional layout.
  let cls = stageClass(s)
  let active = realField(lo,float32)
  let features = [realField(lo,float32),realField(lo,float32),realField(lo,float32)]
  threads:
    for m in features: m := 0
    for x in lo.sites:
      let row = lo.coords[0][x] mod 2
      let col = lo.coords[1][x] mod 2
      let keep = row == cls.p0 and col == cls.p1
      active{x} := (if keep: 1'f32 else: 0'f32)
      features[0]{x} := (if (if cls.dir == 0: row != cls.p0 else: col != cls.p1): 1'f32 else: 0'f32)
      features[2-cls.dir]{x} := (if keep: 0'f32 else: 1'f32)
  (active,features)

proc newNnftStage*[T: SomeFloat](g: seq[DLatticeColorMatrixV]; p: NnftNet[T]; s: int): NnftStage[T] =
  g.requireNnftGauge
  p.requireNet
  let lo = g[0].l
  let masks = nnftMasks(lo,s)
  let t = NnftStage[T](stage:s,geom:newGeom(g),active:masks.active,featureMasks:masks.featureMasks)
  t.input = gauges(lo,2)
  t.output = gauges(lo,2)
  t.loops = gauges(lo,3)
  t.openStaples = gauges(lo,6)
  t.ds = gauges(lo,2)
  t.expa = lo.ColorMatrixD()
  t.m = lo.ColorMatrixD()
  t.features = reals(lo,nnftFeatures,T)
  for l, q in p.layers:
    t.convSum.add reals(lo,q.cout,T)
    t.conv.add reals(lo,q.cout,T)
    if l < p.layers.high: t.gelu.add reals(lo,q.cout,T)
    t.net.add convWorkspace(t.features[0],q)
    t.grad.add reals(lo,q.cin,T)
  t.scaled = reals(lo,nnftCoefs,T)
  t.atan = reals(lo,nnftCoefs,T)
  t.divPi = reals(lo,nnftCoefs,T)
  t.coef = reals(lo,nnftCoefs,T)
  t.dsGrad = gauges(lo,2)
  t.loopGrad = gauges(lo,3)
  t.stapleGrad = gauges(lo,6)
  t.coefGrad = reals(lo,nnftCoefs,T)
  t.tmpGrad1 = reals(lo,nnftCoefs,T)
  t.tmpGrad2 = reals(lo,nnftCoefs,T)
  t

proc evalNet*[T: SomeFloat](t: NnftStage[T]; p: NnftNet[T]; features: seq[RealField[T]]) =
  var x = features
  for l, q in p.layers:
    conv(t.convSum[l], x, q, t.net[l], addBias=false)
    bias(t.conv[l], t.convSum[l], q.bias)
    x = t.conv[l]
    if l < p.layers.high:
      gelu(t.gelu[l], x)
      x = t.gelu[l]
  scale(t.scaled, x, p.scale)
  arctan(t.atan, t.scaled)
  divide(t.divPi, t.atan, T(PI))
  divide(t.coef, t.divPi, T(3))

proc learnedLogJac*[T: SomeFloat](j: T): T =
  ln(max(T(1)+j,T(nnftJacFloor)))

proc evalStage*[T: SomeFloat](t: NnftStage[T]; p: NnftNet[T]; g: seq[DLatticeColorMatrixV]): float =
  let geo = t.geom
  let lo = geo.layout.lo
  if g.len != 2 or g[0].l != lo or g[1].l != lo:
    raise newException(ValueError, "NNFT input gauge does not match workspace")
  let cls = stageClass(t.stage)
  let d = cls.dir
  let sub = lo.getSubset(if cls.parity == 0: "even" else: "odd")
  threads:
    for mu in 0..1:
      t.input[mu] := g[mu]
      t.ds[mu] := 0
  for mu in 0..1:
    geo.halo[mu].field = t.input[mu]
    geo.halo[mu].update(geo.map, lo.comm)
  let nt = geo.offsets.len
  # The masks are exactly zero or one, so multiplying by them selects exactly.
  threads:
    for k in 0..<3:
      for x in t.loops[k]:
        var q: GaugeElem
        q := 1
        for a in geo.loops[k]:
          let v = geo.halo[a.mu][geo.index[x*nt+a.tap]]
          if a.sign > 0: q := q*v
          else: q := q*v.adj
        t.loops[k][x] := q
    for k in 0..<6:
      for x in t.openStaples[k]:
        var q: GaugeElem
        q := 1
        for a in geo.staples[d][k]:
          let v = geo.halo[a.mu][geo.index[x*nt+a.tap]]
          if a.sign > 0: q := q*v
          else: q := q*v.adj
        var a: RealElem
        a := t.active[x]
        q[0,0] := a*q[0,0]
        t.openStaples[k][x] := q
    for k in 0..<3:
      let si = if k == 0: 0 else: k+1
      let ci = if k == 0: 1 else: k+3
      for x in t.features[si]:
        var m: evalType(t.features[si][x])
        m := t.featureMasks[k][x]
        t.features[si][x] := m*t.loops[k][x][0,0].im
        t.features[ci][x] := m*t.loops[k][x][0,0].re+(T(1)-m)
  evalNet(t,p,t.features)
  threads:
    for x in sub:
      var v: evalType(t.ds[d][x][0,0])
      v := 0
      for k in 0..<6:
        var c: RealElem
        c := t.coef[coefSlots[d][k]][x]
        v := v-c*t.openStaples[k][x][0,0]
      t.ds[d][x][0,0] := v
  stoutStepKernel(t.output, t.input, t.ds, t.expa, t.m, 1.0, cls.parity, d, t.active)
  var ld = 0.0
  threads:
    # SIMD locals must stay out of the closure environment with Nim 2.0 refc.
    var floor: RealElem
    floor := nnftJacFloor
    var v = 0.0
    for x in sub:
      var a, j: RealElem
      a := t.active[x]
      j := t.m[x][0,0].re
      v += simdSum(a*ln(max(1.0+j,floor)))
    lo.threadRankSum(v)
    threadSingle: ld = v
  t.logdet = ld
  ld

proc netVjp[T: SomeFloat](t: NnftStage[T]; p: NnftNet[T]) =
  divideVjp(t.tmpGrad1,t.coefGrad,T(3))
  divideVjp(t.tmpGrad2,t.tmpGrad1,T(PI))
  arctanVjp(t.tmpGrad1,t.scaled,t.tmpGrad2)
  scaleVjp(t.tmpGrad2,t.tmpGrad1,p.scale)
  var up = t.tmpGrad2
  for l in countdown(p.layers.high,0):
    convVjp(t.grad[l],up,p.layers[l],t.net[l])
    if l > 0: geluVjp(t.grad[l],t.conv[l-1],t.grad[l])
    up = t.grad[l]

proc stageVjp*[T: SomeFloat](t: NnftStage[T]; p: NnftNet[T]; dy: seq[DLatticeColorMatrixV]; dl: float; dx: seq[DLatticeColorMatrixV]) =
  ## Re sum conj(dx) dW = Re sum conj(dy) dWnew + dl d(logdet).
  ## The tape comes from evalStage at the matching input and parameters.
  let geo = t.geom
  let lo = geo.layout.lo
  if dy.len != 2 or dx.len != 2:
    raise newException(ValueError, "NNFT pullback requires two input/output cotangent fields")
  for mu in 0..1:
    if dy[mu].l != lo or dx[mu].l != lo:
      raise newException(ValueError, "NNFT pullback layouts differ from the tape")
    for nu in 0..1:
      if dx[mu] == dy[nu]:
        raise newException(ValueError, "NNFT pullback destination must not alias its seed")
  let cls = stageClass(t.stage)
  let d = cls.dir
  let sub = lo.getSubset(if cls.parity == 0: "even" else: "odd")
  threads:
    for mu in 0..1:
      dx[mu] := dy[mu]
      t.dsGrad[mu] := 0
    for k in 0..<nnftCoefs: t.coefGrad[k] := 0
    for k in 0..<3: t.loopGrad[k] := 0
    for k in 0..<6: t.stapleGrad[k] := 0
  # Lanes of one SIMD site are handled separately: each takes the clipped log
  # determinant's slope on its own side of the floor.
  threads:
    for x in sub.sites:
      if selected(t.active,x):
        var w, v, e, m, u, rw, rd: evalType(t.input[d]{x})
        w := t.input[d]{x}
        v := t.ds[d]{x}
        e := t.expa{x}
        m := t.m{x}
        u := dy[d]{x}
        var j: float64
        j := m[0,0].re
        let diag = 1.0+j
        let seed = if diag < nnftJacFloor: 0.0 elif diag == nnftJacFloor: 0.5*dl else: dl
        stoutPullbackSite(rw,rd,w,v,e,m,u,1.0,seed)
        dx[d]{x} := rw
        t.dsGrad[d]{x} := rd
        for k in 0..<6:
          var path: evalType(w)
          path := t.openStaples[k]{x}
          t.coefGrad[coefSlots[d][k]]{x} := -redot(rd,path)
          var c: float64
          c := t.coef[coefSlots[d][k]]{x}
          t.stapleGrad[k]{x} := -c*rd
  netVjp(t,p)
  threads:
    for k in 0..<3:
      let si = if k == 0: 0 else: k+1
      let ci = if k == 0: 1 else: k+3
      for x in t.loopGrad[k]:
        var m: RealElem
        m := t.featureMasks[k][x]
        t.loopGrad[k][x][0,0].im := m*t.grad[0][si][x]
        t.loopGrad[k][x][0,0].re := m*t.grad[0][ci][x]
  for mu in 0..1:
    geo.halo[mu].field = t.input[mu]
    geo.halo[mu].update(geo.map,lo.comm)
    geo.back[mu].field = dx[mu]
    geo.back[mu].halo := 0
  let nt = geo.offsets.len
  template pathVjp(path: seq[LinkTerm]; bar: DLatticeColorMatrixV; mu: int) =
    for x in 0..<geo.layout.nOut:
      let b = bar[x]
      for j, a in path:
        if a.mu == mu:
          var q, v: GaugeElem
          q := 1
          for k, z in path:
            if k != j:
              let h = geo.halo[z.mu][geo.index[x*nt+z.tap]]
              if z.sign > 0: q := q*h
              else: q := q*h.adj
          if a.sign > 0: v := b*q.adj
          else: v := b.adj*q
          geo.back[mu][geo.index[x*nt+a.tap]] += v
  # Each thread owns complete input-direction fields, including their halos.
  threads:
    tfor mu, 0..<2:
      for k in 0..<3: pathVjp(geo.loops[k],t.loopGrad[k],mu)
      for k in 0..<6: pathVjp(geo.staples[d][k],t.stapleGrad[k],mu)
  for mu in 0..1: geo.back[mu].updateRev(geo.map,lo.comm)
