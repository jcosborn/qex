import base, layout, field, maths, base/hyper, comms/gather, tables
import physics/color
from physics/qcdTypes import perm
getOptimPragmas()

type
  HaloLayout*[L] = ref object
    lo*: L  # layout
    outerExt*: seq[int32]  # extended outer lattice size
    offset*: seq[int32]  # offset of outer lattice in extended outer
    lex*: seq[int32]  # outerExt lex index of extended sites
    index*: seq[int32]  # index of extended site for given lex index
    neighborFwd*: seq[seq[int32]]  # fwd neighbor for extended outer lattice
    neighborBck*: seq[seq[int32]]  # bck neighbor for extended outer lattice
    nOut*: int  # sites in outer lattice
    nExt*: int  # sites in extended outer lattice
    nOutPar*: array[2,int]  # sites in outer lattice by parity
    nExtPar*: array[2,int]  # sites in extended outer lattice by parity
  Halo*[L,F,T] = ref object
    layout*: HaloLayout[L]
    field*: F
    halo*: alignedMem[T]
    nOut*: int  # sites in outer lattice
    nExt*: int  # sites in extended outer lattice
    nOutPar*: array[2,int]  # sites in outer lattice by parity
    nExtPar*: array[2,int]  # sites in extended outer lattice by parity
  HaloMap*[L] = ref object
    layout*: HaloLayout[L]
    offsets*: seq[seq[int32]]  # offsets of outer sites needed
    gather*: GatherMap  # lanes of the needed halo sites that sites does not fill
    gatherRev*: GatherMap  # gather map to send halo sites back
    sites*: GatherMap  # whole halo sites, lane k from lane k xor perm of one source site
    perm: seq[int32]  # that xor of each halo site in sites, -1 elsewhere
    nsites*, nlanes*: int  # entries of sites and gather on all ranks
  GatherHalo*[F,T;Rev:static bool] = object
    src: F
    dest: ptr UncheckedArray[T]
    elemsize: int
    vlen: int
    perm: ptr UncheckedArray[int32]  # whole sites, see copy

proc norm2*(h: Halo): auto =
  var s = h.halo[0].norm2
  for i in 1..<h.halo.len:
    s += h.halo[i].norm2
  var t = s.simdReduce
  h.field.l.threadRankSum(t)
  t += h.field.norm2
  t

proc sum*(h: Halo): auto =
  var s = h.halo[0]
  for i in 1..<h.halo.len:
    s += h.halo[i]
  var t = s.simdSum
  h.field.l.threadRankSum(t)
  t += h.field.sum
  t

proc inside[T,U,V:openArray[auto]](x: T, lo: U, hi: V): bool {.inline.} =
  result = true
  for i in 0..<x.len:
    if x[i] < lo[i] or x[i] >= hi[i]:
      result = false
      break

proc init[T](s: var seq[T], x: openArray[SomeNumber]) =
  s.newSeq(x.len)
  for i in 0..<x.len: s[i] = T x[i]

# TODO: order halo by layers
proc makeHaloLayout*[L:Layout](lo: L, fwdOffset,bckOffset: openarray[SomeInteger]): HaloLayout[L] =
  tic("makeHaloLayout")
  let outerHigh = lo.outerGeom + bckOffset
  let outerExt = outerHigh + fwdOffset
  let nExt = outerExt.product
  let nOut = lo.nSitesOuter
  var lex = newSeq[int32](nExt)
  var index = newSeq[int32](nExt)
  var idxp = [newSeq[int32](0), newSeq[int32](0)]  # store indices by parity
  let paro = (lo.coordmin.sum + bckOffset.sum) mod 2  # global parity of extended origin
  let nd = lo.nDim
  var x = newSeq[int32](nd)
  for i in 0..<nExt:  # loop over lex indices
    x.lexCoord(i, outerExt)
    if x.inside(bckOffset, outerHigh):
      var y = x - bckOffset
      let loc = lo.localIndex(y) div L.V
      lex[loc] = int32 i
      index[i] = int32 loc
    else:
      let par = (paro + x.sum) mod 2
      idxp[par].add int32 i
  var k = int32 nOut
  for p in 0..1:
    for i in idxp[p]:
      lex[k] = int32 i
      index[i] = int32 k
      inc k
  toc("lex/index")
  result.new
  result.lo = lo
  result.outerExt.init outerExt
  result.offset.init bckOffset
  result.lex = lex
  result.index = index
  result.nOut = nOut
  result.nExt = nExt
  result.nOutPar[0] = lo.nEvenOuter
  result.nOutPar[1] = lo.nOddOuter
  result.nExtPar[0] = result.nOutPar[0] + idxp[0].len
  result.nExtPar[1] = result.nOutPar[1] + idxp[1].len
  result.neighborFwd = newSeq[seq[int32]](nd)
  result.neighborBck = newSeq[seq[int32]](nd)
  for mu in 0..<nd:
    result.neighborFwd[mu] = newSeq[int32](nExt)
    result.neighborBck[mu] = newSeq[int32](nExt)
  for i in 0..<nExt:  # loop over lex indices
    let ii = index[i]
    x.lexCoord(i, outerExt)
    for mu in 0..<nd:
      if x[mu] == 0:
        result.neighborBck[mu][ii] = -1
      x[mu] += 1
      if x[mu] < outerExt[mu]:
        let j = x.lexIndex(outerExt)
        let ij = index[j]
        result.neighborFwd[mu][ii] = ij
        result.neighborBck[mu][ij] = ii
      else:
        #echo x
        result.neighborFwd[mu][ii] = -1
      x[mu] -= 1
  toc("neighbors")

proc makeHaloMap*[L](hl: HaloLayout[L], c: Comm, offsets: seq[seq[int32]]): HaloMap[L] =
  tic("makeHaloMap")
  result.new
  result.layout = hl
  result.offsets = offsets
  let outerHigh = hl.offset + hl.lo.outerGeom
  let nd = hl.lo.nDim
  var x = newSeq[int32](nd)
  var rl = newSeq[RecvList](0)
  var sl = newSeq[SendList](0)
  let vlen = L.V
  var vecOffset = newSeq[seq[int32]](vlen)
  for k in 0 ..< vlen:
    vecOffset[k] = newSeq[int32](nd)
    let li0 = hl.offset.lexIndex(hl.outerExt)
    let i0 = hl.index[li0]
    let ik = vlen*i0 + k
    hl.lo.coord(x, ik)
    vecOffset[k] = x - hl.offset
  # A halo site whose lane k holds lane k xor p of one source site moves as a
  # whole site, permuted by p, and so does the local part of a site with
  # lanes on other ranks: the site gather writes every lane before the lane
  # gather overwrites the remote ones. Other lanes move one by one. Inner
  # extents are 1 or 2 (layoutX), so lanes from one site differ by an xor.
  let rank = c.commrank
  var ql = newSeq[RecvList](0)
  var ri = newSeq[tuple[rank,index:int]](vlen)
  result.perm = newSeq[int32](hl.nExt - hl.nOut)
  for i in hl.nOut ..< hl.nExt:  # loop over halo sites
    result.perm[i-hl.nOut] = -1
    let li = hl.lex[i]
    x.lexCoord(li, hl.outerExt)
    var keep = false
    for j in 0 ..< offsets.len:
      let y = x - offsets[j]
      if y.inside(hl.offset, outerHigh):
        keep = true
        break
    if keep:
      var k0 = 0  # the first local lane, else lane 0, names the source site
      for k in countdown(vlen-1, 0):
        var y = x + vecOffset[k]
        for l in 0..<nd:
          if y[l] < 0: y[l] += int32 hl.lo.physGeom[l]
          if y[l] >= hl.lo.physGeom[l]: y[l] -= int32 hl.lo.physGeom[l]
        ri[k] = hl.lo.rankIndex(y)
        let didx = vlen*(i-hl.nOut)+k
        sl.add SendList(sidx: int32 didx, drank: int32 ri[k].rank, didx: int32 ri[k].index)
        if ri[k].rank == rank: k0 = k
      let
        r = ri[k0].rank
        o = ri[k0].index div vlen
        p = (ri[k0].index mod vlen) xor k0
      template covered(k: int): bool =
        ri[k].rank == r and ri[k].index div vlen == o and ri[k].index mod vlen == (k xor p)
      var all, mine = true
      for k in 0 ..< vlen:
        if not covered(k):
          all = false
          if ri[k].rank == r: mine = false
      if all or (r == rank and mine):
        result.perm[i-hl.nOut] = int32 p
        ql.add RecvList(didx: int32(i-hl.nOut), srank: int32 r, sidx: int32 o)
      for k in 0 ..< vlen:
        if result.perm[i-hl.nOut] < 0 or not covered(k):
          rl.add RecvList(didx: int32(vlen*(i-hl.nOut)+k), srank: int32 ri[k].rank, sidx: int32 ri[k].index)
  toc("RecvList")
  var n = [float rl.len, float ql.len]
  c.allReduce(n[0].addr, 2)
  result.nlanes = int n[0]
  result.nsites = int n[1]
  result.gather = c.makeGatherMap(rl)
  result.sites = c.makeGatherMap(ql)
  result.gatherRev = c.makeGatherMap(sl)
  toc("makeGatherMap")
proc makeHaloMap*[L,N](hl: HaloLayout[L], c: Comm, offsets: seq[array[N,SomeInteger]]): HaloMap[L] =
  var o = newSeq[seq[int32]](offsets.len)
  for i in 0..<offsets.len:
    let t = offsets[i]
    o[i] = @t
  makeHaloMap(hl, c, o)

# Memoized constructors, argument order as makeHaloLayout/makeHaloMap: one
# table per L instantiation, keyed by the address of lo (hl) plus the widths
# (offsets).  The cached value holds lo (hl), so its address cannot be reused
# while cached.  Not thread safe: build outside `threads:`.
proc haloLayout*[L:Layout](lo: L, fwd,bck: seq[int32]): HaloLayout[L] =
  var cache {.global.}: Table[(pointer,seq[int32],seq[int32]), HaloLayout[L]]
  let key = (cast[pointer](lo), fwd, bck)
  result = cache.getOrDefault(key)
  if result.isNil:
    result = makeHaloLayout(lo, fwd, bck)
    cache[key] = result
proc haloLayout*[L:Layout](lo: L, fwdOffset,bckOffset: openarray[SomeInteger]): HaloLayout[L] =
  var f,b: seq[int32]
  f.init fwdOffset
  b.init bckOffset
  haloLayout(lo, f, b)

proc haloLayout*[L:Layout](lo: L, offsets: seq[seq[int32]]): HaloLayout[L] =
  ## Widths from the extreme offsets, through the memoized constructor above.
  ## Offsets must stay inside one period: the map folds a coordinate once.
  let nd = lo.nDim
  var fw = newSeq[int32](nd)
  var bw = newSeq[int32](nd)
  for off in offsets:
    if off.len != nd:
      raise newException(ValueError, "halo offset dimension differs from layout")
    for d in 0..<nd:
      if abs(off[d]) >= lo.physGeom[d]:
        raise newException(ValueError, "halo offset exceeds the periodic extent")
      fw[d] = max(fw[d], off[d])
      bw[d] = max(bw[d], -off[d])
  haloLayout(lo, fw, bw)

proc haloIndex*[L](hl: HaloLayout[L], offsets: seq[seq[int32]]): seq[int32] =
  ## Extended index of outer site i shifted by offsets[t], at i*offsets.len+t.
  ## Not cached; each caller owns its table. Offsets must fit the widths.
  let nd = hl.lo.nDim
  let nt = offsets.len
  for off in offsets:
    for d in 0..<nd:
      let fwd = int(hl.outerExt[d]) - int(hl.offset[d]) - hl.lo.outerGeom[d]
      if off[d] > fwd or -off[d] > hl.offset[d]:
        raise newException(ValueError, "halo offset exceeds the layout widths")
  result = newSeq[int32](hl.nOut*nt)
  var x = newSeq[int32](nd)
  for i in 0..<hl.nOut:
    x.lexCoord(hl.lex[i], hl.outerExt)
    for t, off in offsets:
      result[i*nt+t] = hl.index[(x+off).lexIndex(hl.outerExt)]

proc haloMap*[L](hl: HaloLayout[L], c: Comm, offsets: seq[seq[int32]]): HaloMap[L] =
  var cache {.global.}: Table[(pointer,seq[seq[int32]]), HaloMap[L]]
  let key = (cast[pointer](hl), offsets)
  result = cache.getOrDefault(key)
  if result.isNil:
    result = makeHaloMap(hl, c, offsets)
    cache[key] = result
proc haloMap*[L,N](hl: HaloLayout[L], c: Comm, offsets: seq[array[N,SomeInteger]]): HaloMap[L] =
  var o = newSeq[seq[int32]](offsets.len)
  for i in 0..<offsets.len:
    o[i].init offsets[i]
  haloMap(hl, c, o)

# T is the vectorized site type. The buffer is not zeroed: only cells a map
# fills are defined, and the reverse exchange reads every cell.
proc makeHalo*[L,F,T](hl: HaloLayout[L], f: F, t: typedesc[T]): Halo[L,F,T] =
  result.new
  result.layout = hl
  result.field = f
  result.nOut = hl.nOut
  result.nExt = hl.nExt
  result.nOutPar = hl.nOutPar
  result.nExtPar = hl.nExtPar
  let nhalo = result.nExt - result.nOut
  result.halo.newU(nhalo)

template makeHalo*[L,F](hl: HaloLayout[L], f: F): auto =
  #makeHalo(hl, f, eval(F.type[0]))
  makeHalo(hl, f, eval(F.type.index(int)))

proc setLanes[T](r: var T, x: T, p: int32) {.inline.} =
  ## r = x with lane k of r taken from lane k xor p of x.
  if p == 0: r = x
  else:
    var t = x
    var b = 1'i32
    while b <= p:
      if (p and b) != 0:
        var u {.noinit.}: T
        perm(u, int b, t)
        t = u
      b = b shl 1
    r = t

template copyLanes[F,T;Rev:static bool](gh: GatherHalo[F,T,Rev], d, s: pointer|SomeInteger) =
  when d is pointer:
    type E = eval(index(type T, type asSimd(0)))
    let p = cast[ptr E](d)
    when Rev:
      let o = s div gh.vlen
      let i = s mod gh.vlen
      p[] := gh.dest[o][asSimd(i)]
    else:
      p[] := gh.src{s}
  elif s is pointer:
    type E = eval(index(type T,type asSimd(0)))
    let p = cast[ptr E](s)
    when Rev:
      gh.src{d} += p[]
    else:
      let o = d div gh.vlen
      let i = d mod gh.vlen
      gh.dest[o][asSimd(i)] = p[]
  else:
    when Rev:
      let o = s div gh.vlen
      let i = s mod gh.vlen
      gh.src{d} += gh.dest[o][asSimd(i)]
    else:
      let o = d div gh.vlen
      let i = d mod gh.vlen
      gh.dest[o][asSimd(i)] = gh.src{s}
template copy*[F,T;Rev:static bool](gh: GatherHalo[F,T,Rev], d, s: pointer|SomeInteger) =
  # Keep all three gather transfers in one binding for indirect generic callers.
  # With perm (forward only), elements are whole sites and lane k of halo site d
  # is lane k xor perm[d] of its source; message buffers are not aligned for T.
  when not Rev:
    if gh.perm != nil:
      when d is pointer:
        copyMem(d, gh.src[s].addr, sizeof(T))
      elif s is pointer:
        var t {.noinit.}: T
        copyMem(t.addr, s, sizeof(T))
        setLanes(gh.dest[d], t, gh.perm[d])
      else:
        setLanes(gh.dest[d], gh.src[s], gh.perm[d])
    else:
      copyLanes(gh, d, s)
  else:
    copyLanes(gh, d, s)

proc update*[L,F,T](h: Halo[L,F,T], hm: HaloMap[L], c: Comm) =
  bind copy
  tic("Halo update")
  # The gathers borrow the halo and permutation storage through raw pointers.
  GC_ref(h)
  GC_ref(hm)
  defer:
    GC_unref(hm)
    GC_unref(h)
  if hm.nsites > 0:
    var gs: GatherHalo[F,T,false]
    gs.src = h.field
    gs.dest = cast[ptr UncheckedArray[T]](unsafeAddr h.halo[0])
    gs.perm = cast[ptr UncheckedArray[int32]](unsafeAddr hm.perm[0])
    gs.elemsize = sizeof(T)
    gs.vlen = 1
    c.gather(hm.sites, gs, keep = true)
    toc("sites")
  if hm.nlanes > 0:
    let elemSize = sizeof(T) div L.V
    var gh: GatherHalo[F,T,false]
    gh.src = h.field
    gh.dest = cast[ptr UncheckedArray[T]](unsafeAddr h.halo[0])
    gh.elemsize = elemSize
    gh.vlen = L.V
    c.gather(hm.gather, gh, keep = true)
    toc("gather")

proc updateRev*[L,F,T](h: Halo[L,F,T], hm: HaloMap[L], c: Comm) =
  bind copy
  tic("Halo updateRev")
  let elemSize = sizeof(T) div L.V
  var gh: GatherHalo[F,T,true]
  gh.src = h.field
  gh.dest = cast[ptr UncheckedArray[T]](unsafeAddr h.halo[0])
  gh.elemsize = elemSize
  gh.vlen = L.V
  GC_ref(h)  # as in update
  defer: GC_unref(h)
  c.gather(hm.gatherRev, gh)
  toc("gather")


proc indexPtr*[L,F,T](h: Halo[L,F,T], i: SomeInteger): ptr T {.alwaysInline.} =
  let k = i - h.nOut
  result = if k<0: addr h.field[i] else: addr h.halo[k]
template `[]`*(h: Halo, i: SomeInteger): auto = indexPtr(h,i)[]

proc `[]=`*(h: Halo, i: SomeInteger, x: auto) {.alwaysInline.} =
  let k = i - h.nOut
  if k < 0:
    h.field[i] = x
  else:
    h.halo[k] := x

proc step*[I: SomeInteger](hl: HaloLayout, d: int, i: I, n: SomeInteger): I =
  ## The extended index n sites from i along d; n < 0 walks backward.
  result = i
  if n > 0:
    for _ in 0..<int(n): result = I(hl.neighborFwd[d][result])
  else:
    for _ in 0..<int(-n): result = I(hl.neighborBck[d][result])

proc neighbor*(h: Halo, i: SomeInteger, mu: SomeInteger, fb: SomeInteger): int32 =
  if fb > 0:
    h.layout.neighborFwd[mu][i]
  else:
    h.layout.neighborBck[mu][i]

when isMainModule:
  import qex
  qexInit()
  tic("main")
  var defaultLat = @[4,4,4,4]
  defaultSetup()
  let nd = lo.nDim
  var seed = 987654321'u
  #var rng = newRngField(lo, RngMilc6, seed)
  var rng = newRngField(lo, MRG32k3a, seed)
  g.gaussian rng
  toc "gaussian"
  #var r0 = lo.Real()
  #var cv1 = lo.ColorVector()
  #var cv2 = lo.ColorVector()
  echo 6.0 * g.plaq
  toc "plaq"
  #cv0.gaussian rng
  resetTimers()

  proc testPlaq =
    tic("testPlaq")
    #let hl = lo.makeHaloLayout([1,1,1,1],[1,1,1,1])
    let hl = lo.makeHaloLayout([1,1,1,1],[0,0,0,0])
    toc "makeHaloLayout"
    type HM = HaloMap[type lo]
    let comm {.cursor.} = getDefaultComm()
    var hm = newSeq[HM](nd)
    for mu in 0..<nd:
      var offsets = newSeq[seq[int32]](0)
      for nu in 0..<nd:
        if nu == mu: continue
        var t = newSeq[int32](nd)
        t[nu] = 1
        offsets.add t
      hm[mu] = hl.makeHaloMap(comm, offsets)
    toc "makeHaloMap"
    type H = type makeHalo(hl, g[0])
    var h = newSeq[H](nd)
    for d in 0..<nd:
      h[d] = makeHalo(hl, g[d])
    toc "makeHalo"
    for d in 0..<nd:
      h[d].update hm[d], comm
    toc "update"
    type PT = evalType(g[0][0][0,0].re)
    var pl = newSeq[float](6)
    threads:
      var spl: array[6,PT]
      for i in g[0]:
        var k = 0
        for mu in 1..<4:
          let n0 = hl.neighborFwd[mu][i]
          for nu in 0..<mu:
            let n1 = hl.neighborFwd[nu][i]
            let a = g[mu][i] * h[nu][n0]
            let b = g[nu][i] * h[mu][n1]
            let p = redot(a, b)
            spl[k] += p
            inc k
      for k in 0..<6:
        let tpl = simdReduce spl[k]
        threadCritical:
          pl[k] += tpl[k]
      #var tpl: array[6,float]
      #for k in 0..<6: tpl[k] = simdReduce spl[k]
      #threadSum tpl
      #threadSingle:
      #    for k in 0..<6: pl[k] = tpl[k]
    rankSum pl
    toc "pl"
    let vf = 1.0/(g[0][0].nRows*lo.physVol)
    echo pl * vf
    echo 6.0 * g.plaq
  testPlaq()

  proc testStaple =
    var st = lo.newGauge
    var gf = lo.newGauge
    for mu in 0..<nd: st[mu] := 0
    tic("testStaple")
    let hl = lo.makeHaloLayout([1,1,1,1],[1,1,1,1])
    toc "makeHaloLayout"
    type HM = HaloMap[type lo]
    let comm {.cursor.} = getDefaultComm()
    var hm = newSeq[HM](nd)
    for mu in 0..<nd:
      var offsets = newSeq[seq[int32]](0)
      for nu in 0..<nd:
        var t = newSeq[int32](nd)
        t[nu] = -1
        offsets.add @t
        if nu == mu: continue
        t[nu] = 1
        offsets.add @t
        t[mu] = -1
        offsets.add @t
      hm[mu] = hl.makeHaloMap(comm, offsets)
      #echo offsets
    toc "makeHaloMap"
    type H = type makeHalo(hl, g[0])
    var h = newSeq[H](nd)
    for mu in 0..<nd:
      h[mu] = makeHalo(hl, g[mu])
    toc "makeHalo"
    for mu in 0..<nd:
      h[mu].update hm[mu], comm
    toc "update"
    var cp = 1.0/3.0
    threads:
      for i in g[0]:
        for mu in 1..<4:
          let fmu = hl.neighborFwd[mu][i]
          let bmu = hl.neighborBck[mu][i]
          for nu in 0..<mu:
            let fnu = hl.neighborFwd[nu][i]
            let bnu = hl.neighborBck[nu][i]
            let fmubnu = hl.neighborBck[nu][fmu]
            let fnubmu = hl.neighborFwd[nu][bmu]
            let t0 = h[mu][fnu] * h[nu][fmu].adj
            let t1 = g[nu][i] * t0
            st[mu][i] += cp*t1
            let t2 = g[mu][i] * t0.adj
            st[nu][i] += cp*t2
            let t3 = h[mu][bnu] * h[nu][fmubnu]
            let t4 = h[nu][bnu].adj * t3
            st[mu][i] += cp*t4
            let t5 = h[nu][bmu] * h[mu][fnubmu]
            let t6 = h[mu][bmu].adj * t5
            st[nu][i] += cp*t6
    toc "st"
    GaugeActionCoeffs(plaq: 1.0).gaugeForce(g, gf)
    toc "gaugeForce"
    for mu in 0..<nd:
      contractProjectTAH(st[mu], g[mu], st[mu])
      echo gf[mu].norm2, " ", st[mu].norm2, " ", (gf[mu]-st[mu]).norm2
  testStaple()

  echoTimers()
  qexFinalize()


# u[mu=0] -> [0,1,0,0],[0,0,1,0],[0,0,0,1],
#            [-1,0,0,0],[0,-1,0,0],[0,0,-1,0],[0,0,0,-1],
#            [-1,1,0,0],[-1,0,1,0],[-1,0,0,1]

# [1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],
# [-1,0,0,0],[0,-1,0,0],[0,0,-1,0],[0,0,0,-1],
# [-1,1,0,0],[-1,0,1,0],[-1,0,0,1],
# [1,-1,0,0],[0,-1,1,0],[0,-1,0,1],
# [1,0,-1,0],[0,1,-1,0],[0,0,-1,1],
# [1,0,0,-1],[0,1,0,-1],[0,0,1,-1]
