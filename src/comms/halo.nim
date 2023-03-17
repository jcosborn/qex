import qex, base/hyper, comms/gather
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
  Halo*[L,F,T] = ref object
    layout*: HaloLayout[L]
    field*: F
    halo*: alignedMem[T]
    nOut*: int  # sites in outer lattice
    nExt*: int  # sites in extended outer lattice
  HaloMap*[L] = ref object
    layout*: HaloLayout[L]
    offsets*: seq[seq[int32]]  # offsets of outer sites needed
    gather*: GatherMap  # gather map to fetch needed halo sites
    gatherRev*: GatherMap  # gather map to send halo sites back
  GatherHalo*[F,T;Rev:static bool] = object
    src: F
    dest: ptr UncheckedArray[T]
    elemsize: int
    vlen: int

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
  let nd = lo.nDim
  var x = newSeq[int32](nd)
  var k = int32 nOut
  for i in 0..<nExt:  # loop over lex indices
    x.lexCoord(i, outerExt)
    if x.inside(bckOffset, outerHigh):
      var y = x - bckOffset
      let loc = lo.localIndex(y) div L.V
      lex[loc] = int32 i
      index[i] = int32 loc
    else:
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
  for i in hl.nOut ..< hl.nExt:  # loop over halo sites
    let li = hl.lex[i]
    x.lexCoord(li, hl.outerExt)
    var keep = false
    for j in 0 ..< offsets.len:
      let y = x - offsets[j]
      if y.inside(hl.offset, outerHigh):
        keep = true
        break
    if keep:
      for k in 0 ..< vlen:
        var y = x + vecOffset[k]
        for l in 0..<nd:
          if y[l] < 0: y[l] += int32 hl.lo.physGeom[l]
          if y[l] >= hl.lo.physGeom[l]: y[l] -= int32 hl.lo.physGeom[l]
        let ri = hl.lo.rankIndex(y)
        let didx = vlen*(i-hl.nOut)+k
        rl.add RecvList(didx: int32 didx, srank: int32 ri.rank, sidx: int32 ri.index)
        sl.add SendList(sidx: int32 didx, drank: int32 ri.rank, didx: int32 ri.index)
  toc("RecvList")
  result.gather = c.makeGatherMap(rl)
  result.gatherRev = c.makeGatherMap(sl)
  toc("makeGatherMap")
proc makeHaloMap*[L,N](hl: HaloLayout[L], c: Comm, offsets: seq[array[N,SomeInteger]]): HaloMap[L] =
  var o = newSeq[seq[int32]](offsets.len)
  for i in 0..<offsets.len:
    let t = offsets[i]
    o[i] = @t
  makeHaloMap(hl, c, o)

# T is vectorized type
proc makeHalo*[L,F,T](hl: HaloLayout[L], f: F, t: typedesc[T]): Halo[L,F,T] =
  result.new
  result.layout = hl
  result.field = f
  result.nOut = hl.nOut
  result.nExt = hl.nExt
  let nhalo = result.nExt - result.nOut
  result.halo.newU(nhalo)

template makeHalo*[L,F](hl: HaloLayout[L], f: F): auto =
  makeHalo(hl, f, eval(F.type[0]))

template copy*[F,T;Rev:static bool](gh: GatherHalo[F,T,Rev], d: pointer, s: SomeInteger) =
  type E = eval(index(type T, type asSimd(0)))
  let p = cast[ptr E](d)
  when Rev:
    let o = s div gh.vlen
    let i = s mod gh.vlen
    p[] := gh.dest[o][asSimd(i)]
  else:
    p[] := gh.src{s}
template copy*[F,T;Rev:static bool](gh: GatherHalo[F,T,Rev], d: SomeInteger, s: pointer) =
  type E = eval(index(type T,type asSimd(0)))
  let p = cast[ptr E](s)
  when Rev:
    #threadCritical:
    gh.src{d} += p[]
  else:
    let o = d div gh.vlen
    let i = d mod gh.vlen
    gh.dest[o][asSimd(i)] = p[]
template copy*[F,T;Rev:static bool](gh: GatherHalo[F,T,Rev], d: SomeInteger, s: SomeInteger) =
  when Rev:
    let o = s div gh.vlen
    let i = s mod gh.vlen
    #threadCritical:
    gh.src{d} += gh.dest[o][asSimd(i)]
  else:
    let o = d div gh.vlen
    let i = d mod gh.vlen
    gh.dest[o][asSimd(i)] = gh.src{s}
#[
template copy*[F,T](gh: GatherHalo[F,T], dv: SomeInteger, di: array,
                    s: array, n: SomeInteger) =
  echo n, " ", dv, " ", di, " ", s
  #for i in 0..<n:
  #  gh.dest[dv][asSimd(di[i])] = gh.src{s[i]}
  let t = gh.dest[dv]
  for i in 0..<n:
    t[asSimd(di[i])] = gh.src{s[i]}
  gh.dest[dv] = t
]#

proc update*[L,F,T](h: Halo[L,F,T], hm: HaloMap[L], c: Comm) =
  tic("Halo update")
  let elemSize = sizeof(T) div L.V
  var gh: GatherHalo[F,T,false]
  gh.src = h.field
  gh.dest = cast[ptr UncheckedArray[T]](unsafeAddr h.halo[0])
  gh.elemsize = elemSize
  gh.vlen = L.V
  c.gather(hm.gather, gh)
  toc("gather")

proc updateRev*[L,F,T](h: Halo[L,F,T], hm: HaloMap[L], c: Comm) =
  tic("Halo update")
  let elemSize = sizeof(T) div L.V
  var gh: GatherHalo[F,T,true]
  gh.src = h.field
  gh.dest = cast[ptr UncheckedArray[T]](unsafeAddr h.halo[0])
  gh.elemsize = elemSize
  gh.vlen = L.V
  #let g = hm.gather.reverse
  c.gather(hm.gatherRev, gh)
  toc("gather")

#[
proc `[]`*(h: Halo, i: SomeInteger): auto {.inline,noInit.} =
  let k = i - h.nOut
  if k < 0:
    h.field[i]
  else:
    h.halo[k]

proc `[]`*(h: var Halo, i: SomeInteger): var auto {.inline.} =
  var t {.noInit.}: ptr type(h.field[i])
  let k = i - h.nOut
  if k < 0:
    t = addr h.field[i]
  else:
    t = addr h.halo[k]
  t[]
template `[]`*(h: Halo, i: SomeInteger): auto =
  var t {.noInit.}: ptr type(h.field[i])
  let k = i - h.nOut
  if k < 0:
    t = addr h.field[i]
  else:
    t = addr h.halo[k]
  t[]
]#

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

proc neighbor*(h: Halo, i: SomeInteger, mu: SomeInteger, fb: SomeInteger): int32 =
  if fb > 0:
    h.map.neighborFwd[mu][i]
  else:
    h.map.neighborBck[mu][i]

when isMainModule:
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
    let comm = getDefaultComm()
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
    let comm = getDefaultComm()
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
    gf.gaugeForce(g)
    #gf.gaugeForce2(g)
    #gf.gaugeForce3(g)
    toc "gaugeForce"
    for mu in 0..<nd:
      for e in st[mu]:
        let s = g[mu][e] * st[mu][e].adj
        st[mu][e].projectTAH s
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

