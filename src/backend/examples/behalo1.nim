import qex
import physics/qcdTypes
import backend/[accel,cpugpu,cgfield]
import bench/commonBench
#import strformat
import macros
import base/metaUtils
import parseUtils
import sequtils, strutils
import comms/halo

type
  gpuHaloLayout*[V:static int] = object
    #lo*: L  # layout
    #outerExt*: seq[int32]  # extended outer lattice size
    #offset*: seq[int32]  # offset of outer lattice in extended outer
    #lex*: seq[int32]  # outerExt lex index of extended sites
    #index*: seq[int32]  # index of extended site for given lex index
    #neighborFwd*: gpuSeq[gpuSeq[int32]]  # fwd neighbor for extended outer lattice
    #neighborBck*: gpuSeq[gpuSeq[int32]]  # bck neighbor for extended outer lattice
    neighborFwd*: gpuSeq[ptr UncheckedArray[int32]]  # fwd neighbor for extended outer lattice
    neighborBck*: gpuSeq[ptr UncheckedArray[int32]]  # bck neighbor for extended outer lattice
    nOut*: int  # sites in outer lattice
    nExt*: int  # sites in extended outer lattice
    #nOutPar*: array[2,int]  # sites in outer lattice by parity
    #nExtPar*: array[2,int]  # sites in extended outer lattice by parity

proc nbrFwd[T](ghl: gpuHaloLayout, mu: int, s: T): T =
  let i = (s.V*s[]) div ghl.V
  let j = (s.V*s[]) mod ghl.V
  let n0 = ghl.neighborFwd[mu][i]
  let s0 = T((n0*ghl.V+j)div(s.V))
  result = s0

proc nbrBck[T](ghl: gpuHaloLayout, mu: int, s: T): T =
  let i = (s.V*s[]) div ghl.V
  let j = (s.V*s[]) mod ghl.V
  let n0 = ghl.neighborBck[mu][i]
  let s0 = T((n0*ghl.V+j)div(s.V))
  result = s0

proc getGhl(V:static int): ptr gpuHaloLayout[V] =  # FIXME: need to make pointer map
  var ghl {.global.}: gpuHaloLayout[V]
  result = addr ghl

#proc checkNeighbors(g: gpuHaloLayout) =
#  for i in 0..<g.neighborFwd.n:
#    for j in 0..<g.neighborFwd[i].n:
#      if g.neighborFwd[i][j] >= g.nExt: echo i, " ", j, " ", g.neighborFwd[i][j]
#  for i in 0..<g.neighborBck.n:
#    for j in 0..<g.neighborBck[i].n:
#      if g.neighborBck[i][j] >= g.nExt: echo i, " ", j, " ", g.neighborBck[i][j]

proc toGpu*(x: HaloLayout): auto =  # TODO make CPU not copy inner seq
  tic("toGpuHaloLayout")
  var g: gpuHaloLayout[x.lo.V]
  g.nOut = x.nOut
  g.nExt = x.nExt
  g.neighborFwd.newGpuSeq x.neighborFwd.len
  for i in 0..<g.neighborFwd.n:
    #g.neighborFwd[i].newGpuSeq x.neighborFwd[i].len
    #for j in 0..<g.neighborFwd[i].n:
    #  g.neighborFwd[i][j] = x.neighborFwd[i][j]
    let n = x.neighborFwd[i].len
    let b = n*sizeof(int32)
    #g.neighborFwd[i] = cast[typeof g.neighborFwd[i]](gpuMalloc(b))
    let p = gpuMalloc(b)
    gpuMemCpyToGpu(p, addr x.neighborFwd[i][0], b)  # copy data to GPU
    gpuMemCpyToGpu(addr g.neighborFwd[i], addr p, sizeof(p))  # copy ptr to GPU
  toc("neighborFwd")
  g.neighborBck.newGpuSeq x.neighborBck.len
  for i in 0..<g.neighborBck.n:
    #g.neighborBck[i].newGpuSeq x.neighborBck[i].len
    #for j in 0..<g.neighborBck[i].n:
    #  g.neighborBck[i][j] = x.neighborBck[i][j]
    let n = x.neighborBck[i].len
    let b = n*sizeof(int32)
    #g.neighborBck[i] = cast[typeof g.neighborBck[i]](gpuMalloc(b))
    let p = gpuMalloc(b)
    gpuMemCpyToGpu(p, addr x.neighborBck[i][0], b)  # copy data to GPU
    gpuMemCpyToGpu(addr g.neighborBck[i], addr p, sizeof(p))  # copy ptr to GPU
  #g.checkNeighbors
  toc("neighborBck")
  getGhl(x.lo.V)[] = g
  g

template getGpu*(x: HaloLayout, g: gpuHaloLayout): auto = g

template fromGpu*(x: HaloLayout, g: gpuHaloLayout) = discard

type
  gpuHalo*[V:static int,F,T] = object
    layout*: gpuHaloLayout[V]
    field*: F
    halo*: gpuSeq[T]
    nOut*: int  # sites in outer lattice
    nExt*: int  # sites in extended outer lattice
    #nOutPar*: array[2,int]  # sites in outer lattice by parity
    #nExtPar*: array[2,int]  # sites in extended outer lattice by parity

proc indexPtr*[V:static int,F,T](h: gpuHalo[V,F,T], i: SomeInteger): ptr T =
  #doAssert(i>=0)
  #if i>=h.nExt: echo "i: ", i, "  nExt: ", h.nExt
  #doAssert(i<h.nExt)
  let k = i - h.nOut
  result = if k<0: addr h.field.p[i] else: addr h.halo[k]
template `[]`*(h: gpuHalo, i: SomeInteger): auto = indexPtr(h,i)[]
template `[]`*[F,T;VV,L:static int](h: gpuHalo[VV,F,T], i: SiteV[L]): auto =
  when F.V == L:
    h[i[]]
  else:
    let s = i[] div F.V
    let v = i[] mod F.V
    h[s][asSimd(v)]

proc `[]=`*(h: gpuHalo, i: SomeInteger, x: auto) =
  let k = i - h.nOut
  if k < 0:
    h.field[i] = x
  else:
    h.halo[k] := x

template gpuType*[L,F,T](x: typedesc[Halo[L,F,T]]): typedesc =
  gpuHalo[L.V, gpuType F, gpuType T]

proc toGpu*[L,F,T](x: Halo[L,F,T]): auto {.noInit.} =
  tic("toGpuHalo")
  var g {.noInit.}: gpuHalo[L.V, gpuType F, gpuType T]
  g.layout = getGhl(x.layout.L.V)[]
  g.field = toGpu(x.field)
  toc("toGpuField")
  g.nOut = x.nOut
  g.nExt = x.nExt
  g.halo.newGpuSeq x.halo.len
  toc("newGpuSeq")
  gpuMemCpyToGPU(g.halo.p, addr x.halo[0], g.halo.bytes)
  toc("gpuMemCpyToGPU")
  #ghl.checkNeighbors
  #toc("end")
  g

proc fromGpu*(x: var Halo, g: gpuHalo) =
  tic("fromGpuHalo")
  gpuMemCpyToCPU(addr x.halo[0], g.halo.p, g.halo.bytes)
  toc("end")

proc toGpu*(g: var gpuSeq[gpuHalo], x: seq[Halo]) =
  tic("toGpuSeqHalo")
  for i in 0..<g.n:
    #g[i] = toGpu(x[i])
    let t = toGpu(x[i])
    gpuMemCpyToGpu(addr g[i], addr t, sizeof(t))
  toc("end")

proc copyFromGpu*(x: var seq[Halo], g: gpuSeq[gpuHalo]) =
  tic("fromGpuSeqHalo")
  for i in 0..<x.len:
    var t {.noInit.}: typeof g[i]
    gpuMemCpyToCpu(addr t, addr g[i], sizeof(t))
    x[i].fromGpu(t)
  toc("end")

proc testPlaq(g:auto) =
  tic "testPlaq"
  let lo = g[0].l
  let nd = lo.nDim
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
  #var p = newSeq[typeof g[0]](6)
  #for i in 0..<6:
  #  p[i] = g[0].newOneOf
  #threads:
  #  for i in 0..<6:
  #    p[i] := 0
  var pl = newSeq[float](6)
  var gs = newGpuSum[array[6,float]](lo.nSites)
  toc "create fields"
  proc gpuSite(x: GpuField): auto =
    for i in gpuSites(x): return x[i]
  for nreps in [2,10]:
    resetTimers()
    for rep in 0..<nreps:
      tic "rep"
      for d in 0..<nd:
        h[d].update hm[d], comm
      toc "update"
      when false:
        threads:
          for i in g[0]:
            var k = 0
            for mu in 1..<4:
              let n0 = hl.neighborFwd[mu][i]
              for nu in 0..<mu:
                let n1 = hl.neighborFwd[nu][i]
                let a = g[mu][i] * h[nu][n0]
                let b = g[nu][i] * h[mu][n1]
                p[k][i] += a.adj * b
                inc k
      else:
        onGpu(lo):
          template g(i:int):auto = h[i].field
          var tpl: array[6,float]
          #var tpl: array[6,typeof redot(g(0).gpuSite,g(0).gpuSite)]
          for s in gpuSites(g(0)):
            var k = 0
            for mu in 1..<4:
              let smu = hl.nbrFwd(mu, s)
              for nu in 0..<mu:
                let snu = hl.nbrFwd(nu, s)
                let a = g(mu)[s] * h[nu][smu]
                let b = g(nu)[s] * h[mu][snu]
                tpl[k] += redot(a, b).simdSum
                #tpl[k] += redot(a, b)
                inc k
          gs.reduce tpl
          #var tplf: array[6,float]
          #for k in 0..<6: tplf[k] = tpl[k].simdSum
          #gs.reduce tplf
      let nc = g[0][0].getNc
      let mm = nc*nc*(8*nc-2)
      let rd = 8*nc*nc
      toc("plaq",flops=lo.nSites*6*(2*mm+rd))
  #toc "p"
  #threads:
  for k in 0..<6:
    #pl[k] = p[k].trace.re
    pl[k] = gs.value[k]
  #toc "pl"
  rankSum pl
  let vf = 1.0/(g[0][0].nRows*lo.physVol)
  let ph = pl * vf
  let pp = 6.0 * g.plaq
  let d = ph - pp
  echo ph
  echo d
  echo "norm2 diff: ", sum(d*d)
  #echo pl * vf
  #echo 6.0 * g.plaq

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
  #g.unit
  toc "gaussian"
  #var r0 = lo.Real()
  #var cv1 = lo.ColorVector()
  #var cv2 = lo.ColorVector()
  echo 6.0 * g.plaq
  toc "plaq"
  #cv0.gaussian rng
  resetTimers()
  testPlaq(g)
  echoProf()
  qexFinalize()
