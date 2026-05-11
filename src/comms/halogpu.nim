import qex
import physics/qcdTypes
import backend/[accel,cpugpu,cgfield]
import bench/commonBench
import parseUtils
import sequtils, strutils
import comms/halo

type
  GpuHaloLayout*[V:static int] = object
    #lo*: L  # layout
    #outerExt*: seq[int32]  # extended outer lattice size
    #offset*: seq[int32]  # offset of outer lattice in extended outer
    #lex*: seq[int32]  # outerExt lex index of extended sites
    #index*: seq[int32]  # index of extended site for given lex index
    #neighborFwd*: gpuSeq[gpuSeq[int32]]  # fwd neighbor for extended outer lattice
    #neighborBck*: gpuSeq[gpuSeq[int32]]  # bck neighbor for extended outer lattice
    neighborFwd*: GpuSeq[ptr UncheckedArray[int32]]  # fwd neighbor for extended outer lattice
    neighborBck*: GpuSeq[ptr UncheckedArray[int32]]  # bck neighbor for extended outer lattice
    nOut*: int  # sites in outer lattice
    nExt*: int  # sites in extended outer lattice
    #nOutPar*: array[2,int]  # sites in outer lattice by parity
    #nExtPar*: array[2,int]  # sites in extended outer lattice by parity

proc nbrFwd*[T](ghl: GpuHaloLayout, mu: int, s: T): T =
  let i = (s.V*s[]) div ghl.V
  let j = (s.V*s[]) mod ghl.V
  let n0 = ghl.neighborFwd[mu][i]
  let s0 = T((n0*ghl.V+j)div(s.V))
  result = s0

proc nbrBck*[T](ghl: GpuHaloLayout, mu: int, s: T): T =
  let i = (s.V*s[]) div ghl.V
  let j = (s.V*s[]) mod ghl.V
  let n0 = ghl.neighborBck[mu][i]
  let s0 = T((n0*ghl.V+j)div(s.V))
  result = s0

proc setRO(x: seq[seq]) =
  gpuMemFlagsExcl(addr x[0], {gmCpuWrite,gmGpuWrite}) # set read only
  for i in 0..<x.len:
    gpuMemFlagsExcl(addr x[i][0], {gmCpuWrite,gmGpuWrite}) # set read only

proc toGpu*(g: var GpuHaloLayout, c: HaloLayout) =  # TODO implement copyIn once
  tic("toGpuHaloLayout")
  g.nOut = c.nOut
  g.nExt = c.nExt
  pushGpuMemTag("nbrFwd")
  g.neighborFwd.toGpu(c.neighborFwd)
  setRO(c.neighborFwd)
  popGpuMemTag()
  toc("nbrFwd")
  pushGpuMemTag("nbrBck")
  g.neighborBck.toGpu(c.neighborBck)
  setRO(c.neighborBck)
  popGpuMemTag()
  toc("nbrBck")

proc toGpu*(c: HaloLayout): auto =
  var g: GpuHaloLayout[c.lo.V]
  g.toGpu(c)
  g

template getGpu*(c: HaloLayout, g: GpuHaloLayout): auto = g

proc fromGpu*(c: HaloLayout, g: GpuHaloLayout) =
  tic("fromGpuHaloLayout")
  c.neighborFwd.fromGpu(g.neighborFwd)
  toc("nbrFwd")
  c.neighborBck.fromGpu(g.neighborBck)
  toc("nbrBck")

proc freeGpuMem*(c: HaloLayout) =
  freeGpuMem addr c.neighborFwd[0]
  for i in 0..<c.neighborFwd.len:
    freeGpuMem addr c.neighborFwd[i][0]
  freeGpuMem addr c.neighborBck[0]
  for i in 0..<c.neighborBck.len:
    freeGpuMem addr c.neighborBck[i][0]

type
  GpuHalo*[V:static int,F,T] = object
    layout*: GpuHaloLayout[V]
    field*: F
    halo*: GpuSeq[T]
    nOut*: int  # sites in outer lattice
    nExt*: int  # sites in extended outer lattice
    #nOutPar*: array[2,int]  # sites in outer lattice by parity
    #nExtPar*: array[2,int]  # sites in extended outer lattice by parity
proc displayName*(x: typedesc[GpuHalo]): string =
  result = "GpuHalo"

proc indexPtr*[V:static int,F,T](h: GpuHalo[V,F,T], i: SomeInteger): ptr T =
  #doAssert(i>=0)
  #if i>=h.nExt: echo "i: ", i, "  nExt: ", h.nExt
  #doAssert(i<h.nExt)
  let k = i - h.nOut
  result = if k<0: addr h.field.p[i] else: addr h.halo[k]
template `[]`*(h: GpuHalo, i: SomeInteger): auto = indexPtr(h,i)[]
template `[]`*[F,T;VV,L:static int](h: GpuHalo[VV,F,T], i: SiteV[L]): auto =
  when F.V == L:
    h[i[]]
  else:
    let s = i[] div F.V
    let v = i[] mod F.V
    h[s][asSimd(v)]

proc `[]=`*(h: GpuHalo, i: SomeInteger, x: auto) =
  let k = i - h.nOut
  if k < 0:
    h.field[i] = x
  else:
    h.halo[k] := x

template gpuType*[L,F,T](c: typedesc[Halo[L,F,T]]): typedesc =
  GpuHalo[L.V, gpuType F, gpuType T]

proc gpuFlagsExcl*(x: Halo, f: set[gmFlags]) =
  gpuMemFlagsExcl(addr x.halo[0], f)
proc gpuFlagsIncl*(x: Halo, f: set[gmFlags]) =
  gpuMemFlagsIncl(addr x.halo[0], f)

proc toGpu*(g: var GpuHalo, c: Halo) =
  tic("toGpuHalo")
  g.nOut = c.nOut
  g.nExt = c.nExt
  pushGpuMemTag("Halo")
  g.layout.toGpu(c.layout)
  toc("Layout")
  g.field.toGpu(c.field)
  toc("Field")
  g.halo.toGpu(c.halo)
  popGpuMemTag()
  toc("Halo",flops=g.halo.bytes)

proc toGpu*[L,F,T](c: Halo[L,F,T]): auto {.noInit.} =
  var g {.noInit.}: GpuHalo[L.V, gpuType F, gpuType T]
  g.toGpu(c)
  g

template getGpu*(c: Halo, g: GpuHalo): auto = g

proc fromGpu*(c: var Halo, g: GpuHalo) =
  tic("fromGpuHalo")
  c.field.fromGpu(g.field)
  toc("Field")
  c.halo.fromGpu(g.halo)
  toc("Halo")

proc fromGpu*(c: var Halo) =
  tic("fromGpuHalo")
  c.field.fromGpu()
  toc("Field")
  c.halo.fromGpu()
  toc("Halo")

proc toGpu*(g: var GpuSeq[GpuHalo], c: seq[Halo], pgm: ptr GpuMem) =
  if pgm.needsCopyIn:
    tic("toGpuSeqHaloCopyIn")
    var t = newSeq[typeof g[0]](g.n)
    for i in 0..<g.n:
      t[i].toGpu(c[i])
    toc("loopToGpu")
    pgm.copyIn(addr t[0])
    pgm.flags.excl {gmCpuWrite,gmGpuWrite}  # set seq container read only (not halo data)
    toc("copyIn")
  else:
    tic("toGpuSeqHalo")
    for i in 0..<g.n:
      discard toGpu(c[i])
    toc("loopToGpu")

proc fromGpu*(c: var seq[Halo], g: GpuSeq[GpuHalo], pgm: ptr GpuMem) =
  if pgm.needsCopyOut:
    tic("fromGpuSeqHaloCopyOut")
    var t = newSeq[typeof g[0]](g.n)
    pgm.copyOut(addr t[0])
    toc("copyOut")
    for i in 0..<g.n:
      c[i].fromGpu(t[i])
    toc("loopToGpu")
  else:
    tic("fromGpuSeqHalo")
    for i in 0..<g.n:
      c[i].fromGpu()
    toc("loopFromGpu")

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
    h[d].gpuFlagsExcl {gmGpuWrite}
    h[d].gpuFlagsIncl {gmCpuWriteOnce}
    g[d].gpuFlagsExcl {gmGpuWrite}
    g[d].gpuFlagsIncl {gmCpuWriteOnce}
  toc "makeHalo"
  #var p = newSeq[typeof g[0]](6)
  #for i in 0..<6:
  #  p[i] = g[0].newOneOf
  #threads:
  #  for i in 0..<6:
  #    p[i] := 0
  var pl = newSeq[float](6)
  var gs = newGpuSum[array[6,float]](lo.nSites)
  pushGpuMemTag("testPlaq")
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
                #let a = h[1][s]
                #let b = h[1][s]
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
  echo dumpGpuMem()

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
