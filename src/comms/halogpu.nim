import qex
import physics/qcdTypes
import backend/[accel,cpugpu,cgfield]
import bench/commonBench
import parseUtils
import std/[nativesockets, hashes]
import sequtils, strutils
import comms/[halo,gather,qmp,commsQmp]

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
  toc("Halo")

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

proc toDevice*[T](x: seq[T]): ptr UncheckedArray[T] =
  ## Device copy of x, nil if x is empty.
  if x.len > 0:
    result = cast[ptr UncheckedArray[T]](gpuMalloc(x.len*sizeof(T)))
    gpuMemCpyToGpu(result, unsafeAddr x[0], x.len*sizeof(T))

proc haloSource*[L](hl: HaloLayout[L], gm: GatherMap): seq[int32] =
  ## For each halo lane V*(i-nOut)+l filled by gm: the local site copied into
  ## it, or V*nOut + its position in the receive buffer; -1 for lanes gm
  ## leaves out.
  result = newSeq[int32](L.V*(hl.nExt - hl.nOut))
  for i in 0..<result.len: result[i] = -1
  for k in 0..<gm.ldest.len: result[gm.ldest[k]] = gm.lidx[k]
  for k in 0..<gm.rdest.len: result[gm.rdest[k]] = int32(L.V*hl.nOut + k)

{.passL: "-lze_loader".}
{.emit: """/*INCLUDESECTION*/
#include <omp.h>
""".}
type
  ZeIpcMemHandle {.importc: "ze_ipc_mem_handle_t", header: "level_zero/ze_api.h".} = object
    data: array[64, char]
proc zeMemGetIpcHandle(ctx, p: pointer; h: var ZeIpcMemHandle): cint {.importc, header: "level_zero/ze_api.h".}
proc zeMemOpenIpcHandle(ctx, dev: pointer; h: ZeIpcMemHandle; flags: uint32; p: var pointer): cint {.importc, header: "level_zero/ze_api.h".}
proc zeMemCloseIpcHandle(ctx, p: pointer): cint {.importc, header: "level_zero/ze_api.h".}
proc zeMemGetAddressRange(ctx, p: pointer; base: var pointer; bytes: var csize_t): cint {.importc, header: "level_zero/ze_api.h".}
proc syscall(n: clong): clong {.importc, header: "<unistd.h>", varargs.}
proc getpid(): cint {.importc, header: "<unistd.h>".}
var SYS_pidfd_open {.importc, header: "<sys/syscall.h>".}: clong
var SYS_pidfd_getfd {.importc, header: "<sys/syscall.h>".}: clong

proc zeContext(): tuple[ctx, dev: pointer] =
  ## Level Zero context and device of the default OpenMP device.
  var ctx, dev: pointer
  {.emit: """
  omp_interop_t obj = omp_interop_none;
  int dn = omp_get_default_device(), err;
  #pragma omp interop init(targetsync: obj) device(dn)
  `ctx` = omp_get_interop_ptr(obj, omp_ipr_device_context, &err);
  `dev` = omp_get_interop_ptr(obj, omp_ipr_device, &err);
  #pragma omp interop destroy(obj)
  """.}
  (ctx, dev)

var haloIpc* = true  ## peers on the same host store into each other's device memory

type
  GpuHaloEx*[T] = ref object
    ## Halo exchange of a field in device memory in the qex SIMD layout: real
    ## c of site i at ((i div v)*ne + c)*v + i mod v for v lanes.  The
    ## buffers hold message m at ne*m.start, with
    ## component c of its k-th site at ne*m.start + c*m.count + k.  Send
    ## slot k takes component c at sdst[k][c*sstr[k]]: in sbuf, or, for peers
    ## on the same host, directly in their rbuf through a Level Zero IPC
    ## mapping, announced with a one byte message.  A kernel producing the
    ## field stores site i in its slots sslot[s*n + i], s < nslot, with
    ## sendSite; pack does it for a whole field.  After wait, receive position
    ## p has component c at rbuf[rofs[p] + c*rstr[p]], see recvSite.  Local
    ## halo sites are left to the caller, see haloSource.
    ## A peer rewrites rbuf in its next exchange; callers order exchanges of
    ## the same object with a global sum or barrier in between.
    ne*, n*, v*: int  # reals per site, local sites, lanes
    nsend*, nrecv*, nslot*: int
    sidx*: ptr UncheckedArray[int32]  # local site of each send slot
    sslot*: ptr UncheckedArray[int32]  # [nslot][n]: send slots of each site, -1 for none
    sdst*: ptr UncheckedArray[ptr UncheckedArray[T]]  # component 0 of each send slot
    sstr*: ptr UncheckedArray[int32]  # component stride of each send slot
    rofs*, rstr*: ptr UncheckedArray[int32]  # component 0 and stride of each receive position
    sbuf*, rbuf*: ptr UncheckedArray[T]
    smsg*, rmsg*: seq[MsgInfo]
    speer*, rpeer*: seq[bool]  # messages stored directly in the peer rbuf
    peers: seq[pointer]  # opened IPC mappings
    flags: seq[char]
    mems: seq[QMP_msgmem_t]
    msg: QMP_msghandle_t  # all receives and sends, declared once

proc newGpuHaloEx*[T](gm: GatherMap, ne, n, v: int, c: Comm): GpuHaloEx[T] =
  type Info = object
    pid, fd, off: int  # off: bytes from the exported allocation to the message
    h: ZeIpcMemHandle
  tic("newGpuHaloEx")
  var ex = GpuHaloEx[T](ne: ne, n: n, v: v, nsend: gm.sidx.len, nrecv: gm.rdest.len)
  ex.sidx = gm.sidx.toDevice
  if ex.nsend > 0:
    ex.sbuf = cast[ptr UncheckedArray[T]](gpuMalloc(ne*ex.nsend*sizeof(T)))
  if ex.nrecv > 0:
    ex.rbuf = cast[ptr UncheckedArray[T]](gpuMalloc(ne*ex.nrecv*sizeof(T)))
  ex.smsg = gm.smsginfo
  ex.rmsg = gm.rmsginfo
  ex.flags.newSeq(1 + ex.rmsg.len)
  var cnt = newSeq[int32](n)
  for j in gm.sidx: inc cnt[j]
  for j in 0..<n: ex.nslot = max(ex.nslot, int cnt[j])
  var sl = newSeq[int32](ex.nslot*n)
  for k in 0..<sl.len: sl[k] = -1
  for j in 0..<n: cnt[j] = 0
  for k, j in gm.sidx:
    sl[cnt[j]*n + j] = int32 k
    inc cnt[j]
  ex.sslot = sl.toDevice
  var ro = newSeq[int32](ex.nrecv)
  var rs = newSeq[int32](ex.nrecv)
  for m in ex.rmsg:
    for k in 0..<m.count:
      ro[m.start+k] = int32(ne*m.start + k)
      rs[m.start+k] = int32 m.count
  ex.rofs = ro.toDevice
  ex.rstr = rs.toDevice
  var host = newSeq[float](c.size)  # host name hash of each rank
  host[c.rank] = float(hash(getHostname()) and 0xffffff)
  c.allReduce(addr host[0], c.size)
  ex.speer.newSeq(ex.smsg.len)
  ex.rpeer.newSeq(ex.rmsg.len)
  for i, m in ex.smsg: ex.speer[i] = haloIpc and host[m.rank] == host[c.rank]
  for i, m in ex.rmsg: ex.rpeer[i] = haloIpc and host[m.rank] == host[c.rank]
  toc("hosts")
  let (ctx, dev) = zeContext()
  var rout = newSeq[Info](ex.rmsg.len)
  var sin = newSeq[Info](ex.smsg.len)
  var ns, nr = 0
  for i, m in ex.rmsg:
    if ex.rpeer[i]:
      # the handle maps the whole allocation holding rbuf, e.g. a pool block
      var base: pointer
      var bytes: csize_t
      if zeMemGetAddressRange(ctx, ex.rbuf, base, bytes) != 0 or
         zeMemGetIpcHandle(ctx, base, rout[i].h) != 0:
        qexError("zeMemGetIpcHandle failed")
      rout[i].pid = getpid()
      rout[i].off = cast[int](addr ex.rbuf[ne*m.start]) - cast[int](base)
      copyMem(addr rout[i].fd, addr rout[i].h.data[0], sizeof(cint))
      c.pushSend(m.rank, addr rout[i], sizeof(Info))
      inc ns
  for i, m in ex.smsg:
    if ex.speer[i]:
      c.pushRecv(m.rank, addr sin[i], sizeof(Info))
      inc nr
  if nr > 0: c.waitRecvs(nr)
  if ns > 0: c.waitSends(ns)
  toc("handles")
  var sd = newSeq[ptr UncheckedArray[T]](ex.nsend)
  var ss = newSeq[int32](ex.nsend)
  for i, m in ex.smsg:
    var base = cast[ptr UncheckedArray[T]](addr ex.sbuf[ne*m.start])
    if ex.speer[i]:
      let pfd = syscall(SYS_pidfd_open, sin[i].pid, 0)
      let fd = cint syscall(SYS_pidfd_getfd, pfd, sin[i].fd, 0)
      if pfd < 0 or fd < 0: qexError("pidfd_getfd failed")
      var h = sin[i].h
      copyMem(addr h.data[0], unsafeAddr fd, sizeof(cint))
      var p: pointer
      if zeMemOpenIpcHandle(ctx, dev, h, 0, p) != 0:
        qexError("zeMemOpenIpcHandle failed")
      ex.peers.add p
      base = cast[ptr UncheckedArray[T]](cast[int](p) + sin[i].off)
    for k in 0..<m.count:
      sd[m.start+k] = cast[ptr UncheckedArray[T]](addr base[k])
      ss[m.start+k] = int32 m.count
  ex.sdst = sd.toDevice
  ex.sstr = ss.toDevice
  toc("peers")
  let qc = CommQmp(c).comm
  var hs: seq[QMP_msghandle_t]
  for i, m in ex.rmsg:
    let mm = if ex.rpeer[i]: QMP_declare_msgmem(addr ex.flags[1+i], 1)
             else: QMP_declare_msgmem(addr ex.rbuf[ne*m.start], csize_t(ne*m.count*sizeof(T)))
    ex.mems.add mm
    hs.add QMP_comm_declare_receive_from(qc, mm, cint m.rank, 0)
  for i, m in ex.smsg:
    let mm = if ex.speer[i]: QMP_declare_msgmem(addr ex.flags[0], 1)
             else: QMP_declare_msgmem(addr ex.sbuf[ne*m.start], csize_t(ne*m.count*sizeof(T)))
    ex.mems.add mm
    hs.add QMP_comm_declare_send_to(qc, mm, cint m.rank, 0)
  if hs.len > 0: ex.msg = QMP_declare_multiple(addr hs[0], cint hs.len)
  ex

template sendSite*(sl, sd, st: untyped; nslot, n, i, c0: int; v: untyped) =
  ## Stores v as reals c0, c0+1, ... of site i in its send slots, inside a
  ## kernel.
  for s in 0..<nslot:
    let k = int sl[s*n + i]
    if k >= 0:
      let dp = sd[k]
      let sk = int st[k]
      for c in 0..<v.len: dp[(c0+c)*sk] = v[c]

template recvSite*(ro, rs, rb: untyped; p: int; v: untyped) =
  ## Loads receive position p into v, inside a kernel.
  let o = int ro[p]
  let sk = int rs[p]
  for c in 0..<v.len: v[c] = rb[o + c*sk]

var gpuPackDep {.exportc.}: int  # orders pack before the sends in start

proc pack*[T](ex: GpuHaloEx[T], f: ptr UncheckedArray[T]) =
  ## Stores the send slots of f; returns before the kernel completes, start
  ## waits for it.
  let ne = ex.ne
  let v = ex.v
  let ns = ex.nsend
  let si = ex.sidx
  let sd = ex.sdst
  let st = ex.sstr
  for t in `||`(0, ne*ns-1, "target teams distribute parallel for " & gpuForClause & " nowait depend(out:gpuPackDep)"):
    let c = t div ns
    let k = t - c*ns
    let j = int si[k]
    sd[k][c*int st[k]] = f[((j div v)*ne + c)*v + j mod v]

proc start*[T](ex: GpuHaloEx[T]) =
  ## Waits for pack, then starts the receives and sends; MPI reads the device
  ## buffers directly.  Kernels storing slots must have completed.
  {.emit: "#pragma omp taskwait depend(inout:gpuPackDep)".}
  if ex.mems.len > 0: discard QMP_start(ex.msg)

proc wait*[T](ex: GpuHaloEx[T]) =
  if ex.mems.len > 0: discard QMP_wait(ex.msg)

proc free*[T](ex: GpuHaloEx[T]) =
  if ex.mems.len > 0:
    QMP_free_msghandle(ex.msg)
    for m in ex.mems: QMP_free_msgmem(m)
    ex.mems.setLen(0)
  if ex.peers.len > 0:
    let (ctx, _) = zeContext()
    for p in ex.peers: discard zeMemCloseIpcHandle(ctx, p)
  ex.peers.setLen(0)
  for p in [pointer ex.sidx, ex.sslot, ex.sdst, ex.sstr, ex.rofs, ex.rstr, ex.sbuf, ex.rbuf]:
    if p != nil: gpuFree(p)
  ex.sidx = nil
  ex.sslot = nil
  ex.sdst = nil
  ex.sstr = nil
  ex.rofs = nil
  ex.rstr = nil
  ex.sbuf = nil
  ex.rbuf = nil

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
