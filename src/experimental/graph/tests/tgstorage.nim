#RUNCMD env OMP_NUM_THREADS=1 $RUNJOB
## Joint graph evaluation with plan-owned reusable buffers.
import base/globals
setVLENmax(4)

import math, unittest, std/[tables, sets]
import qex except epsilon
import base/alignedMem
import helpers
import ../[core, scalar, gauge, functional, multi, plan]
import ../gauge/[types, field_ops, transport, matrix, stencil]

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))
qexInit()
letParam:
  expectRanks = nRanks
check nRanks == expectRanks
letParam:
  lat = latticeFromLocalLattice(@[4,4], nRanks)
let lo = lat.newLayout
var rng = lo.newRNGField(Philox4x64, 921260913u64)
let g = lo.newGauge
let u = lo.newGauge
let m = lo.newGauge
threads:
  g.random rng
  u.random rng
  m.randomTAH rng
  for f in m:
    f *= 0.1

proc sameRaw(a, b: auto) =
  let d = a.newOneOf
  threads: d := a - b
  check norm2(d) < 1e-19 * (1.0 + norm2(b))

proc same(a, b: Ggauge) =
  let want = b.eval.gval
  for mu in 0..<want.len:
    sameRaw(a.gval[mu], want[mu])

proc same[F](a, b: GfieldOf[F]) =
  sameRaw(a.fval, b.eval.fval)

proc same(a, b: Gscalar) =
  let want = b.eval.sval
  check abs(a.sval-want) < 1e-11 * (1.0 + abs(want))

proc snapshot(a: auto): auto =
  let z = a.newOneOf
  threads: z := a
  z

proc forwards(rt: GraphRuntime): int =
  for s in rt.runStatsByNode.values:
    result += s.count

proc counts(rt: GraphRuntime): Table[NodeId, int] =
  for id, s in rt.runStatsByNode:
    result[id] = s.count

proc once(rt: GraphRuntime, before: Table[NodeId, int]) =
  for id, s in rt.runStatsByNode:
    check s.count-before.getOrDefault(id) in 0..1

proc countOps(roots: openArray[Gvalue]): int =
  var seen = initHashSet[NodeKey]()
  var n = 0
  proc walk(v: Gvalue) =
    if not seen.markSeenNode(v):
      return
    for x in v.inputs:
      walk(x)
    if v.gfunc != nil and v.gfunc.forward != nil:
      inc n
  for v in roots:
    walk(v)
  n

proc stepF[F](v: Gvalue) =
  let z = GfieldOf[F](v)
  let x = GfieldOf[F](v.inputs[0])
  threads: z.fval := 1.125 * x.fval

proc step[F](x: GfieldOf[F], overwrite = false): GfieldOf[F] =
  # The same pointwise kernel tests both a conservative and an explicit proof.
  let fn = Gfunc(forward: stepF[F], name: "storage step",
                 bufferMode: bmFull, inplace: (if overwrite: @[0] else: @[]))
  graphNode(x.fieldNodeLike, @[Gvalue(x)], fn, "storage step")

proc joinF(v: Gvalue) =
  let z = Gfield(v)
  let a = Gfield(v.inputs[0])
  let b = Gfield(v.inputs[1])
  threads: z.fval := a.fval + b.fval

proc join(a, b: Gfield): Gfield =
  graphNode(a.sameShapeFieldNodeLike(b, "storage join"), @[Gvalue(a), Gvalue(b)],
    Gfunc(forward: joinF, name: "storage join", bufferMode: bmFull), "storage join")

proc failF(v: Gvalue) =
  let z = Gfield(v)
  let x = Gfield(v.inputs[0])
  # A failed forward can dirty its assigned buffer before raising.
  threads: z.fval := 7.0 * x.fval
  if Gint(v.inputs[1]).ival != 0:
    raiseError("requested storage forward failure")
  threads: z.fval := 1.125 * x.fval

proc mayFail(x: Gfield, flag: Gint): Gfield =
  graphNode(x.fieldNodeLike, @[Gvalue(x), Gvalue(flag)],
    Gfunc(forward: failF, name: "storage failure", bufferMode: bmFull), "storage failure")

proc failScalarF(v: Gvalue) =
  if Gint(v.inputs[1]).ival != 0:
    raiseError("requested scalar storage failure")
  Gscalar(v).sval = Gscalar(v.inputs[0]).sval

proc mayFail(x: Gscalar, flag: Gint): Gscalar =
  graphNode(x.scalarNodeLike, @[Gvalue(x), Gvalue(flag)],
    Gfunc(forward: failScalarF, name: "storage scalar failure"), "storage scalar failure")

proc rawIdF(v: Gvalue) =
  Gscalar(v).sval = Gscalar(v.inputs[0]).sval

proc rawIdView(v: Gvalue, mode: InputWalkMode, visit: GnodeVisit) =
  if mode == iwmEval:
    visit v.inputs[0]

proc rawId(x: Gscalar): Gscalar =
  graphNode(x.scalarNodeLike, @[Gvalue(x)],
    Gfunc(forward: rawIdF, inputView: rawIdView, name: "storage raw identity",
      bufferMode: bmFull), "storage raw identity")

type Gpublication = ref object of Gscalar
  fail: Gint

method newOneOf(x: Gpublication): Gvalue =
  Gpublication(runtime: x.runtime, fail: x.fail).assignStableNodeId

method valueLike(x: Gpublication): Gvalue = x.newOneOf

method valAlias(z: Gpublication, x: Gvalue) =
  if z.fail.ival != 0:
    raiseValueError("requested storage publication failure")
  z.sval = Gscalar(x).sval

proc invalid(v: Gvalue) =
  check v.stale
  check not v.valueReady
  expect(GraphError):
    discard v.eval
  check v.stale
  check not v.valueReady

type GnestedCopy = ref object of Gfield
  cache: Gfield

method newOneOf(x: GnestedCopy): Gvalue =
  GnestedCopy(runtime: x.runtime, fval: x.fval.newShape).assignStableNodeId

proc nestedCopyF(v: Gvalue) =
  let z = GnestedCopy(v)
  let x = Gfield(v.inputs[0])
  if z.cache == nil:
    z.cache = z.runtime.toGvalue(x.fval)
  else:
    z.cache.update(x.fval)
  discard z.cache.eval
  z.valCopy(z.cache)

proc nestedCopy(x: Gfield): Gfield =
  graphNode(GnestedCopy(runtime: x.runtime, fval: x.fval.newShape),
    @[Gvalue(x)], Gfunc(forward: nestedCopyF, name: "storage nested copy", bufferMode: bmFull),
    "storage nested copy")

type GdedicatedSubset = ref object of Gfield

method newOneOf(x: GdedicatedSubset): Gvalue =
  GdedicatedSubset(runtime: x.runtime, fval: x.fval.newShape).assignStableNodeId

method bufferProto(x: GdedicatedSubset): Gvalue = nil

proc dedicatedSubsetF(v: Gvalue) =
  let z = GdedicatedSubset(v)
  let x = Gfield(v.inputs[0])
  if Gint(v.inputs[1]).ival != 0:
    threads: z.fval := 9.0
    raiseError("requested dedicated subset failure")
  let parity = Gint(v.inputs[2]).ival
  threads:
    if parity == 0:
      for e in z.fval.even:
        z.fval[e] := x.fval[e]
    else:
      for e in z.fval.odd:
        z.fval[e] := x.fval[e]

proc dedicatedSubset(x: Gfield, flag, parity: Gint): Gfield =
  graphNode(GdedicatedSubset(runtime: x.runtime, fval: x.fval.newShape),
    @[Gvalue(x), Gvalue(flag), Gvalue(parity)],
    Gfunc(forward: dedicatedSubsetF, name: "storage dedicated subset", bufferMode: bmZero),
    "storage dedicated subset")

type
  PoolData = ref object
    kind: int
    vals: array[2, float]
  PoolAudit = ref object
    checks, allocs: int
    homes: seq[tuple[want, got: int]]
    held: PoolData
  GpoolValue = ref object of Gvalue
    family, kind, allocKind: int
    proto: bool
    data: PoolData
    audit: PoolAudit
  GpoolNested = ref object of GpoolValue
    cache: GpoolValue

method newOneOf(x: GpoolValue): Gvalue =
  GpoolValue(runtime: x.runtime, family: x.family, kind: x.kind,
    allocKind: x.allocKind, audit: x.audit).assignStableNodeId

method bufferProto(x: GpoolValue): Gvalue =
  GpoolValue(runtime: x.runtime, family: x.family, kind: x.kind,
    allocKind: x.allocKind, proto: true, audit: x.audit).assignStableNodeId

method hasStorage(x: GpoolValue): bool = x.data != nil

method ensureStorage(x: GpoolValue) =
  if x.data == nil:
    inc x.audit.allocs
    x.data = PoolData(kind: (if x.proto: x.allocKind else: x.kind))

method releaseStorage(x: GpoolValue) = x.data = nil
method bufferBytes(x: GpoolValue): int = 2 * sizeof(float)
method clearBuffer(x: GpoolValue) = x.data.vals = [0.0, 0.0]

method bufferCompatible(x: GpoolValue, y: Gvalue): bool =
  inc x.audit.checks
  if not (y of GpoolValue): return false
  let other = GpoolValue(y)
  let a = if x.data == nil: x.kind else: x.data.kind
  let b = if other.data == nil: other.kind else: other.data.kind
  # 0 ~ 1 and 1 ~ 2, but 0 !~ 2.
  x.family == other.family and abs(a-b) <= 1

method bindBuffer(x: GpoolValue, buffer: Gvalue) =
  let b = GpoolValue(buffer)
  doAssert b.data == nil or
    (x.family == b.family and abs(x.kind-b.data.kind) <= 1)
  x.data = b.data

method valCopy(z: GpoolValue, x: Gvalue) =
  z.ensureStorage
  z.data.vals = GpoolValue(x).data.vals

method valAlias(z: GpoolValue, x: Gvalue) = z.data = GpoolValue(x).data

proc update(x: GpoolValue, value: float) =
  x.ensureStorage
  x.data.vals = [value, 2.0*value]
  x.updated

proc poolInput(rt: GraphRuntime, audit: PoolAudit, family, kind: int,
               value: float): GpoolValue =
  result = GpoolValue(runtime: rt, family: family, kind: kind,
    allocKind: kind, audit: audit).assignStableNodeId
  result.update(value)

proc poolStepF(v: Gvalue) =
  let z = GpoolValue(v)
  let x = GpoolValue(v.inputs[0])
  z.audit.homes.add (z.kind, z.data.kind)
  if v.inputs.len > 1 and Gint(v.inputs[1]).ival != 0:
    z.data.vals = [99.0, 99.0]
    raiseError("requested pool forward failure")
  z.data.vals[0] = x.data.vals[0]+1.0
  if v.gfunc.bufferMode != bmZero:
    z.data.vals[1] = x.data.vals[1]+2.0

proc poolStep(x: GpoolValue, kind = -1, alloc = -1, overwrite = false,
              zero = false, flag: Gint = nil): GpoolValue =
  let tag = if kind < 0: x.kind else: kind
  let dst = if alloc < 0: tag else: alloc
  var inputs = @[Gvalue(x)]
  if flag != nil: inputs.add flag
  graphNode(GpoolValue(runtime: x.runtime, family: x.family, kind: tag,
    allocKind: dst, audit: x.audit), inputs,
    Gfunc(forward: poolStepF, name: "storage pool step",
      bufferMode: (if zero: bmZero else: bmFull),
      inplace: (if overwrite: @[0] else: @[])), "storage pool step")

proc poolSumF(v: Gvalue) =
  Gscalar(v).sval = 0.0
  for x in v.inputs:
    Gscalar(v).sval += GpoolValue(x).data.vals[0]

proc poolSum(values: varargs[GpoolValue]): Gscalar =
  graphNode(Gscalar(runtime: values[0].runtime), values,
    Gfunc(forward: poolSumF, name: "storage pool sum", bufferMode: bmFull),
    "storage pool sum")

method newOneOf(x: GpoolNested): Gvalue =
  GpoolNested(runtime: x.runtime, family: x.family, kind: x.kind,
    allocKind: x.allocKind, audit: x.audit).assignStableNodeId

proc poolHoldF(v: Gvalue) =
  let z = GpoolValue(v)
  let x = GpoolValue(v.inputs[0])
  z.audit.held = x.data
  z.valAlias(x)

proc poolNestedF(v: Gvalue) =
  let z = GpoolNested(v)
  let x = GpoolValue(v.inputs[0])
  if z.cache == nil:
    z.cache = graphNode(GpoolValue(runtime: z.runtime, family: x.family,
      kind: x.kind, allocKind: x.allocKind, audit: x.audit), @[Gvalue(x)],
      Gfunc(forward: poolHoldF, name: "storage opaque hold", bufferMode: bmOpaque),
      "storage opaque hold")
  discard z.cache.eval
  z.valCopy(z.cache)

proc poolNested(x: GpoolValue): GpoolValue =
  graphNode(GpoolNested(runtime: x.runtime, family: x.family, kind: x.kind,
    allocKind: x.allocKind, audit: x.audit), @[Gvalue(x)],
    Gfunc(forward: poolNestedF, name: "storage nested hold", bufferMode: bmFull),
    "storage nested hold")

proc poolWork(n: int): tuple[first, warm: int] =
  let rt = initGraphRuntime()
  let audit = PoolAudit()
  let x = poolInput(rt, audit, 0, 0, 1.0)
  let y = poolInput(rt, audit, 1, 0, 2.0)
  var wide: seq[GpoolValue]
  for _ in 0..<n: wide.add poolStep(x)
  var z = y
  for _ in 0..<n: z = poolStep(z)
  let p = plan(poolSum(wide), z)
  audit.checks = 0
  discard p.eval()
  result.first = audit.checks
  check result.first <= 6*(2*n)
  check p.stats.buffers == n+2
  check Gscalar(p[0]).sval == 2.0*float(n)
  check GpoolValue(p[1]).data.vals == [2.0+float(n), 4.0+2.0*float(n)]
  let allocs = audit.allocs
  let bytes = p.stats.arenaBytes
  x.update(3.0)
  y.update(4.0)
  audit.checks = 0
  discard p.eval()
  result.warm = audit.checks
  check result.warm <= 4*(2*n)
  check p.stats.buffers == n+2
  check p.stats.arenaBytes == bytes
  check audit.allocs == allocs
  check Gscalar(p[0]).sval == 4.0*float(n)
  check GpoolValue(p[1]).data.vals == [4.0+float(n), 8.0+2.0*float(n)]
  audit.checks = 0
  discard p.eval()
  check audit.checks == 0
  p.clear()

proc packedF(v: Gvalue) =
  let z = Gmulti(v)
  for i in 0..<z.len:
    z.storedSlot(i).valCopy(v.inputs[i])

proc packed(values: varargs[Gvalue]): Gmulti =
  newMultiOutputNode(values, values,
    Gfunc(forward: packedF, name: "storage packed copy", bufferMode: bmFull),
    "storage packed copy")

proc chain(overwrite: bool) =
  let rt = initGraphRuntime()
  let keep = initGraphRuntime()
  let x = rt.toGvalue(g[0])
  let kx = keep.toGvalue(g[0])
  let input = x.fval.s.data
  var z = x
  var kz = kx
  for _ in 0..<6:
    z = step(z, overwrite)
    kz = 1.125 * kz
  let raw = getRawMemAllocated()
  let p = plan(z)
  let audits = p.stats.sourceAudits
  let rev = rt.boundaryRevision
  check getRawMemAllocated() == raw
  discard p.eval()
  check p.stats.sourceAudits == audits
  check rt.boundaryRevision == rev
  check p.stats.runs == 1
  check p.stats.forwards == 6
  check forwards(rt) == 6
  check p.stats.buffers == (if overwrite: 1 else: 2)
  check p.stats.reuses >= (if overwrite: 5 else: 4)
  check p.stats.peakLiveBytes > 0
  check p.stats.peakLiveBytes <= p.stats.arenaBytes
  check p[0].gfunc == nil
  check p[0].inputs.len == 0
  check x.fval.s.data == input
  check Gfield(p[0]).fval.s.data != input
  same(Gfield(p[0]), kz)
  sameRaw(x.fval, g[0])
  let n = forwards(rt)
  let f = p.stats.forwards
  let mem = getRawMemAllocated()
  discard p.eval()
  check p.stats.runs == 2
  check p.stats.sourceAudits == audits
  check p.stats.forwards == f
  check forwards(rt) == n
  check getRawMemAllocated() == mem
  x.update(u[0])
  kx.update(u[0])
  let before = counts(rt)
  let bytes = p.stats.arenaBytes
  let used = getRawMemAllocated()
  discard p.eval()
  check p.stats.runs == 3
  check p.stats.sourceAudits == audits
  check rt.boundaryRevision == rev
  check p.stats.forwards-f == 6
  check forwards(rt)-n == 6
  once(rt, before)
  check getRawMemAllocated() == used
  check p.stats.arenaBytes == bytes
  check p.stats.buffers == (if overwrite: 1 else: 2)
  same(Gfield(p[0]), kz)
  sameRaw(x.fval, u[0])

proc moves(g: Ggauge, f: Gfield, n: int): Gfield =
  result = f
  let nd = f.fval.l.nDim
  for i in 0..<n:
    let mu = i mod nd
    let sgn = if ((i div nd) and 1) == 0: 1 else: -1
    if ((i div (2*nd)) and 1) == 0:
      result = shift(result, mu, sgn)
    else:
      result = hop(g, result, mu, sgn)

suite "graph shared storage plans":
  setup:
    let rt = initGraphRuntime()
    let keep = initGraphRuntime()
    let x = rt.toGvalue(g)
    let y = rt.toGvalue(u)
    let b = rt.toGvalue(m)
    let kx = keep.toGvalue(g)
    let ky = keep.toGvalue(u)
    let kb = keep.toGvalue(m)
    defer:
      rt.resetGradCache
      rt.resetApplyCache
      rt.resetLdjCache
      keep.resetGradCache
      keep.resetApplyCache
      keep.resetLdjCache

  test "a preserved input chain reuses two arena buffers":
    chain(false)

  test "an explicit in-place chain reuses one arena buffer":
    chain(true)

  test "shift and hop destinations require the prototype layout":
    let lo2 = lo.physGeom.newLayout
    let dst = lo2.newGauge
    expect(ValueError):
      discard newShifter(g[0], 0, 1, dest = dst[0])
    expect(ValueError):
      discard newTransporter(g[0], g[0], 0, -1, dest = dst[0])

  test "shift and hop requests stay inactive and stable across output rebinding":
    let ds = [g[0].newOneOf, g[0].newOneOf]
    let dt = [g[0].newOneOf, g[0].newOneOf]
    for mu in 0..<lo.nDim:
      for sgn in [-1, 1]:
        let ks = newShifter(g[0], mu, sgn)
        var kt = newTransporter(g[mu], g[0], mu, sgn)
        check ks.field.s.data != g[0].s.data
        check kt.field.s.data != g[0].s.data
        check kt.field.s.data != g[mu].s.data
        let raw = getRawMemAllocated()
        var sh = newShifter(g[0], mu, sgn, dest = ds[0])
        var tr = newTransporter(g[mu], g[0], mu, sgn, dest = dt[0])
        check getRawMemAllocated() == raw
        check sh.field.s.data == ds[0].s.data
        check tr.field.s.data == dt[0].s.data
        let ss = sh.sb.sb
        let ts = tr.sb.sb
        let req = [cast[pointer](ss.sq.smsg), cast[pointer](ss.sq.rmsg),
                   cast[pointer](ss.sq.pairmsg), cast[pointer](ts.sq.smsg),
                   cast[pointer](ts.sq.rmsg), cast[pointer](ts.sq.pairmsg)]
        let buf = [ss.sq.sbuf, ss.sq.rbuf, ss.lbuf, ts.sq.sbuf, ts.sq.rbuf, ts.lbuf]
        check sh.sb.si.sq.nSendSites1 == sh.sb.si.sq.nRecvSites1
        check sh.sb.si.nSendRanks == sh.sb.si.nRecvRanks
        for pass in 0..<4:
          let src = if (pass and 1) == 0: g[0] else: u[0]
          let lnk = if pass < 2: g[mu] else: u[mu]
          sh.field = ds[pass and 1]
          tr.field = dt[pass and 1]
          tr.setLink(lnk)
          kt.setLink(lnk)
          threads:
            discard sh ^* src
            discard tr ^* src
            discard ks ^* src
            discard kt ^* src
          sameRaw(sh.field, ks.field)
          sameRaw(tr.field, kt.field)
          check sh.sb.sb == ss
          check tr.sb.sb == ts
          for i, sb in [sh.sb.sb, tr.sb.sb]:
            check not sb.activeRecv
            check not sb.activeSend
            check cast[pointer](sb.sq.smsg) == req[3*i]
            check cast[pointer](sb.sq.rmsg) == req[3*i+1]
            check cast[pointer](sb.sq.pairmsg) == req[3*i+2]
            check sb.sq.sbuf == buf[3*i]
            check sb.sq.rbuf == buf[3*i+1]
            check sb.lbuf == buf[3*i+2]
            if sh.sb.si.nSendRanks > 0:
              check req[3*i+2] != nil
          sh.field = g[0].newShape
          tr.field = g[0].newShape
          tr.clearLink

  test "mixed shift and hop chains use two arena fields at every length":
    let f = rt.toGvalue(g[0])
    let kf = keep.toGvalue(g[0])
    let fp = f.fval.s.data
    let gp = x.gval
    for n in [2, 4*lo.nDim, 12*lo.nDim]:
      let z = moves(x, f, n)
      let kz = moves(kx, kf, n)
      check z.bufferProto != nil
      check z.gfunc.inplace.len == 0
      let raw = getRawMemAllocated()
      let p = plan(z)
      check getRawMemAllocated() == raw
      for pass in 0..2:
        let src = if pass == 1: u[0] else: g[0]
        let lnk = if pass == 2: u else: g
        f.update(src)
        kf.update(src)
        x.update(lnk)
        kx.update(lnk)
        let used = getRawMemAllocated()
        let before = p.stats.forwards
        discard p.eval()
        check p.stats.forwards-before == n
        check p.stats.buffers == 2
        check p.stats.arenaBytes == 2*f.bufferBytes
        if pass > 0:
          check getRawMemAllocated() == used
        check f.fval.s.data == fp
        check Gfield(p[0]).fval.s.data != fp
        for mu in 0..<gp.len:
          check x.gval[mu].s.data == gp[mu].s.data
          sameRaw(x.gval[mu], lnk[mu])
        sameRaw(f.fval, src)
        same(Gfield(p[0]), kz)
      let before = p.stats.forwards
      discard p.eval()
      check p.stats.forwards == before
      p.clear()

  test "move outputs preserve source aliases through failure retry and double clear":
    let flag = rt.toGvalue(1)
    let f = rt.toGvalue(g[0])
    let kf = keep.toGvalue(g[0])
    let sh = shift(f, 0, -1)
    let z = moves(x, sh, 4*lo.nDim)
    let kz = moves(kx, shift(kf, 0, -1), 4*lo.nDim)
    discard z.eval
    let sf = sh.fval
    let zf = z.fval
    let ss = snapshot(sf)
    let zs = snapshot(zf)
    let p = plan(mayFail(z, flag))
    expect(GraphError):
      discard p.eval()
    check not p[0].valueReady
    sameRaw(sf, ss)
    sameRaw(zf, zs)
    flag.update(0)
    f.update(u[0])
    kf.update(u[0])
    x.update(u)
    kx.update(u)
    let raw = getRawMemAllocated()
    discard p.eval()
    check getRawMemAllocated() == raw
    check p.stats.buffers == 2
    same(Gfield(p[0]), 1.125*kz)
    let outp = Gfield(p[0]).fval
    let saved = snapshot(outp)
    for _ in 0..<2:
      p.clear()
      check p.stats.buffers == 0
      check p.stats.arenaBytes == 0
      GC_fullCollect()
      sameRaw(outp, saved)
      sameRaw(sf, ss)
      sameRaw(zf, zs)
    f.update(g[0])
    kf.update(g[0])
    x.update(g)
    kx.update(g)
    discard p.eval()
    check p.stats.buffers == 2
    check Gfield(p[0]).fval.s.data != outp.s.data
    check sh.fval.s.data == sf.s.data
    check z.fval.s.data == zf.s.data
    same(Gfield(p[0]), 1.125*kz)
    sameRaw(outp, saved)
    sameRaw(sf, ss)
    sameRaw(zf, zs)
    p.clear()
    GC_fullCollect()
    sameRaw(outp, saved)

  test "pooled shifts and hops preserve field link and mixed gradients after updates":
    let f = rt.toGvalue(g[0])
    let kf = keep.toGvalue(g[0])
    let z = moves(x, f, 4*lo.nDim)
    let kz = moves(kx, kf, 4*lo.nDim)
    let score = redot(z, linkField(b, 0))
    let ks = redot(kz, linkField(kb, 0))
    let df = Gfield(grad(score, f))
    let kdf = Gfield(grad(ks, kf))
    let dg = Ggauge(grad(score, x))
    let kdg = Ggauge(grad(ks, kx))
    let mixed = Gfield(grad(redot(dg, b), f))
    let km = Gfield(grad(redot(kdg, kb), kf))
    let p = plan(z, score, df, dg, mixed)
    var bytes = 0
    for pass in 0..2:
      if pass == 1:
        f.update(u[0])
        kf.update(u[0])
        x.update(u)
        kx.update(u)
        b.update(g)
        kb.update(g)
      elif pass == 2:
        f.update(m[0])
        kf.update(m[0])
        x.update(g)
        kx.update(g)
        b.update(u)
        kb.update(u)
      let before = counts(rt)
      let raw = getRawMemAllocated()
      discard p.eval()
      once(rt, before)
      if pass == 0:
        bytes = p.stats.arenaBytes
      else:
        check p.stats.arenaBytes == bytes
        check getRawMemAllocated() == raw
      same(Gfield(p[0]), kz)
      same(Gscalar(p[1]), ks)
      same(Gfield(p[2]), kdf)
      same(Ggauge(p[3]), kdg)
      same(Gfield(p[4]), km)
    check p.stats.reuses > 0

  test "nontransitive requesters use homogeneous prototype storage":
    let audit = PoolAudit()
    let f = poolInput(rt, audit, 0, 1, 3.0)
    let a = poolStep(f, 1)
    let mid = poolStep(f, 0)
    let c = poolStep(f, 2)
    let p = plan(a, poolSum(mid), c)
    audit.allocs = 0
    for value in [3.0, 5.0, -2.0]:
      f.update(value)
      discard p.eval()
      check p.stats.buffers == 2
      check audit.allocs == 2
      check GpoolValue(p[0]).data.kind == 1
      check GpoolValue(p[2]).data.kind == 1
      check GpoolValue(p[0]).data.vals == [value+1.0, 2.0*value+2.0]
      check Gscalar(p[1]).sval == value+1.0
      check GpoolValue(p[2]).data.vals == [value+1.0, 2.0*value+2.0]
      check GpoolValue(p[0]).data != GpoolValue(p[2]).data
    check (want: 0, got: 1) in audit.homes
    check (want: 2, got: 1) in audit.homes

  test "in-place transfers return each allocation to its original pool once":
    let audit = PoolAudit()
    let f = poolInput(rt, audit, 0, 0, 2.0)
    let h = poolInput(rt, audit, 0, 2, 12.0)
    let a = poolStep(f)
    let c = poolStep(poolStep(h), 1, overwrite = true)
    let d = poolStep(h)
    let e = poolStep(d)
    let p = plan(a, poolSum(c), d, e)
    audit.allocs = 0
    for value in [2.0, 5.0, -1.0]:
      f.update(value)
      h.update(value+10.0)
      discard p.eval()
      check p.stats.buffers == 3
      check audit.allocs == 3
      check GpoolValue(p[0]).data.kind == 0
      check Gscalar(p[1]).sval == value+12.0
      check GpoolValue(p[2]).data.vals == [value+11.0, 2.0*value+22.0]
      check GpoolValue(p[3]).data.vals == [value+12.0, 2.0*value+24.0]
      check GpoolValue(p[2]).data != GpoolValue(p[3]).data
      check p.stats.peakLiveBytes <= p.stats.arenaBytes
    check (want: 1, got: 2) in audit.homes
    check p.stats.reuses >= 6

  test "rejected concrete slots stay available while dedicated outputs retry and clear":
    let audit = PoolAudit()
    let f = poolInput(rt, audit, 0, 1, 3.0)
    let flag = rt.toGvalue(1)
    let a = poolStep(f, 1, alloc = 2)
    let mid = poolStep(f, 0, zero = true, flag = flag)
    let c = poolStep(f, 2)
    let p = plan(poolSum(a), mid, c)
    audit.allocs = 0
    expect(GraphError):
      discard p.eval()
    check p.stats.buffers == 1
    check audit.allocs == 2
    flag.update(0)
    for value in [3.0, 5.0]:
      f.update(value)
      discard p.eval()
      check p.stats.buffers == 1
      check audit.allocs == 2
      check GpoolValue(p[1]).data.kind == 0
      check GpoolValue(p[1]).data.vals == [value+1.0, 0.0]
      check GpoolValue(p[2]).data.kind == 2
      check GpoolValue(p[2]).data.vals == [value+1.0, 2.0*value+2.0]
    let left = GpoolValue(p[1]).data
    let right = GpoolValue(p[2]).data
    let aold = left.vals
    let bold = right.vals
    for _ in 0..<2:
      p.clear()
      check p.stats.buffers == 0
      check p.stats.arenaBytes == 0
      check left.vals == aold
      check right.vals == bold
    f.update(7.0)
    discard p.eval()
    check p.stats.buffers == 1
    check GpoolValue(p[1]).data != left
    check GpoolValue(p[2]).data != right
    check GpoolValue(p[1]).data.vals == [8.0, 0.0]
    check GpoolValue(p[2]).data.vals == [8.0, 16.0]
    check left.vals == aold
    check right.vals == bold

  test "incompatible new storage disables pooling for the private generation":
    let audit = PoolAudit()
    let f = poolInput(rt, audit, 0, 0, 2.0)
    let flag = rt.toGvalue(1)
    let p = plan(poolStep(f, alloc = 3, zero = true, flag = flag))
    audit.allocs = 0
    expect(GraphError):
      discard p.eval()
    check p.stats.buffers == 0
    check audit.allocs == 2
    flag.update(0)
    audit.checks = 0
    discard p.eval()
    check audit.checks == 0
    check audit.allocs == 2
    check p.stats.buffers == 0
    check GpoolValue(p[0]).data.kind == 0
    check GpoolValue(p[0]).data.vals == [3.0, 0.0]
    f.update(4.0)
    discard p.eval()
    check audit.checks == 0
    check audit.allocs == 2
    check GpoolValue(p[0]).data.vals == [5.0, 0.0]

  test "a nested opaque consumer fixes its live input allocation":
    let audit = PoolAudit()
    let f = poolInput(rt, audit, 0, 0, 2.0)
    let a = poolNested(poolStep(f))
    let z = poolStep(poolStep(f))
    let p = plan(poolSum(a), z)
    discard p.eval()
    let held = audit.held
    let allocs = audit.allocs
    check held.vals == [3.0, 6.0]
    check p.stats.buffers == 3
    for value in [4.0, -1.0, 3.0]:
      f.update(value)
      discard p.eval()
      check audit.held == held
      check held.vals == [value+1.0, 2.0*value+2.0]
      check Gscalar(p[0]).sval == value+1.0
      check GpoolValue(p[1]).data.vals == [value+2.0, 2.0*value+4.0]
      check GpoolValue(p[1]).data != held
      check p.stats.buffers == 3
      check audit.allocs == allocs
      check p.stats.peakLiveBytes <= p.stats.arenaBytes

  test "selection hooks grow linearly with two pools and many irrelevant free slots":
    let small = poolWork(24)
    let large = poolWork(48)
    check large.first <= 2*small.first+4
    check large.warm <= 2*small.warm+4

  test "plan construction rejects structural and nested function results":
    let a = rt.toGvalue(2.0)
    let fn = lambda(a, a*a)
    var calls = 0
    proc fwd(v: Gvalue) =
      inc calls
    let carrier = newMultiStructureNode([Gvalue(a)], [Gvalue(a)],
      Gfunc(forward: fwd, name: "structural output"), "structural output")
    let nested = multiValues("outer", a, multiValues("inner", a, fn))
    for v in [Gvalue(fn), Gvalue(slotVar(fn)), Gvalue(carrier),
              Gvalue(multiValues("nested carrier", a, carrier)), Gvalue(nested)]:
      expect(GraphValueError):
        discard plan(v)
    check calls == 0

  test "nested numeric bundles preserve stock families and equal-size layout identities":
    let lo2 = lo.physGeom.newLayout
    let g2 = lo2.newGauge
    let r = lo.RealMatrix(1)
    let c = lo.ColorMatrix(1)
    let m8 = lo.RealMatrix(8)
    threads:
      for f in g2: f := 1.0
      r := 1.25
      c := 0.75
      m8 := 2.0
    let f = rt.toGvalue(g[0])
    let h = rt.toGvalue(g2[0])
    let gr = rt.toGvalue(r)
    let gc = rt.toGvalue(c)
    let ga = rt.toGvalue(m8)
    let short = rt.toGvalue(@[g[0]])
    let kf = keep.toGvalue(g[0])
    let kh = keep.toGvalue(g2[0])
    let kr = keep.toGvalue(r)
    let kc = keep.toGvalue(c)
    let ka = keep.toGvalue(m8)
    let ks = keep.toGvalue(@[g[0]])
    let scalar = rt.toGvalue(2.0)
    let inner = packed(step(gr), step(gc), step(ga), short*short, scalar*scalar)
    let z = packed(x*y, step(f), step(h), inner)
    let p = plan(z)
    check f.bufferBytes == h.bufferBytes
    check f.fval.l != h.fval.l
    var bytes = 0
    for pass in 0..1:
      let raw = getRawMemAllocated()
      discard p.eval()
      if pass > 0:
        check getRawMemAllocated() == raw
        check p.stats.arenaBytes == bytes
      let outer = Gmulti(p[0])
      let inside = Gmulti(outer.storedSlot(3))
      same(Ggauge(outer.storedSlot(0)), kx*ky)
      same(Gfield(outer.storedSlot(1)), step(kf))
      same(Gfield(outer.storedSlot(2)), step(kh))
      check Gfield(outer.storedSlot(1)).fval.l == lo
      check Gfield(outer.storedSlot(2)).fval.l == lo2
      same(Grfield(inside.storedSlot(0)), step(kr))
      same(Gcfield(inside.storedSlot(1)), step(kc))
      same(Grmat8(inside.storedSlot(2)), step(ka))
      same(Ggauge(inside.storedSlot(3)), ks*ks)
      check Ggauge(inside.storedSlot(3)).gval.len == 1
      check Gscalar(inside.storedSlot(4)).sval == 4.0
      bytes = p.stats.arenaBytes
      if pass == 0:
        threads:
          for q in g2: q *= 0.8
          r *= 1.1
          c *= 0.9
          m8 *= 1.2
        f.update(u[0])
        kf.update(u[0])
        h.update(g2[0])
        kh.update(g2[0])
        gr.update(r)
        kr.update(r)
        gc.update(c)
        kc.update(c)
        ga.update(m8)
        ka.update(m8)
        x.update(u)
        kx.update(u)
        short.update(@[u[0]])
        ks.update(@[u[0]])

  test "raw input overrides survive initial source discovery":
    let s = rt.toGvalue(2.0)
    let a = s*s
    let z = rawId(a)
    a.update(13.0)
    check z.eval.sval == 13.0
    let p = plan(z)
    discard p.eval()
    check Gscalar(p[0]).sval == 13.0
    check a.valueOverride
    check a.runCount == 0
    a.update(17.0)
    discard p.eval()
    check Gscalar(p[0]).sval == 17.0
    check z.eval.sval == 17.0
    let n = p.stats.forwards
    discard p.eval()
    check p.stats.forwards == n

  test "raw generated constants refresh after a warm mutation":
    let zero = Gscalar(rt.localScalar().zeroLike)
    let z = rawId(zero)
    let p = plan(z)
    discard p.eval()
    check Gscalar(p[0]).sval == 0.0
    let audits = p.stats.sourceAudits
    let n = p.stats.forwards
    zero.update(3.0)
    discard p.eval()
    check Gscalar(p[0]).sval == 3.0
    check z.eval.sval == 3.0
    check p.stats.sourceAudits > audits
    check p.stats.forwards == n+1
    let warm = p.stats.sourceAudits
    discard p.eval()
    check p.stats.sourceAudits == warm
    check p.stats.forwards == n+1

  test "unrelated overrides audit once without rerunning current results":
    let a = rt.toGvalue(2.0)
    let other = rt.toGvalue(3.0)
    let changed = other*other
    let p = plan(a*a)
    discard p.eval()
    let sym = rt.symbolicRevision
    let rev = rt.boundaryRevision
    let audits = p.stats.sourceAudits
    let f = p.stats.forwards
    changed.update(11.0)
    check rt.boundaryRevision == rev+1
    check rt.symbolicRevision == sym
    discard p.eval()
    check Gscalar(p[0]).sval == 4.0
    check p.stats.sourceAudits == audits+1
    check p.stats.forwards == f
    check changed.runCount == 0
    changed.update(12.0)
    check rt.boundaryRevision == rev+1
    discard p.eval()
    check p.stats.sourceAudits == audits+1
    check p.stats.forwards == f
    let ep = changed.epoch
    changed.valueReady = false
    discard changed.eval
    check changed.sval == 9.0
    check changed.epoch > ep
    check not changed.valueOverride
    check rt.boundaryRevision == rev+2
    check rt.symbolicRevision == sym
    discard p.eval()
    check p.stats.sourceAudits == audits+2
    check p.stats.forwards == f
    discard p.eval()
    check p.stats.sourceAudits == audits+2
    check p.stats.forwards == f

  test "static zero membership and constant epochs refresh copied sources":
    let a = rt.toGvalue(2.0)
    let zero = rt.toGvalue(0.0)
    let p = plan(a+zero)
    discard p.eval()
    let sym = rt.symbolicRevision
    let rev = rt.boundaryRevision
    let ep = zero.epoch
    let audits = p.stats.sourceAudits
    let f = p.stats.forwards
    zero.markStaticZeroLeaf
    check zero.epoch == ep
    check rt.boundaryRevision == rev+1
    discard p.eval()
    check Gscalar(p[0]).sval == 2.0
    check p.stats.sourceAudits > audits
    check p.stats.forwards == f+1
    let warm = p.stats.sourceAudits
    zero.markStaticZeroLeaf
    check rt.boundaryRevision == rev+1
    discard p.eval()
    check p.stats.sourceAudits == warm
    check p.stats.forwards == f+1
    zero.update(0.0)
    zero.markStaticZeroLeaf
    check zero.epoch > ep
    check rt.boundaryRevision == rev+3
    discard p.eval()
    check Gscalar(p[0]).sval == 2.0
    check p.stats.sourceAudits > warm
    check p.stats.forwards == f+2
    zero.update(3.0)
    check not zero.staticZeroLeaf
    check rt.boundaryRevision == rev+4
    discard p.eval()
    check Gscalar(p[0]).sval == 5.0
    check p.stats.forwards == f+3
    check rt.symbolicRevision == sym

  test "copied generated constants follow original mutations and updates":
    let unit = x.unitGaugeLike
    let ku = kx.unitGaugeLike
    let zero = Ggauge(x.zeroLike)
    let kz = Ggauge(kx.zeroLike)
    let us = norm2(ku)
    let zs = norm2(kz)
    let pu = plan(norm2(unit))
    let pz = plan(norm2(zero))
    discard pu.eval()
    discard pz.eval()
    same(Gscalar(pu[0]), us)
    same(Gscalar(pz[0]), zs)
    let initial = Gscalar(pu[0]).sval
    check initial > 0.0
    check Gscalar(pz[0]).sval == 0.0
    unit.mutateGauge data:
      threads:
        for f in data:
          f *= 2.0
    ku.mutateGauge data:
      threads:
        for f in data:
          f *= 2.0
    zero.update(u)
    kz.update(u)
    discard pu.eval()
    discard pz.eval()
    same(Gscalar(pu[0]), us)
    same(Gscalar(pz[0]), zs)
    check abs(Gscalar(pu[0]).sval - 4.0*initial) < 1e-11
    check Gscalar(pz[0]).sval > 0.0
    let n = forwards(rt)
    let raw = getRawMemAllocated()
    discard pu.eval()
    discard pz.eval()
    check forwards(rt) == n
    check getRawMemAllocated() == raw

  test "diamonds and repeated inputs execute once across joint roots":
    let f = rt.toGvalue(g[0])
    let kf = keep.toGvalue(g[0])
    let a = step(f)
    let left = join(a, a)
    let right = step(a)
    let z = join(left, right)
    let ka = 1.125*kf
    let kz = (ka+ka) + 1.125*ka
    let p = plan(z, a)
    discard p.eval()
    check p.stats.forwards == 4
    check forwards(rt) == 4
    same(Gfield(p[0]), kz)
    same(Gfield(p[1]), ka)
    let n = forwards(rt)
    let mem = getRawMemAllocated()
    discard p.eval()
    check forwards(rt) == n
    check getRawMemAllocated() == mem
    f.update(u[0])
    kf.update(u[0])
    let before = counts(rt)
    let used = getRawMemAllocated()
    discard p.eval()
    check forwards(rt)-n == 4
    check p.stats.forwards == forwards(rt)
    once(rt, before)
    check getRawMemAllocated() == used
    same(Gfield(p[0]), kz)
    same(Gfield(p[1]), ka)

  test "nested evaluation preserves an initialized ordinary field leaf":
    let f = rt.toGvalue(g[0])
    let kf = keep.toGvalue(g[0])
    let z = step(nestedCopy(step(f)))
    let kz = 1.125*(1.125*kf)
    let p = plan(z)
    discard p.eval()
    same(Gfield(p[0]), kz)
    for pass in 0..2:
      if (pass and 1) == 0:
        f.update(u[0])
        kf.update(u[0])
      else:
        f.update(g[0])
        kf.update(g[0])
      let raw = getRawMemAllocated()
      discard p.eval()
      check getRawMemAllocated() == raw
      same(Gfield(p[0]), kz)
    let n = forwards(rt)
    let raw = getRawMemAllocated()
    discard p.eval()
    check forwards(rt) == n
    check getRawMemAllocated() == raw

  test "matrix8 and scalar fields keep separate compatible buffers":
    let a = lo.RealMatrix(8)
    let r = lo.RealMatrix(1)
    threads:
      a := 2.0
      r := 1.2
      for e in a:
        a[e][0,1] := 0.1
        a[e][2,0] := -0.07
    let ga = rt.toGvalue(a)
    let gr = rt.toGvalue(r)
    let ka = keep.toGvalue(a)
    let kr = keep.toGvalue(r)
    let prod = ga*ga
    let f = inverse(prod)
    let z = exp(gr)+sin(gr)+ln(gr)
    let score = sumLogDet(prod)+sum(z)
    let ks = sumLogDet(ka*ka)+sum(exp(kr)+sin(kr)+ln(kr))
    let p = plan(f, z, score, grad(score, ga), grad(score, gr))
    for pass in 0..1:
      let used = getRawMemAllocated()
      let before = counts(rt)
      discard p.eval()
      if pass > 0:
        check getRawMemAllocated() == used
        once(rt, before)
      same(Grmat8(p[0]), inverse(ka*ka))
      same(Grfield(p[1]), exp(kr)+sin(kr)+ln(kr))
      same(Gscalar(p[2]), ks)
      same(Grmat8(p[3]), Grmat8(grad(ks, ka)))
      same(Grfield(p[4]), Grfield(grad(ks, kr)))
      check Grmat8(p[0]).fval.s.bytes != Grfield(p[1]).fval.s.bytes
      let n = forwards(rt)
      let raw = getRawMemAllocated()
      discard p.eval()
      check forwards(rt) == n
      check getRawMemAllocated() == raw
      if pass == 0:
        threads:
          a *= 0.9
          r *= 1.1
        ga.update(a)
        ka.update(a)
        gr.update(r)
        kr.update(r)

  test "joint primal gradient and mixed derivative roots share one traversal":
    let a = x*y + x*x
    let ka = kx*ky + kx*kx
    let score = norm2(a)
    let ks = norm2(ka)
    let d = Ggauge(grad(score, x))
    let kd = Ggauge(grad(ks, kx))
    let mixed = Ggauge(grad(redot(d, b), y))
    let km = Ggauge(grad(redot(kd, kb), ky))
    let p = plan(a, score, d, mixed)
    let ops = countOps([Gvalue(a), Gvalue(score), Gvalue(d), Gvalue(mixed)])
    discard p.eval()
    check p.stats.forwards == ops
    check forwards(rt) == ops
    same(Ggauge(p[0]), ka)
    same(Gscalar(p[1]), ks)
    same(Ggauge(p[2]), kd)
    same(Ggauge(p[3]), km)
    let n = forwards(rt)
    let mem = getRawMemAllocated()
    discard p.eval()
    check forwards(rt) == n
    check getRawMemAllocated() == mem
    x.update(u)
    kx.update(u)
    b.update(g)
    kb.update(g)
    let before = counts(rt)
    let used = getRawMemAllocated()
    discard p.eval()
    check forwards(rt) > n
    check forwards(rt)-n <= ops
    check p.stats.forwards == forwards(rt)
    once(rt, before)
    check getRawMemAllocated() == used
    same(Ggauge(p[0]), ka)
    same(Gscalar(p[1]), ks)
    same(Ggauge(p[2]), kd)
    same(Ggauge(p[3]), km)

  test "plan execution leaves escaped source gauge and field storage intact":
    let a = x*y
    let f = shift(linkField(a, 0), 1, 1)
    discard a.eval
    discard f.eval
    let oldg = a.gval
    let oldf = f.fval
    var savedg = newSeq[typeof(oldg[0])](oldg.len)
    for mu in 0..<oldg.len:
      savedg[mu] = snapshot(oldg[mu])
    let savedf = snapshot(oldf)
    let ga = a.gval[0].s.data
    let fp = f.fval.s.data
    let p = plan(a, f)
    discard p.eval()
    x.update(u)
    kx.update(u)
    discard p.eval()
    same(Ggauge(p[0]), kx*ky)
    same(Gfield(p[1]), shift(linkField(kx*ky, 0), 1, 1))
    check f.fval.s.data == fp
    for mu in 0..<oldg.len:
      check a.gval[mu].s.data == oldg[mu].s.data
      sameRaw(oldg[mu], savedg[mu])
    sameRaw(oldf, savedf)
    check Ggauge(p[0]).gval[0].s.data != ga
    check Gfield(p[1]).fval.s.data != fp

  test "injected and subset outputs clear dirty reused complements":
    let a = x*y
    let inj = injectLink(linkField(a, 0), 1, x)
    let sub = maskSubset(0, 0, inj+x)
    let ki = injectLink(linkField(kx*ky, 0), 1, kx)
    let ks = maskSubset(0, 0, ki+kx)
    let p = plan(inj, sub)
    for pass in 0..2:
      discard p.eval()
      same(Ggauge(p[0]), ki)
      same(Ggauge(p[1]), ks)
      check norm2(Ggauge(p[0]).gval[0]) == 0.0
      check norm2(Ggauge(p[1]).gval[1]) == 0.0
      if pass == 0:
        x.update(u)
        kx.update(u)
      elif pass == 1:
        x.update(g)
        kx.update(g)
    check p.stats.reuses > 0

  test "changing condition selections evaluates only active branches":
    let sel = rt.toGvalue(1)
    let flag = rt.toGvalue(1)
    let f = rt.toGvalue(g[0])
    let kf = keep.toGvalue(g[0])
    let good = step(f)
    let bad = mayFail(f, flag)
    let z = cond(sel, good, bad)
    let p = plan(z)
    discard p.eval()
    same(Gfield(p[0]), 1.125*kf)
    for s in rt.runStatsByNode.values:
      check s.name != "storage failure"
    flag.update(0)
    sel.update(0)
    discard p.eval()
    same(Gfield(p[0]), 1.125*kf)
    f.update(u[0])
    kf.update(u[0])
    sel.update(1)
    flag.update(1)
    discard p.eval()
    same(Gfield(p[0]), 1.125*kf)
    let raw = getRawMemAllocated()
    flag.update(0)
    sel.update(0)
    discard p.eval()
    check getRawMemAllocated() == raw
    same(Gfield(p[0]), 1.125*kf)

  test "zero scalar scaling preserves its lazy failure guard":
    let a = rt.toGvalue(3.0)
    let scale = rt.toGvalue(0.0)
    let flag = rt.toGvalue(1)
    let contribution = mayFail(a, flag)
    let z = Gscalar(contribution.scaleLike(scale))
    let p = plan(z)
    discard p.eval()
    check Gscalar(p[0]).sval == 0.0
    for s in rt.runStatsByNode.values:
      check s.name != "storage scalar failure"
    scale.update(2.0)
    expect(GraphError):
      discard p.eval()
    flag.update(0)
    discard p.eval()
    check Gscalar(p[0]).sval == 6.0
    scale.update(0.0)
    flag.update(1)
    discard p.eval()
    check Gscalar(p[0]).sval == 0.0

  test "nested identity applications keep shared captures live":
    let v = Ggauge(x.newOneOf)
    let id = lambda(v, v)
    let a = x*y
    let copied = Ggauge(apply(id, Ggauge(apply(id, a))))
    let fn = lambda(v, v*v + y*v + y*v)
    let z = Ggauge(apply(fn, copied))
    let ka = kx*ky
    let kz = ka*ka + ky*ka + ky*ka
    let score = norm2(z)+norm2(copied)
    let ks = norm2(kz)+norm2(ka)
    let p = plan(z, copied, grad(score, x))
    for pass in 0..1:
      discard p.eval()
      same(Ggauge(p[0]), kz)
      same(Ggauge(p[1]), ka)
      same(Ggauge(p[2]), Ggauge(grad(ks, kx)))
      if pass == 0:
        y.update(g)
        ky.update(g)
    let n = forwards(rt)
    let raw = getRawMemAllocated()
    discard p.eval()
    check forwards(rt) == n
    check getRawMemAllocated() == raw

  test "joint scalar roots keep lambda binders separate in either root order":
    let v = rt.toGvalue(2.0)
    let arg = rt.toGvalue(3.0)
    let body = v*v
    let fn = lambda(v, body)
    let call = Gscalar(apply(fn, arg))
    let p = plan(body, call)
    let rev = plan(call, body)
    for vals in [(2.0, 3.0), (5.0, -2.0), (0.5, 4.0)]:
      v.update(vals[0])
      arg.update(vals[1])
      discard p.eval()
      check abs(Gscalar(p[0]).sval-vals[0]*vals[0]) < 1e-12
      check abs(Gscalar(p[1]).sval-vals[1]*vals[1]) < 1e-12
      discard rev.eval()
      check abs(Gscalar(rev[0]).sval-vals[1]*vals[1]) < 1e-12
      check abs(Gscalar(rev[1]).sval-vals[0]*vals[0]) < 1e-12
      check v.sval == vals[0]
      check arg.sval == vals[1]
      check body.runCount == 0
      check call.runCount == 0

  test "numeric root overrides do not replace an abstract lambda body":
    let v = rt.toGvalue(2.0)
    let arg = rt.toGvalue(3.0)
    let body = v*v
    let fn = lambda(v, body)
    let call = Gscalar(apply(fn, arg))
    body.update(13.0)
    let p = plan(body, call)
    let rev = plan(call, body)
    for value in [3.0, 4.0, -2.0]:
      arg.update(value)
      discard rev.eval()
      discard p.eval()
      check Gscalar(p[0]).sval == 13.0
      check Gscalar(p[1]).sval == value*value
      check Gscalar(rev[0]).sval == value*value
      check Gscalar(rev[1]).sval == 13.0
      check body.valueOverride
      check body.sval == 13.0
      check v.sval == 2.0
      check arg.sval == value
    v.update(5.0)
    discard rev.eval()
    discard p.eval()
    check Gscalar(p[0]).sval == 25.0
    check Gscalar(p[1]).sval == 4.0
    check Gscalar(rev[0]).sval == 4.0
    check Gscalar(rev[1]).sval == 25.0
    check not body.valueOverride
    check v.sval == 5.0
    check arg.sval == -2.0

  test "raw source overrides preserve nested binder shadowing in either root order":
    let v = rt.toGvalue(2.0)
    let arg = rt.toGvalue(3.0)
    let body = v*v
    let inner = lambda(v, body)
    let outer = lambda(v, inner)
    let call = Gscalar(apply(apply(outer, 11.0), arg))
    let raw = rawId(body)
    body.update(13.0)
    let p = plan(raw, call)
    let rev = plan(call, raw)
    for value in [3.0, 4.0, -2.0]:
      arg.update(value)
      discard p.eval()
      discard rev.eval()
      check Gscalar(p[0]).sval == 13.0
      check Gscalar(p[1]).sval == value*value
      check Gscalar(rev[0]).sval == value*value
      check Gscalar(rev[1]).sval == 13.0
      check body.valueOverride
      check body.sval == 13.0
      check v.sval == 2.0
      check arg.sval == value
      check body.runCount == 0
      check call.runCount == 0

  test "joint field roots preserve binder scope and original leaf storage":
    let v = rt.toGvalue(g[0])
    let arg = rt.toGvalue(u[0])
    let kv = keep.toGvalue(g[0])
    let ka = keep.toGvalue(u[0])
    let body = v*v
    let fn = lambda(v, body)
    let call = Gfield(apply(fn, arg))
    let p = plan(body, call)
    let rev = plan(call, body)
    let vb = v.fval.s.data
    let ab = arg.fval.s.data
    let want = kv*kv
    let kw = ka*ka
    for pass in 0..2:
      if pass == 1:
        v.update(u[0])
        kv.update(u[0])
        arg.update(g[0])
        ka.update(g[0])
      elif pass == 2:
        v.update(m[0])
        kv.update(m[0])
        arg.update(u[0])
        ka.update(u[0])
      discard p.eval()
      same(Gfield(p[0]), want)
      same(Gfield(p[1]), kw)
      discard rev.eval()
      same(Gfield(rev[0]), kw)
      same(Gfield(rev[1]), want)
      check v.fval.s.data == vb
      check arg.fval.s.data == ab
      sameRaw(v.fval, kv.fval)
      sameRaw(arg.fval, ka.fval)
      check body.runCount == 0
      check call.runCount == 0

  test "conditional lambda changes select the current value and derivative":
    let sel = rt.toGvalue(1)
    let ksel = keep.toGvalue(1)
    let v = Ggauge(x.newOneOf)
    let kv = Ggauge(kx.newOneOf)
    let fn = cond(sel, lambda(v, v+y), lambda(v, v*v+y))
    let kfn = cond(ksel, lambda(kv, kv+ky), lambda(kv, kv*kv+ky))
    let z = Ggauge(apply(fn, x))
    let kz = Ggauge(apply(kfn, kx))
    let p = plan(z, grad(norm2(z), x))
    for flag in [1, 0, 1, 0]:
      sel.update(flag)
      ksel.update(flag)
      discard p.eval()
      same(Ggauge(p[0]), kz)
      same(Ggauge(p[1]), Ggauge(grad(norm2(kz), kx)))

  test "computed lambda selectors expose transitive evaluation inputs":
    let aa = rt.toGvalue(2.0)
    let bb = rt.toGvalue(3.0)
    let xx = rt.toGvalue(4.0)
    let ak = keep.toGvalue(2.0)
    let bk = keep.toGvalue(3.0)
    let xk = keep.toGvalue(4.0)
    let sel = aa*bb
    let v = rt.localScalar()
    let vk = keep.localScalar()
    let fn = cond(sel, lambda(v, 2.0*v), lambda(v, 3.0*v))
    let fk = cond(ak*bk, lambda(vk, 2.0*vk), lambda(vk, 3.0*vk))
    let z = Gscalar(apply(fn, xx))
    let zk = Gscalar(apply(fk, xk))
    let p = plan(z)
    let rev = rt.boundaryRevision
    let audits = p.stats.sourceAudits
    for ab in [(2.0, 3.0), (0.0, 3.0), (2.0, 0.0), (-1.0, 1.0)]:
      aa.update(ab[0])
      bb.update(ab[1])
      ak.update(ab[0])
      bk.update(ab[1])
      discard p.eval()
      same(Gscalar(p[0]), zk)
      check rt.boundaryRevision == rev
      check p.stats.sourceAudits == audits
      check sel.runCount == 0
      check fn.runCount == 0
    let n = forwards(rt)
    let raw = getRawMemAllocated()
    discard p.eval()
    check forwards(rt) == n
    check getRawMemAllocated() == raw
    check p.stats.sourceAudits == audits

  test "alternating function bodies reuse warmed shift and hop work":
    let sel = rt.toGvalue(1)
    let f = rt.toGvalue(g[0])
    let kf = keep.toGvalue(g[0])
    let v = Gfield(f.newOneOf)
    let fn = cond(sel,
      lambda(v, shift(v, 0, 1)), lambda(v, hop(y, v, 1, 1)))
    let z = Gfield(apply(fn, f))
    let p = plan(z)
    let shifted = shift(kf, 0, 1)
    let hopped = hop(ky, kf, 1, 1)
    var nodes, misses: int
    for pass in 0..<8:
      let first = (pass and 1) == 0
      sel.update(if first: 1 else: 0)
      if pass mod 3 == 0:
        f.update(u[0])
        kf.update(u[0])
        y.update(g)
        ky.update(g)
      else:
        f.update(g[0])
        kf.update(g[0])
        y.update(u)
        ky.update(u)
      let raw = getRawMemAllocated()
      discard p.eval()
      if pass == 1:
        nodes = rt.runStatsByNode.len
        misses = rt.functional.applyCacheStats.instantiationMisses
        check misses > 0
      elif pass > 1:
        check rt.runStatsByNode.len == nodes
        check rt.functional.applyCacheStats.instantiationMisses == misses
        if nRanks == 1:
          check getRawMemAllocated() == raw
      same(Gfield(p[0]), if first: shifted else: hopped)

  test "manual nonleaf updates remain authoritative until an input changes":
    let a = x*x
    let ka = kx*kx
    let z = a+y
    let kz = ka+ky
    let p = plan(z)
    discard p.eval()
    same(Ggauge(p[0]), kz)
    a.update(u)
    ka.update(u)
    discard p.eval()
    same(Ggauge(p[0]), kz)
    check a.valueOverride
    let rev = rt.boundaryRevision
    let audits = p.stats.sourceAudits
    a.mutateGauge data:
      threads:
        for f in data:
          f *= 0.5
    ka.mutateGauge data:
      threads:
        for f in data:
          f *= 0.5
    discard p.eval()
    same(Ggauge(p[0]), kz)
    check rt.boundaryRevision == rev
    check p.stats.sourceAudits == audits
    x.update(u)
    kx.update(u)
    discard p.eval()
    same(Ggauge(p[0]), kz)
    check not a.valueOverride
    check rt.boundaryRevision == rev+1
    check p.stats.sourceAudits > audits
    let n = forwards(rt)
    let f = p.stats.forwards
    let raw = getRawMemAllocated()
    discard p.eval()
    check forwards(rt) == n
    check p.stats.forwards == f
    check getRawMemAllocated() == raw

  test "published result updates and lost readiness or storage force reevaluation":
    let f = rt.toGvalue(g[0])
    let kf = keep.toGvalue(g[0])
    let kz = 1.125*(1.125*kf)
    let p = plan(step(step(f)))
    discard p.eval()
    let audits = p.stats.sourceAudits
    let ep = f.epoch
    for change in 0..2:
      let published = Gfield(p[0])
      let old = published.epoch
      case change
      of 0:
        published.update(u[0])
        check published.epoch > old
        sameRaw(published.fval, u[0])
      of 1:
        published.valueReady = false
        check published.epoch == old
        check published.hasStorage
      else:
        published.releaseStorage
        check published.epoch == old
        check published.valueReady
        check not published.hasStorage
      let before = p.stats.forwards
      let n = forwards(rt)
      discard p.eval()
      check p.stats.forwards-before == 2
      check p.stats.sourceAudits == audits
      check forwards(rt)-n == 2
      check p[0].valueReady
      check p[0].hasStorage
      check f.epoch == ep
      same(Gfield(p[0]), kz)
      sameRaw(f.fval, g[0])
    let before = p.stats.forwards
    let n = forwards(rt)
    let raw = getRawMemAllocated()
    discard p.eval()
    check p.stats.forwards == before
    check forwards(rt) == n
    check getRawMemAllocated() == raw

  test "ordinary input readiness and explicit storage replacement refresh the plan":
    let f = rt.toGvalue(g[0])
    let kf = keep.toGvalue(g[0])
    let kz = 1.125*kf
    let p = plan(step(f))
    discard p.eval()
    let audits = p.stats.sourceAudits
    let data = f.fval
    let ep = f.epoch
    let before = p.stats.forwards
    f.valueReady = false
    discard p.eval()
    check f.valueReady
    check f.hasStorage
    check f.epoch == ep
    check f.fval.s.data == data.s.data
    check p.stats.forwards-before == 1
    check p.stats.sourceAudits == audits
    same(Gfield(p[0]), kz)
    f.releaseStorage
    check not f.hasStorage
    check f.epoch == ep
    sameRaw(data, g[0])
    # An authoritative leaf has no producer; its owner supplies replacement data.
    f.update(u[0])
    kf.update(u[0])
    check f.hasStorage
    check f.valueReady
    check f.epoch > ep
    check f.fval.s.data != data.s.data
    let next = p.stats.forwards
    discard p.eval()
    check p.stats.forwards-next == 1
    check p.stats.sourceAudits == audits
    same(Gfield(p[0]), kz)
    sameRaw(data, g[0])
    let n = forwards(rt)
    discard p.eval()
    check forwards(rt) == n

  test "lost computed override storage expires its boundary and advances freshness":
    let f = rt.toGvalue(g[0])
    let kf = keep.toGvalue(g[0])
    let a = step(f)
    a.update(u[0])
    let p = plan(step(a))
    let kz = 1.125*(1.125*kf)
    discard p.eval()
    same(Gfield(p[0]), 1.125*keep.toGvalue(u[0]))
    check a.valueOverride
    check a.runCount == 0
    let ep = a.epoch
    let rev = rt.boundaryRevision
    let audits = p.stats.sourceAudits
    a.releaseStorage
    check a.epoch == ep
    check a.valueReady
    check a.valueOverride
    check not a.hasStorage
    let before = p.stats.forwards
    let n = forwards(rt)
    discard p.eval()
    check a.epoch > ep
    check a.hasStorage
    check not a.valueOverride
    check a.runCount == 1
    check rt.boundaryRevision == rev+1
    check p.stats.sourceAudits > audits
    check p.stats.forwards-before == 2
    check forwards(rt)-n == 3
    same(Gfield(p[0]), kz)
    let warm = forwards(rt)
    let checked = p.stats.sourceAudits
    discard p.eval()
    check forwards(rt) == warm
    check p.stats.sourceAudits == checked
    f.update(u[0])
    kf.update(u[0])
    let next = p.stats.forwards
    discard p.eval()
    check p.stats.forwards-next == 2
    check a.runCount == 1
    same(Gfield(p[0]), kz)

  test "cold override expiry rebuilds on the next call and keeps inactive guards lazy":
    let f = rt.toGvalue(g[0])
    let kf = keep.toGvalue(m[0])
    let sel = rt.toGvalue(1)
    let flag = rt.toGvalue(1)
    let a = step(f)
    let hidden = mayFail(f, flag)
    a.update(u[0])
    hidden.update(u[0])
    let p = plan(cond(sel, step(a), hidden))
    let kz = 1.125*(1.125*kf)
    let rev = rt.boundaryRevision
    let audits = p.stats.sourceAudits
    f.update(m[0])
    discard p.eval()
    same(Gfield(p[0]), kz)
    check not a.valueOverride
    check a.runCount == 1
    check hidden.valueOverride
    check hidden.runCount == 0
    check p.stats.forwards == 2
    check rt.boundaryRevision == rev+1
    check p.stats.sourceAudits == audits
    let before = p.stats.forwards
    let n = forwards(rt)
    # Cold validation reaches the active override during execution; audit it next time.
    discard p.eval()
    check p.stats.forwards-before == 3
    check forwards(rt)-n == 3
    check p.stats.sourceAudits > audits
    same(Gfield(p[0]), kz)
    check a.runCount == 1
    check hidden.valueOverride
    check hidden.runCount == 0
    let warm = forwards(rt)
    let checked = p.stats.sourceAudits
    discard p.eval()
    check forwards(rt) == warm
    check p.stats.sourceAudits == checked

  test "override expiry after a failed pass rebuilds once after recovery":
    let f = rt.toGvalue(g[0])
    let kf = keep.toGvalue(m[0])
    let sel = rt.toGvalue(1)
    let guard = rt.toGvalue(1)
    let flag = rt.toGvalue(0)
    let a = step(f)
    let hidden = mayFail(f, guard)
    a.update(u[0])
    hidden.update(u[0])
    let p = plan(cond(sel, step(a), hidden), mayFail(f, flag))
    let kz = 1.125*(1.125*kf)
    let ks = 1.125*kf
    discard p.eval()
    let rev = rt.boundaryRevision
    let audits = p.stats.sourceAudits
    flag.update(1)
    expect(GraphError):
      discard p.eval()
    check not p[0].valueReady
    check not p[1].valueReady
    check a.valueOverride
    check a.runCount == 0
    check rt.boundaryRevision == rev
    check p.stats.sourceAudits == audits
    check rt.evalFrame == nil
    check rt.applyFrame == nil
    check rt.workFrame == nil
    f.update(m[0])
    flag.update(0)
    let before = p.stats.forwards
    let n = forwards(rt)
    discard p.eval()
    check p.stats.forwards-before == 3
    check forwards(rt)-n == 4
    check not a.valueOverride
    check a.runCount == 1
    check rt.boundaryRevision == rev+1
    check p.stats.sourceAudits == audits
    check hidden.valueOverride
    check hidden.runCount == 0
    same(Gfield(p[0]), kz)
    same(Gfield(p[1]), ks)
    let next = p.stats.forwards
    let runs = forwards(rt)
    discard p.eval()
    check p.stats.forwards-next == 4
    check forwards(rt)-runs == 4
    check p.stats.sourceAudits > audits
    check a.runCount == 1
    check hidden.valueOverride
    check hidden.runCount == 0
    same(Gfield(p[0]), kz)
    same(Gfield(p[1]), ks)
    let warm = forwards(rt)
    let checked = p.stats.sourceAudits
    discard p.eval()
    check forwards(rt) == warm
    check p.stats.sourceAudits == checked

  test "late-bound produced lambdas leave original captures and selectors untouched":
    let f = rt.toGvalue(g[0])
    let c = rt.toGvalue(u[0])
    let kf = keep.toGvalue(g[0])
    let kc = keep.toGvalue(u[0])
    let capture = f*c
    let cap = kf*kc
    discard capture.eval
    let original = capture.fval
    let saved = snapshot(original)
    let capRuns = capture.runCount
    let sel = rt.toGvalue(1)
    let ksel = keep.toGvalue(1)
    let v = Gfield(f.newOneOf)
    let fn = lambdaParam(v, Gfield(v.newOneOf))
    let produced = cond(sel,
      lambda(v, v*v+capture), lambda(v, v+capture))
    fn.valCopy(produced)
    let oldRuns = produced.runCount
    let z = Gfield(apply(fn, f))
    let p = plan(z)
    let want = cond(ksel, kf*kf+cap, kf+cap)
    for flag in [1, 0]:
      sel.update(flag)
      ksel.update(flag)
      discard p.eval()
      same(Gfield(p[0]), want)
      check produced.runCount == oldRuns
      check capture.runCount == capRuns
      check capture.fval.s.data == original.s.data
      sameRaw(original, saved)
    let rebound = cond(sel,
      lambda(v, capture*v), lambda(v, capture+v*v))
    fn.valCopy(rebound)
    let newRuns = rebound.runCount
    let next = cond(ksel, cap*kf, cap+kf*kf)
    f.update(u[0])
    kf.update(u[0])
    c.update(g[0])
    kc.update(g[0])
    for flag in [1, 0]:
      sel.update(flag)
      ksel.update(flag)
      discard p.eval()
      same(Gfield(p[0]), next)
      check rebound.runCount == newRuns
      check produced.runCount == oldRuns
      check capture.runCount == capRuns
      check capture.fval.s.data == original.s.data
      sameRaw(original, saved)

  test "cloned sources and separate plans own independent output buffers":
    let z = x*y+x
    let kz = kx*ky+kx
    let cloned = cloneValues([Gvalue(z)])
    let p = plan(z)
    let q = plan(cloned[0])
    discard p.eval()
    discard q.eval()
    let output = Ggauge(p[0]).gval
    let saved = snapshot(output[0])
    check output[0].s.data != Ggauge(q[0]).gval[0].s.data
    x.update(u)
    kx.update(u)
    discard q.eval()
    same(Ggauge(q[0]), kz)
    sameRaw(output[0], saved)
    discard p.eval()
    same(Ggauge(p[0]), kz)
    check Ggauge(p[0]).gval[0].s.data != Ggauge(q[0]).gval[0].s.data

  test "published leaves and packed selections refuse reads before the first execution":
    let a = rt.toGvalue(2.0)
    let f = rt.toGvalue(g[0])
    let kf = keep.toGvalue(g[0])
    let p = plan(2.0*a, x*x, packed(step(f), a*a))
    let s = Gscalar(p[0])
    let outg = Ggauge(p[1])
    let pair = Gmulti(p[2])
    let outf = Gfield(pair[0])
    let sq = Gscalar(pair[1])
    for v in [p[0], p[1], p[2]]:
      check v.gfunc == nil
      check v.inputs.len == 0
      invalid(v)
    for v in [Gvalue(s*2.0), Gvalue(norm2(outg)), Gvalue(outf), Gvalue(sq)]:
      expect(GraphError):
        discard v.eval
    discard p.eval()
    check s.eval.sval == 4.0
    same(outg.eval, kx*kx)
    same(outf.eval, 1.125*kf)
    check sq.eval.sval == 4.0
    let n = p.stats.forwards
    discard p.eval()
    check p.stats.forwards == n
    for v in [p[0], p[1], p[2]]:
      check not v.stale
      check v.valueReady

  test "failed scalar publications refuse existing cloned captured and planned consumers":
    let a = rt.toGvalue(2.0)
    let flag = rt.toGvalue(0)
    let p = plan(3.0*mayFail(a, flag))
    discard p.eval()
    let pub = Gscalar(p[0])
    let old = 2.0*pub
    let copy = Gscalar(cloneValues([Gvalue(old)])[0])
    let arg = rt.localScalar
    let call = Gscalar(apply(lambda(arg, arg*old), 3.0))
    let dep = plan(old)
    check old.eval.sval == 12.0
    check copy.eval.sval == 12.0
    check call.eval.sval == 36.0
    check Gscalar(dep.eval()[0]).sval == 12.0
    flag.update(1)
    expect(GraphError):
      discard p.eval()
    expect(GraphError):
      discard dep.eval()
    invalid(pub)
    for v in [Gvalue(old), Gvalue(copy), Gvalue(call), Gvalue(3.0*pub)]:
      expect(GraphError):
        discard v.eval
    expect(GraphError):
      discard dep.eval()
    invalid(pub)
    flag.update(0)
    a.update(4.0)
    discard p.eval()
    check p[0].nodeKey == pub.nodeKey
    check pub.eval.sval == 12.0
    check old.eval.sval == 24.0
    check copy.eval.sval == 24.0
    check call.eval.sval == 72.0
    check Gscalar(dep.eval()[0]).sval == 24.0
    p.clear()
    check old.eval.sval == 24.0
    check copy.eval.sval == 24.0
    check call.eval.sval == 72.0
    check Gscalar(dep.eval()[0]).sval == 24.0

  test "a partial publication failure invalidates every result until retry":
    let a = rt.toGvalue(2.0)
    let flag = rt.toGvalue(0)
    let bad = graphNode(Gpublication(runtime: rt, fail: flag), [Gvalue(a)],
      Gfunc(forward: rawIdF, name: "storage publication", bufferMode: bmFull))
    let p = plan(2.0*a, bad)
    discard p.eval()
    let first = Gscalar(p[0])
    let old = first+1.0
    check old.eval.sval == 5.0
    check Gscalar(p[1]).eval.sval == 2.0
    let ep = first.epoch
    flag.update(1)
    a.update(3.0)
    expect(GraphError):
      discard p.eval()
    # The first result was published before the second alias raised.
    check first.epoch > ep
    check first.sval == 6.0
    for v in [p[0], p[1]]:
      invalid(v)
    expect(GraphError):
      discard old.eval
    flag.update(0)
    discard p.eval()
    check p[0].nodeKey == first.nodeKey
    check first.eval.sval == 6.0
    check old.eval.sval == 7.0
    check Gscalar(p[1]).eval.sval == 3.0

  test "rebuilds retire saved results and consumers even after successful recovery":
    let a = rt.toGvalue(2.0)
    let flag = rt.toGvalue(0)
    let zero = Gscalar(a.zeroLike)
    let p = plan(mayFail(a, flag)+zero)
    discard p.eval()
    let old = Gscalar(p[0])
    let use = old*2.0
    let dep = plan(use)
    check use.eval.sval == 4.0
    check Gscalar(dep.eval()[0]).sval == 4.0
    zero.update(1.0)
    flag.update(1)
    expect(GraphError):
      discard p.eval()
    check p[0].nodeKey != old.nodeKey
    expect(GraphError):
      discard dep.eval()
    invalid(old)
    invalid(p[0])
    expect(GraphError):
      discard use.eval
    expect(GraphError):
      discard dep.eval()
    for _ in 0..<2:
      p.clear()
      invalid(old)
      invalid(p[0])
    flag.update(0)
    discard p.eval()
    let fresh = Gscalar(p[0])
    check fresh.eval.sval == 3.0
    invalid(old)
    expect(GraphError):
      discard use.eval
    expect(GraphError):
      discard dep.eval()
    let next = 2.0*fresh
    let newer = plan(next)
    check next.eval.sval == 6.0
    check Gscalar(newer.eval()[0]).sval == 6.0
    p.clear()
    check fresh.eval.sval == 3.0
    check next.eval.sval == 6.0
    check Gscalar(newer.eval()[0]).sval == 6.0
    discard p.eval()
    invalid(fresh)
    expect(GraphError):
      discard next.eval
    expect(GraphError):
      discard newer.eval()
    check Gscalar(p[0]).eval.sval == 3.0

  test "rejected nested plan evaluation preserves valid publications":
    let a = rt.toGvalue(2.0)
    let p = plan(2.0*a)
    discard p.eval()
    let pub = Gscalar(p[0])
    let runs = p.stats.runs
    var hit = false
    proc nested(v: Gvalue) =
      expect(GraphValueError):
        discard p.eval()
      hit = true
      Gscalar(v).sval = Gscalar(v.inputs[0]).sval+1.0
    let z = graphNode(a.scalarNodeLike, [Gvalue(a)],
      Gfunc(forward: nested, name: "storage nested misuse", bufferMode: bmFull))
    let outer = plan(z)
    discard outer.eval()
    check hit
    check p.stats.runs == runs
    check pub.valueReady
    check not pub.stale
    check pub.eval.sval == 4.0
    check Gscalar(outer[0]).sval == 3.0

  test "packed selections refuse failed publications and recover together":
    let a = rt.toGvalue(2.0)
    let flag = rt.toGvalue(0)
    let f = rt.toGvalue(g[0])
    let kf = keep.toGvalue(g[0])
    let p = plan(packed(x*x, mayFail(f, flag), a*a))
    discard p.eval()
    let pub = Gmulti(p[0])
    let outg = Ggauge(pub[0])
    let outf = Gfield(pub[1])
    let outs = Gscalar(pub[2])
    same(outg.eval, kx*kx)
    same(outf.eval, 1.125*kf)
    check outs.eval.sval == 4.0
    flag.update(1)
    expect(GraphError):
      discard p.eval()
    invalid(pub)
    for v in [Gvalue(outg), Gvalue(outf), Gvalue(outs), pub[1]]:
      expect(GraphError):
        discard v.eval
    flag.update(0)
    a.update(3.0)
    f.update(u[0])
    kf.update(u[0])
    discard p.eval()
    same(outg.eval, kx*kx)
    same(outf.eval, 1.125*kf)
    check outs.eval.sval == 9.0

  test "a failed forward can retry after dirtying a reusable buffer":
    let flag = rt.toGvalue(1)
    let f = rt.toGvalue(g[0])
    let kf = keep.toGvalue(g[0])
    let z = step(mayFail(step(f), flag))
    let p = plan(z)
    expect(GraphError):
      discard p.eval()
    invalid(p[0])
    flag.update(0)
    discard p.eval()
    same(Gfield(p[0]), 1.125*(1.125*(1.125*kf)))
    f.update(u[0])
    kf.update(u[0])
    let raw = getRawMemAllocated()
    discard p.eval()
    check getRawMemAllocated() == raw
    same(Gfield(p[0]), 1.125*(1.125*(1.125*kf)))

  test "dedicated partial-write buffers clear complements after failure and updates":
    let f = rt.toGvalue(g[0])
    let kf = keep.toGvalue(g[0])
    let flag = rt.toGvalue(1)
    let parity = rt.toGvalue(0)
    let p = plan(dedicatedSubset(f, flag, parity))
    expect(GraphError):
      discard p.eval()
    flag.update(0)
    discard p.eval()
    check p.stats.buffers == 0
    same(Gfield(p[0]), maskSubset(0, kf))
    check norm2(Gfield(p[0]).fval.odd) == 0.0
    parity.update(1)
    f.update(u[0])
    kf.update(u[0])
    discard p.eval()
    check p.stats.buffers == 0
    same(Gfield(p[0]), maskSubset(1, kf))
    check norm2(Gfield(p[0]).fval.even) == 0.0

  test "failed override validation invalidates previously published outputs":
    let flag = rt.toGvalue(0)
    let f = rt.toGvalue(g[0])
    let kf = keep.toGvalue(g[0])
    let overridden = mayFail(f, flag)
    overridden.update(u[0])
    let p = plan(step(overridden))
    discard p.eval()
    let published = Gfield(p[0])
    same(published, 1.125*keep.toGvalue(u[0]))
    let use = norm2(published)
    discard use.eval
    check published.valueReady
    let n = p.stats.forwards
    flag.update(1)
    expect(GraphError):
      discard p.eval()
    check p.stats.forwards == n
    invalid(published)
    invalid(p[0])
    expect(GraphError):
      discard use.eval
    expect(GraphError):
      discard norm2(published).eval
    check rt.evalFrame == nil
    check rt.applyFrame == nil
    check rt.workFrame == nil
    flag.update(0)
    discard p.eval()
    check p[0].valueReady
    check not p[0].stale
    check p[0].nodeKey == published.nodeKey
    same(Gfield(p[0]), 1.125*(1.125*kf))
    same(use.eval, norm2(1.125*(1.125*kf)))

  test "cache clearing preserves published outputs and plan reuse":
    let v = Ggauge(x.newOneOf)
    let z = Ggauge(apply(lambda(v, v*v+y), x))
    let kz = kx*kx+ky
    let score = norm2(z)
    let ks = norm2(kz)
    let p = plan(z, grad(score, x))
    discard p.eval()
    let output = Ggauge(p[0])
    let d = Ggauge(p[1])
    rt.resetGradCache
    rt.resetApplyCache
    rt.resetLdjCache
    GC_fullCollect()
    same(output, kz)
    same(d, Ggauge(grad(ks, kx)))
    x.update(u)
    kx.update(u)
    discard p.eval()
    same(Ggauge(p[0]), kz)
    same(Ggauge(p[1]), Ggauge(grad(ks, kx)))

  test "double plan clear preserves raw outputs across a new private generation":
    let v = Ggauge(x.newOneOf)
    let st = stapleSum(x, [y])
    let z = Ggauge(apply(lambda(v, v*v+st), x))
    let f = shift(linkField(z, 0), 1, 1)
    let kz = kx*kx+stapleSum(kx, [ky])
    let kf = shift(linkField(kz, 0), 1, 1)
    let p = plan(z, f)
    discard p.eval()
    same(Ggauge(p[0]), kz)
    same(Gfield(p[1]), kf)
    check p.stats.buffers > 0
    check p.stats.arenaBytes > 0
    check p.stats.workspaces > 0
    check p.stats.workspaceBytes > 0
    let pub = Ggauge(p[0])
    let output = pub.gval
    let field = Gfield(p[1]).fval
    var saved = newSeq[typeof(output[0])](output.len)
    for mu in 0..<output.len:
      saved[mu] = snapshot(output[mu])
    let sf = snapshot(field)
    let before = p.stats.forwards
    let audits = p.stats.sourceAudits
    for _ in 0..<2:
      p.clear()
      check p.stats.sourceAudits == audits
      check p.stats.buffers == 0
      check p.stats.arenaBytes == 0
      check p.stats.peakLiveBytes == 0
      check p.stats.workspaces == 0
      check p.stats.workspaceBytes == 0
      check pub.valueReady
      check not pub.stale
      discard pub.eval
      GC_fullCollect()
      for mu in 0..<output.len:
        sameRaw(output[mu], saved[mu])
      sameRaw(field, sf)
    discard p.eval()
    invalid(pub)
    check p.stats.forwards > before
    check p.stats.sourceAudits > audits
    check p.stats.buffers > 0
    check p.stats.arenaBytes > 0
    check p.stats.workspaces > 0
    check p.stats.workspaceBytes > 0
    same(Ggauge(p[0]), kz)
    same(Gfield(p[1]), kf)
    for mu in 0..<output.len:
      check Ggauge(p[0]).gval[mu].s.data != output[mu].s.data
    check Gfield(p[1]).fval.s.data != field.s.data
    GC_fullCollect()
    for mu in 0..<output.len:
      sameRaw(output[mu], saved[mu])
    sameRaw(field, sf)
    x.update(u)
    kx.update(u)
    y.update(g)
    ky.update(g)
    discard p.eval()
    same(Ggauge(p[0]), kz)
    same(Gfield(p[1]), kf)
    for mu in 0..<output.len:
      sameRaw(output[mu], saved[mu])
    sameRaw(field, sf)
    check z.runCount == 0
    check f.runCount == 0
    let n = forwards(rt)
    discard p.eval()
    check forwards(rt) == n

  test "a failed plan can clear twice before a successful new generation":
    let flag = rt.toGvalue(1)
    let f = rt.toGvalue(g[0])
    let kf = keep.toGvalue(g[0])
    let p = plan(step(mayFail(step(f), flag)))
    expect(GraphError):
      discard p.eval()
    invalid(p[0])
    for _ in 0..<2:
      p.clear()
      check p.stats.buffers == 0
      check p.stats.arenaBytes == 0
      check p.stats.workspaces == 0
      check p.stats.workspaceBytes == 0
      invalid(p[0])
    flag.update(0)
    let before = p.stats.forwards
    discard p.eval()
    check p.stats.forwards-before == 3
    check p[0].valueReady
    check not p[0].stale
    check p[0].hasStorage
    same(Gfield(p[0]), 1.125*(1.125*(1.125*kf)))
    let n = forwards(rt)
    discard p.eval()
    check forwards(rt) == n

  test "grouped stout caches retain mixed alpha derivatives across updates":
    let a = rt.toGvalue(0.04)
    let ka = keep.toGvalue(0.04)
    let both = stoutUpdateLogDetJ(x, y, a, 0, 0)
    let sep = stoutUpdate(x, y, a, 0, 0)
    let score = redot(both.Wnew, b) + 0.3*both.lj
    let ksep = stoutUpdate(kx, ky, ka, 0, 0)
    let klj = stoutLogDetJ(kx, ky, ka, 0, 0)
    let ks = redot(ksep, kb) + 0.3*klj
    let da = Gscalar(grad(score, a))
    let kda = Gscalar(grad(ks, ka))
    let dw = Ggauge(grad(score, x))
    let kdw = Ggauge(grad(ks, kx))
    let mixed = Ggauge(grad(da, x))
    let kmixed = Ggauge(grad(kda, kx))
    let p = plan(both.Wnew, both.lj, sep, dw, da, mixed)
    var bytes = 0
    for pass, rho in [0.04, 0.065, 0.04]:
      a.update(rho)
      ka.update(rho)
      if pass == 1:
        x.update(u)
        kx.update(u)
        y.update(g)
        ky.update(g)
      let before = counts(rt)
      discard p.eval()
      once(rt, before)
      same(Ggauge(p[0]), ksep)
      same(Gscalar(p[1]), klj)
      same(Ggauge(p[2]), ksep)
      same(Ggauge(p[3]), kdw)
      same(Gscalar(p[4]), kda)
      same(Ggauge(p[5]), kmixed)
      if pass == 0:
        bytes = p.stats.arenaBytes
      else:
        check p.stats.arenaBytes == bytes
    check p.stats.reuses > 0

  test "stout conditional and application results carry plain gauge values":
    let a = rt.toGvalue(0.04)
    let ka = keep.toGvalue(0.04)
    let sel = rt.toGvalue(1)
    let ksel = keep.toGvalue(1)
    let v = Ggauge(x.newOneOf)
    let kv = Ggauge(kx.newOneOf)
    let cached = stoutUpdate(x, y, a, 1, 1)
    let kc = stoutUpdate(kx, ky, ka, 1, 1)
    let call = Ggauge(apply(lambda(v, stoutUpdate(v, y, a, 1, 1)), x))
    let kcall = Ggauge(apply(lambda(kv, stoutUpdate(kv, ky, ka, 1, 1)), kx))
    let choice = cond(sel, cached, x)
    let kchoice = cond(ksel, kc, kx)
    let composed = stoutUpdate(cached, y, a, 0, 0)
    let kc2 = stoutUpdate(kc, ky, ka, 0, 0)
    let score = redot(call + choice + composed, b)
    let ks = redot(kcall + kchoice + kc2, kb)
    let da = Gscalar(grad(score, a))
    let kda = Gscalar(grad(ks, ka))
    let dw = Ggauge(grad(score, x))
    let kdw = Ggauge(grad(ks, kx))
    check cached.bufferBytes > x.bufferBytes
    check call.bufferBytes == x.bufferBytes
    check choice.bufferBytes == x.bufferBytes
    let p = plan(call, choice, composed, da, dw)
    var bytes = 0
    for pass, flag in [1, 0, 1, 0]:
      sel.update(flag)
      ksel.update(flag)
      a.update(0.04 + 0.005*float(pass))
      ka.update(0.04 + 0.005*float(pass))
      if pass == 1:
        y.update(g)
        ky.update(g)
      discard p.eval()
      same(Ggauge(p[0]), kcall)
      same(Ggauge(p[1]), kchoice)
      same(Ggauge(p[2]), kc2)
      same(Gscalar(p[3]), kda)
      same(Ggauge(p[4]), kdw)
      check p[2].bufferBytes == x.bufferBytes
      if pass == 1:
        bytes = p.stats.arenaBytes
      elif pass > 1:
        check p.stats.arenaBytes == bytes

  test "warm stencil work rebinds gathered and staple input buffers":
    let f = linkField(x, 0)
    let kf = linkField(kx, 0)
    let moved = gather(f, @[1, 0])
    let km = gather(kf, @[1, 0])
    let back = scatter(moved, @[1, 0])
    let kback = scatter(km, @[1, 0])
    let prod = gp(moved, linkField(x, 1), @[0, -1], true)
    let kp = gp(km, linkField(kx, 1), @[0, -1], true)
    let staple = stapleSum(x, [y])
    let kst = stapleSum(kx, [ky])
    let action = plaqSum(x)
    let kact = plaqSum(kx)
    let p = plan(moved, back, prod, staple, action)
    var bytes = 0
    for pass in 0..2:
      if pass == 1:
        x.update(u)
        kx.update(u)
        y.update(m)
        ky.update(m)
      elif pass == 2:
        x.update(g)
        kx.update(g)
        y.update(u)
        ky.update(u)
      let before = counts(rt)
      discard p.eval()
      once(rt, before)
      same(Gfield(p[0]), km)
      same(Gfield(p[1]), kback)
      same(Gfield(p[2]), kp)
      same(Ggauge(p[3]), kst)
      same(Gscalar(p[4]), kact)
      if pass == 0:
        bytes = p.stats.arenaBytes
      else:
        check p.stats.arenaBytes == bytes

  test "distinct staple nodes share one compatible plan workspace":
    let a = stapleSum(x, [y])
    let c = stapleSum(y, [b])
    let d = stapleSum(x+y, [x+b])
    let ka = stapleSum(kx, [ky])
    let kc = stapleSum(ky, [kb])
    let kd = stapleSum(kx+ky, [kx+kb])
    let p = plan(a, c, d)
    var bytes: int
    for pass in 0..2:
      if pass == 1:
        x.update(u)
        kx.update(u)
        y.update(m)
        ky.update(m)
      elif pass == 2:
        x.update(g)
        kx.update(g)
        b.update(u)
        kb.update(u)
      let raw = getRawMemAllocated()
      discard p.eval()
      check p.stats.workspaces == 1
      check p.stats.workspaceBytes > 0
      if pass == 0:
        bytes = p.stats.workspaceBytes
      else:
        check p.stats.workspaceBytes == bytes
        if nRanks == 1:
          check getRawMemAllocated() == raw
      same(Ggauge(p[0]), ka)
      same(Ggauge(p[1]), kc)
      same(Ggauge(p[2]), kd)
    let n = forwards(rt)
    let raw = getRawMemAllocated()
    discard p.eval()
    check forwards(rt) == n
    check getRawMemAllocated() == raw
    check p.stats.workspaces == 1
    check p.stats.workspaceBytes == bytes

qexFinalize()
