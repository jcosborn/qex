## Checks the stopping rules of the solvers of physics/stagGpu against the
## CPU operators, stag.D and stagD2ee, on random SU(3) links, for masses of
## both signs and sources on all sites (full) or on the even ones, whose
## odd sites hold 1e30 to show they are not read:
##   solveM   |b - M x|^2 <= r2req |b|^2, M = m + D/2, of one system and of
##            several (nBatch per hop), in double and mixed precision
##   solveEE  |b_e - A x_e|^2 <= r2req |b_e|^2
## with fixed and atomic dot products.  solveM records the residual it
## reaches, which the CPU one matches, leaves the source unchanged, stops
## within maxits at a residual above the request when maxits is too small,
## and gives x = 0 for b = 0.  Exits with 1 if a check fails.
##   -lat: lattice (default 4^4 per rank), -reals: links of 18 or 12 reals,
##   -ipc:0 MPI for peers on the node too, -split: second hop split around
##   the exchange, 1, 0, -1 with off-node neighbors, -verb: verbosity;
##   -cpu:0 skips the comparisons with the production double/mixed solvers;
##   -inner:0 skips the raw single precision CG comparisons; -innerOnly:1
##   runs only those checks. -r2in sets the mixed inner tolerance; -innerReq
##   overrides the raw CG sweep, -r2req overrides final fixture tolerances.
##   -innerMax and -innerMatch bound true inner r2 and batch/scalar dx2/x2.
##   -refRound and -refFactor bound FP64 CPU-reference roundoff separately
##   from the solver's strict stopping target.
##   -fwd selects link caching; -seed selects the gauge and sources.
##   -recon64:1 uses FP64 intermediates to reconstruct compressed FP32 links.
import qex
import physics/[qcdTypes, stagSolve]
import solvers/cg as cpuCg
# Exercise the same private CGs used by mixed solveM and solveEE.
include physics/stagGpu
import comms/halogpu
import backend/accel
import strformat

qexInit()
let lat = intSeqParam("lat", latticeFromLocalLattice(@[4,4,4,4], nRanks))
let lo = lat.newLayout
echo "lattice ", lo.physGeom, "  ranks ", lo.rankGeom, "  lanes ", lo.innerGeom, "  nBatch ", nBatch
var g = lo.newGauge
let seed = intParam("seed", 987654321).uint
var rng = newRngField(lo, MRG32k3a, seed)
g.random rng
threads:
  g.projectSU
  threadBarrier()
  g.setBC
  threadBarrier()
  g.stagPhase
var s = newStag(g)
haloIpc = intParam("ipc", 1) != 0
hopSplit = intParam("split", -1)
var sd = newStagGpu(g, float64, intParam("reals", 18), intParam("fwd", -1), batch = true)
var sf = newStagGpu(g, float32, intParam("reals", 18), intParam("fwd", -1), batch = true,
                    recon64 = intParam("recon64", 0) != 0)
echo "GPU links: ", sd.nl, " reals", "  split hop: ", sd.split, "  peer IPC: ", haloIpc,
  "  FP32 recon64: ", sf.recon64
let n6 = 6*sd.n
let ne6 = 6*sd.ne
const big = 100000  # maxits of the solves that converge
let verb = intParam("verb", 0)
let cpu = intParam("cpu", 1) != 0
let inner = intParam("inner", 1) != 0
let innerReq = floatParam("r2in", 1e-6)
let innerMax = floatParam("innerMax", 0.25) # hard cases must reduce the true residual norm by at least two
let innerMatch = floatParam("innerMatch", 1e-4) # squared vector difference, batch versus scalar
let innerGoal = floatParam("innerReq", 0.0) # 0 selects the regression sweep
let outerGoal = floatParam("r2req", 0.0) # 0 selects each fixture's final tolerance
let refRound = floatParam("refRound", 128.0)
let refFactor = floatParam("refFactor", 4.0) # at most twice the requested residual norm
doAssert refRound >= 0 and refFactor >= 1
proc requested(r: float): float = (if outerGoal > 0: outerGoal else: r)
echo "solver tests seed ", seed, " r2in ", innerReq, " innerMax ", innerMax,
  " innerMatch ", innerMatch, " innerReq ", innerGoal, " r2req ", outerGoal,
  " refRound ", refRound, " refFactor ", refFactor
var fails = 0

proc referenceLimit(req, opNorm, x2b2: float): float =
  ## The unitary-link fixtures have ||M|| <= |m|+4 and ||A|| <= 4m^2+64.
  ## Compare different FP64 representations/evaluation orders with a roundoff
  ## allowance scaled by these operator bounds and ||x||/||b||. It is independent
  ## of measured residuals and capped by refFactor; the own-operator target is
  ## always req. The conservative refRound covers reconstruction and both hops.
  let rounding = refRound*epsilon(float64)*opNorm*sqrt(x2b2)
  min(refFactor*req, (sqrt(req)+rounding)^2)

proc effective[V: static int](sg: StagGpu[V,float32]): auto =
  ## Decode the actual FP32 forward links to the host FP64 reference.
  ## Cache adjoint consistency is tested separately by tstaglinks.
  result = sg.lo.newGauge()
  let n = sg.n
  let lf = sg.lf
  let nb = sg.nbr
  let tmp = cast[ptr UncheckedArray[float]](gpuMalloc(4*18*n*sizeof(float)))
  template expand(N: static int; r64: static bool) =
    gpuFor(i, 4*n):
      let mu = i div n
      let k = i mod n
      let e = if nb[mu*n+k] < 0: -1'f32 else: 1'f32
      var a: array[18,float32]
      load(a, lf, N*mu*n + uo(V,k,N), V, e, N, r64)
      let o = 18*mu*n + uo(V,k,18)
      forStatic c, 0, 17: tmp[o+c*V] = float(a[c])
  template comp(r64: static bool) =
    if sg.nl == 12: expand(12,r64) else: expand(14,r64)
  if sg.nl == 18: expand(18,false)
  elif sg.recon64: comp(true)
  else: comp(false)
  for mu in 0..3:
    gpuMemCpyToCpu(addr result[mu][0], addr tmp[18*mu*n], 18*n*sizeof(float))
  gpuFree(tmp)

proc checkInner(m, r2req: float) =
  ## Every batch slot gets a distinct source and mass.  Accuracy is checked
  ## independently of the recursive stopping statistic, on each output.
  var ss = toSingle(s)
  var gd = lo.newGauge()
  for mu in 0..<gd.len:
    threads: gd[mu] := ss.g[mu]
  var dc = newStag(gd)
  var de = newStag(effective(sf))
  type F = type(lo.ColorVectorS())
  var bs, xc, xg, xb = newSeq[F](nBatch)
  var b2, rc = newSeq[float](nBatch)
  var ic, ig = newSeq[int](nBatch)
  var rg = newSeq[float](nBatch)
  var ms = newSeq[float](nBatch)
  var ri = newRngField(lo, MRG32k3a, seed xor 192837465'u)
  for j in 0..<nBatch:
    bs[j] = lo.ColorVectorS()
    xc[j] = lo.ColorVectorS()
    xg[j] = lo.ColorVectorS()
    xb[j] = lo.ColorVectorS()
    threads:
      bs[j].gaussian ri
      threadBarrier()
      bs[j].odd := 0
      xc[j] := 0
      xg[j] := 0
      xb[j] := 0
      threadBarrier()
      let bb = bs[j].even.norm2
      threadMaster: b2[j] = bb
    ms[j] = m*(1.0+0.2*float(j))
    let mj = ms[j]
    var sp = initSolverParams()
    sp.backend = sbQex
    sp.r2req = r2req
    sp.maxits = big
    sp.verbosity = 0
    sp.sloppySolve = SloppyNone
    sp.subset.layoutSubset(lo, "even")
    proc op(a,b: F) =
      threadBarrier()
      stagD2ee(ss.se, ss.so, a, ss.g, b, mj*mj)
    var st = cpuCg.newCgState(xc[j], bs[j])
    cpuCg.solve(st, (apply: op, precon: cpuCg.cpNone), sp)
    ic[j] = sp.iterations
    threads:
      let rr = st.r.even.norm2
      threadMaster: rc[j] = rr/b2[j]
    sf.upload(sf.vec[5], bs[j])
    let z = sf.cg(sf.vec[6], sf.vec[5], mj, r2req*b2[j], big, verb, verify=false)
    ig[j] = z.its
    rg[j] = z.r2/b2[j]
    sf.download(xg[j], sf.vec[6])
  let bd = sysB(sf.vecB[5], 6*sf.ne)
  let xd = sysB(sf.vecB[6], 6*sf.ne)
  var st = sf.initB
  for j in 0..<nBatch:
    sf.upload(bd[j], bs[j])
    sf.startB(st, j, xd[j], bd[j], ms[j], r2req*b2[j], big)
  while st.act != 0: discard sf.stepB(st, verb)
  for j in 0..<nBatch: sf.download(xb[j], xd[j])

  proc truth(x,b: F; op: type(s); m: float): float =
    var xx, bb, rr = lo.ColorVectorD()
    var r2 = 0.0
    threads:
      xx := x
      bb := b
      threadBarrier()
      stagD2ee(op.se, op.so, rr, op.g, xx, m*m)
      threadBarrier()
      rr.even := bb - rr
      let r = rr.even.norm2
      let b = bb.even.norm2
      threadMaster: r2 = r/b
    r2
  proc distance(x,y: F): float =
    var d = 0.0
    threads:
      let a = norm2diff(x.even,y.even)
      let b = y.even.norm2
      threadMaster: d = a/b
    d
  proc accepts(r, dx, m: float): bool =
    let lim = if abs(m) >= 0.1: max(4*r2req, 1e-12) else: max(4*r2req, innerMax)
    r <= lim and dx <= innerMatch
  var zero = lo.ColorVectorS()
  threads: zero := 0
  for j in 0..<nBatch:
    let mj = ms[j]
    let rb = st.r2[j]/b2[j]
    let tc = truth(xc[j],bs[j],dc,mj)
    let tg = truth(xg[j],bs[j],de,mj)
    let tb = truth(xb[j],bs[j],de,mj)
    let orig = truth(xb[j],bs[j],s,mj)
    let dx = distance(xb[j],xg[j])
    var ok = rc[j] <= r2req and rg[j] <= r2req and rb <= r2req
    ok = ok and ic[j] <= big and ig[j] <= big and st.its[j] <= big
    ok = ok and accepts(tc,0.0,mj) and accepts(tg,0.0,mj) and accepts(tb,dx,mj)
    # These controls must fail the same output checks as a real damaged copy.
    ok = ok and not accepts(truth(zero,bs[j],de,mj),distance(zero,xg[j]),mj)
    let wrong = if nBatch > 1: xb[(j+1) mod nBatch] else: bs[j]
    ok = ok and not accepts(truth(wrong,bs[j],de,mj),distance(wrong,xg[j]),mj)
    echo &"inner slot {j} m {mj:.6g} requested {r2req:.3e} iterations CPU/GPU/batch {ic[j]}/{ig[j]}/{st.its[j]}",
         if ok: "  ok" else: "  FAILED"
    echo &"  recursive CPU/GPU/batch {rc[j]:.6e}/{rg[j]:.6e}/{rb:.6e} effective FP64 {tc:.6e}/{tg:.6e}/{tb:.6e}",
         &" original FP64 {orig:.6e} batch/scalar dx2/x2 {dx:.6e}"
    if not ok: inc fails
  free(ss)
  free(dc)
  free(de)

proc dvec(): ptr UncheckedArray[float] = cast[ptr UncheckedArray[float]](gpuMalloc(n6*sizeof(float)))

type Src = object
  ## a source on the device and its host copy
  d: ptr UncheckedArray[float]
  h: seq[float]
  full: bool

proc newSrc(scale: float; full: bool; odd = false): Src =
  ## scale times a gaussian source, b_o = 1e30 unless full
  var f = lo.ColorVector()
  threads: f.gaussian rng
  result.h = newSeq[float](n6)
  copyMem(addr result.h[0], addr f[0], n6*sizeof(float))
  for i in 0..<n6:
    result.h[i] = if odd and i < ne6: 0.0 elif full or i < ne6: scale*result.h[i] else: 1e30
  result.d = dvec()
  gpuMemCpyToGpu(result.d, addr result.h[0], n6*sizeof(float))  # not last: a cint with OpenMP
  result.full = full

proc resid(x: ptr UncheckedArray[float]; b: Src; m: float): tuple[r2, x2b2: float] =
  ## |b - M x|^2/|b|^2 by stag.D, b_o = 0 unless full; |M x|^2 for b = 0
  var hx, hb, r = lo.ColorVector()
  gpuMemCpyToCpu(addr hx[0], x, n6*sizeof(float))
  copyMem(addr hb[0], unsafeAddr b.h[0], n6*sizeof(float))
  var rr, bb, xx = 0.0
  threads:
    if not b.full: hb.odd := 0
    threadBarrier()
    s.D(r, hx, m)
    threadBarrier()
    r := hb - r
    let t = r.norm2
    let u = hb.norm2
    let v = hx.norm2
    threadMaster:
      rr = t
      bb = u
      xx = v
  if bb > 0: (rr/bb, xx/bb) else: (rr, xx)

proc unchanged(b: Src): bool =
  ## the device source still holds its host copy
  var h = newSeq[float](n6)
  gpuMemCpyToCpu(addr h[0], b.d, n6*sizeof(float))
  var d = 0.0
  for i in 0..<n6:
    if h[i] != b.h[i]: d += 1
  getDefaultComm().allReduce(d)
  d == 0

proc actualM(x: ptr UncheckedArray[float]; b: Src; m: float): float =
  ## Re-evaluate the returned solution with the operator used by solveM.
  let r = sd.vec[7]
  let src = b.d
  let full = b.full
  let ne = 6*sd.ne
  sd.applyMfull(r,x,m)
  gpuFor(i,n6): r[i] = (if full or i < ne: src[i] else: 0.0) - r[i]
  let rr = sd.norm2(r)
  let bb = sd.norm2EO(src,full)
  let b2 = bb[0]+bb[1]
  if b2 > 0: rr/b2 else: rr

proc compare(label: string; x: ptr UncheckedArray[float]; b: Src; m, r2req, gpuRefLimit: float; mixed: bool) =
  ## M is normal with sigma_min >= |m|. The residual limits of the two
  ## solutions give their distance bound by the triangle inequality.
  var hb, xc, xg, r = lo.ColorVector()
  copyMem(addr hb[0], unsafeAddr b.h[0], n6*sizeof(float))
  gpuMemCpyToCpu(addr xg[0], x, n6*sizeof(float))
  threads:
    if not b.full: hb.odd := 0
  var sp = initSolverParams()
  sp.backend = sbQex
  sp.r2req = requested(r2req)
  sp.maxits = big
  sp.verbosity = 0
  sp.sloppySolve = if mixed: SloppySingle else: SloppyNone
  s.solve(xc, hb, m, sp)
  var r2, dx, bb = 0.0
  threads:
    s.D(r, xc, m)
    threadBarrier()
    r := hb - r
    let rr = r.norm2
    let b2 = hb.norm2
    r := xg - xc
    let d2 = r.norm2
    threadMaster:
      bb = b2
      r2 = if b2 > 0: rr/b2 else: rr
      dx = if b2 > 0: d2/b2 else: d2
  let dxLimit = (sqrt(r2req)+sqrt(gpuRefLimit))^2
  let ok = r2 <= 1.001*r2req and m*m*dx <= 1.001*dxLimit and sp.calls == 1 and sp.iterations <= big and sp.r2.mean <= r2req
  echo &"{label:<28} production its {sp.iterations:5}  |b-Mx|^2/|b|^2 {r2:9.3e}  |xg-xc|^2/|b|^2 {dx:9.3e}",
       if ok: "  ok" else: "  FAILED"
  if not ok: inc fails

proc check(label: string; x: ptr UncheckedArray[float]; b: Src; m: float; sp: SolverParams;
           r2req: float; maxits: int; zero = false; mixed = false) =
  ## the solve of sp, its only one: with maxits = big it meets r2req, with
  ## fewer it misses it within maxits. Check the statistic against the same
  ## operator, and final accuracy against the independent CPU operator.
  let (r2, x2b2) = resid(x, b, m)
  let refLimit = referenceLimit(r2req, abs(m)+4.0, x2b2)
  let own = actualM(x,b,m)
  let rec = sp.r2.max
  var ok = unchanged(b) and sp.calls == 1 and sp.r2.n == 1 and sp.iterations <= maxits and sp.iterationsMax == sp.iterations
  if zero: ok = ok and r2 == 0 and own == 0 and rec == 0 and sp.iterations == 0
  else:
    ok = ok and abs(rec-own) <= 1e-6*max(rec,own) + 1e-10*r2req and (rec <= r2req) == (maxits == big)
    if maxits == big: ok = ok and own <= r2req and r2 <= refLimit
  echo &"{label:<28} m {m:6} its {sp.iterations:5}  |b-Mx|^2/|b|^2 recorded {rec:9.3e} recomputed {own:9.3e} CPU {r2:9.3e}  requested {r2req:7.1e} CPU limit {refLimit:.3e}",
       if ok: "  ok" else: "  FAILED"
  if not ok: inc fails
  if cpu and maxits == big: compare(label, x, b, m, r2req, refLimit, mixed)

proc solveOne(label: string; m, r2req: float; maxits: int; full, mixed: bool; scale = 1.0;
              r2in = innerReq; odd = false) =
  let b = newSrc(scale, full, odd)
  let x = dvec()
  var sp = initSolverParams()
  sp.r2req = requested(r2req)
  sp.maxits = maxits
  sp.verbosity = verb
  sd.solveM(x, b.d, m, sp, if mixed: addr sf else: nil, r2in = r2in, full = full)
  check(label, x, b, m, sp, sp.r2req, maxits, scale == 0.0, mixed)
  gpuFree(x)
  gpuFree(b.d)

proc solveMany(label: string; ns: int; full, mixed: bool; r2in = innerReq) =
  ## ns systems of different masses and requests, system 3 with 5
  ## iterations only, system 5 with b = 0
  const masses = [0.1, -0.02, 0.05, 0.01, -0.2, 0.03, 0.005, 0.07, -0.04]
  const requests = [1e-12, 1e-10, 1e-14, 1e-12, 1e-16, 1e-12, 1e-20, 1e-12, 1e-14]
  var ms, rq = newSeq[float](ns)
  var bs: seq[Src]
  var xs, bd: seq[ptr UncheckedArray[float]]
  var sps = newSeq[SolverParams](ns)
  for j in 0..<ns:
    ms[j] = masses[j mod masses.len]*(1.0+0.03*float(j div masses.len))
    rq[j] = requests[j mod requests.len]
    bs.add newSrc(if j == 5: 0.0 else: 1.0, full)
    bd.add bs[j].d
    xs.add dvec()
    sps[j] = initSolverParams()
    sps[j].r2req = requested(rq[j])
    sps[j].maxits = if j == 3: 5 else: big
    sps[j].verbosity = verb
  sd.solveM(xs, bd, ms, sps, if mixed: addr sf else: nil, r2in = r2in, full = full)
  for j in 0..<ns:
    check(&"{label} {j}/{ns}", xs[j], bs[j], ms[j], sps[j], sps[j].r2req, sps[j].maxits, j == 5, mixed)
    gpuFree(xs[j])
    gpuFree(bs[j].d)

proc solveEEcheck(label: string; m, r2req: float; mixed: bool; maxits = big; reps = 1; zero = false;
                  r2in = innerReq) =
  ## the true residual of solveEE by stagD2ee
  let r2req = requested(r2req)
  var src = lo.ColorVector()
  var xg = lo.ColorVector()
  var d = lo.ColorVector()
  threads:
    src.gaussian rng
    threadBarrier()
    if zero: src := 0
    src.odd := 0
    xg := 0
  var sp = initSolverParams()
  sp.r2req = requested(r2req)
  sp.maxits = maxits
  sp.verbosity = verb
  for k in 0..<reps:
    if mixed: solveEE(sd, sf, xg, src, m, sp, r2in)
    else: sd.solveEE(xg, src, m, sp)
  sd.upload(sd.vec[5],src)
  sd.upload(sd.vec[6],xg)
  sd.applyD2ee(sd.vec[1],sd.vec[6],sd.vec[2],m*m)
  sd.resid(sd.vec[0],sd.vec[5],sd.vec[1])
  var rr = sd.redot(sd.vec[0],sd.vec[0])
  var bb = sd.redot(sd.vec[5],sd.vec[5])
  getDefaultComm().allReduce(rr)
  getDefaultComm().allReduce(bb)
  let own = if bb > 0: rr/bb else: rr
  var r2, x2b2 = 0.0
  threads:
    stagD2ee(s.se, s.so, d, s.g, xg, m*m)
    threadBarrier()
    d.even -= src
    let t = d.even.norm2
    let u = src.even.norm2
    let v = xg.even.norm2
    threadMaster:
      r2 = if u > 0: t/u else: t
      x2b2 = if u > 0: v/u else: v
  let refLimit = referenceLimit(r2req, 4*m*m+64.0, x2b2)
  var ok = sp.calls == reps and sp.r2.n == reps and sp.iterations <= reps*maxits and sp.iterationsMax <= maxits
  ok = ok and abs(sp.r2.max - own) <= 1e-6*max(sp.r2.max, own) + 1e-10*r2req
  if zero: ok = ok and r2 == 0 and sp.r2.max == 0 and sp.iterations == 0
  elif maxits == big: ok = ok and own <= r2req and sp.r2.max <= r2req and r2 <= refLimit
  else: ok = ok and r2 > r2req
  echo &"{label:<28} m {m:6} its {sp.iterations:5}  |b-Ax|^2/|b|^2 recorded {sp.r2.max:9.3e} recomputed {own:9.3e} CPU {r2:9.3e}  requested {r2req:7.1e} CPU limit {refLimit:.3e}",
       if ok: "  ok" else: "  FAILED"
  if not ok: inc fails
  if cpu and maxits == big:
    var xc = lo.ColorVector()
    var spc = initSolverParams()
    spc.backend = sbQex
    spc.r2req = r2req
    spc.maxits = big
    spc.verbosity = 0
    spc.sloppySolve = if mixed: SloppySingle else: SloppyNone
    s.solveEE(xc, src, m, spc)
    var rc, dx = 0.0
    threads:
      stagD2ee(s.se, s.so, d, s.g, xc, m*m)
      threadBarrier()
      d.even -= src
      let rr = d.even.norm2
      let bb = src.even.norm2
      d.even := xg - xc
      let dd = d.even.norm2
      threadMaster:
        rc = if bb > 0: rr/bb else: rr
        dx = if bb > 0: dd/bb else: dd
    let dxLimit = (sqrt(r2req)+sqrt(refLimit))^2
    let ok = rc <= 1.001*r2req and 16*m*m*m*m*dx <= 1.001*dxLimit and spc.calls == 1 and spc.iterations <= big and spc.r2.mean <= r2req
    echo &"{label:<28} production its {spc.iterations:5}  |b-Ax|^2/|b|^2 {rc:9.3e}  |xg-xc|^2/|b|^2 {dx:9.3e}",
         if ok: "  ok" else: "  FAILED"
    if not ok: inc fails

for fixed in [true, false]:
  sd.fixed = fixed
  sf.fixed = fixed
  let f = if fixed: "fixed" else: "atomic"
  if inner:
    for m in [0.1, 0.01, 0.001]:
      for rq in [1e-3, 1e-6, float(epsilon(float32))]:
        if innerGoal == 0 or rq == 1e-3: checkInner(m, if innerGoal > 0: innerGoal else: rq)
  if intParam("innerOnly",0) != 0: continue
  for mixed in [false, true]:
    let p = if mixed: "mixed" else: "double"
    for full in [true, false]:
      let e = if full: "full" else: "even"
      for m in [0.1, 0.01, 0.001, -0.01]:
        solveOne(&"{f} {p} {e}", m, 1e-12, big, full, mixed)
      solveOne(&"{f} {p} {e} tight", 0.001, 1e-24, big, full, mixed)
      solveOne(&"{f} {p} {e} maxits 5", 0.01, 1e-12, 5, full, mixed)
      solveOne(&"{f} {p} {e} zero", 0.01, 1e-12, big, full, mixed, 0.0)
      for ns in [2, nBatch, nBatch + 1, 2*nBatch + 1]:
        solveMany(&"{f} {p} {e} batch", ns, full, mixed)
    solveOne(&"{f} {p} odd source", 0.01, 1e-16, big, true, mixed, odd = true)
    solveOne(&"{f} {p} maxits 0", 0.01, 1e-12, 0, true, mixed)
    solveEEcheck(&"{f} {p} solveEE", 0.01, 1e-14, mixed)
    solveEEcheck(&"{f} {p} solveEE", 0.001, 1e-20, mixed)
    solveEEcheck(&"{f} {p} solveEE stats", 0.01, 1e-14, mixed, maxits = 5, reps = 2)
    solveEEcheck(&"{f} {p} solveEE zero", 0.01, 1e-14, mixed, reps = 2, zero = true)
    solveEEcheck(&"{f} {p} solveEE maxits 0", 0.01, 1e-14, mixed, maxits = 0)
  for r2in in [0.0, 1e-12, 1e-3]:
    solveOne(&"{f} mixed r2in {r2in}", 0.01, 1e-16, big, true, true, r2in = r2in)
    solveMany(&"{f} mixed r2in {r2in}", 2*nBatch + 1, true, true, r2in)
    solveEEcheck(&"{f} mixed r2in {r2in}", 0.01, 1e-16, true, r2in = r2in)

echo if fails == 0: "all solves ok" else: $fails & " solves FAILED"
qexExit(if fails == 0: 0 else: 1)
