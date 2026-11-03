## Complete M solves, scalar or batched, with independent CPU residuals.
## -lat -rg -mass or -masses -nb -r2req -r2in -maxits -mixed -full
## -ncpu -ngpu -reals -recon64 -reconFma -fwd -fixed -ipc -split -seed -verb -prof
## -cpuCg/-gpuCg: 0 existing, 1 reliable; CPU -1 automatic, GPU 2 FP64 accumulation
## -cpuR2in: CPU inner squared target
## -cpuDelta/-gpuDelta -cpuAcc64/-gpuAcc64 -cpuBeta/-gpuBeta -cpuKeep/-gpuKeep
## -d:cpuOnly builds the production reference without GPU modules.
## -d:stagWorkCount counts local stencil work and solver reductions.
import qex, physics/[qcdTypes, stagSolve]
import std/[sequtils, times, strformat]
when not defined(cpuOnly):
  import physics/stagGpu, backend/accel, comms/halogpu

qexInit()
let lat = intSeqParam("lat", @[8,8,8,8])
let rg = intSeqParam("rg", newSeq[int]())
let ig = intSeqParam("ig", newSeq[int]())
let lo = newLayout(lat, VLEN, rg, ig)
let mass = floatParam("mass", 0.1)
let ns = intParam("nb", 1)
let ms = floatSeqParam("masses", toSeq(0..<ns).mapIt(mass*(1.0+0.2*float(it))))
doAssert ns > 0 and ms.len == ns
let mixed = intParam("mixed",0) != 0
let full = intParam("full",1) != 0
let req = floatParam("r2req",1e-16)
let r2in = floatParam("r2in",1e-6)
let maxits = intParam("maxits",100000)
let seed = intParam("seed",987654321).uint
let verb = intParam("verb",0)
let ncpu = intParam("ncpu",1)
let ngpu = intParam("ngpu",3)
var rng = newRngField(lo,MRG32k3a,seed)
var g = lo.newGauge()
g.random rng
threads:
  g.projectSU
  threadBarrier()
  g.setBC
  threadBarrier()
  g.stagPhase
var s = newStag(g)
type F = type(lo.ColorVector())
var bs, xc, xg = newSeq[F](ns)
for j in 0..<ns:
  bs[j] = lo.ColorVector()
  xc[j] = lo.ColorVector()
  xg[j] = lo.ColorVector()
  threads:
    bs[j].gaussian rng
    threadBarrier()
    if not full: bs[j].odd := 0
echo "M benchmark lattice ", lat, " ranks ", lo.rankGeom, " lanes ", lo.innerGeom,
  " seed ", seed, " systems ", ns, " masses ", ms, " mixed ", mixed, " full ", full,
  " r2req ", req, " r2in ", r2in, " maxits ", maxits
var fails = 0

proc params(): SolverParams =
  result = initSolverParams()
  result.backend = sbQex
  result.r2req = req
  result.maxits = maxits
  result.verbosity = verb
  result.sloppySolve = if mixed: SloppySingle else: SloppyNone

proc residual(x,b: F; m: float): float =
  var r = lo.ColorVector()
  var rr = 0.0
  threads:
    s.D(r,x,m)
    threadBarrier()
    r := b-r
    let r2 = r.norm2
    let b2 = b.norm2
    threadMaster: rr = r2/b2
  rr

echo "CPU CG configuration ", params().cg, " inner r2 ", params().r2in
for rep in 0..<ncpu:
  var sp = newSeq[SolverParams](ns)
  for j in 0..<ns: sp[j] = params()
  resetTimers()
  getDefaultComm().barrier
  let t0 = epochTime()
  for j in 0..<ns: s.solve(xc[j],bs[j],ms[j],sp[j])
  var dt = epochTime()-t0
  getDefaultComm().max(dt)
  var its = 0
  var nres = 0
  var rmax = 0.0
  for j in 0..<ns:
    its += sp[j].iterations
    nres += sp[j].reliable
    let r = residual(xc[j],bs[j],ms[j])
    rmax = max(rmax,r)
    if not (r <= 1.01*req) or sp[j].iterations > maxits:
      inc fails
      echo "CPU slot ", j, " residual ", r, " iterations ", sp[j].iterations, " FAILED"
  echo &"BENCH cpu rep {rep} seconds {dt:.9g} iterations {its} true_r2max {rmax:.9e} reliable {nres}"
  echoProf()

when not defined(cpuOnly):
  haloIpc = intParam("ipc",1) != 0
  hopSplit = intParam("split",-1)
  let nl = intParam("reals",18)
  let r64 = intParam("recon64",0) != 0
  let rfm = intParam("reconFma",0) != 0
  let fw = intParam("fwd",-1)
  let fixed = intParam("fixed",1) != 0
  var sd = newStagGpu(g,float64,nl,fw,batch=ns>1)
  var ss: StagGpu[VLEN,float32]
  var sptr: ptr StagGpu[VLEN,float32]
  if mixed:
    ss = newStagGpu(g,float32,nl,fw,batch=ns>1,recon64=r64,reconFma=rfm)
    ss.fixed = fixed
    sptr = addr ss
  sd.fixed = fixed
  echo "GPU configuration reals ", sd.nl, " forward ", sd.lb == nil,
    " fixed ", fixed, " split ", sd.split, " IPC ", haloIpc, " FP32 recon64 ", r64, " reconFma ", rfm,
    " CG ", sd.ctrl
  let n6 = 6*sd.n
  var xd, bd = newSeq[ptr UncheckedArray[float]](ns)
  for j in 0..<ns:
    xd[j] = cast[ptr UncheckedArray[float]](gpuMalloc(n6*sizeof(float)))
    bd[j] = cast[ptr UncheckedArray[float]](gpuMalloc(n6*sizeof(float)))
    gpuMemCpyToGpu(bd[j],addr bs[j][0],n6*sizeof(float))
  for rep in 0..<ngpu:
    var sp = newSeq[SolverParams](ns)
    for j in 0..<ns: sp[j] = params()
    resetTimers()
    when declared(stagWork): stagWork = default(StagWork)
    getDefaultComm().barrier
    let t0 = epochTime()
    sd.solveM(xd,bd,ms,sp,sptr,r2in=r2in,full=full)
    var dt = epochTime()-t0
    getDefaultComm().max(dt)
    var its = 0
    var nres = 0
    var rmax = 0.0
    for j in 0..<ns:
      gpuMemCpyToCpu(addr xg[j][0],xd[j],n6*sizeof(float))
      its += sp[j].iterations
      nres += sp[j].reliable
      let r = residual(xg[j],bs[j],ms[j])
      rmax = max(rmax,r)
      if not (r <= 1.01*req) or sp[j].iterations > maxits or sp[j].calls != 1:
        inc fails
        echo "GPU slot ", j, " residual ", r, " iterations ", sp[j].iterations, " FAILED"
    echo &"BENCH gpu rep {rep} seconds {dt:.9g} iterations {its} true_r2max {rmax:.9e} reliable {nres}"
    when defined(stagWorkCount) and declared(stagWork):
      echo "WORK local ", stagWork
    echoProf()
  for j in 0..<ns:
    gpuFree(xd[j])
    gpuFree(bd[j])
  sd.free()
  if mixed: ss.free()
echo if fails == 0: "all benchmark residuals ok" else: $fails & " benchmark residuals FAILED"
qexExit(if fails == 0: 0 else: 1)
