#import base,
import layout, strformat, stats, base/params
export stats

type
  SolverBackend* = enum
    sbQex, sbQuda, sbGrid
  SloppyType* = enum
    SloppyNone, SloppySingle, SloppyHalf
  CgControl* = object
    kind*: int  # CPU -1 automatic; 0 existing, 1 conventional mixed; GPU 2 wide accumulator
    delta*: float  # reliable update factor in the residual norm; 0 disables the norm trigger
    acc64*, beta*, keep*: bool
    maxInc*, maxTotal*: int
    period*: int  # maximum iterations between precise residual updates; 0 uses only delta
  SolverParams* = object
    # inputs
    r2req*: float
    maxits*: int
    backend*: SolverBackend
    sloppySolve*: SloppyType
    cg*: CgControl
    r2in*: float  # minimum relative squared target of a CPU inner correction
    usePrevSoln*: bool
    verbosity*: int
    subset*: Subset
    subsetName*: string
    # outputs
    calls*: int
    iterations*: int
    iterationsMax*: int
    reliable*: int
    seconds*: float
    flops*: float
    #r2sum*: float
    #r2max*: float
    r2*: RunningStat
# ?split into SolverRequest, SolverFlags(?), SolverPerf?

template finalIterations*(sp: SolverParams): untyped = sp.iterations
template `finalIterations=`*(sp: var SolverParams, x: int): untyped =
  sp.iterations = x

proc cgControl*(prefix: string): CgControl =
  let cpu = prefix=="cpu"
  result.kind = intParam(prefix & "Cg", if cpu: -1 else: 0)
  result.delta = floatParam(prefix & "Delta", if cpu: 0.01 else: 0.1)
  result.acc64 = intParam(prefix & "Acc64", if cpu: 0 else: 1) != 0
  result.beta = intParam(prefix & "Beta", if cpu: 0 else: 1) != 0
  result.keep = intParam(prefix & "Keep", 1) != 0
  result.maxInc = intParam(prefix & "MaxInc", 1)
  result.maxTotal = intParam(prefix & "MaxTotal", 10)
  result.period = intParam(prefix & "Period", 1000)
  doAssert result.kind >= (if prefix=="cpu": -1 else: 0) and result.kind <= (if prefix=="cpu": 1 else: 2)
  doAssert result.delta >= 0.0 and result.delta < 1.0
  doAssert result.maxInc >= 0 and result.maxTotal >= 0
  doAssert result.period >= 0

proc init*(x: var RunningStat, y: SomeNumber) =
  x.clear()
  x.push y

proc resetStats*(sp: var SolverParams) =
  sp.calls = 0
  sp.iterations = 0
  sp.iterationsMax = 0
  sp.reliable = 0
  sp.seconds = 0.0
  sp.flops = 0.0
  #sp.r2sum = 0.0
  #sp.r2max = 0.0
  sp.r2.clear()

proc init*(sp: var SolverParams) =
  sp.r2req = 1e-6
  sp.maxits = 50000
  sp.backend = sbQex
  if defined(qudaDir): sp.backend = sbQuda
  if defined(gridDir): sp.backend = sbGrid
  sp.sloppySolve = intParam("sloppySolve", 0).SloppyType
  sp.cg = cgControl("cpu")
  sp.r2in = floatParam("cpuR2in", 0.0)
  doAssert sp.r2in >= 0.0 and sp.r2in < 1.0
  sp.usePrevSoln = false
  sp.verbosity = 1
  sp.subsetName = "all"
  sp.resetStats()

proc newSolverParams*: SolverParams =
  result.init()
template initSolverParams*: SolverParams = newSolverParams()

proc copyStats*(sp0: var SolverParams, sp1: SolverParams) =
  sp0.calls = sp1.calls
  sp0.iterations = sp1.iterations
  sp0.iterationsMax = sp1.iterationsMax
  sp0.reliable = sp1.reliable
  sp0.seconds = sp1.seconds
  sp0.flops = sp1.flops
  #sp0.r2sum = sp1.r2sum
  #sp0.r2max = sp1.r2max
  sp0.r2 = sp1.r2

proc addStats*(sp0: var SolverParams, sp1: SolverParams) =
  sp0.calls += sp1.calls
  sp0.iterations += sp1.iterations
  sp0.iterationsMax = max(sp0.iterationsMax, sp1.iterationsMax)
  sp0.reliable += sp1.reliable
  sp0.seconds += sp1.seconds
  sp0.flops += sp1.flops
  #sp0.r2sum += sp1.r2sum
  #sp0.r2max = max(sp0.r2max, sp1.r2max)
  sp0.r2 += sp1.r2

proc getStats*(sp: SolverParams, typ0= -1): string =
  let c = max(1,sp.calls)
  let its = sp.iterations
  let ic = its div c
  let im = sp.iterationsMax
  let s = max(1e-12, sp.seconds)
  let sc = s/c.float
  let f = sp.flops
  let gf = 1e-9*f/s
  #let gfc = gf/c.float
  #let r2 = sp.r2sum/c.float
  #let r2m = sp.r2max
  let r2 = sp.r2.mean
  let r2m = sp.r2.max
  #echo "op time: ", top
  var typ = typ0
  if typ<0: typ = if c==1: 3 else: 0
  case typ
  of 0:  # all
    result = &"{c}: {its}({ic}:{im})  {s:.4g}({sc:.4g})s" &
             &"  ({gf:.4g})Gf/s  ({r2:.4g}:{r2m:.4g})"
  of 1:  # total
    result = &"{c}: {its}  {s:.4g}s  {gf:.4g}Gf/s"
  of 2:  # ave
    result = &"{c}: {ic}:{im}  {sc:.4g}s  {gf:.4g}Gf/s  {r2:.4g}:{r2m:.4g}"
  of 3:  # single
    result = &"{ic}  {sc:.4g}s  {gf:.4g}Gf/s  {r2:.4g}"
  else:
    discard

template getTotalStats*(sp: SolverParams): untyped = getStats(sp,1)
template getAveStats*(sp: SolverParams): untyped = getStats(sp,2)
#template getStats*(sp: SolverParams): untyped = getStats(sp,2)

when isMainModule:
  var x = newSolverParams()
  x.r2.push 2
  echo x.r2.mean()
