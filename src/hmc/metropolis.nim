import strUtils
import comms/commsUtils, math, strformat, stats

type
  MetropolisStats* = object
    hOld*: float
    hNew*: float
    rnd*: float
  MetropolisRootObj* {.inheritable.} = object
    verbosity*: int  # 0:quiet, 1:info
    stats*: seq[MetropolisStats]
    hOld*: float
    hNew*: float
    rnd*: float
    deltaH*: float
    hReverse*: float
    expmDeltaH*: float
    pAccept*: float
    accepted*: bool
    nUpdates*: int
    nAccepts*: int
    nRejects*: int
    acceptRatio*: float
    avgDeltaH*: float
    avgDeltaH2*: float
    avgPAccept*: float
    pAcceptStats*: RunningStat
  MetropolisRoot* = ref MetropolisRootObj
  MetropolisWrapper*[T] = ref object of MetropolisRootObj
    state*: T
# required routines:
#   start, logWeight, generate, globalRand, accept, reject

# optional routines
proc finish*[M:MetropolisRoot](m: var M) = discard
proc checkReverse*[M:MetropolisRoot](m: var M): bool = false
proc generateReverse*[M:MetropolisRoot](m: var M) = discard
proc finishReverse*[M:MetropolisRoot](m: var M) = discard


proc clearStats*[M:MetropolisRoot](m: var M) =
  m.stats.setLen(0)
  m.nUpdates = 0
  m.nAccepts = 0
  m.nRejects = 0
  m.acceptRatio = 0
  m.avgDeltaH = 0
  m.avgDeltaH2 = 0
  m.avgPAccept = 0
  clear m.pAcceptStats

proc updateStats*[M:MetropolisRoot](m: var M) =
  m.stats.add MetropolisStats(hOld:m.hOld,hNew:m.hNew,rnd:m.rnd)
  let n = m.nUpdates.float
  inc m.nUpdates
  if m.accepted:
    inc m.nAccepts
  else:
    inc m.nRejects
  m.acceptRatio = m.nAccepts / m.nUpdates
  m.avgDeltaH = (n*m.avgDeltaH + m.deltaH) / (n+1)
  m.avgDeltaH2 = (n*m.avgDeltaH2 + m.deltaH*m.deltaH) / (n+1)
  m.avgPAccept = (n*m.avgPAccept + m.pAccept) / (n+1)
  m.pAcceptStats.push m.pAccept

proc init*[M:MetropolisRoot](m: var M) =
  m.verbosity = 0
  m.stats.newSeq(0)
  m.clearStats

proc update*[T:MetropolisRoot](m: var T) =
  mixin finish
  template ff(x: float): string =
    formatFloat(x, ffDecimal, precision=6)

  m.start

  m.hOld = m.getH
  m.generate
  m.hNew = m.getH
  m.deltaH = m.hNew - m.hOld
  if m.verbosity>0:
    echo &"hOld: {m.hOld:.6f}  hNew: {m.hNew:.6f}"

  m.finish

  if m.checkReverse:
    m.generateReverse
    m.hReverse = m.getH
    m.finishReverse
    # echo?

  m.rnd = m.globalRand
  m.expmDeltaH = exp(-m.deltaH)
  m.pAccept = min(1.0, m.expmDeltaH)
  if m.rnd <= m.pAccept:
    m.accepted = true
    m.updateStats
    if m.verbosity>0:
      echo "ACCEPT deltaH: $1  pAccept: $2  rnd: $3"%
        [m.deltaH.ff, m.pAccept.ff, m.rnd.ff]
    m.accept
  else:
    m.accepted = false
    m.updateStats
    if m.verbosity>0:
      echo "REJECT deltaH: $1  pAccept: $2  rnd: $3"%
        [m.deltaH.ff, m.pAccept.ff, m.rnd.ff]
    m.reject


when isMainModule:
  import random

  type
    Met = ref object of MetropolisRoot
      x: float
      xSave: float
      p: float
      nSteps: int

  proc init*(m: var Met) =
    m.new
    var r = MetropolisRoot m
    init(r)

  proc start(m: var Met) =
    m.xSave = m.x
    m.p = 2.0*rand(1.0) - 1.0  # should really be Gaussian
    echo "start: x: ", m.x, "  p: ", m.p

  proc finish(m: var Met) =
    echo "finish: x: ", m.x, "  p: ", m.p

  proc accept(m: var Met) = discard

  proc reject(m: var Met) =
    m.x = m.xSave

  proc globalRand(m: Met): float =
    rand(1.0)

  proc getH(m: Met): float =
    let x = m.x
    let p = m.p
    1e10 + 0.5*p*p + x*x

  proc updateX(m: var Met, e: float) =
    m.x += e * m.p
  proc updateP(m: var Met, e: float) =
    m.p -= e * 2.0 * m.x
  proc generate(m: var Met) =
    let n = m.nSteps
    let tau = 1.0
    let eps = tau / n.float
    m.updateX(0.5*eps)
    m.updateP(eps)
    for i in 2..n:
      m.updateX(eps)
      m.updateP(eps)
    m.updateX(0.5*eps)

  var m: Met
  m.init
  m.verbosity = 1
  m.x = -1.4
  m.nSteps = 2
  randomize(987654321)

  while m.nUpdates < 10:
    m.update
    echo "nSteps: ", m.nSteps, "  accRatio: ", m.acceptRatio
    if m.accepted:
      m.nSteps = max(2,m.nSteps-2)
    else:
      m.nSteps += 2
