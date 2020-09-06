import strUtils
import math

type
  MetropolisRoot* = object {.inheritable.}
    lwOld*: float    # log weight of starting configuration
    lwNew*: float    # log weight of ending configuration
    dlw*: float      # difference of log weights (lwOld-lwNew)
    pAccept*: float  # acceptance probability
    rnd*: float      # random value chosen for accept/reject
    accepted*: bool  # true if accepted
    nUpdates*: int   # number of updates so far
    nAccepts*: int   # number of accepts so far
    nRejects*: int   # number of rejects so far
    accRatio*: float # acceptance ratio so far
    verbosity*: int  # 0:quiet, 1:info
# required methods:
#   start, logWeight, generate, globalRand, accept, reject
# optional: finish

proc finish*(m: var MetropolisRoot) = discard

proc clear*(m: var MetropolisRoot) =
  m.accepted = false
  m.nUpdates = 0
  m.nAccepts = 0
  m.nRejects = 0
  m.accRatio = 0

proc update*[T:MetropolisRoot](m: var T) =
  mixin finish
  template ff(x: float): string =
    formatFloat(x, ffDecimal, precision=6)

  m.start

  m.lwOld = m.logWeight
  m.generate
  m.lwNew = m.logWeight
  m.dlw = m.lwOld - m.lwNew
  if m.verbosity>0:
    echo "ln(wOld): $1  ln(wNew): $2"%[m.lwOld.ff, m.lwNew.ff]

  m.finish

  m.rnd = m.globalRand
  m.pAccept = exp(m.dlw)
  inc m.nUpdates
  if m.rnd <= m.pAccept:
    m.accepted = true
    inc m.nAccepts
    m.accRatio = m.nAccepts.float / m.nUpdates.float
    if m.verbosity>0:
      echo "ACCEPT ln(pAccept): $1  pAccept: $2  rnd: $3"%
        [m.dlw.ff, m.pAccept.ff, m.rnd.ff]
    m.accept
  else:
    m.accepted = false
    inc m.nRejects
    m.accRatio = m.nAccepts.float / m.nUpdates.float
    if m.verbosity>0:
      echo "REJECT ln(pAccept): $1  pAccept: $2  rnd: $3"%
        [m.dlw.ff, m.pAccept.ff, m.rnd.ff]
    m.reject


when isMainModule:
  import random

  type
    Met = object of MetropolisRoot
      x: float
      xSave: float
      p: float
      nSteps: int

  proc init*(m: var Met) =
    clear(MetropolisRoot(m))

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

  proc logWeight(m: Met): float =
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
    echo "nSteps: ", m.nSteps, "  accRatio: ", m.accRatio
    if m.accepted:
      m.nSteps = max(2,m.nSteps-2)
    else:
      m.nSteps += 2
