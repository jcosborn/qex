import metropolis
import strUtils

type
  HmcRoot* = object of MetropolisRoot
    lwRev: float

proc checkReverse*[T:HmcRoot](h: var T) =
  template ff(x: float): string =
    formatFloat(x, ffDecimal, precision=6)
  h.startReverse
  h.generate
  h.lwRev = h.logWeight
  if h.verbosity>0:
    echo "ln(wOld): $1  ln(wRev): $2  diff: $3"%
        [h.lwOld.ff, h.lwRev.ff, (h.lwRev-h.lwOld).ff]
    h.finishReverse

when isMainModule:
  import random
  type
    HMC = object of HmcRoot
      x: float
      xSave: float
      xRev: float
      p: float
      nSteps: int

  proc init*(h: var HMC) =
    clear(HmcRoot(h))

  proc start(h: var HMC) =
    h.xSave = h.x
    h.p = 2.0*rand(1.0) - 1.0  # should really be Gaussian
    echo "x: ", h.x, "  p: ", h.p

  template startReverse(h: HMC) =
    h.p = -h.p
    h.xRev = h.x
  template finishReverse(h: HMC) =
    h.x = h.xRev
  template finish(h: HMC) =
    let v = h.verbosity
    h.verbosity = 1
    h.checkReverse
    h.verbosity = v

  proc logWeight(h: HMC): float =
    let x = h.x
    let p = h.p
    1e10 + 0.5*p*p + x*x

  proc accept(h: var HMC) =
    echo "x: ", h.x, "  p: ", h.p
  proc reject(h: var HMC) =
    h.x = h.xSave

  proc globalRand(h: HMC): float =
    rand(1.0)

  proc updateX(h: var HMC, e: float) =
    h.x += e * h.p
  proc updateP(h: var HMC, e: float) =
    h.p -= e * 2.0 * h.x
  proc generate(h: var HMC) =
    let n = h.nSteps
    let tau = 1.0
    let eps = tau / n.float
    h.updateX(0.5*eps)
    h.updateP(eps)
    for i in 2..n:
      h.updateX(eps)
      h.updateP(eps)
    h.updateX(0.5*eps)

  var h: HMC
  h.init
  h.x = 1.4
  h.nSteps = 2
  randomize(987654321)

  while h.nUpdates < 10:
    h.update
    echo "nSteps: ", h.nSteps, "  accRatio: ", h.accRatio
    if h.accepted:
      h.nSteps = max(2,h.nSteps-2)
    else:
      h.nSteps += 2
