import ../mcmcTypes
import ../mcmc/randomNumberGeneration

proc read*(self: var LatticeFieldTheory; fn: string; onlyGauge: bool = false) =
  let
    gfn = fn & ".lat"
    pfn = fn & ".rng"
    sfn = fn & ".global_rng"
  if fileExists(gfn):
    if 0 != self.u[].loadGauge(gfn): qexError "unable to read " & gfn
    reunit(self.u[])
  case self.algorithm:
    of HamiltonianMonteCarlo:
      if not onlyGauge:
        if fileExists(pfn): self.pRNG.readRNG(pfn)
        else: qexError "unable to read " & pfn
        if fileExists(sfn): 
          self.sRNG.readRNG(sfn)
          echo "read" & sfn
        else: qexError "unable to read " & sfn
    of HeatbathOverrelax: discard

proc write*(self: var LatticeFieldTheory; fn: string; onlyGauge: bool = false) =
  let
    gfn = fn & ".lat"
    pfn = fn & ".rng"
    sfn = fn & ".global_rng"
  if 0 != self.u[].saveGauge(gfn): qexError "unable to write " & gfn
  case self.algorithm:
    of HamiltonianMonteCarlo:
      if not onlyGauge:
        self.pRNG.writeRNG(pfn)
        self.sRNG.writeRNG(sfn)
    of HeatbathOverrelax: discard