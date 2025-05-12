import ../mcmcTypes
import ../mcmc/randomNumberGeneration
import ../utilities/gaugeUtils

import json
import parseopt
import streams
import strutils
import os

export json

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
  self.start = ReadGauge

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

proc write*[T](u: T; fn: string) =
  let gfn = fn & ".lat"
  if 0 != u.saveGauge(gfn): qexError "unable to write " & gfn

proc readJSON*(fn: string): JsonNode = fn.parseFile

proc readCMD*: JsonNode = 
  var cmd = initOptParser()
  result = parseJson("{}")
  while true:
    cmd.next()
    case cmd.kind:
      of cmdShortOption,cmdLongOption,cmdArgument:
        try: result[cmd.key] = %* parseInt(cmd.val)
        except ValueError:
          try: result[cmd.key] = %* parseFloat(cmd.val)
          except ValueError: result[cmd.key] = %* cmd.val
      of cmdEnd: break

proc setFilename*(info: auto; fn: string) =
  info["filename"] = %* fn 