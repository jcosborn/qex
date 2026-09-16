#RUNCMD env OMP_NUM_THREADS=1 $RUNJOB

import base/globals
setVLENmax(4)

import math, unittest

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))

# basicOps.epsilon collides with fenv.epsilon
import qex except epsilon
import algorithms/numdiff, gauge/stoutsmear
import helpers
import ../[core, scalar, gauge]
from ../gauge/matfun import expPolyGraph, expJet, expTopReplica
import ../gauge/[types, field_ops, transport, cfield, stencil]

let grt = initGraphRuntime()

include gauge/gaugehelpers

# Path equivalences and higher derivatives build reference graphs behind the
# fused kernels; keep the lattice small. Run this alongside tggauge and
# tgtoweru1 (Nc=1) to cover both the path comparisons and derivative towers.
proc runTowerTests*(localLat: seq[int], seed: uint64, subDir: int,
                    smearSteps: int) =
  qexInit()
  letParam:
    expectRanks = nRanks
  check nRanks == expectRanks
  echo "tower ranks: ", nRanks
  letParam:
    lat = latticeFromLocalLattice(localLat,nRanks)
  let lo = lat.newLayout
  var
    r = lo.newRNGField(Philox4x64, seed)
    g = lo.newgauge
    u = lo.newgauge
    q = lo.newgauge
    m = lo.newgauge
    ss = lo.newStoutSmear(0.1)
  threads:
    g.random r
    u.random r
    q.randomTAH r
    m.randomTAH r
  for _ in 0..<smearSteps:
    ss.smear(g, g)
    ss.smear(u, u)
  threads:
    for t in m:
      t *= 0.1
    for t in q:
      t *= 0.1

  include gauge/transport_basic
  include gauge/higher

  qexFinalize()

when isMainModule:
  # Smearing matches the tggauge fixtures, so both test the same kind of configuration.
  runTowerTests(@[4,4,4,4], 987654321u64, 2, 5)
