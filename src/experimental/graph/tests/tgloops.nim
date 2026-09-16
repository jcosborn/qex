#RUNCMD env OMP_NUM_THREADS=1 $RUNJOB

import base/globals
setVLENmax(4)

import math, unittest
import qex except epsilon
import algorithms/numdiff
import helpers
import ../[core, scalar, gauge, functional]
import ../gauge/[types, field_ops, transport, cfield, stencil]

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))
qexInit()
let grt = initGraphRuntime()
include gauge/gaugehelpers
letParam:
  expectRanks = nRanks
check nRanks == expectRanks
echo "graph loop ranks: ", nRanks
# Three directions exercise every fundamental loop. The extent-two direction
# pins repeated physical links; the other extents keep SIMD outer sites even.
# MPI runs override geometry to keep each local outer extent even.
letParam:
  lat = latticeFromLocalLattice(@[2,4,4], nRanks)
let lo = lat.newLayout
let subDir = min(1, lat.len - 1)
const loopTestsOnly = true
var
  rng = lo.newRNGField(Philox4x64, 90911071u64)
  g = lo.newGauge
  u = lo.newGauge
  m = lo.newGauge
  q = lo.newGauge
threads:
  g.random rng
  u.random rng
  m.randomTAH rng
  q.randomTAH rng
  for f in m:
    f *= 0.1
  for f in q:
    f *= 0.1

include gauge/transport_basic
include gauge/higher
qexFinalize()
