#RUNCMD env OMP_NUM_THREADS=2 $RUN1

import math, strutils, unittest

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))

# basicOps.epsilon collides with fenv.epsilon
import qex except epsilon
import algorithms/numdiff, gauge/stoutsmear
import helpers
import ../[core, scalar, multi, gauge]
import ../functional
import ../gauge/types as graphGaugeShared
import ../gauge/basic_ops as graphGaugeBasic
import ../gauge/[types, field_ops, transport, cfield]
import ../hmcgauge/optimizer, ../hmcgauge/integrator
import ../hmcgauge/trajectory
import ../hmcgauge/training
import ../hmcgauge/config
import ../hmcgauge/gauge_io
import ../hmcgauge/params
import ../hmcgauge/rng

let grt = initGraphRuntime()

include gauge/gaugehelpers

qexInit()

let
  lat = @[8,8,8,16]
  lo = lat.newLayout
  seed = 1234567891u64
  vol = lo.physVol
var
  r = lo.newRNGField(Philox4x64, seed)
  g = lo.newgauge
  u = lo.newgauge
  p = lo.newgauge
  q = lo.newgauge
  m = lo.newgauge
  ss = lo.newStoutSmear(0.1)
const nc = g[0][0].nrows
threads:
  g.random r
  u.random r
  p.randomTAH r
  q.randomTAH r
  m.randomTAH r
for i in 0..4:
  ss.smear(g, g)
  ss.smear(u, u)
threads:
  for t in m:
    t *= 0.01

let scalarValues = sampleScalarValues()
let a = scalarValues.a
let b = scalarValues.b

proc zeroGaugeLike(source: graphGaugeShared.Gauge): graphGaugeShared.Gauge =
  result = source.newOneOf
  graphGaugeShared.zeroGaugeStorage(result)

include gauge/coeffs
include gauge/basic
include gauge/field
include gauge/transport
include gauge/fused_basic
include gauge/fused
include gauge/action
include hmcgauge/basic

qexFinalize()
