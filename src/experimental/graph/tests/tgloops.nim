import math, unittest
import qex except epsilon
import base/alignedMem
import algorithms/numdiff
import helpers
import ../[core, scalar, gauge, functional]
import ../gauge/[types, field_ops, transport, cfield, stencil]

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))

proc runCase(f: proc()) {.noinline.} = f()

proc runLoopTests() =
  qexInit()
  defer: qexFinalize()
  let threshold = getRawMemGcThreshold()
  defer: setRawMemGcThreshold(threshold)
  # One coefficient check retains several GB, beyond the raw allocator's trigger.
  setRawMemGcThreshold(int.high)
  let grt = initGraphRuntime()
  template test(name, body: untyped) =
    block:
      # Keep each case's graph on a separate frame for collection.
      let fn = proc() =
        unittest.test name:
          body
      runCase(fn)
      grt.resetGradCache(false)
      grt.resetApplyCache(false)
      grt.resetLdjCache
      GC_fullCollect()
  include gauge/gaugehelpers
  echo "graph loop ranks: ", nRanks
  # Override lat with 2,4,8,8 to also exercise repeated physical links.
  letParam:
    lat = latticeFromLocalLattice(@[4,4,8,8], nRanks)
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

when isMainModule:
  runLoopTests()
