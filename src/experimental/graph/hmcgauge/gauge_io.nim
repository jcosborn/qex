import qex
from os import fileExists
from strformat import `&`
import config
from ../core/base import raiseValueError
from ../gauge/shared import reunitGauge, checkUnitary

proc loadOrInitGauge*(g: var auto, gaugefile: string) =
  if gaugefile.len == 0:
    g.unit
    return
  if not fileExists(gaugefile):
    raiseValueError("gauge file does not exist: " & gaugefile)
  tic("load")
  if 0 != g.loadGauge gaugefile:
    qexError "failed to load gauge file: ", gaugefile
  qexLog "loaded gauge from file: ", gaugefile, " secs: ", getElapsedTime()
  toc("read")
  block:
    tic("reunit")
    let before = g.checkUnitary
    echo "unitary deviation avg: ", before.avg, " max: ", before.max
    g.reunitGauge
    let after = g.checkUnitary
    echo "new unitary deviation avg: ", after.avg, " max: ", after.max
    toc("reunit")

proc maxGaugeDiff2*(a, b: auto): float =
  ## Max over links of the squared per-link difference |a−b|², reduced across
  ## threads/ranks. Used to report the load-time round-trip residual |f(f⁻¹(U))−U|
  ## for the field-transformation samplers.
  var m = 0.0
  threads:
    var mm = 0.0
    for mu in 0..<a.len:
      for x in a[mu]:
        let e = (a[mu][x] - b[mu][x]).norm2.simdMax
        if mm < e: mm = e
    mm.threadRankMax
    threadSingle: m = mm
  m

proc maybeSaveGauge*(g: auto,
                     config: RunConfig,
                     traj: int) =
  if config.savefreq <= 0 or traj mod config.savefreq != 0:
    return
  tic("save")
  let fn = config.savefile & &".{traj:05}.lime"
  if 0 != g.saveGauge(fn):
    qexError "Failed to save gauge to file: ", fn
  qexLog "saved gauge to file: ", fn, " secs: ", getElapsedTime()
  toc("done")
