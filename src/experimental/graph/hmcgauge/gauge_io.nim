import qex
from os import fileExists
from strformat import `&`
import config

proc reunitAndReportGauge(g: auto) =
  threads:
    let d = g.checkSU
    threadBarrier()
    echo "unitary deviation avg: ", d.avg, " max: ", d.max
    g.projectSU
    threadBarrier()
    let dd = g.checkSU
    echo "new unitary deviation avg: ", dd.avg, " max: ", dd.max

proc loadOrInitGauge*(g: var auto, gaugefile: string) =
  if fileExists(gaugefile):
    tic("load")
    if 0 != g.loadGauge gaugefile:
      qexError "failed to load gauge file: ", gaugefile
    qexLog "loaded gauge from file: ", gaugefile, " secs: ", getElapsedTime()
    toc("read")
    block:
      tic("reunit")
      g.reunitAndReportGauge
      toc("reunit")
  else:
    g.unit

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
