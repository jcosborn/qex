import qex, gauge/wflow
import os

proc reunit(g:auto) =
  tic()
  threads:
    let d = g.checkSU
    threadBarrier()
    echo "unitary deviation avg: ",d.avg," max: ",d.max
    g.projectSU
    threadBarrier()
    let dd = g.checkSU
    echo "new unitary deviation avg: ",dd.avg," max: ",dd.max
  toc("reunit")

proc EQ(g:auto,loop:int):auto =
  let
    f = g.fmunu loop
    (es,et) = f.densityE
    q = f.topoQ
  return (es,et,q)

qexinit()

tic()

letParam:
  gaugefile = ""
  lat =
    if fileExists(gaugefile):
      getFileLattice gaugefile
    else:
      if gaugefile.len > 0:
        qexWarn "Nonexistent gauge file: ", gaugefile
      @[8,8,8,8]
  dt = 0.02
  tmax = 2.0
  t2Emax = 0.45
  tdt2Emax = 0.35
  fmunuloop = 5
  showTimers:bool = 0

installHelpParam()
echoParams()
echo "rank ", myRank, "/", nRanks
threads: echo "thread ", threadNum, "/", numThreads

let lo = lat.newLayout
var g = lo.newgauge

if fileExists(gaugefile):
  tic("load")
  if 0 != g.loadGauge gaugefile:
    qexError "failed to load gauge file: ", gaugefile
  qexLog "loaded gauge from file: ", gaugefile," secs: ",getElapsedTime()
  toc("read")
else:
  var r = lo.newRNGField(RngMilc6)
  g.warm(0.5,r)
  toc("random warm gauge")

g.echoPlaq
toc("plaq")
g.reunit
toc("reunit")
g.echoPlaq
toc("plaq")

let (es,et,q) = g.EQ fmunuloop
var t2E = 0.0
var ot2E,dt2E,tdt2E:float

echo "# WFLOW\tt\tEt\tEs\tQ\tt^2E\ttd/dt[t^2E]"
echo "WFLOW 0.0 ",et," ",es," ",q

g.gaugeFlow(dt):
  let (es,et,q) = g.EQ fmunuloop
  ot2E = t2E
  t2E = wflowT*wflowT*(es+et)
  dt2E = (t2E-ot2E)/dt
  tdt2E = wflowT*dt2E
  echo "WFLOW ",wflowT," ",et," ",es," ",q," ",t2E," ",tdt2E
  if (tmax>0 and wflowT>tmax) or (t2Emax>0 and t2E>t2Emax) or (tdt2Emax>0 and tdt2E>tdt2Emax):
    break

if showTimers: echoTimers()
qexFinalize()
