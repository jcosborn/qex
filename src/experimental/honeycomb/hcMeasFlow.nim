## hcMeasFlow -- gradient flow + E + topological charge on the 16-cell
## honeycomb (task **W3**; pattern of src/examples/wflow_topo.nim).
##
## Loads a configuration written by task M's hcio (`-gaugefile:...lime`), or
## generates a warm one in place (`-geom -warms -seed`) so the tool is usable
## without an ensemble.  Flows with the calibrated triangle-action flow
## (hcflow; `t` is continuum flow time in units a^2, directly comparable to
## cubic flow time) and prints per measurement
##
##   HCFLOW  t   E   t^2E   t d/dt(t^2E)   Q
##
## with E = <E(x)> from the hexagon clover (intensive; no a^4/2) and
## Q = (1/2) sum_x q(x) (hctopo).  Stops on tmax / t2Emax / tdt2Emax.  At the
## end interpolates t0 (t^2E = t0target) and w0^2 (t d/dt(t^2E) = w0target)
## and prints `T0 <t0> <Q(t0)>` and `W0SQ <w0^2>`.
##
##   export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
##   cd build_mac && make src/experimental/honeycomb/hcMeasFlow.nim
##   OMP_NUM_THREADS=4 ./bin/hcMeasFlow -geom:8,8,8,8 -warms:0.35 -tmax:2
##   OMP_NUM_THREADS=4 ./bin/hcMeasFlow -gaugefile:cfg.00100.lime -eps:0.02
##
## -simdlen:1 for cell geometries the VLEN=4 layout rejects (6^4, 9^4, ...).
##
import std/[math, os, monotimes, times]
import qex
import physics/qcdTypes
import hcgeom, hclayout, hcgauge, hcaction, hcflow, hctopo, hcio
import hcanalysis

qexInit()
tic()

letParam:
  gaugefile = ""         ## honeycomb .lime configuration (hcio); "" = warm start
  geom = @[8, 8, 8, 8]   ## cell geometry for the warm start
  warms = 0.35           ## warm(s) strength for the generated configuration
  seed = 987654321       ## RNG seed for the warm start
  eps = 0.02             ## RK3 step, continuum flow-time units (a^2)
  cflow = 6.0            ## flow normalisation; hcFlowCflow = 6 is calibrated
  tmax = 2.0             ## stop at this flow time (0 = no limit)
  t2Emax = 0.45          ## stop once t^2 E exceeds this (0 = no limit)
  tdt2Emax = 0.35        ## stop once t d/dt(t^2E) exceeds this (0 = no limit)
  measevery = 1          ## measure E and Q every this many RK steps
  t0target = 0.3         ## t^2 E = t0target defines t0
  w0target = 0.3         ## t d/dt(t^2E) = w0target defines w0^2
  reunitevery = 0        ## reunitarise every N steps (0 = never; drift ~1e-15)
  simdlen = 0            ## 0 = default VLEN; 1 (or 2) for odd geometries
  showTimers: bool = 0

installHelpParam()
echoParams()
echo "rank ", myRank, "/", nRanks
threads: echo "thread ", threadNum, "/", numThreads

if cflow != hcFlowCflow:
  qexWarn "cflow = ", cflow, " != calibrated hcFlowCflow = ", hcFlowCflow,
          ": t is NOT continuum-normalised flow time"

let haveFile = gaugefile.len > 0 and fileExists(gaugefile)
if gaugefile.len > 0 and not haveFile:
  qexError "gauge file not found: ", gaugefile
let cgeom = if haveFile: hcFileGeom(gaugefile) else: geom
echo "cell geometry: ", cgeom, "  (sites ", 2*cgeom[0]*cgeom[1]*cgeom[2]*cgeom[3], ")"

template runAll(hlx: typed) =
  let hl = hlx
  var g = newHcGauge(hl)
  if haveFile:
    let (st, meta) = loadHcGauge(g, gaugefile)
    if st != 0:
      qexError "failed to load ", gaugefile, " (status ", st, ")"
    echo "loaded ", gaugefile, "  beta ", meta.beta, "  traj ", meta.traj,
         "  info '", meta.info, "'"
    threads:
      g.reunit
  else:
    var r = hl.lo.newRNGField(RngMilc6, uint64 seed)
    threads:
      g.warm(warms, r)
    echo "generated warm(", warms, ") configuration, seed ", seed
  block:
    let su = g.checkSU
    echo "checkSU: avg ", su.avg, "  max ", su.max
  toc("setup")

  var wt = newHcTopoWork(g)
  var
    ts = @[0.0]
    t2Es = @[0.0]
    qs: seq[float]
  block:
    let (e0, q0) = hcEQ(wt, g)
    qs.add q0
    echo "# HCFLOW  t  E  t^2E  t*d/dt(t^2E)  Q"
    echo "HCFLOW 0.0 ", e0, " 0.0 0.0 ", q0
  toc("t=0 measurement")

  var nstep = 0
  var nmeas = 0
  let tstart = getMonoTime()
  g.hcGaugeFlow(eps, cflow):
    inc nstep
    if reunitevery > 0 and nstep mod reunitevery == 0:
      threads:
        g.reunit
    if nstep mod measevery == 0:
      inc nmeas
      let (e, q) = hcEQ(wt, g)
      let
        t2E = wflowT*wflowT*e
        # Endpoint estimate for live output and the stopping condition.
        dt2E = (t2E - t2Es[^1])/(eps*measevery.float)
        tdt2E = wflowT*dt2E
      ts.add wflowT
      t2Es.add t2E
      qs.add q
      echo "HCFLOW ", wflowT, " ", e, " ", t2E, " ", tdt2E, " ", q
      if (tmax > 0 and wflowT >= tmax - 1e-9) or
         (t2Emax > 0 and t2E > t2Emax) or
         (tdt2Emax > 0 and tdt2E > tdt2Emax):
        break
  let secs = (getMonoTime()-tstart).inMicroseconds.float*1e-6
  toc("flow")

  let t0 = findT0(ts, t2Es, t0target)
  let w0sq = findW0(ts, derivT2E(ts, t2Es), w0target)
  if t0 > 0.0:
    let q0 = interpAt(ts, qs, t0)
    echo "T0 ", t0, " ", q0, "   # t0/a^2, Q(t0)  [t^2E = ", t0target, "]"
  else:
    echo "T0 -1 0   # t^2E never reached ", t0target,
         " (volume too small or flow too short, cf. STATUS: L/sqrt(t0) >~ 9)"
  if w0sq > 0.0:
    echo "W0SQ ", w0sq, "   # w0^2/a^2  [t d/dt(t^2E) = ", w0target, "]"
  echo "FLOWTIME ", secs, " s for ", nstep, " RK steps + ", nmeas,
       " (E,Q) measurements  (", 1e3*secs/nstep.float, " ms/step incl. meas)"

case simdlen
of 0: runAll(newHcLayout(cgeom))
of 1: runAll(newHcLayoutX(cgeom, 1))
of 2: runAll(newHcLayoutX(cgeom, 2))
else: qexError "unsupported -simdlen:", simdlen, " (0 = default, 1, 2)"

if showTimers: echoTimers()
qexFinalize()
