## measflow -- gradient flow, t0, w0 and topological charge on saved
## configurations, -lattice:hc (hcpuregauge files) or -lattice:cubic
## (refcubicgen files), with the ensemble analysis of the topological
## susceptibility.
##
## Per configuration the flow history (t, t^2E, t d/dt(t^2E), Q) is recorded
## until tmax or t^2E > t2Estop; t0 (t^2E = t0target), w0^2 (t d/dt(t^2E) =
## w0target) and Q(t0) are interpolated and printed as
##   CFG name t0 w0sq Q(t0) Q(tend) tend roundDist
## Ensemble: jackknife (bin -bin) of t0, <Q^2> and 10^4 t0^2 chi with
## chi = <Q^2>/V, V = number of cells (honeycomb, a^4/2 per site) or sites
## (cubic); printed as an ENSEMBLE line.  E and Q: hexagon clover (hctopo) on
## the honeycomb, 1x1 clover fmunu/densityE/topoQ on the cubic lattice; t in
## units of a^2 on both.
##
## -gaugefile: a single configuration instead of -cfgs; -warms: a warm
## honeycomb start without any file; -flowdump: write the first history;
## -chifit: only fit a table of "a^2/t0  10^4 t0^2 chi  err".

import std/[math, os, algorithm, strformat, strutils, monotimes, times]
import qex
import physics/qcdTypes
import gauge, gauge/wflow
import honeycomb, hcanalysis

qexInit()
tic()

letParam:
  lattice = "hc"         ## hc | cubic
  cfgs = ""              ## glob pattern of configuration files
  gaugefile = ""         ## a single configuration file
  geom = @[8, 8, 8, 8]   ## honeycomb cell geometry for a warm start (no file)
  warms = 0.0            ## warm(s) strength of the generated start (0 = none)
  seed = 987654321       ## RNG seed of the warm start
  eps = 0.05             ## RK3 step in continuum flow time (a^2)
  tmax = 0.0             ## stop at this flow time (0 = t2Estop only)
  t2Estop = 0.45         ## stop once t^2E exceeds this (0 = no limit)
  measevery = 2          ## measure E and Q every this many RK steps
  fmunuloop = 1          ## cubic clover size for fmunu: 1, 3, 4 or 5
  t0target = 0.3
  w0target = 0.3
  interpOrder = 1        ## 1 = linear, 3 = cubic interpolation of the crossings
  bin = 1                ## jackknife block size (saved configurations)
  beta = 0.0             ## label only
  flowdump = ""          ## write the flow history of the first configuration
  outfile = ""           ## CFG lines here as well as to stdout
  chifit = ""            ## fit this table instead of measuring
  simdlen = 0            ## honeycomb: 0 = default VLEN, 1 or 2 for odd geometries
  showTimers: bool = 0

installHelpParam()
echoParams()

if chifit.len > 0:
  var x, y, dy: seq[float]
  for r in readColumns chifit:
    if r.len >= 3:
      x.add r[0]
      y.add r[1]
      dy.add r[2]
  echo &"{x.len} points, x in [{x.min:.4f}, {x.max:.4f}]"
  for (name, powers) in [("c0 + c2 x", @[0, 1]), ("c0 + c4 x^2", @[0, 2]),
                         ("c0 + c2 x + c4 x^2", @[0, 1, 2])]:
    if x.len <= powers.len: continue
    let (co, er, cd) = fitPoly(x, y, dy, powers)
    var s = ""
    for k in 0..<co.len: s &= &"  c{2*powers[k]} = {co[k]:8.4f} +- {er[k]:.4f}"
    echo &"FIT {name}:{s}  chi^2/dof = {cd:.3f}"
  qexFinalize()
  quit(0)

var files: seq[string]
if gaugefile.len > 0:
  files.add gaugefile
else:
  for f in walkPattern(cfgs):
    files.add f
  files.sort
let warmStart = files.len == 0 and lattice == "hc" and warms > 0.0
if files.len == 0 and not warmStart:
  qexError "no configuration files matched '", cfgs, "'"
let cgeom = if files.len > 0: getFileLattice files[0] else: geom
echo "lattice ", lattice, " geometry ", cgeom, " files ", files.len
threads: echo "thread ", threadNum, "/", numThreads

var fr: FlowRec
var
  t0s, w0sqs, qs, q2s, qends: seq[float]
  outLines: seq[string]

proc finishConfig(name: string, first: bool) =
  ## crossings of the current history, CFG line, optional dump
  let
    w = derivT2E(fr.ts, fr.t2Es)
    t0 = findT0(fr.ts, fr.t2Es, t0target, interpOrder)
    w0sq = findW0(fr.ts, w, w0target, interpOrder)
  if flowdump.len > 0 and first and myRank == 0:
    var fh = open(flowdump, fmWrite)
    fh.write("# t  t^2E  t d/dt(t^2E)  Q\n")
    for i in 0..<fr.ts.len:
      fh.write(&"{fr.ts[i]:.6f} {fr.t2Es[i]:.8f} {w[i]:.8f} {fr.qs[i]:.8f}\n")
    fh.close
  if t0 <= 0.0:
    echo "WARNING: t^2E never reached ", t0target, " on ", name, ", skipped"
    return
  let
    q0 = interpAt(fr.ts, fr.qs, t0)
    qend = fr.qs[^1]
    tend = fr.ts[^1]
  t0s.add t0
  w0sqs.add w0sq
  qs.add q0
  q2s.add q0*q0
  qends.add qend
  let line = &"CFG {name} {t0:.6f} {w0sq:.6f} {q0:.6f} {qend:.6f} {tend:.4f} {abs(qend - round(qend)):.6f}"
  echo line
  outLines.add line

proc record(t, e, q: float): bool =
  ## add a flow point; true once the stop condition is met
  fr.add(t, e, q)
  (tmax > 0.0 and t >= tmax - 1e-9) or (t2Estop > 0.0 and t*t*e > t2Estop)

proc runHc(lo: Layout) =
  var g = newHcGauge(lo)
  var wt = newTopoWork(g)
  proc eq(): tuple[e, q: float] = EQ(wt, g)
  toc("setup")
  proc history() =
    fr.start(eq().q)
    var nstep = 0
    g.flow(0, eps, cflow):
      inc nstep
      if nstep mod measevery == 0:
        let (e, q) = eq()
        if record(wflowT, e, q): break
  if warmStart:
    var r = lo.newRNGField(RngMilc6, uint64 seed)
    threads:
      g.warm(warms, r)
    history()
    finishConfig(&"warm{warms}", true)
  for i, fn in files:
    let (st, meta) = load(g, fn)
    if st != 0: qexError "failed to load ", fn
    threads:
      g.reunit
    echo &"# {fn.extractFilename}: beta {meta.beta} traj {meta.traj} triangleSum {g.triangleSum:.8f}"
    history()
    finishConfig(fn.extractFilename, i == 0)

proc runCubic(lo: Layout) =
  var g = lo.newGauge
  proc eq(): tuple[e, q: float] =
    let f = g.fmunu fmunuloop
    let (es, et) = f.densityE
    (es + et, f.topoQ)
  toc("setup")
  for i, fn in files:
    if 0 != g.loadGauge fn: qexError "failed to load ", fn
    threads:
      g.projectSU
    fr.start(eq().q)
    var nstep = 0
    g.gaugeFlow(0, eps):
      inc nstep
      if nstep mod measevery == 0:
        let (e, q) = eq()
        if record(wflowT, e, q): break
    finishConfig(fn.extractFilename, i == 0)

echo "# CFG name t0 w0sq Q(t0) Q(tend) tend roundDist"
case lattice
of "hc":
  withLayout(simdlen, cgeom, lo):
    runHc(lo)
of "cubic":
  runCubic(cgeom.newLayout)
else:
  qexError "unknown -lattice:", lattice, " (hc | cubic)"
toc("measure")

if outfile.len > 0 and myRank == 0:
  var fh = open(outfile, fmWrite)
  fh.write("# CFG name t0 w0sq Q(t0) Q(tend) tend roundDist\n")
  for l in outLines: fh.write(l & "\n")
  fh.close

let n = t0s.len
if n == 0:
  qexError "no usable measurements"
var vol = 1.0
for x in cgeom: vol *= x.float
let
  (t0m, t0e) = jackknifeMean(t0s, bin)
  (q2m, q2e) = jackknifeMean(q2s, bin)
  (chiM, chiE, xM, xE) = chiTop(t0s, q2s, vol, bin)
  tauT0 = autocorrTime(t0s)
  tauQ = autocorrTime(qs)
  tauQ2 = autocorrTime(q2s)
var rdAvg, rdMax = 0.0
for q in qends:
  let d = abs(q - round(q))
  rdAvg += d/qends.len.float
  rdMax = max(rdMax, d)
echo ""
echo &"configs        {n}   (jackknife bin = {bin}; V = {vol})"
echo &"t0/a^2         {t0m:.5f} +- {t0e:.5f}"
echo &"w0^2/a^2       {w0sqs.mean:.5f} +- {jackknifeMean(w0sqs, bin).err:.5f}"
echo &"a^2/t0         {xM:.5f} +- {xE:.5f}"
echo &"<Q>            {qs.mean:.4f} +- {jackknifeMean(qs, bin).err:.4f}"
echo &"<Q^2>(t0)      {q2m:.4f} +- {q2e:.4f}"
echo &"10^4 t0^2 chi  {chiM:.4f} +- {chiE:.4f}"
echo &"tau_int        t0 {tauT0:.2f}   Q {tauQ:.2f}   Q^2 {tauQ2:.2f}   (in saved configurations)"
echo &"L/sqrt(t0)     {cgeom[0].float/sqrt(t0m):.2f}"
echo &"Q at t_end     mean dist to nearest integer {rdAvg:.4f}, max {rdMax:.4f}"
echo &"ENSEMBLE {beta} {cgeom.join(\"x\")} {n} {t0m:.6f} {t0e:.6f} {xM:.6f} {xE:.6f} {q2m:.6f} {q2e:.6f} {chiM:.6f} {chiE:.6f} {tauQ:.3f} {tauQ2:.3f}"

if showTimers: echoTimers()
qexFinalize()
