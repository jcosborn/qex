## cubicSpectrum -- the cubic-lattice leg of task **D4**: low-lying clover
## Wilson-Dirac spectrum, chirality, Q_Dirac and Q_flow on saved cubic SU(3)
## configurations (refCubicGen ensembles), mirroring hcSpectrum.nim line for
## line so the two lattices can be compared file against file.
##
## Per configuration (loaded from -cfgs glob):
##   a. gradient flow on the PERIODIC, UNsmeared links (QEX gaugeFlow +
##      fmunu/densityE/topoQ, task C conventions):
##      Q(t) and t^2E(t) on a fixed grid to -flowtmax;
##      Q_flow = Q(t = -t0use), plus the per-config/ensemble t0 crossings;
##   b. 6 stout steps rho = 0.05 (QEX StoutSmear, the plain MP step),
##      antiperiodic time BC, massless r=1 cSW=1 clover
##      operator (task D2 CubicCloverWilson), shift-invert Arnoldi
##      (hcSpectrum machinery: CGNR inner solves, direct residuals with the
##      exact operator, chirality, real-mode classification, Q_Dirac).
##
## Output lines: identical to hcSpectrum (EIG / FLOWQ / CFG / CFGX / ENST2E /
## ENST0), so doc/plots/fermions/harvest.py consumes both.
##
##   export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
##   cd build_mac && make src/experimental/honeycomb/cubicSpectrum.nim
##   OMP_NUM_THREADS=4 ./bin/cubicSpectrum -cfgs:'/tmp/cub586/cfg.*.lime' \
##     -nev:32 -sigma:-0.25 -outfile:cub586.log

import std/[math, complex, os, algorithm, strformat, strutils, monotimes, times]
import qex except epsilon
import physics/qcdTypes
import gauge, gauge/wflow
import hcanalysis
import hcSpectrum

qexInit()
tic()

letParam:
  cfgs = ""                ## glob pattern for refCubicGen .lime configurations
  beta = 0.0               ## label only (recorded in the output header)
  maxcfg = 0               ## measure at most this many (0 = all)
  firstcfg = 1             ## skip files before this index (restart aid)
  nstout = 6               ## stout steps
  rho = 0.05               ## stout rho (paper value; MP alpha)
  mass = 0.0               ## valence quark mass
  cSW = 1.0                ## clover coefficient
  sigma = -0.25            ## shift-invert point (real, left of the spectrum)
  nev = 32
  ncv = 0                  ## 0 = 3*nev
  tol = 1e-7
  maxRestarts = 60
  residcut = 2e-5          ## exclude modes with direct residual above this
  innerR2 = 1e-12
  innerMaxIts = 4000
  epsreal = 1e-6
  doflow: bool = 1
  doeigs: bool = 1
  floweps = 0.05
  flowmeas = 2
  flowtmax = 2.6
  t0use = 1.917            ## Q_flow = Q(t0use); nominal matched t0/a^2
  verb = 0
  outfile = ""
  showTimers: bool = 0

installHelpParam()
echoParams()
echo "rank ", myRank, "/", nRanks
threads: echo "thread ", threadNum, "/", numThreads

var files: seq[string] = @[]
for f in walkPattern(cfgs):
  files.add f
files.sort
if files.len == 0:
  qexError "no configuration files matched: '", cfgs, "'"
echo "found ", files.len, " configurations"

let ncvUse = if ncv > 0: ncv else: 3*nev
var fh: File
let haveOut = outfile.len > 0 and myRank == 0
if haveOut: fh = open(outfile, fmWrite)
proc emit(line: string) =
  echo line
  if haveOut:
    fh.writeLine line
    fh.flushFile

let
  clat = getFileLattice files[0]
  lo = clat.newLayout

emit &"# cubicSpectrum lat {clat} beta {beta} nfiles {files.len}"
emit &"# nstout {nstout} rho {rho} cSW {cSW} mass {mass}"
emit &"# sigma {sigma} nev {nev} ncv {ncvUse} tol {tol} innerR2 {innerR2} epsreal {epsreal}"
emit &"# floweps {floweps} flowmeas {flowmeas} flowtmax {flowtmax} t0use {t0use}"
emit &"# qDiracSign {qDiracSign}"

var g = lo.newGauge          # loaded configuration; the flow runs on it in place
var gferm = lo.newGauge      # smeared + BC copy (the operator lives on this)
var st = newStoutSmear(lo, rho)
var cw = newCubicCloverWilson(gferm, cSW)

var proto = lo.DiracFermion()
type DF = typeof(proto)

var startCount = 0'u64
proc startVec(v: var DF) =
  inc startCount
  let salt = sm64(startCount)*0x10000'u64
  for i in lo.sites:
    for sp in 0..3:
      for c in 0..2:
        let k = salt + uint64(i)*24 + uint64(sp)*6 + uint64(c)*2
        setC(v{i}[sp][c], u01(k), u01(k+1))

var stats = SiStats()
let msolve = mass - sigma
var op = newShiftInvertOp[DF](
  applyM = proc (rr: var DF; x: DF) = cw.D(rr, x, msolve),
  applyMdag = proc (rr: var DF; x: DF) = cw.Ddag(rr, x, msolve),
  newVec = proc (): DF =
    result = newOneOf(proto)
    result := 0,
  startVec = startVec,
  r2req = innerR2, maxits = innerMaxIts, stats = stats)

proc eqm(gg: auto): tuple[e, q: float] =
  ## (E, Q) from the 1x1 clover field strength (task C conventions)
  let f = gg.fmunu 1
  let (es, et) = f.densityE
  (es + et, f.topoQ)

# flow bookkeeping
let nFlowSteps = int(round(flowtmax/floweps))
var flowT: seq[float]
var sumT2E: seq[float]
var nFlowCfg = 0
var lastQflow = 0.0

proc measureFlow(icfg: int) =
  var ts = @[0.0]
  var t2Es = @[0.0]
  var qs: seq[float]
  block:
    let (e0, q0) = eqm(g)
    discard e0
    qs.add q0
  var nstepF = 0
  g.gaugeFlow(nFlowSteps, floweps):
    inc nstepF
    if nstepF mod flowmeas == 0:
      let (e, q) = eqm(g)
      ts.add wflowT
      t2Es.add wflowT*wflowT*e
      qs.add q
  if flowT.len == 0:
    flowT = ts
    sumT2E = newSeq[float](ts.len)
  inc nFlowCfg
  for j in 0..<min(t2Es.len, sumT2E.len):
    sumT2E[j] += t2Es[j]
  let t0cfg = findT0(ts, t2Es, 0.3, 1)
  lastQflow = interpAt(ts, qs, t0use)
  let q15 = interpAt(ts, qs, 1.5)
  let q24 = interpAt(ts, qs, min(2.4, ts[^1]))
  emit &"FLOWQ {icfg} {t0cfg:.6f} {lastQflow:.6f} {q15:.6f} {q24:.6f} {qs[0]:.6f}"

proc measureSpec(icfg: int) =
  let tw0 = getMonoTime()
  threads:
    for mu in 0..<gferm.len:
      gferm[mu] := g[mu]
  st.smearN(gferm, gferm, nstout)
  gferm.setBC
  cw.gaugeRefresh
  stats.reset
  let (mus, vecs, _, napply) =
    arnoldi(op, nev, ncvUse, tol, maxRestarts, "LM", verb)
  let modes = measureModes(mus, vecs, sigma,
    proc(r: var DF; x: DF) = cw.D(r, x, mass), residcut)
  let secs = (getMonoTime() - tw0).inMicroseconds.float*1e-6
  for line in specLines(modes, epsreal, icfg, rho, lastQflow,
                        mus.len - modes.len, napply, stats, secs):
    emit line

toc("setup")
let trun0 = getMonoTime()
var nmeas = 0
for icfg0, fn in files:
  let icfg = icfg0 + 1
  if icfg < firstcfg: continue
  if maxcfg > 0 and nmeas >= maxcfg: break
  inc nmeas
  if 0 != g.loadGauge fn:
    qexError "failed to load ", fn
  threads:
    g.projectSU
  block:
    let pl = g.plaq
    var s = 0.0
    for v in pl: s += v
    emit &"TS {icfg} {s:.8f}   # plaq, {fn.extractFilename}"
  lastQflow = 0.0
  # order matters: copy for the fermion leg BEFORE the flow destroys g
  threads:
    for mu in 0..<gferm.len:
      gferm[mu] := g[mu]
  if doflow: measureFlow(icfg)
  if doeigs:
    threads:                     # measureSpec re-copies from g; restore it
      for mu in 0..<g.len:
        g[mu] := gferm[mu]
    measureSpec(icfg)
  GC_fullCollect()
  let el = (getMonoTime() - trun0).inMicroseconds.float*1e-6
  emit &"# done cfg {icfg} ({nmeas} measured), elapsed {el:.1f} s ({el/float(nmeas):.1f} s/cfg)"

if nFlowCfg > 0:
  for j in 0..<flowT.len:
    emit &"ENST2E {flowT[j]:.4f} {sumT2E[j]/float(nFlowCfg):.8f}"
  var meanT2E = newSeq[float](flowT.len)
  for j in 0..<flowT.len: meanT2E[j] = sumT2E[j]/float(nFlowCfg)
  let t0ens = findT0(flowT, meanT2E, 0.3, 1)
  emit &"ENST0 {t0ens:.6f}   # ensemble <t^2E> = 0.3 crossing; " &
       "finite-volume biased at L/sqrt(t0) < 9"
toc("run")
if haveOut: fh.close
if showTimers: echoTimers()
qexFinalize()
