## spectrum -- low-lying clover Wilson-Dirac spectrum, chirality, Q_Dirac and
## Q_flow on saved configurations: -lattice:hc (files of hcpuregauge) or
## -lattice:cubic (files of refcubicgen), one output format for both.
##
## Per configuration:
##   flow on the periodic, unsmeared links: Q(t), t^2E(t) on a fixed grid to
##   -flowtmax; Q_flow at the per-config t0 (t^2E = 0.3) and at -t0use;
##   nstout stout steps at rho (plain Morningstar-Peardon convention on both
##   lattices; the honeycomb step is rho/3 in flow time, hcflow.stoutKappa),
##   antiperiodic time, clover operator (cSW, mass, r), shift-invert Arnoldi
##   for the nev eigenvalues nearest sigma; per mode lambda, chirality and
##   direct residual; real modes counted for Q_Dirac only below -recut and
##   the reach of the converged window is reported (hcspec).
## Lines (also copied to -outfile):
##   TS icfg loopavg / FLOWQ icfg t0cfg q(t0use) q(t0cfg) q(1.5) q(2.4) q(0) /
##   EIG icfg rho k re im chi resid /
##   CFG icfg rho lam0re lam0im nconv nreal nplus nminus qdirac q(t0use) q(t0cfg) napply cgits secs /
##   CFGX icfg rho reach recut realabove conj worstresid imgapC imgapR sumchi nbad cgmax cghitmax sidirect /
##   ENST2E t <t2E> / ENST0 t0ens

import std/[math, complex, os, algorithm, strformat, strutils, monotimes, times]
import qex except epsilon
import physics/qcdTypes
import gauge, gauge/wflow
import honeycomb, cubic, hcspec, hcanalysis

qexInit()
tic()

letParam:
  lattice = "hc"           ## hc | cubic
  cfgs = ""                ## glob pattern of configuration files
  maxcfg = 0               ## measure at most this many (0 = all)
  firstcfg = 1             ## skip files before this index
  nstout = 6
  rho = 0.05
  mass = 0.0
  rw = 1.0                 ## Wilson r (honeycomb; QEX wilsonD has r = 1)
  cSW = 1.0
  sigma = -0.25            ## shift-invert point, left of the spectrum
  nev = 32
  ncv = 0                  ## Krylov size (0 = 3 nev)
  tol = 1e-7
  maxRestarts = 60
  residcut = 2e-5          ## drop modes with direct residual above this
  innerR2 = 1e-12          ## CGNR |r|^2/|rhs|^2 stop
  innerMaxIts = 4000
  epsreal = 1e-5           ## |Im lambda| below this = real mode
  recut = 0.15             ## real modes counted for Q_Dirac if Re lambda < recut
  doflow: bool = 1
  doeigs: bool = 1
  floweps = 0.05
  flowmeas = 2             ## measure (E,Q) every this many RK steps
  flowtmax = 2.6
  t0use = 1.917            ## nominal matched t0/a^2 for the second Q_flow
  verb = 0
  outfile = ""
  simdlen = 0              ## honeycomb: 0 = default VLEN, 1 or 2 for odd geometries
  showTimers: bool = 0

installHelpParam()
echoParams()
echo "rank ", myRank, "/", nRanks
threads: echo "thread ", threadNum, "/", numThreads

var files: seq[string]
for f in walkPattern(cfgs):
  files.add f
files.sort
if files.len == 0:
  qexError "no configuration files matched '", cfgs, "'"
let cgeom = getFileLattice files[0]
echo "found ", files.len, " configurations, geometry ", cgeom

let ncvUse = if ncv > 0: ncv else: 3*nev
var fh: File
let haveOut = outfile.len > 0 and myRank == 0
if haveOut: fh = open(outfile, fmWrite)
proc emit(line: string) =
  echo line
  if haveOut:
    fh.writeLine line
    fh.flushFile

emit &"# spectrum {lattice} geom {cgeom} nfiles {files.len}"
emit &"# nstout {nstout} rho {rho} cSW {cSW} mass {mass} rw {rw}"
emit &"# sigma {sigma} nev {nev} ncv {ncvUse} tol {tol} innerR2 {innerR2} epsreal {epsreal} recut {recut}"
emit &"# floweps {floweps} flowmeas {flowmeas} flowtmax {flowtmax} t0use {t0use}"
emit &"# qDiracSign {qDiracSign}"

let nFlowSteps = int(round(flowtmax/floweps))
let msolve = mass - sigma
var fr: FlowRec

proc flowResult(icfg: int): tuple[qflow, qflowT0: float] =
  ## FLOWQ line and the ensemble accumulation for the current history
  fr.accumulate
  let t0cfg = findT0(fr.ts, fr.t2Es, 0.3, 1)
  let qflow = interpAt(fr.ts, fr.qs, t0use)
  let qt0 = if t0cfg > 0.0: interpAt(fr.ts, fr.qs, t0cfg) else: 0.0
  let q15 = interpAt(fr.ts, fr.qs, 1.5)
  let q24 = interpAt(fr.ts, fr.qs, min(2.4, fr.ts[^1]))
  emit &"FLOWQ {icfg} {t0cfg:.6f} {qflow:.6f} {qt0:.6f} {q15:.6f} {q24:.6f} {fr.qs[0]:.6f}"
  (qflow, qt0)

proc eigsConfig[F](op: ArnoldiOp[F]; applyD: proc(r: var F; x: F);
                   stats: SiStats; icfg: int; qflow, qflowT0: float;
                   tw0: MonoTime) =
  stats.reset
  let (mus, vecs, _, napply) =
    arnoldi(op, nev, ncvUse, tol, maxRestarts, "LM", verb)
  let modes = measureModes(mus, vecs, sigma, applyD, residcut)
  let secs = (getMonoTime() - tw0).inMicroseconds.float*1e-6
  for line in specLines(modes, epsreal, sigma, recut, icfg, rho, qflow,
                        qflowT0, mus.len - modes.len, napply, stats, secs):
    emit line

proc finish() =
  if fr.ncfg > 0:
    let m = fr.meanT2E
    for j in 0..<fr.gridT.len:
      emit &"ENST2E {fr.gridT[j]:.4f} {m[j]:.8f}"
    emit &"ENST0 {findT0(fr.gridT, m, 0.3, 1):.6f}"
  toc("run")
  if haveOut: fh.close
  if showTimers: echoTimers()

template configLoop(body: untyped) =
  ## `icfg` and `fn` injected; runs body for the selected files
  let trun0 = getMonoTime()
  var nmeas {.inject.} = 0
  for icfg0, fn {.inject.} in files:
    let icfg {.inject.} = icfg0 + 1
    if icfg < firstcfg: continue
    if maxcfg > 0 and nmeas >= maxcfg: break
    inc nmeas
    body
    GC_fullCollect()
    let el {.inject.} = (getMonoTime() - trun0).inMicroseconds.float*1e-6
    emit &"# done cfg {icfg} ({nmeas} measured), elapsed {el:.1f} s ({el/float(nmeas):.1f} s/cfg)"

proc runHc(lo: Layout) =
  var g = newHcGauge(lo)         # loaded configuration, periodic
  var gf = newOneOf(g)           # flow copy
  var gs = newOneOf(g)           # smeared + BC copy, the operator lives here
  var sc = newOneOf(g)           # stout scratch
  var w = newActionWork(g)
  var wt = newTopoWork(gf)
  var cw = newHcWilson(gs, cSW)
  var proto = newHcFermion(lo)
  type HF = typeof(proto)
  var startCount = 0'u64
  proc startVec(v: var HF) =
    inc startCount
    let salt = sm64(startCount)*0x10000'u64
    for i in lo.sites:
      for sp in 0..3:
        for c in 0..2:
          let k = salt + uint64(i)*48 + uint64(sp)*12 + uint64(c)*4
          setC(v.a{i}[sp][c], u01(k), u01(k+1))
          setC(v.b{i}[sp][c], u01(k+2), u01(k+3))
  var stats = SiStats()
  var op = newShiftInvertOp[HF](
    applyM = proc (rr: var HF; x: HF) = cw.D(rr, x, msolve, rw),
    applyMdag = proc (rr: var HF; x: HF) = cw.Ddag(rr, x, msolve, rw),
    newVec = proc (): HF = newOneOf(proto),
    startVec = startVec,
    r2req = innerR2, maxits = innerMaxIts, stats = stats)
  let applyD = proc(r: var HF; x: HF) = cw.D(r, x, mass, rw)
  toc("setup")
  configLoop:
    let (st, meta) = load(g, fn)
    if st != 0: qexError "failed to load ", fn
    threads:
      g.reunit
    emit &"TS {icfg} {g.triangleSum:.8f}   # triangleSum, beta {meta.beta} traj {meta.traj}"
    var qflow, qflowT0 = 0.0
    if doflow:
      threads:
        gf := g
      fr.start(EQ(wt, gf).q)
      var nstep = 0
      gf.flow(nFlowSteps, floweps, cflow):
        inc nstep
        if nstep mod flowmeas == 0:
          let (e, q) = EQ(wt, gf)
          fr.add(wflowT, e, q)
      (qflow, qflowT0) = flowResult(icfg)
    if doeigs:
      let tw0 = getMonoTime()
      threads:
        gs := g
      stout(w, gs, rho, sc, nstout)
      threads:
        gs.setBC
      cw.gaugeRefresh
      eigsConfig(op, applyD, stats, icfg, qflow, qflowT0, tw0)
  finish()

proc runCubic(lo: Layout) =
  var g = lo.newGauge
  var gs = lo.newGauge
  var st = newStoutSmear(lo, rho)
  var cw = newCubicWilson(gs, cSW)
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
  var op = newShiftInvertOp[DF](
    applyM = proc (rr: var DF; x: DF) = cw.D(rr, x, msolve),
    applyMdag = proc (rr: var DF; x: DF) = cw.Ddag(rr, x, msolve),
    newVec = proc (): DF =
      result = newOneOf(proto)
      result := 0,
    startVec = startVec,
    r2req = innerR2, maxits = innerMaxIts, stats = stats)
  let applyD = proc(r: var DF; x: DF) = cw.D(r, x, mass)
  proc eqm(gg: auto): tuple[e, q: float] =
    let f = gg.fmunu 1
    let (es, et) = f.densityE
    (es + et, f.topoQ)
  toc("setup")
  configLoop:
    if 0 != g.loadGauge fn: qexError "failed to load ", fn
    threads:
      g.projectSU
      for mu in 0..<gs.len:
        gs[mu] := g[mu]
    block:
      var s = 0.0
      for v in g.plaq: s += v
      emit &"TS {icfg} {s:.8f}   # plaq, {fn.extractFilename}"
    var qflow, qflowT0 = 0.0
    if doflow:
      fr.start(eqm(g).q)
      var nstep = 0
      g.gaugeFlow(nFlowSteps, floweps):
        inc nstep
        if nstep mod flowmeas == 0:
          let (e, q) = eqm(g)
          fr.add(wflowT, e, q)
      (qflow, qflowT0) = flowResult(icfg)
    if doeigs:
      let tw0 = getMonoTime()
      st.smearN(gs, gs, nstout)
      gs.setBC
      cw.gaugeRefresh
      eigsConfig(op, applyD, stats, icfg, qflow, qflowT0, tw0)
  finish()

case lattice
of "hc":
  withLayout(simdlen, cgeom, lo):
    runHc(lo)
of "cubic":
  runCubic(cgeom.newLayout)
else:
  qexError "unknown -lattice:", lattice, " (hc | cubic)"
qexFinalize()
