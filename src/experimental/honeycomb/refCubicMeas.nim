## refCubicMeas -- gradient flow, t0, and topological charge on saved cubic
## SU(3) configurations, plus the ensemble analysis and the continuum
## extrapolation of slide 10.
##
## Three modes:
##
## 1. measurement (default):  flow every configuration matched by `-cfgs`,
##    record the flow history `(t, t^2E, t d/dt[t^2E], Q)`, extract `t0` and
##    `Q(t0)`, and print one `CFG` line per configuration followed by an
##    `ENSEMBLE` line with the jackknifed `t0`, `<Q^2>`, `chi` and
##    `10^4 t0^2 chi`.
##
## 2. `-flowdump:<file>`:  in addition, write the full flow history of the
##    first configuration (used for the `topoQ` normalisation study and for
##    the flow figure).
##
## 3. `-chifit:<file>`:  skip the lattice entirely, read a table of
##    `a^2/t0  10^4 t0^2 chi  err` and do the O(a^2) and O(a^4) continuum
##    extrapolations with `hcanalysis.fitPoly`.
##
##   export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
##   cd build_mac
##   make src/experimental/honeycomb/refCubicMeas.nim
##   OMP_NUM_THREADS=6 ./bin/refCubicMeas -cfgs:'/tmp/b590/cfg.*.lime' -eps:0.08
##
## Conventions (FORMULATION.md sec. 4):
##   E from the clover F_munu (`gauge/gaugeUtils.fmunu`, `densityE`),
##   t0 from `t^2<E> = 0.3`, w0^2 from `t d/dt[t^2<E>] = 0.3`,
##   Q = `topoQ(f)` -- normalisation derived in `doc/RESULTS_CUBIC.md`,
##   chi = <Q^2>/V with V = prod(lat) in lattice units.

import qex, gauge, gauge/wflow
import std/[math, os, algorithm, strformat, strutils]
import hcanalysis

qexInit()
tic()

letParam:
  abeliantest: bool = 0  ## exact check of the `topoQ` normalisation, see below
  n1 = 1                 ## flux quantum in the (0,1) plane   (abeliantest)
  n2 = 1                 ## flux quantum in the (2,3) plane   (abeliantest)
  lat = @[8, 8, 8, 8]    ## lattice for the abelian test (configs set their own)
  cfgs = ""              ## glob pattern for the configuration files
  eps = 0.05             ## Wilson-flow RK3 step
  tmax = 0.0             ## stop flowing at this t (0 = use t2Estop only)
  t2Estop = 0.45         ## stop flowing once t^2<E> exceeds this
  measevery = 2          ## measure E and Q every this many flow steps
  fmunuloop = 1          ## clover size for fmunu: 1, 3, 4 or 5
  t0target = 0.3         ## t^2<E> = t0target defines t0
  w0target = 0.3         ## t d/dt[t^2<E>] = w0target defines w0^2
  interpOrder = 3        ## 1 = linear, 3 = cubic interpolation for the crossings
  beta = 0.0             ## for labelling only
  bin = 1                ## jackknife bin size (in units of saved configurations)
  flowdump = ""          ## write the full flow history of the first config here
  outfile = ""           ## write the CFG lines here as well as to stdout
  chifit = ""            ## continuum-extrapolation mode: fit this table instead
  showTimers: bool = 0

installHelpParam()
echoParams()

# ---------------------------------------------------------------------------
# mode 0: exact check of QEX's topoQ normalisation (`-abeliantest`)
# ---------------------------------------------------------------------------
# Build the constant-field-strength SU(3) configuration in the Cartan direction
# T = diag(1,-1,0):
#
#   U_1(x) = exp(i phi1 x0 T)                       phi1 = 2 pi n1/(L0 L1)
#   U_0(x) = 1  except at x0 = L0-1 where it is exp(-i phi1 L0 x1 T)
#   U_3(x) = exp(i phi2 x2 T)                       phi2 = 2 pi n2/(L2 L3)
#   U_2(x) = 1  except at x2 = L2-1 where it is exp(-i phi2 L2 x3 T)
#
# Every (0,1) plaquette is then exactly exp(i phi1 T) and every (2,3) plaquette
# exactly exp(i phi2 T); all other plaquettes are 1.  Periodicity requires
# exp(-2 pi i n_k T) = 1, i.e. integer n_k -- that is the flux quantisation.
#
# The gauge field is a direct sum of three U(1) line bundles with charges
# q = (1,-1,0), so by Atiyah-Singer the exact topological charge is
#
#   Q_exact = sum_i q_i^2 n1 n2 = 2 n1 n2 .
#
# The clover (loop=1) is exact here: fmunu gives F_01 = i sin(phi1) T and
# F_23 = i sin(phi2) T, so QEX's expression must return
#
#   topoQ = (V/(2 pi^2)) sin(phi1) sin(phi2)   ->   2 n1 n2 as phi -> 0.
#
# and likewise  E = 2 sin^2(phi1) + 2 sin^2(phi2).
if abeliantest:
  let
    lo = lat.newLayout
    vol = lo.physVol
  var g = lo.newGauge
  const nc = g[0][0].nrows
  when nc != 3:
    qexError "abeliantest assumes Nc = 3"
  let
    phi1 = 2.0*PI*n1.float/float(lat[0]*lat[1])
    phi2 = 2.0*PI*n2.float/float(lat[2]*lat[3])
  g.unit
  for i in lo.sites:
    let
      x0 = lo.coords[0][i]
      x1 = lo.coords[1][i]
      x2 = lo.coords[2][i]
      x3 = lo.coords[3][i]
    template setPhase(m: untyped, th: float) =
      ## m := diag(e^{i th}, e^{-i th}, 1)
      m[0, 0].re := cos(th)
      m[0, 0].im := sin(th)
      m[1, 1].re := cos(th)
      m[1, 1].im := -sin(th)
    setPhase(g[1]{i}, phi1*x0.float)
    if x0 == lat[0]-1:
      setPhase(g[0]{i}, -phi1*float(lat[0]*x1))
    setPhase(g[3]{i}, phi2*x2.float)
    if x2 == lat[2]-1:
      setPhase(g[2]{i}, -phi2*float(lat[2]*x3))

  let
    f = g.fmunu 1
    (es, et) = f.densityE
    q = f.topoQ
    qexact = 2.0*float(n1*n2)
    qpredict = float(vol)/(2.0*PI*PI)*sin(phi1)*sin(phi2)
    epredict = 2.0*sin(phi1)*sin(phi1) + 2.0*sin(phi2)*sin(phi2)
  echo ""
  echo "==== topoQ normalisation: exact abelian test ===="
  echo &"lattice {lat}  V = {vol}   fluxes n1 = {n1}, n2 = {n2}"
  echo &"phi1 = {phi1:.10f}   phi2 = {phi2:.10f}"
  echo &"plaquettes (should be 1 except the 01 and 23 planes): {g.plaq}"
  echo &"E  measured {es+et:.12f}   predicted {epredict:.12f}   diff {abs(es+et-epredict):.3e}"
  echo &"QEX topoQ            {q:.12f}"
  echo &"closed-form lattice  {qpredict:.12f}   diff {abs(q-qpredict):.3e}"
  echo &"exact (index thm)    {qexact:.12f}   diff {abs(q-qexact):.3e}   ratio {q/qexact:.10f}"
  echo &"2 x topoQ            {2.0*q:.12f}   ratio to exact {2.0*q/qexact:.10f}"
  let ok = abs(q/qexact - 1.0) < 0.02
  echo(if ok: "VERDICT: QEX topoQ is correctly normalised (factor 1)."
       else: "VERDICT: QEX topoQ is NOT correctly normalised -- investigate.")
  qexFinalize()
  quit(if ok: 0 else: 1)

# ---------------------------------------------------------------------------
# mode 3: continuum extrapolation of an existing table
# ---------------------------------------------------------------------------
if chifit.len > 0:
  let rows = readColumns chifit
  var x, y, dy: seq[float]
  for r in rows:
    if r.len >= 3:
      x.add r[0]
      y.add r[1]
      dy.add r[2]
  echo ""
  echo "== continuum extrapolation of 10^4 t0^2 chi vs a^2/t0 =="
  echo &"{x.len} points, x in [{x.min:.4f}, {x.max:.4f}]"
  for (name, powers) in [("O(a^2)  c0 + c1 x", @[0, 1]),
                         ("O(a^4)  c0 + c1 x + c2 x^2", @[0, 1, 2])]:
    if x.len <= powers.len: continue
    let (co, er, cd) = fitPoly(x, y, dy, powers)
    var s = ""
    for k in 0..<co.len: s &= &"  c{k} = {co[k]:8.4f} +- {er[k]:.4f}"
    echo &"FIT {name}"
    echo &"    {s}"
    echo &"    chi^2/dof = {cd:.3f}"
    echo &"CONTINUUM {powers.len-1} {co[0]:.5f} {er[0]:.5f} {cd:.4f}"
  qexFinalize()
  quit(0)

# ---------------------------------------------------------------------------
# modes 1 and 2: measure
# ---------------------------------------------------------------------------
var files: seq[string] = @[]
for f in walkPattern(cfgs):
  files.add f
files.sort
if files.len == 0:
  qexError "no configuration files matched: '", cfgs, "'"
echo "found ", files.len, " configurations"

let
  clat = getFileLattice files[0]
  lo = clat.newLayout
  vol = lo.physVol
echo "lattice ", clat, "  vol ", vol
threads: echo "thread ", threadNum, "/", numThreads

var g = lo.newGauge

# flow history of the current configuration
var
  ft: seq[float] = @[]     # flow time
  ft2E: seq[float] = @[]   # t^2 <E>
  fq: seq[float] = @[]     # Q

proc eqm(g: auto, loop: int): auto =
  ## (E_s, E_t, Q) from the clover field strength.
  let
    f = g.fmunu loop
    (es, et) = f.densityE
    q = f.topoQ
  (es, et, q)

var
  t0s: seq[float] = @[]
  w0sqs: seq[float] = @[]
  qs: seq[float] = @[]       # Q at t0
  q2s: seq[float] = @[]
  qends: seq[float] = @[]    # Q at the last flow time
  names: seq[string] = @[]
  outLines: seq[string] = @[]

echo "# CFG name t0 w0sq Q(t0) Q(tend) tend roundDist"
toc("setup")

for icfg, fn in files:
  ft.setLen 0
  ft2E.setLen 0
  fq.setLen 0
  if 0 != g.loadGauge fn:
    qexError "failed to load ", fn
  threads: g.projectSU
  block:
    let (es, et, q) = g.eqm fmunuloop
    ft.add 0.0
    ft2E.add 0.0
    fq.add q
    discard es
    discard et
  var nstep = 0
  g.gaugeFlow(0, eps):
    inc nstep
    if nstep mod measevery == 0:
      let (es, et, q) = g.eqm fmunuloop
      let t2e = wflowT*wflowT*(es+et)
      ft.add wflowT
      ft2E.add t2e
      fq.add q
      if (tmax > 0.0 and wflowT >= tmax) or
         (t2Estop > 0.0 and t2e > t2Estop):
        break

  let
    w = derivT2E(ft, ft2E)
    t0 = findT0(ft, ft2E, t0target, interpOrder)
    w0sq = findW0(ft, w, w0target, interpOrder)
  if t0 <= 0.0:
    echo "WARNING: t^2E never reached ", t0target, " on ", fn, " -- skipped"
    continue
  let
    q0 = interpAt(ft, fq, t0)
    qend = fq[^1]
    tend = ft[^1]
    rd = abs(qend - round(qend))
  t0s.add t0
  w0sqs.add w0sq
  qs.add q0
  q2s.add q0*q0
  qends.add qend
  names.add fn.extractFilename
  let line = &"CFG {fn.extractFilename} {t0:.6f} {w0sq:.6f} {q0:.6f} {qend:.6f} {tend:.4f} {rd:.6f}"
  echo line
  outLines.add line

  if flowdump.len > 0 and icfg == 0 and myRank == 0:
    var fh = open(flowdump, fmWrite)
    fh.write("# flow history of " & fn & "\n")
    # Preserve the legacy duplicate charge columns for existing consumers.
    fh.write("# t  t^2E  W=t d/dt[t^2E]  Q_corrected  Q_qexTopoQ  dist(Q_corr,int)  dist(Q_raw,int)\n")
    for i in 0..<ft.len:
      let a = abs(fq[i]-round(fq[i]))
      fh.write(&"{ft[i]:.6f} {ft2E[i]:.8f} {w[i]:.8f} {fq[i]:.8f} {fq[i]:.8f} {a:.6f} {a:.6f}\n")
    fh.close
    echo "wrote flow history to ", flowdump

toc("measure")

if outfile.len > 0 and myRank == 0:
  var fh = open(outfile, fmWrite)
  fh.write("# CFG name t0 w0sq Q(t0) Q(tend) tend roundDist\n")
  for l in outLines: fh.write(l & "\n")
  fh.close

# ---------------------------------------------------------------------------
# ensemble analysis
# ---------------------------------------------------------------------------
let n = t0s.len
if n == 0:
  qexError "no usable measurements"

let
  (t0m, t0e0) = jackknifeMean(t0s, bin)
  tauT0 = autocorrTime(t0s)
  tauQ = autocorrTime(qs)
  tauQ2 = autocorrTime(q2s)
  (q2m, q2e0) = jackknifeMean(q2s, bin)

# combined jackknife of the two derived quantities, which share configurations
proc chiT0sq(idx: openArray[int]): float =
  ## 10^4 * t0^2 * <Q^2>/V for the subsample `idx`
  var st = 0.0
  var sq = 0.0
  for i in idx:
    st += t0s[i]
    sq += q2s[i]
  let t0b = st/idx.len.float
  let q2b = sq/idx.len.float
  1.0e4*t0b*t0b*q2b/vol.float

proc a2overT0(idx: openArray[int]): float =
  var st = 0.0
  for i in idx: st += t0s[i]
  idx.len.float/st

var allIdx = newSeq[int](n)
for i in 0..<n: allIdx[i] = i

let
  (chiM, chiE) = jackknife(allIdx, chiT0sq, bin)
  (xM, xE) = jackknife(allIdx, a2overT0, bin)

# integrated-autocorrelation-corrected naive errors, for comparison
let
  t0e = t0e0*sqrt(2.0*tauT0)
  q2e = q2e0*sqrt(2.0*tauQ2)

var qhist: array[-12..12, int]
for q in qends:
  let k = int(round(q))
  if k >= -12 and k <= 12: inc qhist[k]

var roundDistMax = 0.0
var roundDistAvg = 0.0
for q in qends:
  let d = abs(q-round(q))
  roundDistAvg += d
  if d > roundDistMax: roundDistMax = d
roundDistAvg = roundDistAvg/qends.len.float

echo ""
echo "==== refCubicMeas ensemble summary ===="
echo &"lattice        {clat}   V = {vol}"
echo &"beta           {beta}"
echo &"configs        {n}   (jackknife bin = {bin})"
echo &"t0/a^2         {t0m:.5f} +- {jackknifeMean(t0s, bin).err:.5f}  (binned jk)"
echo &"               naive jk {t0e0:.5f}, x sqrt(2 tau_int) = {t0e:.5f}"
echo &"w0^2/a^2       {w0sqs.mean:.5f} +- {jackknifeMean(w0sqs, bin).err:.5f}"
echo &"a^2/t0         {xM:.5f} +- {xE:.5f}"
echo &"<Q>            {qs.mean:.4f} +- {jackknifeMean(qs, bin).err:.4f}"
echo &"<Q^2>(t0)      {q2m:.4f} +- {jackknifeMean(q2s, bin).err:.4f}  (naive jk {q2e0:.4f}, tau-corrected {q2e:.4f})"
echo &"chi = <Q^2>/V  {q2m/vol.float:.6e}"
echo &"10^4 t0^2 chi  {chiM:.4f} +- {chiE:.4f}"
echo &"tau_int        t0 {tauT0:.2f}   Q {tauQ:.2f}   Q^2 {tauQ2:.2f}   (in units of the save interval)"
echo &"L/sqrt(t0)     {clat[0].float/sqrt(t0m):.2f}   (needs >~ 9 for t0 to exist, see RESULTS_CUBIC.md)"
echo &"Q at t_end     mean dist to nearest integer {roundDistAvg:.4f}, max {roundDistMax:.4f}"
var hs = ""
for k in -12..12:
  if qhist[k] > 0: hs &= &" {k}:{qhist[k]}"
echo &"Q histogram   {hs}"
echo &"ENSEMBLE {beta} {clat.join(\"x\")} {n} {t0m:.6f} {jackknifeMean(t0s, bin).err:.6f} {xM:.6f} {xE:.6f} {q2m:.6f} {jackknifeMean(q2s, bin).err:.6f} {chiM:.6f} {chiE:.6f} {tauQ:.3f} {tauQ2:.3f}"

if showTimers: echoTimers()
qexFinalize()
