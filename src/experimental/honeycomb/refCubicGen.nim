## refCubicGen -- quenched SU(3) configuration generation on the ordinary
## hypercubic lattice, Wilson plaquette action.
##
## This is the *cubic reference* leg of the 16-cell honeycomb project (PLAN.md
## task R, step 5): everything here uses committed QEX machinery only, so the
## resulting ensembles and the analysis chain built on them can be trusted as
## the yardstick for the honeycomb numbers.
##
## Structure follows `src/examples/puregaugehmc.nim`, stripped to the Wilson
## plaquette action and with the HMC diagnostics that PLAN.md asks for:
## acceptance rate and `<exp(-dH)>` (which must be 1 within errors -- this is
## the standard correctness check on the whole MD/Metropolis chain).
##
##   export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
##   cd build_mac
##   make src/experimental/honeycomb/refCubicGen.nim
##   ./bin/refCubicGen -lat:8,8,8,8 -beta:5.9 -ntraj:200 -nwarm:100 \
##                     -outdir:/tmp/b590 -savefreq:2
##
## Output: `<outdir>/cfg.NNNNN.lime` plus one `HMC` line per trajectory on
## stdout and a final `SUMMARY` line.

import qex, gauge, algorithms/integrator
import std/[math, os, strformat, strutils]
import hcanalysis

qexInit()
tic()

letParam:
  lat = @[8, 8, 8, 8]
  beta = 5.9
  ntraj = 100            ## measured trajectories (after nwarm)
  nwarm = 50             ## thermalisation trajectories, discarded
  tau = 1.0              ## MD trajectory length
  nsteps = 8             ## MD steps per trajectory
  seed: uint64 = 987654321'u64
  savefreq = 0           ## save a configuration every this many measured trajs (0 = never)
  outdir = ""            ## directory for saved configurations
  cfgprefix = "cfg"
  start = "hot"          ## "hot" (random) or "cold" (unit) initial configuration
  gintalg: IntegratorProc = "2MN"
  revCheckFreq = 0       ## reversibility check every N trajectories (0 = never)
  showTimers: bool = 0

installHelpParam()
echoParams()
echo "rank ", myRank, "/", nRanks
threads: echo "thread ", threadNum, "/", numThreads

let
  gc = GaugeActionCoeffs(plaq: beta)   # Wilson plaquette action
  lo = lat.newLayout
  vol = lo.physVol

var r = lo.newRNGField(MRG32k3a, seed)
var R: MRG32k3a
R.seed(seed, 987654321)

var
  g = lo.newGauge
  p = lo.newGauge
  f = lo.newGauge
  g0 = lo.newGauge

proc reunit(g: auto) =
  threads:
    g.projectSU

template pnorm2(p2: float) =
  threads:
    var p2t = 0.0
    for i in 0..<p.len:
      p2t += p[i].norm2
    threadMaster: p2 = p2t

proc gaction(g: auto, p2: float): auto =
  let
    ga = gc.gaugeAction1 g
    t = 0.5*p2 - float(16*vol)     # <T> subtraction: nd*(nc^2-1)/2 per site
    h = ga + t
  (ga, t, h)

proc mdt(dtau: float) =
  threads:
    for mu in 0..<g.len:
      for s in g[mu]:
        g[mu][s] := exp(dtau*p[mu][s])*g[mu][s]

proc mdv(dtau: float) =
  gc.gaugeForce(g, f)
  threads:
    for mu in 0..<f.len:
      p[mu] -= dtau*f[mu]

let
  (V, T) = newIntegratorPair(mdv, mdt)
  H = gintalg(T = T, V = V, steps = nsteps)

echo H

if start == "cold":
  g.unit
else:
  g.random r
g.reunit

if outdir.len > 0: createDir outdir

proc plaqAvg(g: auto): float =
  ## Average plaquette Re tr U_plaq / Nc, 1 for the unit gauge.
  ## QEX's `plaq` returns the 6 planes each already divided by the number of
  ## planes, so the average is their *sum* (cf. `mplaq` in puregaugehmc.nim).
  let pl = g.plaq
  var s = 0.0
  for v in pl: s += v
  s

var
  nacc = 0
  nprop = 0
  sumExpDH = 0.0
  dhs: seq[float] = @[]
  expdhs: seq[float] = @[]
  plaqs: seq[float] = @[]
  nsaved = 0

echo "# HMC traj dH exp(-dH) acc plaq"
toc("setup")

for n in 1-nwarm..ntraj:
  let measuring = n > 0
  var p2 = 0.0
  threads:
    p.randomTAH r
    for i in 0..<g.len:
      g0[i] := g[i]
  p2.pnorm2
  let (ga0, t0, h0) = g.gaction p2

  H.evolve tau
  H.finish

  p2.pnorm2
  let (ga1, t1, h1) = g.gaction p2

  if revCheckFreq > 0 and n mod revCheckFreq == 0:
    var gs = lo.newGauge
    var ps = lo.newGauge
    threads:
      for i in 0..<g.len:
        gs[i] := g[i]
        ps[i] := p[i]
        p[i] := -1*p[i]
    H.evolve tau
    H.finish
    var p2r = 0.0
    p2r.pnorm2
    let (_, _, hr) = g.gaction p2r
    echo &"REVCHECK traj {n} |dH_fwd+dH_bwd| = {abs(hr-h0):.3e}"
    threads:
      for i in 0..<g.len:
        g[i] := gs[i]
        p[i] := ps[i]

  let
    dH = h1 - h0
    acc = exp(-dH)
    accr = R.uniform
    accepted = accr <= acc
  if accepted:
    g.reunit
  else:
    threads:
      for i in 0..<g.len:
        g[i] := g0[i]

  let pl = g.plaqAvg
  if measuring:
    inc nprop
    if accepted: inc nacc
    sumExpDH += acc
    dhs.add dH
    expdhs.add acc
    plaqs.add pl
    echo &"HMC {n} {dH:.8f} {acc:.8f} {(if accepted: 1 else: 0)} {pl:.10f}"
    if savefreq > 0 and outdir.len > 0 and n mod savefreq == 0:
      let fn = outdir / &"{cfgprefix}.{n:05}.lime"
      if 0 != g.saveGauge fn:
        qexError "failed to save gauge to ", fn
      inc nsaved
  else:
    echo &"WARM {n} {dH:.8f} {acc:.8f} {(if accepted: 1 else: 0)} {pl:.10f}"

toc("hmc")

let
  accRate = if nprop > 0: nacc.float/nprop.float else: 0.0
  (edhM, edhE) = jackknifeMean(expdhs)
  tauPlaq = autocorrTime(plaqs)
  (plM, plE0) = jackknifeMean(plaqs)
  plE = plE0*sqrt(2.0*tauPlaq)

echo ""
echo "==== refCubicGen summary ===="
echo &"lat            {lat}   vol = {vol}"
echo &"beta           {beta}"
echo &"tau/nsteps     {tau} / {nsteps}   ({gintalg.repr.len > 0})"
echo &"trajectories   {nprop} measured ({nwarm} warm-up)"
echo &"acceptance     {accRate*100.0:.2f} %   ({nacc}/{nprop})"
echo &"<dH>           {dhs.mean:.6f} +- {stderrMean(dhs):.6f}"
echo &"<exp(-dH)>     {edhM:.6f} +- {edhE:.6f}    (must be 1 within errors)"
echo &"<plaq>         {plM:.8f} +- {plE:.8f}   tau_int = {tauPlaq:.2f}"
echo &"saved          {nsaved} configurations to '{outdir}'"
echo &"SUMMARY {beta} {lat.join(\"x\")} {nprop} {accRate:.4f} {edhM:.6f} {edhE:.6f} {plM:.8f} {plE:.8f} {nsaved}"

if showTimers: echoTimers()
qexFinalize()
