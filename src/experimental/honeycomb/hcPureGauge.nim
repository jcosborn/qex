## hcPureGauge -- quenched SU(3) configuration generation on the 16-cell
## honeycomb with the triangle action (task **M**).
##
## Algorithms:
##   -algo:hb    Cabibbo-Marinari heatbath + overrelaxation (hcheatbath) --
##               what the paper used; one "update" = 1 HB + norsweeps OR
##               sweeps.
##   -algo:hmc   HMC (hchmc); one "update" = one trajectory of length tau
##               with nsteps MD steps of the -gintalg integrator
##               (leapfrog | 2MN | 4MN5FV).
##   -algo:scan  quick <triangleSum> vs beta scan with the heatbath (used
##               once to produce doc/plots/hb_scan.dat for task R); writes
##               "beta <ts> err tau_int" lines to -scanout.
##
## Per measured update one parseable line is printed:
##   HMC <n> <dH> <exp(-dH)> <acc01> <triangleSum>
##   HB <n> <triangleSum>
## plus WARM lines during thermalisation and a final SUMMARY line.
## `triangleSum` here is computed allocation free from the action work space:
## S = (beta/2) 32 nSites (1 - triangleSum)  (cross-checked against Task L's
## `triangleSum` in tests/taction.nim to 4e-16).
##
## Configurations are saved with hcio every -savefreq measured updates:
##   <outdir>/<cfgprefix>.NNNNN.lime      (format: see hcio.nim)
## and a run can be started from a saved configuration with -loadcfg.
##
## -simdlen: QEX's default SIMD layout (VLEN = 4) rejects some cell
## geometries (e.g. 6^4: no valid inner split).  Pass -simdlen:1 (or 2) to
## fall back to a shorter SIMD length for such lattices; 0 means the build
## default.  Results are statistically equivalent but not bit identical
## across simdlen (the site -> RNG-stream pairing changes).
##
## Example:
##   export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
##   cd build_mac && make src/experimental/honeycomb/hcPureGauge.nim
##   OMP_NUM_THREADS=4 ./bin/hcPureGauge -geom:8,8,8,8 -beta:8.0 -algo:hb \
##     -nwarm:200 -ntraj:1000 -norsweeps:3 -savefreq:10 -outdir:/tmp/hc80

import qex
import std/[math, os, strformat, strutils]
import hcgeom, hclayout, hcgauge, hcaction, hchmc, hcheatbath, hcio
import hcanalysis

qexInit()
tic()

letParam:
  geom = @[4, 4, 4, 4]   ## cell geometry (sites = 2x, links = 24x)
  beta = 8.0
  algo = "hb"            ## hb | hmc | scan
  ntraj = 100            ## measured updates (after nwarm)
  nwarm = 100            ## thermalisation updates, discarded
  norsweeps = 3          ## OR sweeps per heatbath update (hb)
  tau = 1.0              ## HMC trajectory length
  nsteps = 10            ## HMC MD steps per trajectory
  gintalg = "2MN"        ## HMC integrator: leapfrog | 2MN | 4MN5FV
  seed: uint64 = 987654321'u64
  start = "hot"          ## hot (random) | cold (unit) | warm (warmSU 0.5)
  loadcfg = ""           ## start from a saved configuration (overrides start)
  savefreq = 0           ## save every this many measured updates (0 = never)
  outdir = ""
  cfgprefix = "hc"
  measfreq = 1           ## measure triangleSum every this many updates
  revCheckFreq = 0       ## HMC reversibility check every N trajs (0 = never)
  simdlen = 0            ## 0 = build default (VLEN); 1 or 2 for odd geometries
  scanbeta0 = 1.0        ## scan: first beta
  scanbeta1 = 12.0       ## scan: last beta
  scannbeta = 23         ## scan: number of beta points
  scanout = "hb_scan.dat"
  showTimers: bool = 0

installHelpParam()
echoParams()
echo "rank ", myRank, "/", nRanks
threads: echo "thread ", threadNum, "/", numThreads

if outdir.len > 0: createDir outdir

# triangleSum from the action (allocation free given a work object)
template triSumOf(w, beta, g: untyped): float =
  1.0 - hcAction(w, beta, g)/(0.5*beta*float(nTriPerSite*g.hl.nSites))

proc saveCfg(g: auto, n: int) =
  if savefreq > 0 and outdir.len > 0 and n mod savefreq == 0:
    let fn = outdir / &"{cfgprefix}.{n:05}.lime"
    if 0 != g.saveHcGauge(fn, beta = beta, traj = n,
                          info = "hcPureGauge " & algo):
      qexError "failed to save configuration to ", fn

template setStart(g, r: untyped) =
  if loadcfg.len > 0:
    let lr = g.loadHcGauge(loadcfg)
    if lr.status != 0:
      qexError "failed to load configuration from ", loadcfg
    echo "# loaded ", loadcfg, ": beta ", lr.meta.beta, " traj ", lr.meta.traj,
         " info '", lr.meta.info, "'"
  elif start == "cold":
    g.unit
  elif start == "warm":
    threads:
      g.warm(0.5, r)
  else:
    threads:
      g.random r

proc runAll(hl: auto) =
  var r = hl.lo.newRNGField(MRG32k3a, seed)
  var R: MRG32k3a          # lattice-global RNG for the Metropolis decision
  R.seed(seed, 987654321)
  var g = newHcGauge(hl)

  case algo
  of "hmc":
    var h = newHcHmc(g, beta, tau, nsteps, gintalg)
    echo &"# HMC integrator {gintalg} nsteps {nsteps}: ",
         &"{nForcePerTraj(h.sched)} force calls/traj"
    setStart(g, r)
    var
      nacc, nprop = 0
      dhs, expdhs, tss: seq[float]
    echo "# HMC n dH exp(-dH) acc triangleSum"
    toc("setup")
    for n in 1-nwarm..ntraj:
      let measuring = n > 0
      let (dH, acc, accepted) = h.trajectory(g, r, R)
      if revCheckFreq > 0 and n mod revCheckFreq == 0:
        let rc = h.revCheck(g, r)
        echo &"REVCHECK {n} |dHf+dHb| = {abs(rc.sumdH):.3e}  linkdiff = {rc.linkDiff:.3e}"
      if measuring:
        inc nprop
        if accepted: inc nacc
        dhs.add dH
        expdhs.add acc
        if n mod measfreq == 0:
          let ts = triSumOf(h.w, beta, g)
          tss.add ts
          echo &"HMC {n} {dH:.8f} {acc:.8f} {(if accepted: 1 else: 0)} {ts:.10f}"
        saveCfg(g, n)
      else:
        echo &"WARM {n} {dH:.8f} {acc:.8f} {(if accepted: 1 else: 0)}"
    toc("hmc")
    let
      accRate = if nprop > 0: nacc.float/nprop.float else: 0.0
      (edhM, edhE) = jackknifeMean(expdhs)
      tauTs = if tss.len > 4: autocorrTime(tss) else: 1.0
      (tsM, tsE0) = jackknifeMean(tss)
      tsE = tsE0*sqrt(max(1.0, 2.0*tauTs))
    echo ""
    echo "==== hcPureGauge HMC summary ===="
    echo &"geom           {geom}   cells = {hl.nCells}  links = {hl.nLinks}"
    echo &"beta           {beta}   tau/nsteps = {tau}/{nsteps}  ({gintalg})"
    echo &"updates        {nprop} measured ({nwarm} warm-up)"
    echo &"acceptance     {accRate*100.0:.2f} %   ({nacc}/{nprop})"
    echo &"<dH>           {dhs.mean:.6f} +- {stderrMean(dhs):.6f}"
    echo &"<exp(-dH)>     {edhM:.6f} +- {edhE:.6f}    (must be 1 within errors)"
    echo &"<triangleSum>  {tsM:.8f} +- {tsE:.8f}   tau_int = {tauTs:.2f}"
    echo &"force calls    {h.nForce}"
    echo &"SUMMARY hmc {beta} {geom.join(\"x\")} {nprop} {accRate:.4f} {edhM:.6f} {edhE:.6f} {tsM:.8f} {tsE:.8f}"

  of "hb":
    var hb = newHcHeatbath(g, beta)
    setStart(g, r)
    var tss: seq[float]
    echo "# HB n triangleSum"
    toc("setup")
    for n in 1-nwarm..ntraj:
      hb.update(g, r, norsweeps)
      if n > 0:
        if n mod measfreq == 0:
          let ts = triSumOf(hb.w, beta, g)
          tss.add ts
          echo &"HB {n} {ts:.10f}"
        saveCfg(g, n)
      elif n mod 10 == 0:
        echo &"WARM {n} {triSumOf(hb.w, beta, g):.10f}"
    # keep unitarity drift down on long streams
    threads:
      g.reunit
    toc("hb")
    let
      tauTs = if tss.len > 4: autocorrTime(tss) else: 1.0
      (tsM, tsE0) = jackknifeMean(tss)
      tsE = tsE0*sqrt(max(1.0, 2.0*tauTs))
    echo ""
    echo "==== hcPureGauge heatbath summary ===="
    echo &"geom           {geom}   cells = {hl.nCells}  links = {hl.nLinks}"
    echo &"beta           {beta}   1 HB + {norsweeps} OR sweeps per update"
    echo &"updates        {tss.len} measured ({nwarm} warm-up)"
    echo &"<triangleSum>  {tsM:.8f} +- {tsE:.8f}   tau_int = {tauTs:.2f}"
    echo &"SUMMARY hb {beta} {geom.join(\"x\")} {tss.len} {tsM:.8f} {tsE:.8f} {tauTs:.3f}"

  of "scan":
    ## <triangleSum> vs beta, heatbath, fresh -start (default hot) per beta
    var hb = newHcHeatbath(g, beta)
    var fh: File
    let haveScanOut = scanout.len > 0 and myRank == 0
    if haveScanOut:
      fh = open(scanout, fmWrite)
      fh.writeLine "# 16-cell honeycomb triangle action: <triangleSum> vs beta"
      fh.writeLine &"# heatbath 1 HB + {norsweeps} OR per update, geom {geom}, " &
        &"{nwarm} warm + {ntraj} measured updates per point, seed {seed}"
      fh.writeLine "# beta <triangleSum> err tau_int"
    echo "# SCAN beta triangleSum err tau_int"
    for ib in 0..<scannbeta:
      let b = scanbeta0 + (scanbeta1 - scanbeta0)*ib.float/(scannbeta-1).float
      hb.beta = b
      setStart(g, r)   # fresh start at each beta: no hysteresis
      var tss: seq[float]
      for n in 1-nwarm..ntraj:
        hb.update(g, r, norsweeps)
        if n > 0 and n mod measfreq == 0:
          tss.add triSumOf(hb.w, b, g)
      let
        tauTs = if tss.len > 4: autocorrTime(tss) else: 1.0
        (tsM, tsE0) = jackknifeMean(tss)
        tsE = tsE0*sqrt(max(1.0, 2.0*tauTs))
      echo &"SCAN {b:.4f} {tsM:.8f} {tsE:.8f} {tauTs:.2f}"
      if haveScanOut:
        fh.writeLine &"{b:.4f} {tsM:.8f} {tsE:.8f} {tauTs:.2f}"
        fh.flushFile
    if haveScanOut: fh.close
    toc("scan")

  else:
    qexError "unknown -algo:", algo, " (hb | hmc | scan)"

case simdlen
of 0: runAll(newHcLayout(geom))
of 1: runAll(newHcLayoutX(geom, 1))
of 2: runAll(newHcLayoutX(geom, 2))
else: qexError "unsupported -simdlen:", simdlen, " (0 = default, 1, 2)"

if showTimers: echoTimers()
qexFinalize()
