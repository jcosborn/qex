## hcpuregauge -- quenched configuration generation on the 16-cell honeycomb
## with the triangle action.
##
##   -algo:hb    Cabibbo-Marinari heatbath + norsweeps overrelaxation sweeps per update
##   -algo:hmc   HMC, one trajectory of length tau with nsteps steps of -gintalg
##   -algo:scan  <triangleSum> vs beta with the heatbath, fresh start per beta
##
## Per measured update:  HMC n dH exp(-dH) acc triangleSum  |  HB n triangleSum
## (WARM lines during thermalisation, SUMMARY line at the end).  triangleSum
## is taken from the action, S = (beta/2) 32 N_sites (1 - triangleSum).
## Configurations are saved every -savefreq measured updates as
## <outdir>/<cfgprefix>.NNNNN.lime (hcgauge.save); -loadcfg restarts from one.
## -simdlen:1 (or 2) for cell geometries the default SIMD layout rejects.

import qex
import std/[math, os, strformat, strutils]
import honeycomb, hcanalysis

qexInit()
tic()

letParam:
  geom = @[4, 4, 4, 4]   ## cell geometry (2 sites, 24 links per cell)
  beta = 8.0
  algo = "hb"            ## hb | hmc | scan
  ntraj = 100            ## measured updates (after nwarm)
  nwarm = 100            ## thermalisation updates, discarded
  norsweeps = 3          ## OR sweeps per heatbath update
  tau = 1.0              ## HMC trajectory length
  nsteps = 10            ## HMC steps per trajectory
  gintalg = "2MN"        ## HMC integrator: leapfrog | 2MN | 4MN5FV
  seed: uint64 = 987654321'u64
  start = "hot"          ## hot (random) | cold (unit) | warm (warmSU 0.5)
  loadcfg = ""           ## start from a saved configuration (overrides start)
  savefreq = 0           ## save every this many measured updates (0 = never)
  outdir = ""
  cfgprefix = "hc"
  measfreq = 1           ## measure triangleSum every this many updates
  revCheckFreq = 0       ## HMC reversibility check every N trajectories (0 = never)
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

proc saveCfg(g: auto, n: int) =
  if savefreq > 0 and outdir.len > 0 and n mod savefreq == 0:
    let fn = outdir / &"{cfgprefix}.{n:05}.lime"
    if 0 != g.save(fn, beta = beta, traj = n, info = "hcpuregauge " & algo):
      qexError "failed to save configuration to ", fn

template setStart(g, r: untyped) =
  if loadcfg.len > 0:
    let lr = g.load(loadcfg)
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

proc summary(tss: seq[float]): tuple[m, e, tau: float] =
  let tau = if tss.len > 4: autocorrTime(tss) else: 1.0
  let (m, e0) = jackknifeMean(tss)
  (m, e0*sqrt(max(1.0, 2.0*tau)), tau)

proc runAll(lo: Layout) =
  var r = lo.newRNGField(MRG32k3a, seed)
  var R: MRG32k3a          # lattice-global RNG for the accept step
  R.seed(seed, 987654321)
  var g = newHcGauge(lo)
  var w = newActionWork(g)
  let ncells = lo.physVol

  case algo
  of "hmc":
    var h = newHcHmc(g, w, beta, tau, nsteps, r, R, gintalg)
    setStart(g, r)
    var dhs, tss: seq[float]
    echo "# HMC n dH exp(-dH) acc triangleSum"
    toc("setup")
    for n in 1-nwarm..ntraj:
      h.update
      let accepted = if h.accepted: 1 else: 0
      if revCheckFreq > 0 and n mod revCheckFreq == 0:
        let rc = h.revCheck
        echo &"REVCHECK {n} |dHf+dHb| = {abs(rc.sumdH):.3e}  linkdiff = {rc.linkDiff:.3e}"
      if n > 0:
        dhs.add h.deltaH
        if n mod measfreq == 0:
          let ts = triSum(w, beta, g)
          tss.add ts
          echo &"HMC {n} {h.deltaH:.8f} {h.expmDeltaH:.8f} {accepted} {ts:.10f}"
        saveCfg(g, n)
      else:
        echo &"WARM {n} {h.deltaH:.8f} {h.expmDeltaH:.8f} {accepted}"
        if n == 0: h.clearStats
    toc("hmc")
    let (tsM, tsE, tauTs) = summary(tss)
    var edh = 0.0
    var edh2 = 0.0
    for s in h.stats:
      let e = exp(s.hOld - s.hNew)
      edh += e
      edh2 += e*e
    let nprop = h.nUpdates
    edh /= nprop.float
    let edhE = sqrt(max(0.0, edh2/nprop.float - edh*edh)/nprop.float)
    echo ""
    echo "==== hcpuregauge HMC summary ===="
    echo &"geom           {geom}   cells = {ncells}  links = {nDirs*ncells}"
    echo &"beta           {beta}   tau/nsteps = {tau}/{nsteps}  ({gintalg})"
    echo &"updates        {nprop} measured ({nwarm} warm-up)"
    echo &"acceptance     {h.acceptRatio*100.0:.2f} %   ({h.nAccepts}/{nprop})"
    echo &"<dH>           {h.avgDeltaH:.6f} +- {stderrMean(dhs):.6f}"
    echo &"<exp(-dH)>     {edh:.6f} +- {edhE:.6f}    (must be 1 within errors)"
    echo &"<triangleSum>  {tsM:.8f} +- {tsE:.8f}   tau_int = {tauTs:.2f}"
    echo &"force calls    {h.nForce}"
    echo &"SUMMARY hmc {beta} {geom.join(\"x\")} {nprop} {h.acceptRatio:.4f} {edh:.6f} {edhE:.6f} {tsM:.8f} {tsE:.8f}"

  of "hb":
    var hb = newHcHeatbath(g, w, beta)
    setStart(g, r)
    var tss: seq[float]
    echo "# HB n triangleSum"
    toc("setup")
    for n in 1-nwarm..ntraj:
      hb.update(g, r, norsweeps)
      if n > 0:
        if n mod measfreq == 0:
          let ts = triSum(w, beta, g)
          tss.add ts
          echo &"HB {n} {ts:.10f}"
        saveCfg(g, n)
      elif n mod 10 == 0:
        echo &"WARM {n} {triSum(w, beta, g):.10f}"
    threads:
      g.reunit
    toc("hb")
    let (tsM, tsE, tauTs) = summary(tss)
    echo ""
    echo "==== hcpuregauge heatbath summary ===="
    echo &"geom           {geom}   cells = {ncells}  links = {nDirs*ncells}"
    echo &"beta           {beta}   1 HB + {norsweeps} OR sweeps per update"
    echo &"updates        {tss.len} measured ({nwarm} warm-up)"
    echo &"<triangleSum>  {tsM:.8f} +- {tsE:.8f}   tau_int = {tauTs:.2f}"
    echo &"SUMMARY hb {beta} {geom.join(\"x\")} {tss.len} {tsM:.8f} {tsE:.8f} {tauTs:.3f}"

  of "scan":
    var hb = newHcHeatbath(g, w, beta)
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
      setStart(g, r)
      var tss: seq[float]
      for n in 1-nwarm..ntraj:
        hb.update(g, r, norsweeps)
        if n > 0 and n mod measfreq == 0:
          tss.add triSum(w, b, g)
      let (tsM, tsE, tauTs) = summary(tss)
      echo &"SCAN {b:.4f} {tsM:.8f} {tsE:.8f} {tauTs:.2f}"
      if haveScanOut:
        fh.writeLine &"{b:.4f} {tsM:.8f} {tsE:.8f} {tauTs:.2f}"
        fh.flushFile
    if haveScanOut: fh.close
    toc("scan")

  else:
    qexError "unknown -algo:", algo, " (hb | hmc | scan)"

withLayout(simdlen, geom, lo):
  runAll(lo)

if showTimers: echoTimers()
qexFinalize()
