## Staggered HMC with Hasenbusch mass ratios as examples/staghmc_sh, with the
## links, momenta and fermion fields on the GPU and every step of the
## molecular dynamics in kernels: the solves and one link forces
## (physics/stagGpu), the gauge action and force (gauge/gaugeGpu), the HYP
## smearing and its force (gauge/hypGpu).  It takes the parameters of
## staghmc_sh, with its species, integrators and force gradient steps,
## reversibility check and pbp, draws the random numbers of its generators
## in its order (rng/rngGpu), and prints its measurements, so that
## tests/extra/tstaghmc_sh runs it as staghmc_sh.  The links go to the host
## for saving only.  The smeared links are U(3), so with smearing the solvers keep
## the determinant of each link with rows 0 and 1, 14 reals.  Beyond the
## parameters of staghmc_sh:
##   -rg: ranks per dimension, lanes as in bestagcg
##   -mixed:1 mixed precision solves; -reals: reals per link of the solvers,
##   12 for SU(3) links, 14 by default with smearing, 18
##   -batch:0 solves the fermion terms that share a force step one by one
##   instead of together (solveM of several systems, nBatch per hop)
##   -check:1 first compares the smeared links and the smearing force with
##   gauge/hypsmear2 on random links
import qex, gauge, gauge/[hypsmear, hypsmear2, gaugeGpu, hypGpu], physics/[qcdTypes, stagSolve, stagGpu]
import backend/accel, rng/rngGpu
import mdevolve
import algorithm, math, os, sequtils, strutils, times

type IntProc = proc(T,V:Integrator; steps:int):Integrator
converter toIntProc(s:string):IntProc =
  template mkProc1(s:untyped):IntProc =
    proc mkInt(T,V:Integrator; steps:int):Integrator {.gensym.} =
      `mk s`(T = T, V = V, steps = steps)
    mkInt
  template mkProc2(s:untyped):IntProc =
    proc mkInt(T,V:Integrator; steps:int):Integrator {.gensym.} =
      `mk s`(T = T, V = V, steps = steps, ss[1].parseFloat)
    mkInt
  template mkProc3(s:untyped):IntProc =
    proc mkInt(T,V:Integrator; steps:int):Integrator {.gensym.} =
      `mk s`(T = T, V = V, steps = steps, ss[1].parseFloat, ss[2].parseFloat)
    mkInt
  template mkProc4(s:untyped):IntProc =
    proc mkInt(T,V:Integrator; steps:int):Integrator {.gensym.} =
      `mk s`(T = T, V = V, steps = steps, ss[1].parseFloat, ss[2].parseFloat, ss[3].parseFloat)
    mkInt
  template mkProc5(s:untyped):IntProc =
    proc mkInt(T,V:Integrator; steps:int):Integrator {.gensym.} =
      `mk s`(T = T, V = V, steps = steps, ss[1].parseFloat, ss[2].parseFloat, ss[3].parseFloat, ss[4].parseFloat)
    mkInt
  let ss = s.split(',')
  # Omelyan's triple star integrators, see Omelyan et. al. (2003)
  case ss[0]:
  of "2MN":
    if ss.len == 1: return mkProc1(Omelyan2MN)
    else: return mkProc2(Omelyan2MN)
  of "4MN5FP":
    if ss.len == 1: return mkProc1(Omelyan4MN5FP)
    elif ss.len == 2: return mkProc2(Omelyan4MN5FP)
    elif ss.len == 3: return mkProc3(Omelyan4MN5FP)
    elif ss.len == 4: return mkProc4(Omelyan4MN5FP)
    elif ss.len == 5: return mkProc5(Omelyan4MN5FP)
    else: return mkProc2(Omelyan4MN5FP)
  of "4MN5FV":
    if ss.len == 1: return mkProc1(Omelyan4MN5FV)
    elif ss.len == 2: return mkProc2(Omelyan4MN5FV)
    elif ss.len == 3: return mkProc3(Omelyan4MN5FV)
    elif ss.len == 4: return mkProc4(Omelyan4MN5FV)
    elif ss.len == 5: return mkProc5(Omelyan4MN5FV)
    else: return mkProc2(Omelyan4MN5FV)
  of "6MN7FV": return mkProc1(Omelyan6MN7FV)
  of "4MN3F1GP":  # lambda = 0.2725431326761773  is  FUEL f3g a0=0.109
    if ss.len == 1: return mkProc1(Omelyan4MN3F1GP)
    else: return mkProc2(Omelyan4MN3F1GP)
  of "4MN4F2GVG": return mkProc1(Omelyan4MN4F2GVG)
  of "4MN4F2GV": return mkProc1(Omelyan4MN4F2GV)
  of "4MN5F1GV": return mkProc1(Omelyan4MN5F1GV)
  of "4MN5F1GP": return mkProc1(Omelyan4MN5F1GP)
  of "4MN5F2GV": return mkProc1(Omelyan4MN5F2GV)
  of "4MN5F2GP": return mkProc1(Omelyan4MN5F2GP)
  of "6MN5F3GP": return mkProc1(Omelyan6MN5F3GP)
  else:
    qexError "Cannot parse integrator: '", s, "'\n",
      """Available integrators (with default parameters):
      2MN,0.1931833275037836
      4MN5FP,0.2750081212332419,-0.1347950099106792,-0.08442961950707149,0.3549000571574260
      4MN5FV,0.2539785108410595,-0.03230286765269967,0.08398315262876693,0.6822365335719091
      6MN7FV
      4MN3F1GP,0.2470939580390842
      4MN4F2GVG
      4MN4F2GV
      4MN5F1GV
      4MN5F1GP
      4MN5F2GV
      4MN5F2GP
      6MN5F3GP"""

qexinit()

tic()

letParam:
  gaugefile = ""
  savefile = "config"
  savefreq = 10
  lat =
    if fileExists(gaugefile):
      getFileLattice gaugefile
    else:
      if gaugefile.len > 0:
        qexWarn "Nonexistent gauge file: ", gaugefile
      @[8,8,8,8]
  beta = 6.0
  adjFac = -0.25
  tau = 2.0
  inittraj = 0
  trajs = 10
  seed:uint64 = int(1000*epochTime())
  gintalg:IntProc = "4MN5F2GP"
  gsteps = 4
  mass = @[0.1]  # mass for each staggered species
  hmasses0 = @[0.2,0.4]  # Hasenbusch masses for mass[0]
  hmasses1 = if mass.len>1: hmasses0 else: @[]  # Hasenbusch masses for mass[1]
  hmasses2 = if mass.len>2: hmasses0 else: @[]  # Hasenbusch masses for mass[2]
  hmasses3 = if mass.len>3: hmasses0 else: @[]  # Hasenbusch masses for mass[3]
  hmasses4 = if mass.len>4: hmasses0 else: @[]  # Hasenbusch masses for mass[4]
  fintalg:IntProc = "4MN5F2GP"
  fsteps = repeat(4, mass.len)  # nsteps for each mass
  hfsteps0 = fsteps[0].repeat hmasses0.len  # nsteps for Hasenbusch masses 0
  hfsteps1 = fsteps[0].repeat hmasses1.len  # nsteps for Hasenbusch masses 1
  hfsteps2 = fsteps[0].repeat hmasses2.len  # nsteps for Hasenbusch masses 2
  hfsteps3 = fsteps[0].repeat hmasses3.len  # nsteps for Hasenbusch masses 3
  hfsteps4 = fsteps[0].repeat hmasses4.len  # nsteps for Hasenbusch masses 4
  arsq = 1e-20  # CG r^2 for fermion action
  frsq = repeat(1e-12, mass.len)  # CG r^2 for fermion force for each mass
  hfrsq0 = frsq[0].repeat hmasses0.len  # frsq for Hasenbusch masses 0
  hfrsq1 = frsq[0].repeat hmasses1.len  # frsq for Hasenbusch masses 1
  hfrsq2 = frsq[0].repeat hmasses2.len  # frsq for Hasenbusch masses 2
  hfrsq3 = frsq[0].repeat hmasses3.len  # frsq for Hasenbusch masses 3
  hfrsq4 = frsq[0].repeat hmasses4.len  # frsq for Hasenbusch masses 4
  alwaysAccept:bool = 0
  revCheckFreq = savefreq
  pbpmass = mass
  pbpreps = repeat(1, pbpmass.len)
  pbprsq = arsq
  maxits = 1000000
  useFG2:bool = 0
  hot:bool = 0  # random links instead of unit ones, as staghmc -hot:1
  alpha1 = 0.4  # HYP smearing coefficients
  alpha2 = 0.5
  alpha3 = 0.5
  showTimers:bool = 1
  timerWasteRatio = 0.05
  timerEchoDropped:bool = 0
  timerExpandRatio = 0.05
  verboseGCStats:bool = 0
  verboseTimer:bool = 0
  rg = newSeq[int](0)  # ranks per dimension, lanes as in bestagcg
  mixed:bool = 0  # mixed precision solves
  reals = 0  # reals per link of the solvers, 12, 14 or 18; 0: 14 with smearing, 12 without
  batch:bool = 1  # solve the fermion terms of a force step together
  check:bool = 0  # first compare the smearing and its force with gauge/hypsmear2

installStandardParams()
echoParams()
echo "rank ", myRank, "/", nRanks
threads: echo "thread ", threadNum, "/", numThreads
processHelpParam()

DropWasteTimerRatio = timerWasteRatio
VerboseGCStats = verboseGCStats
VerboseTimer = verboseTimer

if mass.len > 5:
  qexError "Unimlemented for mass: ", mass
if mass.len != fsteps.len or
    mass.len != frsq.len:
  qexError "Parameters for staggered species mismatch."

let
  hmasses = @[hmasses0, hmasses1, hmasses2, hmasses3, hmasses4][0..<mass.len]
  hfsteps = @[hfsteps0, hfsteps1, hfsteps2, hfsteps3, hfsteps4][0..<mass.len]
  hfrsq = @[hfrsq0, hfrsq1, hfrsq2, hfrsq3, hfrsq4][0..<mass.len]

for k in 0..<mass.len:
  if hmasses[k].len != hfsteps[k].len or
      hmasses[k].len != hfrsq[k].len:
    qexError "Hasenbusch parameters lengths mismatch."

if pbpmass.len != pbpreps.len:
  qexError "The lengths of pbpmass and pbpreps differ."

let
  lo = if rg.len > 0: newLayout(lat, VLEN, rg, innerGeom(lat, rg, static(VLEN))) else: lat.newLayout
  gc = GaugeActionCoeffs(plaq: beta, adjplaq: beta*adjFac)
  vol = lo.physVol
  smear = alpha1 != 0 or alpha2 != 0 or alpha3 != 0
  coef = HypCoefs(alpha1: alpha1, alpha2: alpha2, alpha3: alpha3)
echo "smear = ", coef

var r = lo.newRNGFieldV(RngMilc6, seed)
var R: RngMilc6  # global RNG
R.seed(seed, 987654321)

var g = lo.newgauge  # the host copy of the links
var gs = lo.newgauge  # unit links with the phases of the fermions
gs.unit
threads:
  gs.setBC
  threadBarrier()
  gs.stagPhase

# the fermion terms in the order of the integrators: species k, then its
# Hasenbusch masses, term (k, i) solving with mass mt, applying mh
var terms: seq[(int, int)]
for k in 0..<mass.len:
  for i in 0..hmasses[k].len: terms.add (k, i)
let nt = terms.len

proc mt(j: int): float =
  ## the mass of the solves of term j
  let (k, i) = terms[j]
  if i == 0: mass[k] else: hmasses[k][i-1]
proc mh(j: int): float =
  ## the heavier mass of term j, the last term of a species has mt only
  let (k, i) = terms[j]
  if i < hmasses[k].len: hmasses[k][i] else: mt(j)

func sq(x: float): float = x*x

proc fscale(j: int; t: float): float =
  ## as staghmc_sh: the force of term j is d/dU of |M(mt)^-1 M(mh) phi|^2/2
  let (k, i) = terms[j]
  if hmasses[k].len == 0: 0.5*t/mass[k]
  elif i == 0: 0.5*t*(hmasses[k][0].sq - mass[k].sq)/mass[k]
  elif i < hmasses[k].len: 0.5*t*(hmasses[k][i].sq - hmasses[k][i-1].sq)/hmasses[k][i-1]
  else: 0.5*t/hmasses[k][i-1]

var pbpsp = initSolverParams()
pbpsp.r2req = pbprsq
pbpsp.maxits = maxits
var spa, spf = newSeq[SolverParams](nt)  # action and force solves of each term
for j in 0..<nt:
  let (k, i) = terms[j]
  spa[j] = initSolverParams()
  spa[j].r2req = arsq
  spa[j].maxits = maxits
  spa[j].verbosity = 0
  spf[j] = initSolverParams()
  spf[j].r2req = if i == 0: frsq[k] else: hfrsq[k][i-1]
  spf[j].maxits = maxits
  spf[j].verbosity = 0

proc checkStats(label: string; sp: var SolverParams) =
  echo label, sp.getAveStats
  if sp.r2.max > sp.r2req:
    qexError &"Max r2 ({sp.r2.max}) larger than requested ({sp.r2req})"
  sp.resetStats

proc checkSolvers =
  ## as staghmc_sh
  checkStats("Solver[pbp]: ", pbpsp)
  echo "Solver[action]:"
  for j in 0..<nt: checkStats("  A m=" & $mt(j) & " ", spa[j])
  echo "Solver[force]:"
  for j in 0..<nt: checkStats("  F m=" & $mt(j) & " ", spf[j])

var gg = newGpuGauge(lo)
# the smeared links are U(3), projectU without a determinant condition, so
# the solvers keep their determinant; unsmeared links are SU(3) and 12 do
let nl = if reals != 0: reals elif smear: 14 else: 12
if smear and nl < 14: qexError "the smeared links are U(3): -reals:14 or 18"
let useBatch = batch and nt > 1
var s = newStagGpu(gs, float64, nl, batch = useBatch)
var ss: StagGpu[VLEN,float32]
let ssp = if mixed: (ss = newStagGpu(gs, float32, nl, batch = useBatch); addr ss) else: nil
echo "GPU links: ", s.nl, " reals", if s.lb == nil: ", forward only" else: ""
let sg = stagSigns(gs)
let
  mom = gg.newLinks
  g0 = gg.newLinks  # the links at the start of a trajectory
  gfg = gg.newLinks  # the links before a force gradient step
  fbg = gg.newLinks  # the gauge force there
  fbf = gg.newLinks  # the fermion force there
  grev = if revCheckFreq > 0: gg.newLinks else: nil  # links and momenta at the end of a trajectory
  prev = if revCheckFreq > 0: gg.newLinks else: nil
# with smearing the solvers take the smeared links sgo
var sgo: GpuGauge[VLEN]
var hg: HypGpu[VLEN]
if smear:
  sgo = newGpuGauge(lo)
  hg = newHypGpu(lo)
let sgg = if smear: addr sgo else: addr gg  # the links of the solver
let fF = if smear: gg.newLinks else: nil  # sum of the phased one link forces of a step
let n6 = 6*s.n
proc newVec(): ptr UncheckedArray[float] = cast[ptr UncheckedArray[float]](gpuMalloc(n6*sizeof(float)))
var eta, phi, xs: seq[ptr UncheckedArray[float]]  # per term: gaussian draw, pseudofermion, force solution
for j in 0..<nt:
  eta.add newVec()
  phi.add newVec()
  xs.add newVec()
let psi = newVec()
let ftmp = newVec()

proc zero(v: ptr UncheckedArray[float]; n: int) =
  gpuFor(i, n): v[i] = 0.0

var linksSet = false
proc setLinks() =
  ## the smeared, phased links of the solvers from gg.u
  if linksSet: return
  tic()
  if smear:
    hg.smear(coef, gg, sgo.u)
    sgo.fresh = false
  s.setLinks(sgg[], sg)
  if ssp != nil: ssp[].setLinks(sgg[], sg)
  linksSet = true
  toc("smear & rephase")

proc faction(): seq[seq[float]] =
  ## |psi|^2 per species and term, psi = M(mt)^-1 M(mh) phi, the last
  ## M(mt)^-1 phi, as staghmc_sh; psi in xs, M(mh) phi in eta, together
  ## with batch
  tic("faction")
  result = newSeq[seq[float]](mass.len)
  var b = newSeq[ptr UncheckedArray[float]](nt)
  for j in 0..<nt:
    let (k, i) = terms[j]
    if i != hmasses[k].len:
      s.applyMfull(eta[j], phi[j], mh(j))
      b[j] = eta[j]
    else: b[j] = phi[j]
  if useBatch:
    s.solveM(xs, b, toSeq(0..<nt).mapIt(mt(it)), spa, ssp, full = true)
  else:
    for j in 0..<nt: s.solveM(xs[j], b[j], mt(j), spa[j], ssp, full = true)
  toc("solve")
  for j in 0..<nt: result[terms[j][0]].add s.norm2(xs[j])
  toc("norm")

proc gaction(f2: seq[seq[float]]; p2: float): auto =
  let
    ga = gg.actionA(gc)
    fa = f2.mapit(0.5*it)
    t = 0.5*p2 - float(16*vol)
    h = ga + fa.mapit(sum it).sum + t
  (ga, fa, t, h)

proc mdt(t: float) =
  tic()
  gg.expUpdate(mom, t)
  linksSet = false
  toc("mdt")

proc mdv(t: float) =
  tic()
  gg.forceA(gc, mom, t)
  toc("mdv")

proc solveTerms(ids: seq[int]) =
  ## xs_j = M(mt_j)^-1 phi_j on the current links for the terms ids,
  ## together with batch
  setLinks()
  tic()
  if useBatch and ids.len > 1:
    var sp = ids.mapIt(spf[it])
    s.solveM(ids.mapIt(xs[it]), ids.mapIt(phi[it]), ids.mapIt(mt(it)), sp, ssp)
    for k, j in ids: spf[j] = sp[k]
  else:
    for j in ids: s.solveM(xs[j], phi[j], mt(j), spf[j], ssp)
  toc("fforce solve")

proc force(p: ptr UncheckedArray[float]; ids: seq[int]; ts: openArray[float]) =
  ## p += the fermion force of the terms ids with steps ts, from xs on the
  ## current links, as fforce and smearedOneLinkForce of staghmc_sh
  tic()
  if smear:
    zero(fF, 4*18*s.n)
    for j in ids: s.outerM(fF, xs[j], fscale(j, ts[j]))
    toc("outer")
    hg.force(coef, gg, sgo.u, fF, sg, s.ne, p)
    toc("pullback")
  else:
    for j in ids: s.forceM(sgg[], sg, p, xs[j], -fscale(j, ts[j]))
    toc("force")

proc mdvAll(ts, gs: openarray[float]) =
  ## the combined update of mdvAllfga of staghmc_sh: the momenta from the
  ## forces at the links, then the force gradient steps from the forces at
  ## links displaced by the forces there, as approximateFGcoeff; the terms
  ## with a force at the same links solve together
  tic("mdvAll")
  var
    updateG, updateGG = false
    fts, fgts: seq[int]  # the terms with a force, with a force gradient
    ggs: array[2, tuple[t, g: float]]
    fgs: array[2, tuple[t, g: seq[float]]]
  let order = if useFG2: 2 else: 1
  for o in 0..1:
    fgs[o].t = newSeq[float](nt)
    fgs[o].g = newSeq[float](nt)
  if gs[0] != 0:
    updateGG = true
    if ts[0] == 0:
      qexError "Force gradient without the force update."
    if useFG2:
      let (tf, tg) = approximateFGcoeff2(ts[0], gs[0])
      for o in 0..1: ggs[o] = (t: tf[o], g: tg[o])
    else:
      let (tf, tg) = approximateFGcoeff(ts[0], gs[0])
      ggs[0] = (t: tf, g: tg)
  elif ts[0] != 0:
    updateG = true
  for j in 0..<nt:
    if gs[j+1] != 0:
      fgts.add j
      if ts[j+1] == 0:
        qexError "Force gradient without the force update."
      if useFG2:
        let (tf, tg) = approximateFGcoeff2(ts[j+1], gs[j+1])
        for o in 0..1:
          fgs[o].t[j] = tf[o]
          fgs[o].g[j] = tg[o]
      else:
        let (tf, tg) = approximateFGcoeff(ts[j+1], gs[j+1])
        fgs[0].t[j] = tf
        fgs[0].g[j] = tg
    elif ts[j+1] != 0:
      fts.add j
  let fg = updateGG or fgts.len > 0
  if fg: gg.copy(gfg, gg.u)
  if fts.len + fgts.len > 0: solveTerms(sorted(fts & fgts))
  # MD
  if updateG: mdv ts[0]
  if fts.len > 0: force(mom, fts, ts[1..^1])
  # FG: the forces at gfg displace the links, the forces there update the
  # momenta; each force at gfg before any displacement, as staghmc_sh
  # computes them from its copy gg
  if fg:
    if updateGG:
      zero(fbg, 4*18*s.n)
      gg.forceA(gc, fbg, -1.0)  # fbg = F
    for o in 0..<order:
      if fgts.len > 0:
        if o > 0: solveTerms(fgts)
        zero(fbf, 4*18*s.n)
        force(fbf, fgts, fgs[o].g)
      if updateGG: gg.expUpdate(fbg, -ggs[o].g)
      if fgts.len > 0: gg.expUpdate(fbf, 1.0)
      linksSet = false
      if updateGG: mdv ggs[o].t
      if fgts.len > 0:
        solveTerms(fgts)
        force(mom, fgts, fgs[o].t)
      gg.copy(gg.u, gfg)
      gg.fresh = false
      linksSet = false
  toc("done")

let
  (V, T) = newIntegratorPair(mdvAll, mdt)
  H = newParallelEvolution gintalg(T = T, V = V[0], steps = gsteps)
block:
  var j = 0
  for k in 0..<mass.len:
    inc j
    H.add fintalg(T = T, V = V[j], steps = fsteps[k])
    for i in 0..<hfsteps[k].len:
      inc j
      H.add fintalg(T = T, V = V[j], steps = hfsteps[k][i])

proc revCheck(evo: auto; h0, ga0, t0: float; fa0: seq[seq[float]]) =
  ## as staghmc_sh: evolve back from the end with -p and print H
  tic("reversibility")
  gg.copy(grev, gg.u)
  gg.copy(prev, mom)
  let pm = mom
  gpuFor(i, 4*18*s.n): pm[i] = -pm[i]
  evo.evolve tau
  evo.finish
  let p2 = gg.norm2(mom)
  toc("p norm2 2")
  setLinks()
  let f2 = faction()
  toc("fa solve 2")
  let (ga1, fa1, t1, h1) = gaction(f2, p2)
  var dsf = newseq[seq[float]](fa1.len)
  for k in 0..<fa1.len:
    dsf[k] = newseq[float](fa1[k].len)
    for i in 0..<fa1[k].len:
      dsf[k][i] = fa1[k][i] - fa0[k][i]
  qexLog "Reversed H: ",h1,"  Sg: ",ga1,"  Sf: ",fa1,"  T: ",t1
  echo "Reversibility: dH: ",h1-h0,"  dSg: ",ga1-ga0,"  dSf: ",dsf,"  dT: ",t1-t0
  gg.copy(gg.u, grev)
  gg.copy(mom, prev)
  gg.fresh = false
  linksSet = false
  toc("done")

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

proc mplaq() =
  ## as mplaq of staghmc_sh, on the GPU
  tic()
  let pl = gg.plaqs
  echo "MEASplaq ss: ",pl[0],"  st: ",pl[1],"  tot: ",0.5*(pl[0]+pl[1])
  toc("plaq")

var gw = newGpuGauge(lo)  # scratch for the Polyakov loop
proc ploop() =
  ## as ploop of staghmc_sh, on the GPU
  tic()
  let pl = gg.ploop(gw, gfg)
  let pls = (re: (pl[0].re + pl[1].re + pl[2].re)/3.0, im: (pl[0].im + pl[1].im + pl[2].im)/3.0)
  echo "MEASploop spatial: ",pls.re," ",pls.im," temporal: ",pl[3].re," ",pl[3].im
  toc("ploop")

let rgpu = newRngGpu(r)  # after the host draws of the links
proc pbp() =
  ## as staghmc_sh: m |M(m)^-1 u1|^2/vol for each mass, on the current links
  tic()
  setLinks()
  for k in 0..<pbpmass.len:
    let m = pbpmass[k]
    for i in 0..<pbpreps[k]:
      rgpu.u1(ftmp, 6)
      s.solveM(psi, ftmp, m, pbpsp, ssp, full = true)
      let pbp = s.norm2(psi)
      echo "MEASpbp mass ",m," : ",m*pbp/vol.float
  toc("pbp")

proc rel(a, b: auto): float =
  ## |a-b|^2/|b|^2 summed over the fields
  var d, n = 0.0
  threads:
    var dt, nt = 0.0
    for mu in 0..<a.len:
      dt += norm2(a[mu] - b[mu])
      nt += b[mu].norm2
    threadMaster:
      d = dt
      n = nt
  d/n

if check and smear:
  # random links and chain: the smeared links and p = TAH(F U^+) for the
  # force F of the phased chain against gauge/hypsmear2 and contractProjectTAH
  var rc = lo.newRNGField(RngMilc6, seed + 1)
  var gr = lo.newgauge
  var ch = lo.newgauge
  var cp = lo.newgauge
  var fl = lo.newgauge
  var f = lo.newgauge
  var fg = lo.newgauge
  var z = lo.newgauge
  gr.random rc
  ch.gaussian rc
  threads:
    for mu in 0..<cp.len:
      cp[mu] := ch[mu]
      z[mu] := 0
    cp.setBC
    threadBarrier()
    cp.stagPhase
    threadBarrier()
    for mu in 0..<cp.len:
      for i in cp[mu].odd:
        cp[mu][i] *= -1
  let htc = newHypTemps(gr)
  htc.smear(coef, fl)
  htc.force(coef, f, cp)
  threads: contractProjectTAH(f, f, gr)
  gg.upload(gg.u, gr)
  hg.smear(coef, gg, sgo.u)
  gg.download(fg, sgo.u)
  echo "check smear  |GPU-CPU|^2/|CPU|^2: ", rel(fg, fl)
  gg.upload(fF, ch)
  gg.upload(mom, z)
  hg.force(coef, gg, sgo.u, fF, sg, s.ne, mom)
  gg.download(fg, mom)
  echo "check force  |GPU-CPU|^2/|CPU|^2: ", rel(fg, f)

if fileExists(gaugefile):
  tic("load")
  if 0 != g.loadGauge gaugefile:
    qexError "failed to load gauge file: ", gaugefile
  qexLog "loaded gauge from file: ", gaugefile," secs: ",getElapsedTime()
  toc("read")
  g.reunit
  toc("reunit")
elif hot:
  g.random r
else:
  g.unit
rgpu.upload(r)  # the host generators after the draws of the links
gg.upload(gg.u, g)
gg.fresh = false

mplaq()

echo H

toc("prep")

for n in inittraj+1..inittraj+trajs:
  tic("traj")
  let tt = epochTime()
  rgpu.randomTAH mom
  gg.copy(g0, gg.u)
  toc("p refresh, save g")
  let p2 = gg.norm2(mom)
  toc("p norm2 1")
  setLinks()
  # the gaussian draws in the order of staghmc_sh (bsm.lua): term i of each
  # species, then term i+1
  block:
    var i = 0
    var running = true
    while running:
      running = false
      var j0 = 0
      for k in 0..<mass.len:
        if i <= hmasses[k].len:
          rgpu.gaussian(eta[j0+i], 6)
          running = true
        j0 += hmasses[k].len + 1
      inc i
  # phi = M(-mh)^-1 M(-mt) eta on the even sites, the last M(-mt) eta, the
  # minus signs as staghmc_sh (bsm.lua convention); M(-mt) eta in xs
  var hb: seq[int]  # the terms with a solve
  for j in 0..<nt:
    let (k, i) = terms[j]
    if i != hmasses[k].len:
      s.applyMfull(xs[j], eta[j], -mt(j))
      hb.add j
    else:
      s.applyMfull(phi[j], eta[j], -mt(j))
  if useBatch and hb.len > 1:
    var sp = hb.mapIt(spa[it+1])
    s.solveM(hb.mapIt(phi[it]), hb.mapIt(xs[it]), hb.mapIt(-mh(it)), sp, ssp, full = true)
    for k, j in hb: spa[j+1] = sp[k]
  else:
    for j in hb: s.solveM(phi[j], xs[j], -mh(j), spa[j+1], ssp, full = true)
  for j in 0..<nt: zero(cast[ptr UncheckedArray[float]](addr phi[j][6*s.ne]), n6 - 6*s.ne)
  toc("init")
  var f2 = faction()
  toc("fa solve 1")
  let (ga0, fa0, t0, h0) = gaction(f2, p2)
  toc("init gauge action")
  qexLog "Begin H: ",h0,"  Sg: ",ga0,"  Sf: ",fa0,"  T: ",t0

  H.evolve tau
  H.finish
  toc("evolve")

  let p2e = gg.norm2(mom)
  toc("p norm2 2")
  setLinks()
  f2 = faction()
  toc("fa solve 2")
  let (ga1, fa1, t1, h1) = gaction(f2, p2e)
  toc("final gauge action")
  qexLog "End H: ",h1,"  Sg: ",ga1,"  Sf: ",fa1,"  T: ",t1
  toc("end evolve")

  if revCheckFreq > 0 and n mod revCheckFreq == 0:
    H.revCheck(h0, ga0, t0, fa0)

  let
    dH = h1 - h0
    acc = exp(-dH)
    accr = R.uniform
  if accr <= acc or alwaysAccept:  # accept
    echo "ACCEPT:  dH: ",dH,"  exp(-dH): ",acc,"  r: ",accr,(if alwaysAccept:" (ignored)" else:"")
    tic()
    let du = gg.reunit
    echo "unitary deviation avg: ",du[0].avg," max: ",du[0].max
    echo "new unitary deviation avg: ",du[1].avg," max: ",du[1].max
    toc("reunit")
  else:  # reject
    echo "REJECT:  dH: ",dH,"  exp(-dH): ",acc,"  r: ",accr
    gg.copy(gg.u, g0)
  gg.fresh = false
  linksSet = false
  pbp()

  mplaq()
  ploop()
  toc("measure")

  if savefreq > 0 and n mod savefreq == 0:
    tic("save")
    let fn = savefile & &".{n:05}.lime"
    gg.download(g, gg.u)
    if 0 != g.saveGauge(fn):
      qexError "Failed to save gauge to file: ",fn
    qexLog "saved gauge to file: ",fn," secs: ",getElapsedTime()
    toc("done")

  checkSolvers()

  echo "traj secs: ", epochTime() - tt
  qexLog "traj ",n," secs: ",getElapsedTime()
  toc("traj end")

toc("hmc")

echoProf()
processSaveParams()
writeParamFile()
qexfinalize()
