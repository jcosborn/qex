## staghmc_s.nim (hypsmear) with Hasenbusch masses.

import qex, gauge, gauge/hypsmear, physics/qcdTypes, physics/stagSolve
import mdevolve
import math, os, sequtils, strutils, times

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
  showTimers:bool = 1
  timerWasteRatio = 0.05
  timerEchoDropped:bool = 0
  timerExpandRatio = 0.05
  verboseGCStats:bool = 0
  verboseTimer:bool = 0

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
  lo = lat.newLayout
  #gc = GaugeActionCoeffs(plaq:6)
  gc = GaugeActionCoeffs(plaq: beta, adjplaq: beta*adjFac)
  vol = lo.physVol

var r = lo.newRNGField(RngMilc6, seed)
var R:RngMilc6  # global RNG
R.seed(seed, 987654321)

var
  g = lo.newgauge
  p = lo.newgauge
  f = lo.newgauge
  g0 = lo.newgauge
  sg0 = lo.newgauge
  sg = lo.newGauge  # smeared gauge
  gg = lo.newgauge  # FG backup gauge
  sgg = lo.newgauge  # FG backup smeared gauge
let
  stag = newStag(sg)
  stag0 = newStag(sg0)
  stagg = newStag(sgg)

var
  info: PerfInfo
  coef = HypCoefs(alpha1:0.4, alpha2:0.5, alpha3:0.5)
echo "smear = ",coef

var ftmp = lo.ColorVector()
type CV = typeof(ftmp)
var
  phi = newseq[seq[CV]](mass.len)
  psi = newseq[seq[CV]](mass.len)
for k in 0..<mass.len:
  phi[k] = newseq[typeof(ftmp)](hmasses[k].len+1)
  psi[k] = newseq[typeof(ftmp)](hmasses[k].len+1)
  for i in 0..<phi[k].len:
    phi[k][i] = lo.ColorVector()
    psi[k][i] = lo.ColorVector()

var pbpsp = initSolverParams()
pbpsp.r2req = pbprsq
pbpsp.maxits = maxits

var spa = newseq[typeof(pbpsp)](mass.len)
for k in 0..<mass.len:
  spa[k] = initSolverParams()
  #spa[k].subsetName = "even"
  spa[k].r2req = arsq
  spa[k].maxits = maxits
  spa[k].verbosity = 0

var spah = newseq[typeof(spa)](mass.len)
for k in 0..<mass.len:
  spah[k] = newseq[typeof(pbpsp)](hmasses[k].len)
  for i in 0..<spah[k].len:
    spah[k][i] = initSolverParams()
    spah[k][i].r2req = arsq
    spah[k][i].maxits = maxits
    spah[k][i].verbosity = 0

var spf = newseq[typeof(pbpsp)](mass.len)
for k in 0..<mass.len:
  spf[k] = initSolverParams()
  #spf[k].subsetName = "even"
  spf[k].r2req = frsq[k]
  spf[k].maxits = maxits
  spf[k].verbosity = 0

var spfh = newseq[typeof(spf)](mass.len)
for k in 0..<mass.len:
  spfh[k] = newseq[typeof(pbpsp)](hmasses[k].len)
  for i in 0..<spfh[k].len:
    spfh[k][i] = initSolverParams()
    spfh[k][i].r2req = hfrsq[k][i]
    spfh[k][i].maxits = maxits
    spfh[k][i].verbosity = 0

proc checkStats(label:string, sp:var SolverParams) =
  echo label,sp.getAveStats
  if sp.r2.max > sp.r2req:
    qexError &"Max r2 ({sp.r2.max}) larger than requested ({sp.r2req})"
  sp.resetStats

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

proc pbp(stag:auto) =
  tic()
  var ftmp2 = lo.ColorVector()
  for k in 0..<pbpmass.len:
    let m = pbpmass[k]
    for i in 0..<pbpreps[k]:
      threads:
        ftmp.u1 r
      stag.solve(ftmp2, ftmp, m, pbpsp)
      threads:
        var pbp = ftmp2.norm2
        threadMaster:
          echo "MEASpbp mass ",m," : ",m*pbp/vol.float
  toc("pbp")

proc mplaq(g:auto) =
  tic()
  let
    pl = g.plaq
    nl = pl.len div 2
    ps = pl[0..<nl].sum * 2.0
    pt = pl[nl..^1].sum * 2.0
  echo "MEASplaq ss: ",ps,"  st: ",pt,"  tot: ",0.5*(ps+pt)
  toc("plaq")

proc ploop(g:auto) =
  tic()
  let pg = g[0].l.physGeom
  var pl = newseq[typeof(g.wline @[1])](pg.len)
  for i in 0..<pg.len:
    pl[i] = g.wline repeat(i+1, pg[i])
  let
    pls = pl[0..^2].sum / float(pl.len-1)
    plt = pl[^1]
  echo "MEASploop spatial: ",pls.re," ",pls.im," temporal: ",plt.re," ",plt.im
  toc("ploop")

proc fgsave =
  threads:
    for mu in 0..<g.len:
      gg[mu] := g[mu]
proc fgload =
  threads:
    for mu in 0..<g.len:
      g[mu] := gg[mu]

proc smearRephase(g: auto, sg: auto):auto =
  tic()
  let smearedForce = coef.smearGetForce(g, sg, info)
  toc("smear")
  threads:
    sg.setBC
    threadBarrier()
    sg.stagPhase
  toc("BC & Phase")
  smearedForce

proc smearRephaseDiscardForce(g: auto, sg: auto) =
  tic()
  coef.smear(g, sg, info)
  qexGC "smear done"
  toc("smear w/o force")
  threads:
    sg.setBC
    threadBarrier()
    sg.stagPhase
  toc("BC & Phase")

template pnorm2(p2:float) =
  threads:
    var p2t = 0.0
    for i in 0..<p.len:
      p2t += p[i].norm2
    threadMaster: p2 = p2t

proc gaction(g:auto, f2:seq[seq[float]], p2:float):auto =
  tic()
  let
    ga = gc.actionA g
    fa = f2.mapit(0.5*it)
    t = 0.5*p2 - float(16*vol)
    h = ga + fa.mapit(sum it).sum + t
  toc("gaction")
  (ga, fa, t, h)

template faction(fa:seq[seq[float]]) =
  for k in 0..<fa.len:
    tic("faction")
    for i in 0..<phi[k].len-1:
      threads:
        stag.D(ftmp, phi[k][i], hmasses[k][i])
      if i == 0:
        tic()
        stag.solve(psi[k][i], ftmp, mass[k], spa[k])
        toc("asolve " & $mass[k])
      else:
        tic()
        stag.solve(psi[k][i], ftmp, hmasses[k][i-1], spah[k][i-1])
        toc("asolve " & $hmasses[k][i-1])
    if hmasses[k].len>0:
      tic()
      stag.solve(psi[k][^1], phi[k][^1], hmasses[k][^1], spah[k][^1])
      toc("asolve " & $hmasses[k][^1])
    else:
      tic()
      stag.solve(psi[k][^1], phi[k][^1], mass[k], spa[k])
      toc("asolve " & $mass[k])
    threads:
      for i in 0..<psi[k].len:
        var psi2 = psi[k][i].norm2
        threadMaster: fa[k][i] = psi2
    toc("done")

proc massIndex(i:int):(int,int) =
  var
   i = i
   k = 0
  for j in 0..<mass.len:
    let h = hmasses[j].len + 1
    if i >= h:
      i -= h
    else:
      k = j
      break
  (k,i)

func sq(x:float):float = x*x

proc fscale(k,i:int, t:float):float =
  if hmasses[k].len == 0: 0.5*t/mass[k]
  elif i == 0: 0.5*t*(hmasses[k][0].sq-mass[k].sq)/mass[k]
  elif i < hmasses[k].len: 0.5*t*(hmasses[k][i].sq-hmasses[k][i-1].sq)/hmasses[k][i-1]
  else: 0.5*t/hmasses[k][i-1]

proc smearedOneLinkForce(f: auto, smearedForce: proc, g:auto) =
  tic("olf")
  # Reverse accumulation of the smearing derivatives
  # 1. correcting phase
  threads:
    f.setBC
    threadBarrier()
    f.stagPhase
    threadBarrier()
    for mu in 0..<f.len:
      for i in f[mu].odd:
        f[mu][i] *= -1
  toc("phase")

  # 2. smearing
  f.smearedForce f
  toc("smear")

  # 3. Tₐ ReTr( Tₐ U F† )
  threads:
    for mu in 0..<f.len:
      for i in f[mu]:
        var s {.noinit.}: typeof(f[0][0])
        s := f[mu][i] * g[mu][i].adj
        projectTAH(f[mu][i], s)
  qexGC "combine"
  toc("combine")

proc fforce(stag: auto, f: auto, sf: proc, g: auto, ix:openarray[int], ts:openarray[float]) =
  tic("fforce")
  var t: array[4,Shifter[typeof(ftmp), typeof(ftmp[0])]]
  for mu in 0..<f.len:
    t[mu] = newShifter(ftmp, mu, 1)
  var first = true
  for j in ix:
    tic("floop")
    let (k,i) = massIndex j
    if i == 0:
      tic()
      stag.solve(ftmp, phi[k][i], mass[k], spf[k])
      toc("fsolve " & $mass[k])
    else:
      tic()
      stag.solve(ftmp, phi[k][i], hmasses[k][i-1], spfh[k][i-1])
      toc("fsolve " & $hmasses[k][i-1])
    qexGC "solve"
    toc("solve")

    for mu in 0..<f.len:
      discard t[mu] ^* ftmp
    const n = ftmp[0].len
    let s = fscale(k, i, ts[j])
    if first:
      first = false
      threads:
        for mu in 0..<f.len:
            for i in f[mu]:
              forO a, 0, n-1:
                forO b, 0, n-1:
                  f[mu][i][a,b] := s * ftmp[i][a] * t[mu].field[i][b].adj
    else:
      threads:
        for mu in 0..<f.len:
          for i in f[mu]:
            forO a, 0, n-1:
              forO b, 0, n-1:
                f[mu][i][a,b] += s * ftmp[i][a] * t[mu].field[i][b].adj
    toc("outer")
  toc("solves")
  f.smearedOneLinkForce(sf, g)
  toc("olf")

proc mdt(t: float) =
  tic()
  threads:
    for mu in 0..<g.len:
      for s in g[mu]:
        g[mu][s] := exp(t*p[mu][s])*g[mu][s]
  toc("mdt")
proc mdv(t: float) =
  tic()
  gc.forceA(g, f)
  qexGC "mdv forceA"
  threads:
    for mu in 0..<f.len:
      p[mu] -= t*f[mu]
  toc("mdv")

proc mdvf(ix:openarray[int], sf:proc, ts:openarray[float]) =
  tic()
  stag.fforce(f, sf, g, ix, ts)
  qexGC "mdvf fforce"
  threads:
    for mu in 0..<f.len:
      p[mu] += f[mu]
  toc("mdvf")

# FG update g from backup, gg and sgg
proc fgv(t: float) =
  tic()
  gc.forceA(gg, f)
  qexGC "fgv forceA"
  threads:
    for mu in 0..<g.len:
      for s in g[mu]:
        g[mu][s] := exp((-t)*f[mu][s])*g[mu][s]
  toc("fgv")
proc fgvf(ix:openarray[int], sf:proc, ts:openarray[float]) =
  tic()
  stagg.fforce(f, sf, gg, ix, ts)
  qexGC "fgvf fforce"
  threads:
    for mu in 0..<g.len:
      for s in g[mu]:
        g[mu][s] := exp(f[mu][s])*g[mu][s]
  toc("fgvf")

# Combined update for sharing computations
proc mdvAllfga(ts,gs:openarray[float]) =
  tic("mdvAllfga")
  var
    sforceShared:typeof(g.smearRephase sg) = nil
    updateG = false
    updateGG = false
    updateF = newseq[int](0)
    updateFG = newseq[int](0)
  let approxOrder = if useFG2: 2 else: 1
  type
    GGstep = tuple[t,g:float]
    FGstep = tuple[t,g:seq[float]]
  var
    ggs: array[2,GGstep]
    fgs: array[2,FGstep]  # 2nd item for 2nd order approximation
  for o in 0..1:
    fgs[o].t = newseq[float](gs.len-1)
    fgs[o].g = newseq[float](gs.len-1)

  # gauge
  if gs[0] != 0:
    updateGG = true
    if ts[0] == 0:
      qexError "Force gradient without the force update."
    if useFG2:
      let (tf,tg) = approximateFGcoeff2(ts[0],gs[0])
      for o in 0..1:
        ggs[o] = (t:tf[o], g:tg[o])
    else:
      let (tf,tg) = approximateFGcoeff(ts[0],gs[0])
      ggs[0] = (t:tf, g:tg)
  elif ts[0] != 0:
    updateG = true

  # fermions
  for k in 0..<gs.len-1:
    let i = k+1
    if gs[i] != 0:
      updateFG.add k
      if ts[i] == 0:
        qexError "Force gradient without the force update."
      if useFG2:
        let (tf,tg) = approximateFGcoeff2(ts[i],gs[i])
        for o in 0..1:
          fgs[o].t[k] = tf[o]
          fgs[o].g[k] = tg[o]
      else:
        let (tf,tg) = approximateFGcoeff(ts[i],gs[i])
        fgs[0].t[k] = tf
        fgs[0].g[k] = tg
    elif ts[i] != 0:
      updateF.add k

  # Note that the smeared force proc retains a reference to the input
  # gauge field.  In order to reuse the force proc, we need to use the
  # correct gauge field that would remain the same.

  if updateGG or updateFG.len > 0:
    tic()
    fgsave()
    toc("FG save")
    if updateFG.len > 0:  # fermions in FG
      tic("FG prep")
      sforceShared = gg.smearRephase sgg
      toc("FG smear rephase")
      if updateF.len > 0:  # reuse in mdvf requires sg for stag
        tic("FG copy smeared")
        threads:
          for mu in 0..<sg.len:
            sg[mu] := sgg[mu]
        toc("done")
      toc("done")
    toc("done")

  # MD
  if updateG:
    tic("MD mdv")
    mdv ts[0]
    toc("done")
  if updateF.len > 0:
    tic("MD mdvf")
    if updateFG.len == 0:
      tic()
      var sforcef = g.smearRephase sg
      toc("MD smear rephase")
      mdvf(updateF, sforcef, ts[1..^1])
      sforcef = nil
      qexGC "MD mdvf w/smear done"
    else:
      mdvf(updateF, sforceShared, ts[1..^1])
      qexGC "MD mdvf done"
    toc("done")

  # FG
  if updateGG or updateFG.len > 0:
    tic("updateFG")
    for o in 0..<approxOrder:
      tic()
      if updateGG:
        tic("fgv")
        fgv ggs[o].g
        toc("done")
      if updateFG.len > 0:
        tic("fgvf")
        fgvf(updateFG, sforceShared, fgs[o].g)
        if o+1 == approxOrder: sforceShared = nil
        qexGC "fgvf done"
        toc("done")
      toc("FG")
      if updateGG:
        tic("FG mdv")
        mdv ggs[o].t
        toc("done")
      if updateFG.len > 0:
        tic("FG mdvf")
        var sforceg = g.smearRephase sg
        toc("FG smear rephase temp")
        mdvf(updateFG, sforceg, fgs[o].t)
        sforceg = nil
        qexGC "FG mdvf done"
        toc("done")
      toc("FG MD")
      fgload()
      toc("load")
    toc("done")
  toc("done")

proc revCheck(evo:auto; h0,ga0,t0:float, fa0:seq[seq[float]]) =
  tic("reversibility")
  var
    g1 = lo.newgauge
    p1 = lo.newgauge
    p2 = 0.0
    f2 = newseq[seq[float]](phi.len)
  for k in 0..<phi.len: f2[k] = newseq[float](phi[k].len)
  threads:
    for i in 0..<g1.len:
      g1[i] := g[i]
      p1[i] := p[i]
      p[i] := -1*p[i]
  evo.evolve tau
  evo.finish
  p2.pnorm2
  toc("p norm2 2")
  g.smearRephaseDiscardForce sg
  toc("smear & rephase 2")
  f2.faction
  toc("fa solve 2")
  let (ga1,fa1,t1,h1) = g.gaction(f2,p2)
  var dsf = newseq[seq[float]](fa1.len)
  for k in 0..<fa1.len:
    dsf[k] = newseq[float](fa1[k].len)
    for i in 0..<fa1[k].len:
      dsf[k][i] = fa1[k][i] - fa0[k][i]
  qexLog "Reversed H: ",h1,"  Sg: ",ga1,"  Sf: ",fa1,"  T: ",t1
  echo "Reversibility: dH: ",h1-h0,"  dSg: ",ga1-ga0,"  dSf: ",dsf,"  dT: ",t1-t0
  for i in 0..<g1.len:
    g[i] := g1[i]
    p[i] := p1[i]
  qexGC "revCheck done"
  toc("done")

proc checkSolvers =
  checkStats("Solver[pbp]: ", pbpsp)
  echo "Solver[action]:"
  for k in 0..<spa.len:
    checkStats("  A m=" & $mass[k] & " ", spa[k])
    for i in 0..<spah[k].len:
      checkStats("  A m=" & $hmasses[k][i] & " ", spah[k][i])
  echo "Solver[force]:"
  for k in 0..<spf.len:
    checkStats("  F m=" & $mass[k] & " ", spf[k])
    for i in 0..<spfh[k].len:
      checkStats("  F m=" & $hmasses[k][i] & " ", spfh[k][i])

let
  (V,T) = newIntegratorPair(mdvAllfga, mdt)
  H = newParallelEvolution gintalg(T = T, V = V[0], steps = gsteps)
block:
  var j = 0
  for k in 0..<mass.len:
    inc j
    H.add fintalg(T = T, V = V[j], steps = fsteps[k])
    for i in 0..<hfsteps[k].len:
      inc j
      H.add fintalg(T = T, V = V[j], steps = hfsteps[k][i])

if fileExists(gaugefile):
  tic("load")
  if 0 != g.loadGauge gaugefile:
    qexError "failed to load gauge file: ", gaugefile
  qexLog "loaded gauge from file: ", gaugefile," secs: ",getElapsedTime()
  toc("read")
  g.reunit
  toc("reunit")
else:
  #g.random r
  g.unit

g.mplaq

echo H

toc("prep")

for n in inittraj+1..inittraj+trajs:
  tic("traj")
  var p2 = 0.0
  var f2 = newseq[seq[float]](phi.len)
  for k in 0..<phi.len: f2[k] = newseq[float](phi[k].len)
  threads:
    p.randomTAH r
    for i in 0..<g.len:
      g0[i] := g[i]
  toc("p refresh, save g")
  p2.pnorm2
  toc("p norm2 1")
  g.smearRephaseDiscardForce sg
  toc("smear & rephase 1")
  threads:
    for i in 0..<sg.len:
      sg0[i] := sg[i]
  # conforms to the initialization order in bsm.lua
  threads:
    var i = 0
    var running = true
    while running:
      running = false
      for k in 0..<psi.len:
        if i >= psi[k].len: continue
        psi[k][i].gaussian r
        running = true
      inc i
  # phi = D(m2)^{-1} D(m1) psi
  for k in 0..<phi.len:
    for i in 0..<phi[k].len:
      threads:
        if i != phi[k].len-1:
          stag.D(ftmp, psi[k][i], if i==0: -mass[k] else: -hmasses[k][i-1])  # `-` for bsm.lua convention giving -D⁺
        else:
          stag.D(phi[k][i], psi[k][i], if hmasses[k].len>0: -hmasses[k][i-1] else: -mass[k])
      if i != phi[k].len-1:
        stag.solve(phi[k][i], ftmp, -hmasses[k][i], spah[k][i])
      threads:
        phi[k][i].odd := 0
  toc("init")
  f2.faction
  toc("fa solve 1")
  let (ga0,fa0,t0,h0) = g0.gaction(f2,p2)
  qexGC "init action"
  toc("init gauge action")
  qexLog "Begin H: ",h0,"  Sg: ",ga0,"  Sf: ",fa0,"  T: ",t0

  H.evolve tau
  H.finish
  toc("evolve")

  p2.pnorm2
  toc("p norm2 2")
  g.smearRephaseDiscardForce sg
  toc("smear & rephase 2")
  f2.faction
  toc("fa solve 2")
  let (ga1,fa1,t1,h1) = g.gaction(f2,p2)
  qexGC "end action"
  toc("final gauge action")
  qexLog "End H: ",h1,"  Sg: ",ga1,"  Sf: ",fa1,"  T: ",t1
  toc("end evolve")

  if revCheckFreq > 0 and n mod revCheckFreq == 0:
    H.revCheck(h0,ga0,t0,fa0)

  let
    dH = h1 - h0
    acc = exp(-dH)
    accr = R.uniform
  if accr <= acc or alwaysAccept:  # accept
    echo "ACCEPT:  dH: ",dH,"  exp(-dH): ",acc,"  r: ",accr,(if alwaysAccept:" (ignored)" else:"")
    g.reunit
    g.smearRephaseDiscardForce sg
    stag.pbp
  else:  # reject
    echo "REJECT:  dH: ",dH,"  exp(-dH): ",acc,"  r: ",accr
    threads:
      for i in 0..<g.len:
        g[i] := g0[i]
    stag0.pbp
  qexGC "traj done"

  g.mplaq
  g.ploop
  qexGC "measure done"
  toc("measure")

  if savefreq > 0 and n mod savefreq == 0:
    tic("save")
    let fn = savefile & &".{n:05}.lime"
    if 0 != g.saveGauge(fn):
      qexError "Failed to save gauge to file: ",fn
    qexLog "saved gauge to file: ",fn," secs: ",getElapsedTime()
    toc("done")

  checkSolvers()

  qexLog "traj ",n," secs: ",getElapsedTime()
  toc("traj end")

toc("hmc")

#if showTimers: echoTimers(timerExpandRatio, timerEchoDropped)
echoProf()
processSaveParams()
writeParamFile()
qexfinalize()
