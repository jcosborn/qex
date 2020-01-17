## staghmc_s.nim (hypsmear) with Hasenbusch masses.

import qex, gauge, gauge/hypsmear, physics/qcdTypes, physics/stagSolve
import mdevolve
import math, sequtils, strutils, times

const ReversibilityCheck {.booldefine.} = false

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
    echo "Error: cannot parse integrator: '", ss, "'"
    echo """Available integrators (with default parameters):
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
    qexAbort()

qexinit()

letParam:
  lat = @[8,8,8,8]
  beta = 6.0
  adjFac = -0.25
  tau = 2.0
  gsteps = 4
  fsteps = 4
  trajs = 10
  seed:uint64 = int(1000*epochTime())
  mass = 0.1
  hmasses = @[0.2,0.4]  # Hasenbusch masses
  hfsteps = @[fsteps,fsteps]  # nsteps for Hasenbusch masses
  gintalg:IntProc = "4MN5F2GP"
  fintalg:IntProc = "4MN5F2GP"
  useFG2:bool = 0
  arsq = 1e-20
  frsq = 1e-12
  hfrsq = @[frsq,frsq]
  maxits = 10000

echoParams()

if hmasses.len != hfsteps.len or
    hmasses.len != hfrsq.len:
  echo "Error: Hasenbusch parameters lengths mismatch."
  qexAbort()

let
  lo = lat.newLayout
  #gc = GaugeActionCoeffs(plaq:6)
  gc = GaugeActionCoeffs(plaq: beta, adjplaq: beta*adjFac)
  vol = lo.physVol

var r = lo.newRNGField(RngMilc6, seed)
var R:RngMilc6  # global RNG
R.seed(seed, 987654321)

var g = lo.newgauge
#g.random r
g.unit

echo 6.0*g.plaq
echo gc.actionA g

var
  p = lo.newgauge
  f = lo.newgauge
  g0 = lo.newgauge

var
  ftmp = lo.ColorVector()
  phi = newseq[typeof(lo.ColorVector())](hmasses.len+1)
  psi = newseq[typeof(lo.ColorVector())](hmasses.len+1)
for i in 0..<phi.len:
  phi[i] = lo.ColorVector()
  psi[i] = lo.ColorVector()

var spa = initSolverParams()
#spa.subsetName = "even"
spa.r2req = arsq
spa.maxits = maxits
var spf = initSolverParams()
#spf.subsetName = "even"
spf.r2req = frsq
spf.maxits = maxits
spf.verbosity = 0
var spfh = newseq[typeof(spf)](hfrsq.len)
for i in 0..<spfh.len:
  spfh[i] = initSolverParams()
  spfh[i].r2req = hfrsq[i]
  spfh[i].maxits = maxits
  spfh[i].verbosity = 0

var
  info: PerfInfo
  coef = HypCoefs(alpha1:0.4, alpha2:0.5, alpha3:0.5)
echo "smear = ",coef
var sg = lo.newGauge
let stag = newStag(sg)

proc smearRephase(g: any, sg: any):auto {.discardable.} =
  tic()
  let smearedForce = coef.smear(g, sg, info)
  toc("smear")
  threads:
    sg.setBC
    threadBarrier()
    sg.stagPhase
  toc("BC & Phase")
  smearedForce

template pnorm2(p2:float) =
  threads:
    var p2t = 0.0
    for i in 0..<p.len:
      p2t += p[i].norm2
    threadMaster: p2 = p2t

proc gaction(g:any, f2:seq[float], p2:float):auto =
  let
    ga = gc.actionA g
    fa = f2.mapit(0.5*it)
    t = 0.5*p2 - float(16*vol)
    h = ga + fa.sum + t
  (ga, fa, t, h)

template faction(fa:seq[float]) =
  for i in 0..<phi.len-1:
    threads:
      stag.D(ftmp, phi[i], hmasses[i])
    stag.solve(psi[i], ftmp, if i==0: mass else: hmasses[i-1], spa)
  stag.solve(psi[^1], phi[^1], hmasses[^1], spa)
  threads:
    for i in 0..<psi.len:
      var psi2 = psi[i].norm2()
      threadMaster: fa[i] = psi2

proc smearedOneLinkForce(f: any, smearedForce: proc, p: any, g:any) =
  # reverse accumulation of the derivative
  # 1. Dslash
  var t: array[4,Shifter[typeof(p), typeof(p[0])]]
  for mu in 0..<f.len:
    t[mu] = newShifter(p, mu, 1)
    discard t[mu] ^* p
  const n = p[0].len
  threads:
    for mu in 0..<f.len:
      for i in f[mu]:
        forO a, 0, n-1:
          forO b, 0, n-1:
            f[mu][i][a,b] := p[i][a] * t[mu].field[i][b].adj

  # 2. correcting phase
  threads:
    f.setBC
    threadBarrier()
    f.stagPhase
    threadBarrier()
    for mu in 0..<f.len:
      for i in f[mu].odd:
        f[mu][i] *= -1

  # 3. smearing
  f.smearedForce f

  # 4. Tₐ ReTr( Tₐ U F† )
  threads:
    for mu in 0..<f.len:
      for i in f[mu]:
        var s {.noinit.}: typeof(f[0][0])
        s := f[mu][i] * g[mu][i].adj
        projectTAH(f[mu][i], s)


proc fforce(f: any, sf: proc, i: int) =
  tic()
  if i == 0:
    stag.solve(ftmp, phi[i], mass, spf)
  else:
    stag.solve(ftmp, phi[i], hmasses[i-1], spfh[i-1])
  toc("fforce solve")
  f.smearedOneLinkForce(sf, ftmp, g)
  toc("fforce olf")

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
  threads:
    for mu in 0..<f.len:
      p[mu] -= t*f[mu]
  toc("mdv")

func sq(x:float):float = x*x

proc mdvf(i:int, sf:proc, t:float) =
  tic()
  let s =
    if i == 0: -0.5*t*(hmasses[0].sq-mass.sq)/mass
    elif i < hmasses.len: -0.5*t*(hmasses[i].sq-hmasses[i-1].sq)/hmasses[i-1]
    else: -0.5*t/hmasses[i-1]
  f.fforce sf, i
  threads:
    for mu in 0..<f.len:
      p[mu] -= s*f[mu]
  toc("mdvf")

proc fgv(t: float) =
  gc.forceA(g, f)
  threads:
    for mu in 0..<g.len:
      for s in g[mu]:
        g[mu][s] := exp((-t)*f[mu][s])*g[mu][s]
proc fgvf(i:int, sf:proc, t:float) =
  tic()
  let t =
    if i == 0: -0.5*t*(hmasses[0].sq-mass.sq)/mass
    elif i < hmasses.len: -0.5*t*(hmasses[i].sq-hmasses[i-1].sq)/hmasses[i-1]
    else: -0.5*t/hmasses[i-1]
  f.fforce sf, i
  threads:
    for mu in 0..<g.len:
      for s in g[mu]:
        g[mu][s] := exp((-t)*f[mu][s])*g[mu][s]
  toc("fgvf")
var gg = lo.newgauge
var sgg = lo.newgauge
proc fgsave =
  threads:
    for mu in 0..<g.len:
      gg[mu] := g[mu]
proc fgsaves =
  threads:
    for mu in 0..<g.len:
      sgg[mu] := sg[mu]
proc fgload =
  threads:
    for mu in 0..<g.len:
      g[mu] := gg[mu]
proc fgloads =
  threads:
    for mu in 0..<g.len:
      sg[mu] := sgg[mu]

# Compined update for sharing computations
proc mdvAll(t: openarray[float]) =
  if t[0] != 0: mdv t[0]
  var ff = false
  for i in 0..hmasses.len:
    if t[i+1] != 0:
      ff = true
      break
  if ff:
    tic()
    let smearedForce = g.smearRephase sg
    toc("mdvAll smear rephase")
    for i in 0..hmasses.len:
      if t[i+1] != 0:
        mdvf(i, smearedForce, t[i+1])
    toc("mdvAll ff")
proc mdvAllfga(ts,gs:openarray[float]) =
  var
    saved = false
    saveds = false
  proc save1 =
    if not saved:
      fgsave()
      saved = true
  proc saves1 =
    if not saveds:
      fgsaves()
      saveds = true
  let
    gt = ts[0] # 0 for gauge
    gg = gs[0]
  # For gauge
  if gg != 0:
    if gt != 0:
      save1()
      if useFG2:
        # Approximate the force gradient update with two Taylor expansions.
        let (tf,tg) = approximateFGcoeff2(gt,gg)
        fgv tg[0]
        mdv tf[0]
        fgload()
        fgv tg[1]
        mdv tf[1]
      else:
        # Approximate the force gradient update with a Taylor expansion.
        let (tf,tg) = approximateFGcoeff(gt,gg)
        # echo "gauge fg: ",tf," ",tg
        fgv tg
        mdv tf
      fgload()
    else:
      quit("Force gradient without the force update.")
  elif gt != 0:
    mdv gt
  # For fermion
  var ff = false
  for i in 0..hmasses.len:
    if ts[i+1] != 0:
      ff = true
      break
  if ff:
    tic()
    let smearedForce = g.smearRephase sg
    toc("mdvAllfga smear rephase")
    for i in 0..hmasses.len:
      let
        ft = ts[i+1]
        fg = gs[i+1]
      if fg != 0:
        if ft != 0:
          save1()
          saves1()
          if useFG2:
            # Approximate the force gradient update with two Taylor expansions.
            let (tf,tg) = approximateFGcoeff2(ft,fg)
            fgvf i,smearedForce,tg[0]
            block:
              let smearedForce2 = g.smearRephase sg
              mdvf i,smearedForce2,tf[0]
            fgload()
            fgloads()
            fgvf i,smearedForce,tg[1]
            block:
              let smearedForce2 = g.smearRephase sg
              mdvf i,smearedForce2,tf[1]
          else:
            # Approximate the force gradient update with a Taylor expansion.
            let (tf,tg) = approximateFGcoeff(ft,fg)
            # echo "fermion fg: ",tf," ",tg
            fgvf i,smearedForce,tg
            let smearedForce2 = g.smearRephase sg
            mdvf i,smearedForce2,tf
          fgload()
          fgloads()
        else:
          quit("Force gradient without the force update.")
      elif ft != 0:
        mdvf i,smearedForce,ft
    toc("mdvAllfga ff")

#[ Nested integrator
let (V, T) = newIntegratorPair(mdvAll, mdt)
var H = mkOmelyan2MN(steps = gsteps div fsteps, V = V[0], T = T)
for i in 0..<hmasses.len:
  H = mkOmelyan2MN(steps = 1, V = V[i+1], T = H)
H = mkOmelyan2MN(steps = fsteps, V = V[hmasses.len+1], T = H)
]#

let
  (V,T) = newIntegratorPair(mdvAllfga, mdt)
  Hg = gintalg(T = T, V = V[0], steps = gsteps)
  Hf = fintalg(T = T, V = V[1], steps = fsteps)
  H = newParallelEvolution(Hg, Hf)
for i in 0..<hmasses.len:
  H.add fintalg(T = T, V = V[i+2], steps = hfsteps[i])

echo H

for n in 1..trajs:
  tic()
  var p2 = 0.0
  var f2 = newseq[float](phi.len)
  threads:
    p.randomTAH r
    for i in 0..<p.len:
      g0[i] := g[i]
  toc("p refresh, save g")
  p2.pnorm2
  toc("p norm2 1")
  g.smearRephase sg
  toc("smear & rephase 1")
  # phi = D(m2)^{-1} D(m1) psi
  for i in 0..<phi.len:
    threads:
      psi[i].gaussian r
      threadBarrier()
      if i != phi.len-1:
        stag.D(ftmp, psi[i], if i==0: -mass else: -hmasses[i-1])  # `-` for bsm.lua convention giving -D⁺
      else:
        stag.D(phi[i], psi[i], -hmasses[i-1])
    if i != phi.len-1:
      stag.solve(phi[i], ftmp, -hmasses[i], spa)
    threads:
      phi[i].odd := 0
  toc("init traj")
  f2.faction
  toc("fa solve 1")
  let (ga0,fa0,t0,h0) = g0.gaction(f2,p2)
  toc("init gauge action")
  echo "Begin H: ",h0,"  Sg: ",ga0,"  Sf: ",fa0,"  T: ",t0

  H.evolve tau
  H.finish
  toc("evolve")

  p2.pnorm2
  toc("p norm2 2")
  g.smearRephase sg
  toc("smear & rephase 2")
  f2.faction
  toc("fa solve 2")
  let (ga1,fa1,t1,h1) = g.gaction(f2,p2)
  toc("final gauge action")
  echo "End H: ",h1,"  Sg: ",ga1,"  Sf: ",fa1,"  T: ",t1

  when ReversibilityCheck:
    block:
      var g1 = lo.newgauge
      var p1 = lo.newgauge
      threads:
        for i in 0..<g1.len:
          g1[i] := g[i]
          p1[i] := p[i]
          p[i] := -1*p[i]
      H.evolve tau
      H.finish
      p2.pnorm2
      toc("p norm2 2")
      g.smearRephase sg
      toc("smear & rephase 2")
      f2.faction
      toc("fa solve 2")
      let (ga1,fa1,t1,h1) = g.gaction(f2,p2)
      var dsf = newseq[float](fa1.len)
      for i in 0..<fa1.len: dsf[i] = fa1[i] - fa0[i]
      echo "Reversed H: ",h1,"  Sg: ",ga1,"  Sf: ",fa1,"  T: ",t1
      echo "Reversibility: dH: ",h1-h0,"  dSg: ",ga1-ga0,"  dSf: ",dsf,"  dT: ",t1-t0
      for i in 0..<g1.len:
        g[i] := g1[i]
        p[i] := p1[i]
    toc("reversibility")

  let
    dH = h1 - h0
    acc = exp(-dH)
    accr = R.uniform
  if accr <= acc:  # accept
    echo "ACCEPT:  dH: ",dH,"  exp(-dH): ",acc,"  r: ",accr
  else:  # reject
    echo "REJECT:  dH: ",dH,"  exp(-dH): ",acc,"  r: ",accr
    threads:
      for i in 0..<g.len:
        g[i] := g0[i]

  echo 6.0*g.plaq

echoTimers()
qexfinalize()
