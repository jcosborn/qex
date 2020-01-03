## staghmc_s.nim (hypsmear) with Hasenbusch masses.

import qex
import gauge, gauge/hypsmear, physics/qcdTypes
import physics/stagSolve
import mdevolve
import math, sequtils

const ReversibilityCheck {.booldefine.} = false

qexinit()

let
  lat = intSeqParam("lat", @[8,8,8,8])
  #lat = @[8,8,8]
  #lat = @[32,32]
  #lat = @[1024,1024]
  lo = lat.newLayout
  #gc = GaugeActionCoeffs(plaq:6)
  gc = GaugeActionCoeffs(plaq:6,adjplaq:1)
  seed = intParam("seed", 987654321).uint
var r = lo.newRNGField(RngMilc6, seed)
var R:RngMilc6  # global RNG
R.seed(987654321u+seed, 987654321)

var g = lo.newgauge
#g.random r
g.unit

echo 6.0*g.plaq
echo gc.actionA g

var
  p = lo.newgauge
  f = lo.newgauge
  g0 = lo.newgauge

let mass = floatParam("mass", 0.1)
let hmasses = floatSeqParam("hmasses", @[0.2,0.4])  # Hasenbusch masses

var
  ftmp = lo.ColorVector()
  phi = newseq[typeof(lo.ColorVector())](hmasses.len+1)
  psi = newseq[typeof(lo.ColorVector())](hmasses.len+1)
for i in 0..<phi.len:
  phi[i] = lo.ColorVector()
  psi[i] = lo.ColorVector()

var spa = initSolverParams()
#spa.subsetName = "even"
spa.r2req = floatParam("arsq", 1e-20)
spa.maxits = 10000
var spf = initSolverParams()  # TODO: separate for hmasses
#spf.subsetName = "even"
spf.r2req = floatParam("frsq", 1e-12)
spf.maxits = 10000
spf.verbosity = 0

let
  tau = floatParam("tau", 1.0)
  gsteps = intParam("gsteps", 4)
  fsteps = intParam("fsteps", 4)  # TODO: separate for hmasses
  trajs = intParam("ntraj", 10)

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
    t = 0.5*p2
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
  stag.solve(ftmp, phi[i], if i==0: mass else: hmasses[i-1], spf)
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

# For force gradient update
const useFG = true
const useApproxFG2 = false
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
proc fgsave =
  threads:
    for mu in 0..<g.len:
      gg[mu] := g[mu]
proc fgload =
  threads:
    for mu in 0..<g.len:
      g[mu] := gg[mu]

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
  let
    gt = ts[0] # 0 for gauge
    gg = gs[0]
  # For gauge
  if gg != 0:
    if gt != 0:
      fgsave()
      if useApproxFG2:
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
          fgsave()
          if useApproxFG2:
            # Approximate the force gradient update with two Taylor expansions.
            let (tf,tg) = approximateFGcoeff2(ft,fg)
            fgvf i,smearedForce,tg[0]
            mdvf i,smearedForce,tf[0]
            fgload()
            fgvf i,smearedForce,tg[1]
            mdvf i,smearedForce,tf[1]
          else:
            # Approximate the force gradient update with a Taylor expansion.
            let (tf,tg) = approximateFGcoeff(ft,fg)
            # echo "fermion fg: ",tf," ",tg
            fgvf i,smearedForce,tg
            mdvf i,smearedForce,tf
          fgload()
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

# Omelyan's triple star integrators, see Omelyan et. al. (2003)
when useFG:
  let
    (V,T) = newIntegratorPair(mdvAllfga, mdt)
    # mkOmelyan4MN4F2GVG(steps = gsteps, V = V[0], T = T),
    # mkOmelyan4MN4F2GV(steps = gsteps, V = V[0], T = T),
    # mkOmelyan4MN5F1GV(steps = gsteps, V = V[0], T = T),
    # mkOmelyan4MN5F1GP(steps = gsteps, V = V[0], T = T),
    # mkOmelyan4MN5F2GV(steps = gsteps, V = V[0], T = T),
    # mkOmelyan4MN5F2GP(steps = gsteps, V = V[0], T = T),
    # mkOmelyan6MN5F3GP(steps = gsteps, V = V[0], T = T),
    Hg = mkOmelyan4MN5F2GP(steps = gsteps, V = V[0], T = T)
    Hf = mkOmelyan4MN5F2GP(steps = fsteps, V = V[1], T = T)
    H = newParallelEvolution(Hg, Hf)
  for i in 0..<hmasses.len:
    H.add mkOmelyan4MN5F2GP(steps = fsteps, V = V[i+2], T = T)
else:
  let
    (V,T) = newIntegratorPair(mdvAll, mdt)
    # mkOmelyan2MN(steps = gsteps, V = V[0], T = T),
    # mkOmelyan4MN5FP(steps = gsteps, V = V[0], T = T),
    # mkOmelyan4MN5FV(steps = gsteps, V = V[0], T = T),
    # mkOmelyan6MN7FV(steps = gsteps, V = V[0], T = T),
    Hg = mkOmelyan6MN7FV(steps = gsteps, V = V[0], T = T)
    Hf = mkOmelyan6MN7FV(steps = fsteps, V = V[1], T = T)
    H = newParallelEvolution(Hg, Hf)
  for i in 0..<hmasses.len:
    H.add mkOmelyan6MN7FV(steps = fsteps, V = V[i+2], T = T)

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
        stag.D(ftmp, psi[i], if i==0: mass else: hmasses[i-1])
      else:
        stag.D(phi[i], psi[i], hmasses[i-1])
    if i != phi.len-1:
      stag.solve(phi[i], ftmp, hmasses[i], spa)
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
