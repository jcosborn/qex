## staghmc.nim with Hasenbusch masses.

import qex, gauge, physics/qcdTypes, physics/stagSolve
import mdevolve
import math, sequtils, times

const ReversibilityCheck {.booldefine.} = false

qexinit()

letParam:
  lat = @[8,8,8,8]
  beta = 6.0
  adjFac = -0.25
  tau = 2.0
  gsteps = 64
  fsteps = 32
  trajs = 10
  seed:uint64 = int(1000*epochTime())
  mass = 0.1
  hmasses = @[0.2,0.4]  # Hasenbusch masses
  hfsteps = @[fsteps,fsteps]  # nsteps for Hasenbusch masses
  arsq = 1e-20
  frsq = 1e-12
  maxits = 10000

echoParams()

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

let stag = newStag(g)
var spa = initSolverParams()
#spa.subsetName = "even"
spa.r2req = arsq
spa.maxits = maxits
var spf = initSolverParams()  # TODO: separate for hmasses
#spf.subsetName = "even"
spf.r2req = frsq
spf.maxits = maxits
spf.verbosity = 0

template rephase(g: typed) =
  g.setBC
  threadBarrier()
  g.stagPhase

template pnorm2(p2:float) =
  threads:
    var p2t = 0.0
    for i in 0..<p.len:
      p2t += p[i].norm2
    threadMaster: p2 = p2t
    g.rephase

proc gaction(g:auto, f2:seq[float], p2:float):auto =
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
    g.rephase

proc olf(f: var auto, v1: auto, v2: auto) =
  var t {.noInit.}: type(f)
  for i in 0..<v1.len:
    for j in 0..<v2.len:
      t[i,j] := v1[i] * v2[j].adj
  projectTAH(f, t)

proc oneLinkForce(f: auto, p: auto, g: auto) =
  let t = newTransporters(g, p, 1)
  for mu in 0..<g.len:
    discard t[mu] ^* p
  for mu in 0..<g.len:
    for i in f[mu]:
      olf(f[mu][i], p[i], t[mu].field[i])
    for i in f[mu].odd:
      f[mu][i] *= -1

proc fforce(f: auto, i: int) =
  tic()
  threads:
    g.rephase
  toc("fforce rephase")
  stag.solve(ftmp, phi[i], if i==0: mass else: hmasses[i-1], spf)
  toc("fforce solve")
  f.oneLinkForce(ftmp, g)
  toc("fforce olf")
  threads:
    g.rephase
  toc("fforce rephase 2")

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

proc mdvf(i:int, t:float) =
  tic()
  let s =
    if i == 0: -0.5*t*(hmasses[0].sq-mass.sq)/mass
    elif i < hmasses.len: -0.5*t*(hmasses[i].sq-hmasses[i-1].sq)/hmasses[i-1]
    else: -0.5*t/hmasses[i-1]
  f.fforce i
  threads:
    for mu in 0..<f.len:
      p[mu] -= s*f[mu]
  toc("mdvf")

# For force gradient update
const useFG = false
const useApproxFG2 = false
proc fgv(t: float) =
  gc.forceA(g, f)
  threads:
    for mu in 0..<g.len:
      for s in g[mu]:
        g[mu][s] := exp((-t)*f[mu][s])*g[mu][s]
proc fgvf(i:int, t:float) =
  tic()
  let t =
    if i == 0: -0.5*t*(hmasses[0].sq-mass.sq)/mass
    elif i < hmasses.len: -0.5*t*(hmasses[i].sq-hmasses[i-1].sq)/hmasses[i-1]
    else: -0.5*t/hmasses[i-1]
  f.fforce i
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
  # TODO: actually share computation.
  # For now, just do it separately.
  if t[0] != 0: mdv t[0]
  for i in 0..hmasses.len:
    if t[i+1] != 0:
      mdvf(i, t[i+1])
proc mdvAllfga(ts,gs:openarray[float]) =
  # TODO: actually share computation.
  # For now, just do it separately.
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
          fgvf i,tg[0]
          mdvf i,tf[0]
          fgload()
          fgvf i,tg[1]
          mdvf i,tf[1]
        else:
          # Approximate the force gradient update with a Taylor expansion.
          let (tf,tg) = approximateFGcoeff(ft,fg)
          # echo "fermion fg: ",tf," ",tg
          fgvf i,tg
          mdvf i,tf
        fgload()
      else:
        quit("Force gradient without the force update.")
    elif ft != 0:
      mdvf i,ft

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
    H.add mkOmelyan4MN5F2GP(steps = hfsteps[i], V = V[i+2], T = T)
else:
  let
    (V,T) = newIntegratorPair(mdvAll, mdt)
    # mkOmelyan2MN(steps = gsteps, V = V[0], T = T),
    # mkOmelyan4MN5FP(steps = gsteps, V = V[0], T = T),
    # mkOmelyan4MN5FV(steps = gsteps, V = V[0], T = T),
    # mkOmelyan6MN7FV(steps = gsteps, V = V[0], T = T),
    Hg = mkOmelyan2MN(steps = gsteps, V = V[0], T = T)
    Hf = mkOmelyan2MN(steps = fsteps, V = V[1], T = T)
    H = newParallelEvolution(Hg, Hf)
  for i in 0..<hmasses.len:
    H.add mkOmelyan2MN(steps = hfsteps[i], V = V[i+2], T = T)

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
  toc("p norm2 1, rephase")
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
  toc("fa solve 1, rephase")
  let (ga0,fa0,t0,h0) = g0.gaction(f2,p2)
  toc("init gauge action")
  echo "Begin H: ",h0,"  Sg: ",ga0,"  Sf: ",fa0,"  T: ",t0

  H.evolve tau
  H.finish
  toc("evolve")

  p2.pnorm2
  toc("p norm2 2, rephase")
  f2.faction
  toc("fa solve 2, rephase")
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
      f2.faction
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
