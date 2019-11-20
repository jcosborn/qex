## staghmc.nim with Hasenbusch masses.

import qex
import gauge, physics/qcdTypes
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
echo g.gaugeAction2 gc

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

let stag = newStag(g)
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
  gsteps = intParam("gsteps", 100)
  fsteps = intParam("fsteps", 20)  # TODO: separate for hmasses
  trajs = intParam("ntraj", 10)

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

proc gaction(g:any, f2:seq[float], p2:float):auto =
  let
    ga = g.gaugeAction2 gc
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
    g.rephase

proc olf(f: var any, v1: any, v2: any) =
  var t {.noInit.}: type(f)
  for i in 0..<v1.len:
    for j in 0..<v2.len:
      t[i,j] := v1[i] * v2[j].adj
  projectTAH(f, t)

proc oneLinkForce(f: any, p: any, g: any) =
  let t = newTransporters(g, p, 1)
  for mu in 0..<g.len:
    discard t[mu] ^* p
  for mu in 0..<g.len:
    for i in f[mu]:
      olf(f[mu][i], p[i], t[mu].field[i])
    for i in f[mu].odd:
      f[mu][i] *= -1

proc fforce(f: any, i: int) =
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
  f.gaugeforce2(g, gc)
  threads:
    for mu in 0..<f.len:
      p[mu] -= t*f[mu]
  toc("mdv")

proc mdvf(i:int, t:float) =
  proc sq(x:float):float = x*x
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

#proc mdvf2(t: float) =
#  mdv(t)
#  mdvf(t)

# For FGYin11
proc fgv(t: float) =
  f.gaugeforce2(g, gc)
  threads:
    for mu in 0..<g.len:
      for s in g[mu]:
        g[mu][s] := exp((-t)*f[mu][s])*g[mu][s]
var gg = lo.newgauge
proc fgsave =
  threads:
    for mu in 0..<g.len:
      gg[mu] := g[mu]
proc fgload =
  threads:
    for mu in 0..<g.len:
      g[mu] := gg[mu]

# Compined update for sharing computations, requires mdevolve v1
proc mdvAll(t: openarray[float]) =
  # TODO: actually share computation.
  # For now, just do it separately.
  if t[0] != 0: mdv t[0]
  for i in 0..hmasses.len:
    if t[i+1] != 0:
      mdvf(i, t[i+1])

#[ mdevolve v0
let
  # H = mkLeapfrog(steps = steps, V = mdv, T = mdt)
  # H = mkSW92(steps = steps, V = mdv, T = mdt)
  #H = mkOmelyan2MN(steps = gsteps, V = mdvf2, T = mdt)
  # H = mkOmelyan4MN4FP(steps = steps, V = mdv, T = mdt)
  #H = mkOmelyan4MN5FV(steps = gsteps, V = mdvf2, T = mdt)
  #H = mkFGYin11(steps = steps, V = mdv, T = mdt, Vfg = fgv, save = fgsave(), load = fgload())
  #Hg = mkLeapfrog(steps = gsteps, V = mdv, T = mdt, shared=0)
  #Hf = mkLeapfrog(steps = fsteps, V = mdvf, T = mdt, shared=0)
  Hg = mkOmelyan2MN(steps = gsteps, V = mdv, T = mdt, shared=0)
  Hf = mkOmelyan2MN(steps = fsteps, V = mdvf, T = mdt, shared=0)
  #Hg = mkOmelyan4MN5FV(steps = gsteps, V = mdv, T = mdt, shared=1)
  #Hf = mkOmelyan4MN5FV(steps = fsteps, V = mdvf, T = mdt, shared=1)
  H = mkSharedEvolution(Hg, Hf)
]#

# Nested integrator
let (V, T) = wrap(mdvAll, mdt)
var H = mkOmelyan2MN(steps = gsteps div fsteps, V = V[0], T = T)
for i in 0..<hmasses.len:
  H = mkOmelyan2MN(steps = 1, V = V[i+1], T = H)
H = mkOmelyan2MN(steps = fsteps, V = V[hmasses.len+1], T = H)

#[ TODO: test parallel evolution
let
  (V, T) = wrap(mdvAll, mdt)
  Hg = mkOmelyan2MN(steps = gsteps, V = V[0], T = mdt)
  Hf = mkOmelyan2MN(steps = fsteps, V = V[1], T = mdt)
var H = newParallelEvolution(Hg, Hf)
for i in 0..<hmasses.len:
  H.add mkOmelyan2MN(steps = fsteps, V = V[i+2], T = mdt)
]#

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
        stag.D(ftmp, psi[i], if i==0: mass else: hmasses[i-1])
      else:
        stag.D(phi[i], psi[i], hmasses[i-1])
    if i != phi.len-1:
      stag.solve(phi[i], ftmp, hmasses[i], spa)
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
