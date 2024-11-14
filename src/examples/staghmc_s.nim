import qex, gauge, gauge/hypsmear, physics/qcdTypes, physics/stagSolve
import mdevolve
import times, macros

const ReversibilityCheck {.booldefine.} = false

qexinit()

let
  lat = intSeqParam("lat", @[8,8,8,8])
  #lat = @[8,8,8]
  #lat = @[32,32]
  #lat = @[1024,1024]
  beta = floatParam("beta", 6.0)
  adjFac = floatParam("adjFac", -0.25)
  tau = floatParam("tau", 2.0)
  gsteps = intParam("gsteps", 64)
  fsteps = intParam("fsteps", 32)
  trajs = intParam("trajs", 10)
  seed = intParam("seed", int(1000*epochTime())).uint64
  mass = floatParam("mass", 0.1)
  arsq = floatParam("arsq", 1e-20)
  frsq = floatParam("frsq", 1e-12)

macro echoparam(x: typed): untyped =
  let n = x.repr
  result = quote do:
    echo `n`, ": ", `x`

echoparam(beta)
echoparam(adjFac)
echoparam(tau)
echoparam(gsteps)
echoparam(fsteps)
echoparam(trajs)
echoparam(seed)
echoparam(mass)
echoparam(arsq)
echoparam(frsq)

let
  gc = GaugeActionCoeffs(plaq: beta, adjplaq: beta*adjFac)
  lo = lat.newLayout
  vol = lo.physVol

var r = lo.newRNGField(RngMilc6, seed)
var R:RngMilc6  # global RNG
R.seed(seed, 987654321)

var g = lo.newgauge
#g.random r
g.unit

echo 6.0*g.plaq
echo g.gaugeAction2 gc
echo gc.actionA g

var
  p = lo.newgauge
  f = lo.newgauge
  g0 = lo.newgauge
  phi = lo.ColorVector()
  psi = lo.ColorVector()

var spa = initSolverParams()
#spa.subsetName = "even"
spa.r2req = arsq
spa.maxits = 10000
var spf = initSolverParams()
#spf.subsetName = "even"
spf.r2req = frsq
spf.maxits = 10000
spf.verbosity = 0

var
  info: PerfInfo
  coef = HypCoefs(alpha1:0.4, alpha2:0.5, alpha3:0.5)
echo "smear = ",coef
var sg = lo.newGauge
let stag = newStag(sg)

proc smearRephase(g: auto, sg: auto):auto {.discardable.} =
  tic()
  #let smearedForce = coef.smear(g, sg, info)
  let smearedForce = coef.smearGetForce(g, sg, info)
  toc("smear")
  threads:
    sg.setBC
    threadBarrier()
    sg.stagPhase
  toc("BC & Phase")
  smearedForce

proc smearedOneLinkForce(f: auto, smearedForce: proc, p: auto, g:auto) =
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

proc fforce(f: auto) =
  tic()
  let smearedForce = g.smearRephase sg
  toc("fforce smear rephase")
  stag.solve(psi, phi, mass, spf)
  toc("fforce solve")
  f.smearedOneLinkForce(smearedForce, psi, g)
  toc("fforce olf")

proc mdt(t: float) =
  tic()
  threads:
    for mu in 0..<g.len:
      for s in g[mu]:
        g[mu][s] := exp(t*p[mu][s])*g[mu][s]
  toc("mdt")
  GC_fullCollect()
  toc("GC")
proc mdv(t: float) =
  tic()
  gc.forceA(g, f)
  threads:
    for mu in 0..<f.len:
      p[mu] -= t*f[mu]
  toc("mdv")

proc mdvf(t: float) =
  tic()
  #let s = t*floatParam("s", 1.0)
  let s = -0.5*t/mass
  f.fforce()
  threads:
    for mu in 0..<f.len:
      p[mu] -= s*f[mu]
  toc("mdvf")
  GC_fullCollect()
  toc("GC")

# For force gradient update
#const useFG = true
const useFG = false
const useApproxFG2 = false
proc fgv(t: float) =
  tic()
  gc.forceA(g, f)
  threads:
    for mu in 0..<g.len:
      for s in g[mu]:
        g[mu][s] := exp((-t)*f[mu][s])*g[mu][s]
  toc("fgv")
proc fgvf(t: float) =
  tic()
  let t = -0.5*t/mass
  f.fforce()
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
  if t[1] != 0: mdvf t[1]
proc mdvAllfga(ts,gs:openarray[float]) =
  # TODO: actually share computation.
  # For now, just do it separately.
  let
    gt = ts[0] # 0 for gauge
    gg = gs[0]
    ft = ts[1] # 1 for fermion
    fg = gs[1]
  # echo "mdvAll: gauge: ",gt," ",gg,"  fermion: ",ft," ",fg
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
  if fg != 0:
    if ft != 0:
      fgsave()
      if useApproxFG2:
        # Approximate the force gradient update with two Taylor expansions.
        let (tf,tg) = approximateFGcoeff2(ft,fg)
        fgvf tg[0]
        mdvf tf[0]
        fgload()
        fgvf tg[1]
        mdvf tf[1]
      else:
        # Approximate the force gradient update with a Taylor expansion.
        let (tf,tg) = approximateFGcoeff(ft,fg)
        # echo "fermion fg: ",tf," ",tg
        fgvf tg
        mdvf tf
      fgload()
    else:
      quit("Force gradient without the force update.")
  elif ft != 0:
    mdvf ft

#[ Nested integrators
let
  (V,T) = newIntegratorPair(mdvAll, mdt)
  Hg = mkOmelyan2MN(steps = gsteps div fsteps, V = V[0], T = T)
  H = mkOmelyan2MN(steps = fsteps, V = V[1], T = Hg)
]#
let
  # Omelyan's triple star integrators, see Omelyan et. al. (2003)
  H =
    when useFG:
      let
        (VAllG,T) = newIntegratorPair(mdvAllfga, mdt)
        V = VAllG[0]
        Vf = VAllG[1]
      newParallelEvolution(
        # mkOmelyan4MN4F2GVG(steps = gsteps, V = V, T = T),
        # mkOmelyan4MN4F2GV(steps = gsteps, V = V, T = T),
        # mkOmelyan4MN5F1GV(steps = gsteps, V = V, T = T),
        # mkOmelyan4MN5F1GP(steps = gsteps, V = V, T = T),
        # mkOmelyan4MN5F2GV(steps = gsteps, V = V, T = T),
        mkOmelyan4MN5F2GP(steps = gsteps, V = V, T = T),
        # mkOmelyan6MN5F3GP(steps = gsteps, V = V, T = T),
        mkOmelyan4MN5F2GP(steps = fsteps, V = Vf, T = T))
    else:
      let
        (VAll,T) = newIntegratorPair(mdvAll, mdt)
        V = VAll[0]
        Vf = VAll[1]
      newParallelEvolution(
        mkOmelyan2MN(steps = gsteps, V = V, T = T),
        mkOmelyan2MN(steps = fsteps, V = Vf, T = T))
        # mkOmelyan4MN5FP(steps = gsteps, V = V, T = T),
        # mkOmelyan4MN5FV(steps = gsteps, V = V, T = T),
        #mkOmelyan6MN7FV(steps = gsteps, V = V, T = T),
        #mkOmelyan6MN7FV(steps = fsteps, V = Vf, T = T))

for n in 1..trajs:
  tic()
  var p2 = 0.0
  var f2 = 0.0
  threads:
    p.randomTAH r
    var p2t = 0.0
    for i in 0..<p.len:
      p2t += p[i].norm2
      g0[i] := g[i]
    threadMaster: p2 = p2t
    psi.gaussian r
  toc("init traj")
  g.smearRephase sg
  toc("init traj smear & rephase")
  threads:
    stag.D(phi, psi, -mass)   # match bsm.lua convention
    threadBarrier()
    phi.odd := 0
  toc("init traj D")
  stag.solve(psi, phi, mass, spa)
  toc("fa solve 1")
  threads:
    var psi2 = psi.norm2()
    threadMaster: f2 = psi2
  let
    ga0 = gc.actionA g0
    fa0 = 0.5*f2
    t0 = 0.5*p2 - (16*vol).float
    h0 = ga0 + fa0 + t0
  toc("init gauge action")
  echo "Begin H: ",h0,"  Sg: ",ga0,"  Sf: ",fa0,"  T: ",t0

  H.evolve tau
  H.finish
  toc("evolve")

  threads:
    var p2t = 0.0
    for i in 0..<p.len:
      p2t += p[i].norm2
    threadMaster: p2 = p2t
  toc("p norm2")
  g.smearRephase sg
  toc("smear & rephase")
  stag.solve(psi, phi, mass, spa)
  toc("fa solve 2")
  threads:
    var psi2 = psi.norm2()
    threadMaster: f2 = psi2
  let
    ga1 = gc.actionA g
    fa1 = 0.5*f2
    t1 = 0.5*p2 - (16*vol).float
    h1 = ga1 + fa1 + t1
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
      threads:
        var p2t = 0.0
        for i in 0..<p.len:
          p2t += p[i].norm2
        threadMaster: p2 = p2t
      g.smearRephase sg
      toc("p norm2, rephase")
      stag.solve(psi, phi, mass, spa)
      toc("fa solve 2")
      threads:
        var psi2 = psi.norm2()
        threadMaster: f2 = psi2
      let
        ga1 = gc.actionA g
        fa1 = 0.5*f2
        t1 = 0.5*p2
        h1 = ga1 + fa1 + t1
      echo "Reversed H: ",h1,"  Sg: ",ga1,"  Sf: ",fa1,"  T: ",t1
      echo "Reversibility: dH: ",h1-h0,"  dSg: ",ga1-ga0,"  dSf: ",fa1-fa0,"  dT: ",t1-t0
      #echo p[0][0]
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
