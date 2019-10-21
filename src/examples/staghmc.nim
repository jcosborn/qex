import qex
import gauge, physics/qcdTypes
import physics/stagSolve
import mdevolve

qexinit()

let
  lat = intSeqParam("lat", @[8,8,8,8])
  #lat = @[8,8,8]
  #lat = @[32,32]
  #lat = @[1024,1024]
  lo = lat.newLayout
  #gc = GaugeActionCoeffs(plaq:6)
  gc = GaugeActionCoeffs(plaq:6,adjplaq:1)
var r = lo.newRNGField(RngMilc6, 987654321)
var R:RngMilc6  # global RNG
R.seed(987654321, 987654321)

var g = lo.newgauge
#g.random r
g.unit

echo 6.0*g.plaq
echo g.gaugeAction2 gc

var
  p = lo.newgauge
  f = lo.newgauge
  g0 = lo.newgauge
  phi = lo.ColorVector()
  psi = lo.ColorVector()

let mass = floatParam("mass", 0.1)
let stag = newStag(g)
var spa = initSolverParams()
#spa.subsetName = "even"
spa.r2req = floatParam("arsq", 1e-20)
spa.maxits = 10000
var spf = initSolverParams()
#spf.subsetName = "even"
spf.r2req = floatParam("frsq", 1e-12)
spf.maxits = 10000
spf.verbosity = 0

let
  tau = floatParam("tau", 1.0)
  gsteps = intParam("gsteps", 100)
  fsteps = intParam("fsteps", 100)
  trajs = intParam("ntraj", 10)

template rephase(g: typed) =
  g.setBC
  threadBarrier()
  g.stagPhase

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

proc fforce(f: any) =
  tic()
  threads:
    g.rephase
  toc("fforce rephase")
  stag.solve(psi, phi, mass, spf)
  toc("fforce solve")
  #stagD(stag.so, psi, g, psi, 0.0)
  f.oneLinkForce(psi, g)
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

proc mdvf(t: float) =
  tic()
  #let s = t*floatParam("s", 1.0)
  let s = -0.5*t/mass
  f.fforce()
  threads:
    for mu in 0..<f.len:
      p[mu] -= s*f[mu]
  toc("mdvf")

proc mdvf2(t: float) =
  mdv(t)
  mdvf(t)

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
    g.rephase
    threadBarrier()
    stag.D(phi, psi, mass)
    threadBarrier()
    phi.odd := 0
  toc("init traj")
  stag.solve(psi, phi, mass, spa)
  toc("fa solve 1")
  threads:
    var psi2 = psi.norm2()
    threadMaster: f2 = psi2
    g.rephase
  let
    ga0 = g0.gaugeAction2 gc
    fa0 = 0.5*f2
    t0 = 0.5*p2
    h0 = ga0 + fa0 + t0
  toc("init gauge action")
  echo "Begin H: ",h0,"  Sg: ",ga0,"  Sf: ",fa0,"  T: ",t0

  H.evolve tau
  toc("evolve")

  threads:
    var p2t = 0.0
    for i in 0..<p.len:
      p2t += p[i].norm2
    threadMaster: p2 = p2t
    g.rephase
  toc("p norm2, rephase")
  stag.solve(psi, phi, mass, spa)
  toc("fa solve 2")
  threads:
    var psi2 = psi.norm2()
    threadMaster: f2 = psi2
    g.rephase
  let
    ga1 = g.gaugeAction2 gc
    fa1 = 0.5*f2
    t1 = 0.5*p2
    h1 = ga1 + fa1 + t1
  toc("final gauge action")
  echo "End H: ",h1,"  Sg: ",ga1,"  Sf: ",fa1,"  T: ",t1

  #when true:
  when false:
    block:
      var g1 = lo.newgauge
      var p1 = lo.newgauge
      threads:
        for i in 0..<g1.len:
          g1[i] := g[i]
          p1[i] := p[i]
          p[i] := -1*p[i]
      H.evolve tau
      threads:
        var p2t = 0.0
        for i in 0..<p.len:
          p2t += p[i].norm2
        threadMaster: p2 = p2t
      let
        ga1 = g.gaugeAction2 gc
        t1 = 0.5*p2
        h1 = ga1 + t1
      echo "Reversed H: ",h1,"  Sg: ",ga1,"  T: ",t1
      echo "Reversibility: dH: ",h1-h0,"  dSg: ",ga1-ga0,"  dT: ",t1-t0
      #echo p[0][0]
      for i in 0..<g1.len:
        g[i] := g1[i]
        p[i] := p1[i]

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
