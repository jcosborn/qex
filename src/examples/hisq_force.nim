# C.T. Peterson: force test inspired from conversation with Peter Boyle
# See Grid implementation here:
#   -https://github.com/paboyle/Grid/blob/develop/tests/forces/Test_bdy.cc
import qex
import gauge/[hisqsmear]
import physics/[stagD,stagSolve]

qexInit()

defaultSetup()
var
  sg = lo.newGauge()
  sgl = lo.newGauge()
  f = lo.newGauge()
  ff = lo.newGauge()
  p = lo.newGauge()
  phi = lo.ColorVector()
  psi = lo.ColorVector()
  r = lo.newRNGField(RngMilc6,123456789)
  mass = 0.1
  eps = floatParam("eps", 1e-4)
  spa = initSolverParams()
  spf = initSolverParams()
  info: PerfInfo
  hisq = newHisq()
let
  stag = newStag3(sg,sgl)
  arsq = 1e-20
  frsq = 1e-12

spa.r2req = arsq
spa.maxits = 10000
spf.r2req = frsq
spf.maxits = 10000
spf.verbosity = 1

# -- Generic

template rephase(g: auto) =
  g.setBC
  threadBarrier()
  g.stagPhase

proc smearRephase(g: auto, sg,sgl: auto): auto {.discardable.} =
  tic()
  threads: g.rephase
  let smearedForce = hisq.smearGetForce(g,sg,sgl)
  threads: g.rephase
  smearedForce

proc reTrMul(x,y:auto):auto =
  var d: type(eval(toDouble(redot(x[0],y[0]))))
  for ir in x: d += redot(x[ir].adj, y[ir])
  result = simdSum(d)
  x.l.threadRankSum(result)

# -- Action calculation

proc action(): float =
  var s: float
  stag.solve(psi, phi, mass, spa)
  threads:
    var st = psi.norm2
    threadMaster: s = st
  result = 0.5*s

# -- Force calculation & momentum update

proc smearedOneAndThreeLinkForce(f: auto, smearedForce: proc, p: auto, g:auto) =
  # reverse accumulation of the derivative
  # 1. Dslash
  var
    f1 = f.newOneOf()
    f3 = f.newOneOf()
    ff = f.newOneOf()
    t,t3: array[4,Shifter[typeof(p),typeof(p[0])]]
  for mu in 0..<f.len:
    t[mu] = newShifter(p,mu,1)
    discard t[mu] ^* p
    t3[mu] = newShifter(p,mu,3)
    discard t3[mu] ^* p
  const n = p[0].len
  threads:
    for mu in 0..<f.len:
      for i in f[mu]:
        forO a, 0, n-1:
          forO b, 0, n-1:
            f1[mu][i][a,b] := p[i][a] * t[mu].field[i][b].adj
            f3[mu][i][a,b] := p[i][a] * t3[mu].field[i][b].adj

  # 2. correcting phase
  threads:
    g.rephase
    for mu in 0..<f.len:
      for i in f[mu].odd:
        f1[mu][i] *= -1
        f3[mu][i] *= -1

  # 3. smearing
  ff.smearedForce(f1,f3)

  # 4. Tₐ ReTr( Tₐ U F† )
  threads:
    for mu in 0..<f.len:
      for i in f[mu]:
        var s {.noinit.}: typeof(f[0][0])
        #s := g[mu][i]*ff[mu][i]
        s := ff[mu][i] * g[mu][i].adj
        f[mu][i].projectTAH(s)
    threadBarrier()
    g.rephase

proc fforce(f: auto) =
  tic()
  let smearedForce = g.smearRephase(sg,sgl)
  toc("fforce smear rephase")
  stag.solve(psi, phi, mass, spf)
  toc("fforce solve")
  f.smearedOneAndThreeLinkForce(smearedForce, psi, g)
  toc("fforce olf")

proc mdt(delta: float) =
  tic()
  threads:
    for mu in 0..<g.len:
      for s in g[mu]:
        g[mu][s] := exp(delta*p[mu][s])*g[mu][s]

proc mdv() =
  let s = -0.5/mass
  f.fforce()
  threads:
    for mu in 0..<f.len: f[mu] *= s

# -- Test

var p1: float
#g.unit
g.random
threads:
  p.randomTAH r
  psi.gaussian r
  var p2t = 0.0
  for i in 0..<p.len: p2t += p[i].norm2
  threadMaster: p1 = 0.5*p2t
discard g.smearRephase(sg,sgl)
threads:
  stag.D(phi, psi, -mass)
  threadBarrier()
  phi.odd := 0
  psi := 0

# Calculate initial action
let s1 = action()
echo "ACTION 1: ", s1

# Update g and calculate force in middle (don't update p)
mdt(0.5*eps); mdv(); mdt(0.5*eps);

# Calculate final action
var p2: float
discard g.smearRephase(sg,sgl)
let s2 = action()
echo "ACTION 2: ", s2
threads:
  var p2t = 0.0
  for i in 0..<p.len:
    p2t += p[i].norm2
  threadMaster: p2 = 0.5*p2t

# Calculate dS = P U dSdU
var dS: float
threads:
  var dSt = 0.0
  for mu in 0..<p.len:
    dSt = dSt - reTrMul(p[mu],f[mu])
  threadMaster: dS = dSt

# Compare differences
let dH = s2+p2-s1-p1
let (dSdt1,dSdt2) = (dS*eps,s2-s1)
echo "dt*dS/dt, dS, difference = ", dSdt1,", ", dSdt2, ", ", dSdt1-dSdt2
