import qex, gauge, physics/[qcdTypes,stagSolve]
import times, macros, algorithm, sequtils
import hmc/metropolis
import observables/sources
import quda/qudaWrapper
import hmc/agradOps

qexInit()
echo "rank ", myRank, "/", nRanks
threads:
  echo "thread ", threadNum, "/", numThreads

proc `:=`*(r: var seq, x: seq) =
  for i in 0..<r.len: r[i] := x[i]

let
  lat = intSeqParam("lat", @[8,8,8,8])
  beta = floatParam("beta", 5.6)
  tau = floatParam("tau", 0.04)
  fixtau = (intParam("fixtau",0) != 0)
  fixparams = (intParam("fixparams",0) != 0)
  fixhmasses = (intParam("fixhmasses",0) != 0)
  md = stringParam("md", "aba")
  upit = intParam("upit", 1)
  anneal = floatParam("anneal", 0.99)
  alpha = floatParam("alpha", 0)
  checkg = (intParam("checkg", 0) != 0)
  nsteps = intParam("nsteps", 1)
  nwarm = intParam("nwarm", 10)
  warmmd = (intParam("warmmd", 1) != 0)
  ntrain = intParam("ntrain", 10)
  trajs = intParam("trajs", 10)
  nf = intParam("nf", 1)
  mass = floatParam("mass", 0.1)
  hmasses = floatSeqParam("hmasses")  # Hasenbusch masses
  arsq = floatParam("arsq", 1e-20)
  frsq = floatParam("frsq", 1e-12)
  seed0 = defaultComm.broadcast(int(1000*epochTime()))
  seed = uint64 intParam("seed", seed0)
  infn = stringParam("infn", "")
  outfn = stringParam("outfn", "")
var
  lrate = floatParam("lrate", 0.001)
  lrateh = floatParam("lrateh", lrate)
  pt0 = floatParam("t0", 0)
  pg0 = floatParam("g0", 0)
  pf0 = floatParam("f0", 0)
  pgf0 = floatParam("gf0", 0)
  pff0 = floatParam("ff0", 0)
  pt1 = floatParam("t1", 0)
  pg1 = floatParam("g1", 0)
  pf1 = floatParam("f1", 0)
  pgf1 = floatParam("gf1", 0)
  pff1 = floatParam("ff1", 0)
  pt2 = floatParam("t2", 0)
  pg2 = floatParam("g2", 0)
  pf2 = floatParam("f2", 0)
  pgf2 = floatParam("gf2", 0)
  pff2 = floatParam("ff2", 0)
  pt3 = floatParam("t3", 0)
  pg3 = floatParam("g3", 0)
  pf3 = floatParam("f3", 0)
  pgf3 = floatParam("gf3", 0)
  pff3 = floatParam("ff3", 0)
  pt4 = floatParam("t4", 0)
  pg4 = floatParam("g4", 0)
  pf4 = floatParam("f4", 0)
  pgf4 = floatParam("gf4", 0)
  pff4 = floatParam("ff4", 0)

macro echoparam(x: auto): auto =
  let n = x.repr
  result = quote do:
    echo `n`, ": ", `x`
macro echoparam(x: auto, y: untyped): auto =
  let n = y.repr
  result = quote do:
    echo `n`, ": ", `x`

echoparam(beta)
echoparam(tau)
echoparam(fixtau)
echoparam(fixparams)
echoparam(fixhmasses)
echoparam(md)
echoparam(upit)
echoparam(lrate)
echoparam(lrateh)
echoparam(anneal)
echoparam(alpha)
echoparam(checkg)
echoparam(nsteps)
echoparam(nwarm)
echoparam(warmmd)
echoparam(ntrain)
echoparam(trajs)
echoparam(nf)
echoparam(mass)
echoparam(hmasses)
echoparam(arsq)
echoparam(frsq)
echoparam(seed)
echoparam(infn)
echoparam(outfn)
echoparam(pt0, t0)
echoparam(pg0, g0)
echoparam(pf0, f0)
echoparam(pgf0, gf0)
echoparam(pff0, ff0)
echoparam(pt1, t1)
echoparam(pg1, g1)
echoparam(pf1, f1)
echoparam(pgf1, gf1)
echoparam(pff1, ff1)
echoparam(pt2, t2)
echoparam(pg2, g2)
echoparam(pf2, f2)
echoparam(pgf2, gf2)
echoparam(pff2, ff2)

let
  gc = GaugeActionCoeffs(plaq: beta, adjplaq: 0)
  lo = lat.newLayout
  vol = lo.physVol

var r = lo.newRNGField(RngMilc6, seed)
var R:RngMilc6  # global RNG
R.seed(seed, 987654321)

var
  g = lo.newgauge
  p = lo.newgauge
  f = lo.newgauge
  g0 = lo.newgauge
  eta = newseq[typeof(lo.ColorVector())](1+hmasses.len) # random gaussian
  phi = newseq[typeof(lo.ColorVector())](1+hmasses.len) # pseudofermion
  psi = newseq[typeof(lo.ColorVector())](1+hmasses.len) # intermediate
  ftmp = lo.ColorVector()
for i in 0..<phi.len:
  eta[i] = lo.ColorVector()
  phi[i] = lo.ColorVector()
  psi[i] = lo.ColorVector()
type
  Gauge = typeof(g)
  GaugeV = AgVar[Gauge]
  Cvec = typeof(phi[0])
  CvecV = AgVar[Cvec]
template newGaugeV(c: AgTape, x: Gauge): auto = newGaugeFV(c, x)
case infn
of "":
  g.random r
  #g.unit
of "hot":
  g.random r
of "cold":
  g.unit
else:
  echo "Loading gauge field from file: ", infn
  let err = g.loadGauge infn

echo "plaq: ", 6.0*g.plaq
echo "gaugeAction2: ", g.gaugeAction2 gc
echo "actionA: ", gc.actionA g

let stag = newStag(g)
var spa = initSolverParams()
#spa.subsetName = "even"
spa.r2req = arsq
spa.maxits = 10000
#spa.backend = sbQex
var spf = newSeq[type spa](hmasses.len+1)  # fermion force forward
var spfb = newSeq[type spa](hmasses.len+1)  # fermion force backward
for i in 0..<spf.len:
  spf[i] = initSolverParams()
  spf[i].r2req = frsq
  spf[i].maxits = 10000
  spf[i].verbosity = 0
  #spf.backend = sbQex
  spfb[i] = initSolverParams()
  spfb[i].r2req = frsq
  spfb[i].maxits = 10000
  spfb[i].verbosity = 0
  #spf.backend = sbQex

proc norm2*(x: seq): float =
  for i in 0..<x.len:
    result += x[i].norm2
    #echo result
  #result /= (x.len * x[0].l.physVol).float

proc norm2v*(x: seq): float =
  for i in 0..<x.len:
    result += x[i].norm2
    #echo result
  result /= (x.len * x[0].l.physVol).float

proc norm2subtract*(x: seq, y: float): float =
  for i in 0..<x.len:
    result += x[i].norm2subtract(y)

proc rmsStats*(f: Field, scale = 1.0): RunningStat =
  let presult = addr result
  threads:
    let ne = f[0].numNumbers / f[0].simdLength
    var s: RunningStat
    for e in f:
      let t = scale*sqrt(f[e].norm2/ne)
      for i in 0..<simdLength(t):
        s.push t[i]
    if threadNum == 0:
      presult[] += s
    for i in 1..<numThreads:
      threadBarrier()
      if threadNum == i:
        presult[] += s
  var c = getDefaultComm()
  var offset = 1
  while offset < c.size:
    let other = c.rank xor offset
    if other < c.rank:
      c.pushSend(other, result)
      c.waitSend()
    elif other < c.size:
      var t {.noInit.}: RunningStat
      c.pushRecv(other, t)
      c.waitRecv()
      result += t
    offset *= 2

proc rmsStats*(f: seq, scale = 1.0): RunningStat =
  for i in 0..<f.len:
    result += rmsStats(f[i], scale)

var
  gauges = newSeq[Gauge](0)
  moms = newSeq[Gauge](0)
  gaugevs = newSeq[GaugeV](0)
  momvs = newSeq[GaugeV](0)
  params = newSeq[FloatV](0)
  ptemps = newSeq[FloatV](0)
  cvecs = newSeq[Cvec](0)
  cvecvs = newSeq[CvecV](0)
var tape = newAgTape()
proc pushParam(x: float): FloatV =
  params.add tape.newFloatV(x)
  result = params[^1]
proc echoParams =
  var s = "params:"
  for p in params:
    s &= " " & $p.obj
  echo s

gauges.add g0
gaugevs.add tape.newGaugeV(g0)
var p0 = lo.newGauge
moms.add p0
momvs.add tape.newGaugeV(p0)
#params.add tape.newFloatV(tau)
#let vtau = params[0]
let vtau = pushParam(tau)
#vtau.doGrad = false
var nff = newSeq[int](1+hmasses.len)
var vhmasses = newSeq[type pushParam(0.0)](hmasses.len)
for i in 0..<vhmasses.len: vhmasses[i] = pushParam(hmasses[i])
var phiv = newSeq[type tape.newAgVar(phi[0])](phi.len)
for i in 0..<phi.len: phiv[i] = tape.newAgVar(phi[i])

proc pushTemp =
  ptemps.add tape.newFloatV
proc pushGauge =
  gauges.add lo.newGauge
  gaugevs.add tape.newGaugeV(gauges[^1])
proc pushMom =
  moms.add lo.newGauge
  momvs.add tape.newGaugeV(moms[^1])
proc pushCvec =
  cvecs.add lo.ColorVector
  cvecvs.add tape.newAgVar(cvecs[^1])
proc `+`(x: FloatV, y: FloatV): FloatV =
  pushTemp()
  result = ptemps[^1]
  add(result, x, y)
proc `-`(x: FloatV, y: float): FloatV =
  pushTemp()
  result = ptemps[^1]
  sub(result, x, y)
proc `-`(x: FloatV, y: FloatV): FloatV =
  pushTemp()
  result = ptemps[^1]
  sub(result, x, y)
proc `*`(x: float, y: FloatV): FloatV =
  pushTemp()
  result = ptemps[^1]
  mul(result, x, y)
proc `*`(x: FloatV, y: FloatV): FloatV =
  pushTemp()
  result = ptemps[^1]
  mul(result, x, y)
proc `/`(x: FloatV, y: SomeNumber): FloatV =
  pushTemp()
  result = ptemps[^1]
  divd(result, x, float y)
proc `/`(x: SomeNumber, y: FloatV): FloatV =
  pushTemp()
  result = ptemps[^1]
  divd(result, float x, y)
proc `/`(x: FloatV, y: FloatV): FloatV =
  pushTemp()
  result = ptemps[^1]
  divd(result, x, y)

proc addT(veps: FloatV) =
  pushGauge()
  exp(gaugevs[^1], veps, momvs[^1])
  pushGauge()
  mul(gaugevs[^1], gaugevs[^2], gaugevs[^3])

proc addGf(g: GaugeV) =
  if nf == 0: nff[0] += 1
  pushMom()
  gc.gderiv(momvs[^1], g)
  pushMom()
  mulna(momvs[^1], g, momvs[^2])
  pushMom()
  projtah(momvs[^1], momvs[^2])

var gfList = newSeq[int](0)
var gfeList = newSeq[FloatV](0)
proc addGx(veps: FloatV, p,g: GaugeV) =
  addGf(g)
  pushMom()
  xpay(momvs[^1], p, veps, momvs[^2])
  gfList.add momvs.len-2
  gfeList.add veps

proc addG(veps: FloatV) =
  addGx(veps, momvs[^1], gaugevs[^1])

proc addGF(va, vb: FloatV) =
  let p = momvs[^1]
  addGf(gaugevs[^1])
  pushMom()
  exp(momvs[^1], vb, momvs[^2])
  pushMom()
  mul(momvs[^1], momvs[^2], gaugevs[^1])
  addGx(va, p, momvs[^1])

proc addFf(g: GaugeV, i = 0) =
  if nf == 0: return
  pushCvec()
  let cv = cvecvs[^1]
  if i == 0:
    stag.agradSolve(g, cv, phiv[i], mass, spf[i], spfb[i])
  else:
    stag.agradSolve(g, cv, phiv[i], vhmasses[i-1], spf[i], spfb[i])
  pushMom()
  stag.agradStagDeriv(momvs[^1], cv)
  pushMom()
  mulna(momvs[^1], g, momvs[^2])
  pushMom()
  projtah(momvs[^1], momvs[^2])
  nff[i] += 1

proc massFactor0(veps: FloatV, i: int): FloatV =
  if i == 0:
    result = (-0.5/mass) * veps
  else:
    result = -0.5 * (veps/vhmasses[i-1])

proc massFactor(veps: FloatV, i: int): FloatV =
  result = massFactor0(veps, i)
  if i < hmasses.len: # last term is just inverse (no ratio)
    if i == 0:
      result = result * (vhmasses[0]*vhmasses[0]-mass*mass)
    else:
      result = result * (vhmasses[i]*vhmasses[i]-vhmasses[i-1]*vhmasses[i-1])

var ffList = newSeqWith(1+hmasses.len, newSeq[int](0))
var ffeList = newSeqWith(1+hmasses.len, newSeq[FloatV](0))
proc addFx(veps: FloatV; p0,g: GaugeV; i0,i1: int) =
  var p = p0
  if nf == 0: return
  for i in i0..i1:
    addFf(g, i)
    pushMom()
    let va = massFactor(veps, i)
    xpay(momvs[^1], p, va, momvs[^2])
    p = momvs[^1]
    ffList[i].add momvs.len-2
    ffeList[i].add va

proc addF(veps: FloatV) =
  if nf == 0: return
  addFx(veps, momvs[^1], gaugevs[^1], 0, hmasses.len)

proc addF(veps: FloatV, i: int) =
  if nf == 0: return
  addFx(veps, momvs[^1], gaugevs[^1], i, i)

proc addF(veps: FloatV; i0,i1: int) =
  if nf == 0: return
  addFx(veps, momvs[^1], gaugevs[^1], i0, i1)

proc addFF(va, vb: FloatV; i0,i1: int) =
  if nf == 0: return
  var p = momvs[^1]
  let g = gaugevs[^1]
  for i in i0..i1:
    addFf(g, i)
    pushMom()
    let vbx = massFactor(vb, i)
    exp(momvs[^1], vbx, momvs[^2])
    pushMom()
    mul(momvs[^1], momvs[^2], gaugevs[^1])
    addFx(va, p, momvs[^1], i, i)
    p = momvs[^1]

template addFF(va, vb: FloatV) = addFF(va, vb, 0, hmasses.len)

proc addGFF(vag, vbg, vaf, vbf: FloatV; i0,i1: int) =
  var p = momvs[^1]
  let g = gaugevs[^1]
  addGf(g)
  pushMom()
  var s = momvs[^1]
  mul(s, vbg, momvs[^2])
  for i in i0..i1:
    addFf(g, i)
    let vbx = massFactor(vbf, i)
    pushMom()
    xpay(momvs[^1], s, vbx, momvs[^2])
    s = momvs[^1]
  pushMom()
  exp(momvs[^1], s)
  pushMom()
  mul(momvs[^1], momvs[^2], g)
  s = momvs[^1]
  addGx(vag, p, s)
  addFx(vaf, momvs[^1], s, i0, i1)

template addGFF(vag, vbg, vaf, vbf: FloatV) = addGFF(vag, vbg, vaf, vbf, 0, hmasses.len)
template addGFF(va, vb: FloatV) = addGFF(va, vb, va, vb, 0, hmasses.len)

proc setupMDx =
  pushTemp()
  let vtau2 = ptemps[^1]
  mul(vtau2, 0.5, vtau)
  addT(vtau2)
  addG(vtau)
  addT(vtau2)

proc setupMDa =
  pushTemp()
  let vtau2 = ptemps[^1]
  mul(vtau2, 0.5, vtau)
  pushTemp()
  let vtau4 = ptemps[^1]
  mul(vtau4, 0.25, vtau)
  addT(vtau4)
  addG(vtau2)
  addT(vtau2)
  addG(vtau2)
  addT(vtau4)

# fixed integrators

# Order 2

proc setupMDaba =
  let t0 = 0.5 * vtau
  let g0 = vtau
  let f0 = vtau
  addT(t0)
  addG(g0)
  addF(f0)
  addT(t0)

proc setupMDaca =
  if pgf0 == 0: pgf0 = 1.0/6.0
  if pff0 == 0: pff0 = 1.0/6.0
  let t0 = 0.5 * vtau
  let t02 = vtau
  let g0 = vtau
  let f0 = vtau
  let g0g = vtau*vtau*pushParam(pgf0)
  let f0g = if nf==0: g0g else: vtau*vtau*pushParam(pff0)
  addT(t0)
  for i in 0..<nsteps:
    if i!=0: addT(t02)
    addGF(g0, g0g)
    addFF(f0, f0g)
  addT(t0)

proc setupMDababa =
  if pt0 == 0: pt0 = 0.1931833275037836
  #let eps = vtau / nsteps
  let t0 = vtau * pushParam(pt0)
  let t02 = 2 * t0
  let t1 = vtau - t02
  let g0 = 0.5 * vtau
  let f0 = 0.5 * vtau
  addT(t0)
  for i in 0..<nsteps:
    if i!=0: addT(t02)
    addG(g0)
    addF(f0)
    addT(t1)
    addF(f0)
    addG(g0)
  addT(t0)

#proc setupMDagabagabaga =
proc setupMDg5f2 =
  if pt0 == 0: pt0 = 0.1
  if pt1 == 0: pt1 = 0.1
  if pg0 == 0: pg0 = 0.14
  if pg1 == 0: pg1 = 0.24
  let t0 = vtau * pushParam(pt0)
  let t02 = 2 * t0
  let t1 = vtau * pushParam(pt1)
  let t2 = 0.5*vtau - t0 - t1
  let g0 = vtau * pushParam(pg0)
  let g1 = vtau * pushParam(pg1)
  let g2 = vtau - 2*( g0 + g1 )
  let f0 = 0.5*vtau
  addT(t0)
  for i in 0..<nsteps:
    if i!=0: addT(t02)
    addG(g0)
    addT(t1)
    addG(g1)
    addF(f0)
    addT(t2)
    addG(g2)
    addT(t2)
    addF(f0)
    addG(g1)
    addT(t1)
    addG(g0)
  addT(t0)

proc setupMDg5f23 =
  if pt0 == 0: pt0 = 0.09
  if pt1 == 0: pt1 = 0.13
  if pg0 == 0: pg0 = 0.13
  if pg1 == 0: pg1 = 0.23
  if pf0 == 0: pf0 = 0.26
  let t0 = vtau * pushParam(pt0)
  let t02 = 2 * t0
  let t1 = vtau * pushParam(pt1)
  let t2 = 0.5*vtau - t0 - t1
  let g0 = vtau * pushParam(pg0)
  let g1 = vtau * pushParam(pg1)
  let g2 = vtau - 2*( g0 + g1 )
  let f00 = vtau * pushParam(pf0)
  let f01 = vtau - 2*f00
  let f10 = 0.5*vtau
  let i0 = hmasses.len
  let i1 = hmasses.len - 1
  addT(t0)
  for i in 0..<nsteps:
    if i!=0: addT(t02)
    addG(g0)
    addF(f00, i0)
    addT(t1)
    addG(g1)
    addF(f10, 0, i1)
    addT(t2)
    addG(g2)
    addF(f01, i0)
    addT(t2)
    addF(f10, 0, i1)
    addG(g1)
    addT(t1)
    addF(f00, i0)
    addG(g0)
  addT(t0)

proc setupMDg6f24 =
  if pt0 == 0: pt0 = 0.1
  if pt1 == 0: pt1 = 0.11
  if pt2 == 0: pt2 = 0.13
  if pg0 == 0: pg0 = 0.12
  if pg1 == 0: pg1 = 0.22
  if pf0 == 0: pf0 = 0.22
  let t0 = vtau * pushParam(pt0)
  let t02 = 2 * t0
  let t1 = vtau * pushParam(pt1)
  let t2 = vtau * pushParam(pt2)
  let t3 = vtau - 2*(t0 + t1 + t2)
  let g0 = vtau * pushParam(pg0)
  let g1 = vtau * pushParam(pg1)
  let g2 = 0.5*vtau - g0 - g1
  let f00 = vtau * pushParam(pf0)
  let f01 = 0.5*vtau - f00
  let f10 = 0.5*vtau
  let i0 = hmasses.len
  let i1 = hmasses.len - 1
  addT(t0)
  for i in 0..<nsteps:
    if i!=0: addT(t02)
    addG(g0)
    addF(f00, i0)
    addT(t1)
    addG(g1)
    addF(f10, 0, i1)
    addT(t2)
    addG(g2)
    addF(f01, i0)
    addT(t3)
    addF(f01, i0)
    addG(g2)
    addT(t2)
    addF(f10, 0, i1)
    addG(g1)
    addT(t1)
    addF(f00, i0)
    addG(g0)
  addT(t0)

proc setupMDg10f2 =
  if pt0 == 0: pt0 = 0.075
  if pt1 == 0: pt1 = 0.07
  if pt2 == 0: pt2 = 0.06
  if pt3 == 0: pt3 = 0.1
  if pt4 == 0: pt4 = 0.1
  if pg0 == 0: pg0 = 0.095
  if pg1 == 0: pg1 = 0.095
  if pg2 == 0: pg2 = 0.09
  if pg3 == 0: pg3 = 0.09
  #let eps = vtau / nsteps
  let t0 = vtau * pushParam(pt0)
  let t02 = 2 * t0
  let t1 = vtau * pushParam(pt1)
  let t2 = vtau * pushParam(pt2)
  let t3 = vtau * pushParam(pt3)
  let t4 = vtau * pushParam(pt4)
  let t5 = vtau - t02 - 2 * (t1 + t2 + t3 + t4)
  let g0 = vtau * pushParam(pg0)
  let g1 = vtau * pushParam(pg1)
  let g2 = vtau * pushParam(pg2)
  let g3 = vtau * pushParam(pg3)
  let g4 = 0.5 * vtau - (g0 + g1 + g2 + g3)
  let f0 = 0.5 * vtau
  addT(t0)
  for i in 0..<nsteps:
    if i!=0: addT(t02)
    addG(g0)
    addT(t1)
    addG(g1)
    addT(t2)
    addG(g2)
    addF(f0)
    addT(t3)
    addG(g3)
    addT(t4)
    addG(g4)
    addT(t5)
    addG(g4)
    addT(t4)
    addG(g3)
    addT(t3)
    addF(f0)
    addG(g2)
    addT(t2)
    addG(g1)
    addT(t1)
    addG(g0)
  addT(t0)

proc setupMDg10f28 =
  if pt0 == 0: pt0 = 0.075
  if pt1 == 0: pt1 = 0.07
  if pt2 == 0: pt2 = 0.06
  if pt3 == 0: pt3 = 0.09
  if pt4 == 0: pt4 = 0.09
  if pg0 == 0: pg0 = 0.095
  if pg1 == 0: pg1 = 0.095
  if pg2 == 0: pg2 = 0.09
  if pg3 == 0: pg3 = 0.09
  if pf0 == 0: pf0 = 0.11
  if pf1 == 0: pf1 = 0.11
  if pf2 == 0: pf2 = 0.16
  #let eps = vtau / nsteps
  let t0 = vtau * pushParam(pt0)
  let t02 = 2 * t0
  let t1 = vtau * pushParam(pt1)
  let t2 = vtau * pushParam(pt2)
  let t3 = vtau * pushParam(pt3)
  let t4 = vtau * pushParam(pt4)
  let t5 = vtau - t02 - 2 * (t1 + t2 + t3 + t4)
  let g0 = vtau * pushParam(pg0)
  let g1 = vtau * pushParam(pg1)
  let g2 = vtau * pushParam(pg2)
  let g3 = vtau * pushParam(pg3)
  let g4 = 0.5 * vtau - (g0 + g1 + g2 + g3)
  let f00 = vtau * pushParam(pf0)
  let f01 = vtau * pushParam(pf1)
  let f02 = vtau * pushParam(pf2)
  let f03 = 0.5*vtau - f00 - f01 - f02
  let f10 = 0.5*vtau
  let i0 = hmasses.len
  let i1 = hmasses.len - 1
  addT(t0)
  for i in 0..<nsteps:
    if i!=0: addT(t02)
    addG(g0)
    addF(f00, i0)
    addT(t1)
    addG(g1)
    addF(f01, i0)
    addT(t2)
    addG(g2)
    addF(f10, 0, i1)
    addT(t3)
    addG(g3)
    addF(f02, i0)
    addT(t4)
    addG(g4)
    addF(f03, i0)
    addT(t5)
    addF(f03, i0)
    addG(g4)
    addT(t4)
    addF(f02, i0)
    addG(g3)
    addT(t3)
    addF(f10, 0, i1)
    addG(g2)
    addT(t2)
    addF(f01, i0)
    addG(g1)
    addT(t1)
    addF(f00, i0)
    addG(g0)
  addT(t0)

proc setupMDbabab =
  if pg0 == 0: pg0 = 0.1931833275037836
  if pf0 == 0: pf0 = 0.1931833275037836
  let g0 = vtau * pushParam(pg0)
  let f0 = if nf==0: g0 else: vtau * pushParam(pf0)
  let t0 = 0.5 * vtau
  let g02 = 2 * g0
  let f02 = 2 * f0
  let g1 = vtau - g02
  let f1 = vtau - f02
  addG(g0)
  addF(f0)
  for i in 0..<nsteps:
    if i!=0:
      addG(g02)
      addF(f02)
    addT(t0)
    addG(g1)
    addF(f1)
    addT(t0)
  addF(f0)
  addG(g0)

# Order 4

proc setupMDacaca =
  if pt0 == 0: pt0 = 0.211324865
  if pgf0 == 0: pgf0 = 0.022329099
  if pff0 == 0: pff0 = 0.022329099
  let t0 = vtau * pushParam(pt0)
  let t1 = vtau - 2 * t0
  let g0 = 0.5 * vtau
  let f0 = 0.5 * vtau
  let g0g = vtau*vtau*pushParam(pgf0)
  let f0g = if nf==0: g0g else: vtau*vtau*pushParam(pff0)
  addT(t0)
  addGF(g0, g0g)
  addFF(f0, f0g)
  addT(t1)
  addFF(f0, f0g)
  addGF(g0, g0g)
  addT(t0)

proc setupMDabababa =
  if pt0 == 0: pt0 = 0.6756035959798288
  if pg0 == 0: pg0 = 1.351207191959658
  if pf0 == 0: pf0 = 1.351207191959658
  let t0 = vtau * pushParam(pt0)
  let g0 = vtau * pushParam(pg0)
  let f0 = if nf==0: g0 else: vtau * pushParam(pf0)
  let t1 = 0.5 * vtau - t0
  let g1 = vtau - 2 * g0
  let f1 = vtau - 2 * f0
  addT(t0)
  addG(g0)
  addF(f0)
  addT(t1)
  addG(g1)
  addF(f1)
  addT(t1)
  addF(f0)
  addG(g0)
  addT(t0)

proc setupMDababababa =
  if pt0 == 0: pt0 = 0.1720865590295143
  if pg0 == 0: pg0 = 0.5915620307551568
  if pf0 == 0: pf0 = 0.5915620307551568
  if pt1 == 0: pt1 = -0.1616217622107222
  let t0 = vtau * pushParam(pt0)
  let g0 = vtau * pushParam(pg0)
  let f0 = if nf==0: g0 else: vtau * pushParam(pf0)
  let t1 = vtau * pushParam(pt1)
  let g1 = 0.5 * vtau - g0
  let f1 = 0.5 * vtau - f0
  let t2 = vtau - 2*t0 - 2*t1
  addT(t0)
  addG(g0)
  addF(f0)
  addT(t1)
  addG(g1)
  addF(f1)
  addT(t2)
  addF(f1)
  addG(g1)
  addT(t1)
  addF(f0)
  addG(g0)
  addT(t0)

proc setupMDbacab =
  if pg0 == 0: pg0 = 1.0/6.0
  if pf0 == 0: pf0 = 1.0/6.0
  if pgf1 == 0: pgf1 = 1.0/24.0
  if pff1 == 0: pff1 = 1.0/24.0
  let g0 = vtau * pushParam(pg0)
  let f0 = if nf==0: g0 else: vtau * pushParam(pf0)
  let t0 = 0.5 * vtau
  let g1 = vtau - 2 * g0
  let f1 = vtau - 2 * f0
  let g1g = vtau*vtau*pushParam(pgf1)
  let f1g = if nf==0: g1g else: vtau*vtau*pushParam(pff1)
  addG(g0)
  addF(f0)
  addT(t0)
  addGF(g1, g1g)
  addFF(f1, f1g)
  addT(t0)
  addF(f0)
  addG(g0)

proc setupMDabacaba =
  if pt0 == 0: pt0 = 0.08935804763220157
  if pg0 == 0: pg0 = 0.2470939580390842
  if pf0 == 0: pf0 = 0.2470939580390842
  if pgf1 == 0: pgf1 = 0.006938106540706989*(2.0/(1-2*pg0))
  if pff1 == 0: pff1 = 0.006938106540706989*(2.0/(1-2*pf0))
  let t0 = vtau * pushParam(pt0)
  let t02 = 2 * t0
  let g0 = vtau * pushParam(pg0)
  let f0 = if nf==0: g0 else: vtau * pushParam(pf0)
  let t1 = 0.5 * vtau - t0
  let g1 = vtau - 2 * g0
  let f1 = vtau - 2 * f0
  let gf1 = vtau*vtau*pushParam(pgf1)
  let ff1 = if nf==0: gf1 else: vtau*vtau*pushParam(pff1)
  addT(t0)
  for i in 0..<nsteps:
    if i!=0: addT(t02)
    addG(g0)
    addF(f0)
    addT(t1)
    #addGF(g1, gf1)
    #addFF(f1, ff1)
    addGFF(g1, gf1, f1, ff1)
    addT(t1)
    addF(f0)
    addG(g0)
  addT(t0)

proc setupMDacacaca =
  if pt0 == 0: pt0 = 0.116438749543126
  if pg0 == 0: pg0 = 0.283216992495952
  if pf0 == 0: pf0 = 0.283216992495952
  if pgf0 == 0: pgf0 = 0.001247201195115*(2.0/pg0)
  if pff0 == 0: pff0 = 0.001247201195115*(2.0/pf0)
  if pgf1 == 0: pgf1 = 0.002974030329635*(2.0/(1-2*pg0))
  if pff1 == 0: pff1 = 0.002974030329635*(2.0/(1-2*pf0))
  let t0 = vtau * pushParam(pt0)
  let t02 = 2 * t0
  let g0 = vtau * pushParam(pg0)
  let f0 = if nf==0: g0 else: vtau * pushParam(pf0)
  let t1 = 0.5 * vtau - t0
  let g1 = vtau - 2 * g0
  let f1 = vtau - 2 * f0
  let gf0 = vtau*vtau*pushParam(pgf0)
  let ff0 = if nf==0: gf0 else: vtau*vtau*pushParam(pff0)
  let gf1 = vtau*vtau*pushParam(pgf1)
  let ff1 = if nf==0: gf1 else: vtau*vtau*pushParam(pff1)
  addT(t0)
  for i in 0..<nsteps:
    if i!=0: addT(t02)
    #addG(g0)
    #addF(f0)
    addGFF(g0, gf0, f0, ff0)
    addT(t1)
    #addGF(g1, gf1)
    #addFF(f1, ff1)
    addGFF(g1, gf1, f1, ff1)
    addT(t1)
    addGFF(g0, gf0, f0, ff0)
    #addF(f0)
    #addG(g0)
  addT(t0)

# Order 6

proc setupMDacabacabaca =
  if pt0 == 0: pt0 = 0.1097059723948682
  if pg0 == 0: pg0 = 0.2693315848935301
  if pf0 == 0: pf0 = 0.2693315848935301
  if pgf0 == 0: pgf0 = 0.0008642161339706166*(2.0/pg0)
  if pff0 == 0: pff0 = 0.0008642161339706166*(2.0/pf0)
  if pt1 == 0: pt1 = 0.4140632267310831
  if pg1 == 0: pg1 = 1.131980348651556
  if pf1 == 0: pf1 = 1.131980348651556
  if pgf2 == 0: pgf2 = -0.01324638643416052*(2.0/(1-2*(pg0+pg1)))
  if pff2 == 0: pff2 = -0.01324638643416052*(2.0/(1-2*(pf0+pf1)))
  let t0 = vtau * pushParam(pt0)
  let g0 = vtau * pushParam(pg0)
  let f0 = if nf==0: g0 else: vtau * pushParam(pf0)
  let g0f = vtau*vtau*pushParam(pgf0)
  let f0f = if nf==0: g0f else: vtau*vtau*pushParam(pff0)
  let t1 = vtau * pushParam(pt1)
  let g1 = vtau * pushParam(pg1)
  let f1 = if nf==0: g1 else: vtau * pushParam(pf1)
  let t2 = 0.5 * vtau - t0 - t1
  let g2 = vtau - 2*g0 - 2*g1
  let f2 = vtau - 2*f0 - 2*f1
  let g2f = vtau*vtau*pushParam(pgf2)
  let f2f = if nf==0: g2f else: vtau*vtau*pushParam(pff2)
  addT(t0)
  addGF(g0, g0f)
  addFF(f0, f0f)
  addT(t1)
  addG(g1)
  addF(f1)
  addT(t2)
  addGF(g2, g2f)
  addFF(f2, f2f)
  addT(t2)
  addF(f1)
  addG(g1)
  addT(t1)
  addFF(f0, f0f)
  addGF(g0, g0f)
  addT(t0)

# Other

proc setupMD2xn(n: int) =
  if vtau.obj==0.04: vtau.obj = (if nf==0: 0.1*n else: 0.04*n)
  let eps = vtau / n
  let lamt = eps * pushParam(pt0)
  let lamt2 = 2*lamt
  let alp = eps - lamt2
  let eps2 = 0.5*eps
  addT(lamt)
  addG(eps2)
  addF(eps2)
  addT(alp)
  addG(eps2)
  addF(eps2)
  for i in 1..<n:
    addT(lamt2)
    addG(eps2)
    addF(eps2)
    addT(alp)
    addG(eps2)
    addF(eps2)
  addT(lamt)

proc setupMD2pn(n: int) =
  if vtau.obj==0.04: vtau.obj = (if nf==0: 0.1*n else: 0.04*n)
  let eps = vtau / n
  let lamt = eps * pushParam(pt0)
  let lamt2 = 2*lamt
  let alp = eps - lamt2
  let eps2 = 0.5*eps
  addG(lamt)
  addF(lamt)
  addT(eps2)
  addG(alp)
  addF(alp)
  addT(eps2)
  for i in 1..<n:
    addG(lamt2)
    addF(lamt2)
    addT(eps2)
    addG(alp)
    addF(alp)
    addT(eps2)
  addG(lamt)
  addF(lamt)

proc setupMD2pgn(n: int) =
  if vtau.obj==0.04: vtau.obj = (if nf==0: 0.1*n else: 0.04*n)
  if pt0==0: pt0 = 1.0/6.0
  let n2 = (n-1) div 2
  let nr = n - 2*n2
  let eps = vtau / n
  let lamt = eps * pushParam(pt0)
  let lamt2 = 2*lamt
  let alp = eps - lamt2
  let eps2 = 0.5*eps
  let rhot = vtau * vtau * vtau * pushParam(0.2/n)
  addG(lamt)
  addF(lamt)
  addT(eps2)
  for i in 0..<n2:
    addG(alp)
    addF(alp)
    addT(eps2)
    addG(lamt2)
    addF(lamt2)
    addT(eps2)
  if nr == 1:
    addGF(alp, rhot)
    addFF(alp, rhot)
  else:
    addG(alp)
    addF(alp)
    addT(eps2)
    addGF(lamt2, rhot)
    addFF(lamt2, rhot)
    addT(eps2)
    addF(alp)
    addG(alp)
  for i in 0..<n2:
    addT(eps2)
    addF(lamt2)
    addG(lamt2)
    addT(eps2)
    addF(alp)
    addG(alp)
  addT(eps2)
  addF(lamt)
  addG(lamt)

proc setupMD2xn2(n: int) =
  if vtau.obj==0.04: vtau.obj = (if nf==0: 0.1*n else: 0.04*n)
  if pt0==0: pt0 = 1.0/6.0
  let n2 = (n-1) div 2
  let nr = n - 2*n2
  let eps = vtau / n
  var ts = newSeq[FloatV]()
  var gs = newSeq[FloatV]()
  var fs = newSeq[FloatV]()
  let t0 = eps * pushParam(pt0)
  var tc = t0
  addT(t0)
  for i in 0..<n2:
    let g0 = eps * pushParam(0.5)
    let f0 = if nf==0: g0 else: eps * pushParam(0.5)
    let t1 = eps * pushParam(1-2*pt0)
    let g1 = eps * pushParam(0.5)
    let f1 = if nf==0: g1 else: eps * pushParam(0.5)
    let t2 = eps * pushParam(2*pt0)
    ts.add t1
    ts.add t2
    gs.add g0
    gs.add g1
    fs.add f0
    fs.add f1
    tc = tc + t1 + t2
    addG(g0)
    addF(f0)
    addT(t1)
    addG(g1)
    addF(f1)
    addT(t2)
  if nr == 1:
    tc = vtau - 2*tc
    var gc = 0.5*vtau
    var fc = 0.5*vtau
    for i in 0..<gs.len:
      gc = gc - gs[i]
      fc = fc - fs[i]
    addG(gc)
    addF(fc)
    addT(tc)
    addF(fc)
    addG(gc)
  else: # nr == 2
    let g0 = eps * pushParam(0.5)
    let f0 = if nf==0: g0 else: eps * pushParam(0.5)
    let t1 = eps * pushParam(1-2*pt0)
    gs.add g0
    fs.add f0
    tc = tc + t1
    tc = vtau - 2*tc
    var gc = 0.5*vtau
    var fc = 0.5*vtau
    for i in 0..<gs.len:
      gc = gc - gs[i]
      fc = fc - fs[i]
    addG(g0)
    addF(f0)
    addT(t1)
    addG(gc)
    addF(fc)
    addT(tc)
    addF(fc)
    addG(gc)
    addT(t1)
    addF(f0)
    addG(g0)
  for i in countdown(n2-1,0):
    addT(ts[2*i+1])
    addF(fs[2*i+1])
    addG(gs[2*i+1])
    addT(ts[2*i])
    addF(fs[2*i])
    addG(gs[2*i])
  addT(t0)

proc setupMD2pn2(n: int) =
  if vtau.obj==0.04: vtau.obj = (if nf==0: 0.1*n else: 0.04*n)
  let eps = vtau / n
  let n2 = (n-1) div 2
  let nr = n - 2*n2
  var ts = newSeq[FloatV]()
  var gs = newSeq[FloatV]()
  var fs = newSeq[FloatV]()
  let g0 = eps * pushParam(pt0)
  let f0 = if nf==0: g0 else: eps * pushParam(pt0)
  var gc = g0
  var fc = f0
  addG(g0)
  addF(f0)
  for i in 0..<n2:
    let t0 = eps * pushParam(0.5)
    let g1 = eps * pushParam(1-2*pt0)
    let f1 = if nf==0: g1 else: eps * pushParam(1-2*pt0)
    let t1 = eps * pushParam(0.5)
    let g2 = eps * pushParam(2*pt0)
    let f2 = if nf==0: g2 else: eps * pushParam(2*pt0)
    ts.add t0
    ts.add t1
    gs.add g1
    gs.add g2
    fs.add f1
    fs.add f2
    gc = gc + g1 + g2
    fc = fc + f1 + f2
    addT(t0)
    addG(g1)
    addF(f1)
    addT(t1)
    addG(g2)
    addF(f2)
  if nr == 1:
    gc = vtau - 2*gc
    fc = vtau - 2*fc
    var tc = 0.5*vtau
    for i in 0..<ts.len:
      tc = tc - ts[i]
    addT(tc)
    addG(gc)
    addF(fc)
    addT(tc)
  else: # nr == 2
    let t0 = eps * pushParam(0.5)
    let g1 = eps * pushParam(1-2*pt0)
    let f1 = if nf==0: g1 else: eps * pushParam(1-2*pt0)
    ts.add t0
    gc = gc + g1
    gc = vtau - 2*gc
    fc = fc + f1
    fc = vtau - 2*fc
    var tc = 0.5*vtau
    for i in 0..<ts.len:
      tc = tc - ts[i]
    addT(t0)
    addG(g1)
    addF(f1)
    addT(tc)
    addG(gc)
    addF(fc)
    addt(tc)
    addF(f1)
    addG(g1)
    addT(t0)
  for i in countdown(n2-1,0):
    addF(fs[2*i+1])
    addG(gs[2*i+1])
    addT(ts[2*i+1])
    addF(fs[2*i])
    addG(gs[2*i])
    addT(ts[2*i])
  addF(f0)
  addG(g0)

#[
proc setupMD3 =
  vtau.obj = 0.12
  let lamt = vtau * pushParam(0.12)
  let sigt = vtau * pushParam(0.3)
  let sigtf = vtau * pushParam(0.3)
  let alp = 0.5*vtau - lamt
  let bet = vtau - 2*sigt
  let betf = vtau - 2*sigtf
  addT(lamt)
  addG(sigt)
  addF(sigtf)
  addT(alp)
  addG(bet)
  addF(betf)
  addT(alp)
  addG(sigt)
  addF(sigtf)
  addT(lamt)

proc setupMD3g =
  vtau.obj = 0.12
  let lamt = vtau * pushParam(0.11)
  let sigt = vtau * pushParam(0.28)
  let sigtf = vtau * pushParam(0.28)
  let alp = 0.5*vtau - lamt
  let bet = vtau - 2*sigt
  let betf = vtau - 2*sigtf
  let rhot = vtau * vtau * vtau * pushParam(0.12)
  let rhotf = vtau * vtau * vtau * pushParam(0.12)
  addT(lamt)
  addG(sigt)
  addF(sigtf)
  addT(alp)
  addGF(bet, rhot)
  addFF(betf, rhotf)
  addT(alp)
  addG(sigt)
  addF(sigtf)
  addT(lamt)

proc setupMD3p =
  let lam = pushParam(pt0)
  let lamt = lam * vtau
  let sig = pushParam(sig0)
  let sigt = sig * vtau
  let alp = 0.5*vtau - lamt
  let bet = vtau - 2*sigt
  addG(lamt)
  addT(sigt)
  addG(alp)
  addT(bet)
  addG(alp)
  addT(sigt)
  addG(lamt)

# 2x OMF2: 0.1 0.25 0.3
proc setupMD4o =
  if vtau.obj==0.04: vtau.obj = 0.08
  var a0 = if pt0==0.19: 0.1 else: pt0
  let lamt = vtau * pushParam(a0)
  let lam2t = 2 * lamt
  let sigt = 0.5*vtau - lam2t
  let tau4 = 0.25 * vtau
  addT(lamt)
  addG(tau4)
  addF(tau4)
  addT(sigt)
  addG(tau4)
  addF(tau4)
  addT(lam2t)
  addG(tau4)
  addF(tau4)
  addT(sigt)
  addG(tau4)
  addF(tau4)
  addT(lamt)

proc setupMD4 =
  if vtau.obj==0.04: vtau.obj = 0.08
  var a0 = if pt0==0.19: 0.1 else: pt0
  let lamt = vtau * pushParam(a0)
  let sigt = vtau * pushParam(0.25)
  let sigtf = vtau * pushParam(0.25)
  let rhot = vtau * pushParam(0.3)
  #let lamt = vtau * pushParam(0.1786178958448091)
  #let sigt = vtau * pushParam(0.7123418310626056)
  #let sigtf = vtau * pushParam(0.7123418310626056)
  #let rhot = vtau * pushParam(-0.06626458266981843)
  let alp = 0.5*vtau - sigt
  let bet = vtau - 2*lamt - 2*rhot
  #let sigtf = sigt
  #let alpf = alp
  let alpf = 0.5*vtau - sigtf
  addT(lamt)
  addG(sigt)
  addF(sigtf)
  addT(rhot)
  addG(alp)
  addF(alpf)
  addT(bet)
  addG(alp)
  addF(alpf)
  addT(rhot)
  addG(sigt)
  addF(sigtf)
  addT(lamt)

proc setupMD4p1 =
  if vtau.obj==0.04: vtau.obj = 0.08
  var a0 = if pt0==0.19: 0.1 else: pt0
  let lamt = vtau * pushParam(a0)
  let lam2t = 2 * lamt
  let alp = 0.5*vtau - lam2t
  let tau4 = 0.25 * vtau
  addG(lamt)
  addF(lamt)
  addT(tau4)
  addG(alp)
  addF(alp)
  addT(tau4)
  addG(lam2t)
  addF(lam2t)
  addT(tau4)
  addG(alp)
  addF(alp)
  addT(tau4)
  addG(lamt)
  addF(lamt)

proc setupMD4p2 =
  if vtau.obj==0.04: vtau.obj = 0.08
  var a0 = if pt0==0.19: 0.1 else: pt0
  let g0 = vtau * pushParam(a0)
  let f0 = vtau * pushParam(a0)
  let t0 = vtau * pushParam(0.25)
  let g1 = vtau * pushParam(0.3)
  let f1 = vtau * pushParam(0.3)
  let t1 = 0.5*vtau - t0
  let g2 = vtau - 2*g0 - 2*g1
  let f2 = vtau - 2*f0 - 2*f1
  addG(g0)
  addF(f0)
  addT(t0)
  addG(g1)
  addF(f1)
  addT(t1)
  addG(g2)
  addF(f2)
  addT(t1)
  addG(g1)
  addF(f1)
  addT(t0)
  addG(g0)
  addF(f0)

proc setupMD4g =
  let lam = pushParam(pt0)
  let lamt = lam * vtau
  let sig = pushParam(sig0)
  let sigt = sig * vtau
  let rho = pushParam(rho0)
  let rhot = rho * vtau
  let alp = 0.5*vtau - sigt
  let bet = vtau - 2*lamt - 2*rhot
  let theta = pushParam(theta0)
  let thetat = theta * vtau * vtau * vtau
  addG(lamt)
  addT(sigt)
  addG(rhot)
  addT(alp)
  addGF(bet, thetat)
  addT(alp)
  addG(rhot)
  addT(sigt)
  addG(lamt)

proc setupMD5 =
  vtau.obj = 0.12
  let lamt = vtau * pushParam(0.11)
  let sigt = vtau * pushParam(0.27)
  let sigtf = vtau * pushParam(0.27)
  let rhot = vtau * pushParam(0.27)
  let thetat = vtau * pushParam(0.17)
  let thetatf = vtau * pushParam(0.17)
  let alp = 0.5*vtau - lamt - rhot
  let bet = vtau - 2*sigt - 2*thetat
  let betf = vtau - 2*sigtf - 2*thetatf
  addT(lamt)
  addG(sigt)
  addF(sigtf)
  addT(rhot)
  addG(thetat)
  addF(thetatf)
  addT(alp)
  addG(bet)
  addF(betf)
  addT(alp)
  addG(thetat)
  addF(thetatf)
  addT(rhot)
  addG(sigt)
  addF(sigtf)
  addT(lamt)

# 6x LF: 0.08 0.17 0.17 0.17 0.17
# 3x OMF2: 0.06 0.17 0.2 0.17 0.13
proc setupMD6 =
  vtau.obj = 0.14
  let lamt = vtau * pushParam(0.06)
  let sigt = vtau * pushParam(0.17)
  let sigtf = vtau * pushParam(0.17)
  let rhot = vtau * pushParam(0.2)
  let thetat = vtau * pushParam(0.17)
  let thetatf = vtau * pushParam(0.17)
  let upst = vtau * pushParam(0.13)
  let alp = 0.5*vtau - sigt - thetat
  let alpf = 0.5*vtau - sigtf - thetatf
  let bet = vtau - 2*lamt - 2*rhot - 2*upst
  addT(lamt)
  addG(sigt)
  addF(sigtf)
  addT(rhot)
  addG(thetat)
  addF(thetatf)
  addT(upst)
  addG(alp)
  addF(alpf)
  addT(bet)
  addG(alp)
  addF(alpf)
  addT(upst)
  addG(thetat)
  addF(thetatf)
  addT(rhot)
  addG(sigt)
  addF(sigtf)
  addT(lamt)

# G T GF T G T GF T G
proc setupMD52p =
  let g0 = vtau * pushParam(0.5*pt0)
  let g1 = vtau * pushParam(0.5*sig0)
  let g2 = vtau - 2 * g0 - 2 * g1
  let t0 = vtau * pushParam(rho0)
  let t1 = 0.5 * vtau - t0
  let f0 = 0.5 * vtau
  addG(g0)
  addT(t0)
  addG(g1)
  addF(f0)
  addT(t1)
  addG(g2)
  addT(t1)
  addG(g1)
  addF(f0)
  addT(t0)
  addG(g0)
]#

proc updatePseudo =
  tape.setTrack 3
  tape.run
  tape.setTrack 0

proc update =
  tic("update")
  tape.run
  for mu in 0..<g.len:
    g[mu] := gauges[^1][mu]
    p[mu] := moms[^1][mu]
  toc("end update")

proc updateAll =
  updatePseudo()
  update()

var
  p2v = newSeq[type tape.newFloatV](0)
  gav = newSeq[type tape.newFloatV](0)
  fav = newSeq[type tape.newFloatV](0)
  hv = newSeq[type tape.newFloatV](0)
proc addAction(p: GaugeV, g: GaugeV) =
  pushTemp()
  norm2subtract(ptemps[^1], p, 8.0)
  pushTemp()
  mul(ptemps[^1], 0.5, ptemps[^2])
  p2v.add ptemps[^1]
  pushTemp()
  gaction(gc, ptemps[^1], g)
  gav.add ptemps[^1]
  if nf == 0:
    hv.add p2v[^1] + gav[^1]
  else:
    pushTemp()
    let fasum = ptemps[^1]
    fasum := 0
    for i in 0..<phi.len-1:
      pushCvec()
      stag.agradD(g, cvecvs[^1], phiv[i], vhmasses[i])
      pushCvec()
      if i == 0:
        stag.agradSolve(g, cvecvs[^1], cvecvs[^2], mass, spa)
      else:
        stag.agradSolve(g, cvecvs[^1], cvecvs[^2], vhmasses[i-1], spa)
      pushTemp()
      norm2subtract(ptemps[^1], cvecvs[^1], 3.0)
      fasum += ptemps[^1]
    pushCvec()
    if phi.len == 1:
      stag.agradSolve(g, cvecvs[^1], phiv[^1], mass, spa)
    else:
      stag.agradSolve(g, cvecvs[^1], phiv[^1], vhmasses[^1], spa)
    pushTemp()
    norm2subtract(ptemps[^1], cvecvs[^1], 3.0)
    fasum += ptemps[^1]
    fav.add 0.5*fasum
    hv.add p2v[^1] + gav[^1] + fav[^1]

proc addPseudo(g: GaugeV) =
  for i in 0..<hmasses.len:
    pushCvec()
    if i == 0:
      stag.agradD(g, cvecvs[^1], eta[i], mass)
    else:
      stag.agradD(g, cvecvs[^1], eta[i], vhmasses[i-1])
    stag.agradSolve(g, phiv[i], cvecvs[^1], vhmasses[i], spa)
    maskEven phiv[i]
  if vhmasses.len == 0:
    stag.agradD(g, phiv[^1], eta[^1], mass)
  else:
    stag.agradD(g, phiv[^1], eta[^1], vhmasses[^1])
  maskEven phiv[^1]

proc setupAction =
  tape.addTrack
  addAction(momvs[0], gaugevs[0])
  tape.addTrack
  addAction(momvs[^1], gaugevs[^1])
  tape.addTrack
  addPseudo(gaugevs[0])
  tape.setTrack 0
  #echo tape

type
  Met = ref object of MetropolisRoot
    state: int

proc init(m: var Met) =
  m.new
  var r = MetropolisRoot m
  init(r)
  m.verbosity = 1

template masses(i: int): float =
  if i==0: mass
  else: vhmasses[i-1].obj

template masses(bi: BackwardsIndex): float =
  let i = 1 + hmasses.len - int(bi)
  if i==0: mass
  else: vhmasses[i-1].obj

proc start*(m: var Met) =
  tic("start")
  m.state = 0
  threads:
    for i in 0..<g.len:
      g0[i] := g[i]
    p.randomTAH r
    threadBarrier()
    for i in 0..<p.len:
      p0[i] := p[i]
  if nf != 0:
    threads:
      for i in 0..<eta.len:
        eta[i].gaussian r
    updatePseudo()
  toc("init p, phi")

proc getH*(m: Met): float =
  tic("getH")
  var p2 = 0.0
  threads:
    var p2t = 0.0
    for i in 0..<p.len:
      p2t += p[i].norm2subtract(8.0)
    threadMaster: p2 = p2t
  toc("p2")
  var f2 = 0.0
  if nf != 0:
    threads:
      stag.rephase
    for i in 0..<(phi.len-1):
      threads:
        stag.D(ftmp, phi[i], masses(i+1))
        #ftmp := phi[i]
      stag.solve(psi[i], ftmp, masses(i), spa)
    stag.solve(psi[^1], phi[^1], masses(^1), spa)
    threads:
      for i in 0..<psi.len:
        var psi2 = psi[i].norm2subtract(3.0)
        threadMaster: f2 += psi2
      stag.rephase
  toc("fa solve")
  let
    ga0 = gc.actionA g
    fa0 = 0.5*f2
    t0 = 0.5*p2
    #h0 = t0 + ga0
    h0 = t0 + ga0 + fa0
  result = h0
  if m.state == 0:
    inc m.state
    echo &"Begin H: {h0}  T: {t0}  Sg: {ga0}  Sf: {fa0}"
    tape.setTrack 1
    tape.run
    echo &"      H: {hv[0].obj}  T: {p2v[0].obj}  Sg: {gav[0].obj}  Sf: {fav[0].obj}"
    tape.setTrack 0
  else:
    echo &"End H: {h0}  T: {t0}  Sg: {ga0}  Sf: {fa0}"
    tape.setTrack 2
    tape.run
    echo &"    H: {hv[1].obj}  T: {p2v[1].obj}  Sg: {gav[1].obj}  Sf: {fav[1].obj}"
    tape.setTrack 0
  toc("end getH")

# w = sum_i w_i
# sum = sum_i w_i x_i
# sum2 = sum_i w_i x_i^2
# X = sum/w
# var = sum_i w_i (x_i - X)^2 = sum_i w_i x_i^2 - w X^2
type DecayStat = object
  n: int
  fac: float
  fac2: float
  w: float
  sum: float
  w2: float
  sum2: float
proc clear(d: var DecayStat) =
  d.n = 0
  d.w = 0
  d.sum = 0
  d.w2 = 0
  d.sum2 = 0
proc newDecayStat(f = 0.9, f2 = 0.999): DecayStat =
  result.fac = f
  result.fac2 = f2
proc push(d: var DecayStat, x: float) =
  inc d.n
  let f = d.fac
  let o = 1.0 - f
  d.w = f*d.w + o
  d.sum = f*d.sum + o*x
  let f2 = d.fac2
  let o2 = 1.0 - f2
  d.w2 = f2*d.w2 + o2
  d.sum2 = f2*d.sum2 + o2*x*x
proc mean(d: DecayStat): float =
  result = d.sum / d.w
proc rms(d: DecayStat): float =
  result = sqrt( d.sum2 / d.w2 )
proc variance(d: DecayStat): float =
  let m = d.sum / d.w
  let m2 = d.sum2 / d.w2
  result = sqrt(m2 - m*m)

var pacc = newDecayStat()
var paccg = newSeq[DecayStat]()

#[
proc getCost1(m: Met): seq[float] =
  let c0 = 0.0
  result.newSeq(0)
  let ct = c0 + vtau.obj
  #var cost = m.avgPAccept
  #var cost = ct * m.avgPAccept
  #var cost = ct * ct * m.avgPAccept
  #result.add cost
  var cost = 1.0/(ct * ct * m.avgPAccept)
  result.add nff*cost
  for i in 0..<params.len:
    #var costg = 0.0
    #if m.hNew > m.hOld:
    #  costg -= m.pAccept * params[i].grad
    #var costg = (if i==0: 1.0 else: 0.0)
    var costg = (if i==0: 2.0*ct else: 0.0)
    if m.hNew > m.hOld:
      #costg = m.pAccept * (costg - ct * params[i].grad)
      costg = m.pAccept * (costg - ct * ct * params[i].grad)
    costg = costg/(cost*cost)  # extra - to make it minimize
    #if i==0: costg = 0
    #result.add costg
    result.add nff*costg
]#

proc forceCost: float =
  result = nff[0].float
  for i in 0..<vhmasses.len:
    result += nff[i+1].float * (mass/vhmasses[i].obj)

proc forceCostGrad(k: int): float =
  let i = k - 1
  if i>=0 and i<vhmasses.len:
    result = -nff[i+1].float*mass/(vhmasses[i].obj*vhmasses[i].obj)

proc getCost0(m: Met): float =
  let ct = nsteps * vtau.obj
  forceCost()/(ct * ct * m.avgPAccept)
  #nff/(ct * ct * pacc.mean)

proc getCost(m: Met): seq[float] =
  result.newSeq(0)
  #let pm = pacc.mean
  let pm = m.pAccept
  let ct = nsteps * vtau.obj
  #result.add (ct*ct*pm)/nff
  let fc = forceCost()
  result.add fc/(ct*ct*pm)
  #var cost = 1.0/()
  #result.add nff*cost
  for i in 0..<params.len:
    if i >= paccg.len: paccg.add newDecayStat()
    # grad (p = min(1,exp(ho-hn))) -> 0 or - p grad(hn)
    var pg = 0.0
    if m.hNew > m.hOld:
      pg = - pm * params[i].grad
    paccg[i].push pg
    #pg = paccg[i].mean
    var costg = (if i==0: 2.0*nsteps*ct*pm else: 0.0)
    costg = costg + ct * ct * pg
    let d = m.hNew - m.hOld
    #let alp = alpha
    #costg += alp*d*(d*pg + 2*(pm-1)*params[i].grad)
    costg = costg/fc - (ct*ct*pm/(fc*fc))*forceCostGrad(i)
    #costg = nff*costg*(cost*cost)  # extra - to make it minimize
    if i==0:
      if m.deltaH > 10:
        costg = -0.1
    if fixtau and i==0: costg = 0
    if fixhmasses and i>0 and i<=hmasses.len: costg = 0
    if fixparams and i>hmasses.len: costg = 0
    #if i > 0:
    #  costg = (m.hOld-m.hNew)*params[i].grad
    result.add costg

proc checkGrad(m: Met) =
  var gx = g.newOneOf
  gx := g
  #let tg = vtau.grad
  #let f0 = vtau.obj * exp(m.hOld-m.hNew)
  #let h0 = m.hOld
  let eps = 1e-6
  var gs = newSeq[float](0)
  for i in 0..<params.len: gs.add params[i].grad
  for i in 0..<params.len:
    let t = params[i].obj
    params[i].obj += eps
    updatePseudo()
    threads:
      for mu in 0..<g.len:
        g[mu] := g0[mu]
        p[mu] := p0[mu]
    let h0 = m.getH
    update()
    let h1 = m.getH
    let dh = h1-h0
    echo "params[",i,"] grad: ", gs[i]
    echo "params[",i,"] diff: ", (dh-m.deltaH)/eps
    params[i].obj = t
  #update()
  g := gx
  for i in 0..<params.len: params[i].grad = gs[i]

#var cgstat = newSeq[RunningStat](0)
var cgstat = newSeq[DecayStat](0)

proc getGrad(m: Met) =
  hv[0].grad = -1.0
  hv[1].grad = 1.0
  tape.setTrack 1
  tape.grad
  tape.setTrack 2
  tape.grad
  tape.setTrack 0
  #gc.gaugeDeriv2(gaugevs[^1].obj, gaugevs[^1].grad)
  #for mu in 0..<g.len:
  #  momvs[^1].grad[mu] := momvs[^1].obj[mu]
  tape.grad
  tape.setTrack 3
  tape.grad
  if checkg:
    checkGrad(m)
  #let tg = vtau.grad
  let cst = getCost(m)
  echo "cost:  ", cst[0]
  #if cgstat.len < params.len: cgstat.setLen(params.len)
  if cgstat.len < params.len:
    for i in 0..<params.len: cgstat.add newDecayStat()
  for i in 0..<cgstat.len:
    cgstat[i].push cst[i+1]
    let m = cgstat[i].mean
    let v = cgstat[i].rms
    #echo "costg: ", cgstat[i].mean, " ", cgstat[i].standardDeviationS
    #echo "costg: ", cgstat[i].mean, " ", cst[i+1]
    echo &"costg: {m/v:8.6f} {m:10.6f} {cst[i+1]}"

proc updateParams(ratefac: float) =
  let eps = 1e-8
  let n = cgstat.len
  #var s = 0.0
  var p = newSeq[float](n)
  var g = newSeq[float](n)
  for i in 0..<n:
    let m = cgstat[i].mean
    let v = cgstat[i].rms
    let t = m / (v+eps)
    #echo "t: ", t
    g[i] = t
    #s += t*t
    p[i] = params[i].obj
  echo "Params: ", p
  #s = rate*sqrt(n/s)
  for i in 0..<n:
    var rate = ratefac * lrate
    if i>0 and i<=hmasses.len:
      rate = ratefac * lrateh
    #let m = cgstat[i].mean
    #let d = s*g[i]
    let d = rate*g[i]
    params[i].obj += d
    p[i] = params[i].obj
  echo "Params: ", p
  #for i in 0..<n:
  #  clear cgstat[i]

proc finish*(m: var Met) =
  discard
  #for mu in 0..<g.len:
  #  g[mu] := g0[mu]
  #  p[mu] := p0[mu]
  #discard m.getH
  #mdt(0.5*tau)
  #echo "mdt: ", g.norm2, "  ", vg2.obj.norm2
  #mdv(tau)
  #echo "mdv: ", p.norm2, "  ", vp2.obj.norm2
  #mdt(0.5*tau)
  #echo "mdt: ", g.norm2, "  ", vg6.obj.norm2
  # pacc = max(1, exp(hOld-hNew))
  #if m.hNew <= m.hold:
  #gc.gaugeDeriv2(gaugevs[^1].obj, gaugevs[^1].grad)
  #for mu in 0..<g.len:
  #  momvs[^1].grad[mu] := momvs[^1].obj[mu]
  #tape.grad
  ##let tg = vtau.grad
  #let cst = getCost(m)
  #echo "cost:  ", cst[0]
  ##if cgstat.len < params.len: cgstat.setLen(params.len)
  #if cgstat.len < params.len:
  #  for i in 0..<params.len: cgstat.add newDecayStat(decay)
  #for i in 0..<cgstat.len:
  #  cgstat[i].push cst[i+1]
  #  #echo "costg: ", cgstat[i].mean, " ", cgstat[i].standardDeviationS
  #  echo "costg: ", cgstat[i].mean, " ", cst[i+1]
  #checkGrad(m)

proc disp(m: Met) =
  let p = g.plaq
  echo "plaq: ", 6.0*p
  echo "tplaq: ", p.sum
  #echo "gauge force: ", sqrt(gf2s/ngf2)
  #gf2s = 0.0
  #ngf2 = 0
  echo "rmsDeltaH: ", sqrt(m.avgDeltaH2)
  echo "avgPAccept: ", m.avgPAccept
  echo "pacc.mean:  ", pacc.mean

proc accept*(m: var Met) =
  pacc.push m.pAccept
  disp(m)

proc reject*(m: var Met) =
  threads:
    for i in 0..<g.len:
      g[i] := g0[i]
  pacc.push m.pAccept
  disp(m)

var alwaysAccept = false
proc globalRand*(m: Met): float =
  if not alwaysAccept:
    result = R.uniform

proc generate*(m: var Met) =
  update()

var m: Met
m.init

var src = lo.ColorVector
var prop = lo.ColorVector
let nt = lat[3]
let nt2 = (nt+2) div 2
#var picorr = newSeq[float](lat[3])
var picorr = newSeq[RunningStat](nt)
var picorrf = newSeq[RunningStat](nt2)
var sp = newSolverParams()
sp.r2req = 1e-20
block:
  var v: evalType(src{0})
  v := 0
  v[0] = 1
  threads:
    src.wallSource(0, v)
  #echo src.norm2slice(3)

var gfStats: RunningStat
var ffStats = newSeq[RunningStat](1+hmasses.len)

proc resetMeasure =
  gfStats.clear
  for i in 0..<ffStats.len:
    ffStats[i].clear

proc avgabs(x: seq[FloatV]): float =
  for i in 0..<x.len:
    result += abs(x[i].obj)
  result /= x.len.float

proc fstats(s: RunningStat): string =
  let f2 = sqrt(s.variance)
  let f4 = f2*sqrt(sqrt((s.kurtosis+3.0)))
  result = &"{s.mean:8.6f}  max {s.max:8.6f}  std {f2:8.6f}  m4 {f4:8.6f}"

proc measure =
  #threads:
  #  g.rephase
  #stag.solve(prop, src, mass, sp)
  #threads:
  #  g.rephase
  #picorr += prop.norm2slice(3)
  #echo picorr
  #let pic = prop.norm2slice(3)
  #for i in 0..<nt:
  #  picorr[i].push pic[i]
  #  let k = min(i, nt-i)
  #  picorrf[k].push pic[i]
  #for i in 0..<nt:
  #  echo i, " ", picorr[i].mean, " ", picorr[i].standardDeviationS
  #for i in 0..<nt2:
  #  let k = nt2 - 1 - i
  #  echo i, " ", picorr[k].mean, " ", picorr[k].standardDeviationS
  #discard
  #let sf = 1.0/(lat.len.float*lo.physVol.float*8.0)  # 8 SU3 generators
  var gfs: RunningStat
  for i in 0..<gfList.len:
    #let f2 = momvs[gfList[i]].obj.norm2
    #gfStats.push gfeList[i].obj*sqrt(sf*f2)
    #let gfs = momvs[gfList[i]].obj.rmsStats gfeList[i].obj
    #let t = momvs[gfList[i]].obj.rmsStats (1.5/gfList.len)  # 1.5=sqrt(18/8)
    let t = momvs[gfList[i]].obj.rmsStats (1.5*gfeList.avgabs)  # 1.5=sqrt(18/8)
    gfs += t
  echo "GF: ", fstats(gfs)
  gfStats += gfs
  for i in 0..<ffList.len:
    var ffs: RunningStat
    for j in 0..<ffList[i].len:
      #let f2 = momvs[ffList[i][j]].obj.norm2
      #ffStats[i].push abs(ffeList[i][j].obj)*sqrt(sf*f2)
      #let ffs = momvs[ffList[i][j]].obj.rmsStats abs(ffeList[i][j].obj)
      #let t = momvs[ffList[i][j]].obj.rmsStats (1.5/ffList[i].len)  # 1.5=sqrt(18/8)
      let t = momvs[ffList[i][j]].obj.rmsStats (1.5*ffeList[i].avgabs)  # 1.5=sqrt(18/8)
      ffs += t
    echo "FF", i, ": ", fstats(ffs)
    ffStats[i] += ffs
  echo "TGF: ", fstats(gfStats)
  for i in 0..<ffList.len:
    echo "TFF", i, ": ", fstats(ffStats[i])

case md
of "aba": setupMDaba()
of "aca": setupMDaca()
of "ababa": setupMDababa()
of "acaca": setupMDacaca()
of "babab": setupMDbabab()
of "bacab": setupMDbacab()
of "abababa": setupMDabababa()
of "abacaba": setupMDabacaba()
of "acacaca": setupMDacacaca()
of "ababababa": setupMDababababa()
of "acabacabaca": setupMDacabacabaca()
of "g5f2": setupMDg5f2()
of "g5f23": setupMDg5f23()
of "g6f24": setupMDg6f24()
of "g10f2": setupMDg10f2()
of "g10f28": setupMDg10f28()
else:
  echo "unknown MD string: ", md
  qexAbort()

#setupMD2()
#setupMD2p()
#setupMD2p2()
#setupMD2g()
#setupMD3()
#setupMD4o()
#setupMD4()
#setupMD4p1()
#setupMD4p2()
#setupMD5()
#setupMD6()

#setupMD2xn(nsteps)
#setupMD2pn(nsteps)
#setupMD2pgn(nsteps)
#setupMD2xn2(nsteps)
#setupMD2pn2(nsteps)

#setupMD1g()
#setupMD3g()
#setupMD4g()
#setupMD52p()

setupAction()
echo "nFF: ", nff
echo "initial force cost: ", forceCost()

block:
  tic("warmup")
  if nwarm > 0:
    echo "Starting warmups"
    #setupMDx()
    alwaysAccept = warmmd
    for n in 1..nwarm:
      m.update
    m.clearStats
    pacc.clear
  toc("end warmup")

echo "Starting HMC"
#setupMD5()
alwaysAccept = false
#gutime = 0.0
#gftime = 0.0
#fftime = 0.0
for i in 0..<spf.len:
  spf[i].resetStats
  spfb[i].resetStats
block:
  tic("training")
  for n in 1..ntrain:
    echo "Starting trajectory: ", n
    echoParams()
    tic("train step")
    m.update
    getGrad(m)
    let tup = getElapsedTime()
    for i in 0..<spf.len:
      echo &"FFits{i}: ", spf[i].getAveStats
    measure()
    if upit > 0:
      if n mod upit == 0:
        updateParams(sqrt(float upit))
        lrate *= anneal
        lrateh *= anneal
    let ttot = getElapsedTime()
    echo "End trajectory update: ", tup, "  measure: ", ttot-tup, "  total: ", ttot
  let et = getElapsedTime()
  toc("end training")
  echo "HMC time: ", et
  #let at = gutime + gftime + fftime
  #echo &"gu: {gutime}  gf: {gftime}  ff: {fftime}  ot: {et-at}  tt: {et}"

resetMeasure()
for i in 0..<spf.len:
  spf[i].resetStats
  spfb[i].resetStats
if trajs > 0:
  m.clearStats
  pacc.clear
  tic("inference")
  for n in 1..trajs:
    echo "Starting inference: ", n
    echoParams()
    tic()
    m.update
    #echo "cost: ", vtau.obj*m.avgPAccept
    #echo "cost: ", vtau.obj*vtau.obj*m.avgPAccept
    #echo "cost: ", nff/(vtau.obj*vtau.obj*m.avgPAccept)
    echo "cost: ", getCost0(m)
    let tup = getElapsedTime()
    for i in 0..<spf.len:
      echo &"FFits{i}: ", spf[i].getAveStats
    measure()
    let ttot = getElapsedTime()
    echo "End inference update: ", tup, "  measure: ", ttot-tup, "  total: ", ttot
  let et = getElapsedTime()
  toc("end inference")
  echo "Inference time: ", et

if outfn != "":
  echo "Saving gauge field to file: ", outfn
  let err = g.saveGauge outfn

echoProf()
qexFinalize()
