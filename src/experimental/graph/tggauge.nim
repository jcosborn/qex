import math, unittest

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))

# basicOps.epsilon collide with fenv.epsilon
import qex except epsilon
import algorithms/numdiff, gauge/stoutsmear
import core, scalar, gauge

template checkeq(ii: tuple[filename:string, line:int, column:int], sa: string, a: float, sb: string, b: float) =
  if not almostEqual(a, b, unitsInLastPlace = 64):
    checkpoint(ii.filename & ":" & $ii.line & ":" & $ii.column & ": Check failed: " & sa & " :~ " & sb)
    checkpoint("  " & sa & ": " & $a)
    checkpoint("  " & sb & ": " & $b)
    fail()

template `:~`(a:Gvalue, b:Gvalue) =
  checkeq(instantiationInfo(), astToStr a, a.eval.getfloat, astToStr b, b.eval.getfloat)

template `:<`(a:Gvalue, b:float) =
  let av = a.eval.getfloat.abs
  if av >= b:
    let ii = instantiationInfo()
    let sa = astToStr a
    let sb = astToStr b
    checkpoint(ii.filename & ":" & $ii.line & ":" & $ii.column & ": Check failed: " & sa & " :< " & sb)
    checkpoint("  " & sa & ": " & $av)
    checkpoint("  " & sb & ": " & $b)
    fail()

# basic test: y <- f(x), or z = y B† = f(x) B†, with x = x + t A
# d/dt z = d/dt f(x+tA) B†
# d/dt z = (d/dt y) (d/dy z)† = (d/dt x) (d/dx z)† = (d/dx z) A†

proc ndiff(zt: Gvalue, t: Gscalar): (float, float) =
  proc z(v:float):float =
    t.update v
    zt.eval.getfloat
  var dzdt,e: float
  ndiff(dzdt, e, z, 0.0, 0.125, ordMax=3)
  (dzdt, e)

template check(ii: tuple[filename:string, line:int, column:int], ast: string, dzdt, e, gdota: float) =
  if not almostEqual(gdota, dzdt, unitsInLastPlace = 512*1024):
    checkpoint(ii.filename & ":" & $ii.line & ":" & $ii.column & ": Check failed: " & ast)
    checkpoint("  ndiff: " & $dzdt & " +/- " & $e)
    checkpoint("  grad: " & $gdota)
    checkpoint("  reldelta: " & $(abs(dzdt-gdota)/abs(dzdt+gdota)))
    fail()

template ckforce(s: untyped, f: untyped, x: Gvalue, p: Gvalue) =
  let t = Gscalar()
  let (dsdt, e) = ndiff(s(exp(t*p)*x), t)
  let pdotf = eval(redot(p, f(x))).getfloat
  check(instantiationInfo(), astTostr(s(x) -> f(x)), dsdt, e, pdotf)

template ckgrad(f: untyped, x: Gvalue, a: Gvalue) =
  let t = Gscalar()
  let (dzdt, e) = ndiff(f(x+t*a), t)
  let ff = f(x)
  let gdota = eval(redot(grad(ff, x), a)).getfloat
  check(instantiationInfo(), astTostr(f(x)), dzdt, e, gdota)

template ckgrad2(f: untyped, x: Gvalue, y: Gvalue, ax: Gvalue, ay: Gvalue) =
  let t = Gscalar()
  let (dzdt, e) = ndiff(f(x+t*ax, y+t*ay), t)
  let ff = f(x, y)
  let gdota = eval(redot(grad(ff, x), ax) + redot(grad(ff, y), ay)).getfloat
  check(instantiationInfo(), astTostr(f(x,y)), dzdt, e, gdota)

template ckgradm(f: untyped, x: Gvalue, a: Gvalue, b: Gvalue) =
  let t = Gscalar()
  let (dzdt, e) = ndiff(f(x+t*a).redot b, t)
  let ff = f(x).redot b
  let gdota = eval(redot(grad(ff, x), a)).getfloat
  check(instantiationInfo(), astTostr(f(x)), dzdt, e, gdota)

template ckgradm2(f: untyped, x: Gvalue, y: Gvalue, ax: Gvalue, ay: Gvalue, b: Gvalue) =
  let t = Gscalar()
  let (dzdt, e) = ndiff(f(x+t*ax, y+t*ay).redot b, t)
  let ff = f(x, y).redot b
  let gdota = eval(redot(grad(ff, x), ax) + redot(grad(ff, y), ay)).getfloat
  check(instantiationInfo(), astTostr(f(x,y)), dzdt, e, gdota)

template ckgradm3(f: untyped, x: Gvalue, y: Gvalue, u: Gvalue, ax: Gvalue, ay: Gvalue, au: Gvalue, b: Gvalue) =
  let t = Gscalar()
  let (dzdt, e) = ndiff(f(x+t*ax, y+t*ay, u+t*au).redot b, t)
  let ff = f(x, y, u).redot b
  let gdota = eval(redot(grad(ff, x), ax) + redot(grad(ff, y), ay) + redot(grad(ff, u), au)).getfloat
  check(instantiationInfo(), astToStr(f(x,y,u)), dzdt, e, gdota)

qexInit()

let
  lat = @[8,8,8,16]
  lo = lat.newLayout
  seed = 1234567891u64
  vol = lo.physVol
var
  r = lo.newRNGField(MRG32k3a, seed)
  g = lo.newgauge
  u = lo.newgauge
  p = lo.newgauge
  q = lo.newgauge
  m = lo.newgauge
  ss = lo.newStoutSmear(0.1)
const nc = g[0][0].nrows
threads:
  g.random r
  u.random r
  p.randomTAH r
  q.randomTAH r
  m.randomTAH r
for i in 0..4:
  ss.smear(g, g)
  ss.smear(u, u)
threads:
  for t in m:
    t *= 0.01

let a = 0.5 * (sqrt(5.0) - 1.0)
let b = sqrt(2.0) - 1.0

suite "gauge basic":
  setup:
    let gg {.used.} = toGvalue g
    let gu {.used.} = toGvalue u
    let gp {.used.} = toGvalue p
    let gq {.used.} = toGvalue q
    let gm {.used.} = toGvalue m
    let x {.used.} = toGvalue a
    let y {.used.} = toGvalue b

  test "norm2":
    let n2 = gg.norm2
    let p2 = gp.norm2
    let dp = grad(0.5 * p2, gp)
    n2 :~ 4.0*float(nc*vol)
    dp.norm2 :~ p2
    norm2(dp-gp) :~ 0
    ckgrad(norm2, gm, gq)

  test "redot":
    let n2 = gg.redot gg
    let p2 = gp.redot gp
    let dp = grad(0.5 * p2, gp)
    n2.eval :~ 4.0*float(nc*vol)
    dp.norm2 :~ p2
    norm2(dp-gp).eval :~ 0
    let pq = gp.redot gq
    norm2(grad(pq, gp) - gq) :< 1e-26
    norm2(grad(pq, gq) - gp) :< 1e-26
    ckgrad2(redot, gp, gq, gg, gu)

  test "retr":
    let rtp = gp.retr
    let n2 = retr(gg * gg.adj)
    rtp*rtp :< 1e-20
    n2.eval :~ 4.0*float(nc*vol)
    let p2 = retr(gp * gq.adj)
    p2 :~ redot(gp, gq)
    norm2(grad(p2, gp) - gq) :< 1e-26
    norm2(grad(p2, gq) - gp) :< 1e-26
    ckgrad(retr, gp, gq)

  test "adj":
    norm2(gg.adj*gg - 1.0)/float(4*nc*vol) :< 1e-22
    norm2(gg*gg.adj - 1.0)/float(4*nc*vol) :< 1e-22
    norm2(gp.adj + gp) :< 1e-26
    norm2(grad(gp.adj.norm2, gp) - 2.0*gp) :< 1e-26
    ckgradm(adj, gg, gp, gq)

  test "neg":
    norm2(gp.adj - (-gp)) :< 1e-26
    norm2(-gp) :~ gp.norm2
    ckgradm(`-`, gg, gp, gq)

  test "addsg":
    let p2 = norm2(x+gp)
    grad(p2, x) :~ retr(2.0*(a+gp))
    norm2(grad(p2, gp) - 2.0*(a+gp)) :< 1e-26
    ckgradm2(`+`, x, gp, y, gq, gg)

  test "addgg":
    let pq = norm2(gp+gq)
    norm2(grad(pq, gp) - 2.0*(gp+gq)) :< 1e-26
    norm2(grad(pq, gq) - 2.0*(gp+gq)) :< 1e-26
    ckgradm2(`+`, gq, gp, gu, gg, gm)

  test "mulsg":
    let p2 = norm2(x*gp)
    grad(p2, x) :~ 2.0*a*gp.norm2
    norm2(grad(p2, gp) - 2.0*a*a*gp) :< 1e-26
    ckgradm2(`*`, x, gp, y, gq, gg)

  test "mulgg":
    let pq = norm2(gp*gq)
    norm2(grad(pq, gp) - 2.0*gp*gq*gq.adj) :< 1e-24
    norm2(grad(pq, gq) - 2.0*gp.adj*gp*gq) :< 1e-24
    ckgradm2(`*`, gq, gp, gu, gg, gm)

  test "subgs":
    let p2 = norm2(gp-x)
    grad(p2, x) :~ retr(-2.0*(gp-a))
    norm2(grad(p2, gp) - 2.0*(gp-x)) :< 1e-26
    ckgradm2(`-`, gp, x, gq, y, gg)

  test "subgg":
    let pq = norm2(gp-gq)
    norm2(grad(pq, gp) - 2.0*(gp-gq)) :< 1e-26
    norm2(grad(pq, gq) - 2.0*(gq-gp)) :< 1e-26
    ckgradm2(`-`, gq, gp, gu, gg, gm)

  test "exp":
    let egp = exp(gp)
    norm2(egp.adj*egp - 1.0) :< 1e-20
    norm2(egp*egp.adj - 1.0) :< 1e-20
    ckgradm(exp, gm, 0.1*gp, gg)

  test "projTAH":
    let gt = gg.projTAH
    let tgt = gt.retr
    tgt*tgt :< 1e-26
    ckgradm(projTAH, gg, gp, gu)

suite "gauge fused":
  setup:
    let gg {.used.} = toGvalue g
    let gu {.used.} = toGvalue u
    let gp {.used.} = toGvalue p
    let gq {.used.} = toGvalue q
    let gm {.used.} = toGvalue m
    let x {.used.} = toGvalue a
    let y {.used.} = toGvalue b

  test "adjmul":
    let rf = gg.adjmul gu
    let rg = gg.adj * gu
    norm2(rf - rg) :< 1e-26
    let srf = rf.norm2
    let srg = rg.norm2
    norm2(grad(srf, gg) - grad(srg, gg)) :< 1e-25
    norm2(grad(srf, gu) - grad(srg, gu)) :< 1e-25
    ckgradm2(adjmul, gg, gu, gp, gq, gm)

  test "muladj":
    let rf = gg.muladj gu
    let rg = gg * gu.adj
    norm2(rf - rg) :< 1e-26
    let srf = rf.norm2
    let srg = rg.norm2
    norm2(grad(srf, gg) - grad(srg, gg)) :< 1e-25
    norm2(grad(srf, gu) - grad(srg, gu)) :< 1e-25
    ckgradm2(muladj, gg, gu, gp, gq, gm)

  test "contractProjTAH":
    let rf = contractProjTAH(gg, gu)
    let rg = projTAH(gg * gu.adj)
    norm2(rf - rg) :< 1e-26
    let srf = rf.norm2
    let srg = rg.norm2
    norm2(grad(srf, gg) - grad(srg, gg)) :< 1e-26
    norm2(grad(srf, gu) - grad(srg, gu)) :< 1e-26
    ckgradm2(contractProjTAH, gg, gu, gp, gq, gm)

  test "axexp":
    let rf = axexp(x, gm)
    let rg = exp(x*gm)
    norm2(rf - rg) :< 1e-26
    let srf = retr(rf * gu)
    let srg = retr(rg * gu)
    grad(srf, x) :~ grad(srg, x)
    norm2(grad(srf, gm) - grad(srg, gm)) :< 1e-26
    ckgradm2(axexp, x, gm, y, 0.05*gq, gp)

  test "axexpmuly":
    let rf = axexpmuly(x, gm, gg)
    let rg = exp(x*gm)*gg
    norm2(rf - rg) :< 1e-26
    let srf = retr(rf * gu)
    let srg = retr(rg * gu)
    grad(srf, x) :~ grad(srg, x)
    norm2(grad(srf, gm) - grad(srg, gm)) :< 1e-26
    norm2(grad(srf, gg) - grad(srg, gg)) :< 1e-26
    ckgradm3(axexpmuly, x, gm, gu, y, 0.05*gq, gg, gp)

suite "gauge action":
  let gplaq = block:
    var pl = 0.0
    for t in g.plaq:
      pl += t
    pl

  setup:
    let gg {.used.} = toGvalue g
    let gu {.used.} = toGvalue u
    let gm {.used.} = toGvalue m

  test "wilson action":
    let beta = 5.4
    let c = actWilson(beta)
    let s = gaugeAction(c, gg)
    s :~ -gplaq*float(6*vol*beta)
    proc act(x: Gvalue): Gvalue = gaugeAction(c, x)
    ckgrad(act, gg, gu)

  test "wilson force":
    let beta = 5.4
    let c = actWilson(beta)
    proc act(x: Gvalue): Gvalue = gaugeAction(c, x)
    proc force(x: Gvalue): Gvalue = gaugeForce(c, x)
    ckforce(act, force, gg, 10.0*gm)

  test "wilson force gradient":
    let beta = 5.4
    let c = actWilson(beta)
    proc force(x: Gvalue): Gvalue = gaugeForce(c, x)
    ckgradm(force, gg, gu, gm)

  test "wilson force gradient recomp":
    let beta = 5.4
    let c = actWilson(beta)
    let a = gaugeAction(c, gg)
    let f2 = gaugeForce(c, gg).norm2
    let df2 = grad(f2, gg).norm2
    let rs1 = [a.eval.getfloat, f2.eval.getfloat, df2.eval.getfloat]
    c.updated
    gg.updated
    let rs2 = [a.eval.getfloat, f2.eval.getfloat, df2.eval.getfloat]
    c.updated
    gg.updated
    let rs3 = [a.eval.getfloat, f2.eval.getfloat, df2.eval.getfloat]
    check rs1 == rs2
    check rs1 == rs3

qexFinalize()
