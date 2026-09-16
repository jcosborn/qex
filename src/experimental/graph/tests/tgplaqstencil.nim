#RUNCMD env OMP_NUM_THREADS=1 $RUNJOB

import base/globals
setVLENmax(4)

import math, unittest
import qex except epsilon
import algorithms/numdiff
import helpers
import ../[core, scalar, gauge, functional]
import ../gauge/[types, basic_ops, field_ops, transport, stencil]

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))
let grt = initGraphRuntime()

proc plaqRef(x: Ggauge): Gscalar =
  result = grt.toGvalue(0.0)
  for mu in 1..<x.gval.len:
    for nu in 0..<mu:
      result = result + retr(wilsonLine(x, plaqPath(mu, nu)))

proc jetRef(x: Ggauge, ds: openArray[Ggauge]): Ggauge =
  let slot = slotVar(x)
  var d = gradSeeded(plaqRef(slot), slot, grt.toGvalue(1.0))
  for seed in ds:
    d = gradSeeded(d, slot, seed)
  Ggauge(d)

proc compare(a, b: Ggauge, tol = 2.0e-22) =
  let e = norm2(a-b).eval.sval
  let s = norm2(b).eval.sval
  check e/(1.0+s) < tol

proc compare(a, b: Gscalar, tol = 2.0e-11) =
  let av = a.eval.sval
  let bv = b.eval.sval
  check abs(av-bv) < tol*(1.0+abs(bv))

qexInit()
letParam:
  expectRanks = nRanks
check nRanks == expectRanks
letParam:
  lat = latticeFromLocalLattice(@[2,4,4], nRanks)
let
  lo = lat.newLayout
  g = lo.newGauge
  h = lo.newGauge
  k = lo.newGauge
  q = lo.newGauge
  v = lo.newGauge
var r = lo.newRNGField(Philox4x64, 920260911'u64)
threads:
  g.random r
  h.gaussian r
  k.gaussian r
  q.gaussian r
  v.gaussian r
  for mu in 0..<g.len:
    h[mu] *= 0.1
    k[mu] *= 0.1
    q[mu] *= 0.1
    v[mu] *= 0.1
    g[mu] += h[mu]
let
  x = grt.toGvalue(g)
  a = grt.toGvalue(h)
  b = grt.toGvalue(k)
  c = grt.toGvalue(q)
  u = grt.toGvalue(v)
  all = @[a,b,c]

suite "Fused plaquette stencil":
  test "plaquette and staple values match independent hop chains":
    compare(plaqSum(x), plaqRef(x))
    compare(stapleSum(x), jetRef(x, @[]))
    compare(Ggauge(grad(plaqSum(x), x)), stapleSum(x))

  test "every jet's field and seed pullbacks match hop derivatives":
    for order in 0..3:
      let
        ds = all[0..<order]
        fused = stapleSum(x, ds)
        reference = jetRef(x, ds)
        f = redot(u, fused)
        rr = redot(u, reference)
      compare(fused, reference)
      compare(Ggauge(grad(f, x)), Ggauge(grad(rr, x)))
      for seed in ds:
        compare(Ggauge(grad(f, seed)), Ggauge(grad(rr, seed)))

  test "jets above cubic order and their pullbacks are exactly zero":
    let
      z = stapleSum(x, @[a,b,c,a])
      f = redot(u, z)
    check norm2(z).eval.sval == 0.0
    for arg in [x,a,b,c,u]:
      check norm2(grad(f, arg)).eval.sval == 0.0

  test "aliased and field-dependent seeds remain live":
    for ds in [@[x], @[a,a], @[x*x,a,x.adj], @[a,a,a]]:
      let
        fused = redot(x*u, stapleSum(x, ds))
        reference = redot(x*u, jetRef(x, ds))
      compare(fused, reference)
      let
        df = Ggauge(grad(fused, x))
        dr = Ggauge(grad(reference, x))
      compare(df, dr)
      compare(Ggauge(grad(redot(a, df), x)), Ggauge(grad(redot(a, dr), x)))
      compare(Ggauge(grad(fused, a)), Ggauge(grad(reference, a)))

  test "coefficient bases survive zero crossings and cached updates":
    let
      beta = grt.toGvalue(0.0)
      fp = beta*plaqSum(x)
      rp = beta*plaqRef(x)
      fj = redot(u, beta*stapleSum(x, @[a,b]))
      rj = redot(u, beta*jetRef(x, @[a,b]))
      dbp = Gscalar(grad(fp, beta))
      dbr = Gscalar(grad(rp, beta))
      djp = Gscalar(grad(fj, beta))
      djr = Gscalar(grad(rj, beta))
      mixp = Ggauge(grad(djp, x))
      mixr = Ggauge(grad(djr, x))
    for value in [0.0, 0.7, 0.0, -0.3]:
      beta.update(value)
      compare(fp, rp)
      compare(fj, rj)
      compare(dbp, dbr)
      compare(djp, djr)
      compare(mixp, mixr)
    x.update(g)
    a.update(k)
    compare(fj, rj)
    compare(mixp, mixr)
    a.update(h)

  test "cloned lambda applications own and rebind stencil work":
    let
      y = grt.toGvalue(k)
      body = plaqSum(x)+0.2*norm2(stapleSum(x, @[a]))
      fun = lambda(x, body)
    discard body.eval  # clone an already allocated workspace
    let
      applied = Gscalar(apply(fun, y))
      direct = plaqSum(y)+0.2*norm2(stapleSum(y, @[a]))
    compare(applied, direct)
    y.update(q)
    a.update(v)
    compare(applied, direct)
    compare(Ggauge(grad(applied, y)), Ggauge(grad(direct, y)))
    compare(body, plaqSum(x)+0.2*norm2(stapleSum(x, @[a])))
    a.update(h)

  test "third-jet seed pullbacks agree with ambient finite differences":
    let
      t = grt.toGvalue(0.0)
      value = redot(u, stapleSum(x, @[a+t*b, c, x]))
      want = redot(b, Ggauge(grad(redot(u, stapleSum(x, @[a,c,x])), a)))
    proc at(z: float): float =
      t.update(z)
      value.eval.sval
    var num, err: float
    ndiff(num, err, at, 0.0, 0.125, ordMax=3)
    let ana = want.eval.sval
    check abs(num-ana) < max(1.0e-10, 32.0*err)

qexFinalize()
