# Equivalences in definition order: gauge/transport, gauge/action/ops, gauge/stout.
# Run on the small SU(3) and U(1) fixtures in tgtower and tgtoweru1.
# Action derivatives use Wilson coefficients; rectangles cover S and dS only.

suite "gauge transport vs basic":
  setup:
    let gg {.used.} = grt.toGvalue(g)
    let gu {.used.} = grt.toGvalue(u)
    let gm {.used.} = grt.toGvalue(m)

  test "hop = link multiplication and shift, for both signs":
    let f = linkField(gu, 0)
    let b = linkField(gm, 1)
    for mu in 0..<g.len:
      for (rf, rr) in [
          (hop(gg, f, mu, 1), linkField(gg, mu) * shift(f, mu, 1)),
          (hop(gg, f, mu, -1), shift(linkField(gg, mu).adj * f, mu, -1))]:
        let sf = redot(rf, b)
        let sr = redot(rr, b)
        norm2(rf - rr) :< 1e-18
        norm2(grad(sf, gg) - grad(sr, gg)) :< 1e-20
        norm2(grad(sf, f) - grad(sr, f)) :< 1e-20

  test "transport([1, 2, -1]) = hop_1(hop_2(hop_-1(f)))":
    let
      f = linkField(gu, 0)
      b = linkField(gm, 1)
      rf = transport(gg, f, [1, 2, -1])
      rr = hop(gg, hop(gg, hop(gg, f, 0, -1), 1, 1), 0, 1)
      sf = redot(rf, b)
      sr = redot(rr, b)
      df = grad(sf, gg)
      dr = grad(sr, gg)
    for k in 0..1:
      norm2(rf - rr) :< 1e-18
      norm2(df - dr) :< 1e-20
      norm2(grad(sf, f) - grad(sr, f)) :< 1e-20
      if k == 0:
        gg.update u
    norm2(transport(gg, f, []) - f) :< 1e-24

  test "lineProducts = hop chains, shared plaquette pair, rectangles":
    let
      ps = lineProducts(gg, [@(plaqPath(0, 1)), @[1, 1, 2, -1, -1, -2], @[2, 1, -2, -1]])
      f0 = linkField(gm, 0)
      f1 = linkField(gm, 1)
      rf = redot(ps[0], f0) + redot(ps[1], f1) + redot(ps[2], f0)
      rr = redot(wilsonLine(gg, plaqPath(0, 1)), f0) +
           redot(wilsonLine(gg, [1, 1, 2, -1, -1, -2]), f1) +
           redot(wilsonLine(gg, [2, 1, -2, -1]), f0)
    for i in 0..2:
      norm2(ps[i] - wilsonLine(gg, (if i == 0: @[1, 2, -1, -2] elif i == 1: @[1, 1, 2, -1, -1, -2] else: @[2, 1, -2, -1]))) :< 1e-18
    (rf - rr) :< 1e-8
    norm2(grad(rf, gg) - grad(rr, gg)) :< 1e-20
    # the two orientations of one plaquette share their product node
    let shared = lineProducts(gg, [@[1, 2, -1, -2], @[2, 1, -2, -1]])
    check shared[1].inputs[0].nodeKey == shared[0].nodeKey

  test "wilsonLine plaquette = U_mu shift_mu(U_nu) shift_nu(U_mu).adj U_nu.adj":
    let
      a = linkField(gg, 0)
      b = linkField(gg, 1)
      rf = wilsonLine(gg, plaqPath(0, 1))
      rr = a * shift(b, 0, 1) * shift(a, 1, 1).adj * b.adj
      sf = redot(rf, linkField(gm, 0))
      sr = redot(rr, linkField(gm, 0))
    norm2(rf - rr) :< 1e-18
    norm2(grad(sf, gg) - grad(sr, gg)) :< 1e-20
    norm2(wilsonLine(gg, []) - gg.unitFieldLike) :< 1e-24

proc plaqSum(x: Ggauge): Gscalar =
  # sum_{mu>nu} retr[U_mu U_nu(+mu) U_mu(+nu)^dag U_nu^dag].
  result = grt.toGvalue(0.0)
  for mu in 1..<x.gval.len:
    for nu in 0..<mu:
      let a = linkField(x, mu)
      let b = linkField(x, nu)
      result = result + retr(a * shift(b, mu, 1) * shift(a, nu, 1).adj * b.adj)

proc staples(x: Ggauge): Ggauge =
  # C_mu = sum_{nu!=mu} [U_nu U_mu(+nu) U_nu(+mu)^dag
  #                     + U_nu(-nu)^dag U_mu(-nu) U_nu(+mu-nu)].
  result = Ggauge(x.zeroLike)
  for mu in 0..<x.gval.len:
    let a = linkField(x, mu)
    var s = Gfield(a.zeroLike)
    for nu in 0..<x.gval.len:
      if nu == mu:
        continue
      let b = linkField(x, nu)
      s = s + hop(x, a, nu, 1) * shift(b, mu, 1).adj + hop(x, hop(x, b, mu, 1), nu, -1)
    result = result + injectLink(s, mu, x)

suite "gauge action vs paths":
  setup:
    let gg {.used.} = grt.toGvalue(g)
    let gu {.used.} = grt.toGvalue(u)
    let gm {.used.} = grt.toGvalue(m)
    let gq {.used.} = grt.toGvalue(q)
    let beta {.used.} = grt.toGvalue(5.4)
    let c {.used.} = actWilson(beta)
    const nc = g[0][0].nrows

  test "gaugeActionGraph = -beta/Nc times the plaquette trace sum":
    let
      sf = gaugeActionGraph(c, gg)
      sr = (-1.0 / float(nc)) * beta * plaqSum(gg)
      df = grad(sf, gg)
      dr = grad(sr, gg)
      bf = grad(sf, beta)
      br = grad(sr, beta)
    for k in 0..1:
      (sf - sr) :< 1e-8
      norm2(df - dr) :< 1e-16
      bf :~ br
      if k == 0:
        beta.update 6.1
        gg.update u

  test "gaugeAction Wilson = -beta/Nc times the plaquette trace sum":
    let
      sf = gaugeAction(c, gg)
      sr = (-1.0 / float(nc)) * beta * plaqSum(gg)
    (sf - sr) :< 1e-8
    norm2(grad(sf, gg) - grad(sr, gg)) :< 1e-16

  test "gaugeActionGraph with rectangles = gaugeAction, in value and derivative":
    let
      cc = grt.toGvalue(GaugeActionCoeffs(plaq: 1.3, rect: -0.1))
      sf = gaugeActionGraph(cc, gg)
      sr = gaugeAction(cc, gg)
    (sf - sr) :< 1e-8
    norm2(grad(sf, gg) - gaugeActionDeriv(cc, gg)) :< 1e-16
    # linear in beta: beta dS/dbeta = S for a rectangle family
    let sym = gaugeActionGraph(actSymanzik(beta), gg)
    (beta * grad(sym, beta) - sym) :< 1e-6

  test "gaugeAction with rectangles = weighted plaquette and 1x2 Wilson loops":
    let cc = grt.toGvalue(GaugeActionCoeffs(plaq: 1.3, rect: -0.1))
    var r = grt.toGvalue(0.0)
    for mu in 1..<g.len:
      for nu in 0..<mu:
        let a = mu + 1
        let b = nu + 1
        r = r + retr(wilsonLine(gg, [a, a, b, -a, -a, -b]))
        r = r + retr(wilsonLine(gg, [a, b, b, -a, -b, -b]))
    let
      sf = gaugeAction(cc, gg)
      sr = (-1.0 / float(nc)) * (1.3 * plaqSum(gg) - 0.1 * r)
    (sf - sr) :< 1e-8
    norm2(grad(sf, gg) - grad(sr, gg)) :< 1e-16

  test "gaugeActionDeriv Wilson = -beta/Nc times the staple sum":
    let
      rf = gaugeActionDeriv(c, gg)
      rr = (-1.0 / float(nc)) * beta * staples(gg)
      s = (-1.0 / float(nc)) * beta * plaqSum(gg)
      sf = redot(rf, gm)
      sr = redot(rr, gm)
    norm2(rf - rr) :< 1e-16
    norm2(rr - grad(s, gg)) :< 1e-16
    norm2(grad(sf, gg) - grad(sr, gg)) :< 1e-16

  test "gaugeForce Wilson = projTAH((-beta/Nc * staples) * adj(g))":
    let
      d = (-1.0 / float(nc)) * beta * staples(gg)
      rf = gaugeForce(c, gg)
      rr = projTAH(d * gg.adj)
      sf = redot(rf, gm)
      sr = redot(rr, gm)
    norm2(rf - rr) :< 1e-16
    norm2(grad(sf, gg) - grad(sr, gg)) :< 1e-16

  test "gaugeActionDeriv subset = masked staple sum, with neighbour pullback":
    let
      d = (-1.0 / float(nc)) * beta * staples(gg)
      zero = Ggauge(gg.zeroLike)
    for parity in 0..1:
      for dir in 0..<g.len:
        let
          rf = gaugeActionDeriv(c, gg, parity, dir)
          rr = blendSubset(parity, dir, d, zero)
          sf = redot(rf, gm)
          sr = redot(rr, gm)
        norm2(rf - rr) :< 1e-16
        norm2(grad(sf, gg) - grad(sr, gg)) :< 1e-16

  test "gaugeActionDeriv2 Wilson = grad(redot(b, -beta/Nc * staples), g)":
    let
      d = (-1.0 / float(nc)) * beta * staples(gg)
      rf = gaugeActionDeriv2(gm, c, gg)
      rr = grad(redot(gm, d), gg)
      sf = redot(rf, gq)
      sr = redot(rr, gq)
    norm2(rf - rr) :< 1e-16
    norm2(grad(sf, gm) - grad(sr, gm)) :< 1e-16
    norm2(grad(sf, gg) - grad(sr, gg)) :< 1e-16

  test "subset Hessian with summed seeds = staple pullback of the masked sum":
    let
      d = (-1.0 / float(nc)) * beta * staples(gg)
      sub = gaugeActionDeriv(c, gg, 1, subDir)
      seed = blendSubset(1, subDir, gm + gq, Ggauge(gg.zeroLike))
      rf = grad(redot(gm, sub) + redot(gq, sub), gg)
      rr = grad(redot(seed, d), gg)
    norm2(rf - rr) :< 1e-16

suite "stout update vs paths":
  setup:
    let gg {.used.} = grt.toGvalue(g)
    let gm {.used.} = grt.toGvalue(m)
    let beta {.used.} = grt.toGvalue(5.4)
    let c {.used.} = actWilson(beta)
    let alpha {.used.} = grt.toGvalue(0.02)
    const nc = g[0][0].nrows
    let ds = blendSubset(1, subDir, (-1.0 / float(nc)) * beta * staples(gg), Ggauge(gg.zeroLike))

  test "stoutUpdate(g, ds, alpha) = subset exp(alpha * projTAH(g * adj(ds))) * g":
    let
      rf = stoutUpdate(gg, ds, alpha, 1, subDir)
      rr = blendSubset(1, subDir, exp(alpha * projTAH(gg * ds.adj)) * gg, gg)
      sf = redot(rf, gm)
      sr = redot(rr, gm)
    norm2(rf - rr) :< 1e-20
    norm2(grad(sf, gg) - grad(sr, gg)) :< 1e-16
    grad(sf, alpha) :~ grad(sr, alpha)

  test "stoutUpdateLogDetJ with supplied staples: Wnew = basic subset update":
    let
      rf = stoutUpdateLogDetJ(gg, ds, alpha, 1, subDir).Wnew
      rr = blendSubset(1, subDir, exp(alpha * projTAH(gg * ds.adj)) * gg, gg)
      sf = redot(rf, gm)
      sr = redot(rr, gm)
    norm2(rf - rr) :< 1e-20
    norm2(grad(sf, gg) - grad(sr, gg)) :< 1e-16
    grad(sf, alpha) :~ grad(sr, alpha)

  test "stoutUpdateLogDetJ with coefficients: Wnew = basic staple subset update":
    # lj additionally needs a site-local matrix Jacobian, beyond shift/hop algebra.
    let
      rf = stoutUpdateLogDetJ(gg, c, alpha, 1, subDir).Wnew
      rr = blendSubset(1, subDir, exp(alpha * projTAH(gg * ds.adj)) * gg, gg)
      sf = redot(rf, gm)
      sr = redot(rr, gm)
    norm2(rf - rr) :< 1e-20
    norm2(grad(sf, gg) - grad(sr, gg)) :< 1e-16
    grad(sf, alpha) :~ grad(sr, alpha)
