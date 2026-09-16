# Equivalences in definition order: gauge/transport, gauge/action/ops, gauge/stout.
# tgtower/tgtoweru1 cover the Wilson baseline; tgloops adds all fundamental
# loop families, coefficient transitions, and their derivative contractions.

when not declared(loopTestsOnly):
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

    test "gather and scatter = shifts, mutually adjoint":
      let
        f = linkField(gu, 0)
        b = linkField(gm, 1)
      for sh in [@[1], @[0, -2], @[1, -1], @[-1, 2]]:
        var r = f
        for mu, n in sh:
          if n != 0:
            r = shift(r, mu, -n)
        norm2(gather(f, sh) - r) :< 1e-18
        norm2(scatter(gather(f, sh), sh) - f) :< 1e-18
        (redot(gather(f, sh), b) - redot(f, scatter(b, sh))) :< 1e-8
      proc ga(x: Ggauge): Gscalar = redot(gather(linkField(x, 0) * linkField(x, 1), @[1, -1]), b)
      ckgrad(ga, gg, gu)
      proc sc(x: Ggauge): Gscalar = redot(scatter(linkField(x, 1), @[0, 2]), b * linkField(x, 0))
      ckgrad(sc, gg, gu)

    test "gp = flagged product of a and gathered b":
      let
        a = linkField(gg, 0)
        c = linkField(gu, 1)
        b = linkField(gm, 0)
        sh = @[1, -1]
        cs = gather(c, sh)
      for fa in [false, true]:
        for fb in [false, true]:
          let aa = (if fa: a.adj else: a)
          let bb = (if fb: cs.adj else: cs)
          norm2(gp(a, c, sh, fa, fb) - aa * bb) :< 1e-18
      proc g1(x: Ggauge): Gscalar = redot(gp(linkField(x, 0), linkField(x, 1), sh, true, false), b)
      ckgrad(g1, gg, gu)
      proc g2(x: Ggauge): Gscalar = redot(gp(linkField(x, 1), linkField(x, 0), sh, false, true), b)
      ckgrad(g2, gg, gu)
      proc g3(x: Ggauge): Gscalar = redot(gm, grad(g1(x), x))
      ckgrad(g3, gg, gu)

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
      # with origin the outputs start at the base site like hop chains do
      let og = lineProducts(gg, [@[-1, 2, 1, -2]], origin = true)
      norm2(og[0] - wilsonLine(gg, [-1, 2, 1, -2])) :< 1e-18

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

proc pgmSum(x: Ggauge): Gscalar =
  # Four pairs of three-link transports, matching gaugeAction2 independently
  # of the closed-path planner used by gaugeActionGraph.
  result = grt.toGvalue(0.0)
  for mu in 2..<x.gval.len:
    for nu in 1..<mu:
      for sg in 0..<nu:
        let a = linkField(x, mu)
        let b = linkField(x, nu)
        let d = linkField(x, sg)
        result = result + redot(hop(x, hop(x, d, nu, 1), mu, 1), hop(x, hop(x, a, nu, 1), sg, 1))
        result = result + redot(hop(x, hop(x, b, sg, 1), mu, 1), hop(x, hop(x, a, sg, 1), nu, 1))
        result = result + redot(hop(x, hop(x, d, mu, 1), nu, 1), hop(x, hop(x, b, mu, 1), sg, 1))
        result = result + redot(hop(x, hop(x, d, nu, -1), mu, 1), hop(x, hop(x, a, nu, -1), sg, 1))

suite "gauge action vs paths":
  teardown:
    grt.resetGradCache

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

  when declared(loopTestsOnly):
    test "parallelogram reference matches four open path contractions":
      let cc = grt.toGvalue(GaugeActionCoeffs(pgm: 1.0))
      let sf = gaugeActionGraph(cc, gg)
      let sr = (-1.0 / float(nc)) * pgmSum(gg)
      (sf - sr) :< 1e-8
      norm2(grad(sf, gg) - grad(sr, gg)) :< 1e-16

    test "fundamental production action, derivative, force, and Hessian match paths":
      let cc = grt.toGvalue(GaugeActionCoeffs())
      let sf = gaugeAction(cc, gg)
      let sr = gaugeActionGraph(cc, gg)
      let dr = grad(sr, gg)
      let hr = grad(redot(gm, dr), gg)
      let sub = gaugeActionDeriv(cc, gg, 1, subDir)
      let sh = grad(redot(gm, sub) + redot(gq, sub), gg)
      let sdif = norm2(sub - maskSubset(1, subDir, dr))
      let hdif = norm2(sh - grad(redot(maskSubset(1, subDir, gm + gq), dr), gg))
      let bdif = norm2(grad(redot(gq, sh), gm) - maskSubset(1, subDir, gaugeActionDeriv2(gq, cc, gg)))
      for gc in [GaugeActionCoeffs(rect: 0.17), GaugeActionCoeffs(pgm: -0.11),
                 GaugeActionCoeffs(plaq: 0.7, rect: -0.03, pgm: 0.02),
                 GaugeActionCoeffs(plaq: 1.0), GaugeActionCoeffs()]:
        cc.update gc
        (sf - sr) :< 1e-8
        sf :~ gc.gaugeAction2(g)
        norm2(gaugeActionDeriv(cc, gg) - dr) :< 1e-16
        norm2(gaugeForce(cc, gg) - projTAH(dr * gg.adj)) :< 1e-16
        norm2(gaugeActionDeriv2(gm, cc, gg) - hr) :< 1e-16
        sdif :< 1e-16
        hdif :< 1e-16
        bdif :< 1e-16
        (redot(gq, gaugeActionDeriv2(gm, cc, gg)) -
         redot(gm, gaugeActionDeriv2(gq, cc, gg))) :< 1e-8
      cc.update GaugeActionCoeffs(plaq: 0.7, rect: -0.03, pgm: 0.02)
      proc act(x: Ggauge): Gscalar = gaugeAction(cc, x)
      proc deriv(x: Ggauge): Ggauge = gaugeActionDeriv(cc, x)
      proc force(x: Ggauge): Ggauge = gaugeForce(cc, x)
      ckgrad(act, gg, gu)
      ckgradm(deriv, gg, gu, gm)
      ckgradm(force, gg, gu, gm)
      ckforce(act, force, gg, gm)

    test "identity coefficient bases stay nonzero when coefficients vanish":
      let id = gg.unitGaugeLike
      let cc = grt.toGvalue(GaugeActionCoeffs())
      let dc = grad(gaugeAction(cc, id), cc)
      let rc = grad(gaugeActionGraph(cc, id), cc)
      let nd = g.len
      let nv = float(g[0].l.physVol)
      let want = GaugeActionCoeffs(plaq: -0.5*float(nd*(nd-1))*nv,
                                   rect: -float(nd*(nd-1))*nv,
                                   pgm: -(2.0/3.0)*float(nd*(nd-1)*(nd-2))*nv)
      for gc in [GaugeActionCoeffs(), GaugeActionCoeffs(plaq: 1.0, rect: 0.1, pgm: -0.2),
                 GaugeActionCoeffs(plaq: 1.0), GaugeActionCoeffs()]:
        cc.update gc
        discard dc.eval
        discard rc.eval
        for got, expected in fields(dc.cval, want):
          check almostEqual(got, expected)
        for got, expected in fields(rc.cval, want):
          check almostEqual(got, expected)
      cc.update GaugeActionCoeffs(rect: 0.1, adjplaq: 0.2)
      expect GraphValueError:
        discard dc.eval
      expect GraphValueError:
        discard rc.eval

    test "zero coefficient gradients and mixed field derivatives remain live":
      let cc = grt.toGvalue(GaugeActionCoeffs())
      let dc = grad(gaugeAction(cc, gg), cc)
      let rc = grad(gaugeActionGraph(cc, gg), cc)
      let p = 1
      let dir = subDir
      let ds = gaugeActionDeriv(cc, gg, p, dir)
      let hc = grad(redot(gm + gq, ds), gg)
      let dcoef = grad(redot(gm, gaugeActionDeriv(cc, gg)), cc)
      let fcoef = grad(redot(gm, gaugeForce(cc, gg)), cc)
      let scoef = grad(redot(gm, ds), cc)
      let hcoef = grad(redot(gq, gaugeActionDeriv2(gm, cc, gg)), cc)
      let qcoef = grad(redot(gq, hc), cc)
      var checks, norms: seq[Gscalar]
      for basis in [GaugeActionCoeffs(plaq: 1.0), GaugeActionCoeffs(rect: 1.0), GaugeActionCoeffs(pgm: 1.0)]:
        let cb = grt.toGvalue(basis)
        let sb = gaugeActionGraph(cb, gg)
        let db = grad(sb, gg)
        let hb = grad(redot(gm, db), gg)
        checks.add redot(dc, cb) - sb
        checks.add redot(rc, cb) - sb
        norms.add norm2(grad(redot(dc, cb), gg) - db)
        checks.add redot(dcoef, cb) - redot(gm, db)
        checks.add redot(fcoef, cb) - redot(gm, projTAH(db * gg.adj))
        checks.add redot(scoef, cb) - redot(maskSubset(p, dir, gm), db)
        checks.add redot(hcoef, cb) - redot(gq, hb)
        checks.add redot(qcoef, cb) - redot(gq, grad(redot(maskSubset(p, dir, gm + gq), db), gg))
      for gc in [GaugeActionCoeffs(), GaugeActionCoeffs(plaq: 0.7, rect: -0.03, pgm: 0.02),
                 GaugeActionCoeffs(plaq: 1.0), GaugeActionCoeffs()]:
        cc.update gc
        for k, ch in checks:
          checkpoint("coefficient contraction " & $k & ", coefficients " & $gc)
          ch :< 1e-8
        for k, ch in norms:
          checkpoint("coefficient field derivative " & $k & ", coefficients " & $gc)
          ch :< 1e-16
      let cb = grt.toGvalue(GaugeActionCoeffs(rect: 0.3, pgm: -0.2))
      proc mixed(x: Ggauge): Gscalar = redot(grad(gaugeAction(cc, x), cc), cb)
      checkpoint("mixed coefficient/field finite difference")
      ckgrad(mixed, gg, gm)
      let t = grt.toGvalue(0.0)
      let score = gaugeAction(t * cb, gg)
      let dd = grad(grad(score * score, t), t)
      let basisScore = gaugeActionGraph(cb, gg)
      checkpoint("second coefficient derivative with live scalar upstream at zero")
      dd :~ 2.0 * basisScore * basisScore
      checkpoint("cached coefficient derivative after field update")
      gg.update u
      redot(dc, cb) :~ gaugeActionGraph(cb, gg)
      cc.update GaugeActionCoeffs(adjplaq: 0.2)
      expect GraphValueError:
        discard hcoef.eval
      expect GraphValueError:
        discard qcoef.eval

    test "cached coefficient pullbacks follow fundamental and adjoint families":
      let cc = grt.toGvalue(GaugeActionCoeffs(plaq: 1.0))
      let cb = grt.toGvalue(GaugeActionCoeffs(plaq: 0.8, rect: -0.03, pgm: 0.02, adjplaq: 0.4))
      let dc = [grad(gaugeAction(cc, gg), cc),
                grad(redot(gm, gaugeActionDeriv(cc, gg)), cc),
                grad(redot(gm, gaugeForce(cc, gg)), cc)]
      let scores = [redot(dc[0], cb), redot(dc[1], cb), redot(dc[2], cb)]
      let mixed = [grad(scores[0], gg), grad(scores[1], gg), grad(scores[2], gg)]
      let cf = grt.toGvalue(GaugeActionCoeffs(plaq: 0.8, rect: -0.03, pgm: 0.02))
      let ca = grt.toGvalue(GaugeActionCoeffs(plaq: 0.8, adjplaq: 0.4))
      var errs, norms: array[2, array[3, Gscalar]]
      for family, refa in [gaugeActionGraph(cf, gg), adjPlaqAction(ca, gg)]:
        let d = grad(refa, gg)
        let refs = [refa, redot(gm, d), redot(gm, projTAH(d * gg.adj))]
        for k in 0..2:
          errs[family][k] = scores[k] - refs[k]
          norms[family][k] = norm2(mixed[k] - grad(refs[k], gg))
      for step, gc in [GaugeActionCoeffs(plaq: 0.7, rect: -0.03, pgm: 0.02),
                      GaugeActionCoeffs(plaq: 1.2, adjplaq: 0.3),
                      GaugeActionCoeffs(),
                      GaugeActionCoeffs(plaq: 0.5, adjplaq: -0.2),
                      GaugeActionCoeffs(plaq: 1.0)]:
        cc.update gc
        if step == 2:
          gg.update u
        let family = if gc.adjplaq == 0.0: 0 else: 1
        for k in 0..2:
          checkpoint("cached coefficient pullback " & $k & ", coefficients " & $gc)
          errs[family][k] :< 1e-8
          norms[family][k] :< 1e-16
          if family == 0:
            check dc[k].cval.adjplaq == 0.0
          else:
            check dc[k].cval.rect == 0.0
            check dc[k].cval.pgm == 0.0

    test "actAdj at zero beta uses the current fundamental coefficient derivative":
      let bt = grt.toGvalue(0.0)
      let af = grt.toGvalue(0.25)
      let sa = gaugeAction(actAdj(bt, af), gg)
      let db = grad(sa, bt)
      let da = grad(sa, af)
      let sf = (-1.0 / float(nc)) * plaqSum(gg)
      let sr = adjPlaqAction(actAdj(grt.toGvalue(1.0), af), gg)
      let df = norm2(grad(db, gg) - grad(sf, gg))
      let dr = norm2(grad(db, gg) - grad(sr, gg))
      for b in [0.0, 0.8, 0.0]:
        bt.update b
        if b == 0.0:
          sa :~ 0.0
          da :~ 0.0
          (db - sf) :< 1e-8
          df :< 1e-16
        else:
          (db - sr) :< 1e-8
          dr :< 1e-16

  test "adjPlaqAction = gaugeAction for the adjoint-plaquette family":
    let
      cc = grt.toGvalue(GaugeActionCoeffs(plaq: 1.2, adjplaq: 0.3))
      sf = adjPlaqAction(cc, gg)
      sr = gaugeAction(cc, gg)
    (sf - sr) :< 1e-8
    norm2(grad(sf, gg) - gaugeActionDeriv(cc, gg)) :< 1e-16
    # linear in both coefficients
    let adjFac = grt.toGvalue(0.25)
    let sa = adjPlaqAction(actAdj(beta, adjFac), gg)
    (beta * grad(sa, beta) - sa) :< 1e-6
    var adjBasis = grt.toGvalue(0.0)
    for mu in 1..<g.len:
      for nu in 0..<mu:
        let trp = trace(wilsonLine(gg, plaqPath(mu, nu)))
        adjBasis = adjBasis + float(g[0].l.physVol) - (1.0 / float(nc*nc)) * norm2(trp)
    (grad(sa, adjFac) - beta * adjBasis) :< 1e-6
    (grad(gaugeAction(actAdj(beta, adjFac), gg), adjFac) - beta * adjBasis) :< 1e-6

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
  teardown:
    grt.resetGradCache

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

  when declared(loopTestsOnly):
    test "ordinary loop subsets retain the Wilson stout factorization guard":
      let cc = grt.toGvalue(GaugeActionCoeffs(plaq: 1.0))
      let st = stoutUpdateLogDetJ(gg, cc, alpha, 1, subDir)
      discard st.Wnew.eval
      discard st.lj.eval
      for gc in [GaugeActionCoeffs(rect: 0.1), GaugeActionCoeffs(pgm: 0.1)]:
        cc.update gc
        discard gaugeActionDeriv(cc, gg, 1, subDir).eval
        expect GraphValueError:
          discard st.Wnew.eval
        expect GraphValueError:
          discard st.lj.eval

    test "grouped Wilson stout coefficient gradient matches scalar differences":
      let b = grt.toGvalue(1.0)
      let cc = actWilson(b)
      let st = stoutUpdateLogDetJ(gg, cc, alpha, 1, subDir)
      let score = redot(st.Wnew, gm) + st.lj
      let dc = grad(score, cc)
      let db = grad(score, b)
      let (dv, err) = ndiff(score, b)
      b.update 0.0
      check(instantiationInfo(), "grouped stout coefficient gradient", dv, err, db.eval.sval)
      discard dc.eval
      check dc.cval.rect == 0.0
      check dc.cval.pgm == 0.0
      check dc.cval.adjplaq == 0.0
      let raw = grt.toGvalue(GaugeActionCoeffs(plaq: 1.0))
      let rr = stoutUpdateLogDetJ(gg, raw, alpha, 1, subDir)
      let dr = grad(redot(rr.Wnew, gm) + rr.lj, raw)
      discard dr.eval
      raw.update GaugeActionCoeffs(plaq: 1.0, rect: 0.1)
      expect GraphValueError:
        discard dr.eval

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
