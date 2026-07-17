suite "gauge action":
  let gplaq = block:
    var pl = 0.0
    for t in g.plaq:
      pl += t
    pl

  setup:
    let gg {.used.} = grt.toGvalue(g)
    let gu {.used.} = grt.toGvalue(u)
    let gm {.used.} = grt.toGvalue(m)

  test "wilson action":
    let beta = 5.4
    let c = actWilson(scalar.toGvalue(grt, beta))
    let s = gaugeAction(c, gg)
    s :~ -gplaq*float(6*vol*beta)
    proc act(x: Ggauge): Gscalar = gaugeAction(c, x)
    ckgrad(act, gg, gu)

  test "wilson force":
    let beta = 5.4
    let c = actWilson(scalar.toGvalue(grt, beta))
    proc act(x: Ggauge): Gscalar = gaugeAction(c, x)
    proc force(x: Ggauge): Ggauge = gaugeForce(c, x)
    ckforce(act, force, gg, 10.0*gm)

  test "wilson force matches separate derivative and projection":
    let
      c = actWilson(scalar.toGvalue(grt, 5.4))
      fused = gaugeForce(c, gg)
      separate = contractProjTAH(gaugeActionDeriv(c, gg), gg)
    norm2(fused - separate) :< 1e-24

  test "wilson force gradient":
    let beta = 5.4
    let c = actWilson(scalar.toGvalue(grt, beta))
    proc force(x: Ggauge): Ggauge = gaugeForce(c, x)
    ckgradm(force, gg, gu, gm)

  test "gaugeActionDerivSubset matches masked full derivative and ndiff":
    let
      c = actWilson(scalar.toGvalue(grt, 5.4))
      zero = grt.toGvalue(zeroGaugeLike(g))
      full = gaugeActionDeriv(c, gg)
    for parity in 0..1:
      for dir in 0..<g.len:
        let
          sub = gaugeActionDeriv(c, gg, parity, dir)
          reference = blendSubset(parity, dir, full, zero)
        norm2(sub - reference) :< 1e-24

    let
      parity = 1
      dir = 2
      aSub = blendSubset(parity, dir, gu, zero)
      t = grt.toGvalue(0.0)
      action = gaugeAction(c, gg + t*aSub)
      (d, e) = ndiff(action, t)
      analytic = redot(aSub, gaugeActionDeriv(c, gg, parity, dir)).eval.sval
    check(instantiationInfo(), "gaugeActionDeriv subset action derivative", d, e, analytic)

  test "raw subset derivative overwrites only its requested output":
    let c = GaugeActionCoeffs(plaq: 5.4)
    var full = zeroGaugeLike(g)
    c.gaugeDeriv2(g, full)
    for p in 0..1:
      let sub = g[0].l.getSubset(if p == 0: "even" else: "odd")
      for d in 0..<g.len:
        var
          got = zeroGaugeLike(g)
          reference = zeroGaugeLike(g)
        threads:
          for mu in 0..<g.len:
            got[mu] := u[mu]
            reference[mu] := u[mu]
          threadBarrier()
          reference[d] := 0
          threadBarrier()
          for x in sub:
            reference[d][x] := full[d][x]
        c.gaugeDeriv2Subset(g, got, p, d)
        norm2(grt.toGvalue(got) - grt.toGvalue(reference)) :< 1e-24

  test "subset derivative refresh overwrites graph-owned storage":
    let
      c = actWilson(scalar.toGvalue(grt, 5.4))
      x = grt.toGvalue(g)
      zero = grt.toGvalue(zeroGaugeLike(g))
      p = 1
      d = min(2, g.len - 1)
      got = gaugeActionDeriv(c, x, p, d)
      full = gaugeActionDeriv(c, x)
      reference = blendSubset(p, d, full, zero)
    norm2(got - reference) :< 1e-24
    x.update(u)
    norm2(got - reference) :< 1e-24

  test "subset derivative functional clones own shift workspace":
    let
      c = actWilson(scalar.toGvalue(grt, 5.4))
      p = 1
      d = min(2, g.len - 1)
      xg = grt.toGvalue(g)
      xu = grt.toGvalue(u)
      x = Ggauge(xg.newOneOf)
      body = gaugeActionDeriv(c, x, p, d)
      fn = lambda(x, body)
      gotG = Ggauge(apply(fn, xg))
      gotU = Ggauge(apply(fn, xu))
      refG = gaugeActionDeriv(c, xg, p, d)
      refU = gaugeActionDeriv(c, xu, p, d)
    norm2(gotG - refG) :< 1e-24
    norm2(gotU - refU) :< 1e-24
    xg.update(u)
    xu.update(g)
    norm2(gotU - refU) :< 1e-24
    norm2(gotG - refG) :< 1e-24

  test "gaugeActionDeriv2Subset matches the full Hessian and ndiff":
    let
      c = actWilson(scalar.toGvalue(grt, 5.4))
      zero = grt.toGvalue(zeroGaugeLike(g))
    for parity in 0..1:
      for dir in 0..<g.len:
        let
          p = parity
          d = dir
        proc subsetDeriv(x: Ggauge): Ggauge =
          gaugeActionDeriv(c, x, p, d)
        ckgradm(subsetDeriv, gg, gu, gm)

        let
          bSub = blendSubset(p, d, gm, zero)
          got = grad(redot(subsetDeriv(gg), gm), gg)
          reference = gaugeActionDeriv2(bSub, c, gg)
        norm2(got - reference) :< 1e-22

    block:
      let
        p = 1
        d = min(2, g.len - 1)
        sub = gaugeActionDeriv(c, gg, p, d)
        split = redot(sub, gm) + redot(sub, gu)
        got = grad(split, gg)
        bSub = blendSubset(p, d, gm + gu, zero)
        reference = gaugeActionDeriv2(bSub, c, gg)
      norm2(got - reference) :< 1e-22

  test "wilson subset derivative rejects unsupported coefficients":
    let c = actSymanzik(scalar.toGvalue(grt, 5.4))
    expect(GraphValueError):
      discard gaugeActionDeriv(c, gg, 0, 0).eval

  test "raw subset kernels reject unsupported coefficients":
    let c = GaugeActionCoeffs(rect: 1.0)
    var outg = zeroGaugeLike(g)
    let w = g[0].newOneOf
    expect(ValueError):
      c.gaugeDeriv2Subset(g, outg, 0, 0)
    expect(ValueError):
      c.gaugeDerivDeriv2Subset(g, m, outg, 0, 0)
    expect(ValueError):
      c.gaugeDerivDeriv2SubsetSum(g, @[m[0]], w, outg, 0, 0)
    expect(ValueError):
      c.gaugeDerivDeriv2SubsetAdd(g, m[0], outg, 0, 0)
    expect(ValueError):
      c.gaugeDerivDeriv2SubsetAddBase(g, m[0], u, outg, 0, 0)

  test "raw summed subset Hessian matches a materialized sum":
    let c = GaugeActionCoeffs(plaq: 5.4)
    for p in 0..1:
      for d in 0..<g.len:
        let
          hs = @[m[d], u[d], m[d]]
          w = g[d].newOneOf
        var
          h = zeroGaugeLike(g)
          got = zeroGaugeLike(g)
          reference = zeroGaugeLike(g)
        threads:
          h[d] := m[d]
          h[d] += u[d]
          h[d] += m[d]
        c.gaugeDerivDeriv2SubsetSum(g, hs, w, got, p, d)
        c.gaugeDerivDeriv2Subset(g, h, reference, p, d)
        norm2(grt.toGvalue(got) - grt.toGvalue(reference)) :< 1e-24

  test "raw accumulating subset Hessian adds to an arbitrary output":
    let c = GaugeActionCoeffs(plaq: 5.4)
    for p in 0..1:
      for d in 0..<g.len:
        var
          got = zeroGaugeLike(g)
          delta = zeroGaugeLike(g)
          reference = zeroGaugeLike(g)
        threads:
          for mu in 0..<g.len:
            got[mu] := u[mu]
            reference[mu] := u[mu]
        c.gaugeDerivDeriv2Subset(g, m, delta, p, d)
        threads:
          for mu in 0..<g.len:
            reference[mu] += delta[mu]
        c.gaugeDerivDeriv2SubsetAdd(g, m[d], got, p, d)
        norm2(grt.toGvalue(got) - grt.toGvalue(reference)) :< 1e-24

  test "raw subset Hessian initializes its direct base in one pass":
    let c = GaugeActionCoeffs(plaq: 5.4)
    for p in 0..1:
      let sub = g[0].l.getSubset(if p == 0: "even" else: "odd")
      for d in 0..<g.len:
        var
          got = zeroGaugeLike(g)
          delta = zeroGaugeLike(g)
          reference = zeroGaugeLike(g)
        threads:
          for mu in 0..<g.len:
            reference[mu] := u[mu]
          threadBarrier()
          for x in sub:
            got[d][x] := m[d][x]
            reference[d][x] := m[d][x]
        c.gaugeDerivDeriv2Subset(g, m, delta, p, d)
        threads:
          for mu in 0..<g.len:
            reference[mu] += delta[mu]
        c.gaugeDerivDeriv2SubsetAddBase(g, m[d], u, got, p, d)
        norm2(grt.toGvalue(got) - grt.toGvalue(reference)) :< 1e-24

  test "single-term subset Hessian refreshes its direction":
    let
      c = actWilson(scalar.toGvalue(grt, 5.4))
      zero = grt.toGvalue(zeroGaugeLike(g))
      p = 1
      d = min(2, g.len - 1)
      sub = gaugeActionDeriv(c, gg, p, d)
      got = grad(redot(sub, gm), gg)
    norm2(got - gaugeActionDeriv2(blendSubset(p, d, gm, zero), c, gg)) :< 1e-22
    gm.update(u)
    norm2(got - gaugeActionDeriv2(blendSubset(p, d, gm, zero), c, gg)) :< 1e-22

  test "summed subset Hessian functional clone owns inputs and scratch":
    let
      c = actWilson(scalar.toGvalue(grt, 5.4))
      p = 1
      d = min(2, g.len - 1)
      x = Ggauge(gg.newOneOf)
      forceX = gaugeActionDeriv(c, x, p, d)
      body = grad(redot(forceX, x) + redot(forceX, x), x)
      fn = lambda(x, body)
      cloned = Ggauge(apply(fn, gg))
      forceG = gaugeActionDeriv(c, gg, p, d)
      direct = grad(redot(forceG, gg) + redot(forceG, gg), gg)
    x.update(u)
    discard body.eval
    norm2(cloned - direct) :< 1e-22
    gg.update(u)
    x.update(g)
    discard body.eval
    norm2(cloned - direct) :< 1e-22

  test "wilson force gradient recomp":
    let beta = 5.4
    let c = actWilson(scalar.toGvalue(grt, beta))
    let a = gaugeAction(c, gg)
    let f2 = gaugeForce(c, gg).norm2
    let df2 = grad(f2, gg).norm2
    let rs1 = [a.eval.sval, f2.eval.sval, df2.eval.sval]
    c.updated
    gg.updated
    let rs2 = [a.eval.sval, f2.eval.sval, df2.eval.sval]
    c.updated
    gg.updated
    let rs3 = [a.eval.sval, f2.eval.sval, df2.eval.sval]
    check rs1 == rs2
    check rs1 == rs3

  test "gaugeAction rejects coefficient gradients through the action layer":
    let beta = grt.toGvalue(5.4)
    let c = actWilson(beta)
    expect(GraphValueError):
      discard gaugeAction(c, gg).grad beta

  test "gaugeActionDeriv rejects coefficient gradients through the action layer":
    let beta = grt.toGvalue(5.4)
    let c = actWilson(beta)
    expect(GraphValueError):
      discard gaugeActionDeriv(c, gg).grad beta

  test "gaugeForce rejects coefficient gradients through the action layer":
    let beta = grt.toGvalue(5.4)
    let c = actWilson(beta)
    expect(GraphValueError):
      discard gaugeForce(c, gg).grad beta

  test "gaugeActionDeriv backward rejects missing upstream":
    let beta = grt.toGvalue(5.4)
    let c = actWilson(beta)
    let force = gaugeActionDeriv(c, gg)

    expect(GraphValueError):
      discard force.gfunc.backward(nil, force, 1, gg)

  test "gaugeActionDeriv subset backward rejects invalid paths":
    let
      beta = grt.toGvalue(5.4)
      c = actWilson(beta)
      force = gaugeActionDeriv(c, gg, 1, 2)
    expect(GraphValueError):
      discard force.grad beta
    expect(GraphValueError):
      discard force.gfunc.backward(nil, force, 1, gg)
