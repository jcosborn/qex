# Higher-derivative towers over the grad-complete fallbacks. These build
# replica graphs behind the fused kernels, so they run on a small lattice.
# The including program provides the fixtures (grt, g, u, m, q) and
# `subDir`, a direction valid on its lattice for the subset tests.
#
# Every replica backward is pinned one order past the derivative its hook
# builds, over a cotangent slot that aliases the field and over one that
# merely depends on it: that is the order a frozen seed would silently lose,
# and the orders below it stay exact either way. Those three checks dominate
# the runtime of this file at Nc=3 (seconds to about a minute); trim them
# there before anything else if the cost bites, since tgtoweru1 runs the same
# suite cheaply.

const towerNc = g[0][0].nrows

suite "exp derivative tower":
  setup:
    let gg {.used.} = grt.toGvalue(g)
    let gm {.used.} = grt.toGvalue(m)
    let gq {.used.} = grt.toGvalue(q)

  when towerNc == 1:
    test "scalar exp derivative matches its exact graph identity":
      let x = 5.0 * projTAH(gg)
      norm2(expDeriv(gm, x) - exp(x.adj) * gm) :< 1e-20
  else:
    test "expPolyGraph matches the exp kernel":
      let x = projTAH(gg)
      norm2(exp(x) - expPolyGraph(x)) :< 1e-20
      norm2(grad(redot(gm, exp(x)), gg) - grad(redot(gm, expPolyGraph(x)), gg)) :< 1e-20

  test "expJet matches its basic-op replica":
    # Fused jet kernels (m <= 3) against the nested polynomial replica that
    # also serves past m = 3. Directions are general matrices, not algebra.
    let ds = [gm, gq, gg * gm]
    for m in 1..3:
      norm2(expJet(gg, ds[0..<m]) - expTopReplica(gg, ds[0..<m])) :< 1e-20
    # symmetric in the directions
    norm2(expJet(gg, [gm, gq]) - expJet(gg, [gq, gm])) :< 1e-24

  test "exp second derivative, constant cotangent":
    # grad hits only the x input of the expDeriv node (replica branch).
    proc h(x: Ggauge): Gscalar = redot(gm, exp(projTAH(x)))
    ckgrad(h, gg, gm)
    proc h2(x: Ggauge): Gscalar = redot(gq, grad(h(x), x))
    ckgrad(h2, gg, gm)

  test "exp second derivative, field-dependent cotangent":
    # The product makes the expDeriv b input depend on the field, so grad
    # exercises both the analytic b branch and the replica x branch. Use the
    # small TAH direction: this observable is too nonlinear for the ndiff
    # step along an O(1) link direction.
    proc f(x: Ggauge): Gscalar = redot(gm, exp(projTAH(x)) * exp(projTAH(x)))
    ckgrad(f, gg, gm)
    proc f2(x: Ggauge): Gscalar = redot(gq, grad(f(x), x))
    ckgrad(f2, gg, gm)

  test "exp third derivative":
    proc h(x: Ggauge): Gscalar = redot(gm, exp(projTAH(x)))
    proc h2(x: Ggauge): Gscalar = redot(gq, grad(h(x), x))
    proc h3(x: Ggauge): Gscalar = redot(gm, grad(h2(x), x))
    ckgrad(h3, gg, gm)

  test "expDeriv with one node in both slots":
    proc dup(x: Ggauge): Gscalar =
      let px = projTAH(x)
      redot(gm, expDeriv(px, px))
    ckgrad(dup, gg, gm)
    # One order past the slot the backward hook differentiates: the cotangent
    # slot stays live there, so its cross term survives into this derivative.
    proc dup2(x: Ggauge): Gscalar = redot(gq, grad(dup(x), x))
    ckgrad(dup2, gg, gm)

  test "subset expDeriv tower":
    # Field-dependent, non-commuting cotangent, so the masked analytic branch
    # and the masked replica branch are both exercised. (A cotangent that
    # commutes with the exponent, like exp(px), makes the value the identity
    # and the observable degenerate.)
    proc se(x: Ggauge): Gscalar =
      redot(gm, expDeriv(x * gq, projTAH(x), 1, subDir))
    ckgrad(se, gg, gm)

suite "gauge action derivative tower":
  setup:
    let gg {.used.} = grt.toGvalue(g)
    let gu {.used.} = grt.toGvalue(u)
    let gm {.used.} = grt.toGvalue(m)
    let gq {.used.} = grt.toGvalue(q)
    let c {.used.} = actWilson(scalar.toGvalue(grt, 5.4))

  test "gauge action third derivative":
    proc s1(x: Ggauge): Gscalar = gaugeAction(c, x)
    proc s2(x: Ggauge): Gscalar = redot(gm, grad(s1(x), x))
    ckgrad(s2, gg, gu)
    proc s3(x: Ggauge): Gscalar = redot(gq, grad(s2(x), x))
    ckgrad(s3, gg, gu)

  test "rectangle reference action second derivative":
    # The Hessian kernel is plaquette-only; the reference action carries
    # the rectangle family through basic ops.
    let cr = actSymanzik(scalar.toGvalue(grt, 5.4))
    proc r1(x: Ggauge): Gscalar = gaugeActionGraph(cr, x)
    proc r2(x: Ggauge): Gscalar = redot(gm, grad(r1(x), x))
    ckgrad(r2, gg, gu)
    expect(GraphValueError):
      discard gaugeActionDeriv2(gm, cr, gg).eval

  test "second derivative with field-dependent force direction":
    # b input of gaugeActionDeriv2 depends on the field through norm2's
    # backward, exercising the self-adjoint Hessian branch.
    proc n1(x: Ggauge): Gscalar = norm2(gaugeActionDeriv(c, x))
    ckgrad(n1, gg, gu)
    proc n2(x: Ggauge): Gscalar = redot(gm, grad(n1(x), x))
    ckgrad(n2, gg, gu)
    # A cotangent slot that depends on the field but is not the field: the
    # next order needs that dependence differentiated, not frozen.
    proc n3(x: Ggauge): Gscalar = redot(gq, grad(n2(x), x))
    ckgrad(n3, gg, gu)

  test "subset derivative tower":
    proc p1(x: Ggauge): Gscalar = redot(gm, gaugeActionDeriv(c, x, 1, subDir))
    ckgrad(p1, gg, gu)
    proc p2(x: Ggauge): Gscalar = redot(gq, grad(p1(x), x))
    ckgrad(p2, gg, gu)

  test "subset derivative tower, multi-term force direction":
    proc m1(x: Ggauge): Gscalar =
      let d = gaugeActionDeriv(c, x, 0, 1)
      redot(gm, d) + redot(gq, d)
    ckgrad(m1, gg, gu)
    proc m2(x: Ggauge): Gscalar = redot(gm, grad(m1(x), x))
    ckgrad(m2, gg, gu)

  test "subset derivative tower, field-dependent force direction":
    proc q1(x: Ggauge): Gscalar =
      let d = gaugeActionDeriv(c, x, 0, 1)
      redot(d, d)
    ckgrad(q1, gg, gu)
    proc q2(x: Ggauge): Gscalar = redot(gm, grad(q1(x), x))
    ckgrad(q2, gg, gu)

  test "second derivative with the gauge field as force direction":
    # One node in the b and g slots of gaugeActionDeriv2.
    proc hgg(x: Ggauge): Gscalar = redot(gm, gaugeActionDeriv2(x, c, x))
    ckgrad(hgg, gg, gu)
    # Both slots aliased AND one order further: the g-slot contribution has to
    # stay a live function of the b slot for this to come out right.
    proc hgg2(x: Ggauge): Gscalar = redot(gq, grad(hgg(x), x))
    ckgrad(hgg2, gg, gu)

  test "subset derivative with the gauge field as force direction":
    # The force direction is the gauge field itself, so the b slots alias g in
    # the subset second-derivative backward.
    proc r1(x: Ggauge): Gscalar = redot(x, gaugeActionDeriv(c, x, 1, 0))
    ckgrad(r1, gg, gu)
    proc r2(x: Ggauge): Gscalar = redot(gm, grad(r1(x), x))
    ckgrad(r2, gg, gu)

suite "stout update tower":
  setup:
    let gg {.used.} = grt.toGvalue(g)
    let gu {.used.} = grt.toGvalue(u)
    let gm {.used.} = grt.toGvalue(m)
    let gq {.used.} = grt.toGvalue(q)
    let alpha {.used.} = grt.toGvalue(0.1)
    let c {.used.} = actWilson(scalar.toGvalue(grt, 5.4))

  test "stout update second derivative, field-dependent staple":
    # The staple sum makes both pullback slots of the update kernel live.
    # (A staple like x*gu degenerates for U(1): W ds^dag = |x|^2 gu^dag.)
    # Too nonlinear for the ndiff step along an O(1) link direction; use
    # the small TAH direction as the exp tower does.
    proc s1(x: Ggauge): Gscalar = redot(gm, stoutUpdate(x, gaugeActionDeriv(c, x), alpha, 1, subDir))
    ckgrad(s1, gg, gm)
    proc s2(x: Ggauge): Gscalar = redot(gq, grad(s1(x), x))
    ckgrad(s2, gg, gm)

  test "stout update alpha derivative differentiates in the field":
    proc s1(x: Ggauge): Gscalar = redot(gm, stoutUpdate(x, gu, alpha, 0, subDir))
    proc a1(x: Ggauge): Gscalar = grad(s1(x), alpha)
    ckgrad(a1, gg, gm)
