## Stout FTHMC tests: SU(3) here, U(1) through tftstoutu1.nim.
## Check graph gradients against ndiff, the small-rho log-Jacobian, and rho=0.

import qex, algorithms/numdiff, maths/groupOps
import ../[core, scalar, gauge]
import ../functional
import ../hmcgauge/ftstout

proc runFtStoutTests*(lat: seq[int]; beta, rho: float; nsmear = 1) =
  qexInit()
  const eps = 1e-3
  let seed = 1234567891'u
  let
    grt = initGraphRuntime()
    lo = lat.newLayout
    gc = actWilson(scalar.toGvalue(grt, beta))
    c1 = actWilson(scalar.toGvalue(grt, 1.0))
    rhoG = scalar.toGvalue(grt, rho)
  var
    r = lo.newRNGField(Philox4x64, seed)
    V0 = lo.newgauge
    Aconst = lo.newgauge
    Bconst = lo.newgauge
  V0.random r
  Aconst.random r
  Bconst.random r
  let
    Ag = gauge.toGvalue(grt, Aconst)
    Bg = gauge.toGvalue(grt, Bconst)
  const nc = V0[0][0].nrows
  echo "=== ftstout tests: Nc=", nc, " ", lat.len, "D lat=", lat,
    " beta=", beta, " rho=", rho, " ==="

  var nfail = 0
  proc checkGaugeEq(name: string, x, y: Ggauge; tol = 1e-18) =
    let
      d = (x - y).norm2.eval.sval
      n = y.norm2.eval.sval
      ok = d <= tol*(1.0 + n)
    echo "E[", name, "]: |dx|^2=", d, " ref=", n
    if not ok: inc nfail

  proc checkScalarEq(name: string, x, y: Gscalar; tol = 1e-10) =
    let
      xv = x.eval.sval
      yv = y.eval.sval
      rel = abs(xv - yv)/(1.0 + abs(yv))
    echo "E[", name, "]: x=", xv, " ref=", yv, " rel=", rel
    if rel > tol: inc nfail

  proc checkNodeGrad(name: string, f, t: Gscalar, x0: float) =
    let ana = f.grad(t).eval.sval
    proc act(x: float): float =
      t.update x
      f.eval.sval
    var nd, err: float
    ndiff(nd, err, act, x0, eps, ordMax = 3)
    t.update x0
    let
      delta = abs(nd - ana)
      rel = delta/(abs(nd) + abs(ana) + 1e-30)
      ok = rel < 1e-6 or delta < max(1e-9, 32.0*err)
    echo "N[", name, "]: ana=", ana, " num=", nd,
      " +/- ", err, " rel=", rel
    if not ok: inc nfail

  block:
    let Vg = gauge.toGvalue(grt, V0)
    var rejected = false
    try:
      discard smearedField(Vg, rho, -1)
    except GraphValueError:
      rejected = true
    doAssert rejected

    rejected = false
    try:
      discard stoutAction(gc, rho, -1)
    except GraphValueError:
      rejected = true
    doAssert rejected

    var tmp = lo.newgauge
    rejected = false
    try:
      discard invertStoutFlow(tmp, rho, 1, maxiter = 0)
    except GraphValueError:
      rejected = true
    doAssert rejected

    threads:
      for mu in 0..<tmp.len:
        tmp[mu] := V0[mu]
    rejected = false
    try:
      discard invertStoutFlow(tmp, rho, 1, maxiter = 1)
    except GraphValueError:
      rejected = true
    doAssert rejected

  block:
    let
      Vg = gauge.toGvalue(grt, V0)
      u = smearFlow(Vg, c1, rhoG, 0)
      ld = logDetJ(u, Vg)
      seff = stoutAction(gc, rho, 0).action(Vg)
      plain = gaugeAction(gc, Vg)
    doAssert u.nodeKey == Vg.nodeKey
    doAssert ld.isStaticZeroLeaf
    checkScalarEq("stoutAction zero sweep", seff, plain)

  proc checkGrad(name: string, build: proc(Vt: Ggauge): Gscalar) =
    ## Perturb V along a random algebra direction R via V(t)=exp(t R)V0 and compare
    ## the graph gradient d/dt with the numerical derivative of the forward.
    var R = lo.newgauge
    R.randomTAH r
    let
      V0g = gauge.toGvalue(grt, V0)
      Rg = gauge.toGvalue(grt, R)
      gt = scalar.toGvalue(grt, 0.0)
      Vt = axexpmuly(gt, Rg, V0g)
      f = build(Vt)
      dsdt = f.grad gt
    let ana = dsdt.eval.sval
    proc act(t: float): float =
      gt.update t
      f.eval.sval
    var nd, err: float
    ndiff(nd, err, act, 0.0, eps, ordMax = 3)
    gt.update 0.0
    let rel = abs((nd - ana) / (abs(nd) + abs(ana) + 1e-30))
    echo "A[", name, "]: ana=", ana, " num=", nd, " +/- ", err, " rel=", rel
    if rel >= 1e-6: inc nfail

  checkGrad("ax_xgauge", proc(Vt: Ggauge): Gscalar =
    # x = projTAH(ds Vt†) depends on Vt: isolates expDeriv (Nc=1 fix guard)
    gaugeAction(gc, axexpmuly(rhoG, contractProjTAH(Vt, Ag), Ag)))
  checkGrad("blendW", proc(Vt: Ggauge): Gscalar =
    # Vt on the complement of (even, dir0): blendSubset W-slot backward
    gaugeAction(gc, blendSubset(0, 0, Ag, Vt)))
  checkGrad("blendCand", proc(Vt: Ggauge): Gscalar =
    # Vt on the (even, dir0) subset: blendSubset cand-slot backward
    gaugeAction(gc, blendSubset(0, 0, Vt, Ag)))
  checkGrad("lndet", proc(Vt: Ggauge): Gscalar =
    logDetJ(smearFlow(Vt, c1, rhoG, nsmear), Vt))
  checkGrad("seff", proc(Vt: Ggauge): Gscalar =
    # the full S_eff whose grad is the integrator force
    let u = smearFlow(Vt, c1, rhoG, nsmear)
    gaugeAction(gc, u) - logDetJ(u, Vt))

  proc checkAlphaGrad(name: string, build: proc(a: Gscalar): Gscalar) =
    let
      a = scalar.toGvalue(grt, rho)
      f = build(a)
      ana = f.grad(a).eval.sval
    proc act(x: float): float =
      a.update x
      f.eval.sval
    var nd, err: float
    ndiff(nd, err, act, rho, eps, ordMax = 3)
    a.update rho
    let rel = abs((nd - ana) / (abs(nd) + abs(ana) + 1e-30))
    echo "A[", name, "]: ana=", ana, " num=", nd, " +/- ", err, " rel=", rel
    if rel >= 1e-6: inc nfail

  checkAlphaGrad("lndet_rho", proc(a: Gscalar): Gscalar =
    let Vg = gauge.toGvalue(grt, V0)
    logDetJ(smearFlow(Vg, c1, a, nsmear), Vg))
  checkAlphaGrad("seff_rho", proc(a: Gscalar): Gscalar =
    let Vg = gauge.toGvalue(grt, V0)
    let u = smearFlow(Vg, c1, a, nsmear)
    gaugeAction(gc, u) - logDetJ(u, Vg))

  # Direct fused/reference test for one odd-parity subset. The split objective
  # also checks a summed upstream and reevaluation after an alpha update.
  block:
    let
      parity = 1
      dir = min(1, V0.len - 1)
      W = gauge.toGvalue(grt, V0)
      ds = gauge.toGvalue(grt, Aconst)
      a = scalar.toGvalue(grt, rho)
      F = contractProjTAH(W, ds, parity, dir)
      cand = axexpmuly(a, F, W, parity, dir)
      separate = blendSubset(parity, dir, cand, W)
      fused = stoutUpdate(W, ds, a, parity, dir)
      sf = redot(fused, Bg) + redot(fused, Ag)
      sr = redot(separate, Bg) + redot(separate, Ag)
      ddsf = grad(sf, ds)
      ddsr = grad(sr, ds)
      dWf = grad(sf, W)
      dWr = grad(sr, W)
      daf = grad(sf, a)
      dar = grad(sr, a)
    checkGaugeEq("stoutUpdate value", fused, separate)
    checkGaugeEq("stoutUpdate ds grad", ddsf, ddsr)
    checkGaugeEq("stoutUpdate W grad", dWf, dWr)
    checkScalarEq("stoutUpdate alpha grad", daf, dar)
    a.update(0.8*rho)
    checkGaugeEq("stoutUpdate value update", fused, separate)
    checkGaugeEq("stoutUpdate ds grad update", ddsf, ddsr)
    checkGaugeEq("stoutUpdate W grad update", dWf, dWr)
    checkScalarEq("stoutUpdate alpha grad update", daf, dar)
    a.update(rho)

  # The grouped step must match the independent ops while merging gauge and
  # scalar cotangents into one subset-gradient kernel.
  block:
    let
      parity = 1
      dir = min(1, V0.len - 1)
      W = gauge.toGvalue(grt, V0)
      ds = gauge.toGvalue(grt, Aconst)
      a = scalar.toGvalue(grt, rho)
      grouped = stoutUpdateLogDetJ(W, ds, a, parity, dir)
      updateRef = stoutUpdate(W, ds, a, parity, dir)
      logdetRef = stoutLogDetJ(W, ds, a, parity, dir)
      fg = redot(grouped.Wnew, Bg) + redot(grouped.Wnew, Ag) - grouped.lj
      fr = redot(updateRef, Bg) + redot(updateRef, Ag) - logdetRef
      dWg = grad(fg, W)
      dWr = grad(fr, W)
      ddsg = grad(fg, ds)
      ddsr = grad(fr, ds)
      dag = grad(fg, a)
      dar = grad(fr, a)
    checkScalarEq("stoutStep logdet", grouped.lj, logdetRef)
    checkGaugeEq("stoutStep update", grouped.Wnew, updateRef)
    checkScalarEq("stoutStep objective", fg, fr)
    checkGaugeEq("stoutStep W grad", dWg, dWr)
    checkGaugeEq("stoutStep ds grad", ddsg, ddsr)
    checkScalarEq("stoutStep alpha grad", dag, dar)
    a.update(0.8*rho)
    checkGaugeEq("stoutStep update refresh", grouped.Wnew, updateRef)
    checkScalarEq("stoutStep logdet refresh", grouped.lj, logdetRef)
    checkGaugeEq("stoutStep W grad refresh", dWg, dWr)
    checkGaugeEq("stoutStep ds grad refresh", ddsg, ddsr)
    checkScalarEq("stoutStep alpha grad refresh", dag, dar)
    a.update(rho)

  # logDetJ must match the hand-accumulated substep sum and its gradients.
  block:
    let
      W0 = gauge.toGvalue(grt, V0)
      a = scalar.toGvalue(grt, rho)
      nd = V0.len
    var
      W = W0
      steps: seq[tuple[Wnew: Ggauge, lj: Gscalar]]
      manual = scalar.toGvalue(grt, 0.0)
    for parity in 0..1:
      for dir in 0..<nd:
        let st = stoutUpdateLogDetJ(W, c1, a, parity, dir)
        W = st.Wnew
        steps.add st
        manual = manual + st.lj
    let auto = logDetJ(W, W0)
    checkScalarEq("logDetJ fused chain", auto, manual)
    doAssert logDetJ(W, W0).nodeKey == auto.nodeKey   # cached: the same node
    # A one-step chain is the step's own fused logdet view, not a duplicate.
    doAssert logDetJ(steps[0].Wnew, W0).nodeKey == steps[0].lj.nodeKey
    # An intermediate flow node is a valid base.
    checkScalarEq("logDetJ tail chain",
      logDetJ(W, steps[0].Wnew) + steps[0].lj, manual)
    checkGaugeEq("logDetJ fused chain grad", grad(auto, W0), grad(manual, W0))
    checkScalarEq("logDetJ fused chain alpha grad", grad(auto, a), grad(manual, a))
    a.update(0.8*rho)
    checkScalarEq("logDetJ fused chain refresh", auto, manual)
    a.update(rho)
    # cond between two flows: the sum is piecewise and follows the selector.
    let
      k = scalar.toGvalue(grt, 1.0)
      condLd = logDetJ(cond(k, W, steps[0].Wnew), W0)
    checkScalarEq("logDetJ cond flow", condLd, auto)
    k.update 0.0
    checkScalarEq("logDetJ cond flow flip", condLd, steps[0].lj)
    k.update 1.0
    # Plain stoutUpdate steps declare the same factorization.
    let
      ds = gauge.toGvalue(grt, Aconst)
      u1 = stoutUpdate(W0, ds, a, 1, min(1, nd - 1))
      ldRef = stoutLogDetJ(W0, ds, a, 1, min(1, nd - 1))
    checkScalarEq("logDetJ stoutUpdate", logDetJ(u1, W0), ldRef)

  # A local stout determinant holds free auxiliaries fixed. Do not advertise it
  # as a total factorization when those auxiliaries feed back from the flow
  # input; the action-aware staple is the explicit frozen-context exception.
  block:
    let
      parity = 1
      dir = min(1, V0.len - 1)
      W = gauge.toGvalue(grt, V0)
      ds = gauge.toGvalue(grt, Aconst)
      a = scalar.toGvalue(grt, rho)

    proc rejectsLogdet(u: Ggauge): bool =
      try:
        discard logDetJ(u, W)
      except GraphValueError:
        return true
      false

    let aliased = stoutUpdate(W, W, a, parity, dir)
    when nc == 1:
      # projectTAH(W W†) is zero, so the composed U(1) map is the identity.
      checkGaugeEq("reentrant stout identity", aliased, W)
    doAssert rejectsLogdet(aliased)
    doAssert rejectsLogdet(
      stoutUpdateLogDetJ(W, W, a, parity, dir).Wnew)

    let fieldDs = W + ds
    doAssert rejectsLogdet(stoutUpdate(W, fieldDs, a, parity, dir))

    let fieldAlpha = redot(W, Ag)
    doAssert rejectsLogdet(stoutUpdate(W, ds, fieldAlpha, parity, dir))
    doAssert rejectsLogdet(
      stoutUpdateLogDetJ(W, c1, fieldAlpha, parity, dir).Wnew)

    let fieldCoefficients = actWilson(redot(W, Ag))
    doAssert rejectsLogdet(
      stoutUpdateLogDetJ(W, fieldCoefficients, a, parity, dir).Wnew)

  # logDetJ sums participate in functional cloning like ordinary graph values.
  block:
    let
      W = gauge.toGvalue(grt, V0)
      p = Ggauge(W.newOneOf)
      a = scalar.toGvalue(grt, rho)
      stepP0 = stoutUpdateLogDetJ(p, c1, a, 0, 0)
      stepP1 = stoutUpdateLogDetJ(stepP0.Wnew, c1, a, 1, min(1, V0.len - 1))
      body = redot(stepP1.Wnew, Bg) - logDetJ(stepP1.Wnew, p)
      fn = lambda(p, body)
      cloned = Gscalar(apply(fn, W))
      stepW0 = stoutUpdateLogDetJ(W, c1, a, 0, 0)
      stepW1 = stoutUpdateLogDetJ(stepW0.Wnew, c1, a, 1, min(1, V0.len - 1))
      direct = redot(stepW1.Wnew, Bg) - logDetJ(stepW1.Wnew, W)
      clonedGrad = grad(cloned, W)
      directGrad = grad(direct, W)
    p.update(Aconst)
    discard cloned.eval
    discard body.eval
    checkScalarEq("logDetJ functional clone", cloned, direct)
    checkGaugeEq("logDetJ functional clone grad", clonedGrad, directGrad)

  # The action-aware overload folds the internally constructed staple pullback
  # into W while preserving the explicit-staple operation as a reference.
  block:
    let
      parity = 1
      dir = min(1, V0.len - 1)
      W = gauge.toGvalue(grt, V0)
      up = gauge.toGvalue(grt, Bconst)
      a = scalar.toGvalue(grt, rho)
      ds = gaugeActionDeriv(c1, W, parity, dir)
      fused = stoutUpdateLogDetJ(W, c1, a, parity, dir)
      reference = stoutUpdateLogDetJ(W, ds, a, parity, dir)
      ff = redot(fused.Wnew, up) + redot(fused.Wnew, Ag) - fused.lj
      fr = redot(reference.Wnew, up) + redot(reference.Wnew, Ag) - reference.lj
      dWf = grad(ff, W)
      dWr = grad(fr, W)
      daf = grad(ff, a)
      dar = grad(fr, a)
      dWuf = grad(redot(fused.Wnew, up), W)
      dWur = grad(redot(reference.Wnew, up), W)
      dWlf = grad(fused.lj, W)
      dWlr = grad(reference.lj, W)
      dalf = grad(fused.lj, a)
      dalr = grad(reference.lj, a)
    checkScalarEq("stoutActionStep logdet", fused.lj, reference.lj)
    checkGaugeEq("stoutActionStep update", fused.Wnew, reference.Wnew)
    checkGaugeEq("stoutActionStep W grad", dWf, dWr)
    checkScalarEq("stoutActionStep alpha grad", daf, dar)
    checkGaugeEq("stoutActionStep update-only W grad", dWuf, dWur)
    checkGaugeEq("stoutActionStep log-only W grad", dWlf, dWlr)
    checkScalarEq("stoutActionStep log-only alpha grad", dalf, dalr)
    var coeffRejected = false
    try:
      discard grad(ff, c1)
    except GraphValueError:
      coeffRejected = true
    doAssert coeffRejected
    W.update(Aconst)
    up.update(V0)
    a.update(0.8*rho)
    checkGaugeEq(
      "stoutActionStep update refresh", fused.Wnew, reference.Wnew)
    checkScalarEq(
      "stoutActionStep logdet refresh", fused.lj, reference.lj)
    checkGaugeEq("stoutActionStep W grad refresh", dWf, dWr)
    checkGaugeEq(
      "stoutActionStep update-only W grad refresh", dWuf, dWur)
    checkGaugeEq("stoutActionStep log-only W grad refresh", dWlf, dWlr)
    checkScalarEq(
      "stoutActionStep log-only alpha grad refresh", dalf, dalr)
    checkScalarEq("stoutActionStep alpha grad refresh", daf, dar)

  # Functional clones own their update cache and restore zero-copy output views.
  block:
    let
      parity = 1
      dir = min(1, V0.len - 1)
      W = gauge.toGvalue(grt, V0)
      p = Ggauge(W.newOneOf)
      ds = gauge.toGvalue(grt, Bconst)
      a = scalar.toGvalue(grt, rho)
      stepP = stoutUpdateLogDetJ(p, ds, a, parity, dir)
      body = redot(stepP.Wnew, Bg) - stepP.lj
      fn = lambda(p, body)
      cloned = Gscalar(apply(fn, W))
      stepW = stoutUpdateLogDetJ(W, ds, a, parity, dir)
      direct = redot(stepW.Wnew, Bg) - stepW.lj
      clonedGrad = grad(cloned, W)
      directGrad = grad(direct, W)
    p.update(Aconst)
    discard cloned.eval
    discard body.eval
    checkScalarEq("stoutStep functional clone", cloned, direct)
    checkGaugeEq("stoutStep functional clone grad", clonedGrad, directGrad)
    W.update(Aconst)
    p.update(V0)
    discard body.eval
    checkScalarEq("stoutStep functional clone refresh", cloned, direct)
    checkGaugeEq(
      "stoutStep functional clone grad refresh", clonedGrad, directGrad)

  # The action-aware clone also owns its one-direction Hessian scratch.
  block:
    let
      parity = 1
      dir = min(1, V0.len - 1)
      W = gauge.toGvalue(grt, V0)
      p = Ggauge(W.newOneOf)
      a = scalar.toGvalue(grt, rho)
      stepP = stoutUpdateLogDetJ(p, c1, a, parity, dir)
      body = redot(stepP.Wnew, Bg) - stepP.lj
      fn = lambda(p, body)
      cloned = Gscalar(apply(fn, W))
      stepW = stoutUpdateLogDetJ(W, c1, a, parity, dir)
      direct = redot(stepW.Wnew, Bg) - stepW.lj
      clonedGrad = grad(cloned, W)
      directGrad = grad(direct, W)
    p.update(Aconst)
    discard cloned.eval
    discard body.eval
    checkScalarEq("stoutActionStep functional clone", cloned, direct)
    checkGaugeEq(
      "stoutActionStep functional clone grad", clonedGrad, directGrad)
    W.update(Aconst)
    p.update(V0)
    discard body.eval
    checkScalarEq(
      "stoutActionStep functional clone refresh", cloned, direct)
    checkGaugeEq(
      "stoutActionStep functional clone grad refresh", clonedGrad, directGrad)

  # Check stoutLogDetJ independently of the flow composition: exact forward
  # value and numerical derivatives of each free input.
  block:
    let
      parity = 1
      dir = min(1, V0.len - 1)
      sub = lo.getSubset("odd")
      W = gauge.toGvalue(grt, V0)
      ds = gauge.toGvalue(grt, Aconst)
      a = scalar.toGvalue(grt, rho)
      lj = stoutLogDetJ(W, ds, a, parity, dir)
    var exact = 0.0
    threads:
      var s = 0.0
      for e in sub:
        var M {.noinit.}: evalType(V0[dir][e])
        M := rho * (V0[dir][e] * Aconst[dir][e].adj)
        when nc == 1:
          s += simdSum(ln(1.0 + M[0, 0].re))
        elif nc == 3:
          s += simdSum(expProjMulLogJac(M[]))
      s.threadRankSum
      threadSingle: exact = s
    let got = lj.eval.sval
    echo "E[stoutLogDetJ value]: x=", got, " ref=", exact
    if abs(got - exact) > 1e-12*(1.0 + abs(exact)): inc nfail

    var
      Rw = lo.newgauge
      Rd = lo.newgauge
    Rw.random r
    Rd.random r
    let
      Rwg = gauge.toGvalue(grt, Rw)
      Rdg = gauge.toGvalue(grt, Rd)
      tw = scalar.toGvalue(grt, 0.0)
      td = scalar.toGvalue(grt, 0.0)
      fw = stoutLogDetJ(W + tw*Rwg, ds, a, parity, dir)
      fd = stoutLogDetJ(W, ds + td*Rdg, a, parity, dir)
    checkNodeGrad("stoutLogDetJ W", fw, tw, 0.0)
    checkNodeGrad("stoutLogDetJ ds", fd, td, 0.0)
    checkNodeGrad("stoutLogDetJ alpha", lj, a, rho)

    let
      lj2 = stoutLogDetJ(W, ds, a, parity, dir)
      dd0 = grad(lj, ds)
      dW0 = grad(lj, W)
      dW1 = grad(lj2, W)
      dd1 = grad(lj2, ds)
    checkGaugeEq("stoutLogDetJ ds grad order", dd0, dd1)
    checkGaugeEq("stoutLogDetJ W grad order", dW0, dW1)
    a.update(0.8*rho)
    checkGaugeEq("stoutLogDetJ ds grad update", dd0, dd1)
    checkGaugeEq("stoutLogDetJ W grad update", dW0, dW1)
    a.update(rho)

  block:
    var lat2 = lat
    lat2[^1] *= 2
    let
      lo2 = lat2.newLayout
      other = gauge.toGvalue(grt, lo2.newgauge)
      W = gauge.toGvalue(grt, V0)
      a = scalar.toGvalue(grt, rho)
    var updateRejected = false
    try:
      discard stoutUpdate(W, other, a, 0, 0)
    except GraphValueError:
      updateRejected = true
    doAssert updateRejected
    var logdetRejected = false
    try:
      discard stoutLogDetJ(W, other, a, 0, 0)
    except GraphValueError:
      logdetRejected = true
    doAssert logdetRejected

  doAssert nfail == 0

  # The action and measurement APIs must share one flow for the exact same graph
  # node, while distinct graph leaves receive distinct cached flow nodes.
  block:
    let
      Vg = gauge.toGvalue(grt, V0)
      sa = stoutAction(gc, rho, nsmear)
      f0 = sa.flow(Vg)
      f1 = sa.flow(Vg)
      seff = sa.action(Vg)
      f2 = sa.flow(Vg)
      Vother = gauge.toGvalue(grt, V0)
      fo = sa.flow(Vother)
      direct = gaugeAction(gc, f0) - logDetJ(f0, Vg)
    doAssert f0.nodeKey == f1.nodeKey
    doAssert f0.nodeKey == f2.nodeKey
    doAssert f0.nodeKey != fo.nodeKey
    # The action evaluates the same cached logDetJ node returned for the flow.
    doAssert seff.reaches(logDetJ(f0, Vg), iwmEval)
    doAssert abs(seff.eval.sval - direct.eval.sval) < 1e-12
    Vg.update(Aconst)
    let
      refreshed = sa.flow(Vg)
      fresh = smearedField(Vg, rho, nsmear)
    doAssert refreshed.nodeKey == f0.nodeKey
    checkGaugeEq("stoutAction cache refresh", refreshed, fresh)
    checkScalarEq("stoutAction logdet refresh",
      logDetJ(refreshed, Vg), logDetJ(fresh, Vg))
    doAssert nfail == 0

  # --- B. log-Jacobian forward leading order at small rho --------------------
  block:
    let smallRho = 1e-4
    # Leading order of ln det f' for the canonical projTAH(W ds†) direction:
    # U(1) → +rho*redot(V,ds); SU(N) → +rho*(N^2-1)/N*redot(V,ds) (trace part of dF).
    const linCoef = when nc == 1: 1.0 else: float(nc*nc - 1) / float(nc)
    let
      V0g = gauge.toGvalue(grt, V0)
      u = smearFlow(V0g, c1, scalar.toGvalue(grt, smallRho), 1)
      lin = scalar.toGvalue(grt, linCoef*smallRho) * redot(V0g, gaugeActionDeriv(c1, V0g))
      lndetVal = logDetJ(u, V0g).eval.sval
      linVal = lin.eval.sval
    let rel = abs((lndetVal - linVal) / (abs(lndetVal) + abs(linVal) + 1e-30))
    echo "B lndet: actual=", lndetVal, " linear=", linVal, " rel=", rel
    doAssert rel < 1e-2

  # --- C. rho = 0 identity baseline ------------------------------------------
  block:
    let
      V0g = gauge.toGvalue(grt, V0)
      u = smearFlow(V0g, c1, scalar.toGvalue(grt, 0.0), nsmear)
      seff = gaugeAction(gc, u) - logDetJ(u, V0g)
      sPlain = gaugeAction(gc, V0g)
      fForce = contractProjTAH(grad(seff, V0g), V0g)
      gForce = gaugeForce(gc, V0g)
    let
      lndet0 = logDetJ(u, V0g).eval.sval
      dS = abs(seff.eval.sval - sPlain.eval.sval)
      fNorm = gForce.norm2.eval.sval
      fDiff = (fForce - gForce).norm2.eval.sval
    echo "C alpha=0: lndet=", lndet0, " dSeff=", dS, " |dForce|^2=", fDiff
    doAssert abs(lndet0) < 1e-12
    doAssert dS < 1e-8 * (1.0 + abs(sPlain.eval.sval))
    doAssert fDiff < 1e-16 * (1.0 + fNorm)

  # --- D. inverse round trip: invertStoutFlow(smearFlow(V0)) recovers V0 ------
  block:
    let smeared = smearFlow(gauge.toGvalue(grt, V0), c1, rhoG, nsmear)
    discard smeared.eval
    let Usnap = smeared.gaugeSnapshot          # physical U = f(V0)
    var Uinv = lo.newgauge
    threads:
      for mu in 0..<Uinv.len: Uinv[mu] := Usnap[mu]
    let inv = invertStoutFlow(Uinv, rho, nsmear)   # U -> V in place
    var dV, dU = 0.0
    # recovered field vs the original V0 (forward is a diffeomorphism at this rho)
    threads:
      var m = 0.0
      for mu in 0..<Uinv.len:
        for x in Uinv[mu]:
          let d = (Uinv[mu][x] - V0[mu][x]).norm2.simdMax
          if m < d: m = d
      m.threadRankMax
      threadSingle: dV = m
    # forward of the recovered field vs the physical U (independent of injectivity)
    let reSmeared = smearFlow(gauge.toGvalue(grt, Uinv), c1, rhoG, nsmear)
    discard reSmeared.eval
    let reSnap = reSmeared.gaugeSnapshot
    threads:
      var m = 0.0
      for mu in 0..<reSnap.len:
        for x in reSnap[mu]:
          let d = (reSnap[mu][x] - Usnap[mu][x]).norm2.simdMax
          if m < d: m = d
      m.threadRankMax
      threadSingle: dU = m
    echo "D inverse: iter=", inv.iter, " rdf2=", inv.rdf2,
      " max|V'-V0|^2=", dV, " max|f(V')-U|^2=", dU
    doAssert dU < 1e-16
    doAssert dV < 1e-16

  echo "tftstout OK (Nc=", nc, ", ", lat.len, "D)"
  qexFinalize()

when isMainModule:
  runFtStoutTests(@[4, 4, 4, 4], 6.0, 0.02)
