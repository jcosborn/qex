import qex
import maths/groupOps
import testutils, scaledexpRef
import std/math

type
  M = MatrixArray[3,3,ComplexType[float64]]
  A = MatrixArray[8,8,float64]
  V = VectorArray[8,float64]

let refs=buildRefs()

proc rel(a, b: auto): float =
  sqrt(norm2(a-b))/max(1.0, sqrt(norm2(b)))

proc scalarRel(a, b: float): float =
  abs(a-b)/max(1.0, abs(b))

proc finiteJet(p, dp: var A, x, dx: A, order, scale: int) =
  # Direct powers provide a reference for the structured SU(3) evaluation.
  let s = 1.0 / float(1 shl scale)
  let y = s*x
  let dy = s*dx
  var t, dt: A
  p := 1.0
  dp := 0.0
  t := 1.0
  dt := 0.0
  for k in 1..order:
    dt := (1.0/float(k+1))*(dt*y + t*dy)
    t := (1.0/float(k+1))*(t*y)
    p += t
    dp += dt
  for j in 0..<scale:
    let c = 1.0 / float(1 shl (scale+1-j))
    let pp = p*p
    dp := dp + c*(dx*pp + x*(dp*p + p*dp))
    p := p + c*(x*pp)

proc poly12(e, l: var M, f, d: M): int =
  # The branch belongs to the primal. L_next = E L + L E.
  var n2 = norm2(f)
  var s = 1.0
  while n2 > 1.0/16.0:
    n2 *= 0.25
    s *= 0.5
    inc result
  let y = s*f
  let dy = s*d
  var t, dt: M
  e := 1.0
  l := 0.0
  t := 1.0
  dt := 0.0
  for k in 1..12:
    dt := (1.0/float(k))*(dt*y + t*dy)
    t := (1.0/float(k))*(t*y)
    e += t
    l += dt
  for k in 0..<result:
    l := e*l + l*e
    e := e*e

proc checkFinite(scale: static int) =
  var mp, md, mg, ml, ma, mt, mi, mc: float
  for c in refs.cases:
    let
      m = c.m
      dm = c.dm
      cot = c.c
      x = c.x
      dx = c.dx
      b = c.b
      v = c.v
    var bx, bd: A
    bx.su3AdNeg(m)
    bd.su3ProjectDeriv(m)
    check rel(bx, x) < 2e-14
    check rel(bd, c.d) < 2e-14
    # Every CH perturbation is induced by a 3x3 input direction.
    bx.su3AdNeg(dm)
    check rel(bx, dx) < 2e-14
    check norm2(b-b.adj) > 1.0
    for r in c.finite:
      if r.scale != scale: continue
      let order = r.order
      checkpoint c.name & " order=" & $order & " scale=" & $scale
      let
        rp = r.phi
        rd = r.dphi
        rg = r.grad
      var j, p, df, ad, jp, dp: A
      var f, g, gc, gm, pb, pm: M
      var pv: V
      j.diffExpProjectTAHMul(p, df, ad, f, m, order=order, scale=scale)
      finiteJet(jp, dp, x, dx, order, scale)
      g.expProjMulLogJacGrad(m, order=order, scale=scale)
      gc.expProjMulLogJacGrad(pv, m, v, order=order, scale=scale)
      gm.expProjMulLogJacGrad(pm, m, cot, order=order, scale=scale)
      pb.expProjectTAHPullback(m, cot, order=order, scale=scale)
      let
        pe = rel(p, rp)
        de = rel(dp, rd)
        ge = rel(g, rg)
        le = scalarRel(expProjMulLogJac(m, order=order, scale=scale), r.log)
        ae = max(rel(pv, r.apply), rel(pb, r.pullback))
        te = scalarRel(redot(b, dp), r.pair)
      mp = max(mp, pe)
      md = max(md, de)
      mg = max(mg, ge)
      ml = max(ml, le)
      ma = max(ma, ae)
      mt = max(mt, te)
      mi = max(mi, r.invNorm)
      mc = max(mc, r.cond)
      check pe < 2e-13
      check rel(jp, rp) < 2e-13
      check rel(j, r.k) < 2e-13
      check de < 2e-12
      check te < 2e-12
      check ge < 2e-11
      check le < 2e-11
      check ae < 2e-12
      check rel(gc, rg) < 2e-11
      check rel(gm, rg) < 2e-11
      check rel(pm, r.pullback) < 2e-12
      check scalarRel(redot(g, dm), r.dlog) < 2e-12
      check r.invNorm < 60.0
      when scale == 0:
        if order == 13:
          var old: M
          old.expProjMulLogJacGrad(m)
          check norm2(old-g) == 0.0
          check expProjMulLogJac(m) == expProjMulLogJac(m, order=13, scale=0)
  echo "finite Phi scale=", scale, " max normalized value=", mp,
    " derivative=", md, " gradient=", mg, " logdet=", ml,
    " application=", ma, " nonsymmetric cotangent pairing=", mt,
    " max KinvNorm=", mi, " condF=", mc

suite "scaled SU3 differential references":
  test "finite value and gradient at scale zero":
    checkFinite(0)

  test "finite value and gradient at scale five":
    checkFinite(5)

  test "log Jacobian stays finite for a large positive Hermitian part":
    var m: M
    m := 1e40
    # F = 0, D = 1e40 I, K = (1 + 1e40) I.
    let want = 8.0*ln(1e40)
    check scalarRel(expProjMulLogJac(m), want) < 2e-14
    check scalarRel(expProjMulLogJac(m, scale=5), want) < 2e-14

  test "scaled Phi versus analytic Phi and selected Poly12 logdet":
    var mp, md, ml, mj, mi, mc: float
    for c in refs.cases:
      checkpoint c.name
      let m = c.m
      var j, p, df, ad: A
      var f, g: M
      j.diffExpProjectTAHMul(p, df, ad, f, m, order=13, scale=5)
      g.expProjMulLogJacGrad(m, order=13, scale=5)
      let
        pe = rel(p, c.analyticPhi)
        de = scalarRel(redot(g, c.dm), c.analyticDLog)
        ll = expProjMulLogJac(m, order=13, scale=5)
        le = max(scalarRel(ll, c.analyticLog), scalarRel(ll, c.selectedLog))
        je = scalarRel(determinant(j), exp(c.selectedLog))
      mp = max(mp, pe)
      md = max(md, de)
      ml = max(ml, le)
      mj = max(mj, je)
      mi = max(mi, c.analyticInvNorm)
      mc = max(mc, c.analyticCond)
      echo c.name, " Fnorm=", c.fNorm,
        " Dnorm=", c.dNorm, " KinvNorm=", c.analyticInvNorm,
        " condF=", c.analyticCond
      check pe < 2e-13
      check de < 2e-12
      check le < 2e-9
      check je < 2e-9
    echo "analytic Phi max normalized value=", mp, " logdet derivative=", md,
      " analytic/selected Poly12 logdet=", ml, " selected determinant=", mj,
      " max KinvNorm=", mi, " condF=", mc

  test "selected Poly12 derivative follows the primal branch through norm eight":
    var mv, mj, mf, ma: float
    for c in refs.exp:
      checkpoint c.name
      let
        f = c.f
        dm = c.dm
        a = c.alpha
        cot = c.c
      var d, e, l, ea, la, actual, pb: M
      d.projectTAH(dm)
      let ns = poly12(e, l, f, d)
      discard poly12(ea, la, f, a)
      actual := expAH(f)
      pb.expProjectTAHPullback(f, actual.adj*cot, order=13, scale=5)
      let
        ve = rel(actual, c.e)
        je = max(rel(l, c.de), rel(la, c.da))
        fe = scalarRel(redot(pb, dm), c.field)
        ae = scalarRel(redot(pb, a), c.alphaPair)
      mv = max(mv, ve)
      mj = max(mj, je)
      mf = max(mf, fe)
      ma = max(ma, ae)
      # Exact norm thresholds can round to either side under a different sum order.
      # Adjacent threshold fixtures have a 1e-8 relative margin and check the branch.
      if c.checkBranch:
        check ns == c.branch
      check rel(e, c.e) < 2e-13
      check ve < 2e-13
      check je < 2e-13
      check fe < 2e-10
      check ae < 2e-10
    echo "selected Poly12 max normalized value=", mv, " jet=", mj,
      " scaled Lie field pullback=", mf, " alpha pullback=", ma

  test "U1 retains its exact scalar specialization with a scale argument":
    type U = MatrixArray[1,1,ComplexType[float64]]
    var m, g, gc, p, c: U
    m[0,0] = newComplex(0.2, 8.0)
    c[0,0] = newComplex(-0.4, 0.7)
    for order in [1, 3, 7, 11, 13]:
      check abs(expProjMulLogJac(m, order=order, scale=5)-ln(1.2)) < 1e-15
      g.expProjMulLogJacGrad(m, order=order, scale=5)
      gc.expProjMulLogJacGrad(p, m, c, order=order, scale=5)
      check abs(g[0,0].re-1.0/1.2) < 1e-15
      check g[0,0].im == 0.0
      check norm2(g-gc) == 0.0
      check p[0,0].re == 0.0
      check p[0,0].im == 0.7

when declared(SimdD4):
  suite "scaled differential with mixed SIMD thresholds":
    test "all lanes use the branch selected by the largest primal norm":
      type S = MatrixArray[3,3,ComplexType[SimdD4]]
      var f, dm, cot, e, pb: S
      for i in 0..<3:
        for j in 0..<3:
          var fr, fi, dr, di, cr, ci: array[4,float64]
          for k in 0..<4:
            let c = refs.simd[k]
            fr[k] = c.f[i,j].re
            fi[k] = c.f[i,j].im
            dr[k] = c.dm[i,j].re
            di[k] = c.dm[i,j].im
            cr[k] = c.c[i,j].re
            ci[k] = c.c[i,j].im
          f[i,j].re := fr
          f[i,j].im := fi
          dm[i,j].re := dr
          dm[i,j].im := di
          cot[i,j].re := cr
          cot[i,j].im := ci
      e := expAH(f)
      pb.expProjectTAHPullback(f, e.adj*cot, order=13, scale=5)
      let got = redot(pb, dm)
      for k in 0..<4:
        var ek: M
        for i in 0..<3:
          for j in 0..<3:
            ek[i,j].re = e[i,j].re[k]
            ek[i,j].im = e[i,j].im[k]
        check rel(ek, refs.simd[k].e) < 2e-13
        check scalarRel(got[k], refs.simd[k].field) < 2e-10
