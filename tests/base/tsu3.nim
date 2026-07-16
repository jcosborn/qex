import qex
import maths/groupOps
import algorithms/numdiff
import testutils

const
  dim = 8
  nc = 3

proc test(T: typedesc) =
  type
    V = VectorArray[dim,T]
    M = MatrixArray[nc,nc,ComplexType[T]]
    A = MatrixArray[dim,dim,T]
  suite("SU(3) group ops :: " & $T):
    var
      x0,x1,x2:V
      X0,X1,X2:M
      S0,S1,S2:M

    block:
      const
        v0:array[dim,float] = [0.043, -0.06, 0.1, 0.16, -0.26, 0.43, -0.68, 1.1]
        v1:array[dim,float] = [-0.574, -0.56, -0.508, -0.442, -0.324, -0.14, 0.162, 0.648]
        v2:array[dim,float] = [0.481, -0.755, 0.009, 0.773, -0.463, 0.301, -0.916, -0.172]
        sl = simdLength(x0[0])

      for i in 0..<nc:
        if sl==1:
          x0[i] := v0[i]
          x1[i] := v1[i]
          x2[i] := v2[i]
        else:
          var t0,t1,t2: array[sl, float]
          for k in 0..<sl:
            t0[k] = v0[i] + 0.3/(k.float-1.5)
            t1[k] = v1[i] - 0.01/(k.float+0.3)
            t2[k] = v2[i] + 0.31/(k.float-2.113)
          x0[i] := t0
          x1[i] := t1
          x2[i] := t2
      X0.suFromVec x0
      X1.suFromVec x1
      X2.suFromVec x2
      S0 := exp(X0)
      S1 := exp(X1)
      S2 := exp(X2)

    test "X = X^a T^a":
      var X: type(X0)
      for a in 0..<8:
        X += x0[a]*su3gen[a]
      check(X~X0)

    test "Tr[T^a T^b] = -½δ^{ab}":
      for a in 0..<8:
        for b in 0..<8:
          let t = trace(su3gen[a]*su3gen[b])
          if a==b:
            check(t.re ~ -0.5)
          else:
            check(t.re ~ 0)
          check(t.im ~ 0)

    test "suFromVec":
      var T0,T1,T2: M
      T0.suFromVec(x0)
      T1 := suFromVec(x0)
      T2 := suFromVec_mat(x0)
      check(T0 ~ T1)
      check(T0 ~ T2)

    test "suToVec":
      var t0,t1,t2:V
      t0.suToVec(X1)
      t1 := suToVec(X1)
      t2 := suToVec_mat(X1)
      check(t0 ~ x1)
      check(t1 ~ x1)
      check(t2 ~ x1)

    test "suTo(From)Vec generator":
      for a in 0..<8:
        let v = su3gen[a].suToVec
        check(su3gen[a] ~ v.suFromVec)

    test "adx y":
      # adx(y) = [x,y]
      var x,y:M
      x.suadApply(suad X0, X1)
      y := X0*X1 - X1*X0
      check(x ~ y)

    test "AdX":
      # exp(adx) = Ad[exp(x)]
      var x,y:A
      x.SUAd(exp(X0))
      y.suad(X0)
      y := exp(y)
      check(x ~ y)

    test "projectTAH vs suTo(From)Vec convention":
      #[
        projectTAH(M)
            = - T^a tr[T^a (M - M†)]
            = 1/2 { δ_il δ_jk (M - M†)_lk - 1/3 δ_ij δ_kl (M - M†)_lk }
            = 1/2 { (M - M†)_ij - 1/3 δ_ij tr(M - M†) }
      ]#
      var ss,T0,T1: M
      var v:V
      ss := S0+S1
      T0.projectTAH ss
      v.suToVec(ss-ss.adj)
      v *= 0.5
      T1.suFromVec(v)
      check(T0 ~ T1)

    test "projectTAH from derivative of SU(3)":
      # projectTAH(M) = - T^a tr[T^a (M - M†)] = T^a ∂_a (- tr[M + M†]) = T^a ∂_a (-2 ReTr M)
      var ss:M
      ss := S0+S1
      proc f(x:M):T = (-2.0)*trace(x).re
      var vd,ve:V
      ndiffSUtoReal(vd, ve, f, ss)
      var p,m:M
      p.suFromVec vd
      m.projectTAH ss
      check(m ~ p)

    test "diffProjectTAH":
      proc f(x:M):M {.noinit.} = result.projectTAH x
      var vd,ve:A
      ndiffSUtoAlg(vd, ve, f, S0)
      var m:M
      var j:A
      m.projectTAH(S0)
      j.diffProjectTAH(S0, m)
      check(vd ~ j)

    test "diffCrossProjectTAH":
      #[
        ∂_Y^b ∂_X^a (-2) ReTr[ X (Z Y)† ] = - ∂_Y^b ∂_X^a tr[ X Y† Z† + Z Y X† ]
            = - 2 ReTr[T^a (- X Y†) T^b Z†]
        Note the extra negative sign from ∂_Y^b.
      ]#
      proc f(x:M):M {.noinit.} =
        result := S0 * adj(S2 * x)
        result.projectTAH(result)
      var vd,ve:A
      ndiffSUtoAlg(vd, ve, f, S1)
      var s,m:M
      var j,dp:A
      s := S0 * adj(S2 * S1)
      m.projectTAH(s)
      dp.diffProjectTAH(s, m)
      j.diffCrossProjectTAH(SUAd(S0 * S1.adj), dp)
      j := -j
      check(vd ~ j)

    test "diff2ProjectTAH":
      #[
        P^a = -tr[T^a (M - M†)]
        ∂_c P^a = -tr[T^a (T^c M + M† T^c)]
                = -1/2 { d^acb tr[T^b i(M+M†)] - 1/3 δ^ac tr(M+M†) + f^acb F^b }
        ∂_d ∂_c P^a = -tr[T^a T^c T^d M - T^c T^a M† T^d]
        Use diffProjTAH, but use T^d M for M.
        The same goes with ∂_d on a different matrix.
        ∂_Y^d ∂_X^c ∂_X^a (-2) ReTr[ X (Z Y)† ] = - 2 ReTr[ T^a T^c X (Z T^d Y)† ]
      ]#
      proc ff(y:M):A {.noinit.} =
        let ss = S0 * adj(S2 * y)
        proc f(x:M):M {.noinit.} = result.projectTAH x
        var vd,ve:A
        ndiffSUtoAlg(vd, ve, f, ss)
        vd
      var dr,er:A
      var z,d:T
      z := 0.0
      d := 2.0
      for a in 0..<8:
        ndiff(dr, er, proc(l:T):A {.noinit.} = ff(exp(l*su3gen[a])*S1), z, d, scale=5.0, ordMax=4)
        var m:M
        var j:A
        let ss = S0 * adj(S2 * su3gen[a] * S1)
        m.projectTAH(ss)
        j.diffProjectTAH(ss, m)
        withCT(1e-12):
          check(dr ~ j)

    test "diffExp(T)":
      proc f(m:M):M = exp(m)
      var dr,er:A
      for a in 0..<8:
        ndiffAlgtoSU(dr, er, f, su3gen[a])
        var j,adx:A
        adx.suad(su3gen[a])
        j.diffExp adx
        withCT(1e-12):
          check(dr ~ j)

    test "diffExp":
      proc f(m:M):M = exp(m)
      let m = X0
      var dr,er:A
      ndiffAlgtoSU(dr, er, f, m)
      var j,adx:A
      adx.suad(m)
      j.diffExp adx
      check(dr ~ j)

    test "diffExp paired recurrence":
      var a, r, rs, p:A
      var v, vr, ve:V
      a.suad(X0)
      v := x1
      for order in 0..14:
        r.diffExp(a, order=order)
        rs := 1.0
        p := 1.0
        var fac = 1.0
        for k in 1..order:
          p := p*a
          fac *= float(k+1)
          rs += (1.0/fac)*p
        withCT(1e-13):
          check(r ~ rs)
        vr.diffExpApply(a, v, order=order)
        ve := r*v
        withCT(1e-13):
          check(vr ~ ve)

    test "projected expDeriv in adjoint representation":
      for s in [0.0, 0.03, 0.12]:
        var m, c, e, q, qp, qr:M
        var a:A
        var v, p:V
        m := s*X0
        c := S1 + 0.37*S2
        e := exp(m)
        q.projectTAH(e.adj*c)
        a.suad(m)
        v.suToVec(q)
        p.diffExpApply(a, v)
        q.suFromVec(p)
        qp.expProjectTAHPullback(m, e.adj*c)
        qr.projectTAH(expDeriv(m, c))
        withCT(1e-12):
          check(q ~ qr)
          check(qp ~ qr)

    test "diffDiffExp":
      var a, da, r, dr, er:A
      a.suad(X0)
      da.suad(X1)
      for h in 1..8:
        r.diffDiffExp(a, da, halfOrder=h)
        proc f(x:T):A =
          result.diffExp(a + x*da, order=2*h-1)
        ndiff(dr, er, f, 0.0, 1.0, scale=5.0, ordMax=4)
        withCT(1e-10):
          check(r ~ dr)

    test "log det ∂_X [exp(ProjectTAH(X Y†)) X]":
      let eps = 0.12
      var X,Y:M
      X := S0
      Y := S1+S2
      proc f(X:M):M {.noinit.} =
        var r:M
        r := eps * (X * Y.adj)
        r.projectTAH r
        r := exp(r)*X
        r
      var dr,er:A
      ndiffSUtoSU(dr, er, f, X)
      let detj = determinant(dr)
      var j,K,adF,dexpf:A
      var m,F:M
      # combined
      #[
        -2 tr[(∂_c Z) Z† T^a]
            = exp(adF)^ac + J(-F)^ab [∂_c F]^b
            = exp(adF)^ac
              - 1/2 [(exp(adF)-1)/adF]^ab {d^bcd tr[T^d i (M + M†)] - 1/N δ^bc tr[M + M†]}
              - 1/2 [(exp(adF)-1)]^ac
            = 1/2 { [(exp(adF)+1)]^ac - [(exp(adF)-1)/adF]^ab {d^bcd tr[T^d i (M + M†)] - 1/N δ^bc tr[M + M†]} }
            = 1/2 { [(exp(adF)+1)]^ac - J(-F)^ab {d^bcd tr[T^d i (M + M†)] - 1/N δ^bc tr[M + M†]} }
      ]#
      m := eps * (X * Y.adj)
      F.projectTAH(m)
      adF.suad(F)
      let Ms = m + m.adj
      let trMs = trace(Ms).re
      const ii = newComplex(0, -0.5)
      var v:V
      v.suToVec(ii*Ms)
      K.sudabc v
      K += (-1.0/3.0)*trMs
      dexpf.diffExp(adF)
      j = 0.5*(exp(adF) + 1.0 - dexpf * K)
      check(j ~ dr)
      # alt
      var df,ja:A
      dF.diffProjectTAH(m,F)
      ja = exp(adF) + dexpf * dF
      check(ja ~ dr)
      # simplified detJ
      var ff:M
      var jj,aF,jF,dd:A
      jj.diffExpProjectTAHMul(jF, dd, aF, ff, m)
      check(ff ~ F)
      check(dd ~ dF)
      check(aF ~ adF)
      dexpf.diffExp(-adF)
      check(jF ~ dexpf)
      check(determinant(jj) ~ detj)
      check(expProjMulLogJac(m) ~ ln(detj))

    test "free-matrix gradient of projected-exponential multiplication log-Jacobian":
      var m0, m1, m2, m3, mz:M
      m0 := 0.04 * (S0 + S1 * S2.adj)
      m1 := 0.03 * (S2 + S0 * S1.adj)
      m2 := 0.12 * (S0 + S1 * S2.adj)
      m3 := 0.60 * (S0 + S1 * S2.adj)
      mz := 0

      template checkFreeGrad(m: untyped, o=13) =
        block:
          var g:M
          g.expProjMulLogJacGrad(m, order=o)
          for i in 0..<nc:
            for j in 0..<nc:
              for k in 0..1:
                var b:M
                b := 0
                if k == 0:
                  b[i, j] := 1.0
                else:
                  b[i, j] := newComplex(0.0, 1.0)
                proc f(x:T):T =
                  expProjMulLogJac(m + x*b, order=o)
                var d, e:T
                ndiff(d, e, f, 0.0, 0.1, ordMax=4)
                check abs(redot(g, b) - d) < 1e-10

      checkFreeGrad(m0)
      checkFreeGrad(m1)
      checkFreeGrad(m2)
      checkFreeGrad(m3)
      checkFreeGrad(mz)
      checkFreeGrad(m0, 1)
      checkFreeGrad(m0, 3)
      checkFreeGrad(m0, 7)
      checkFreeGrad(m0, 11)
      checkFreeGrad(m0, 13)

    test "∂_X log det ∂_X [exp(ProjectTAH(X Y†)) X]":
      let eps = 0.12
      var X,Y:M
      X := S0
      Y := S1+S2
      proc ff(X:M):T =
        proc f(X:M):M {.noinit.} =
          var r:M
          r := eps * (X * Y.adj)
          r.projectTAH r
          r := exp(r)*X
          r
        var dr,er:A
        ndiffSUtoSU(dr, er, f, X, dx=1.0)
        ln(determinant(dr))
      var vd,ve:V
      ndiffSUtoReal(vd, ve, ff, X)
      var r:V
      let m = eps * (X * Y.adj)
      r.diffLnDetDiffExpProjectTAHMul(m)
      withCT 1e-11:
        check(r ~ vd)

    test "∂_Y log det ∂_X [exp(ProjectTAH(X Y)) X]":
      # let eps = 0.12
      let eps = 1.0    # 0.01
      var X,Y:M
      # X := 1.0    # this makes it correct ?!
      # X := S0
      # Y := S1
      # JXY WAS HERE: debug, let's work out explicitly with T⁰
      X := exp(su3gen[0])
      Y := exp(su3gen[0])
      proc ff(Y:M):T =
        proc f(X:M):M {.noinit.} =
          var r:M
          r := eps * (X * Y)
          r.projectTAH r
          r := exp(r)*X
          r
        var dr,er:A
        ndiffSUtoSU(dr, er, f, X, dx=1.0)
        ln(determinant(dr))
      var vd,ve:V
      ndiffSUtoReal(vd, ve, ff, Y)
      var r:V
      let mx = X
      let my = eps * Y
      r.diffCrossLnDetDiffExpProjectTAHMul(mx, my, order=63)
      withCT 1e-9:
        check(r ~ vd)

    test "∂_Y log det ∂_X [exp(ProjectTAH(X Y Z)) X]":
      let eps = 0.12
      var X,Y,Z:M
      X := S0
      Y := S1
      Z := S2
      proc ff(Y:M):T =
        proc f(X:M):M {.noinit.} =
          var r:M
          r := eps * (X * Y * Z)
          r.projectTAH r
          r := exp(r)*X
          r
        var dr,er:A
        ndiffSUtoSU(dr, er, f, X, dx=1.0)
        ln(determinant(dr))
      var vd,ve:V
      ndiffSUtoReal(vd, ve, ff, Y)
      var r:V
      let mx = eps * X
      let my = Y * Z
      r.diffCrossLnDetDiffExpProjectTAHMul(mx, my)
      withCT 1e-11:
        check(r ~ vd)

    test "∂_Y log det ∂_X [exp(ProjectTAH(X (Z Y)†)) X]":
      let eps = 0.12
      var X,Y,Z:M
      X := S0
      Y := S1
      Z := S2
      proc ff(Y:M):T =
        proc f(X:M):M {.noinit.} =
          var r:M
          r := eps * (X * (Z * Y).adj)
          r.projectTAH r
          r := exp(r)*X
          r
        var dr,er:A
        ndiffSUtoSU(dr, er, f, X, dx=1.0)
        ln(determinant(dr))
      var vd,ve:V
      ndiffSUtoReal(vd, ve, ff, Y)
      var r:V
      let mx = eps * (X * Y.adj)
      let my = Z.adj
      r.diffCrossAdjLnDetDiffExpProjectTAHMul(mx, my)
      withCT 1e-11:
        check(r ~ vd)

    test "∂_Y log det ∂_X [exp(ProjectTAH((X Y) Z)) X] (cross on middle)":
      let eps = 0.23
      var X,Y,Z:M
      X := S0
      Y := S1
      Z := S2
      proc ff(Y:M):T =
        proc f(X:M):M {.noinit.} =
          var r:M
          r := eps * ((X * Y) * Z)
          r.projectTAH r
          r := exp(r)*X
          r
        var dr,er:A
        ndiffSUtoSU(dr, er, f, X, dx=1.0)
        ln(determinant(dr))
      var vd,ve:V
      ndiffSUtoReal(vd, ve, ff, Y)
      var r:V
      let mx = eps * X
      let my = Y * Z
      r.diffCrossLnDetDiffExpProjectTAHMul(mx, my, order=63)
      withCT 1e-11:
        check(r ~ vd)

    test "∂_Z log det ∂_X [exp(ProjectTAH(X (Y Z))) X] (cross on right)":
      let eps = 0.07
      var X,Y,Z:M
      X := S0
      Y := S1
      Z := S2
      proc ff(Z:M):T =
        proc f(X:M):M {.noinit.} =
          var r:M
          r := eps * (X * (Y * Z))
          r.projectTAH r
          r := exp(r)*X
          r
        var dr,er:A
        ndiffSUtoSU(dr, er, f, X, dx=1.0)
        ln(determinant(dr))
      var vd,ve:V
      ndiffSUtoReal(vd, ve, ff, Z)
      var r:V
      # Factor M = X (Y Z) as (X Y) · Z so the cross routine applies to the right factor
      let mx = eps * (X * Y)
      let my = Z
      r.diffCrossLnDetDiffExpProjectTAHMul(mx, my, order=63)
      withCT 1e-11:
        check(r ~ vd)

template doTest(t:untyped) =
  when declared(t):
    test(t)
doTest(float64)
# doTest(SimdD1)
# doTest(SimdD2)
# doTest(SimdD4)
# doTest(SimdD8)
# doTest(SimdD16)
#[
doTest(float32)
doTest(SimdS1)
doTest(SimdS2)
doTest(SimdS4)
doTest(SimdS8)
doTest(SimdS16)
]#

proc testDiffExpSu3Ad(T: typedesc) =
  type
    V = VectorArray[dim, T]
    M = MatrixArray[nc, nc, ComplexType[T]]
    A = MatrixArray[dim, dim, T]
  suite("SU(3) adjoint Cayley-Hamilton :: " & $T):
    test "matches the dense degree-13 series":
      const av: array[dim, float] = [0.043, -0.06, 0.1, 0.16, -0.26, 0.43, -0.68, 1.1]
      var a: V
      for i in 0..<dim:
        when T is SomeFloat:
          a[i] := av[i]
        else:
          var aa: array[simdLength(T), float]
          for k in 0..<aa.len:
            aa[k] = av[i] + 0.017*float(k-i)
          a[i] := aa
      var m0: M
      m0.suFromVec(a)
      for scale in [0.0, 0.1, 0.7, 1.5]:
        var m, f: M
        var j, jf, df, ad, expected: A
        m := scale*m0
        j.diffExpProjectTAHMul(jf, df, ad, f, m)
        expected.diffExp(-ad)
        let
          d2 = simdSum(norm2(jf - expected))
          e2 = simdSum(norm2(expected))
        check d2 <= 2e-25*max(1.0, e2)

template doCayleyHamiltonTest(t: untyped) =
  when declared(t):
    testDiffExpSu3Ad(t)
doCayleyHamiltonTest(float64)
doCayleyHamiltonTest(SimdD4)

proc testDiffExpSuApply(T: typedesc) =
  type
    V = VectorArray[dim, T]
    M = MatrixArray[nc, nc, ComplexType[T]]
    A = MatrixArray[dim, dim, T]
  suite("SU(3) adjoint apply :: " & $T):
    test "matches dense adjoint series":
      const
        av: array[dim, float] = [0.043, -0.06, 0.1, 0.16, -0.26, 0.43, -0.68, 1.1]
        xv: array[dim, float] = [-0.574, -0.56, -0.508, -0.442, -0.324, -0.14, 0.162, 0.648]
        sl = simdLength(T)
      var a, x: V
      for i in 0..<dim:
        when sl == 1:
          a[i] := av[i]
          x[i] := xv[i]
        else:
          var aa, xx: array[sl, float]
          for k in 0..<sl:
            aa[k] = av[i] + 0.017*float(k-i)
            xx[k] = xv[i] - 0.013*float(k+i)
          a[i] := aa
          x[i] := xx
      var m: M
      var ad: A
      m.suFromVec(a)
      ad.suad(m)
      for order in 0..14:
        var got, expected: V
        got.diffExpSuApply(m, x, order)
        expected.diffExpApply(ad, x, order)
        let
          d2 = simdSum(norm2(got - expected))
          e2 = simdSum(norm2(expected))
        check(d2 <= 4e-26*e2)
        if order == 0:
          check simdSum(norm2(got - x)) == 0
      var q: T
      when sl == 1:
        q := 0.4
      else:
        var qq: array[sl, float]
        for k in 0..<sl:
          qq[k] = 0.2 + 0.1*float(k)
        q := qq
      m := 0
      m[0, 0].im := q
      m[1, 1].im := q
      m[2, 2].im := -2.0*q
      ad.suad(m)
      var got, expected: V
      got.diffExpSuApply(m, x)
      expected.diffExpApply(ad, x)
      let
        d2 = simdSum(norm2(got - expected))
        e2 = simdSum(norm2(expected))
      check d2 <= 4e-26*e2

    test "logdet gradient shares its adjoint application":
      const
        av: array[dim, float] = [0.043, -0.06, 0.1, 0.16, -0.26, 0.43, -0.68, 1.1]
        xv: array[dim, float] = [-0.574, -0.56, -0.508, -0.442, -0.324, -0.14, 0.162, 0.648]
        sl = simdLength(T)
      var a, x: V
      for i in 0..<dim:
        when sl == 1:
          a[i] := av[i]
          x[i] := xv[i]
        else:
          var aa, xx: array[sl, float]
          for k in 0..<sl:
            aa[k] = av[i] + 0.017*float(k-i)
            xx[k] = xv[i] - 0.013*float(k+i)
          a[i] := aa
          x[i] := xx
      var m, f, g, gm, pmat, xmat, expectedG, expectedP: M
      var got, expected: V
      m.suFromVec(a)
      m *= 0.05
      f.projectTAH(m)
      expectedG.expProjMulLogJacGrad(m)
      g.expProjMulLogJacGrad(got, m, x)
      expected.diffExpSuApply(f, x)
      xmat.suFromVec(x)
      gm.expProjMulLogJacGrad(pmat, m, xmat)
      expectedP.suFromVec(expected)
      let
        dg2 = simdSum(norm2(g - expectedG))
        g2 = simdSum(norm2(expectedG))
        dp2 = simdSum(norm2(got - expected))
        p2 = simdSum(norm2(expected))
        dgm2 = simdSum(norm2(gm - expectedG))
        dpm2 = simdSum(norm2(pmat - expectedP))
      check dg2 <= 4e-26*max(g2, 1.0)
      check dp2 <= 4e-26*p2
      check dgm2 <= 4e-26*max(g2, 1.0)
      check dpm2 <= 4e-26*p2

    test "contracted reverse logdet gradient matches directional series":
      const sl = simdLength(T)
      var m0: M
      for i in 0..<nc:
        for j in 0..<nc:
          when sl == 1:
            m0[i, j].re := 0.11*float(1 + 2*i - j)
            m0[i, j].im := -0.07*float(2 - i + 3*j)
          else:
            var rr, ii: array[sl, float]
            for k in 0..<sl:
              rr[k] = 0.11*float(1 + 2*i - j) + 0.013*float(k)
              ii[k] = -0.07*float(2 - i + 3*j) - 0.009*float(k)
            m0[i, j].re := rr
            m0[i, j].im := ii
      for s in [0.03, 0.08, 0.14, 0.60, 1.00]:
        let m = s*m0
        var g: M
        var got, expected: V
        g.expProjMulLogJacGrad(m)
        expected.diffLnDetDiffExpProjectTAHMul(m)
        for d in 0..<dim:
          got[d] := redot(g, su3gen[d]*m)
        let
          d2 = simdSum(norm2(got - expected))
          e2 = simdSum(norm2(expected))
        check d2 <= 2e-18*max(1.0, e2)

template doApplyTest(t: untyped) =
  when declared(t):
    testDiffExpSuApply(t)
doApplyTest(float64)
doApplyTest(SimdD4)
