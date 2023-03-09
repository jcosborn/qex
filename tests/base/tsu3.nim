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
      A0,A1,A2:A

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

    test "log det ∂/∂X [exp(F) X]":
      let eps = 0.12
      var X,Y:M
      X := S0
      Y := S1+S2
      proc f(X:M):M {.noinit.} =
        var r:M
        r := eps * (X * Y.adj)
        r.projectTAH r
        r := -r
        r := exp(r)*X
        r
      var dr,er:A
      ndiffSUtoSU(dr, er, f, X)
      var detj = determinant(dr)
      var j,K,adF,dexpf:A
      var m,F:M
      # combined
      m := eps * (X * Y.adj)
      F.projectTAH(m)
      F := -F
      adF.suad(F)
      let Ms = m + m.adj
      let trMs = trace(Ms).re
      const ii = newComplex(0, -0.5)
      var v:V
      v.suToVec(ii*Ms)
      K.sudabc v
      K += (-1.0/3.0)*trMs
      dexpf.diffExp(adF)
      j = 0.5*(exp(adF) + 1.0 + dexpf * K)
      check(j ~ dr)
      # alt
      var df,ja:A
      dF.diffProjectTAH(-m,F)
      ja = exp(adF) + dexpf * dF
      check(ja ~ dr)
      # simplified detJ
      var ff:M
      var ldj:T
      ldj = smearIndepLogDetJacobian(ff, m)
      check(ff ~ F)
      check(ldj ~ ln(detj))

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
