import qex
import maths/groupOps
import algorithms/numdiff
import testutils

type
  M = MatrixArray[3,3,ComplexType[float64]]
  A = MatrixArray[8,8,float64]
  V = VectorArray[8,float64]

var a, d, b: M
block:
  var av, dv, bv: V
  for i in 0..7:
    av[i] = float(i+1)/13.0
    dv[i] = float((i*3) mod 7 - 3)/11.0
    bv[i] = float((i*5) mod 9 - 4)/17.0
  a.suFromVec av
  d.suFromVec dv
  b.suFromVec bv
  a := (1.0 / sqrt(norm2(a))) * a
  d := (1.0 / sqrt(norm2(d))) * d
  b := (1.0 / sqrt(norm2(b))) * b

suite "SU3 exponential approximation contract":
  test "field and scale pullbacks around scaling thresholds":
    # expAH scales at ||A|| = 1/4, then at successive factors of two.
    for amp in [0.0, 0.03, 0.12, 0.249999, 0.250001,
                0.499999, 0.500001, 0.999999, 1.000001, 2.0, 4.0, 8.0]:
      var x, e, p, q, q5: M
      x := amp*a
      e := expAH(x)
      p.projectTAH(expDeriv(x, b))
      q.expProjectTAHPullback(x, e.adj*b)
      q5.expProjectTAHPullback(x, e.adj*b, scale=expProjectTAHScale)
      proc fieldAt(t: float): float = redot(b, expAH(x+t*d))
      proc scaleAt(t: float): float = redot(b, expAH((amp+t)*a))
      var nf, ef, ns, es: float
      ndiff(nf, ef, fieldAt, 0.0, 0.02, ordMax=5)
      ndiff(ns, es, scaleAt, 0.0, 0.02, ordMax=5)
      let pe = abs(redot(p,d)-nf)
      let qe = abs(redot(q,d)-nf)
      let q5e = abs(redot(q5,d)-nf)
      let se = abs(redot(p,a)-ns)
      echo "exp norm=", amp, " poly4 field=", pe,
        " adjoint13 field=", qe, " scaled field=", q5e, " alpha=", se,
        " finite_difference_error=", max(ef,es)
      check pe < 2e-10
      check se < 2e-10
      check q5e < 2e-10
      if amp <= 1.000001:
        check qe < 2e-10

  test "local determinant versus the actual expAH link update":
    var w: M
    w := 1.0
    for amp in [0.0, 0.03, 0.12, 0.25, 0.5, 1.0, 2.0, 4.0]:
      var m, f: M
      var j, j5, jf, df, ad, nj, err: A
      m := 0.1 + amp*a
      j.diffExpProjectTAHMul(jf, df, ad, f, m)
      j5.diffExpProjectTAHMul(jf, df, ad, f, m, scale=expProjectTAHScale)
      proc update(u: M): M =
        var x: M
        x.projectTAH(u*m)
        result := expAH(x)*u
      ndiffSUtoSU(nj, err, update, w, dx=0.02, ordMax=5)
      let rel = abs(determinant(j)-determinant(nj))/max(1.0, abs(determinant(nj)))
      let rel5 = abs(determinant(j5)-determinant(nj))/max(1.0, abs(determinant(nj)))
      echo "Jacobian norm=", amp, " determinant relative error=", rel, " scaled=", rel5
      if amp <= 1.0:
        check rel < 2e-9
      check rel5 < 2e-9
