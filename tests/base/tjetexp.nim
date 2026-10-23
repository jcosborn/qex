import testutils
import base
import maths/complexNumbers
import maths/matrixConcept
import maths/matrixFunctions
import maths/types
from maths/matexp import expAH, ExpAHWork, newExpAHWork, bytes
import algorithms/numdiff
import math
import simd
import base/alignedMem

# expTop(y; d_1..d_M) = d/de_1 ... d/de_M exp(y + sum e_i d_i) for the exp kernel.
#   expTop(x^dag; c)        == expDeriv(x, c)
#   expTop(x^dag; u^dag, b) == d/dt expDeriv(x + t u, b)      (the order-two backward)
#   expTop(y; d, w)         == d/dt expTop(y + t w; d)         (one more direction)

proc fillGen(m: var auto, s: float) =
  const nc = m.nrows
  for i in 0..<nc:
    for j in 0..<nc:
      m[i,j].re := 0.3*sin(s + 1.1*i.float - 0.7*j.float)
      m[i,j].im := 0.3*cos(s + 0.6*i.float + 1.2*j.float)

proc runJetTests(nc: static int) =
  type M = MatrixArray[nc, nc, ComplexType[float]]
  var x, u, b, w, e, xa, ua, wa: M
  fillGen(x, 4.0)
  fillGen(u, 2.0)
  fillGen(b, 3.0)
  fillGen(w, 5.0)
  fillGen(e, 1.0)
  xa := x.adj
  ua := u.adj
  wa := w.adj
  proc nd(f: proc(t: float): float): float =
    var err: float
    ndiff(result, err, f, 0.0, 0.1)
  suite "expTop Nc=" & $nc:
    test "one direction is expDeriv":
      withCT(1e-12): check redot(e, expTop(xa, [b])) ~ redot(e, expDeriv(x, b))
    test "two directions is the derivative of expDeriv":
      let fd = nd(proc(t: float): float = redot(e, expDeriv(x + t*u, b)))
      withCT(1e-8): check fd ~ redot(e, expTop(xa, [ua, b]))
    test "two directions are symmetric":
      withCT(1e-12): check redot(e, expTop(xa, [ua, b])) ~ redot(e, expTop(xa, [b, ua]))
    test "three directions is the derivative of two":
      let fd = nd(proc(t: float): float = redot(e, expTop(xa + t*wa, [ua, b])))
      withCT(1e-8): check fd ~ redot(e, expTop(xa, [ua, b, wa]))

  var ah: M
  ah.projectTAH(x)
  ah := (1.0/sqrt(norm2(ah)))*ah
  suite "expAH jets Nc=" & $nc:
    test "directional derivatives on both sides of scaling thresholds":
      for amp in [0.03, 0.249999, 0.250001, 0.499999, 0.500001, 2.0, 8.0]:
        var y, ya: M
        y := amp*ah
        ya := y.adj
        let fd = nd(proc(t: float): float = redot(e, expAH(y+t*b)))
        let j1 = redot(e, expAH(y, [b]))
        withCT(1e-8): check fd ~ j1
        withCT(1e-12): check j1 ~ redot(expAH(ya, [e]), b)
        let fd2 = nd(proc(t: float): float = redot(e, expAH(y+t*u, [b])))
        withCT(1e-8): check fd2 ~ redot(e, expAH(y, [u,b]))
        let fd3 = nd(proc(t: float): float = redot(e, expAH(y+t*w, [u,b])))
        withCT(1e-8): check fd3 ~ redot(e, expAH(y, [w,u,b]))
        withCT(1e-12): check redot(e, expAH(y, [u,b,w])) ~ redot(e, expAH(y, [w,b,u]))

    test "runtime scratch agrees with static jets and can change order":
      var work: ExpAHWork[M]
      var y: M
      y := 0.8*ah
      for ds in [@[u,b,w,e], @[b], @[u,b], @[u,b,w,e,u]]:
        var expected: M
        case ds.len
        of 1: expected := expAH(y, [b])
        of 2: expected := expAH(y, [u,b])
        of 4: expected := expAH(y, [u,b,w,e])
        else: expected := expAH(y, [u,b,w,e,u])
        check norm2(expAH(y, ds, work)-expected) < 1e-26
      let fd4 = nd(proc(t: float): float = redot(e, expAH(y+t*w, [u,b,e])))
      withCT(1e-8): check fd4 ~ redot(e, expAH(y, [u,b,e,w]))

    test "presized scratch matches lazily sized scratch exactly":
      var y: M
      y := 0.8*ah
      var sized = newExpAHWork[M](4)
      var lazy: ExpAHWork[M]
      let capacity = sized.bytes
      check capacity <= 6*(16*sizeof(M)+64)
      check norm2(expAH(y, @[u,b,w,e], sized) - expAH(y, @[u,b,w,e], lazy)) == 0
      check norm2(expAH(y, @[u,b,w,e], sized) - expAH(y, [u,b,w,e])) == 0
      check sized.bytes == capacity

    test "runtime scratch reuses capacity when order changes":
      var y: M
      y := 0.8*ah
      var work = newExpAHWork[M](4)
      let bytes = work.bytes
      let allocated = getRawMemAllocated()
      for ds in [@[u,b,w,e], @[b], @[u,b], @[u,b,w,e]]:
        var expected: M
        case ds.len
        of 1: expected := expAH(y, [b])
        of 2: expected := expAH(y, [u,b])
        else: expected := expAH(y, [u,b,w,e])
        check norm2(expAH(y, ds, work)-expected) < 1e-26
        check work.bytes == bytes
        check getRawMemAllocated() == allocated

    test "the first jet also returns the exponential":
      for amp in [0.03, 0.249999, 0.250001, 0.8, 8.0]:
        var y, ev, dv: M
        y := amp*ah
        dv := expAH(y, b, ev)
        # The same jet; compilers may contract the two procs differently.
        check norm2(dv-expAH(y, [b])) <= (4*epsilon(float))^2*(1+norm2(dv))
        # Plain squaring against expm1 squaring: the same E to rounding.
        check norm2(ev-expAH(y)) < 1e-24

    test "large diagonal anti-Hermitian argument has an accurate pullback":
      for t in [10000.0, 100000.0]:
        var y, d, expected: M
        y := 0
        d := 0
        expected := 0
        y[0,0].im := t
        d[0,0].im := 1.0
        expected[0,0].re := -sin(t)
        expected[0,0].im := cos(t)
        when nc > 1:
          y[1,1].im := -t
          d[1,1].im := -1.0
          expected[1,1].re := -sin(t)
          expected[1,1].im := -cos(t)
        check norm2(expAH(y, [d])-expected) < 1e-16

proc pairedViews(nc: static int, S: typedesc) =
  type M = MatrixArray[nc,nc,ComplexType[S]]
  var x, d, ah: M
  when S is SomeFloat:
    fillGen(x, 4.0)
    fillGen(d, 3.0)
  else:
    for k in 0..<simdLength(S):
      var a, b: MatrixArray[nc,nc,ComplexType[float64]]
      fillGen(a, 4.0+0.1*k.float)
      fillGen(b, 3.0+0.2*k.float)
      x[asSimd(k)] := a
      d[asSimd(k)] := b
  ah.projectTAH(x)
  ah := (1.0/sqrt(norm2(ah).simdMax))*ah
  suite "paired expAH views Nc=" & $nc & " " & $S:
    test "mixed views and owning matrices agree for both results":
      for amp in [0.0, 0.249999, 0.250001, 0.8, 8.0]:
        var y, ya, da: M
        y := amp*ah
        ya := y.adj
        da := d.adj
        template checkPair(a, b, ar, br: untyped) =
          block:
            var e, v: M
            v := expAH(a, b, e)
            let jr = expAH(ar, [br])
            let er = expAH(ar)
            check norm2(v-jr).simdMax <= (4*epsilon(float))^2*(1+norm2(jr).simdMax)
            check norm2(e-er).simdMax < 1e-24
        checkPair(y.adj, d, ya, d)
        checkPair(y, d.adj, y, da)
        checkPair(y.adj, d.adj, ya, da)

when isMainModule:
  runJetTests(1)
  runJetTests(3)
  pairedViews(1, float64)
  pairedViews(3, float64)
  template simdViews(n: static int) =
    pairedViews(1, `SimdD n`)
    pairedViews(3, `SimdD n`)
  simdViews(VLEN)
  when 2*VLEN <= 16:
    simdViews(2*VLEN)
