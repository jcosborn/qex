import testutils
import base
import maths/complexNumbers
import maths/matrixConcept
import maths/matrixFunctions
import maths/types
import algorithms/numdiff
import math

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

when isMainModule:
  runJetTests(1)
  runJetTests(3)
