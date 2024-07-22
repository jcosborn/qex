import math, unittest

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))

import core, scalar

template checkeq(ii: tuple[filename:string, line:int, column:int], sa: string, a: float, sb: string, b: float) =
  if not almostEqual(a, b, unitsInLastPlace = 64):
    checkpoint(ii.filename & ":" & $ii.line & ":" & $ii.column & ": Check failed: " & sa & " :~ " & sb)
    checkpoint("  " & sa & ": " & $a)
    checkpoint("  " & sb & ": " & $b)
    fail()

template checkeq(ii: tuple[filename:string, line:int, column:int], sa: string, a: int, sb: string, b: int) =
  if a != b:
    checkpoint(ii.filename & ":" & $ii.line & ":" & $ii.column & ": Check failed: " & sa & " :~ " & sb)
    checkpoint("  " & sa & ": " & $a)
    checkpoint("  " & sb & ": " & $b)
    fail()

template `:~`(a:Gvalue, b:float) =
  checkeq(instantiationInfo(), astToStr a, a.eval.getfloat, astToStr b, b)

template `:~`(a:Gvalue, b:int) =
  checkeq(instantiationInfo(), astToStr a, a.eval.getint, astToStr b, b)

suite "scalar basic":
  # run once before
  setup:
    # before each test
    let a = 0.5 * (sqrt(5.0) - 1.0)
    let b = sqrt(2.0) - 1.0
    let x = toGvalue(a)
    let y = toGvalue(b)
  #teardown:
    # after each test
  # run once after

  test "assign":
    x :~ a
    y :~ b

  test "n":
    let z = -x
    let dx = z.grad x
    z :~ -a
    dx :~ -1.0

  test "a":
    let z = x+y
    let dx = z.grad x
    let dy = z.grad y
    z :~ a+b
    dx :~ 1.0
    dy :~ 1.0

  test "m":
    let z = x*y
    let dx = z.grad x
    let dy = z.grad y
    z :~ a*b
    dx :~ b
    dy :~ a

  test "s":
    let z = x-y
    let dx = z.grad x
    let dy = z.grad y
    z :~ a-b
    dx :~ 1.0
    dy :~ -1.0

  test "d":
    let z = x/y
    let dx = z.grad x
    let dy = z.grad y
    z :~ a/b
    dx :~ 1.0/b
    dy :~ -a/(b*b)

  test "exp":
    let z = exp(x)
    let dx = z.grad x
    let ddx = dx.grad x
    let dddx = ddx.grad x
    let e = exp(a)
    z :~ e
    dx :~ e
    ddx :~ e
    dddx :~ e

  test "nm":
    let z = (-x)*x
    let dx = z.grad x
    z :~ -a*a
    dx :~ -2.0*a

  test "nm exp":
    let z = (-exp(x))*exp(x)
    let dx = z.grad x
    z :~ -exp(2.0*a)
    dx :~ -2.0*exp(2.0*a)

  test "am":
    let z = (x+y)*x
    let dx = z.grad x
    let dy = z.grad y
    z :~ (a+b)*a
    dx :~ 2.0*a+b
    dy :~ a

  test "ama":
    let w = x
    let v = w+y
    let z = v*v
    let dy = z.grad y
    z :~ (a+b)*(a+b)
    dy :~ 2.0*(a+b)

  test "amd":
    let w = x
    let v = w+y
    let z = v*v/w
    let dy = z.grad y
    z :~ (a+b)*(a+b)/a
    dy :~ 2.0*(a+b)/a

  test "amnd":
    let w = x
    let v = w+y
    let z = v*(-v)/w
    let dy = z.grad y
    z :~ (a+b)*(-a-b)/a
    dy :~ -2.0*(a+b)/a

  test "samnd":
    let w = x-2.0
    let v = w+y
    let z = v*(-v)/w
    let dy = z.grad y
    z :~ (a+b-2.0)*(2.0-a-b)/(a-2.0)
    dy :~ -2.0*(a+b-2.0)/(a-2.0)

suite "scalar d2":
  setup:
    let a = 0.5 * (sqrt(5.0) - 1.0)
    let b = sqrt(2.0) - 1.0
    let c = 2.0 * a - 1.0
    let d = a + 3.0 * b - 1.0
    let x = Gscalar()
    let y = Gscalar()
    x.update a
    y.update b

  test "samnd dx dy":
    let w = x-2.0
    let v = w+y
    let z = v*(-v)/w
    let dy = z.grad y
    let dxy = dy.grad x
    z :~ (a+b-2.0)*(2.0-a-b)/(a-2.0)
    dy :~ -2.0*(a+b-2.0)/(a-2.0)
    dxy :~ 2.0*b/((a-2.0)*(a-2.0))

  test "samnd dx dy repeat":
    let w = x-2.0
    let v = w+y
    let z = v*(-v)/w
    let dy = z.grad y
    let dxy = dy.grad x
    z :~ (a+b-2.0)*(2.0-a-b)/(a-2.0)
    dy :~ -2.0*(a+b-2.0)/(a-2.0)
    dxy :~ 2.0*b/((a-2.0)*(a-2.0))
    y.update c
    dy :~ -2.0*(a+c-2.0)/(a-2.0)
    x.update d
    dxy :~ 2.0*c/((d-2.0)*(d-2.0))
    y.update a
    z :~ (d+a-2.0)*(2.0-d-a)/(d-2.0)
    dy :~ -2.0*(d+a-2.0)/(d-2.0)
    dxy :~ 2.0*a/((d-2.0)*(d-2.0))

  test "samndpdy dx":
    let w = x-2.0
    let v = w+y
    let z = v*(-v)/w
    let dy = z.grad y
    let u = z+0.1*dy
    let dx = (u*u).grad x
    z :~ (a+b-2.0)*(2.0-a-b)/(a-2.0)
    dx :~ -2.0*(b+a-2.0)*(5.0*b+5.0*a-9.0)*(5.0*b*b+b-5.0*a*a+20.0*a-20.0)/(25.0*(a-2.0)*(a-2.0)*(a-2.0))
    y.update c
    dx :~ -2.0*(c+a-2.0)*(5.0*c+5.0*a-9.0)*(5.0*c*c+c-5.0*a*a+20.0*a-20.0)/(25.0*(a-2.0)*(a-2.0)*(a-2.0))
    x.update d
    y.update a
    u :~ (d+a-2.0)*(2.0-d-a)/(d-2.0) - 0.1*2.0*(d+a-2.0)/(d-2.0)
    dx :~ -2.0*(a+d-2.0)*(5.0*a+5.0*d-9.0)*(5.0*a*a+a-5.0*d*d+20.0*d-20.0)/(25.0*(d-2.0)*(d-2.0)*(d-2.0))

suite "bool and cond":
  setup:
    let a = 0.5 * (sqrt(5.0) - 1.0)
    let b = sqrt(2.0) - 1.0
    let c = 2.0 * a - 1.0
    let d = a + 3.0 * b - 1.0
    let x = toGvalue a
    let y = toGvalue b

  test "not":
    let f = toGvalue 0
    not(f) :~ 1
    not(not f) :~ 0
    let t = toGvalue 1.0
    not(t) :~ 0.0
    not(not t) :~ 1.0

  test "and":
    let fi = toGvalue 0
    let ti = toGvalue 1
    let t = toGvalue 1.0
    let f = toGvalue 0.0
    fi and t :~ 0.0
    t and fi :~ 0
    ti and t :~ 1.0
    t and ti :~ 1
    f and fi :~ 0
    fi and f :~ 0.0

  test "or":
    let fi = toGvalue 0
    let ti = toGvalue 1
    let t = toGvalue 1.0
    let f = toGvalue 0.0
    fi or t :~ 1.0
    t or fi :~ 1
    ti or t :~ 1.0
    t or ti :~ 1
    f or fi :~ 0
    fi or f :~ 0.0

  test "xor":
    let fi = toGvalue 0
    let ti = toGvalue 1
    let t = toGvalue 1.0
    let f = toGvalue 0.0
    fi xor t :~ 1.0
    t xor fi :~ 1
    ti xor t :~ 0.0
    t xor ti :~ 0
    f xor fi :~ 0
    fi xor f :~ 0.0

  test "condi":
    let k = toGvalue 0
    let z = cond(k, x, y)
    let dx = z.grad x
    let dy = z.grad y
    z :~ b
    dx :~ 0.0
    dy :~ 1.0
    k.update 1
    z :~ a
    dx :~ 1.0
    dy :~ 0.0

  test "conds":
    let k = toGvalue 1.0
    let z = cond(k, x, y)
    let dx = z.grad x
    let dy = z.grad y
    z :~ a
    dx :~ 1.0
    dy :~ 0.0
    k.update 0.0
    z :~ b
    dx :~ 0.0
    dy :~ 1.0

  test "condi 2":
    let k = toGvalue 0
    let z = cond(k, x, y)
    let z2 = z*z
    let dx = z2.grad x
    let dy = z2.grad y
    z2 :~ b*b
    dx :~ 0.0
    dy :~ 2.0*b
    k.update 1
    y.update c
    z2 :~ a*a
    dx :~ 2.0*a
    dy :~ 0.0
    k.update 0
    z2 :~ c*c
    dx :~ 0.0
    dy :~ 2.0*c

  test "conds 2":
    let k = toGvalue 1.0
    let z = cond(k, x, y)
    let z2 = z*z
    let dx = z2.grad x
    let dy = z2.grad y
    z2 :~ a*a
    dx :~ 2.0*a
    dy :~ 0.0
    k.update 0.0
    x.update d
    z2 :~ b*b
    dx :~ 0.0
    dy :~ 2.0*b
    k.update 1.0
    x.update c
    z2 :~ c*c
    dx :~ 2.0*c
    dy :~ 0.0

  test "cond eval shortcut":
    let t = toGvalue 2.0
    let f = toGvalue 0.0
    let t2 = t*t
    let t3 = t*t*t
    check t2.getfloat == 0.0  # should be zero before eval
    check t3.getfloat == 0.0  # ditto
    var tt = cond(t, t2, t3)
    tt :~ 4.0
    check t2.getfloat == 4.0
    check t3.getfloat == 0.0  # should remain zero after eval
    tt = cond(f, t3, t2)
    tt :~ 4.0
    check t2.getfloat == 4.0
    check t3.getfloat == 0.0  # should remain zero after eval
