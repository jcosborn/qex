import math, unittest

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))

import core, scalar
proc `:~`(a:Gvalue, b:float):bool = a.getfloat.almostEqual b

suite "scalar basic":
  # run once before
  setup:
    # before each test
    let a = 0.5 * (sqrt(5.0) - 1.0)
    let b = sqrt(2.0) - 1.0
    let x = Gscalar()
    let y = Gscalar()
    x := a
    y := b
  #teardown:
    # after each test
  # run once after

  test "assign":
    require x :~ a
    require y :~ b

  test "n":
    let z = -x
    let dx = z.grad x
    z.eval
    dx.eval
    check z :~ -a
    check dx :~ -1.0

  test "a":
    let z = x+y
    let dx = z.grad x
    let dy = z.grad y
    z.eval
    dx.eval
    dy.eval
    check z :~ a+b
    check dx :~ 1.0
    check dy :~ 1.0

  test "m":
    let z = x*y
    let dx = z.grad x
    let dy = z.grad y
    z.eval
    dx.eval
    dy.eval
    check z :~ a*b
    check dx :~ b
    check dy :~ a

  test "s":
    let z = x-y
    let dx = z.grad x
    let dy = z.grad y
    z.eval
    dx.eval
    dy.eval
    check z :~ a-b
    check dx :~ 1.0
    check dy :~ -1.0

  test "d":
    let z = x/y
    let dx = z.grad x
    let dy = z.grad y
    z.eval
    dx.eval
    dy.eval
    check z :~ a/b
    check dx :~ 1.0/b
    check dy :~ -a/(b*b)

  test "nm":
    let z = (-x)*x
    let dx = z.grad x
    z.eval
    dx.eval
    check z :~ -a*a
    check dx :~ -2.0*a

  test "am":
    let z = (x+y)*x
    let dx = z.grad x
    let dy = z.grad y
    z.eval
    dx.eval
    dy.eval
    check z :~ (a+b)*a
    check dx :~ 2.0*a+b
    check dy :~ a

  test "ama":
    let w = x
    let v = w+y
    let z = v*v
    let dy = z.grad y
    z.eval
    dy.eval
    check z :~ (a+b)*(a+b)
    check dy :~ 2.0*(a+b)

  test "amd":
    let w = x
    let v = w+y
    let z = v*v/w
    let dy = z.grad y
    z.eval
    dy.eval
    check z :~ (a+b)*(a+b)/a
    check dy :~ 2.0*(a+b)/a

  test "amnd":
    let w = x
    let v = w+y
    let z = v*(-v)/w
    let dy = z.grad y
    z.eval
    dy.eval
    check z :~ (a+b)*(-a-b)/a
    check dy :~ -2.0*(a+b)/a

  test "samnd":
    let w = x-2.0
    let v = w+y
    let z = v*(-v)/w
    let dy = z.grad y
    z.eval
    dy.eval
    check z :~ (a+b-2.0)*(2.0-a-b)/(a-2.0)
    check dy :~ -2.0*(a+b-2.0)/(a-2.0)

suite "scalar d2":
  setup:
    let a = 0.5 * (sqrt(5.0) - 1.0)
    let b = sqrt(2.0) - 1.0
    let c = 2.0 * a - 1.0
    let d = a + 3.0 * b - 1.0
    let x = Gscalar()
    let y = Gscalar()
    x := a
    y := b

  test "samnd dx dy":
    let w = x-2.0
    let v = w+y
    let z = v*(-v)/w
    let dy = z.grad y
    let dxy = dy.grad x
    z.eval
    dy.eval
    dxy.eval
    check z :~ (a+b-2.0)*(2.0-a-b)/(a-2.0)
    check dy :~ -2.0*(a+b-2.0)/(a-2.0)
    check dxy :~ 2.0*b/((a-2.0)*(a-2.0))

  test "samnd dx dy repeat":
    let w = x-2.0
    let v = w+y
    let z = v*(-v)/w
    let dy = z.grad y
    let dxy = dy.grad x
    z.eval
    dy.eval
    dxy.eval
    check z :~ (a+b-2.0)*(2.0-a-b)/(a-2.0)
    check dy :~ -2.0*(a+b-2.0)/(a-2.0)
    check dxy :~ 2.0*b/((a-2.0)*(a-2.0))
    y := c
    dy.eval
    check dy :~ -2.0*(a+c-2.0)/(a-2.0)
    x := d
    dxy.eval
    check dxy :~ 2.0*c/((d-2.0)*(d-2.0))
    y := a
    z.eval
    dy.eval
    dxy.eval
    check z :~ (d+a-2.0)*(2.0-d-a)/(d-2.0)
    check dy :~ -2.0*(d+a-2.0)/(d-2.0)
    check dxy :~ 2.0*a/((d-2.0)*(d-2.0))

  test "samndpdy dx":
    let w = x-2.0
    let v = w+y
    let z = v*(-v)/w
    let dy = z.grad y
    let u = z+0.1*dy
    let dx = (u*u).grad x
    z.eval
    dx.eval
    check z :~ (a+b-2.0)*(2.0-a-b)/(a-2.0)
    check dx :~ -2.0*(b+a-2.0)*(5.0*b+5.0*a-9.0)*(5.0*b*b+b-5.0*a*a+20.0*a-20.0)/(25.0*(a-2.0)*(a-2.0)*(a-2.0))
    y := c
    dx.eval
    check dx :~ -2.0*(c+a-2.0)*(5.0*c+5.0*a-9.0)*(5.0*c*c+c-5.0*a*a+20.0*a-20.0)/(25.0*(a-2.0)*(a-2.0)*(a-2.0))
    x := d
    y := a
    dx.eval
    u.eval
    check u :~ (d+a-2.0)*(2.0-d-a)/(d-2.0) - 0.1*2.0*(d+a-2.0)/(d-2.0)
    check dx :~ -2.0*(a+d-2.0)*(5.0*a+5.0*d-9.0)*(5.0*a*a+a-5.0*d*d+20.0*d-20.0)/(25.0*(d-2.0)*(d-2.0)*(d-2.0))
