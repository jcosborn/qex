suite "scalar d2":
  setup:
    let fixture = initMutableScalarPair(grt)
    let a {.used.} = fixture.a
    let b {.used.} = fixture.b
    let c {.used.} = fixture.c
    let d {.used.} = fixture.d
    let x {.used.} = fixture.x
    let y {.used.} = fixture.y

  test "zero upstream scalar keeps higher order gradient live":
    let s = grt.toGvalue(0.0)
    let z = s * (x * x)
    let dzdx = z.grad x
    let d2zdxds = dzdx.grad s

    dzdx :~ 0.0
    d2zdxds :~ 2.0 * a

    x.update d
    d2zdxds :~ 2.0 * d

    s.update 5.0
    dzdx :~ 10.0 * d
    d2zdxds :~ 2.0 * d

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
