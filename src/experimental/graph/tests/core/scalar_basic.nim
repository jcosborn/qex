suite "scalar basic":
  setup:
    let fixture = initScalarLeafPair(grt)
    let a {.used.} = fixture.a
    let b {.used.} = fixture.b
    let x {.used.} = fixture.x
    let y {.used.} = fixture.y

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

  test "shared intermediate sees distinct upstream gradients across outputs":
    let w = x * y
    let dwsqdx = (w * w).grad x
    let dexpwdx = exp(w).grad x
    let dwdx = w.grad x

    dexpwdx :~ exp(a * b) * b
    dwdx :~ b
    dwsqdx :~ 2.0 * a * b * b

    y.update 0.25
    dwdx :~ 0.25
    dwsqdx :~ 2.0 * a * 0.25 * 0.25
    dexpwdx :~ exp(a * 0.25) * 0.25

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

  test "division by zero fails fast":
    let x2 = grt.toGvalue(2.0)
    let zero = grt.toGvalue(0.0)
    let z = x2 / zero
    expect(GraphValueError):
      discard z.eval

    let dzdx = z.grad x2
    expect(GraphValueError):
      discard dzdx.eval

  test "scalar getters reject wrong value type":
    expect(GraphValueError):
      discard grt.toGvalue(1).getfloat
    expect(GraphValueError):
      discard grt.toGvalue(1.0).getint

  test "separate runtimes isolate node ids and reject mixed graphs":
    let leftGrt = initGraphRuntime()
    let rightGrt = initGraphRuntime()

    let xLeft = scalarNodeIn(leftGrt)
    let yLeft = scalarNodeIn(leftGrt)
    xLeft.update 2.0
    yLeft.update 3.0

    let xRight = scalarNodeIn(rightGrt)
    xRight.update 5.0

    let zLeft = xLeft + yLeft
    zLeft :~ 5.0
    check xLeft.stableNodeId == 1
    check yLeft.stableNodeId == 2
    check xRight.stableNodeId == 1

    expect(GraphValueError):
      discard xLeft + xRight

  test "findGrad rejects runtime-local node id collisions":
    let leftGrt = initGraphRuntime()
    let rightGrt = initGraphRuntime()

    let xLeft = leftGrt.toGvalue(2.0)
    let yLeft = leftGrt.toGvalue(3.0)
    let xRight = rightGrt.toGvalue(5.0)
    let zLeft = xLeft + yLeft

    discard zLeft.grad(xLeft).eval
    check xLeft.stableNodeId == xRight.stableNodeId

    expect(GraphValueError):
      discard findGrad(xRight, zLeft)

  test "mixed numeric literals stay inside the anchored runtime":
    let grt = initGraphRuntime()
    let xRuntime = grt.toGvalue(2.0)
    let kRuntime = grt.toGvalue(1)
    let sum = xRuntime + 3.0
    let scaled = 0.5 * sum
    let chosen: Gscalar = cond(kRuntime, scaled, 0.0)
    let dx = chosen.grad xRuntime

    chosen :~ 2.5
    dx :~ 0.5
    check sum.runtime == grt
    check scaled.runtime == grt
    check chosen.runtime == grt

  test "literal and erased-right scalar overloads preserve concrete type":
    let erasedY: Gvalue = y

    let addRight: Gscalar = x + 2.0
    let addLeft: Gscalar = 2.0 + x
    let subRight: Gscalar = x - 2
    let subLeft: Gscalar = 2 - x
    let mulRight: Gscalar = x * 2.0
    let mulLeft: Gscalar = 2.0 * x
    let divRight: Gscalar = x / 2
    let divLeft: Gscalar = 2 / x

    addRight :~ a + 2.0
    addLeft :~ 2.0 + a
    subRight :~ a - 2.0
    subLeft :~ 2.0 - a
    mulRight :~ 2.0 * a
    mulLeft :~ 2.0 * a
    divRight :~ a / 2.0
    divLeft :~ 2.0 / a

    let addErased: Gscalar = x + erasedY
    let subErased: Gscalar = x - erasedY
    let mulErased: Gscalar = x * erasedY
    let divErased: Gscalar = x / erasedY
    let chained: Gscalar = (x + erasedY) * 2.0

    addErased :~ a + b
    subErased :~ a - b
    mulErased :~ a * b
    divErased :~ a / b
    chained :~ 2.0 * (a + b)
    chained.grad(x) :~ 2.0

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

  test "backwardTarget takes precedence over backward":
    proc targetPriorityf(v: Gvalue) =
      v.valCopy v.inputs[0]

    proc targetPriorityb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
      discard z
      discard dep
      if i != 0:
        raiseValueError("targetPriority backward expects input 0, got " & $i)
      let plain = grt.toGvalue(7.0)
      if zb == nil:
        return plain
      plain.scaleLike zb

    proc targetPriorityBackwardTarget(zb: Gvalue,
                                      z: Gvalue,
                                      target: Gvalue,
                                      dep: Gvalue): Gvalue =
      discard z
      discard target
      discard dep
      let targeted = grt.toGvalue(11.0)
      if zb == nil:
        return targeted
      targeted.scaleLike zb

    let gtargetPriority = newGfunc(
      forward = targetPriorityf,
      backward = targetPriorityb,
      backwardTarget = targetPriorityBackwardTarget,
      name = "targetPriority")

    let z = graphNode(scalarNodeLike(x), @[x], gtargetPriority)
    let dx = z.grad x
    z :~ a
    dx :~ 11.0

  test "backwardTarget works without backward hook":
    proc targetOnlyf(v: Gvalue) =
      v.valCopy v.inputs[0]

    proc targetOnlyBackwardTarget(zb: Gvalue,
                                  z: Gvalue,
                                  target: Gvalue,
                                  dep: Gvalue): Gvalue =
      discard z
      discard target
      discard dep
      let targeted = grt.toGvalue(13.0)
      if zb == nil:
        return targeted
      targeted.scaleLike zb

    let gtargetOnly = newGfunc(
      forward = targetOnlyf,
      backwardTarget = targetOnlyBackwardTarget,
      name = "targetOnly")

    let z = graphNode(scalarNodeLike(x), @[x], gtargetOnly)
    let dx = z.grad x
    z :~ a
    dx :~ 13.0

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
