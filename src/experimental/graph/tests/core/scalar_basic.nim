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

  test "update refreshes cached scalar values and gradients":
    let mutable = grt.toGvalue(2.0)
    let z = mutable * mutable
    let dzdx = z.grad mutable

    z :~ 4.0
    dzdx :~ 4.0

    mutable.update 5.0
    z :~ 25.0
    dzdx :~ 10.0

  test "direct scalar field writes do not mark freshness":
    let mutable = grt.toGvalue(2.0)
    let z = mutable + 1.0

    z :~ 3.0
    check mutable.sval == 2.0

    mutable.sval = 5.0
    check mutable.sval == 5.0
    z :~ 3.0

    mutable.update 5.0
    z :~ 6.0

  test "direct int field writes do not mark freshness":
    let selector = grt.toGvalue(0)
    let whenTrue = grt.toGvalue(10.0)
    let whenFalse = grt.toGvalue(20.0)
    let z = cond(selector, whenTrue, whenFalse)

    z :~ 20.0
    check selector.ival == 0

    selector.ival = 1
    check selector.ival == 1
    z :~ 20.0

    selector.update 1
    z :~ 10.0

  test "scalar and int erased copy compatibility stays direct":
    let scalarDst = scalarNodeLike(x)
    scalarDst.valCopy(x)
    check scalarDst.sval == x.sval
    check scalarDst.copyCompatible(x)

    let intSrc = grt.toGvalue(7)
    let intDst = intNodeLike(intSrc)
    intDst.valCopy(intSrc)
    check intDst.ival == intSrc.ival
    check intDst.copyCompatible(intSrc)

    check not scalarDst.copyCompatible(intSrc)
    check not intDst.copyCompatible(x)

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

  test "run counts are per node, not per graph function":
    let first = x + y
    let second = x + y
    let firstId = first.stableNodeId
    let secondId = second.stableNodeId

    check firstId != secondId

    check first.runCount == 0
    check second.runCount == 0

    discard first.eval
    check first.runCount == 1
    check second.runCount == 0

    discard second.eval
    check first.runCount == 1
    check second.runCount == 1

    x.update a + 1.0
    discard first.eval
    check first.runCount == 2
    check second.runCount == 1

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

  test "scalar validators reject wrong value type":
    let intValue: Gvalue = grt.toGvalue(1)
    let scalarValue: Gvalue = grt.toGvalue(1.0)
    expect(GraphValueError):
      discard intValue.requireScalar("requireScalar")
    expect(GraphValueError):
      discard scalarValue.requireInt("requireInt")

  test "separate runtimes isolate node ids and reject mixed graphs":
    let leftGrt = initGraphRuntime()
    let rightGrt = initGraphRuntime()

    let xLeft = scalarNodeIn(leftGrt)
    let yLeft = scalarNodeIn(leftGrt)
    xLeft.update 2.0
    yLeft.update 3.0

    let xRight = scalarNodeIn(rightGrt)
    xRight.update 5.0

    let xLeftId = xLeft.stableNodeId
    let yLeftId = yLeft.stableNodeId
    let xRightId = xRight.stableNodeId
    check xLeftId != yLeftId
    check xLeftId == xRightId

    let zLeft = xLeft + yLeft
    zLeft :~ 5.0

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

  test "findGrad returns cached gradient for matching runtime":
    let localX = grt.toGvalue(2.0)
    let localY = grt.toGvalue(3.0)
    let z = localX * localY
    let dzdx = z.grad localX

    discard dzdx.eval
    let found = findGrad(localX, z)
    check found != nil
    if found != nil:
      check sameNode(found, dzdx)

  test "graph construction rejects values without runtimes":
    let rawNode = Gvalue()
    let rawInput = Gvalue()
    let gnoop = newGfunc(name = "runtime-less")

    expect(GraphValueError):
      discard graphNode(rawNode, @[x], gnoop, "runtime-less result")
    expect(GraphValueError):
      discard graphNode(rawGraphValueIn(grt), @[rawInput], gnoop, "runtime-less input")

  test "runtime attachment rejects nil values":
    let missing: Gvalue = nil

    try:
      discard missing.attachRuntime(grt)
      check false
    except GraphValueError as e:
      check e.msg.contains("cannot attach graph runtime to nil value")

  test "graph construction rejects nil and runtime-less inputs":
    let gnoop = newGfunc(name = "input precheck")

    block:
      let missing: Gvalue = nil
      try:
        discard graphNode(rawGraphValueIn(grt), @[missing], gnoop, "nil input precheck")
        check false
      except GraphValueError as e:
        check e.msg.contains("nil input precheck input 0 cannot be nil")

    block:
      let rawInput = Gvalue()
      try:
        discard graphNode(
          rawGraphValueIn(grt),
          @[rawInput],
          gnoop,
          "runtime-less input precheck")
        check false
      except GraphValueError as e:
        check e.msg.contains("runtime-less input precheck input 0 has no graph runtime")

  test "graph construction rejects inputful nodes without graph functions":
    let missingFunc: Gfunc = nil

    try:
      discard graphNode(
        scalarNodeLike(x),
        @[Gvalue(x)],
        missingFunc,
        "nil graph function")
      check false
    except GraphValueError as e:
      check e.msg.contains("nil graph function with inputs requires a graph function")

  test "graph construction exposes checked node fields directly":
    proc copyForward(v: Gvalue) =
      v.valCopy v.inputs[0]

    let gcopy = newGfunc(forward = copyForward, name = "field copy")
    let z = graphNode(scalarNodeLike(x), @[Gvalue(x)], gcopy, "field copy")

    check z.inputs.len == 1
    check sameNode(z.inputs[0], x)
    check z.gfunc == gcopy
    check z.runtime == grt
    check z.epoch == 0
    z :~ a

  test "typed node validators return simple tuples":
    let z = x + y
    let view = z.requireBinaryNodeView(Gscalar, Gscalar, "binary tuple")

    check sameNode(view.x, x)
    check sameNode(view.y, y)

    expect(GraphValueError):
      discard x.requireBinaryNodeView(Gscalar, Gscalar, "wrong arity tuple")

    let intValue = grt.toGvalue(1)
    let mixed = graphNode(
      scalarNodeLike(x),
      @[Gvalue(x), Gvalue(intValue)],
      newGfunc(name = "mixed validator"),
      "mixed validator")

    expect(GraphValueError):
      discard mixed.requireBinaryNodeView(Gscalar, Gscalar, "wrong type tuple")

    try:
      discard mixed.requireNodeInput(3, "labeled input", "right")
      check false
    except GraphValueError as e:
      check e.msg.contains("right index out of range")

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

  test "literal scalar overloads and explicit erased scalar casts preserve concrete type":
    let erasedY: Gvalue = y
    let castY = erasedY.requireScalar("erased scalar test right")

    check not compiles(x + erasedY)
    check not compiles(x - erasedY)
    check not compiles(x * erasedY)
    check not compiles(x / erasedY)

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

    let addErased: Gscalar = x + castY
    let subErased: Gscalar = x - castY
    let mulErased: Gscalar = x * castY
    let divErased: Gscalar = x / castY
    let chained: Gscalar = (x + castY) * 2.0

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

  test "newGfunc rejects conflicting gradient hooks":
    proc rawBackward(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
      discard zb
      discard z
      discard i
      discard dep
      nil

    proc targetBackward(zb: Gvalue, z: Gvalue, target: Gvalue, dep: Gvalue): Gvalue =
      discard zb
      discard z
      discard target
      discard dep
      nil

    expect(GraphValueError):
      discard newGfunc(
        backward = rawBackward,
        backwardTarget = targetBackward,
        name = "badGradientHooks")

  test "newGfunc rejects blank names":
    for name in ["", "   "]:
      var failed = false
      try:
        discard newGfunc(name = name)
      except GraphValueError as e:
        failed = true
        check e.msg.contains("graph function name")
      check failed

  test "raw backward returning nil fails fast":
    proc copyForward(v: Gvalue) =
      v.valCopy v.inputs[0]

    proc nilBackward(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
      discard zb
      discard z
      discard i
      discard dep
      nil

    let gnilBackward = newGfunc(
      forward = copyForward,
      backward = nilBackward,
      name = "nilBackward")
    let z = graphNode(scalarNodeLike(x), @[x], gnilBackward, "nilBackward")

    expect(GraphValueError):
      discard z.grad x

  test "gradient planning skips irrelevant raw inputs":
    proc leftForward(v: Gvalue) =
      v.valCopy v.inputs[0]

    proc guardedBackward(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
      discard zb
      discard dep
      if i == 0:
        return z.oneLike
      raiseValueError("irrelevant input backward was called")

    let guarded = newGfunc(
      forward = leftForward,
      backward = guardedBackward,
      name = "guardedLeft")
    let z = graphNode(
      scalarNodeLike(x),
      @[Gvalue(x), Gvalue(y)],
      guarded,
      "guardedLeft")

    z.grad(x) :~ 1.0

  test "gradient planning rejects dependency cycles":
    let left = scalarNodeLike(x)
    let right = scalarNodeLike(x)
    let gcycle = newGfunc(name = "cycle")
    left.inputs = @[Gvalue(right)]
    left.gfunc = gcycle
    right.inputs = @[Gvalue(left)]
    right.gfunc = gcycle

    expect(GraphError):
      discard left.grad x

  test "eval rejects dependency cycles":
    let left = scalarNodeLike(x)
    let right = scalarNodeLike(x)
    let gcycle = newGfunc(name = "eval cycle")
    left.inputs = @[Gvalue(right)]
    left.gfunc = gcycle
    right.inputs = @[Gvalue(left)]
    right.gfunc = gcycle

    expect(GraphError):
      discard left.eval

  test "dependency walks reject nil dependencies at traversal boundary":
    proc walkNil(node: Gvalue, visit: GnodeVisit) =
      discard node
      let missing: Gvalue = nil
      visit missing

    proc expectNilDependency(mode: InputWalkMode,
                             walks: GdepWalks,
                             label: string) =
      let node = graphNode(
        scalarNodeLike(x),
        @[Gvalue(x)],
        newGfunc(depWalks = walks, name = label),
        label)
      try:
        discard node.collectNodeInputs(mode)
        check false
      except GraphValueError as e:
        check e.msg.contains("dependency walk produced nil dependency")

    expectNilDependency(iwmEval, GdepWalks(eval: walkNil), "nil eval dependency")
    expectNilDependency(
      iwmDepend,
      GdepWalks(depend: walkNil),
      "nil depend dependency")
    expectNilDependency(
      iwmGradSignature,
      GdepWalks(gradSignature: walkNil),
      "nil signature dependency")

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

  test "custom depend walk does not make hidden deps raw backward inputs":
    let hidden = grt.toGvalue(5.0)

    proc copyForward(v: Gvalue) =
      v.valCopy v.inputs[0]

    proc walkHiddenDepend(node: Gvalue, visit: GnodeVisit) =
      visit node.inputs[0]
      visit hidden

    proc rawBackward(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
      discard dep
      if i != 0:
        raiseValueError("custom depend raw backward input index")
      if zb == nil:
        return z.oneLike
      zb

    let rawFunc = newGfunc(
      forward = copyForward,
      backward = rawBackward,
      depWalks = GdepWalks(depend: walkHiddenDepend),
      name = "custom depend raw backward")
    let rawNode = graphNode(scalarNodeLike(x), @[Gvalue(x)], rawFunc, "custom depend raw")

    rawNode :~ a
    rawNode.grad(x) :~ 1.0
    rawNode.grad(hidden) :~ 0.0

    proc targetBackward(zb: Gvalue,
                        z: Gvalue,
                        target: Gvalue,
                        dep: Gvalue): Gvalue =
      discard z
      discard dep
      if not sameNode(target, hidden):
        return nil
      let contribution = hidden.oneLike
      if zb == nil:
        return contribution
      contribution.scaleLike zb

    let targetFunc = newGfunc(
      forward = copyForward,
      backwardTarget = targetBackward,
      depWalks = GdepWalks(depend: walkHiddenDepend),
      name = "custom depend target backward")
    let targetNode = graphNode(scalarNodeLike(x), @[Gvalue(x)], targetFunc, "custom depend target")

    targetNode :~ a
    targetNode.grad(hidden) :~ 1.0

  test "nil backwardTarget contribution is treated as zero":
    proc targetNilf(v: Gvalue) =
      v.valCopy v.inputs[0]

    proc targetNilBackwardTarget(zb: Gvalue,
                                 z: Gvalue,
                                 target: Gvalue,
                                 dep: Gvalue): Gvalue =
      discard zb
      discard z
      discard target
      discard dep
      nil

    let gtargetNil = newGfunc(
      forward = targetNilf,
      backwardTarget = targetNilBackwardTarget,
      name = "targetNil")

    let z = graphNode(scalarNodeLike(x), @[x], gtargetNil)
    let dx = z.grad x
    z :~ a
    dx :~ 0.0

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
