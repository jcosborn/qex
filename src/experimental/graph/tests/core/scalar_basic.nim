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

  test "mutable zero upstream stays live in cached gradient":
    grt.resetGradCache()

    let s = grt.toGvalue(0.0)
    let x2 = grt.toGvalue(2.0)
    let z = s * (x2 * x2)
    let dzdx = z.grad x2
    let dzdxPointer = cast[pointer](dzdx)

    dzdx :~ 0.0

    s.update 3.0
    dzdx :~ 12.0

    let dzdxAgain = z.grad x2
    dzdxAgain :~ 12.0
    check cast[pointer](dzdxAgain) == dzdxPointer

  test "explicit static zero marker distinguishes mutable zero leaves":
    let mutableScalar = scalarNodeLike(x)
    let mutableInt = intNodeLike(grt.toGvalue(0))
    let scalarZero = x.zeroLike
    let intZero = grt.toGvalue(1).zeroLike

    check mutableScalar.isZero
    check mutableInt.isZero
    check not mutableScalar.isStaticZeroLeaf
    check not mutableInt.isStaticZeroLeaf
    check scalarZero.isStaticZeroLeaf
    check intZero.isStaticZeroLeaf

    Gscalar(scalarZero).update 2.0
    check not scalarZero.isStaticZeroLeaf

  test "computed zero graph nodes are not static zero leaves":
    let z = x - x

    z :~ 0.0
    check z.isZero
    check not z.isStaticZeroLeaf

  test "stable node ids are assigned at construction and read-only":
    let beforeLeaf = grt.nextStableNodeId
    let leaf = Gscalar(runtime: grt).assignStableNodeId
    let leafId = leaf.stableNodeId

    check leafId == beforeLeaf + 1
    check grt.nextStableNodeId == leafId
    check leaf.stableNodeId == leafId
    check grt.nextStableNodeId == leafId

    let beforeNode = grt.nextStableNodeId
    let z = x + y
    let nodeId = z.stableNodeId

    check nodeId == beforeNode + 1
    check grt.nextStableNodeId == nodeId
    check z.stableNodeId == nodeId
    check grt.nextStableNodeId == nodeId

    expect(GraphValueError):
      discard Gscalar().stableNodeId

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

  test "graphNode rejects nil values and nil runtimes":
    let copyFunc = Gfunc(
      forward: proc(v: Gvalue) =
        v.valCopy(v.inputs[0]),
      name: "construction test copy")
    let missing: Gvalue = nil
    let raw = Gscalar()
    let rawWithRuntime = Gscalar(runtime: grt)

    expect(GraphValueError):
      discard graphNode(
        scalarNodeLike(x),
        [missing],
        copyFunc,
        "nil input graph node")

    expect(GraphValueError):
      discard graphNode(
        scalarNodeLike(x),
        [Gvalue(raw)],
        copyFunc,
        "raw input graph node")

    expect(GraphValueError):
      discard graphNode(
        scalarNodeLike(x),
        [Gvalue(rawWithRuntime)],
        copyFunc,
        "unconstructed input graph node")

    expect(GraphValueError):
      discard graphNode(
        Gscalar(nil),
        [Gvalue(x)],
        copyFunc,
        "nil result graph node")

    expect(GraphValueError):
      discard graphNode(
        raw,
        [Gvalue(x)],
        copyFunc,
        "raw result graph node")

  test "gradSeeded returns compatible same-node seed":
    let seed = grt.toGvalue(3.0)
    let result = x.gradSeeded(x, seed)

    check result.nodeKey == seed.nodeKey
    result :~ 3.0

  test "gradSeeded rejects seed incompatible with output":
    let badSeed = grt.toGvalue(1)

    try:
      discard x.gradSeeded(x, badSeed)
      fail()
    except GraphValueError as e:
      check e.msg.contains("gradSeeded seed is incompatible with output")

    let z = x + y
    try:
      discard z.gradSeeded(x, badSeed)
      fail()
    except GraphValueError as e:
      check e.msg.contains("gradSeeded seed is incompatible with output")

  test "gradSeeded with unit seed matches ordinary scalar gradient":
    let z = x * y + x
    let seeded = z.gradSeeded(x, grt.toGvalue(1.0))
    let ordinary = z.grad x

    z :~ a * b + a
    seeded :~ b + 1.0
    ordinary :~ b + 1.0

  test "gradSeeded unreachable target returns zero without populating cache":
    grt.resetGradCache()

    let z = x * x
    let seeded = z.gradSeeded(y, grt.toGvalue(7.0))

    seeded :~ 0.0
    check findGrad(y, z) == nil
    check grt.gradCacheStats == GradCacheStats()

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

  test "eval freshness is per node, not per graph function":
    proc addForward(v: Gvalue) =
      Gscalar(v).sval = Gscalar(v.inputs[0]).sval + Gscalar(v.inputs[1]).sval
    let f = Gfunc(forward: addForward, name: "counted add")

    let first = graphNode(scalarNodeLike(x), @[Gvalue(x), Gvalue(y)], f, "counted first add")
    let second = graphNode(scalarNodeLike(x), @[Gvalue(x), Gvalue(y)], f, "counted second add")
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

  test "separate runtimes isolate node ids and reject mixed graphs":
    let leftGrt = initGraphRuntime()
    let rightGrt = initGraphRuntime()

    let xLeft = leftGrt.toGvalue(2.0)
    let yLeft = leftGrt.toGvalue(3.0)

    let xRight = rightGrt.toGvalue(5.0)

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
    discard xLeft.stableNodeId
    discard xRight.stableNodeId
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
      check found.nodeKey == dzdx.nodeKey

  test "findGrad ignores stale symbolic revisions":
    grt.resetGradCache()
    let localX = grt.toGvalue(2.0)
    let localY = grt.toGvalue(3.0)
    let z = localX * localY
    discard z.grad localX

    check findGrad(localX, z) != nil
    inc grt.symbolicRevision
    check findGrad(localX, z) == nil

  test "same-node gradients keep revision cache current":
    grt.resetGradCache()
    let localX = grt.toGvalue(2.0)

    let dx = localX.grad localX
    dx :~ 1.0
    check grt.gradCacheStats.revisionMisses == 1
    check grt.gradCacheStats.directHits == 0
    check grt.gradCacheStats.invalidations == 0

    let dxAgain = localX.grad localX
    dxAgain :~ 1.0
    check grt.gradCacheStats.revisionHits == 1
    check grt.gradCacheStats.directHits == 1
    check grt.gradCacheStats.revisionMisses == 1

  test "gradient cache reuses intermediate adjoints across targets and deps":
    var fLeftCalls = 0
    var fRightCalls = 0
    var gLeftCalls = 0
    var gRightCalls = 0
    var hCalls = 0

    proc upstreamOrOne(zb: Gvalue, z: Gvalue): Gscalar =
      if zb == nil:
        return toGvalue(z.runtime, 1.0)
      Gscalar(zb)

    proc addForward(v: Gvalue) =
      let left = Gscalar(v.inputs[0])
      let right = Gscalar(v.inputs[1])
      Gscalar(v).sval = left.sval + right.sval

    proc addBackward(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
      let upstream = upstreamOrOne(zb, z)
      case i
      of 0:
        inc fLeftCalls
        upstream
      of 1:
        inc fRightCalls
        upstream
      else:
        raiseValueError("tracked add input index")

    proc mulForward(v: Gvalue) =
      let left = Gscalar(v.inputs[0])
      let right = Gscalar(v.inputs[1])
      Gscalar(v).sval = left.sval * right.sval

    proc mulBackward(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
      let left = Gscalar(z.inputs[0])
      let right = Gscalar(z.inputs[1])
      let upstream = upstreamOrOne(zb, z)
      case i
      of 0:
        inc gLeftCalls
        upstream * right
      of 1:
        inc gRightCalls
        upstream * left
      else:
        raiseValueError("tracked mul input index")

    proc squareForward(v: Gvalue) =
      let inputNode = v.inputs[0]
      let input = Gscalar(inputNode)
      Gscalar(v).sval = input.sval * input.sval

    proc squareBackward(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
      if i != 0:
        raiseValueError("tracked square input index")
      inc hCalls
      let inputNode = z.inputs[0]
      let input = Gscalar(inputNode)
      upstreamOrOne(zb, z) * input * 2.0

    let gFunc = Gfunc(
      forward: mulForward,
      backward: mulBackward,
      name: "tracked g")
    let hFunc = Gfunc(
      forward: squareForward,
      backward: squareBackward,
      name: "tracked h")
    let fFunc = Gfunc(
      forward: addForward,
      backward: addBackward,
      name: "tracked f")

    let gNode = graphNode(
      scalarNodeLike(x),
      @[Gvalue(x), Gvalue(y)],
      gFunc,
      "tracked g")
    let hNode = graphNode(scalarNodeLike(y), @[Gvalue(y)], hFunc, "tracked h")
    let dep = graphNode(
      scalarNodeLike(gNode),
      @[Gvalue(gNode), Gvalue(hNode)],
      fFunc,
      "tracked f")

    let dx = dep.grad x
    dx :~ b
    check fLeftCalls == 1
    check fRightCalls == 0
    check gLeftCalls == 1
    check gRightCalls == 0
    check hCalls == 0

    let cachedGAdjoint = findGrad(gNode, dep)
    check cachedGAdjoint != nil

    let e = dep + dx
    let dedy = e.grad y
    dedy :~ a + 2.0 * b + 1.0

    let cAdjointInE = findGrad(dep, e)
    check cAdjointInE != nil
    let gAdjointInE = findGrad(gNode, e)
    check gAdjointInE != nil

    let callsBeforeDy = (
      fLeft: fLeftCalls,
      fRight: fRightCalls,
      gLeft: gLeftCalls,
      gRight: gRightCalls,
      h: hCalls)

    let dy = dep.grad y
    dy :~ a + 2.0 * b
    check fLeftCalls == callsBeforeDy.fLeft
    check fRightCalls == callsBeforeDy.fRight + 1
    check gLeftCalls == callsBeforeDy.gLeft
    check gRightCalls == callsBeforeDy.gRight + 1
    check hCalls == callsBeforeDy.h + 1

    let reusedGAdjoint = findGrad(gNode, dep)
    check reusedGAdjoint != nil
    if cachedGAdjoint != nil and reusedGAdjoint != nil:
      check cachedGAdjoint.nodeKey == reusedGAdjoint.nodeKey

    let z = dedy + dy
    z :~ 2.0 * a + 4.0 * b + 1.0

  test "failed gradient expansion is not cached as expanded":
    var failRight = true
    var rightCalls = 0

    proc addForward(v: Gvalue) =
      let left = Gscalar(v.inputs[0])
      let right = Gscalar(v.inputs[1])
      Gscalar(v).sval = left.sval + right.sval

    proc addBackward(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
      let upstream =
        if zb == nil:
          toGvalue(z.runtime, 1.0)
        else:
          Gscalar(zb)
      case i
      of 0:
        upstream
      of 1:
        inc rightCalls
        if failRight:
          raiseValueError("flaky right gradient")
        upstream
      else:
        raiseValueError("flaky add input index")

    let flakyAdd = Gfunc(
      forward: addForward,
      backward: addBackward,
      name: "flaky add")
    let z = graphNode(
      scalarNodeLike(x),
      @[Gvalue(x), Gvalue(y)],
      flakyAdd,
      "flaky add")

    z.grad(x) :~ 1.0
    check rightCalls == 0

    expect(GraphValueError):
      discard z.grad(y)
    check rightCalls == 1
    check findGrad(y, z) == nil

    failRight = false
    z.grad(y) :~ 1.0
    check rightCalls == 2

  test "graph values are constructed with explicit runtimes":
    let node = Gvalue(runtime: grt)
    let gnoop = Gfunc(name: "explicit runtime")
    let z = graphNode(node, newSeq[Gvalue](), gnoop, "explicit runtime")

    check z.runtime == grt
    check z.gfunc == gnoop

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

    let gcopy = Gfunc(forward: copyForward, name: "field copy")
    let z: Gscalar = graphNode(scalarNodeLike(x), @[Gvalue(x)], gcopy, "field copy")

    check z.inputs.len == 1
    check z.inputs[0].nodeKey == x.nodeKey
    check z.gfunc == gcopy
    check z.runtime == grt
    check z.epoch == 0
    z :~ a

  test "eval keeps graph node topology for ordinary forward hooks":
    proc copyForward(v: Gvalue) =
      v.valCopy v.inputs[0]

    let copying = Gfunc(
      forward: copyForward,
      name: "copy forward")
    let z = graphNode(
      scalarNodeLike(x),
      @[Gvalue(x)],
      copying,
      "copy forward")

    discard z.eval
    check z.gfunc == copying
    check z.inputs.len == 1
    check z.inputs[0].nodeKey == x.nodeKey

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
    let castY = Gscalar(erasedY)

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

  test "Gfunc accepts only graph operation hooks":
    check not compiles(Gfunc(
      cloneInputs = proc(src: Gvalue, clonedInputs: seq[Gvalue]): seq[Gvalue] =
        discard src
        clonedInputs,
      name: "clone hook"))

  test "backward mode exposes dynamic deps without mutating raw inputs":
    var sawDynamicInput = false

    proc copyForward(v: Gvalue) =
      v.valCopy v.inputs[0]

    proc dynamicInputView(v: Gvalue, mode: InputWalkMode, visit: GnodeVisit) =
      case mode
      of iwmBackward:
        visit y
        visit v.inputs[0]
      of iwmEval, iwmReachable:
        visit v.inputs[0]

    proc dynamicBackward(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
      discard zb
      discard z
      case i
      of 0:
        check input.nodeKey == y.nodeKey
        sawDynamicInput = true
        toGvalue(y.runtime, 7.0)
      of 1:
        check input.nodeKey == x.nodeKey
        toGvalue(x.runtime, 1.0)
      else:
        raiseValueError("dynamic input index")

    let dynamicFunc = Gfunc(
      forward: copyForward,
      backward: dynamicBackward,
      inputView: dynamicInputView,
      name: "dynamic backward deps")
    let z = graphNode(
      scalarNodeLike(x),
      @[Gvalue(x)],
      dynamicFunc,
      "dynamic backward deps")

    check z.inputs.len == 1
    check z.inputs[0].nodeKey == x.nodeKey
    z.grad(y) :~ 7.0
    check sawDynamicInput
    check z.inputs.len == 1
    check z.inputs[0].nodeKey == x.nodeKey

  test "backward mode rejects cross-runtime values":
    proc copyForward(v: Gvalue) =
      v.valCopy v.inputs[0]

    proc passthroughBackward(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
      discard i
      discard input
      if zb == nil:
        return z.oneLike
      zb

    block:
      let other = initGraphRuntime().toGvalue(1.0)

      proc mixedInputView(v: Gvalue, mode: InputWalkMode, visit: GnodeVisit) =
        case mode
        of iwmBackward:
          visit other
        of iwmEval, iwmReachable:
          visit v.inputs[0]

      let mixedFunc = Gfunc(
        forward: copyForward,
        backward: passthroughBackward,
        inputView: mixedInputView,
        name: "mixed backward deps")
      let z = graphNode(
        scalarNodeLike(x),
        @[Gvalue(x)],
        mixedFunc,
        "mixed backward deps")

      expect(GraphValueError):
        discard z.grad x

  test "raw backward returning nil fails fast":
    proc copyForward(v: Gvalue) =
      v.valCopy v.inputs[0]

    proc nilBackward(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
      discard zb
      discard z
      discard i
      nil

    let gnilBackward = Gfunc(
      forward: copyForward,
      backward: nilBackward,
      name: "nilBackward")
    let z = graphNode(scalarNodeLike(x), @[x], gnilBackward, "nilBackward")

    expect(GraphValueError):
      discard z.grad x

  test "malformed intermediate backward shape fails at its edge":
    proc copyForward(v: Gvalue) =
      v.valCopy v.inputs[0]

    proc intBackward(zb: Gvalue, z: Gvalue,
                     i: int, input: Gvalue): Gvalue =
      discard zb
      discard i
      discard input
      Gvalue(toGvalue(z.runtime, 1))

    let
      middle = exp(x)
      malformed = graphNode(
        scalarNodeLike(middle),
        @[Gvalue(middle)],
        Gfunc(
          forward: copyForward,
          backward: intBackward,
          name: "wrong-shaped backward"),
        "wrong-shaped backward")

    expect(GraphValueError):
      discard malformed.grad x

  test "gradient planning skips irrelevant raw inputs":
    proc leftForward(v: Gvalue) =
      v.valCopy v.inputs[0]

    proc guardedBackward(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
      discard zb
      if i == 0:
        return z.oneLike
      raiseValueError("irrelevant input backward was called")

    let guarded = Gfunc(
      forward: leftForward,
      backward: guardedBackward,
      name: "guardedLeft")
    let z = graphNode(
      scalarNodeLike(x),
      @[Gvalue(x), Gvalue(y)],
      guarded,
      "guardedLeft")

    z.grad(x) :~ 1.0

  test "gradient planning rejects dependency cycles":
    var left: Gscalar
    var right: Gscalar
    proc leftInputView(node: Gvalue, mode: InputWalkMode, visit: GnodeVisit) =
      discard node
      discard mode
      visit right
    proc rightInputView(node: Gvalue, mode: InputWalkMode, visit: GnodeVisit) =
      discard node
      discard mode
      visit left

    left = graphNode(
      scalarNodeLike(x),
      newSeq[Gvalue](),
      Gfunc(inputView: leftInputView, name: "cycle left"),
      "cycle left")
    right = graphNode(
      scalarNodeLike(x),
      newSeq[Gvalue](),
      Gfunc(inputView: rightInputView, name: "cycle right"),
      "cycle right")

    expect(GraphError):
      discard left.grad x

  test "eval rejects dependency cycles":
    var left: Gscalar
    var right: Gscalar
    proc leftInputView(node: Gvalue, mode: InputWalkMode, visit: GnodeVisit) =
      discard node
      discard mode
      visit right
    proc rightInputView(node: Gvalue, mode: InputWalkMode, visit: GnodeVisit) =
      discard node
      discard mode
      visit left

    left = graphNode(
      scalarNodeLike(x),
      newSeq[Gvalue](),
      Gfunc(inputView: leftInputView, name: "eval cycle left"),
      "eval cycle left")
    right = graphNode(
      scalarNodeLike(x),
      newSeq[Gvalue](),
      Gfunc(inputView: rightInputView, name: "eval cycle right"),
      "eval cycle right")

    expect(GraphError):
      discard left.eval

  test "input views reject nil dependencies at traversal boundary":
    proc walkNil(node: Gvalue, mode: InputWalkMode, visit: GnodeVisit) =
      discard node
      discard mode
      let missing: Gvalue = nil
      visit missing

    proc expectNilDependency(mode: InputWalkMode,
                             inputView: GinputViewHook,
                             label: string) =
      let node = graphNode(
        scalarNodeLike(x),
        @[Gvalue(x)],
        Gfunc(inputView: inputView, name: label),
        label)
      try:
        discard node.collectInputView(mode)
        check false
      except GraphValueError as e:
        check e.msg.contains("input view produced nil dependency")

    expectNilDependency(iwmEval, walkNil, "nil eval dependency")
    expectNilDependency(iwmReachable, walkNil, "nil reachable dependency")
    expectNilDependency(iwmBackward, walkNil, "nil backward dependency")

  test "input views reject unconstructed dependencies at traversal boundary":
    let raw = Gscalar(runtime: grt)

    proc walkRaw(node: Gvalue, mode: InputWalkMode, visit: GnodeVisit) =
      discard node
      discard mode
      visit raw

    proc expectRawDependency(mode: InputWalkMode,
                             label: string) =
      let node = graphNode(
        scalarNodeLike(x),
        @[Gvalue(x)],
        Gfunc(inputView: walkRaw, name: label),
        label)
      try:
        discard node.collectInputView(mode)
        check false
      except GraphValueError as e:
        check e.msg.contains("graph value has no stable node id")

    expectRawDependency(iwmEval, "raw eval dependency")
    expectRawDependency(iwmReachable, "raw reachable dependency")
    expectRawDependency(iwmBackward, "raw backward dependency")

  test "input views reject cross-runtime dependencies at traversal boundary":
    let other = initGraphRuntime().toGvalue(9.0)

    proc walkOtherRuntime(node: Gvalue,
                          mode: InputWalkMode,
                          visit: GnodeVisit) =
      discard node
      discard mode
      visit other

    proc expectCrossRuntimeDependency(mode: InputWalkMode,
                                      label: string) =
      let node = graphNode(
        scalarNodeLike(x),
        @[Gvalue(x)],
        Gfunc(inputView: walkOtherRuntime, name: label),
        label)
      try:
        discard node.collectInputView(mode)
        check false
      except GraphValueError as e:
        check e.msg.contains("mixes graph runtimes")

    expectCrossRuntimeDependency(iwmEval, "foreign eval dependency")
    expectCrossRuntimeDependency(iwmReachable, "foreign reachable dependency")
    expectCrossRuntimeDependency(iwmBackward, "foreign backward dependency")

  test "custom reachable view does not make extra deps raw backward deps":
    let hidden = grt.toGvalue(5.0)

    proc copyForward(v: Gvalue) =
      v.valCopy v.inputs[0]

    proc walkExtraReachable(node: Gvalue, mode: InputWalkMode, visit: GnodeVisit) =
      visit node.inputs[0]
      if mode == iwmReachable:
        visit hidden

    proc rawBackward(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
      if i != 0:
        raiseValueError("custom reachable raw backward dep index")
      if zb == nil:
        return z.oneLike
      zb

    let rawFunc = Gfunc(
      forward: copyForward,
      backward: rawBackward,
      inputView: walkExtraReachable,
      name: "custom reachable raw backward")
    let rawNode = graphNode(scalarNodeLike(x), @[Gvalue(x)], rawFunc, "custom reachable raw")

    rawNode :~ a
    rawNode.grad(x) :~ 1.0
    rawNode.grad(hidden) :~ 0.0

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
