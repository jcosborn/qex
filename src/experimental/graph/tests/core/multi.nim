proc copyMultiForward(v: Gvalue) =
  let z = v.requireMultiValue("copyMulti forward")
  for i in 0..<z.inputs.len:
    z.storedSlot(i).valCopy(z.inputs[i])

proc copyMultiBackward(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard dep
  let upstream = requireMultiUpstream(zb, "copyMulti backward")
  upstream[i]

let copyMultiFunc = newGfunc(
  forward = copyMultiForward,
  backward = copyMultiBackward,
  name = "copyMulti")

proc newScalarMulti[T: Gvalue](values: openArray[T],
                               label = "scalar multi"): Gmulti =
  var inputs = newseq[Gvalue](values.len)
  for i in 0..<values.len:
    inputs[i] = values[i]
  newMultiOutputNode(inputs, inputs, copyMultiFunc, label)

suite "graph multi":
  test "selection forwards and scatters gradients":
    let x = grt.toGvalue(2.0)
    let y = grt.toGvalue(3.0)
    let pair = newScalarMulti([x, y], "pair")
    let first = pair[0]
    let second = pair[1]

    first :~ 2.0
    second :~ 3.0
    grad(first, x) :~ 1.0
    grad(first, y) :~ 0.0
    grad(second, x) :~ 0.0
    grad(second, y) :~ 1.0

  test "multi add keeps slotwise forward and backward behavior":
    let x = grt.toGvalue(2.0)
    let y = grt.toGvalue(3.0)
    let left = newScalarMulti([x, y], "left")
    let right = newScalarMulti([y, x], "right")
    let added = left + right

    added[0] :~ 5.0
    added[1] :~ 5.0
    grad(added[0], x) :~ 1.0
    grad(added[0], y) :~ 1.0
    grad(added[1], x) :~ 1.0
    grad(added[1], y) :~ 1.0

  test "oneLike and zeroLike preserve multi shape":
    let x = grt.toGvalue(2.0)
    let y = grt.toGvalue(3.0)
    let pair = newScalarMulti([x, y], "pair")
    let onePair = requireMultiValue(pair.oneLike, "one pair")
    let zeroPair = requireMultiValue(pair.zeroLike, "zero pair")

    onePair[0] :~ 1.0
    onePair[1] :~ 1.0
    zeroPair[0] :~ 0.0
    zeroPair[1] :~ 0.0
    grad(onePair[0], x) :~ 0.0
    grad(onePair[1], y) :~ 0.0

  test "selection rejects out-of-range indices early":
    let x = grt.toGvalue(2.0)
    let y = grt.toGvalue(3.0)
    let pair = newScalarMulti([x, y], "pair")

    expect(GraphValueError):
      discard pair[-1]
    expect(GraphValueError):
      discard pair[2]

  test "selection exposes only the base as a graph dependency":
    let x = grt.toGvalue(2.0)
    let y = grt.toGvalue(3.0)
    let pair = newScalarMulti([x, y], "pair")
    let first = pair[0]
    let deps = first.collectNodeInputs(iwmDepend)

    check deps.len == 1
    if deps.len > 0:
      check sameNode(deps[0], pair)

  test "selection slot index is signature metadata":
    let x = grt.toGvalue(2.0)
    let y = grt.toGvalue(3.0)
    let pair = newScalarMulti([x, y], "pair")
    let first = pair[0]
    let firstAgain = pair[0]
    let second = pair[1]

    check signatureKey(first.gfunc) == signatureKey(firstAgain.gfunc)
    check signatureKey(first.gfunc) != signatureKey(second.gfunc)

  test "multi add rejects mismatched arity":
    let x = grt.toGvalue(2.0)
    let y = grt.toGvalue(3.0)
    let left = newScalarMulti([x, y], "left")
    let right = newScalarMulti([x], "right")

    expect(GraphValueError):
      discard left + right

  test "multi constructor rejects slot and input runtime mismatch":
    let leftGrt = initGraphRuntime()
    let rightGrt = initGraphRuntime()
    let slot = leftGrt.toGvalue(0.0)
    let input = rightGrt.toGvalue(1.0)

    expect(GraphValueError):
      discard newMultiOutputNode(
        [Gvalue(slot)],
        [Gvalue(input)],
        copyMultiFunc,
        "mixed runtime multi")

  test "multi constructor rejects invalid slot prototypes":
    let x = grt.toGvalue(2.0)
    let rawSlot = Gvalue()
    let missingSlot: Gvalue = nil

    expect(GraphValueError):
      discard newMultiOutputNode(
        [rawSlot],
        [Gvalue(x)],
        copyMultiFunc,
        "runtime-less slot multi")
    expect(GraphValueError):
      discard newMultiOutputNode(
        [missingSlot],
        [Gvalue(x)],
        copyMultiFunc,
        "nil slot multi")

  test "multi constructor rejects inputful carriers without graph functions":
    let x = grt.toGvalue(2.0)
    let missingFunc: Gfunc = nil

    try:
      discard newMultiOutputNode(
        [Gvalue(x)],
        [Gvalue(x)],
        missingFunc,
        "nil multi function")
      check false
    except GraphValueError as e:
      check e.msg.contains("nil multi function with inputs requires a graph function")

  test "carrier supports heterogeneous slot prototypes and per-slot symbolic access":
    let scalar = grt.toGvalue(2.0)
    let intValue = grt.toGvalue(3)
    let mixed = newMultiOutputNode(
      [scalar, intValue],
      [scalar, intValue],
      copyMultiFunc,
      "mixed")
    let zeroMixed = requireMultiValue(mixed.zeroLike, "mixed zero")
    discard mixed.eval
    discard zeroMixed.eval

    mixed.storedSlot(0) :~ 2.0
    mixed.storedSlot(1) :~ 3
    zeroMixed.storedSlot(0) :~ 0.0
    zeroMixed.storedSlot(1) :~ 0
    mixed[0] :~ 2.0
    mixed[1] :~ 3
    zeroMixed[0] :~ 0.0
    zeroMixed[1] :~ 0

  test "multiValues packs varargs as symbolic slots":
    let x = grt.toGvalue(2.0)
    let y = grt.toGvalue(3.0)
    let pair = multiValues("vararg pair", x, y)

    pair[0] :~ 2.0
    pair[1] :~ 3.0
    grad(pair[0], x) :~ 1.0
    grad(pair[0], y) :~ 0.0
    grad(pair[1], x) :~ 0.0
    grad(pair[1], y) :~ 1.0

  test "multiValues packs mixed varargs":
    let scalar = grt.toGvalue(2.0)
    let intValue = grt.toGvalue(3)
    let other = grt.toGvalue(4.0)
    let mixed = multiValues("mixed varargs", scalar, intValue, other)
    discard mixed.eval

    mixed.storedSlot(0) :~ 2.0
    mixed.storedSlot(1) :~ 3
    mixed.storedSlot(2) :~ 4.0

  test "symbolicSlots returns typed tuple selections":
    let x = grt.toGvalue(2.0)
    let y = grt.toGvalue(3.0)
    let pair = multiValues("symbolic tuple pair", x, y)
    let (left, right) =
      symbolicSlots[tuple[left: Gscalar, right: Gscalar]](pair, "symbolic tuple")

    left :~ 2.0
    right :~ 3.0
    grad(left, x) :~ 1.0
    grad(left, y) :~ 0.0
    grad(right, x) :~ 0.0
    grad(right, y) :~ 1.0

  test "storedSlots returns typed tuple storage":
    let x = grt.toGvalue(2.0)
    let y = grt.toGvalue(3.0)
    let pair = multiValues("stored tuple pair", x, y)
    discard pair.eval
    let (left, right) =
      storedSlots[tuple[left: Gscalar, right: Gscalar]](pair, "stored tuple")

    left :~ 2.0
    right :~ 3.0

  test "tuple slot unpacking reports arity mismatches with caller label":
    let x = grt.toGvalue(2.0)
    let y = grt.toGvalue(3.0)
    let pair = multiValues("arity pair", x, y)

    try:
      discard symbolicSlots[
        tuple[left: Gscalar, middle: Gscalar, right: Gscalar]](pair, "arity tuple")
      check false
    except GraphValueError as e:
      check e.msg.contains("arity tuple expects 3 packed slots, got 2")

  test "tuple slot unpacking reports type mismatches with field label":
    let scalar = grt.toGvalue(2.0)
    let intValue = grt.toGvalue(3)
    let mixed = multiValues("type pair", scalar, intValue)

    try:
      discard symbolicSlots[
        tuple[left: Gscalar, right: Gscalar]](mixed, "typed tuple")
      check false
    except GraphValueError as e:
      check e.msg.contains("typed tuple.right expects Gscalar")

  test "multiValues rejects empty packs":
    expect(GraphValueError):
      discard multiValues("empty varargs")
