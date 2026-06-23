proc copyMultiForward(v: Gvalue) =
  let z = Gmulti(v)
  for i in 0..<z.inputs.len:
    z.storedSlot(i).valCopy(z.inputs[i])

proc copyMultiBackward(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let zb = Gmulti(zb)
  zb[i]

let copyMultiFunc = Gfunc(
  forward: copyMultiForward,
  backward: copyMultiBackward,
  name: "copyMulti")

suite "graph multi":
  test "selection forwards and scatters gradients":
    let x = grt.toGvalue(2.0)
    let y = grt.toGvalue(3.0)
    let pair = multiValues("pair", x, y)
    let first = pair[0]
    let second = pair[1]

    first :~ 2.0
    second :~ 3.0
    grad(first, x) :~ 1.0
    grad(first, y) :~ 0.0
    grad(second, x) :~ 0.0
    grad(second, y) :~ 1.0

  test "multi addLike keeps slotwise forward and backward behavior":
    let x = grt.toGvalue(2.0)
    let y = grt.toGvalue(3.0)
    let left = multiValues("left", x, y)
    let right = multiValues("right", y, x)
    let added = Gmulti(left.addLike(left, right))

    added[0] :~ 5.0
    added[1] :~ 5.0
    grad(added[0], x) :~ 1.0
    grad(added[0], y) :~ 1.0
    grad(added[1], x) :~ 1.0
    grad(added[1], y) :~ 1.0

  test "cond over multi carriers selects slots and gradients":
    let x0 = grt.toGvalue(2.0)
    let x1 = grt.toGvalue(3.0)
    let y0 = grt.toGvalue(5.0)
    let y1 = grt.toGvalue(7.0)
    let k = grt.toGvalue(1)
    let left = multiValues("cond left", x0, x1)
    let right = multiValues("cond right", y0, y1)
    let selected = cond(k, left, right)
    let first = selected[0]
    let second = selected[1]

    first :~ 2.0
    second :~ 3.0
    grad(first, x0) :~ 1.0
    grad(first, y0) :~ 0.0
    grad(second, x1) :~ 1.0
    grad(second, y1) :~ 0.0

    k.update 0
    first :~ 5.0
    second :~ 7.0
    grad(first, x0) :~ 0.0
    grad(first, y0) :~ 1.0
    grad(second, x1) :~ 0.0
    grad(second, y1) :~ 1.0

  test "oneLike and zeroLike preserve multi shape":
    let x = grt.toGvalue(2.0)
    let y = grt.toGvalue(3.0)
    let pair = multiValues("pair", x, y)
    let onePair = Gmulti(pair.oneLike)
    let zeroPair = Gmulti(pair.zeroLike)

    onePair[0] :~ 1.0
    onePair[1] :~ 1.0
    zeroPair[0] :~ 0.0
    zeroPair[1] :~ 0.0
    grad(onePair[0], x) :~ 0.0
    grad(onePair[1], y) :~ 0.0

  test "selection exposes only the base as a graph dependency":
    let x = grt.toGvalue(2.0)
    let y = grt.toGvalue(3.0)
    let pair = multiValues("pair", x, y)
    let first = pair[0]
    let deps = first.collectInputView(iwmReachable)

    check deps.len == 1
    if deps.len > 0:
      check deps[0].nodeKey == pair.nodeKey

  test "selection slot index stays fixed by the selection node":
    let x = grt.toGvalue(2.0)
    let y = grt.toGvalue(3.0)
    let pair = multiValues("pair", x, y)
    let first = pair[0]
    let firstAgain = pair[0]
    let second = pair[1]

    first :~ 2.0
    firstAgain :~ 2.0
    second :~ 3.0
    grad(first, x) :~ 1.0
    grad(firstAgain, y) :~ 0.0
    grad(second, y) :~ 1.0

  test "multi addLike rejects mismatched arity":
    let x = grt.toGvalue(2.0)
    let y = grt.toGvalue(3.0)
    let left = multiValues("left", x, y)
    let right = multiValues("right", x)

    expect(GraphValueError):
      discard left.addLike(left, right)

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

  test "multi constructor rejects nil values and nil runtimes":
    let x = grt.toGvalue(2.0)
    let missing: Gvalue = nil
    let raw = Gscalar()
    let rawWithRuntime = Gscalar(runtime: grt)

    expect(GraphValueError):
      discard newMultiOutputNode(
        [missing],
        [Gvalue(x)],
        copyMultiFunc,
        "nil slot multi")

    expect(GraphValueError):
      discard newMultiOutputNode(
        [Gvalue(x)],
        [missing],
        copyMultiFunc,
        "nil input multi")

    expect(GraphValueError):
      discard newMultiOutputNode(
        [Gvalue(raw)],
        [Gvalue(x)],
        copyMultiFunc,
        "raw slot multi")

    expect(GraphValueError):
      discard newMultiOutputNode(
        [Gvalue(x)],
        [Gvalue(raw)],
        copyMultiFunc,
        "raw input multi")

    expect(GraphValueError):
      discard newMultiOutputNode(
        [Gvalue(rawWithRuntime)],
        [Gvalue(x)],
        copyMultiFunc,
        "unconstructed slot multi")

    expect(GraphValueError):
      discard newMultiOutputNode(
        [Gvalue(x)],
        [Gvalue(rawWithRuntime)],
        copyMultiFunc,
        "unconstructed input multi")

  test "carrier supports heterogeneous slot prototypes and per-slot symbolic access":
    let scalar = grt.toGvalue(2.0)
    let intValue = grt.toGvalue(3)
    let mixed = multiValues("mixed", scalar, intValue)
    let zeroMixed = Gmulti(mixed.zeroLike)
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

  test "multiValues root gradient seeds all slots":
    let x = grt.toGvalue(2.0)
    let y = grt.toGvalue(3.0)
    let pair = multiValues("root pair", x, y)

    grad(pair, x) :~ 1.0
    grad(pair, y) :~ 1.0

    let mixed = multiValues("root mixed", x, x * x)
    grad(mixed, x) :~ 5.0

  test "multiValues packs mixed varargs":
    let scalar = grt.toGvalue(2.0)
    let intValue = grt.toGvalue(3)
    let other = grt.toGvalue(4.0)
    let mixed = multiValues("mixed varargs", scalar, intValue, other)
    discard mixed.eval

    mixed.storedSlot(0) :~ 2.0
    mixed.storedSlot(1) :~ 3
    mixed.storedSlot(2) :~ 4.0

  test "multiValues rejects empty packs":
    expect(GraphValueError):
      discard multiValues("empty varargs")
