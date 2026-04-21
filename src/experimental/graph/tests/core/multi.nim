proc copyMultiForward(v: Gvalue) =
  let z = v.requireMultiValue("copyMulti forward")
  for i in 0..<z.inputs.len:
    z.slotValue(i).valCopy(z.inputs[i])

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
    let added = requireMultiValue(left + right, "multi add")

    added[0] :~ 5.0
    added[1] :~ 5.0
    grad(added[0], x) :~ 1.0
    grad(added[0], y) :~ 1.0
    grad(added[1], x) :~ 1.0
    grad(added[1], y) :~ 1.0

  test "constLike and zeroLike preserve multi shape":
    let x = grt.toGvalue(2.0)
    let y = grt.toGvalue(3.0)
    let pair = newScalarMulti([x, y], "pair")
    let constPair = requireMultiValue(pair.constLike(7), "const pair")
    let zeroPair = requireMultiValue(pair.zeroLike, "zero pair")

    constPair[0] :~ 7.0
    constPair[1] :~ 7.0
    zeroPair[0] :~ 0.0
    zeroPair[1] :~ 0.0
    grad(constPair[0], x) :~ 0.0
    grad(constPair[1], y) :~ 0.0

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

    mixed.slotValue(0) :~ 2.0
    mixed.slotValue(1) :~ 3
    zeroMixed.slotValue(0) :~ 0.0
    zeroMixed.slotValue(1) :~ 0
    mixed[0] :~ 2.0
    mixed[1] :~ 3
    zeroMixed[0] :~ 0.0
    zeroMixed[1] :~ 0
