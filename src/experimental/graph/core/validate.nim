import base_types

type
  UnaryNodeView* = object
    node*: Gvalue
    x*: Gvalue
  BinaryNodeView* = object
    node*: Gvalue
    x*: Gvalue
    y*: Gvalue
  TernaryNodeView* = object
    node*: Gvalue
    x*: Gvalue
    y*: Gvalue
    z*: Gvalue
  UnaryTypedNodeView*[X: Gvalue] = object
    node*: Gvalue
    x*: X
  BinaryTypedNodeView*[X: Gvalue, Y: Gvalue] = object
    node*: Gvalue
    x*: X
    y*: Y
  TernaryTypedNodeView*[X: Gvalue, Y: Gvalue, Z: Gvalue] = object
    node*: Gvalue
    x*: X
    y*: Y
    z*: Z

proc nodeContext*(node: Gvalue): string =
  if node == nil:
    return ""
  ":\n" & node.nodeRepr

proc requireInputCountAtLeast*(node: Gvalue,
                               minimum: int,
                               label: string) =
  if node.inputs.len < minimum:
    raiseValueError(
      label & " requires at least " & $minimum & " inputs, got " &
      $node.inputs.len &
      node.nodeContext)

proc requireInputCountExactly*(node: Gvalue,
                               expected: int,
                               label: string) =
  if node.inputs.len != expected:
    raiseValueError(
      label & " requires " & $expected & " inputs, got " &
      $node.inputs.len &
      node.nodeContext)

proc requireNodeInput*(node: Gvalue,
                       index: int,
                       label: string,
                       inputLabel = "input"): Gvalue =
  if index < 0 or index >= node.inputs.len:
    raiseValueError(
      label & " input index out of range: " & $index &
      " for " & $node.inputs.len &
      node.nodeContext)
  result = node.inputs[index]
  if result == nil:
    raiseError(label & " has nil " & inputLabel & " input" & node.nodeContext)

proc requireNodeInput*[T: Gvalue](node: Gvalue,
                                  index: int,
                                  inputType: typedesc[T],
                                  label: string,
                                  inputLabel = "input"): T =
  let value = node.requireNodeInput(index, label, inputLabel)
  if not (value of T):
    raiseValueError(
      label & " expects " & inputLabel & " input of type " & $inputType &
      ", got:\n" & value.nodeRepr)
  T(value)

proc requireUnaryNodeView*(node: Gvalue,
                           label: string,
                           inputLabel = "input"): UnaryNodeView =
  node.requireInputCountExactly(1, label)
  UnaryNodeView(
    node: node,
    x: node.requireNodeInput(0, label, inputLabel))

proc requireUnaryNodeView*[X: Gvalue](node: Gvalue,
                                      inputType: typedesc[X],
                                      label: string,
                                      inputLabel = "input"): UnaryTypedNodeView[X] =
  node.requireInputCountExactly(1, label)
  result.node = node
  result.x = node.requireNodeInput(0, inputType, label, inputLabel)

proc requireBinaryNodeView*(node: Gvalue,
                            label: string,
                            leftLabel = "left",
                            rightLabel = "right"): BinaryNodeView =
  node.requireInputCountExactly(2, label)
  BinaryNodeView(
    node: node,
    x: node.requireNodeInput(0, label, leftLabel),
    y: node.requireNodeInput(1, label, rightLabel))

proc requireBinaryNodeView*[X: Gvalue, Y: Gvalue](node: Gvalue,
                                                  leftType: typedesc[X],
                                                  rightType: typedesc[Y],
                                                  label: string,
                                                  leftLabel = "left",
                                                  rightLabel = "right"): BinaryTypedNodeView[X, Y] =
  node.requireInputCountExactly(2, label)
  result.node = node
  result.x = node.requireNodeInput(0, leftType, label, leftLabel)
  result.y = node.requireNodeInput(1, rightType, label, rightLabel)

proc requireTernaryNodeView*(node: Gvalue,
                             label: string,
                             firstLabel = "first",
                             secondLabel = "second",
                             thirdLabel = "third"): TernaryNodeView =
  node.requireInputCountExactly(3, label)
  TernaryNodeView(
    node: node,
    x: node.requireNodeInput(0, label, firstLabel),
    y: node.requireNodeInput(1, label, secondLabel),
    z: node.requireNodeInput(2, label, thirdLabel))

proc requireTernaryNodeView*[X: Gvalue, Y: Gvalue, Z: Gvalue](
    node: Gvalue,
    firstType: typedesc[X],
    secondType: typedesc[Y],
    thirdType: typedesc[Z],
    label: string,
    firstLabel = "first",
    secondLabel = "second",
    thirdLabel = "third"): TernaryTypedNodeView[X, Y, Z] =
  node.requireInputCountExactly(3, label)
  result.node = node
  result.x = node.requireNodeInput(0, firstType, label, firstLabel)
  result.y = node.requireNodeInput(1, secondType, label, secondLabel)
  result.z = node.requireNodeInput(2, thirdType, label, thirdLabel)

proc checkedInputValues*[T: Gvalue](values: openArray[T],
                                    label: string): seq[Gvalue] =
  result = newseq[Gvalue](values.len)
  for i in 0..<values.len:
    let value = values[i]
    if value == nil:
      raiseValueError(label & " element cannot be nil")
    result[i] = value
