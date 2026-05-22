import base_types

proc nodeContext*(node: Gvalue): string =
  if node == nil:
    return ""
  ":\n" & node.nodeRepr

proc requireGraphValue*(value: Gvalue,
                        label: string): Gvalue =
  if value == nil:
    raiseValueError(label & " cannot be nil")
  if value.runtime == nil:
    raiseValueError(label & " has no graph runtime")
  value

proc requireSameGraphRuntime*(left: Gvalue,
                              right: Gvalue,
                              label: string,
                              leftLabel = "left",
                              rightLabel = "right"): GraphRuntime =
  let checkedLeft = left.requireGraphValue(label & " " & leftLabel)
  let checkedRight = right.requireGraphValue(label & " " & rightLabel)
  result = checkedLeft.runtime
  if result != checkedRight.runtime:
    raiseValueError(label & " mixes multiple graph runtimes")

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
      label & " " & inputLabel & " index out of range: " & $index &
      " for " & $node.inputs.len &
      node.nodeContext)
  node.inputs[index]

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

proc requireUnaryNodeView*[X: Gvalue](node: Gvalue,
                                      inputType: typedesc[X],
                                      label: string,
                                      inputLabel = "input"): tuple[x: X] =
  node.requireInputCountExactly(1, label)
  result.x = node.requireNodeInput(0, inputType, label, inputLabel)

proc requireBinaryNodeView*[X: Gvalue, Y: Gvalue](node: Gvalue,
                                                  leftType: typedesc[X],
                                                  rightType: typedesc[Y],
                                                  label: string,
                                                  leftLabel = "left",
                                                  rightLabel = "right"): tuple[x: X, y: Y] =
  node.requireInputCountExactly(2, label)
  result.x = node.requireNodeInput(0, leftType, label, leftLabel)
  result.y = node.requireNodeInput(1, rightType, label, rightLabel)

proc requireTernaryNodeView*[X: Gvalue, Y: Gvalue, Z: Gvalue](
    node: Gvalue,
    firstType: typedesc[X],
    secondType: typedesc[Y],
    thirdType: typedesc[Z],
    label: string,
    firstLabel = "first",
    secondLabel = "second",
    thirdLabel = "third"): tuple[x: X, y: Y, z: Z] =
  node.requireInputCountExactly(3, label)
  result.x = node.requireNodeInput(0, firstType, label, firstLabel)
  result.y = node.requireNodeInput(1, secondType, label, secondLabel)
  result.z = node.requireNodeInput(2, thirdType, label, thirdLabel)

proc checkedInputValues*[T: Gvalue](values: openArray[T],
                                    label: string): seq[Gvalue] =
  result = newseq[Gvalue](values.len)
  for i in 0..<values.len:
    result[i] = values[i].requireGraphValue(label & " input " & $i)
