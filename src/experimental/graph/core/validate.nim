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

proc requireUnaryNodeView*(node: Gvalue,
                           label: string,
                           inputLabel = "input"): UnaryNodeView =
  node.requireInputCountExactly(1, label)
  UnaryNodeView(
    node: node,
    x: node.requireNodeInput(0, label, inputLabel))

proc requireBinaryNodeView*(node: Gvalue,
                            label: string,
                            leftLabel = "left",
                            rightLabel = "right"): BinaryNodeView =
  node.requireInputCountExactly(2, label)
  BinaryNodeView(
    node: node,
    x: node.requireNodeInput(0, label, leftLabel),
    y: node.requireNodeInput(1, label, rightLabel))

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

proc checkedInputValues*(values: openArray[Gvalue],
                         label: string): seq[Gvalue] =
  result = newseq[Gvalue](values.len)
  for i in 0..<values.len:
    let value = values[i]
    if value == nil:
      raiseValueError(label & " element cannot be nil")
    result[i] = value
