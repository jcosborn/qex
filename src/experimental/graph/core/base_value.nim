from strutils import join
import base_types

method newOneOf*(x: Gvalue): Gvalue {.base.} =
  raiseErrorBaseMethod("newOneOf(" & $x & ")")
method valCopy*(z: Gvalue, x: Gvalue) {.base.} =
  raiseErrorBaseMethod("valCopy(" & $z & "," & $x & ")")
method copyCompatible*(prototype: Gvalue, value: Gvalue): bool {.base.} =
  false
method zeroLike*(x: Gvalue): Gvalue {.base.} = x.newOneOf

method isZero*(x: Gvalue): bool {.base.} =
  raiseErrorBaseMethod("isZero(" & $x & ")")
method update*(x: Gvalue, y: int) {.base.} =
  raiseErrorBaseMethod("update(" & $x & "," & $y & ")")
method update*(x: Gvalue, y: float) {.base.} =
  raiseErrorBaseMethod("update(" & $x & "," & $y & ")")

method constLike*(x: Gvalue, value: int): Gvalue {.base.} =
  ## Convenience for scalar-like graph values.
  ## This is not a universal algebraic identity for every Gvalue subtype.
  result = x.newOneOf
  result.update value

proc treeReprImpl*(v: Gvalue,
                   children: proc(node: Gvalue): seq[Gvalue] {.closure.}): string =
  var shared = newseq[Gvalue]()
  var rendered = initNodeSet()

  proc render(x: Gvalue): seq[string] =
    let sharedIndex = shared.find x
    result = @[x.nodeRepr]
    if rendered.containsNode(x):
      if sharedIndex >= 0:
        result[0] &= " #" & $sharedIndex
      return
    if sharedIndex >= 0:
      result[0] &= " #" & $sharedIndex & "#"
    rendered.inclNode x
    for child in children(x):
      for renderedChild in render(child):
        result.add("  " & renderedChild)

  var seen = initNodeSet()
  proc collectShared(node: Gvalue) =
    if node == nil:
      return
    if not seen.markSeenNode(node):
      if shared.find(node) < 0:
        shared.add node
      return
    for input in children(node):
      collectShared(input)

  collectShared(v)

  result = render(v).join "\n"

proc treeRepr*(v: Gvalue): string =
  treeReprImpl(v, proc(node: Gvalue): seq[Gvalue] = node.inputs)

method `-`*(x: Gvalue): Gvalue {.base.} =
  raiseErrorBaseMethod("`-`(" & $x & ")")
method `+`*(x: Gvalue, y: Gvalue): Gvalue {.base.} =
  raiseErrorBaseMethod("`+`(" & $x & ", " & $y & ")")
method `*`*(x: Gvalue, y: Gvalue): Gvalue {.base.} =
  raiseErrorBaseMethod("`*`(" & $x & ", " & $y & ")")
method `-`*(x: Gvalue, y: Gvalue): Gvalue {.base.} =
  raiseErrorBaseMethod("`-`(" & $x & ", " & $y & ")")
method `/`*(x: Gvalue, y: Gvalue): Gvalue {.base.} =
  raiseErrorBaseMethod("`/`(" & $x & ", " & $y & ")")
method exp*(x: Gvalue): Gvalue {.base.} =
  raiseErrorBaseMethod("exp(" & $x & ")")

method walkHiddenDeps*(v: Gvalue,
                       mode: InputWalkMode,
                       visit: GnodeVisit) {.base.} =
  discard
method appendSignatureTokens*(v: Gvalue,
                              tokens: var seq[GradSigToken]) {.base.} =
  discard
