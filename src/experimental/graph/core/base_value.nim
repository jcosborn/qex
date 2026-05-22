from strutils import join
import std/sets
import base_types

method newOneOf*(x: Gvalue): Gvalue {.base.} =
  raiseErrorBaseMethod("newOneOf(" & $x & ")")
method valCopy*(z: Gvalue, x: Gvalue) {.base.} =
  raiseErrorBaseMethod("valCopy(" & $z & "," & $x & ")")
method copyCompatible*(prototype: Gvalue, value: Gvalue): bool {.base.} =
  false
method zeroLike*(x: Gvalue): Gvalue {.base.} = x.newOneOf
method oneLike*(x: Gvalue): Gvalue {.base.} =
  raiseErrorBaseMethod("oneLike(" & $x & ")")

method addLike*(prototype: Gvalue, x: Gvalue, y: Gvalue): Gvalue {.base.} =
  raiseErrorBaseMethod("addLike(" & $prototype & "," & $x & "," & $y & ")")

method isZero*(x: Gvalue): bool {.base.} =
  raiseErrorBaseMethod("isZero(" & $x & ")")

method supportsCondSelection*(x: Gvalue): bool {.base.} =
  false

proc treeReprImpl*(v: Gvalue,
                   children: proc(node: Gvalue): seq[Gvalue] {.closure.}): string =
  var shared = newseq[Gvalue]()
  var rendered = initHashSet[NodeKey]()

  proc render(x: Gvalue): seq[string] =
    let sharedIndex = shared.find x
    result = @[x.nodeRepr]
    if x.nodeKey in rendered:
      if sharedIndex >= 0:
        result[0] &= " #" & $sharedIndex
      return
    if sharedIndex >= 0:
      result[0] &= " #" & $sharedIndex & "#"
    rendered.incl x.nodeKey
    for child in children(x):
      for renderedChild in render(child):
        result.add("  " & renderedChild)

  var seen = initHashSet[NodeKey]()
  proc collectShared(node: Gvalue) =
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

method scaleLike*(contribution: Gvalue, upstream: Gvalue): Gvalue {.base.} =
  raiseErrorBaseMethod("scaleLike(" & $contribution & ", " & $upstream & ")")

method walkHiddenDeps*(v: Gvalue,
                       mode: InputWalkMode,
                       visit: GnodeVisit) {.base.} =
  discard
method appendSignatureTokens*(v: Gvalue,
                              tokens: var seq[GradSigToken]) {.base.} =
  discard
