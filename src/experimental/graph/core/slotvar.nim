import base

## `slotVar(x)` is a transparent alias of `x`: the same value and derivatives,
## but a distinct node. It is the graph spelling of `let slot = x`. It has no
## `inputView` because its raw input is the whole dependency surface. It lets a
## backward hook differentiate a node replica with respect to one slot without
## reaching a sibling slot that holds the same node. See DESIGN.md section 5.

proc slotVarForward(v: Gvalue) =
  v.valCopy v.inputs[0]

proc slotVarBackward(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  rootedUpstream(zb, z)

let slotVarFunc = Gfunc(
  forward: slotVarForward,
  backward: slotVarBackward,
  name: "slotVar")

proc isSlotVarNode*(x: Gvalue): bool =
  x.gfunc == slotVarFunc

proc slotVar*[T: Gvalue](x: T): T =
  ## A fresh differentiation target holding x's value (see the module doc).
  # graphNode clears any static-zero marker when it installs the input edge.
  graphNode(T(x.newOneOf), @[Gvalue(x)], slotVarFunc, "slotVar")
