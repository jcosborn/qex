import ../../core
import ../callable

type
  ApplyPartialView* = object
    node: Gvalue
    base*: Gvalue
    target*: Gvalue

proc requireApplyPartialBase(base: Gvalue,
                             label: string,
                             node: Gvalue = nil): Gvalue =
  if base == nil:
    raiseError(label & " base input cannot be nil" & node.nodeContext)
  base

proc requireApplyPartialTarget*(target: Gvalue): Gvalue =
  if target == nil:
    raiseValueError("apply partial target cannot be nil")
  target

proc needsConcreteApplyPartial*(target: Gvalue): bool =
  target != nil and not isCallableLike(target)

proc acceptsApplyPartialScale*(zb: Gvalue): bool =
  zb == nil or not zb.isCallableLike

proc isApplyPartialNode*(v: Gvalue): bool =
  v != nil and v.gfunc != nil and v.gfunc.deferredApplyKind == dakApplyPartial

proc requireApplyPartialView*(v: Gvalue,
                              label = "applyPartialDeferred node"): ApplyPartialView =
  v.requireInputCountExactly(2, label)
  result.node = v
  result.base = v.requireNodeInput(0, label, "base").requireApplyPartialBase(label, v)
  result.target = v.requireNodeInput(1, label, "target").requireApplyPartialTarget

proc newApplyPartialNode*(partialGfunc: Gfunc,
                          base: Gvalue,
                          target: Gvalue,
                          label: string,
                          node: Gvalue = nil): Gvalue =
  let checkedBase = base.requireApplyPartialBase(label, node)
  let checkedTarget = target.requireApplyPartialTarget
  graphNode(checkedTarget.newOneOf, @[checkedBase, checkedTarget], partialGfunc, label)

proc contributeApplyPartialTarget*(partialGfunc: Gfunc,
                                   zb: Gvalue,
                                   z: Gvalue,
                                   target: Gvalue,
                                   label: string): Gvalue =
  result = newApplyPartialNode(partialGfunc, z, target, label)
  if zb != nil:
    result = result.scaleLike zb

proc applyPartialBase*(v: Gvalue): Gvalue =
  v.requireApplyPartialView.base

proc applyPartialTarget*(v: Gvalue): Gvalue =
  v.requireApplyPartialView.target
