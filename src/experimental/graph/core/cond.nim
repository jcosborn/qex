import base

type
  Gcond = ref object of Gfunc

proc newCondNodeErased(selector: Gvalue,
                       whenTrue: Gvalue,
                       whenFalse: Gvalue): Gvalue
proc isCondNode*(v: Gvalue): bool

proc condInputView(v: Gvalue,
                   mode: InputWalkMode,
                   visit: GnodeVisit) =
  visit(v.inputs[0])
  case mode
  of iwmEval:
    visit(if v.inputs[0].isZero: v.inputs[2] else: v.inputs[1])
  of iwmReachable, iwmBackward:
    visit(v.inputs[1])
    visit(v.inputs[2])

proc condBackward(zb: Gvalue,
                  z: Gvalue,
                  i: int,
                  input: Gvalue): Gvalue =
  let upstream = rootedUpstream(zb, z)
  case i
  of 0:
    result = z.inputs[0].zeroLike
  of 1:
    result = newCondNodeErased(z.inputs[0], upstream, upstream.zeroLike)
  of 2:
    result = newCondNodeErased(z.inputs[0], upstream.zeroLike, upstream)
  else:
    raiseValueError("malformed cond backward input index: " & $i)

proc condf(v: Gvalue) =
  v.valCopy(if v.inputs[0].isZero: v.inputs[2] else: v.inputs[1])

# Eval follows only the selected branch; reachable keeps both branches visible
# for traversal and diagnostics.
let gcond = Gcond(
  forward: condf,
  inputView: condInputView,
  backward: condBackward,
  name: "cond")

proc isCondNode*(v: Gvalue): bool =
  v.gfunc of Gcond

proc condParts*(v: Gvalue): tuple[selector, whenTrue, whenFalse: Gvalue] =
  if not v.isCondNode:
    raiseValueError("expected cond node:\n" & v.nodeRepr)
  (v.inputs[0], v.inputs[1], v.inputs[2])

proc newCondNodeErased(selector: Gvalue,
                       whenTrue: Gvalue,
                       whenFalse: Gvalue): Gvalue =
  if not whenTrue.copyCompatible(whenFalse):
    raiseValueError(
      "cond branches must have compatible result types: " &
      "true-branch prototype cannot copy false-branch" &
      "\nprototype: " & whenTrue.nodeRepr &
      "\nbranch: " & whenFalse.nodeRepr)
  if not whenFalse.copyCompatible(whenTrue):
    raiseValueError(
      "cond branches must have compatible result types: " &
      "false-branch prototype cannot copy true-branch" &
      "\nprototype: " & whenFalse.nodeRepr &
      "\nbranch: " & whenTrue.nodeRepr)
  graphNode(
    whenTrue.newOneOf,
    @[selector, whenTrue, whenFalse],
    Gfunc(gcond),
    "cond")

proc newCondNode*[T: Gvalue](selector: Gvalue,
                             whenTrue: T,
                             whenFalse: T): T =
  when Gvalue is T:
    newCondNodeErased(selector, whenTrue, whenFalse)
  else:
    T(newCondNodeErased(selector, whenTrue, whenFalse))

proc distributeCond*(v: Gvalue,
                     branch: proc(b: Gvalue): Gvalue): Gvalue =
  ## Apply `branch` to both branches of cond node `v` and rebuild the cond;
  ## the selector contributes nothing (the a.e. piecewise convention). A nil
  ## branch result means "no contribution" and is filled with the other side's
  ## zeroLike; returns nil when both are nil.
  let parts = v.condParts
  let whenTrue = branch(parts.whenTrue)
  let whenFalse = branch(parts.whenFalse)
  if whenTrue == nil and whenFalse == nil:
    return nil
  newCondNodeErased(
    parts.selector,
    (if whenTrue == nil: whenFalse.zeroLike else: whenTrue),
    (if whenFalse == nil: whenTrue.zeroLike else: whenFalse))
