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
