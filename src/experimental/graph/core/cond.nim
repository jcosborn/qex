import base, grad_engine

type
  CondView = object
    selector: Gvalue
    whenTrue: Gvalue
    whenFalse: Gvalue

proc newCondNode*(selector: Gvalue,
                  whenTrue: Gvalue,
                  whenFalse: Gvalue): Gvalue
proc newCondNode*[T: Gvalue](selector: Gvalue,
                             whenTrue: T,
                             whenFalse: T): T

proc requireCondOperand(value: Gvalue,
                        label: string) =
  if value == nil:
    raiseValueError("cond " & label & " cannot be nil")

proc requireCondBranchCopyCompatible(protoSource: Gvalue,
                                     branch: Gvalue,
                                     label: string) =
  if not protoSource.copyCompatible(branch):
    raiseValueError(
      "cond branches must have compatible result types: " & label &
      "\nprototype: " & protoSource.nodeRepr &
      "\nbranch: " & branch.nodeRepr)

proc requireCompatibleCondBranches(x: Gvalue,
                                   y: Gvalue) =
  x.requireCondBranchCopyCompatible(y, "true-branch prototype cannot copy false-branch")
  y.requireCondBranchCopyCompatible(x, "false-branch prototype cannot copy true-branch")

proc requireCondView(v: Gvalue): CondView =
  v.requireInputCountExactly(3, "cond node")
  CondView(
    selector: v.requireNodeInput(0, "cond node", "condition"),
    whenTrue: v.requireNodeInput(1, "cond node", "true-branch"),
    whenFalse: v.requireNodeInput(2, "cond node", "false-branch"))

proc selectedBranch(view: CondView): Gvalue =
  if view.selector.isZero:
    return view.whenFalse
  view.whenTrue

proc branchExpr(view: CondView,
                whenTrue: Gvalue,
                whenFalse: Gvalue): Gvalue =
  newCondNode(view.selector, whenTrue, whenFalse)

proc newCondResult(whenTrue: Gvalue): Gvalue =
  whenTrue.newOneOf

proc newCondResult[T: Gvalue](whenTrue: T): T =
  T(whenTrue.newOneOf)

proc walkCondEvalInputs(view: CondView,
                        visit: GnodeVisit) =
  visit(view.selector)
  visit(view.selectedBranch)

proc walkCondDependInputs(view: CondView,
                          visit: GnodeVisit) =
  visit(view.selector)
  visit(view.whenTrue)
  visit(view.whenFalse)

proc condWalkEvalInputs(v: Gvalue,
                        visit: GnodeVisit) =
  v.requireCondView.walkCondEvalInputs(visit)

proc condWalkDependInputs(v: Gvalue,
                          visit: GnodeVisit) =
  v.requireCondView.walkCondDependInputs(visit)

proc condWalkGradSignatureInputs(v: Gvalue,
                                 visit: GnodeVisit) =
  visit(v.requireCondView.selector)

proc condBackwardTarget(zb: Gvalue,
                        z: Gvalue,
                        target: Gvalue,
                        dep: Gvalue): Gvalue =
  discard dep
  let view = z.requireCondView
  result = view.branchExpr(
    gradOrZero(view.whenTrue, target),
    gradOrZero(view.whenFalse, target))
  if zb != nil:
    result = result.scaleLike zb

proc condf(v: Gvalue) =
  v.valCopy v.requireCondView.selectedBranch

let gcond = newGfunc(
  forward = condf,
  depWalks = newDepWalks(
    eval = walkedInputs(condWalkEvalInputs),
    gradSignature = walkedInputs(condWalkGradSignatureInputs),
    depend = walkedInputs(condWalkDependInputs)),
  backwardTarget = condBackwardTarget,
  name = "cond")

proc newCondNode*(selector: Gvalue,
                  whenTrue: Gvalue,
                  whenFalse: Gvalue): Gvalue =
  selector.requireCondOperand("condition")
  whenTrue.requireCondOperand("true-branch")
  whenFalse.requireCondOperand("false-branch")
  requireCompatibleCondBranches(whenTrue, whenFalse)
  graphNode(newCondResult(whenTrue), @[selector, whenTrue, whenFalse], gcond, "cond")

proc newCondNode*[T: Gvalue](selector: Gvalue,
                             whenTrue: T,
                             whenFalse: T): T =
  selector.requireCondOperand("condition")
  whenTrue.requireCondOperand("true-branch")
  whenFalse.requireCondOperand("false-branch")
  requireCompatibleCondBranches(whenTrue, whenFalse)
  graphNode(newCondResult(whenTrue), @[selector, whenTrue, whenFalse], gcond, "cond")
