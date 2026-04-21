import base, grad_engine

type
  CondView = object
    selector: Gvalue
    whenTrue: Gvalue
    whenFalse: Gvalue

proc cond*(c: Gvalue, x: Gvalue, y: Gvalue): Gvalue
proc newCondNode(selector: Gvalue,
                 whenTrue: Gvalue,
                 whenFalse: Gvalue): Gvalue

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
    result = zb * result

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

proc newCondNode(selector: Gvalue,
                 whenTrue: Gvalue,
                 whenFalse: Gvalue): Gvalue =
  graphNode(whenTrue.newOneOf, @[selector, whenTrue, whenFalse], gcond, "cond")

proc cond*(c: Gvalue, x: Gvalue, y: Gvalue): Gvalue =
  c.requireCondOperand("condition")
  x.requireCondOperand("true-branch")
  y.requireCondOperand("false-branch")
  requireCompatibleCondBranches(x, y)
  result = newCondNode(c, x, y)
