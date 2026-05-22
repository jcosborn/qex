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

proc requireCondSelector(selector: Gvalue) =
  if not selector.supportsCondSelection:
    raiseValueError(
      "cond condition must be a scalar or int selector, got:\n" &
      selector.nodeRepr)

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

proc condWalkEvalInputs(v: Gvalue,
                        visit: GnodeVisit) =
  let view = v.requireCondView
  visit(view.selector)
  visit(view.selectedBranch)

proc condWalkDependInputs(v: Gvalue,
                          visit: GnodeVisit) =
  let view = v.requireCondView
  visit(view.selector)
  visit(view.whenTrue)
  visit(view.whenFalse)

proc condWalkGradSignatureInputs(v: Gvalue,
                                 visit: GnodeVisit) =
  visit(v.requireCondView.selector)

proc condBackwardTarget(zb: Gvalue,
                        z: Gvalue,
                        target: Gvalue,
                        dep: Gvalue): Gvalue =
  discard dep
  let view = z.requireCondView
  result = newCondNode(
    view.selector,
    view.whenTrue.gradIsolated(target),
    view.whenFalse.gradIsolated(target))
  if zb != nil:
    result = result.scaleLike zb

proc condf(v: Gvalue) =
  v.valCopy v.requireCondView.selectedBranch

let gcond = newGfunc(
  forward = condf,
  # Eval follows only the selected branch; depend-mode keeps both branches
  # reachable; grad signatures recurse through the selector while raw branch
  # identities still enter the node signature.
  depWalks = GdepWalks(
    eval: condWalkEvalInputs,
    gradSignature: condWalkGradSignatureInputs,
    depend: condWalkDependInputs),
  backwardTarget = condBackwardTarget,
  name = "cond")

proc newCondNode*(selector: Gvalue,
                  whenTrue: Gvalue,
                  whenFalse: Gvalue): Gvalue =
  let checkedSelector = selector.requireGraphValue("cond condition")
  let checkedTrue = whenTrue.requireGraphValue("cond true-branch")
  let checkedFalse = whenFalse.requireGraphValue("cond false-branch")
  checkedSelector.requireCondSelector
  requireCompatibleCondBranches(checkedTrue, checkedFalse)
  graphNode(
    checkedTrue.newOneOf,
    @[checkedSelector, checkedTrue, checkedFalse],
    gcond,
    "cond")

proc newCondNode*[T: Gvalue](selector: Gvalue,
                             whenTrue: T,
                             whenFalse: T): T =
  T(newCondNode(selector, Gvalue(whenTrue), Gvalue(whenFalse)))
