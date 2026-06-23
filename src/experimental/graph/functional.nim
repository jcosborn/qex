import std/[sets, tables]
import core
import scalar

# --- lambda values, normalization, and substitution ---


type
  LambdaRefKind = enum
    lrkLocal, lrkProduced
  GlambdaRef* {.final.} = ref object of Gvalue
    kind: LambdaRefKind
    paramProto: Gvalue
    resultProto: Gvalue
    binding: Gvalue
  Glambda* {.final.} = ref object of Gvalue
    param: Gvalue
    body: Gvalue
    captureParams: seq[Gvalue]
  Bindings = NodeTable[Gvalue]
  LambdaProtoPair = tuple[expected: NodeKey, actual: NodeKey]

proc isResolvedLambda(fn: Glambda): bool {.inline.} =
  fn != nil and fn.body != nil

proc isLambdaShaped(v: Gvalue): bool {.inline.} =
  v of Glambda or v of GlambdaRef

proc isLocalLambdaRef(w: GlambdaRef): bool {.inline.} =
  w != nil and w.kind == lrkLocal

proc isProducedLambdaRef(w: GlambdaRef): bool {.inline.} =
  w != nil and w.kind == lrkProduced

proc isGraphProducedLambdaRef(w: GlambdaRef): bool {.inline.} =
  w.isProducedLambdaRef and Gvalue(w).gfunc != nil

proc lambdaParamProtoOrNil(v: Gvalue): Gvalue =
  if v of Glambda:
    let fn = Glambda(v)
    if fn.isResolvedLambda:
      return fn.param
    return nil
  if v of GlambdaRef:
    return GlambdaRef(v).paramProto
  nil

proc lambdaResultProtoOrNil(v: Gvalue): Gvalue =
  if v of Glambda:
    let fn = Glambda(v)
    if fn.isResolvedLambda:
      return fn.body
    return nil
  if v of GlambdaRef:
    return GlambdaRef(v).resultProto
  nil

proc lambdaParamProto*(v: Gvalue): Gvalue =
  result = v.lambdaParamProtoOrNil
  if result == nil and v.isLambdaShaped:
    raiseValueError("lambda value is missing parameter prototype:\n" & v.nodeRepr)

proc lambdaResultProto*(v: Gvalue): Gvalue =
  result = v.lambdaResultProtoOrNil
  if result == nil and v.isLambdaShaped:
    raiseValueError("lambda value is missing result prototype:\n" & v.nodeRepr)

proc requireLambdaProto(v: Gvalue,
                        label: string): tuple[paramProto, resultProto: Gvalue] =
  result.paramProto = v.lambdaParamProto
  result.resultProto = v.lambdaResultProto
  if result.paramProto == nil or result.resultProto == nil:
    raiseValueError(label & " expects full lambda prototype:\n" & v.nodeRepr)

proc lambdaProtoMatch(expected: Gvalue,
                      actual: Gvalue,
                      seen: var HashSet[LambdaProtoPair]): bool =
  ## Structural prototype compatibility, walking nested lambda prototypes. The
  ## walk is symmetric (every leaf copyCompatible reachable from a lambda
  ## prototype is a type test), so one direction suffices wherever full
  ## equivalence was previously checked.
  if expected == nil or actual == nil:
    return expected == nil and actual == nil

  let pair = (expected: expected.nodeKey, actual: actual.nodeKey)
  if pair in seen:
    return true
  seen.incl pair

  let expectedParam = expected.lambdaParamProtoOrNil
  let expectedResult = expected.lambdaResultProtoOrNil
  let actualParam = actual.lambdaParamProtoOrNil
  let actualResult = actual.lambdaResultProtoOrNil
  if expected of Glambda or expected of GlambdaRef or
      actual of Glambda or actual of GlambdaRef or
      expectedParam != nil or expectedResult != nil or
      actualParam != nil or actualResult != nil:
    if expectedParam == nil or expectedResult == nil or
        actualParam == nil or actualResult == nil:
      return false
    return lambdaProtoMatch(expectedParam, actualParam, seen) and
      lambdaProtoMatch(expectedResult, actualResult, seen)

  expected.copyCompatible(actual)

proc lambdaProtoMatch(expected: Gvalue,
                      actual: Gvalue): bool =
  var seen = initHashSet[LambdaProtoPair]()
  lambdaProtoMatch(expected, actual, seen)

proc producedLambdaRef(paramPrototype: Gvalue,
                       resultPrototype: Gvalue): GlambdaRef =
  let grt = sharedGraphRuntime([paramPrototype, resultPrototype], "produced lambda ref")
  GlambdaRef(
    runtime: grt,
    kind: lrkProduced,
    paramProto: paramPrototype,
    resultProto: resultPrototype).assignStableNodeId

proc lambdaParam*(paramPrototype: Gvalue,
                  resultPrototype: Gvalue): GlambdaRef =
  let grt = sharedGraphRuntime([paramPrototype, resultPrototype], "lambda parameter")
  result = GlambdaRef(
    runtime: grt,
    kind: lrkLocal,
    paramProto: paramPrototype,
    resultProto: resultPrototype).assignStableNodeId
  result.updated

proc requireCompatibleLambdaRefBinding(z: GlambdaRef,
                                     x: Gvalue,
                                     label: string) =
  discard x.requireLambdaProto(label)
  if z.runtime != x.runtime:
    raiseValueError(label & " mixes multiple graph runtimes")
  if not Gvalue(z).lambdaProtoMatch(x):
    raiseValueError(
      label & " binding has incompatible lambda prototype" &
      "\nexpected: " & Gvalue(z).nodeRepr &
      "\nactual: " & x.nodeRepr)

proc hasFirstClassCotangent(v: Gvalue): bool =
  ## A value has a first-class cotangent iff it is not lambda-valued.
  if v.isLambdaShaped:
    discard v.requireLambdaProto("first-class cotangent")
    return false
  true

proc requireFirstClassCotangent(v: Gvalue) =
  ## Whole-lambda values have no first-class cotangent; differentiate applied
  ## lambdas with vjpOf/apply instead. Shared by zeroLike/addLike on lambdas.
  if not v.hasFirstClassCotangent:
    raiseValueError(
      "lambda-valued gradients are not first-class graph values; " &
      "use vjpOf/apply to build structural lambda VJPs instead" &
      "\nlambda: " & v.nodeRepr)

proc cotangentProto(primal: Gvalue): Gvalue =
  requireFirstClassCotangent(primal)
  primal.zeroLike

proc requireNoLambdaBodyBinding(x: Gvalue,
                                label: string) =
  if x of Glambda:
    raiseValueError(label & " does not accept lambda bodies:\n" & x.nodeRepr)
  if x of GlambdaRef and GlambdaRef(x).binding of Glambda:
    raiseValueError(label & " does not copy lambda bodies:\n" & x.nodeRepr)

proc bindLambdaRef(z: GlambdaRef,
                   x: Gvalue,
                   label: string) =
  z.requireCompatibleLambdaRefBinding(x, label)
  z.binding = x
  inc z.runtime.symbolicRevision
  Gvalue(z).updated

method newOneOf*(x: GlambdaRef): Gvalue =
  GlambdaRef(
    runtime: x.runtime,
    kind: x.kind,
    paramProto: x.paramProto,
    resultProto: x.resultProto).assignStableNodeId

method zeroLike*(x: GlambdaRef): Gvalue =
  cotangentProto(x)

method rootGradientSeed*(x: GlambdaRef): Gvalue =
  x.zeroLike

method isZero*(x: GlambdaRef): bool =
  x.gfunc == nil and x.inputs.len == 0 and x.binding == nil

proc zeroLambdaAdd(prototype: Gvalue, x: Gvalue, y: Gvalue): Gvalue =
  requireFirstClassCotangent(prototype)
  if x.isZero:
    return y
  if y.isZero:
    return x
  raiseValueError(
    "cannot accumulate non-zero lambda graph gradients" &
    "\nleft: " & x.nodeRepr &
    "\nright: " & y.nodeRepr)

method addLike*(prototype: GlambdaRef, x: Gvalue, y: Gvalue): Gvalue =
  zeroLambdaAdd(prototype, x, y)

method valCopy*(z: GlambdaRef, x: Gvalue) =
  if z.isGraphProducedLambdaRef:
    # Structural lambda-captures copy: refresh freshness only, no rebinding and
    # no symbolic-revision bump.
    z.requireCompatibleLambdaRefBinding(x, "lambda ref graph valCopy")
    z.updated
    return
  x.requireNoLambdaBodyBinding("lambda ref valCopy")
  z.bindLambdaRef(x, "lambda ref valCopy")

method copyCompatible*(prototype: GlambdaRef, value: Gvalue): bool =
  lambdaProtoMatch(prototype, value)

method `$`*(x: GlambdaRef): string =
  let label =
    case x.kind
    of lrkLocal:
      "local"
    of lrkProduced:
      "lambda"
  if x.binding == nil:
    return label
  label & "(bound)"

method newOneOf*(x: Glambda): Gvalue =
  if x.isResolvedLambda:
    return producedLambdaRef(x.param, x.body)
  Glambda(runtime: x.runtime, param: x.param).assignStableNodeId

method zeroLike*(x: Glambda): Gvalue =
  cotangentProto(x)

method rootGradientSeed*(x: Glambda): Gvalue =
  x.zeroLike

method copyCompatible*(prototype: Glambda, value: Gvalue): bool =
  if not prototype.isResolvedLambda:
    return false
  lambdaProtoMatch(prototype, value)

method addLike*(prototype: Glambda, x: Gvalue, y: Gvalue): Gvalue =
  zeroLambdaAdd(prototype, x, y)

proc grad*(dep: Gvalue, x: Glambda | GlambdaRef): Gvalue =
  raiseValueError(
    "grad with respect to a whole lambda value is not first-class; " &
    "differentiate an applied lambda with vjpOf/apply or differentiate scalar captures" &
    "\nlambda: " & Gvalue(x).nodeRepr)

proc cond*[C: Gscalar | Gint](c: C, x: Glambda, y: Glambda): Gvalue =
  newCondNode(c, Gvalue(x), Gvalue(y))

method `$`*(x: Glambda): string =
  result = "lambda(" & $x.param & " -> " & $x.body & ")"
  if x.captureParams.len > 0:
    result &= "[captureParams:" & $x.captureParams.len & "]"

proc structuralLambdaGraphDeps(v: Gvalue,
                               mode: InputWalkMode): seq[Gvalue] =
  ## Collect traversal-visible deps for a structural lambda-valued graph node.
  ## LambdaRef self-bindings are excluded here; including them would make stale
  ## lambda captures look current.
  var deps: seq[Gvalue] = @[]
  v.walkInputView(mode, proc(input: Gvalue) =
    if input.nodeKey != v.nodeKey:
      deps.add input)
  deps

proc freshLambdaBinding(v: Gvalue): Gvalue =
  ## Follow only symbolic lambda-ref bindings. Graph-produced nodes keep their
  ## lambda structure in ordinary graph inputs such as apply/cond branches.
  if not (v of GlambdaRef):
    return nil
  let w = GlambdaRef(v)
  if not w.isProducedLambdaRef:
    return nil
  w.binding

proc pushReverse(stack: var seq[Gvalue], deps: openArray[Gvalue]) =
  for i in countdown(deps.len - 1, 0):
    stack.add deps[i]

proc walkLambdaGraph(roots: openArray[Gvalue],
                     excludeValues: openArray[Gvalue],
                     mode: InputWalkMode,
                     visitValue: GnodeVisit = nil,
                     visitNode: GnodeVisit = nil,
                     stopValues: openArray[Gvalue] = []) =
  ## Walk graph structure until lambda boundaries, then follow only symbolic
  ## structural lambda bindings.
  var seenLambda = initHashSet[NodeKey]()
  var seenValues = initHashSet[NodeKey]()
  var excludedKeys = initHashSet[NodeKey]()
  var stopKeys = initHashSet[NodeKey]()
  for value in excludeValues:
    if value != nil:
      excludedKeys.incl value.nodeKey
  for value in stopValues:
    if value != nil:
      stopKeys.incl value.nodeKey
  var stack: seq[Gvalue] = @[]
  stack.pushReverse(roots)

  while stack.len > 0:
    let node = stack[^1]
    stack.setLen(stack.len - 1)
    if node.nodeKey in stopKeys:
      continue

    if node of Glambda:
      let fn = Glambda(node)
      if fn.isResolvedLambda:
        if not seenLambda.markSeenNode(node):
          continue
        if visitNode != nil:
          visitNode node
        stack.pushReverse(fn.inputs)
        continue

    let w =
      if node of GlambdaRef:
        GlambdaRef(node)
      else:
        nil
    if w != nil:
      if not seenLambda.markSeenNode(node):
        continue
      if visitNode != nil:
        visitNode node
      var deps: seq[Gvalue] = @[]
      if w.isLocalLambdaRef:
        if mode == iwmEval and w.binding != nil:
          stack.add w.binding
          continue
        if mode == iwmReachable:
          if w.binding != nil:
            deps.add w.binding
        else:
          deps = node.collectInputView(mode)
      else:
        let lambdaBinding = node.freshLambdaBinding
        if lambdaBinding != nil:
          stack.add lambdaBinding
          continue
        deps = node.structuralLambdaGraphDeps(mode)
      stack.pushReverse(deps)
      continue

    if not seenValues.markSeenNode(node):
      continue
    if visitNode != nil:
      visitNode node
    if node.nodeKey notin excludedKeys:
      if visitValue != nil:
        visitValue node
    let deps = node.collectInputView(mode)
    stack.pushReverse(deps)

proc collectLambdaValueDeps(roots: openArray[Gvalue],
                            seedValues: openArray[Gvalue],
                            mode: InputWalkMode = iwmReachable): seq[Gvalue] =
  ## Return ordinary value dependencies under lambda boundaries. Seed values are
  ## excluded from the result; leaf lambda refs expose their symbolic binding,
  ## and structural lambda-valued graph refs expose traversal-visible graph deps.
  var deps: seq[Gvalue] = @[]
  walkLambdaGraph(roots, seedValues, mode, proc(value: Gvalue) =
    deps.add value)
  deps

type
  CloneCtx = object
    subst: Bindings
    memo: Bindings
    preserve: NodeSet
    extraInputs: NodeTable[seq[Gvalue]]

proc clone(v: Gvalue, ctx: var CloneCtx): Gvalue

proc clonedGraphInputs(src: Gvalue,
                       ctx: var CloneCtx): seq[Gvalue] =
  result = newseq[Gvalue](src.inputs.len)
  for i in 0..<src.inputs.len:
    result[i] = clone(src.inputs[i], ctx)
  if ctx.extraInputs.hasKey(src.nodeKey):
    for input in ctx.extraInputs[src.nodeKey]:
      let clonedInput = clone(input, ctx)
      var found = false
      for value in result:
        if value.nodeKey == clonedInput.nodeKey:
          found = true
          break
      if not found:
        result.add clonedInput

proc withoutLambdaBindings(subst: Bindings, fn: Glambda): Bindings =
  result = subst
  if result.hasKey(fn.param.nodeKey):
    result.del(fn.param.nodeKey)
  for binding in fn.captureParams:
    if result.hasKey(binding.nodeKey):
      result.del(binding.nodeKey)

proc cloneLambdaBinder(v: Gvalue, ctx: var CloneCtx): Gvalue =
  # Every GlambdaRef is constructed with non-nil param/result prototypes
  # (producedLambdaRef / lambdaParam / newOneOf all enforce it), so a binder's
  # prototypes clone unconditionally.
  result =
    if v of GlambdaRef:
      let w = GlambdaRef(v)
      Gvalue(GlambdaRef(
        runtime: w.runtime,
        kind: w.kind,
        paramProto: clone(w.paramProto, ctx),
        resultProto: clone(w.resultProto, ctx)).assignStableNodeId)
    else:
      v.newOneOf
  result.updated
  ctx.memo[v.nodeKey] = result
  ctx.subst[v.nodeKey] = result

proc cloneResolvedLambda(fn: Glambda,
                         ctx: var CloneCtx): Gvalue =
  if not fn.isResolvedLambda:
    raiseValueError("lambda clone expects resolved lambda")
  result = Glambda(runtime: fn.runtime).assignStableNodeId
  ctx.memo[fn.nodeKey] = result

  let r = Glambda(result)
  # Inner clone forks the substitution (binder params shadow) but shares the memo
  # so already-cloned nodes stay identity-stable.
  var innerCtx = CloneCtx(
    subst: ctx.subst.withoutLambdaBindings(fn),
    memo: ctx.memo,
    preserve: ctx.preserve,
    extraInputs: ctx.extraInputs)
  r.param = cloneLambdaBinder(fn.param, innerCtx)
  r.captureParams = newseq[Gvalue](fn.captureParams.len)
  for i in 0..<fn.captureParams.len:
    r.captureParams[i] = cloneLambdaBinder(fn.captureParams[i], innerCtx)
  r.body = clone(fn.body, innerCtx)
  r.inputs = fn.clonedGraphInputs(ctx)
  r.gfunc = fn.gfunc

proc cloneWithFreshMemo(v: Gvalue,
                        subst: Bindings,
                        preserve: openArray[Gvalue] = [],
                        extraInputs: NodeTable[seq[Gvalue]] = initTable[NodeKey, seq[Gvalue]]()): Gvalue =
  var preserveSet = initHashSet[NodeKey]()
  for value in preserve:
    if value != nil:
      preserveSet.incl value.nodeKey
  var ctx = CloneCtx(
    subst: subst,
    memo: initTable[NodeKey, Gvalue](),
    preserve: preserveSet,
    extraInputs: extraInputs)
  clone(v, ctx)

proc clone(v: Gvalue, ctx: var CloneCtx): Gvalue =
  let key = v.nodeKey
  if key in ctx.preserve:
    return v
  if ctx.subst.hasKey(key):
    return ctx.subst[key]
  if ctx.memo.hasKey(key):
    return ctx.memo[key]

  if v of Glambda:
    let fn = Glambda(v)
    if fn.isResolvedLambda:
      return fn.cloneResolvedLambda(ctx)

  let lambdaBinding = v.freshLambdaBinding
  if lambdaBinding != nil:
    result = clone(lambdaBinding, ctx)
    ctx.memo[key] = result
    return result

  if v.gfunc == nil and v.inputs.len == 0:
    return v

  # Generic graph-node clone is valid only for nodes whose structural state is
  # fully represented by `inputs` plus `gfunc`. Nodes with extra structural
  # fields need explicit clone handling or prototype state preserved by
  # `newOneOf`.
  result = v.newOneOf
  ctx.memo[key] = result
  result.inputs = v.clonedGraphInputs(ctx)
  result.gfunc = v.gfunc

let lambdaCaptureGraph = Gfunc(
  forward: proc(v: Gvalue) =
    discard v,
  name: "lambda captures")

proc initLambdaSubst(fn: Glambda,
                     arg: Gvalue = nil,
                     extraBindings: Bindings = initTable[NodeKey, Gvalue]()): Bindings =
  if not fn.isResolvedLambda:
    raiseValueError("lambda substitution expects resolved lambda")
  var bindings = extraBindings
  for i in 0..<fn.captureParams.len:
    var captureValue = fn.inputs[i]
    if bindings.hasKey(captureValue.nodeKey):
      captureValue = bindings[captureValue.nodeKey]
    bindings[fn.captureParams[i].nodeKey] = captureValue
  if arg != nil:
    bindings[fn.param.nodeKey] = arg
  result = bindings

proc bindLambdaCaptureParams(fn: Glambda,
                             captures: seq[Gvalue],
                             preserve: openArray[Gvalue] = [],
                             extraInputs: NodeTable[seq[Gvalue]] = initTable[NodeKey, seq[Gvalue]]()) =
  if captures.len == 0:
    fn.captureParams.setLen(0)
    fn.inputs.setLen(0)
    fn.gfunc = nil
    return

  fn.captureParams = newseq[Gvalue](captures.len)
  fn.inputs = newseq[Gvalue](captures.len)
  var subst = initTable[NodeKey, Gvalue]()
  for i in 0..<captures.len:
    let param = captures[i].newOneOf
    param.updated
    fn.captureParams[i] = param
    fn.inputs[i] = captures[i]
    subst[captures[i].nodeKey] = param
  fn.gfunc = lambdaCaptureGraph
  fn.body = fn.body.cloneWithFreshMemo(
    subst,
    preserve,
    extraInputs)

proc normalizeLambda(fn: Glambda,
                                    preserve: openArray[Gvalue] = [],
                                    extraInputs: NodeTable[seq[Gvalue]] = initTable[NodeKey, seq[Gvalue]]()) =
  ## Normalize a resolved lambda into capture-normal form.
  ## After this pass, captures live in `captureParams` and the body can be instantiated by
  ## ordinary substitution without having to rediscover free variables.
  if not fn.isResolvedLambda:
    return

  if fn.captureParams.len > 0:
    fn.body = fn.body.cloneWithFreshMemo(
      fn.initLambdaSubst(),
      preserve,
      extraInputs)
  fn.captureParams.setLen(0)
  fn.inputs.setLen(0)
  fn.gfunc = nil
  var stops = newseq[Gvalue](preserve.len + 1)
  stops[0] = fn.param
  for i in 0..<preserve.len:
    stops[i + 1] = preserve[i]
  var seenCaps = initHashSet[NodeKey]()
  var captures: seq[Gvalue] = @[]
  proc addCapture(value: Gvalue) =
    if value.gfunc != nil or value.inputs.len != 0 or value of Glambda:
      return
    if seenCaps.markSeenNode(value):
      captures.add value

  walkLambdaGraph(
    [fn.body],
    [],
    iwmReachable,
    visitValue = addCapture,
    visitNode = addCapture,
    stopValues = stops)
  fn.bindLambdaCaptureParams(captures, preserve, extraInputs)

proc instantiateNormalizedBody(fn: Glambda,
                               x: Gvalue,
                               extraBindings: Bindings = initTable[NodeKey, Gvalue](),
                               extraInputs: NodeTable[seq[Gvalue]] = initTable[NodeKey, seq[Gvalue]]()): Gvalue =
  result = fn.body.cloneWithFreshMemo(
    fn.initLambdaSubst(x, extraBindings),
    [Gvalue(fn)],
    extraInputs)
  if result of Glambda and Glambda(result).isResolvedLambda:
    normalizeLambda(Glambda(result), [Gvalue(fn)], extraInputs)

proc lambda*(param: Gvalue, body: Gvalue): Glambda =
  # Lambda captures must stay in one runtime for cloning and cache identity.
  let grt = sharedGraphRuntime([param, body], "lambda")
  result = Glambda(runtime: grt, param: param, body: body).assignStableNodeId
  normalizeLambda(result)
  result.updated

# --- apply, symbolic VJPs, and higher-order AD ---


type
  LambdaVjpTargetKind = enum
    lvtkArgument, lvtkValue
  LambdaVjpKind = enum
    lvkCall, lvkResult
  LambdaVjpTarget = object
    case kind: LambdaVjpTargetKind
    of lvtkArgument:
      discard
    of lvtkValue:
      value: Gvalue
  ApplyBackwardDep = object
    input: Gvalue
    targets: seq[LambdaVjpTarget]
  ApplyBackwardDeps = seq[ApplyBackwardDep]
  Gapply = ref object of Gfunc
  GvjpOf = ref object of Gfunc
    kind: LambdaVjpKind
  VjpSpec = object
    ## Both VJP kinds carry the same two prototypes: the lambda argument shape
    ## and the (placeholder) body shape. `kind` selects how they are used, not
    ## which fields exist. `callVjpSpec` requires `argProto`; `resultVjpSpec`
    ## requires `bodyProto`.
    kind: LambdaVjpKind
    target: LambdaVjpTarget
    argProto: Gvalue
    bodyProto: Gvalue
  VjpKey = object
    # The key omits the resolved lambda: directLambda(fun) is determined by funId,
    # so selectedId would never separate two equal-funId keys. bodyProto is omitted
    # for the same reason: it is fixed by funId + target + argProto, so it adds no
    # distinguishing power.
    funId: NodeId
    targetKind: LambdaVjpTargetKind
    targetId: NodeId
    argProtoId: NodeId
  ActiveVjp = object
    kind: LambdaVjpKind
    owner: Glambda
    target: LambdaVjpTarget
    bodyProto: Gvalue
    argProto: Gvalue
    shell: Glambda
  ApplyVjpBuildCtx = ref object
    ## Build-local context threaded through apply/VJP construction. It is never
    ## process-global: independent VJP builds must not share `active`/`memo`.
    ##
    ## `active` also drives the *value-target injection*, the one part of apply
    ## backward that is not "raw inputs plus ordinary backward". A generated apply
    ## that must expose a captured scalar cotangent target appends that target as
    ## an extra raw input (indices >= `FirstApplyValueTargetInput`) so backward
    ## traversal can reach it without reducing the apply node. The pieces:
    ##   - `FirstApplyValueTargetInput` / `extraApplyValueTargets`: the input convention;
    ##   - `activeApplyValueTargetInputs`: maps active value targets onto apply nodes;
    ##   - apply-owned VJP clone paths append them during body instantiation.
    memo: Table[VjpKey, Gvalue]
    active: seq[ActiveVjp]

const
  ApplyFunInput = 0
  ApplyArgInput = 1
  # Apply raw inputs are `[fun, arg, active value targets...]`. The extra value
  # targets are injected while cloning active structural VJP graphs so backward
  # traversal can see captured scalar targets without reducing the apply node.
  FirstApplyValueTargetInput = 2

proc newApplyVjpBuildCtx(): ApplyVjpBuildCtx =
  ApplyVjpBuildCtx(memo: initTable[VjpKey, Gvalue]())

proc makeVjpKey(fun: Gvalue,
                spec: VjpSpec): VjpKey =
  result.funId = fun.stableNodeId
  result.targetKind = spec.target.kind
  case spec.target.kind
  of lvtkArgument:
    result.targetId = 0
  of lvtkValue:
    result.targetId = spec.target.value.stableNodeId
  let argProto = spec.argProto
  result.argProtoId =
    if argProto == nil:
      0
    else:
      argProto.stableNodeId

proc sameVjpTarget(a: LambdaVjpTarget, b: LambdaVjpTarget): bool =
  if a.kind != b.kind:
    return false
  case a.kind
  of lvtkArgument:
    true
  of lvtkValue:
    a.value.nodeKey == b.value.nodeKey

# --- apply result prototype materialization ---

proc materializeApplyResultProto(proto: Gvalue,
                                 seen: var NodeSet): Gvalue =
  if not seen.markSeenNode(proto):
    raiseValueError(
      "apply result prototype contains a cycle:\n" &
      proto.nodeRepr)
  let nestedProto = proto.lambdaResultProto
  if nestedProto != nil:
    let paramProto = proto.lambdaParamProto
    return producedLambdaRef(
      paramProto,
      materializeApplyResultProto(nestedProto, seen))
  proto.newOneOf

proc applyResultProto(fun: Gvalue): Gvalue =
  let resultProto = fun.lambdaResultProto
  if resultProto == nil:
    return nil
  var seen = initHashSet[NodeKey]()
  materializeApplyResultProto(resultProto, seen)

# --- SVJP type algebra (specs + prototype derivations) ---

proc structuralVjpResultProto(bodyProto: Gvalue, gradProto: Gvalue): Gvalue =
  let resultProto = bodyProto.lambdaResultProto
  if resultProto == nil:
    producedLambdaRef(cotangentProto(bodyProto), gradProto)
  else:
    let paramProto = bodyProto.lambdaParamProto
    producedLambdaRef(paramProto, structuralVjpResultProto(resultProto, gradProto))

proc callVjpNodeProto(fun: Gvalue, gradProto: Gvalue): Gvalue =
  ## Result type of a call VJP over `fun : A -> B`:  `A -> SVJP(B, gradProto)`.
  ## Shared by the symbolic `vjpOf` node type, the function-position VJP, and the
  ## type threaded through a lambda-valued apply.
  let proto = fun.requireLambdaProto("lambda VJP prototype")
  producedLambdaRef(
    proto.paramProto,
    structuralVjpResultProto(proto.resultProto, gradProto))

proc targetGradProto(fun: Gvalue,
                     target: LambdaVjpTarget,
                     argumentProto: Gvalue): Gvalue =
  case target.kind
  of lvtkArgument:
    if argumentProto == nil:
      raiseValueError(
        "lambda VJP prototype requires argument prototype:\n" &
        fun.nodeRepr)
    argumentProto.zeroLike
  of lvtkValue:
    target.value.zeroLike

proc callVjpSpec(target: LambdaVjpTarget,
                 argProto: Gvalue,
                 bodyProto: Gvalue = nil): VjpSpec =
  if argProto == nil:
    raiseValueError("lambda call VJP requires argument prototype")
  VjpSpec(
    kind: lvkCall,
    target: target,
    argProto: argProto,
    bodyProto: bodyProto)

proc resultVjpSpec(target: LambdaVjpTarget,
                   bodyProto: Gvalue,
                   argProto: Gvalue = nil): VjpSpec =
  if bodyProto == nil:
    raiseValueError("lambda result VJP requires result prototype")
  VjpSpec(
    kind: lvkResult,
    target: target,
    argProto: argProto,
    bodyProto: bodyProto)

# Structural VJP result types follow:
#
#   SVJP(R, G) =
#     Cot(R) -> G      when R is an ordinary graph value
#     A -> SVJP(B, G)  when R is a lambda A -> B
#
# The helpers below derive the public result type of a symbolic `vjpOf`, the
# placeholder body type of an active recursive shell, and the type threaded
# through a lambda-valued apply.

proc vjpNodeProto(fun: Gvalue, spec: VjpSpec): Gvalue =
  ## Erased result type for a symbolic `vjpOf` node over `fun`.
  case spec.kind
  of lvkCall:
    callVjpNodeProto(fun, targetGradProto(fun, spec.target, fun.lambdaParamProto))
  of lvkResult:
    let proto = fun.requireLambdaProto("lambda VJP prototype")
    producedLambdaRef(proto.paramProto, spec.bodyProto)

proc vjpShellBodyProto(fun: Gvalue, spec: VjpSpec): Gvalue =
  ## Placeholder body type for a generated VJP shell before its body is built.
  case spec.kind
  of lvkCall:
    if spec.bodyProto != nil:
      spec.bodyProto
    else:
      let proto = fun.requireLambdaProto("lambda VJP prototype")
      let gradProto = targetGradProto(fun, spec.target, proto.paramProto)
      producedLambdaRef(gradProto, gradProto)
  of lvkResult:
    spec.bodyProto

proc vjpApplyResultProto(fun: Gvalue, spec: VjpSpec): Gvalue =
  ## Result type when threading a VJP through a lambda-valued apply node.
  case spec.kind
  of lvkCall:
    callVjpNodeProto(fun, targetGradProto(fun, spec.target, spec.argProto))
  of lvkResult:
    spec.bodyProto

# --- symbolic vjpOf node ---

proc vjpOfTarget(v: Gvalue): LambdaVjpTarget =
  case v.inputs.len
  of 1:
    LambdaVjpTarget(kind: lvtkArgument)
  of 2:
    LambdaVjpTarget(kind: lvtkValue, value: v.inputs[1])
  else:
    raiseValueError("malformed vjpOf node:\n" & v.nodeRepr)

proc vjpOfForward(v: Gvalue) =
  v.updated

let gvjpOfCall = GvjpOf(
  kind: lvkCall,
  forward: vjpOfForward,
  name: "vjpOf")
let gvjpOfResult = GvjpOf(
  kind: lvkResult,
  forward: vjpOfForward,
  name: "vjpOfResult")

proc newVjpOfNode(fun: Gvalue,
                  spec: VjpSpec): Gvalue =
  let proto = vjpNodeProto(fun, spec)
  var inputs = @[fun]
  if spec.target.kind == lvtkValue:
    inputs.add spec.target.value
  let gfunc =
    case spec.kind
    of lvkCall:
      gvjpOfCall
    of lvkResult:
      gvjpOfResult
  graphNode(proto.newOneOf, inputs, gfunc, "vjpOf")

# --- apply backward-dependency discovery ---

proc addUniqueNode(dst: var seq[Gvalue],
                   seen: var HashSet[NodeKey],
                   value: Gvalue) =
  if seen.markSeenNode(value):
    dst.add value

proc addUniqueTarget(dst: var seq[LambdaVjpTarget],
                     target: LambdaVjpTarget) =
  for existing in dst:
    if sameVjpTarget(existing, target):
      return
  dst.add target

proc addBackwardTarget(deps: var ApplyBackwardDeps,
                       input: Gvalue,
                       target: LambdaVjpTarget) =
  for i in 0..<deps.len:
    if deps[i].input.nodeKey == input.nodeKey:
      deps[i].targets.addUniqueTarget(target)
      return
  deps.add ApplyBackwardDep(
    input: input,
    targets: @[target])

proc activeVjpValueTargets(ctx: ApplyVjpBuildCtx): seq[Gvalue] =
  if ctx == nil:
    return
  var seen = initHashSet[NodeKey]()
  for active in ctx.active:
    if active.target.kind == lvtkValue and active.target.value != nil:
      result.addUniqueNode(seen, active.target.value)

proc activeApplyValueTargetInputs(ctx: ApplyVjpBuildCtx,
                                  root: Gvalue): NodeTable[seq[Gvalue]] =
  ## Apply raw inputs after index 1 are explicit value targets injected while
  ## cloning structural VJP graphs. Every apply node reached gets all active
  ## value targets: clone dedups them against fun/arg, and a target not reachable
  ## through an apply's captures contributes a zero gradient.
  var extraInputs = initTable[NodeKey, seq[Gvalue]]()
  let targets = activeVjpValueTargets(ctx)
  if targets.len == 0 or root == nil:
    return extraInputs
  walkLambdaGraph([root], [], iwmReachable, visitNode = proc(node: Gvalue) =
    if node.gfunc of Gapply:
      extraInputs[node.nodeKey] = targets)
  extraInputs

proc instantiateBodyForVjp(fn: Glambda,
                           arg: Gvalue,
                           ctx: ApplyVjpBuildCtx): Gvalue =
  instantiateNormalizedBody(
    fn,
    arg,
    extraInputs = activeApplyValueTargetInputs(ctx, fn.body))

proc normalizeLambdaForVjp(fn: Glambda,
                           preserve: openArray[Gvalue],
                           ctx: ApplyVjpBuildCtx,
                           body: Gvalue) =
  normalizeLambda(
    fn,
    preserve,
    activeApplyValueTargetInputs(ctx, body))

proc backwardReachableNodes(root: Gvalue): HashSet[NodeKey] =
  var owned = initHashSet[NodeKey]()
  proc visit(node: Gvalue) =
    if not owned.markSeenNode(node):
      return
    if node.gfunc == nil or node.gfunc.backward == nil:
      return
    node.walkInputView(iwmBackward, proc(child: Gvalue) =
      visit child
    )

  visit root
  owned

proc addMaximalBackwardTargets(dst: var seq[Gvalue],
                               seen: var HashSet[NodeKey],
                               deps: openArray[Gvalue]) =
  ## Keep only targets that are not already exposed through another target's
  ## backward surface.
  var ownerNodes: seq[HashSet[NodeKey]] = @[]
  for owner in deps:
    ownerNodes.add backwardReachableNodes(owner)

  for depIndex, dep in deps:
    var ownedByAnotherDep = false
    for ownerIndex in 0..<deps.len:
      if ownerIndex != depIndex and deps[ownerIndex].nodeKey != dep.nodeKey and
          dep.nodeKey in ownerNodes[ownerIndex]:
        ownedByAnotherDep = true
        break
    if not ownedByAnotherDep:
      dst.addUniqueNode(seen, dep)

proc extraApplyValueTargets(v: Gvalue): seq[Gvalue] =
  for i in FirstApplyValueTargetInput..<v.inputs.len:
    result.add v.inputs[i]

proc collectApplyValueTargets(fun: Gvalue,
                              arg: Gvalue): seq[Gvalue] =
  var seenTargets = initHashSet[NodeKey]()
  for target in collectLambdaValueDeps([fun], [fun]):
    result.addUniqueNode(seenTargets, target)

  if not arg.hasFirstClassCotangent:
    for target in collectLambdaValueDeps([arg], [arg]):
      result.addUniqueNode(seenTargets, target)

proc applyBackwardDeps(v: Gvalue): ApplyBackwardDeps =
  ## Apply raw inputs are `[fun, arg]`. Backward deps are the ordinary argument
  ## when it has a first-class cotangent, plus value targets discovered through
  ## lambda captures and active VJP targets appended as explicit raw inputs.
  let fun = v.inputs[ApplyFunInput]
  let arg = v.inputs[ApplyArgInput]

  let valueTargets = collectApplyValueTargets(fun, arg)
  if arg.hasFirstClassCotangent:
    result.addBackwardTarget(
      arg,
      LambdaVjpTarget(kind: lvtkArgument))

  var seen = initHashSet[NodeKey]()
  seen.incl fun.nodeKey
  var frontier: seq[Gvalue] = @[]
  frontier.addMaximalBackwardTargets(seen, valueTargets)
  for input in frontier:
    result.addBackwardTarget(
      input,
      LambdaVjpTarget(kind: lvtkValue, value: input))

  for input in v.extraApplyValueTargets:
    result.addBackwardTarget(
      input,
      LambdaVjpTarget(kind: lvtkValue, value: input))

proc resetApplyCache*(grt: GraphRuntime, stats = true) =
  grt.functional.applyCacheByNode = initTable[NodeId, ApplyCacheEntry]()
  if stats:
    grt.functional.applyCacheStats = ApplyCacheStats()

proc ensureInstantiation(v: Gvalue,
                         arg: Gvalue,
                         fn: Glambda): ApplyCacheEntry =
  let grt = v.runtime
  let key = ApplyCacheKey(
    selectedLambdaId: fn.stableNodeId,
    revision: grt.symbolicRevision)

  let nodeId = v.stableNodeId
  if not grt.functional.applyCacheByNode.hasKey(nodeId):
    grt.functional.applyCacheByNode[nodeId] = ApplyCacheEntry()
  result = grt.functional.applyCacheByNode[nodeId]
  if result.instantiated != nil and result.key == key:
    grt.functional.applyCacheStats.instantiationHits.inc
    return
  grt.functional.applyCacheStats.instantiationMisses.inc

  let instantiated = instantiateNormalizedBody(fn, arg)
  result.instantiated = instantiated
  result.key = key

proc apply*(fun: Gvalue, x: Gvalue): Gvalue

# --- lambda resolution walkers ---

proc followLambdaRefs(start: Gvalue, label: string): Gvalue =
  ## Follow GlambdaRef.binding chains (cycle-checked). Stops at the first
  ## non-ref node or unbound ref. Does NOT step cond/apply/vjpOf.
  result = start
  var seen = initHashSet[NodeKey]()
  while result of GlambdaRef and GlambdaRef(result).binding != nil:
    if not seen.markSeenNode(result):
      raiseValueError(label & " cyclic lambda binding:\n" & result.nodeRepr)
    result = GlambdaRef(result).binding

proc directLambda(fun: Gvalue, label: string): Glambda =
  let resolved = followLambdaRefs(fun, label)
  if resolved of Glambda:
    let fn = Glambda(resolved)
    if fn.isResolvedLambda:
      return fn
  nil

proc evaluatedLambda(fun: Gvalue,
                      label: string): Glambda
proc evaluatedVjpOfLambda(v: Gvalue,
                          label: string): Gvalue

proc nextEvalLambdaExpr(current: Gvalue,
                        label: string,
                        vjpLabel: string): Gvalue =
  if current.isCondNode and current.lambdaResultProto != nil:
    let parts = current.condParts
    discard parts.selector.eval
    return if parts.selector.isZero: parts.whenFalse else: parts.whenTrue

  if current.gfunc of Gapply and current.lambdaResultProto != nil:
    let fn = evaluatedLambda(
      current.inputs[ApplyFunInput],
      label & " applied function")
    return current.ensureInstantiation(current.inputs[ApplyArgInput], fn).instantiated

  if current.gfunc of GvjpOf and current.lambdaResultProto != nil:
    return current.evaluatedVjpOfLambda(label & " " & vjpLabel)

  if current of GlambdaRef:
    return GlambdaRef(current).binding
  nil

proc evaluatedLambda(fun: Gvalue,
                      label: string): Glambda =
  var current = fun
  var seen = initHashSet[NodeKey]()
  while true:
    if not seen.markSeenNode(current):
      raiseValueError(label & " cyclic lambda resolution:\n" & current.nodeRepr)
    if current of Glambda:
      let fn = Glambda(current)
      if fn.isResolvedLambda:
        return fn

    let next = current.nextEvalLambdaExpr(label, "vjpOf")
    if next == nil:
      raiseValueError(label & " unresolved lambda:\n" & current.nodeRepr)
    current = next

# --- structural VJP builder ---

proc buildStructuralLambdaVjp(fun: Gvalue,
                               spec: VjpSpec,
                               ctx: ApplyVjpBuildCtx): Gvalue

proc reduceVjpOfNode(v: Gvalue,
                     label: string,
                     ctx: ApplyVjpBuildCtx): Gvalue =
  var sourceFun = v.inputs[0]
  let sourceTarget = v.vjpOfTarget
  let sourceKind = GvjpOf(v.gfunc).kind

  proc sourceSpec(argProto: Gvalue): VjpSpec =
    case sourceKind
    of lvkCall:
      callVjpSpec(sourceTarget, argProto, v.lambdaResultProto)
    of lvkResult:
      resultVjpSpec(sourceTarget, v.lambdaResultProto, argProto)

  var seen = initHashSet[NodeKey]()
  while true:
    if not seen.markSeenNode(sourceFun):
      raiseValueError(label & " cyclic vjpOf lambda resolution:\n" & sourceFun.nodeRepr)

    sourceFun = followLambdaRefs(sourceFun, label)

    if sourceFun.lambdaResultProto != nil and sourceFun.isCondNode:
      return buildStructuralLambdaVjp(
        sourceFun, sourceSpec(sourceFun.lambdaParamProto), ctx)

    # Conditional lambdas are handled above, so the eval-step walker never
    # reaches its cond branch here.
    let next = sourceFun.nextEvalLambdaExpr(label, "nested vjpOf")
    if next != nil:
      sourceFun = next
      continue

    let fn = directLambda(sourceFun, label)
    if fn == nil:
      if sourceFun of GlambdaRef:
        return nil
      raiseValueError(label & " unresolved lambda:\n" & sourceFun.nodeRepr)
    return buildStructuralLambdaVjp(sourceFun, sourceSpec(fn.param), ctx)

proc preservedVjpShells(ctx: ApplyVjpBuildCtx, shell: Glambda): seq[Gvalue] =
  var seen = initHashSet[NodeKey]()
  for active in ctx.active:
    if seen.markSeenNode(active.shell):
      result.add Gvalue(active.shell)
  if seen.markSeenNode(shell):
    result.add Gvalue(shell)

proc sumGradientPieces(proto: Gvalue, pieces: seq[Gvalue]): Gvalue =
  if pieces.len == 0:
    return proto.zeroLike
  result = pieces[0]
  for i in 1..<pieces.len:
    result = proto.addLike(result, pieces[i])

proc generatedLambda(param: Gvalue,
                     body: Gvalue,
                     preserve: openArray[Gvalue] = [],
                     ctx: ApplyVjpBuildCtx = nil): Glambda =
  let grt = sharedGraphRuntime([param, body], "generated lambda")
  result = Glambda(runtime: grt, param: param, body: body).assignStableNodeId
  # VJP normalization injects active scalar capture targets onto interior apply
  # nodes so backward traversal can reach them without reducing those apply nodes.
  normalizeLambdaForVjp(result, preserve, ctx, body)
  result.updated

proc lowerValueCotangent(value: Gvalue,
                         target: Gvalue,
                         seed: Gvalue): Gvalue =
  if not value.hasFirstClassCotangent:
    raiseValueError(
      "cannot lower a lambda-valued result with a first-class cotangent; " &
      "build a structural VJP with vjpOf/apply instead" &
      "\nvalue: " & value.nodeRepr)
  value.gradSeeded(target, seed)

proc activeVjpShellFor(fun: Gvalue,
                       spec: VjpSpec,
                       bodyProto: Gvalue,
                       ctx: ApplyVjpBuildCtx): Glambda =
  if ctx == nil:
    return nil
  let fn = directLambda(fun, "active VJP function")
  if fn == nil:
    return nil
  for active in ctx.active:
    if active.kind == spec.kind and
        active.owner != nil and active.owner.nodeKey == fn.nodeKey and
        sameVjpTarget(active.target, spec.target) and
        lambdaProtoMatch(active.bodyProto, bodyProto) and
        lambdaProtoMatch(active.argProto, spec.argProto):
      return active.shell
  nil

proc buildFunctionApplyValueVjp(fun: Gvalue,
                                arg: Gvalue,
                                targetValue: Gvalue,
                                seed: Gvalue,
                                ctx: ApplyVjpBuildCtx): Gvalue =
  let appliedFun = fun.inputs[ApplyFunInput]
  let appliedArg = fun.inputs[ApplyArgInput]
  var pieces: seq[Gvalue] = @[]
  let gradProto = targetValue.zeroLike
  let directResultProto = callVjpNodeProto(fun, gradProto)
  let directVjp = apply(
    buildStructuralLambdaVjp(
      appliedFun,
      resultVjpSpec(
        LambdaVjpTarget(kind: lvtkValue, value: targetValue),
        directResultProto,
        arg),
      ctx),
    appliedArg)
  pieces.add apply(apply(directVjp, arg), seed)

  if appliedArg.hasFirstClassCotangent:
    let argGradProto = appliedArg.zeroLike
    let argVjpResultProto = callVjpNodeProto(fun, argGradProto)
    let argVjp = apply(
      buildStructuralLambdaVjp(
        appliedFun,
        callVjpSpec(
          LambdaVjpTarget(kind: lvtkArgument),
          appliedArg,
          argVjpResultProto),
        ctx),
      appliedArg)
    let argCotangent = apply(apply(argVjp, arg), seed)
    pieces.add lowerValueCotangent(appliedArg, targetValue, argCotangent)

  sumGradientPieces(gradProto, pieces)

proc buildApplyVjp(fun: Gvalue,
                    arg: Gvalue,
                    target: LambdaVjpTarget,
                    seed: Gvalue,
                    ctx: ApplyVjpBuildCtx): Gvalue =
  if target.kind == lvtkValue and fun.gfunc of Gapply and
      fun.lambdaResultProto != nil:
    return buildFunctionApplyValueVjp(fun, arg, target.value, seed, ctx)

  let gradProto = targetGradProto(fun, target, arg)
  let resultProto = fun.lambdaResultProto
  if resultProto == nil:
    raiseValueError("apply VJP expects lambda result prototype:\n" & fun.nodeRepr)
  let vjpResultProto = structuralVjpResultProto(resultProto, gradProto)
  let spec = callVjpSpec(target, arg, vjpResultProto)
  let activeShell = activeVjpShellFor(fun, spec, vjpResultProto, ctx)
  if activeShell != nil:
    return apply(apply(Gvalue(activeShell), arg), seed)

  # A non-lambda vjp or seed function fails inside apply's own construction check;
  # builder consistency is a functional-test contract, not a per-use guard here.
  apply(apply(buildStructuralLambdaVjp(fun, spec, ctx), arg), seed)

proc reusableDirectVjpShell(fun: Gvalue,
                            spec: VjpSpec,
                            ctx: ApplyVjpBuildCtx,
                            bodyProto: Gvalue): Glambda =
  if spec.kind == lvkCall:
    # The key (funId + target + argProto, with the resolved lambda pinned by funId)
    # fully determines the shell shape, so a memo hit is always usable as-is.
    let memoValue = ctx.memo.getOrDefault(makeVjpKey(fun, spec))
    if memoValue != nil:
      return Glambda(memoValue)

  let activeShell = activeVjpShellFor(fun, spec, bodyProto, ctx)
  if activeShell != nil:
    if spec.kind == lvkCall:
      ctx.memo[makeVjpKey(fun, spec)] = Gvalue(activeShell)
    return activeShell

proc buildResultVjpShellBody(bodyAtArg: Gvalue,
                             target: LambdaVjpTarget,
                             ctx: ApplyVjpBuildCtx): Gvalue =
  if bodyAtArg.lambdaResultProto == nil:
    raiseValueError(
      "lambda result VJP expects lambda-valued body, got:\n" &
      bodyAtArg.nodeRepr)
  buildStructuralLambdaVjp(
    bodyAtArg,
    callVjpSpec(target, bodyAtArg.lambdaParamProto),
    ctx)

proc buildCallVjpShellBody(bodyAtArg: Gvalue,
                           bodyTarget: Gvalue,
                           ctx: ApplyVjpBuildCtx,
                           shell: Glambda): Gvalue =
  if bodyAtArg.lambdaResultProto != nil:
    return buildStructuralLambdaVjp(
      bodyAtArg,
      callVjpSpec(
        LambdaVjpTarget(kind: lvtkValue, value: bodyTarget),
        bodyAtArg.lambdaParamProto),
      ctx)

  let seedParam = cotangentProto(bodyAtArg)
  seedParam.updated
  let bodyBar = lowerValueCotangent(bodyAtArg, bodyTarget, seedParam)
  generatedLambda(
    seedParam,
    bodyBar,
    preservedVjpShells(ctx, shell),
    ctx)

proc buildDirectLambdaVjp(fun: Gvalue,
                          fn: Glambda,
                          spec: VjpSpec,
                          ctx: ApplyVjpBuildCtx): Gvalue =
  let bodyProto = vjpShellBodyProto(Gvalue(fn), spec)
  let reusable = reusableDirectVjpShell(
    fun,
    spec,
    ctx,
    bodyProto)
  if reusable != nil:
    return reusable

  let argParam = fn.param.newOneOf
  argParam.updated
  let shell = Glambda(
    runtime: fn.runtime,
    param: argParam,
    body: bodyProto).assignStableNodeId
  shell.updated
  if spec.kind == lvkCall:
    ctx.memo[makeVjpKey(fun, spec)] = Gvalue(shell)

  ctx.active.add ActiveVjp(
    kind: spec.kind,
    owner: fn,
    target: spec.target,
    bodyProto: bodyProto,
    argProto: spec.argProto,
    shell: shell)
  try:
    # VJP instantiation threads active scalar capture targets onto interior apply
    # nodes so the shell's backward can reach them.
    let bodyAtArg = instantiateBodyForVjp(fn, argParam, ctx)
    case spec.kind
    of lvkResult:
      shell.body = buildResultVjpShellBody(
        bodyAtArg,
        spec.target,
        ctx)
    of lvkCall:
      let bodyTarget =
        case spec.target.kind
        of lvtkArgument:
          argParam
        of lvtkValue:
          spec.target.value
      shell.body = buildCallVjpShellBody(
        bodyAtArg,
        bodyTarget,
        ctx,
        shell)
  finally:
    ctx.active.setLen(ctx.active.len - 1)

  normalizeLambdaForVjp(
    shell,
    preservedVjpShells(ctx, shell),
    ctx,
    shell.body)
  shell.updated
  shell

proc buildStructuralLambdaVjp(fun: Gvalue,
                               spec: VjpSpec,
                               ctx: ApplyVjpBuildCtx): Gvalue =
  # Dispatch on the structural form of the function and rewrite the VJP through
  # it. Every caller passes a non-nil ctx. Each branch below is one rewrite rule.
  let label =
    case spec.kind
    of lvkCall:
      "apply VJP"
    of lvkResult:
      "lambda result VJP"

  # vjpOf(g): reduce g's source lambda, then VJP that; stay symbolic if unresolved.
  if fun.gfunc of GvjpOf:
    let sourceVjp = reduceVjpOfNode(fun, label, ctx)
    if sourceVjp == nil:
      return newVjpOfNode(fun, spec)
    return buildStructuralLambdaVjp(sourceVjp, spec, ctx)

  # VJP distributes over a conditional lambda:
  #   vjpOf(cond(k, f1, f2)) = cond(k, vjpOf(f1), vjpOf(f2))
  if fun.isCondNode:
    let parts = fun.condParts
    return newCondNode(
      parts.selector,
      buildStructuralLambdaVjp(parts.whenTrue, spec, ctx),
      buildStructuralLambdaVjp(parts.whenFalse, spec, ctx))

  # VJP of a lambda-valued apply propagates through the apply's function position.
  if fun.gfunc of Gapply and fun.lambdaResultProto != nil:
    return apply(
      buildStructuralLambdaVjp(
        fun.inputs[ApplyFunInput],
        resultVjpSpec(
          spec.target,
          vjpApplyResultProto(fun, spec),
          spec.argProto),
        ctx),
      fun.inputs[ApplyArgInput])

  let fn = directLambda(fun, label)
  if fn == nil:
    if fun of GlambdaRef:
      return newVjpOfNode(fun, spec)
    raiseValueError(label & " cannot resolve lambda:\n" & fun.nodeRepr)

  buildDirectLambdaVjp(fun, fn, spec, ctx)

proc evaluatedVjpOfLambda(v: Gvalue,
                          label: string): Gvalue =
  result = reduceVjpOfNode(v, label, newApplyVjpBuildCtx())
  if result == nil:
    raiseValueError(label & " unresolved lambda:\n" & v.inputs[0].nodeRepr)

# --- apply node (forward / inputView / backward / construction) ---

proc applyForward(v: Gvalue) =
  let fn = evaluatedLambda(v.inputs[ApplyFunInput], "apply")
  # ensureInstantiation owns the apply-cache lookup/miss; forward only evaluates
  # and re-copies the (cached) instantiated body into this node.
  let instantiated = v.ensureInstantiation(v.inputs[ApplyArgInput], fn).instantiated
  if v.lambdaResultProto != nil:
    if instantiated.lambdaResultProto == nil:
      raiseValueError(
        "apply expected lambda-valued result, got:\n" &
        instantiated.nodeRepr)
    v.updated
    return
  discard instantiated.eval
  v.valCopy instantiated

proc vjpOf*(fun: Gvalue): Gvalue =
  let target = LambdaVjpTarget(kind: lvtkArgument)
  let argProto = fun.lambdaParamProto
  let proto = vjpNodeProto(fun, callVjpSpec(target, argProto))
  buildStructuralLambdaVjp(
    fun,
    callVjpSpec(target, argProto, proto.lambdaResultProto),
    newApplyVjpBuildCtx())

proc applyInputView(v: Gvalue,
                    mode: InputWalkMode,
                    visit: GnodeVisit) =
  case mode
  of iwmEval:
    let fun = v.inputs[ApplyFunInput]
    let arg = v.inputs[ApplyArgInput]
    visit fun
    visit arg
    for dep in collectLambdaValueDeps([fun], [fun, arg], iwmEval):
      visit dep
    if arg of GlambdaRef and GlambdaRef(arg).binding != nil:
      for dep in collectLambdaValueDeps([arg], [fun, arg], iwmEval):
        visit dep
  of iwmReachable:
    visit v.inputs[ApplyFunInput]
    visit v.inputs[ApplyArgInput]
    for dep in v.applyBackwardDeps:
      visit dep.input
  of iwmBackward:
    for dep in v.applyBackwardDeps:
      visit dep.input

proc applyBackward(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let fun = z.inputs[ApplyFunInput]
  let arg = z.inputs[ApplyArgInput]
  let plan = z.applyBackwardDeps
  let seed = rootedUpstream(zb, z)

  if not seed.hasFirstClassCotangent:
    raiseValueError(
      "apply VJP cannot scale lambda upstream gradient" &
      "\napply: " & z.nodeRepr &
      "\nupstream: " & seed.nodeRepr)

  if seed.isStaticZeroLeaf:
    return input.zeroLike

  if i < 0 or i >= plan.len or plan[i].input.nodeKey != input.nodeKey:
    raiseValueError("apply backward input does not match dependency plan")
  let targets = plan[i].targets

  var pieces: seq[Gvalue] = @[]
  let ctx = newApplyVjpBuildCtx()
  for target in targets:
    pieces.add buildApplyVjp(
      fun,
      arg,
      target,
      seed,
      ctx)

  sumGradientPieces(input, pieces)

let gapply = block:
  var f: Gapply
  new f
  f.forward = applyForward
  f.inputView = applyInputView
  f.backward = applyBackward
  f.name = "apply"
  f

proc apply*(fun: Gvalue, x: Gvalue): Gvalue =
  ## An apply node's raw inputs are just `[fun, arg]`. Active structural value
  ## targets reach an apply only via apply-owned VJP clone paths, so construction
  ## needs no build context.
  let proto = applyResultProto(fun)
  if proto == nil:
    raiseValueError("apply expects a lambda value or lambda placeholder, got: " & fun.nodeRepr)
  discard sharedGraphRuntime([fun, x], "apply")
  graphNode(proto.newOneOf, @[fun, x], Gfunc(gapply), "apply")

proc applyLiteralParamProto(fun: Gvalue): Gvalue =
  let direct = directLambda(fun, "apply literal function")
  if direct != nil:
    return direct.param
  result = fun.lambdaParamProto
  if result == nil:
    raiseValueError(
      "apply literal requires lambda parameter prototype" &
      "\nfunction: " & fun.nodeRepr)

proc apply*(fun: Gvalue,
            x: float): Gvalue =
  let proto = applyLiteralParamProto(fun)
  apply(fun, toGvalue(proto.runtime, x))

proc apply*(fun: Gvalue,
            x: int): Gvalue =
  apply(fun, numericLeafLike(applyLiteralParamProto(fun), x))
