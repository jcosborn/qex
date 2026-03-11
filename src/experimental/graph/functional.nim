import core, scalar
import std/tables
import std/sets

# Section: Functional Node Types and Runtime State

type
  Glocal* {.final.} = ref object of Gvalue
    retProto: Gvalue
    bound: Gvalue
  Gcallable* {.final.} = ref object of Gvalue
    retProto: Gvalue
    bound: Gvalue
  Glambda* {.final.} = ref object of Gvalue
    param: Gvalue
    body: Gvalue
    envParams: seq[Gvalue]
    envValues: seq[Gvalue]
    resolved: bool
  ApplyInputView = object
    inputs: seq[Gvalue]
  ApplySignature = object
    callableKey: pointer
    inputKeys: seq[pointer]
  ApplyCacheEntry = ref object
    hasSignature: bool
    signature: ApplySignature
    freshEpoch: int
    reduced: Gvalue
    partials: Table[pointer, Gvalue]
  DependencyEpochContext = object
    seenCallable: HashSet[pointer]
    seenValues: HashSet[pointer]
  CallableDepContext = ref object
    seenCallable: HashSet[pointer]
    seenValues: HashSet[pointer]
    deps: seq[Gvalue]
  ApplyCacheStats* = object
    reduceHits*: int
    reduceMisses*: int
    partialHits*: int
    partialMisses*: int

var applyCache = initTable[pointer, ApplyCacheEntry]()
var applyCacheStats*: ApplyCacheStats
var gapplyDeferred: Gfunc
var gapplyPartialDeferred: Gfunc

# Section: Private Forward Declarations

# Lambda cloning and closure conversion.
proc cloneWithSubst(v: Gvalue, subst: Table[pointer, Gvalue],
                    memo: var Table[pointer, Gvalue]): Gvalue
proc cloneResolvedLambdaWithSubst(fn: Glambda,
                                  subst: Table[pointer, Gvalue],
                                  memo: var Table[pointer, Gvalue]): Gvalue
proc recloseLambda(fn: Glambda)
proc instantiateLambdaBody(fn: Glambda, x: Gvalue): Gvalue

# Callable resolution and apply normalization.
proc resolveDirectLambda(fun: Gvalue): Glambda
proc resolveLambda(fun: Gvalue): Glambda
proc buildApplyInputView(v: Gvalue): ApplyInputView
proc appendApplySignature(v: Gvalue, tokens: var seq[GradSigToken])
proc collectCallableValueDeps(roots: openArray[Gvalue], deps: var seq[Gvalue])

# Deferred apply-partial node construction.
proc deferredTargetNode(baseInputs: seq[Gvalue],
                        targets: seq[Gvalue],
                        gfunc: Gfunc,
                        label: string): Gvalue
proc applyPartialDeferredNode(baseInputs: seq[Gvalue], targets: seq[Gvalue]): Gvalue
proc applyPartialDeferred(dep: Gvalue, target: Gvalue): Gvalue
proc applyPartialDeferredAppendTarget(z: Gvalue, target: Gvalue): Gvalue
proc applyPartialDeferredContrib(zb: Gvalue, z: Gvalue, target: Gvalue): Gvalue
proc dependencyEpoch(ctx: var DependencyEpochContext, v: Gvalue): int
proc dependencyEpoch(v: Gvalue): int
proc inputDependencyEpoch(ctx: var DependencyEpochContext,
                          inputs: openArray[Gvalue]): int
proc inputDependencyEpoch(inputs: openArray[Gvalue]): int
proc ensureApplyReduction(v: Gvalue): ApplyCacheEntry
proc requireApplyReduction(v: Gvalue): ApplyCacheEntry
proc setApplyReduction(entry: ApplyCacheEntry,
                       callableKey: pointer,
                       inputs: seq[Gvalue],
                       reduced: Gvalue,
                       freshEpoch: int)
proc freshCallableBound(v: Gvalue): Gvalue

# Section: Runtime Limits and Stats Reset

var lambdaResolveDepthLimit* = 64
var applyGradPrepareDepthLimit* = 4096
var applyGradPrepareDepth = 0

proc resetApplyCacheStats*() =
  applyCache = initTable[pointer, ApplyCacheEntry]()
  applyCacheStats = ApplyCacheStats()

# Section: Small Shared Helpers

proc copySeq[T](src: seq[T]): seq[T] =
  result = newseq[T](src.len)
  for i in 0..<src.len:
    result[i] = src[i]

proc updateMax(value: var int, candidate: int) =
  if value < candidate:
    value = candidate

proc initDependencyEpochContext(): DependencyEpochContext =
  result.seenCallable = initHashSet[pointer]()
  result.seenValues = initHashSet[pointer]()

proc markSeenNode(seen: var HashSet[pointer], v: Gvalue): bool =
  let key = cast[pointer](v)
  if key in seen:
    return false
  seen.incl key
  true

proc scanDependencyInputs(ctx: var DependencyEpochContext,
                          value: var int,
                          inputs: openArray[Gvalue]) =
  for input in inputs:
    value.updateMax ctx.dependencyEpoch(input)

# Section: Callable Resolution and Apply Signatures

proc callableKey(fn: Glambda): pointer =
  if fn == nil:
    return nil
  cast[pointer](fn)

proc applySignatureMatches(signature: ApplySignature,
                           callableKey: pointer,
                           inputs: seq[Gvalue]): bool =
  if signature.callableKey != callableKey:
    return false
  if signature.inputKeys.len != inputs.len:
    return false
  for i in 0..<inputs.len:
    if signature.inputKeys[i] != cast[pointer](inputs[i]):
      return false
  true

proc setApplySignature(entry: ApplyCacheEntry,
                       callableKey: pointer,
                       inputs: seq[Gvalue]) =
  entry.signature.callableKey = callableKey
  entry.signature.inputKeys = newseq[pointer](inputs.len)
  for i in 0..<inputs.len:
    entry.signature.inputKeys[i] = cast[pointer](inputs[i])

proc applyCacheEntry(v: Gvalue): ApplyCacheEntry =
  let key = cast[pointer](v)
  if not applyCache.hasKey(key):
    applyCache[key] = ApplyCacheEntry()
  applyCache[key]

proc applyEntryMatches(entry: ApplyCacheEntry,
                       callableKey: pointer,
                       inputs: seq[Gvalue]): bool =
  if not entry.hasSignature or entry.reduced == nil:
    return false
  if not entry.signature.applySignatureMatches(callableKey, inputs):
    return false
  entry.freshEpoch >= inputs.inputDependencyEpoch

proc maxInputEpoch(v: Gvalue): int =
  if v == nil:
    return 0
  for input in v.inputs:
    if input != nil and result < input.epochOf:
      result = input.epochOf

proc freshCallableBound(v: Gvalue): Gvalue =
  ## Reuse an evaluated callable wrapper only while none of its inputs are newer.
  if v == nil or not (v of Gcallable):
    return nil
  let cv = Gcallable(v)
  if cv.bound == nil:
    return nil
  var maxep = v.maxInputEpoch
  var ctx = initDependencyEpochContext()
  let boundEpoch = ctx.dependencyEpoch(cv.bound)
  if maxep < boundEpoch:
    maxep = boundEpoch
  if v.epochOf < maxep:
    return nil
  cv.bound

proc callableResultProto(v: Gvalue): Gvalue =
  if v == nil:
    return nil
  let fn = resolveDirectLambda(v)
  if fn != nil:
    return fn.body
  if v of Gcallable:
    return Gcallable(v).retProto
  if v of Glocal:
    return Glocal(v).retProto
  if v of Glambda:
    return Glambda(v).body
  nil

proc nextCallableBinding(v: Gvalue): Gvalue =
  if v == nil:
    return nil
  if v of Glocal:
    let lf = Glocal(v)
    if lf.bound != nil:
      return lf.bound
  let boundCallable = v.freshCallableBound
  if boundCallable != nil:
    return boundCallable
  if v of Gcallable:
    let cv = Gcallable(v)
    if cv.bound != nil and v.inputs.len == 0:
      return cv.bound
  nil

proc resolveCallableValue(v: Gvalue, reduceApply: bool): Gvalue =
  var current = v
  var depth = 0
  while depth < lambdaResolveDepthLimit:
    if current == nil:
      return nil
    let next = current.nextCallableBinding
    if next != nil:
      current = next
      inc depth
      continue
    if reduceApply and current.gfunc == gapplyDeferred and current.inputs.len >= 2:
      try:
        let entry = current.ensureApplyReduction
        current = entry.reduced
      except GraphUnresolvedValueError:
        return nil
      inc depth
      continue
    return current
  raiseValueError(
    "callable resolution exceeded depth limit " & $lambdaResolveDepthLimit &
    ":\n" & v.nodeRepr)

proc resolveDirectLambda(fun: Gvalue): Glambda =
  let current = fun.resolveCallableValue(reduceApply = false)
  if current == nil or not (current of Glambda):
    return nil
  let fn = Glambda(current)
  if fn.resolved:
    return fn
  nil

proc symbolicCallableToken(v: Gvalue): pointer =
  let current = v.resolveCallableValue(reduceApply = false)
  if current == nil:
    return nil
  let fn = resolveDirectLambda(current)
  if fn != nil:
    return callableKey(fn)
  if current != v:
    return cast[pointer](current)
  nil

proc appendApplySignature(v: Gvalue, tokens: var seq[GradSigToken]) =
  for j in 0..<v.inputs.len:
    if j >= 2:
      break
    let token = symbolicCallableToken(v.inputs[j])
    if token != nil:
      tokens.add GradSigToken(kind: gstCallable, nodePtr: token)
  let view = v.buildApplyInputView
  for j in 2..<view.inputs.len:
    tokens.add GradSigToken(kind: gstInput, nodePtr: cast[pointer](view.inputs[j]))

proc isCallableLike(v: Gvalue): bool =
  v.callableResultProto != nil

# Section: Public Local and Prototype Constructors

proc local*(): Gvalue =
  result = Glocal()
  result.updated
proc local*(retPrototype: Gvalue): Gvalue =
  result = Glocal(retProto: retPrototype)
  result.updated
proc localValue*(prototype: Gvalue): Gvalue = prototype.newOneOf
proc localScalar*(): Gscalar =
  result = Gscalar()
  result.updated
proc localInt*(): Gint =
  result = Gint()
  result.updated

# Section: Functional Node Methods

method newOneOf*(x: Glocal): Gvalue = Glocal(retProto: x.retProto)
method valCopy*(z: Glocal, x: Glocal) =
  z.retProto = x.retProto
  z.bound = x.bound
  z.updated
method valCopy*(z: Glocal, x: Gvalue) =
  z.bound = x
  z.updated
method `$`*(x: Glocal): string =
  if x.bound == nil: "local" else: "local(" & $x.bound & ")"

method newOneOf*(x: Gcallable): Gvalue = Gcallable(retProto: x.retProto)
method valCopy*(z: Gcallable, x: Gcallable) =
  z.retProto = x.retProto
  z.bound = x.bound
method valCopy*(z: Gcallable, x: Gvalue) = z.bound = x
method `$`*(x: Gcallable): string =
  if x.bound == nil: "callable" else: "callable(" & $x.bound & ")"

method newOneOf*(x: Glambda): Gvalue =
  Glambda(
    param: x.param,
    body: x.body,
    envParams: x.envParams,
    envValues: x.envValues,
    resolved: false)
method `$`*(x: Glambda): string =
  result = "lambda(" & $x.param & " -> " & $x.body & ")"
  if x.envValues.len > 0:
    result &= "[env:" & $x.envValues.len & "]"

method walkSymbolicDeps*(v: Glocal, visit: proc(n: Gvalue) {.closure.}) =
  if v.bound != nil:
    visit v.bound

method appendSignatureTokens*(v: Glocal, tokens: var seq[GradSigToken]) =
  let token = symbolicCallableToken(v)
  if token != nil:
    tokens.add GradSigToken(kind: gstCallable, nodePtr: token)

method walkSymbolicDeps*(v: Gcallable, visit: proc(n: Gvalue) {.closure.}) =
  let boundCallable = freshCallableBound(v)
  if boundCallable != nil:
    visit boundCallable
  elif v.inputs.len == 0 and v.bound != nil:
    visit v.bound

method appendSignatureTokens*(v: Gcallable, tokens: var seq[GradSigToken]) =
  let token = symbolicCallableToken(v)
  if token != nil:
    tokens.add GradSigToken(kind: gstCallable, nodePtr: token)

method walkSymbolicDeps*(v: Glambda, visit: proc(n: Gvalue) {.closure.}) =
  for ev in v.envValues:
    if ev != nil:
      visit ev

# Section: Lambda Closure Conversion and Cloning

proc collectCaptureValues(v: Gvalue, bound: var HashSet[pointer],
                          seenCaps: var HashSet[pointer],
                          seenNodes: var HashSet[pointer],
                          caps: var seq[Gvalue]) =
  if v == nil:
    return
  let key = cast[pointer](v)
  if key in seenNodes:
    return
  seenNodes.incl key
  if key in bound:
    return

  if v of Glambda and Glambda(v).resolved:
    let fn = Glambda(v)
    for ev in fn.envValues:
      collectCaptureValues(ev, bound, seenCaps, seenNodes, caps)
    return

  let boundCallable = v.freshCallableBound
  if boundCallable != nil:
    collectCaptureValues(boundCallable, bound, seenCaps, seenNodes, caps)
    return

  if v.gfunc == nil and v.inputs.len == 0:
    if key notin seenCaps and not (v of Glambda):
      seenCaps.incl key
      caps.add v
    return

  for i in v.inputs:
    collectCaptureValues(i, bound, seenCaps, seenNodes, caps)

proc trimSubstForLambda(subst: Table[pointer, Gvalue], fn: Glambda): Table[pointer, Gvalue] =
  result = subst
  let pkey = cast[pointer](fn.param)
  if result.hasKey(pkey):
    result.del(pkey)
  for ep in fn.envParams:
    let ekey = cast[pointer](ep)
    if result.hasKey(ekey):
      result.del(ekey)

proc cloneResolvedLambdaWithSubst(fn: Glambda,
                                  subst: Table[pointer, Gvalue],
                                  memo: var Table[pointer, Gvalue]): Gvalue =
  let key = cast[pointer](fn)
  result = Glambda(param: fn.param, resolved: fn.resolved)
  memo[key] = result

  let r = Glambda(result)
  if fn.envParams.len > 0:
    r.envParams = fn.envParams.copySeq
  if fn.envValues.len > 0:
    r.envValues = newseq[Gvalue](fn.envValues.len)
    for i in 0..<fn.envValues.len:
      r.envValues[i] = cloneWithSubst(fn.envValues[i], subst, memo)
  if fn.body != nil:
    let innerSubst = subst.trimSubstForLambda(fn)
    var innerMemo = initTable[pointer, Gvalue]()
    r.body = cloneWithSubst(fn.body, innerSubst, innerMemo)

proc bindSubst(subst: var Table[pointer, Gvalue], key, value: Gvalue) =
  subst[cast[pointer](key)] = value

proc cloneWithFreshMemo(v: Gvalue, subst: Table[pointer, Gvalue]): Gvalue =
  var memo = initTable[pointer, Gvalue]()
  cloneWithSubst(v, subst, memo)

proc bindLambdaEnvValues(subst: var Table[pointer, Gvalue], fn: Glambda) =
  if fn.envParams.len != fn.envValues.len:
    raiseValueError("lambda env arity mismatch during instantiation")
  for i in 0..<fn.envParams.len:
    subst.bindSubst(fn.envParams[i], fn.envValues[i])

proc cloneWithSubst(v: Gvalue, subst: Table[pointer, Gvalue],
                    memo: var Table[pointer, Gvalue]): Gvalue =
  if v == nil:
    return nil

  let key = cast[pointer](v)
  if subst.hasKey(key):
    return subst[key]
  if memo.hasKey(key):
    return memo[key]

  if v of Glambda and Glambda(v).resolved:
    return Glambda(v).cloneResolvedLambdaWithSubst(subst, memo)

  let boundCallable = v.freshCallableBound
  if boundCallable != nil:
    result = cloneWithSubst(boundCallable, subst, memo)
    memo[key] = result
    return result

  if v.gfunc == nil and v.inputs.len == 0:
    return v

  result = v.newOneOf
  memo[key] = result
  if v.inputs.len > 0:
    result.inputs = newseq[Gvalue](v.inputs.len)
    for i in 0..<v.inputs.len:
      result.inputs[i] = cloneWithSubst(v.inputs[i], subst, memo)
  result.gfunc = v.gfunc
  if v.locals.len > 0:
    result.locals = newseq[Gvalue](v.locals.len)

proc openLambdaBody(fn: Glambda): Gvalue =
  if fn.envParams.len == 0:
    return fn.body
  var subst = initTable[pointer, Gvalue]()
  subst.bindLambdaEnvValues(fn)
  fn.body.cloneWithFreshMemo(subst)

proc collectLambdaCaptures(fn: Glambda): seq[Gvalue] =
  var bound = initHashSet[pointer]()
  bound.incl cast[pointer](fn.param)
  var seenCaps = initHashSet[pointer]()
  var seenNodes = initHashSet[pointer]()
  collectCaptureValues(fn.body, bound, seenCaps, seenNodes, result)

proc rewriteLambdaCaptures(fn: Glambda, captures: seq[Gvalue]) =
  fn.envValues = captures
  fn.envParams.setLen(0)
  if captures.len == 0:
    return

  fn.envParams = newseq[Gvalue](captures.len)
  var subst = initTable[pointer, Gvalue]()
  for i in 0..<captures.len:
    let ep = localValue(captures[i])
    fn.envParams[i] = ep
    subst.bindSubst(captures[i], ep)
  fn.body = fn.body.cloneWithFreshMemo(subst)

proc recloseLambda(fn: Glambda) =
  if fn == nil or not fn.resolved:
    return
  if fn.body == nil:
    raiseValueError("resolved lambda body cannot be nil")

  fn.body = openLambdaBody(fn)
  let captures = fn.collectLambdaCaptures
  fn.rewriteLambdaCaptures(captures)

proc instantiateLambdaBody(fn: Glambda, x: Gvalue): Gvalue =
  var subst = initTable[pointer, Gvalue]()
  subst.bindSubst(fn.param, x)
  subst.bindLambdaEnvValues(fn)
  fn.body.cloneWithFreshMemo(subst)

proc resolveLambda(fun: Gvalue): Glambda =
  let current = fun.resolveCallableValue(reduceApply = true)
  if current == nil or not (current of Glambda):
    return nil
  let fn = Glambda(current)
  if fn.resolved:
    return fn
  nil

# Section: Dependency Epoch Tracking

proc dependencyEpoch(ctx: var DependencyEpochContext, v: Gvalue): int =
  if v == nil:
    return 0

  if v of Glambda and Glambda(v).resolved:
    if not ctx.seenCallable.markSeenNode(v):
      return 0
    let fn = Glambda(v)
    ctx.scanDependencyInputs(result, fn.envValues)
    return result

  if v of Gcallable:
    if not ctx.seenCallable.markSeenNode(v):
      return 0
    let cv = Gcallable(v)
    if cv.bound != nil:
      result.updateMax ctx.dependencyEpoch(cv.bound)
    ctx.scanDependencyInputs(result, v.inputs)
    return result

  if v of Glocal and Glocal(v).retProto != nil:
    if not ctx.seenCallable.markSeenNode(v):
      return 0
    result.updateMax v.epochOf
    let lv = Glocal(v)
    if lv.bound != nil:
      result.updateMax ctx.dependencyEpoch(lv.bound)
    ctx.scanDependencyInputs(result, v.inputs)
    return result

  if not ctx.seenValues.markSeenNode(v):
    return 0
  result = v.epochOf
  ctx.scanDependencyInputs(result, v.inputs)

proc dependencyEpoch(v: Gvalue): int =
  if v == nil:
    return 0
  var ctx = initDependencyEpochContext()
  result = ctx.dependencyEpoch(v)

proc inputDependencyEpoch(ctx: var DependencyEpochContext,
                          inputs: openArray[Gvalue]): int =
  for input in inputs:
    if input != nil:
      result.updateMax ctx.dependencyEpoch(input)

proc inputDependencyEpoch(inputs: openArray[Gvalue]): int =
  var ctx = initDependencyEpochContext()
  result = ctx.inputDependencyEpoch(inputs)

# Section: Public Lambda Construction and Apply Result Prototypes

proc lambda*(param: Gvalue, body: Gvalue): Gvalue =
  if param == nil:
    raiseValueError("lambda parameter cannot be nil")
  if body == nil:
    raiseValueError("lambda body cannot be nil")
  let fn = Glambda(param: param, body: body, resolved: true)
  recloseLambda(fn)
  result = fn
  result.updated

proc materializeApplyResultProto(proto: Gvalue): Gvalue =
  if proto == nil:
    return nil
  let nestedProto = proto.callableResultProto
  if nestedProto != nil:
    return Gcallable(retProto: materializeApplyResultProto(nestedProto))
  proto.newOneOf

proc applyResultProto(fun: Gvalue): Gvalue =
  let retProto = fun.callableResultProto
  if retProto == nil:
    return nil
  materializeApplyResultProto(retProto)

# Section: Deferred Apply-Partial Node Construction

proc deferredTargetNode(baseInputs: seq[Gvalue],
                        targets: seq[Gvalue],
                        gfunc: Gfunc,
                        label: string): Gvalue =
  if baseInputs.len == 0:
    raiseValueError(label & " requires at least one base input")
  if targets.len == 0:
    raiseValueError(label & " requires at least one target")
  for j in 0..<baseInputs.len:
    if baseInputs[j] == nil:
      raiseError(label & " base input cannot be nil")
  for target in targets:
    if target == nil:
      raiseError(label & " target cannot be nil")
  result = targets[^1].newOneOf
  result.inputs = newseq[Gvalue](baseInputs.len + targets.len)
  for j in 0..<baseInputs.len:
    result.inputs[j] = baseInputs[j]
  for j in 0..<targets.len:
    result.inputs[baseInputs.len + j] = targets[j]
  result.gfunc = gfunc

proc applyPartialDeferredNode(baseInputs: seq[Gvalue], targets: seq[Gvalue]): Gvalue =
  if baseInputs.len != 1:
    raiseValueError("applyPartialDeferred base input length mismatch: expected 1, got: " &
      $baseInputs.len)
  deferredTargetNode(baseInputs, targets, gapplyPartialDeferred, "applyPartialDeferred")

proc applyPartialDeferred(dep: Gvalue, target: Gvalue): Gvalue =
  if dep == nil or target == nil:
    raiseError("applyPartialDeferred has nil input")
  applyPartialDeferredNode(@[dep], @[target])

proc applyPartialDeferredAppendTarget(z: Gvalue, target: Gvalue): Gvalue =
  if target == nil:
    raiseError("applyPartialDeferred append target cannot be nil")
  if z.gfunc != gapplyPartialDeferred:
    raiseValueError("applyPartialDeferred append target expects applyPartialDeferred node")
  let baseLen = 1
  if z.inputs.len < baseLen + 1:
    raiseValueError("applyPartialDeferred node missing target inputs")
  let existingTargetsLen = z.inputs.len - baseLen
  var targets = newseq[Gvalue](existingTargetsLen + 1)
  for j in 0..<existingTargetsLen:
    targets[j] = z.inputs[baseLen + j]
  targets[existingTargetsLen] = target
  applyPartialDeferredNode(z.inputs[0..<baseLen], targets)

proc applyPartialDeferredContrib(zb: Gvalue, z: Gvalue, target: Gvalue): Gvalue =
  result = applyPartialDeferredAppendTarget(z, target)
  if zb != nil:
    result = zb * result

# Section: Callable Dependency Collection and Apply Views

proc collectCallableValueDeps(roots: openArray[Gvalue], deps: var seq[Gvalue]) =
  let ctx = CallableDepContext(
    seenCallable: initHashSet[pointer](),
    seenValues: initHashSet[pointer](),
    deps: deps)

  for dep in deps:
    if dep != nil:
      ctx.seenValues.incl cast[pointer](dep)

  proc collect(v: Gvalue) =
    if v == nil:
      return
    if isCallableLike(v):
      let ckey = cast[pointer](v)
      if ckey in ctx.seenCallable:
        return
      ctx.seenCallable.incl ckey
      if v.gfunc != nil:
        v.prepareNode
      let fn = resolveDirectLambda(v)
      if fn != nil:
        for ev in fn.envValues:
          collect(ev)
        return
      let boundCallable = v.freshCallableBound
      if boundCallable != nil:
        collect(boundCallable)
        return
      v.walkSymbolicDeps(proc(n: Gvalue) =
        collect(n))
      for input in v.inputs:
        collect(input)
      return
    let vkey = cast[pointer](v)
    if vkey in ctx.seenValues:
      return
    ctx.seenValues.incl vkey
    ctx.deps.add v

  for root in roots:
    collect(root)
  deps = ctx.deps

proc buildApplyInputView(v: Gvalue): ApplyInputView =
  if v.inputs.len < 2:
    raiseValueError("apply node requires at least two inputs")
  let f = v.inputs[0]
  let x = v.inputs[1]
  result.inputs = @[f, x]

  let fn = resolveDirectLambda(f)
  if fn != nil:
    for ev in fn.envValues:
      result.inputs.add ev

  collectCallableValueDeps([f, x], result.inputs)

proc walkApplyInputView(v: Gvalue, visit: proc(n: Gvalue) {.closure.}) =
  let view = v.buildApplyInputView
  for input in view.inputs:
    visit input

proc emptyLike(x: Gvalue): Gvalue =
  ## Type-preserving neutral constructor; callable-like nodes may not support update(0).
  x.newOneOf

proc resolvedCallableValue(v: Gvalue): Gvalue =
  ## Collapse freshly-bound callable wrappers to their callable descriptor after eval.
  let direct = v.resolveCallableValue(reduceApply = false)
  if direct == nil:
    return nil
  let resolved = direct.resolveCallableValue(reduceApply = true)
  if resolved != nil:
    return resolved
  direct

# Section: Apply Cache and Lazy Reduction

proc setApplyReduction(entry: ApplyCacheEntry,
                       callableKey: pointer,
                       inputs: seq[Gvalue],
                       reduced: Gvalue,
                       freshEpoch: int) =
  entry.hasSignature = true
  entry.setApplySignature(callableKey, inputs)
  entry.freshEpoch = freshEpoch
  entry.reduced = reduced
  entry.partials = initTable[pointer, Gvalue]()

proc ensureApplyReduction(v: Gvalue): ApplyCacheEntry =
  let view = v.buildApplyInputView
  let fn = resolveLambda(view.inputs[0])
  if fn == nil:
    raiseUnresolvedValueError("deferred apply unresolved at eval: " & view.inputs[0].nodeRepr)
  let key = callableKey(fn)

  result = v.applyCacheEntry
  if result.applyEntryMatches(key, view.inputs):
    applyCacheStats.reduceHits.inc
    return
  applyCacheStats.reduceMisses.inc

  let freshEpoch = view.inputs.inputDependencyEpoch
  let reduced = instantiateLambdaBody(fn, view.inputs[1])
  if reduced of Glambda and Glambda(reduced).resolved:
    recloseLambda(Glambda(reduced))

  result.setApplyReduction(key, view.inputs, reduced, freshEpoch)

proc requireApplyReduction(v: Gvalue): ApplyCacheEntry =
  result = v.ensureApplyReduction
  if result.reduced == nil:
    raiseValueError("deferred apply preparation did not produce reduced body")

# Section: Lazy Apply-Partials and Runtime Hooks

proc ensureApplyPartial(z: Gvalue, target: Gvalue): Gvalue =
  if target == nil:
    raiseValueError("apply partial target cannot be nil")
  let entry = z.requireApplyReduction
  let targetKey = cast[pointer](target)
  if entry.partials.hasKey(targetKey):
    applyCacheStats.partialHits.inc
    return entry.partials[targetKey]
  applyCacheStats.partialMisses.inc
  if isCallableLike(target):
    result = emptyLike(target)
    entry.partials[targetKey] = result
    return
  inc applyGradPrepareDepth
  defer:
    dec applyGradPrepareDepth
  if applyGradPrepareDepth > applyGradPrepareDepthLimit:
    raiseValueError(
      "apply partial materialization exceeded depth limit " &
      $applyGradPrepareDepthLimit &
      "\napply: " & z.nodeRepr &
      "\ntarget: " & target.nodeRepr)
  let reduced = entry.reduced
  let g = reduced.gradIsolated(target)
  if g == nil:
    result = emptyLike(target)
  else:
    result = g
  entry.partials[targetKey] = result

proc applyDeferredf(v: Gvalue) =
  let entry = v.requireApplyReduction
  let reduced = entry.reduced
  let callable = resolveLambda(reduced)
  if callable != nil:
    v.valCopy callable
    return
  discard reduced.eval
  v.valCopy reduced.resolvedCallableValue

proc applyDeferredContribution(zb: Gvalue, z: Gvalue, target: Gvalue): Gvalue =
  if target == nil or isCallableLike(target):
    return nil
  let partial = applyPartialDeferred(z, target)
  if zb != nil and zb.isCallableLike:
    return nil
  if zb == nil:
    return partial
  zb * partial

proc applyDeferredb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard dep
  if i < 0 or i >= z.inputs.len:
    raiseValueError("apply backward input index out of range: " & $i)
  if i == 0:
    return nil
  applyDeferredContribution(zb, z, z.inputs[i])

proc applyDeferredBackwardTarget(zb: Gvalue,
                                 z: Gvalue,
                                 target: Gvalue,
                                 dep: Gvalue): Gvalue =
  discard dep
  applyDeferredContribution(zb, z, target)

proc applyDeferredWalkEvalInputs(v: Gvalue, visit: proc(n: Gvalue) {.closure.}) =
  v.walkApplyInputView(visit)

proc applyDeferredWalkGradSignatureInputs(v: Gvalue, visit: proc(n: Gvalue) {.closure.}) =
  v.walkApplyInputView(visit)

proc walkApplyInputViewWithUnknown(v: Gvalue,
                                   visit: proc(n: Gvalue) {.closure.},
                                   onUnknown: proc(tbranch, fbranch: Gvalue) {.closure.}) =
  let _ = onUnknown
  v.walkApplyInputView(visit)

proc applyDeferredWalkDependInputs(v: Gvalue,
                                   visit: proc(n: Gvalue) {.closure.},
                                   onUnknown: proc(tbranch, fbranch: Gvalue) {.closure.}) =
  v.walkApplyInputViewWithUnknown(visit, onUnknown)

proc applyDeferredWalkGradMarkInputs(v: Gvalue,
                                     visit: proc(n: Gvalue) {.closure.},
                                     onUnknown: proc(tbranch, fbranch: Gvalue) {.closure.}) =
  v.walkApplyInputViewWithUnknown(visit, onUnknown)

proc applyPartialDeferredf(v: Gvalue) =
  if v.inputs.len < 2:
    raiseValueError("applyPartialDeferred requires a base apply and at least one target")
  let base = v.inputs[0]
  if base == nil:
    raiseError("applyPartialDeferred has nil base apply")
  var expr = base.ensureApplyPartial(v.inputs[1])
  for j in 2..<v.inputs.len:
    let target = v.inputs[j]
    if target == nil:
      raiseError("applyPartialDeferred target cannot be nil")
    var g = expr.gradIsolated(target)
    if g == nil:
      g = emptyLike(target)
    expr = g
  discard expr.eval
  v.valCopy expr

proc applyPartialDeferredb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  if i < 0 or i >= z.inputs.len:
    raiseValueError("applyPartialDeferred backward input index out of range: " & $i)
  emptyLike(z.inputs[i])

proc requireApplyPartialBase(v: Gvalue): Gvalue =
  if v.inputs.len < 2:
    raiseValueError("applyPartialDeferred node requires a base apply and at least one target")
  result = v.inputs[0]
  if result == nil:
    raiseError("applyPartialDeferred node has nil base input:\n" & v.nodeRepr)

proc applyPartialDeferredWalkEvalInputs(v: Gvalue, visit: proc(n: Gvalue) {.closure.}) =
  let base = v.requireApplyPartialBase
  var seen = initHashSet[pointer]()
  proc walkEvalDeps(n: Gvalue) =
    if n == nil:
      raiseError("applyPartialDeferred eval dependency walk encountered nil node")
    let key = cast[pointer](n)
    if key in seen:
      return
    seen.incl key
    n.prepareNode
    if n.gfunc == gapplyDeferred or n.gfunc == gapplyPartialDeferred:
      n.walkPreparedDependInputs(walkEvalDeps)
      return
    n.walkPreparedEvalInputs(walkEvalDeps)
    n.walkSymbolicDeps(walkEvalDeps)
    visit n
  walkEvalDeps(base)

proc applyPartialDeferredWalkGradSignatureInputs(v: Gvalue, visit: proc(n: Gvalue) {.closure.}) =
  visit(v.requireApplyPartialBase)

proc visitApplyPartialBaseWithUnknown(v: Gvalue,
                                      visit: proc(n: Gvalue) {.closure.},
                                      onUnknown: proc(tbranch, fbranch: Gvalue) {.closure.}) =
  let _ = onUnknown
  visit(v.requireApplyPartialBase)

proc applyPartialDeferredWalkDependInputs(v: Gvalue,
                                          visit: proc(n: Gvalue) {.closure.},
                                          onUnknown: proc(tbranch, fbranch: Gvalue) {.closure.}) =
  v.visitApplyPartialBaseWithUnknown(visit, onUnknown)

proc applyPartialDeferredWalkGradMarkInputs(v: Gvalue,
                                            visit: proc(n: Gvalue) {.closure.},
                                            onUnknown: proc(tbranch, fbranch: Gvalue) {.closure.}) =
  v.visitApplyPartialBaseWithUnknown(visit, onUnknown)

proc applyPartialDeferredBackwardTarget(zb: Gvalue,
                                        z: Gvalue,
                                        target: Gvalue,
                                        dep: Gvalue): Gvalue =
  discard dep
  applyPartialDeferredContrib(zb, z, target)

gapplyPartialDeferred = newGfunc(
  forward = applyPartialDeferredf,
  backward = applyPartialDeferredb,
  walkEvalInputs = applyPartialDeferredWalkEvalInputs,
  walkGradSignatureInputs = applyPartialDeferredWalkGradSignatureInputs,
  walkDependInputs = applyPartialDeferredWalkDependInputs,
  walkGradMarkInputs = applyPartialDeferredWalkGradMarkInputs,
  backwardTarget = applyPartialDeferredBackwardTarget,
  name = "applyPartialDeferred")
gapplyDeferred = newGfunc(
  forward = applyDeferredf,
  backward = applyDeferredb,
  walkEvalInputs = applyDeferredWalkEvalInputs,
  walkGradSignatureInputs = applyDeferredWalkGradSignatureInputs,
  walkDependInputs = applyDeferredWalkDependInputs,
  walkGradMarkInputs = applyDeferredWalkGradMarkInputs,
  backwardTarget = applyDeferredBackwardTarget,
  signature = appendApplySignature,
  name = "applyDeferred")

# Section: Public Apply Construction

proc apply*(fun: Gvalue, x: Gvalue): Gvalue =
  if fun == nil:
    raiseValueError("apply function cannot be nil")
  if x == nil:
    raiseValueError("apply argument cannot be nil")

  let proto = applyResultProto(fun)
  if proto == nil:
    raiseValueError("apply expects a lambda or local function placeholder, got: " & fun.nodeRepr)

  result = proto.newOneOf
  result.inputs = @[fun, x]
  result.gfunc = gapplyDeferred

# Section: Demonstrations

when isMainModule:
  block:
    let x = toGvalue(3.0)
    let y = toGvalue(2.0)
    let v = localScalar()
    let z = apply(lambda(v, v + v), x * y)
    echo "## z (before eval)"
    echo z.treeRepr
    z.eval
    echo "## z (after eval)"
    echo z.treeRepr
    echo "z = ",z
    let dzdx = z.grad x
    echo "## dzdx (before eval)"
    echo dzdx.treeRepr
    dzdx.eval
    echo "## dzdx (after eval)"
    echo dzdx.treeRepr
    echo "dzdx = ",dzdx
  block:
    let f = local(localScalar())
    let u = localScalar()
    let hof = lambda(f, lambda(u, apply(f, u) + 1.0))
    let v = localScalar()
    let a = toGvalue(2.0)
    let g = lambda(v, a * v)
    let z = apply(apply(hof, g), 3.0)
    echo "## z (before eval)"
    echo z.treeRepr
    z.eval
    echo "## z (after eval)"
    echo z.treeRepr
    echo "z = ",z
    let dzda = z.grad a
    echo "## dzda (before eval)"
    echo dzda.treeRepr
    dzda.eval
    echo "## dzda (after eval)"
    echo dzda.treeRepr
    echo "dzda = ",dzda
  block:
    # Short-depth recursion via Y combinator: F(n) = if n == 0 then base else F(n-1) + step
    let protoArg = localScalar()
    let protoRet = localScalar()
    let fnProto = lambda(protoArg, protoRet)
    let x = local(fnProto)
    let f = local(fnProto)
    let Y = lambda(f, apply(lambda(x, apply(f, apply(x, x))), lambda(x, apply(f, apply(x, x)))))

    let rf = local(localScalar())
    let u = localScalar()
    let v = localScalar()
    let base = toGvalue(1.0)
    let step = toGvalue(1.0)
    let F = lambda(rf, lambda(u,
      cond(equal(u, 0.0), base,
        apply(lambda(v, apply(rf, v) + step), u - 1.0))))

    let z = apply(apply(Y, F), 3.0)
    echo "## Y recursion z (before eval)"
    echo z.treeRepr
    z.eval
    echo "## Y recursion z (after eval)"
    echo z.treeRepr
    echo "Y recursion z = ",z
    let dzdstep = z.grad step
    echo "## Y recursion dzdstep (before eval)"
    echo dzdstep.treeRepr
    dzdstep.eval
    echo "## Y recursion dzdstep (after eval)"
    echo dzdstep.treeRepr
    echo "Y recursion dzdstep = ",dzdstep
    let dzdbase = z.grad base
    echo "## Y recursion dzdbase (before eval)"
    echo dzdbase.treeRepr
    dzdbase.eval
    echo "## Y recursion dzdbase (after eval)"
    echo dzdbase.treeRepr
    echo "Y recursion dzdbase = ",dzdbase
