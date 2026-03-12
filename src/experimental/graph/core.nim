#[

- Graph traversals are not thread-safe.
- Backward functions for scalar outputs may receive a nil upstream gradient.

]#

from strutils import join, toHex, strip
import std/tables
import std/sets

type
  NodeKey* = pointer
  NodeTable*[T] = Table[NodeKey, T]
  NodeSet* = HashSet[NodeKey]
  GnodeVisit = proc(n: Gvalue) {.closure.}
  GbranchVisit = proc(tbranch, fbranch: Gvalue) {.closure.}
  InputWalkMode* = enum
    iwmEval, iwmGradSignature, iwmDepend
  GradSigTokenKind* = enum
    gstNode, gstInput, gstCallable
  GradSigToken* = object
    kind*: GradSigTokenKind
    nodePtr*: pointer
  Gfunc* {.acyclic.} = ref object
    ## Represents a graph operation: inputs -> output.
    forward: proc(z: Gvalue)
    backward: proc(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue  ## builds a symbolic backprop contribution
    prepare: proc(z: Gvalue)  ## optional node normalization hook before traversal/eval/grad work
    signature: proc(z: Gvalue, tokens: var seq[GradSigToken])  ## optional cache-signature hook after prepare
    walkInputs: proc(z: Gvalue, mode: InputWalkMode, visit: GnodeVisit, onUnknown: GbranchVisit)
    backwardTarget: proc(zb: Gvalue, z: Gvalue, target: Gvalue, dep: Gvalue): Gvalue
    runCount: int
    name: string
  Gvalue* {.acyclic.} = ref object of RootObj
    ## A value tracks its dependencies, which enables symbolic backpropagation.
    inputs*: seq[Gvalue]
    gfunc*: Gfunc
    locals*: seq[Gvalue]  ## scratch storage shared by forward/backward helper code
    epoch: int

type
  GraphError* = object of Defect
  GraphValueError* = object of GraphError
  GraphUnresolvedValueError* = object of GraphValueError

template raiseError*(msg: string) =
  raise newException(GraphError, msg)

template raiseValueError*(msg: string) =
  raise newException(GraphValueError, msg)

template raiseUnresolvedValueError*(msg: string) =
  raise newException(GraphUnresolvedValueError, msg)

template raiseErrorBaseMethod*(msg: string) =
  raiseError(
    "Base method invoked: " & msg &
    "\nMake sure to pass `--multimethods:on` and check there is a custom method for each derived type.")

var graphDebug* = false

proc newGfunc*(
    forward: proc(z: Gvalue) = nil,
    backward: proc(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue = nil,
    prepare: proc(z: Gvalue) = nil,
    signature: proc(z: Gvalue, tokens: var seq[GradSigToken]) = nil,
    walkInputs: proc(z: Gvalue, mode: InputWalkMode, visit: GnodeVisit, onUnknown: GbranchVisit) = nil,
    backwardTarget: proc(zb: Gvalue, z: Gvalue, target: Gvalue, dep: Gvalue): Gvalue = nil,
    name: string): Gfunc =
  Gfunc(
    forward: forward,
    backward: backward,
    prepare: prepare,
    signature: signature,
    walkInputs: walkInputs,
    backwardTarget: backwardTarget,
    name: name)

proc runCount*(f: Gfunc): int = f.runCount

proc `$`*(x: Gfunc): string

method `$`*(x: Gvalue): string {.base.} =
  let f = x.gfunc
  result = "Gvalue(" & $x.epoch & ")"
  if f != nil:
    result &= " " & $f

proc `$`*(x: Gfunc): string = x.name & "<" & $x.runCount & ">"

proc nodeRepr*(x: Gvalue): string =
  let f = x.gfunc
  result = $x & " (" & $x.epoch & ")" & "@0X" & strip(toHex(cast[int](x)), trailing = false, chars = {'0'})
  if f != nil:
    result &= " " & $f & "@0X" & strip(toHex(cast[int](f)), trailing = false, chars = {'0'})

proc epochOf*(x: Gvalue): int =
  if x == nil:
    return 0
  x.epoch

proc nodeKey*(x: Gvalue): NodeKey {.inline.} =
  cast[NodeKey](x)

proc nodeFromKey(key: NodeKey): Gvalue {.inline.} =
  cast[Gvalue](key)

proc initNodeTable*[T](): NodeTable[T] {.inline.} =
  initTable[NodeKey, T]()

proc initNodeSet*(): NodeSet {.inline.} =
  initHashSet[NodeKey]()

proc hasNode*[T](t: NodeTable[T], x: Gvalue): bool {.inline.} =
  t.hasKey(x.nodeKey)

proc getNode*[T](t: NodeTable[T], x: Gvalue): T {.inline.} =
  t[x.nodeKey]

proc putNode*[T](t: var NodeTable[T], x: Gvalue, value: T) {.inline.} =
  t[x.nodeKey] = value

proc delNode*[T](t: var NodeTable[T], x: Gvalue) {.inline.} =
  t.del(x.nodeKey)

proc nodeOrDefault*[T](t: NodeTable[T], x: Gvalue, default: T): T {.inline.} =
  t.getOrDefault(x.nodeKey, default)

proc containsNode*(s: NodeSet, x: Gvalue): bool {.inline.} =
  x.nodeKey in s

proc inclNode*(s: var NodeSet, x: Gvalue) {.inline.} =
  s.incl x.nodeKey

proc exclNode*(s: var NodeSet, x: Gvalue) {.inline.} =
  s.excl x.nodeKey

method newOneOf*(x: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("newOneOf(" & $x & ")")  ## Derived types should return a zero-initialized sibling node.
method valCopy*(z: Gvalue, x: Gvalue) {.base.} = raiseErrorBaseMethod("valCopy(" & $z & "," & $x & ")")

method isZero*(x: Gvalue): bool {.base.} = raiseErrorBaseMethod("isZero(" & $x & ")")
method update*(x: Gvalue, y: int) {.base.} = raiseErrorBaseMethod("update(" & $x & "," & $y & ")")
method update*(x: Gvalue, y: float) {.base.} = raiseErrorBaseMethod("update(" & $x & "," & $y & ")")

proc constLike*(x: Gvalue, value: int): Gvalue =
  result = x.newOneOf
  result.update value

proc treeRepr*(v: Gvalue): string =
  var seen = initNodeSet()
  var shared = newseq[Gvalue]()
  proc mark(x: Gvalue) =
    if x == nil:
      return
    if seen.containsNode(x):
      if shared.find(x) < 0:
        shared.add x
      return
    seen.inclNode x
    for i in x.inputs:
      mark(i)
  var rendered = initNodeSet()
  proc render(x: Gvalue): seq[string] =
    let si = shared.find x
    result = @[x.nodeRepr]
    if rendered.containsNode(x):
      if si >= 0:
        result[0] &= " #" & $si
      return
    if si >= 0:
      result[0] &= " #" & $si & "#"
    rendered.inclNode x
    for i in x.inputs:
      for ir in render(i):
        result.add("  " & ir)
  mark(v)
  result = render(v).join "\n"

method `-`*(x: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("`-`(" & $x & ")")
method `+`*(x: Gvalue, y: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("`+`(" & $x & ", " & $y & ")")
method `*`*(x: Gvalue, y: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("`*`(" & $x & ", " & $y & ")")
method `-`*(x: Gvalue, y: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("`-`(" & $x & ", " & $y & ")")
method `/`*(x: Gvalue, y: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("`/`(" & $x & ", " & $y & ")")
method exp*(x: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("exp(" & $x & ")")

proc cond*(c: Gvalue, x: Gvalue, y: Gvalue): Gvalue
proc buildGradExpr(dep: Gvalue, x: Gvalue): Gvalue

var gcond: Gfunc

type
  CondInputs = tuple[c, t, f: Gvalue]

method walkSymbolicDeps*(v: Gvalue, visit: proc(n: Gvalue) {.closure.}) {.base.} = discard
method appendSignatureTokens*(v: Gvalue, tokens: var seq[GradSigToken]) {.base.} = discard

proc condInputs(v: Gvalue): CondInputs =
  if v.inputs.len != 3:
    raiseError("cond node requires 3 inputs, got " & $v.inputs.len & ":\n" & v.nodeRepr)
  let c = v.inputs[0]
  if c == nil:
    raiseError("cond node has nil condition input:\n" & v.nodeRepr)
  let t = v.inputs[1]
  if t == nil:
    raiseError("cond node has nil true-branch input:\n" & v.nodeRepr)
  let f = v.inputs[2]
  if f == nil:
    raiseError("cond node has nil false-branch input:\n" & v.nodeRepr)
  (c, t, f)

proc prepareNode*(v: Gvalue) =
  if v == nil:
    raiseError("node preparation received nil node")
  let f = v.gfunc
  if f != nil and f.prepare != nil:
    f.prepare v

proc walkRawInputs(v: Gvalue, visit: GnodeVisit) =
  for j in 0..<v.inputs.len:
    let input = v.inputs[j]
    if input == nil:
      raiseError("node has nil input at index " & $j & ":\n" & v.nodeRepr)
    visit input

proc includesSymbolicDeps(mode: InputWalkMode): bool {.inline.} =
  mode != iwmEval

proc walkDefaultPreparedInputs(v: Gvalue,
                               mode: InputWalkMode,
                               visit: GnodeVisit) =
  v.walkRawInputs(visit)
  if mode.includesSymbolicDeps:
    v.walkSymbolicDeps(visit)

proc walkPreparedInputs(v: Gvalue,
                        mode: InputWalkMode,
                        visit: GnodeVisit,
                        onUnknown: GbranchVisit = nil) =
  if v == nil:
    case mode
    of iwmEval:
      raiseError("eval input walk received nil node")
    of iwmGradSignature:
      raiseError("signature input walk received nil node")
    of iwmDepend:
      raiseError("dependency input walk received nil node")

  let f = v.gfunc
  if f != nil and f.walkInputs != nil:
    f.walkInputs(v, mode, visit, onUnknown)
    return
  v.walkDefaultPreparedInputs(mode, visit)

proc walkPreparedEvalInputs*(v: Gvalue, visit: GnodeVisit) =
  v.walkPreparedInputs(iwmEval, visit)

proc walkPreparedGradSignatureInputs*(v: Gvalue, visit: GnodeVisit) =
  v.walkPreparedInputs(iwmGradSignature, visit)

proc walkPreparedDependInputs*(v: Gvalue,
                               visit: GnodeVisit,
                               onUnknown: GbranchVisit = nil) =
  v.walkPreparedInputs(iwmDepend, visit, onUnknown)

proc walkEvalInputs(v: Gvalue, visit: GnodeVisit) =
  v.prepareNode
  v.walkPreparedInputs(iwmEval, visit)

proc walkGradSignatureInputs(v: Gvalue, visit: GnodeVisit) =
  v.prepareNode
  v.walkPreparedInputs(iwmGradSignature, visit)

proc walkDependInputs(v: Gvalue,
                      visit: GnodeVisit,
                      onUnknown: GbranchVisit = nil) =
  v.prepareNode
  v.walkPreparedInputs(iwmDepend, visit, onUnknown)

proc sameNode(a: Gvalue, b: Gvalue): bool =
  a.nodeKey == b.nodeKey

proc dependsOnPreparedGraph(root: Gvalue, target: Gvalue): bool =
  var seen = initNodeSet()
  proc d(v: Gvalue): bool =
    if v == nil:
      return false
    if sameNode(v, target):
      return true
    if seen.containsNode(v):
      return false
    seen.inclNode v
    var found = false
    proc search(n: Gvalue) =
      if not found and d(n):
        found = true
    v.walkDependInputs(
      search,
      onUnknown = proc(tbranch: Gvalue, fbranch: Gvalue) =
        if not found:
          search(tbranch)
        if not found:
          search(fbranch))
    found
  d(root)

proc condb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  proc branchGrad(c: Gvalue, branch: Gvalue, takeTrue: bool): Gvalue =
    if zb == nil:
      let one = branch.constLike(1)
      let zero = branch.constLike(0)
      if takeTrue:
        return cond(c, one, zero)
      return cond(c, zero, one)

    let zero = branch.constLike(0)
    if takeTrue:
      return cond(c, zb, zero)
    return cond(c, zero, zb)

  case i
  of 0:
    return z.inputs[0].constLike(0)
  of 1:
    # Keep the selector in the returned graph so cached gradients remain live when it changes.
    return branchGrad(z.inputs[0], z.inputs[1], takeTrue = true)
  of 2:
    return branchGrad(z.inputs[0], z.inputs[2], takeTrue = false)
  else:
    raiseValueError("i must be 0, 1, or 2, got: " & $i)

proc condWalkInputs(v: Gvalue,
                    mode: InputWalkMode,
                    visit: GnodeVisit,
                    onUnknown: GbranchVisit) =
  let ci = v.condInputs
  case mode
  of iwmEval:
    visit(ci.c)
    if ci.c.isZero:
      visit(ci.f)
    else:
      visit(ci.t)
  of iwmGradSignature:
    visit(ci.c)
  of iwmDepend:
    visit(ci.c)
    if onUnknown != nil:
      onUnknown(ci.t, ci.f)

proc condBackwardTarget(zb: Gvalue, z: Gvalue, target: Gvalue, dep: Gvalue): Gvalue =
  discard dep
  let ci = z.condInputs
  result = cond(ci.c, buildGradExpr(ci.t, target), buildGradExpr(ci.f, target))
  if zb != nil:
    result = zb * result

proc condf(v: Gvalue) =
  if v.inputs[0].isZero:
    v.valCopy v.inputs[2]
  else:
    v.valCopy v.inputs[1]

gcond = newGfunc(
  forward = condf,
  backward = condb,
  walkInputs = condWalkInputs,
  backwardTarget = condBackwardTarget,
  name = "cond")

proc cond*(c: Gvalue, x: Gvalue, y: Gvalue): Gvalue =
  ## The result type follows the branch values; mismatched branch types fail later in `valCopy`.
  result = y.newOneOf
  result.inputs = @[c, x, y]
  result.gfunc = gcond

proc updated*(x: Gvalue) =
  var epoch {.global.} = 0
  inc epoch
  x.epoch = epoch

proc evaluated*(x: Gvalue) =
  ## Marks a node current with respect to its inputs.
  ## Useful when a forward hook updates node-local state outside the usual eval path.
  var maxep = 0
  for i in x.inputs:
    if maxep < i.epoch:
      maxep = i.epoch
  x.epoch = maxep

proc eval*(v: Gvalue): Gvalue {.discardable.} =
  var seen = initNodeSet()
  proc r(x: Gvalue) =
    if x == nil:
      raiseError("eval traversal encountered nil node")
    if seen.containsNode(x):
      return
    seen.inclNode x
    var maxep = 0
    proc visit(n: Gvalue) =
      n.r
      if maxep < n.epoch:
        maxep = n.epoch
    x.walkEvalInputs(visit)
    if x.epoch < maxep:
      let f = x.gfunc
      if graphDebug:
        var s = "[graph/core] eval: " & x.nodeRepr
        for c in x.inputs:
          s &= "\n  " & c.nodeRepr
        echo s
      if f.forward != nil:
        f.forward x
        x.epoch = maxep
        f.runCount.inc
      else:
        raiseError("inputs.len: " & $x.inputs.len & ", but no forward function defined for:\n" & x.nodeRepr)
  v.r
  v

type
  GradSignature = object
    tokens: seq[GradSigToken]
  GradCacheEntry = object
    hasSignature: bool
    signature: GradSignature
    grads: NodeTable[Gvalue]
  GradBuildContext = object
    dep: Gvalue
    x: Gvalue
    signature: GradSignature
    cache: GradCacheEntry
    relevant: NodeTable[bool]
    order: seq[Gvalue]
    contribs: NodeTable[seq[Gvalue]]
  GradCacheStats* = object
    signatureHits*: int
    signatureMisses*: int
    directHits*: int
    directMisses*: int
    invalidations*: int

var gradCacheByOutput = initNodeTable[GradCacheEntry]()
var gradCacheStats*: GradCacheStats

proc buildGradSignature(dep: Gvalue): GradSignature
proc findGrad*(input: Gvalue, output: Gvalue): Gvalue

proc resetGradCacheStats*() =
  gradCacheByOutput = initNodeTable[GradCacheEntry]()
  gradCacheStats = GradCacheStats()

proc buildGradSignature(dep: Gvalue): GradSignature =
  var sig: GradSignature
  var seen = initNodeSet()
  proc walk(v: Gvalue) =
    if v == nil:
      raiseError("grad signature traversal encountered nil node")
    if seen.containsNode(v):
      return
    seen.inclNode v
    sig.tokens.add GradSigToken(kind: gstNode, nodePtr: v.nodeKey)
    v.prepareNode
    let f = v.gfunc
    for i in v.inputs:
      sig.tokens.add GradSigToken(kind: gstInput, nodePtr: i.nodeKey)
    if f != nil and f.signature != nil:
      f.signature(v, sig.tokens)
    v.appendSignatureTokens(sig.tokens)
    v.walkPreparedGradSignatureInputs(proc(n: Gvalue) = n.walk)
  dep.walk
  result = sig

proc gradIsolated*(dep: Gvalue, x: Gvalue): Gvalue

proc zeroLikeNode(x: Gvalue): Gvalue =
  ## Numeric zero constructor for differentiable value nodes.
  x.constLike(0)

proc buildGradExpr(dep: Gvalue, x: Gvalue): Gvalue =
  if dep == nil or x == nil:
    raiseError("buildGradExpr has nil input")
  if not dep.dependsOnPreparedGraph(x):
    return zeroLikeNode(x)
  var g = dep.gradIsolated(x)
  if g == nil:
    g = zeroLikeNode(x)
  g

proc dumpGradientList* =
  echo "# Gradient Cache:"
  for outputKey, entry in gradCacheByOutput.pairs:
    let output = nodeFromKey(outputKey)
    echo "## output: ",output.nodeRepr
    for inputKey, grad in entry.grads.pairs:
      let input = nodeFromKey(inputKey)
      echo "### w.r.t.: ",input.nodeRepr
      echo grad.treeRepr

proc findGrad*(input: Gvalue, output: Gvalue): Gvalue =
  if not gradCacheByOutput.hasNode(output):
    return nil
  let entry = gradCacheByOutput.getNode(output)
  if not entry.grads.hasNode(input):
    return nil
  entry.grads.getNode(input)

proc sumGradContributions(parts: seq[Gvalue]): Gvalue =
  if parts.len == 0:
    return nil
  result = parts[0]
  for j in 1..<parts.len:
    result = result + parts[j]

proc initGradBuildContext(dep: Gvalue, x: Gvalue): GradBuildContext =
  result.dep = dep
  result.x = x
  result.relevant = initNodeTable[bool]()
  result.contribs = initNodeTable[seq[Gvalue]]()

proc storeGradCache(ctx: GradBuildContext) =
  gradCacheByOutput.putNode(ctx.dep, ctx.cache)

proc addGradContribution(ctx: var GradBuildContext,
                         input: Gvalue,
                         contrib: Gvalue) =
  if input == nil or contrib == nil:
    return
  if not ctx.contribs.hasNode(input):
    ctx.contribs.putNode(input, @[])
  ctx.contribs[input.nodeKey].add contrib

proc prepareGradCache(ctx: var GradBuildContext): Gvalue =
  ctx.signature = ctx.dep.buildGradSignature
  ctx.cache =
    if gradCacheByOutput.hasNode(ctx.dep):
      gradCacheByOutput.getNode(ctx.dep)
    else:
      GradCacheEntry(grads: initNodeTable[Gvalue]())
  if ctx.cache.grads.len == 0:
    ctx.cache.grads = initNodeTable[Gvalue]()

  let sameSignature = ctx.cache.hasSignature and ctx.cache.signature == ctx.signature
  if sameSignature:
    gradCacheStats.signatureHits.inc
    if ctx.cache.grads.hasNode(ctx.x):
      let direct = ctx.cache.grads.getNode(ctx.x)
      if direct != nil:
        gradCacheStats.directHits.inc
        return direct
    gradCacheStats.directMisses.inc
  else:
    gradCacheStats.signatureMisses.inc
    if ctx.cache.hasSignature:
      gradCacheStats.invalidations.inc
    ctx.cache.hasSignature = true
    ctx.cache.signature = ctx.signature
    ctx.cache.grads = initNodeTable[Gvalue]()
    ctx.storeGradCache

  if sameNode(ctx.dep, ctx.x):
    let one = ctx.x.constLike(1)
    ctx.cache.grads.putNode(ctx.x, one)
    ctx.storeGradCache
    return one

proc collectGradBuildPlan(ctx: var GradBuildContext) =
  let dep = ctx.dep
  let x = ctx.x
  var active = initNodeSet()
  var relevant = ctx.relevant
  var order = ctx.order

  proc mark(v: Gvalue): bool =
    if v == nil:
      return false
    if relevant.hasNode(v):
      return relevant.getNode(v)
    if active.containsNode(v):
      return false
    active.inclNode v
    v.prepareNode
    var need = sameNode(v, x)
    proc visit(n: Gvalue) =
      if mark(n):
        need = true
    v.walkPreparedDependInputs(
      visit,
      onUnknown = proc(tbranch: Gvalue, fbranch: Gvalue) =
        if mark(tbranch):
          need = true
        if mark(fbranch):
          need = true)
    active.exclNode v
    relevant.putNode(v, need)
    if need:
      order.add v
    need

  discard mark(dep)
  ctx.relevant = relevant
  ctx.order = order

proc upstreamGradient(ctx: GradBuildContext, v: Gvalue): Gvalue =
  if sameNode(v, ctx.dep):
    return nil
  if ctx.contribs.hasNode(v):
    return sumGradContributions(ctx.contribs.getNode(v))

proc accumulateGradContributions(ctx: var GradBuildContext): Gvalue =
  for j in countdown(ctx.order.high, 0):
    let v = ctx.order[j]
    let hasUpstream = sameNode(v, ctx.dep) or ctx.contribs.hasNode(v)
    if not hasUpstream:
      continue
    let f = v.gfunc
    let vgr = ctx.upstreamGradient(v)
    if f != nil and f.backwardTarget != nil:
      ctx.addGradContribution(ctx.x, f.backwardTarget(vgr, v, ctx.x, ctx.dep))
      continue
    if f == nil:
      continue
    for i in 0..<v.inputs.len:
      let input = v.inputs[i]
      if input == nil:
        raiseError("node has nil input at index " & $i & ":\n" & v.nodeRepr)
      if not ctx.relevant.nodeOrDefault(input, false):
        continue
      if f.backward == nil:
        raiseError(v.nodeRepr & ":" & $i & ":" & input.nodeRepr & ": backward undefined")
      ctx.addGradContribution(input, f.backward(vgr, v, i, ctx.dep))

  if ctx.contribs.hasNode(ctx.x):
    return sumGradContributions(ctx.contribs.getNode(ctx.x))
  zeroLikeNode(ctx.x)

proc storeGradResult(ctx: GradBuildContext, grad: Gvalue) =
  var cache = ctx.cache
  cache.grads.putNode(ctx.x, grad)
  gradCacheByOutput.putNode(ctx.dep, cache)

proc gradImpl(dep: Gvalue, x: Gvalue): Gvalue =
  var ctx = initGradBuildContext(dep, x)
  result = ctx.prepareGradCache
  if result != nil:
    return

  ctx.collectGradBuildPlan
  if ctx.relevant.nodeOrDefault(dep, false):
    result = ctx.accumulateGradContributions
  else:
    result = zeroLikeNode(x)
  ctx.storeGradResult(result)

proc grad*(dep: Gvalue, x: Gvalue): Gvalue =
  gradImpl(dep, x)

proc gradIsolated*(dep: Gvalue, x: Gvalue): Gvalue =
  let savedGradCacheByOutput = gradCacheByOutput
  let savedGradCacheStats = gradCacheStats
  gradCacheByOutput = initNodeTable[GradCacheEntry]()
  gradCacheStats = GradCacheStats()
  try:
    result = gradImpl(dep, x)
  finally:
    gradCacheByOutput = savedGradCacheByOutput
    gradCacheStats = savedGradCacheStats
