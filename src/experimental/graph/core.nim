#[

- Graph traversals are not thread-safe.
- Backward functions for scalar outputs may receive a nil upstream gradient.

]#

from strutils import join, toHex, strip
import std/tables
import std/sets

type
  GnodeVisit = proc(n: Gvalue) {.closure.}
  GbranchVisit = proc(tbranch, fbranch: Gvalue) {.closure.}
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
    walkEvalInputs: proc(z: Gvalue, visit: GnodeVisit)
    walkGradSignatureInputs: proc(z: Gvalue, visit: GnodeVisit)
    walkDependInputs: proc(z: Gvalue, visit: GnodeVisit, onUnknown: GbranchVisit)
    walkGradMarkInputs: proc(z: Gvalue, visit: GnodeVisit, onUnknown: GbranchVisit)
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
    walkEvalInputs: proc(z: Gvalue, visit: GnodeVisit) = nil,
    walkGradSignatureInputs: proc(z: Gvalue, visit: GnodeVisit) = nil,
    walkDependInputs: proc(z: Gvalue, visit: GnodeVisit, onUnknown: GbranchVisit) = nil,
    walkGradMarkInputs: proc(z: Gvalue, visit: GnodeVisit, onUnknown: GbranchVisit) = nil,
    backwardTarget: proc(zb: Gvalue, z: Gvalue, target: Gvalue, dep: Gvalue): Gvalue = nil,
    name: string): Gfunc =
  Gfunc(
    forward: forward,
    backward: backward,
    prepare: prepare,
    signature: signature,
    walkEvalInputs: walkEvalInputs,
    walkGradSignatureInputs: walkGradSignatureInputs,
    walkDependInputs: walkDependInputs,
    walkGradMarkInputs: walkGradMarkInputs,
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

method newOneOf*(x: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("newOneOf(" & $x & ")")  ## Derived types should return a zero-initialized sibling node.
method valCopy*(z: Gvalue, x: Gvalue) {.base.} = raiseErrorBaseMethod("valCopy(" & $z & "," & $x & ")")

method isZero*(x: Gvalue): bool {.base.} = raiseErrorBaseMethod("isZero(" & $x & ")")
method update*(x: Gvalue, y: int) {.base.} = raiseErrorBaseMethod("update(" & $x & "," & $y & ")")
method update*(x: Gvalue, y: float) {.base.} = raiseErrorBaseMethod("update(" & $x & "," & $y & ")")

proc treeRepr*(v: Gvalue): string =
  var seen = initHashSet[pointer]()
  var shared = newseq[Gvalue]()
  proc mark(x: Gvalue) =
    if x == nil:
      return
    let key = cast[pointer](x)
    if key in seen:
      if shared.find(x) < 0:
        shared.add x
      return
    seen.incl key
    for i in x.inputs:
      mark(i)
  var rendered = initHashSet[pointer]()
  proc render(x: Gvalue): seq[string] =
    let si = shared.find x
    result = @[x.nodeRepr]
    let key = cast[pointer](x)
    if key in rendered:
      if si >= 0:
        result[0] &= " #" & $si
      return
    if si >= 0:
      result[0] &= " #" & $si & "#"
    rendered.incl key
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

proc walkPreparedEvalInputs*(v: Gvalue, visit: GnodeVisit) =
  if v == nil:
    raiseError("eval input walk received nil node")

  let f = v.gfunc
  if f != nil and f.walkEvalInputs != nil:
    f.walkEvalInputs(v, visit)
    return

  v.walkRawInputs(visit)

proc walkPreparedGradSignatureInputs*(v: Gvalue, visit: GnodeVisit) =
  if v == nil:
    raiseError("signature input walk received nil node")

  let f = v.gfunc
  if f != nil and f.walkGradSignatureInputs != nil:
    f.walkGradSignatureInputs(v, visit)
    return

  v.walkRawInputs(visit)
  v.walkSymbolicDeps(visit)

proc walkPreparedDependInputs*(v: Gvalue,
                               visit: GnodeVisit,
                               onUnknown: GbranchVisit = nil) =
  if v == nil:
    raiseError("dependency input walk received nil node")

  let f = v.gfunc
  if f != nil and f.walkDependInputs != nil:
    f.walkDependInputs(v, visit, onUnknown)
    return

  v.walkRawInputs(visit)
  v.walkSymbolicDeps(visit)

proc walkPreparedGradMarkInputs*(v: Gvalue,
                                 visit: GnodeVisit,
                                 onUnknown: GbranchVisit = nil) =
  if v == nil:
    raiseError("gradient-mark input walk received nil node")

  let f = v.gfunc
  if f != nil and f.walkGradMarkInputs != nil:
    f.walkGradMarkInputs(v, visit, onUnknown)
    return

  v.walkRawInputs(visit)
  v.walkSymbolicDeps(visit)

proc walkEvalInputs(v: Gvalue, visit: GnodeVisit) =
  v.prepareNode
  v.walkPreparedEvalInputs(visit)

proc walkGradSignatureInputs(v: Gvalue, visit: GnodeVisit) =
  v.prepareNode
  v.walkPreparedGradSignatureInputs(visit)

proc walkDependInputs(v: Gvalue,
                      visit: GnodeVisit,
                      onUnknown: GbranchVisit = nil) =
  v.prepareNode
  v.walkPreparedDependInputs(visit, onUnknown)

proc walkGradMarkInputs(v: Gvalue,
                        visit: GnodeVisit,
                        onUnknown: GbranchVisit = nil) =
  v.prepareNode
  v.walkPreparedGradMarkInputs(visit, onUnknown)

proc sameNode(a: Gvalue, b: Gvalue): bool =
  cast[pointer](a) == cast[pointer](b)

proc dependsOnPreparedGraph(root: Gvalue, target: Gvalue): bool =
  var seen = initHashSet[pointer]()
  proc d(v: Gvalue): bool =
    if v == nil:
      return false
    if sameNode(v, target):
      return true
    let key = cast[pointer](v)
    if key in seen:
      return false
    seen.incl key
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
  proc zeroLikeInput(input: Gvalue): Gvalue =
    result = input.newOneOf
    result.update 0

  proc oneLikeInput(input: Gvalue): Gvalue =
    result = input.newOneOf
    result.update 1

  proc branchGrad(c: Gvalue, branch: Gvalue, takeTrue: bool): Gvalue =
    if zb == nil:
      let one = oneLikeInput(branch)
      let zero = zeroLikeInput(branch)
      if takeTrue:
        return cond(c, one, zero)
      return cond(c, zero, one)

    let zero = branch.newOneOf
    if takeTrue:
      return cond(c, zb, zero)
    return cond(c, zero, zb)

  case i
  of 0:
    return zeroLikeInput(z.inputs[0])
  of 1:
    # Keep the selector in the returned graph so cached gradients remain live when it changes.
    return branchGrad(z.inputs[0], z.inputs[1], takeTrue = true)
  of 2:
    return branchGrad(z.inputs[0], z.inputs[2], takeTrue = false)
  else:
    raiseValueError("i must be 0, 1, or 2, got: " & $i)

proc condWalkEvalInputs(v: Gvalue, visit: GnodeVisit) =
  let ci = v.condInputs
  visit(ci.c)
  if ci.c.isZero:
    visit(ci.f)
  else:
    visit(ci.t)

proc condWalkGradSignatureInputs(v: Gvalue, visit: GnodeVisit) =
  visit(v.condInputs.c)

proc condWalkDependInputs(v: Gvalue,
                          visit: GnodeVisit,
                          onUnknown: GbranchVisit) =
  let ci = v.condInputs
  visit(ci.c)
  if onUnknown != nil:
    onUnknown(ci.t, ci.f)

proc condWalkGradMarkInputs(v: Gvalue,
                            visit: GnodeVisit,
                            onUnknown: GbranchVisit) =
  condWalkDependInputs(v, visit, onUnknown)

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
  walkEvalInputs = condWalkEvalInputs,
  walkGradSignatureInputs = condWalkGradSignatureInputs,
  walkDependInputs = condWalkDependInputs,
  walkGradMarkInputs = condWalkGradMarkInputs,
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
  var seen = initHashSet[pointer]()
  proc r(x: Gvalue) =
    if x == nil:
      raiseError("eval traversal encountered nil node")
    let key = cast[pointer](x)
    if key in seen:
      return
    seen.incl key
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
    grads: Table[pointer, Gvalue]
  GradBuildPlan = object
    relevant: Table[pointer, bool]
    order: seq[Gvalue]
  GradCacheStats* = object
    signatureHits*: int
    signatureMisses*: int
    directHits*: int
    directMisses*: int
    invalidations*: int

var gradCacheByOutput = initTable[pointer, GradCacheEntry]()
var gradCacheStats*: GradCacheStats

proc buildGradSignature(dep: Gvalue): GradSignature
proc findGrad*(input: Gvalue, output: Gvalue): Gvalue

proc resetGradCacheStats*() =
  gradCacheByOutput = initTable[pointer, GradCacheEntry]()
  gradCacheStats = GradCacheStats()

proc buildGradSignature(dep: Gvalue): GradSignature =
  var sig: GradSignature
  var seen = initHashSet[pointer]()
  proc walk(v: Gvalue) =
    if v == nil:
      raiseError("grad signature traversal encountered nil node")
    let key = cast[pointer](v)
    if key in seen:
      return
    seen.incl key
    sig.tokens.add GradSigToken(kind: gstNode, nodePtr: cast[pointer](v))
    v.prepareNode
    let f = v.gfunc
    for i in v.inputs:
      sig.tokens.add GradSigToken(kind: gstInput, nodePtr: cast[pointer](i))
    if f != nil and f.signature != nil:
      f.signature(v, sig.tokens)
    v.appendSignatureTokens(sig.tokens)
    v.walkPreparedGradSignatureInputs(proc(n: Gvalue) = n.walk)
  dep.walk
  result = sig

proc gradIsolated*(dep: Gvalue, x: Gvalue): Gvalue

proc zeroLikeNode(x: Gvalue): Gvalue =
  ## Numeric zero constructor for differentiable value nodes.
  result = x.newOneOf
  result.update 0

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
    let output = cast[Gvalue](outputKey)
    echo "## output: ",output.nodeRepr
    for inputKey, grad in entry.grads.pairs:
      let input = cast[Gvalue](inputKey)
      echo "### w.r.t.: ",input.nodeRepr
      echo grad.treeRepr

proc findGrad*(input: Gvalue, output: Gvalue): Gvalue =
  let outputKey = cast[pointer](output)
  if not gradCacheByOutput.hasKey(outputKey):
    return nil
  let entry = gradCacheByOutput[outputKey]
  let inputKey = cast[pointer](input)
  if not entry.grads.hasKey(inputKey):
    return nil
  entry.grads[inputKey]

proc sumGradContributions(parts: seq[Gvalue]): Gvalue =
  if parts.len == 0:
    return nil
  result = parts[0]
  for j in 1..<parts.len:
    result = result + parts[j]

proc addGradContribution(contribs: var Table[pointer, seq[Gvalue]],
                         input: Gvalue,
                         contrib: Gvalue) =
  if input == nil or contrib == nil:
    return
  let key = cast[pointer](input)
  if not contribs.hasKey(key):
    contribs[key] = @[]
  contribs[key].add contrib

proc prepareGradCache(dep: Gvalue,
                      x: Gvalue,
                      sig: GradSignature,
                      depKey: var pointer,
                      xKey: var pointer,
                      cache: var GradCacheEntry): Gvalue =
  depKey = cast[pointer](dep)
  cache =
    if gradCacheByOutput.hasKey(depKey):
      gradCacheByOutput[depKey]
    else:
      GradCacheEntry(grads: initTable[pointer, Gvalue]())
  if cache.grads.len == 0:
    cache.grads = initTable[pointer, Gvalue]()

  xKey = cast[pointer](x)
  let sameSignature = cache.hasSignature and cache.signature == sig
  if sameSignature:
    gradCacheStats.signatureHits.inc
    if cache.grads.hasKey(xKey):
      let direct = cache.grads[xKey]
      if direct != nil:
        gradCacheStats.directHits.inc
        return direct
    gradCacheStats.directMisses.inc
  else:
    gradCacheStats.signatureMisses.inc
    if cache.hasSignature:
      gradCacheStats.invalidations.inc
    cache.hasSignature = true
    cache.signature = sig
    cache.grads = initTable[pointer, Gvalue]()
    gradCacheByOutput[depKey] = cache

  if sameNode(dep, x):
    var one = x.newOneOf
    one.update 1
    cache.grads[xKey] = one
    gradCacheByOutput[depKey] = cache
    return one

proc collectGradBuildPlan(dep: Gvalue, x: Gvalue): GradBuildPlan =
  var plan: GradBuildPlan
  plan.relevant = initTable[pointer, bool]()
  var active = initHashSet[pointer]()

  proc mark(v: Gvalue): bool =
    if v == nil:
      return false
    let key = cast[pointer](v)
    if plan.relevant.hasKey(key):
      return plan.relevant[key]
    if key in active:
      return false
    active.incl key
    v.prepareNode
    var need = sameNode(v, x)
    proc visit(n: Gvalue) =
      if mark(n):
        need = true
    v.walkPreparedGradMarkInputs(
      visit,
      onUnknown = proc(tbranch: Gvalue, fbranch: Gvalue) =
        if mark(tbranch):
          need = true
        if mark(fbranch):
          need = true)
    active.excl key
    plan.relevant[key] = need
    if need:
      plan.order.add v
    need

  discard mark(dep)
  result = plan

proc accumulateGradContributions(dep: Gvalue,
                                 x: Gvalue,
                                 plan: GradBuildPlan): Gvalue =
  let xKey = cast[pointer](x)
  var contribs = initTable[pointer, seq[Gvalue]]()
  for j in countdown(plan.order.high, 0):
    let v = plan.order[j]
    let vKey = cast[pointer](v)
    let hasUpstream = sameNode(v, dep) or contribs.hasKey(vKey)
    if not hasUpstream:
      continue
    let f = v.gfunc
    let vgr =
      if sameNode(v, dep):
        nil
      elif contribs.hasKey(vKey):
        sumGradContributions(contribs[vKey])
      else:
        nil
    if f != nil and f.backwardTarget != nil:
      addGradContribution(contribs, x, f.backwardTarget(vgr, v, x, dep))
      continue
    if f == nil:
      continue
    for i in 0..<v.inputs.len:
      let input = v.inputs[i]
      if input == nil:
        raiseError("node has nil input at index " & $i & ":\n" & v.nodeRepr)
      let inputKey = cast[pointer](input)
      if not plan.relevant.getOrDefault(inputKey, false):
        continue
      if f.backward == nil:
        raiseError(v.nodeRepr & ":" & $i & ":" & input.nodeRepr & ": backward undefined")
      addGradContribution(contribs, input, f.backward(vgr, v, i, dep))

  if contribs.hasKey(xKey):
    return sumGradContributions(contribs[xKey])
  zeroLikeNode(x)

proc gradImpl(dep: Gvalue, x: Gvalue): Gvalue =
  let sig = dep.buildGradSignature
  var depKey, xKey: pointer
  var cache: GradCacheEntry
  result = prepareGradCache(dep, x, sig, depKey, xKey, cache)
  if result != nil:
    return

  let plan = collectGradBuildPlan(dep, x)
  if plan.relevant.getOrDefault(depKey, false):
    result = accumulateGradContributions(dep, x, plan)
  else:
    result = zeroLikeNode(x)
  cache.grads[xKey] = result
  gradCacheByOutput[depKey] = cache

proc grad*(dep: Gvalue, x: Gvalue): Gvalue =
  gradImpl(dep, x)

proc gradIsolated*(dep: Gvalue, x: Gvalue): Gvalue =
  let savedGradCacheByOutput = gradCacheByOutput
  let savedGradCacheStats = gradCacheStats
  gradCacheByOutput = initTable[pointer, GradCacheEntry]()
  gradCacheStats = GradCacheStats()
  try:
    result = gradImpl(dep, x)
  finally:
    gradCacheByOutput = savedGradCacheByOutput
    gradCacheStats = savedGradCacheStats
