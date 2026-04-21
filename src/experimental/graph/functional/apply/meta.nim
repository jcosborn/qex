import std/algorithm
import ../../core
import ../callable

type
  ApplyDeps = object
    extra: seq[Gvalue]
    inputKeys: seq[NodeId]
  ApplyAnalysis* = object
    fun*: Gvalue
    arg*: Gvalue
    extraDeps*: seq[Gvalue]
    inputKeys*: seq[NodeId]
    fn*: Glambda
    callableKey*: NodeId

proc addUniqueApplyInput(inputs: var seq[Gvalue],
                         seen: var NodeSet,
                         value: Gvalue) =
  if value != nil and seen.markSeenNode(value):
    inputs.add value

proc requireApplyInputs(v: Gvalue): tuple[fun: Gvalue, arg: Gvalue] =
  v.requireInputCountExactly(2, "apply node")
  (
    fun: v.requireNodeInput(0, "apply node", "function"),
    arg: v.requireNodeInput(1, "apply node", "argument"))

proc collectApplyValueDeps(fun: Gvalue,
                           arg: Gvalue,
                           deps: var seq[Gvalue]) =
  var seed = @[fun, arg]
  seed.add deps
  deps.add collectCallableValueDeps([fun, arg], seed)

proc collectApplyDeps(fun: Gvalue,
                      arg: Gvalue): ApplyDeps =
  fun.collectApplyValueDeps(arg, result.extra)
  var seen = initNodeSet()
  var cacheInputs: seq[Gvalue] = @[]
  cacheInputs.addUniqueApplyInput(seen, fun)
  cacheInputs.addUniqueApplyInput(seen, arg)
  for dep in result.extra:
    cacheInputs.addUniqueApplyInput(seen, dep)
  result.inputKeys = newseq[NodeId](cacheInputs.len)
  for i in 0..<cacheInputs.len:
    result.inputKeys[i] = cacheInputs[i].stableNodeId
  result.inputKeys.sort()

proc analyzeApply*(v: Gvalue): ApplyAnalysis =
  (result.fun, result.arg) = v.requireApplyInputs
  let deps = result.fun.collectApplyDeps(result.arg)
  result.extraDeps = deps.extra
  result.inputKeys = deps.inputKeys

proc prepareApplyReduction*(analysis: var ApplyAnalysis) =
  analysis.fn = resolveLambda(analysis.fun)
  if analysis.fn != nil:
    analysis.callableKey = callableKey(analysis.fn)

proc walkApplyDeps*(analysis: ApplyAnalysis,
                    visit: proc(n: Gvalue) {.closure.}) =
  visit analysis.fun
  visit analysis.arg
  for dep in analysis.extraDeps:
    visit dep

proc walkApplyDeps*(v: Gvalue,
                    visit: proc(n: Gvalue) {.closure.}) =
  v.analyzeApply.walkApplyDeps(visit)

proc appendApplySignature*(analysis: ApplyAnalysis,
                           tokens: var seq[GradSigToken]) =
  for input in [analysis.fun, analysis.arg]:
    let token = symbolicCallableToken(input)
    if token != 0:
      tokens.add GradSigToken(kind: gstCallable, key: token)
  for dep in analysis.extraDeps:
    tokens.add GradSigToken(kind: gstInput, key: dep.stableNodeId)

proc appendApplySignature*(v: Gvalue, tokens: var seq[GradSigToken]) =
  v.analyzeApply.appendApplySignature(tokens)
