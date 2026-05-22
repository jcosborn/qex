import std/[sets, tables]
import base

export GradCacheEntry, GradCacheStats

proc resetGradCacheEntry*(entry: GradCacheEntry, signature: GradSignature) =
  entry.hasSignature = true
  entry.signature = signature
  entry.adjoints = initTable[NodeId, Gvalue]()
  entry.completeAdjoints = initHashSet[NodeId]()
  entry.inputContributions = initTable[NodeId, Table[int, Gvalue]]()
  entry.expandedInputs = initTable[NodeId, HashSet[int]]()
  entry.expandedTargets = initTable[NodeId, HashSet[NodeId]]()

proc ensureGradCacheEntry*(dep: Gvalue): GradCacheEntry =
  let grt = dep.runtime
  let depId = dep.stableNodeId
  if not grt.gradCacheByOutput.hasKey(depId):
    grt.gradCacheByOutput[depId] = GradCacheEntry(
      adjoints: initTable[NodeId, Gvalue](),
      completeAdjoints: initHashSet[NodeId](),
      inputContributions: initTable[NodeId, Table[int, Gvalue]](),
      expandedInputs: initTable[NodeId, HashSet[int]](),
      expandedTargets: initTable[NodeId, HashSet[NodeId]]())
  grt.gradCacheByOutput[depId]

proc findGrad*(input: Gvalue, output: Gvalue): Gvalue =
  let grt = requireSameGraphRuntime(input, output, "findGrad", "input", "output")
  let outputId = output.stableNodeId
  if not grt.gradCacheByOutput.hasKey(outputId):
    return nil
  let entry = grt.gradCacheByOutput[outputId]
  let inputId = input.stableNodeId
  if inputId notin entry.completeAdjoints:
    return nil
  if not entry.adjoints.hasKey(inputId):
    return nil
  entry.adjoints[inputId]

proc withIsolatedGradCache*[T](grt: GraphRuntime,
                               body: proc(): T {.closure.}): T =
  let savedGradCacheByOutput = grt.gradCacheByOutput
  let savedGradCacheStats = grt.gradCacheStats
  resetGradCache(grt)
  try:
    result = body()
  finally:
    grt.gradCacheByOutput = savedGradCacheByOutput
    grt.gradCacheStats = savedGradCacheStats
