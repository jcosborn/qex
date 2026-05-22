import std/tables
import base

export GradCacheEntry, GradCacheStats

proc ensureGradCacheEntry*(dep: Gvalue): GradCacheEntry =
  let grt = dep.runtime
  let depId = dep.stableNodeId
  if not grt.gradCacheByOutput.hasKey(depId):
    grt.gradCacheByOutput[depId] = GradCacheEntry(
      grads: initTable[NodeId, Gvalue]())
  grt.gradCacheByOutput[depId]

proc findGrad*(input: Gvalue, output: Gvalue): Gvalue =
  let grt = requireSameGraphRuntime(input, output, "findGrad", "input", "output")
  let outputId = output.stableNodeId
  if not grt.gradCacheByOutput.hasKey(outputId):
    return nil
  let entry = grt.gradCacheByOutput[outputId]
  let inputId = input.stableNodeId
  if not entry.grads.hasKey(inputId):
    return nil
  entry.grads[inputId]

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
