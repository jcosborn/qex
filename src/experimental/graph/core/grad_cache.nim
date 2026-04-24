import std/tables
import base

export GradCacheEntry, GradCacheStats

proc clearGradCache*(grt: GraphRuntime) =
  grt.gradCacheByOutput = initTable[NodeId, GradCacheEntry]()

proc resetGradCacheStats*(grt: GraphRuntime) =
  grt.gradCacheStats = GradCacheStats()

proc resetGradCache*(grt: GraphRuntime) =
  clearGradCache(grt)
  resetGradCacheStats(grt)

proc initGradCacheEntry(): GradCacheEntry =
  GradCacheEntry(grads: initTable[NodeId, Gvalue]())

proc ensureGradCacheEntry*(dep: Gvalue): GradCacheEntry =
  let grt = dep.runtime
  let depId = dep.stableNodeId
  if not grt.gradCacheByOutput.hasKey(depId):
    grt.gradCacheByOutput[depId] = initGradCacheEntry()
  result = grt.gradCacheByOutput[depId]

proc hasCachedGrad*(entry: GradCacheEntry, input: Gvalue): bool =
  entry != nil and input != nil and entry.grads.hasKey(input.stableNodeId)

proc cachedGrad*(entry: GradCacheEntry, input: Gvalue): Gvalue =
  entry.grads[input.stableNodeId]

proc cacheGrad*(entry: GradCacheEntry,
                input: Gvalue,
                grad: Gvalue) =
  entry.grads[input.stableNodeId] = grad

proc findGrad*(input: Gvalue, output: Gvalue): Gvalue =
  let grt = requireSameGraphRuntime(input, output, "findGrad", "input", "output")
  let outputId = output.stableNodeId
  if not grt.gradCacheByOutput.hasKey(outputId):
    return nil
  let entry = grt.gradCacheByOutput[outputId]
  if not entry.hasCachedGrad(input):
    return nil
  entry.cachedGrad(input)

proc dumpGradientList*(grt: GraphRuntime) =
  echo "# Gradient Cache:"
  for outputId, entry in grt.gradCacheByOutput.pairs:
    echo "## output id: ", outputId
    for inputId, grad in entry.grads.pairs:
      echo "### w.r.t. id: ", inputId
      echo grad.treeRepr

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
