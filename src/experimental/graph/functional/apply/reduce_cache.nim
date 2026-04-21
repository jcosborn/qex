import std/tables
import ../../core
import ../closure/normalize, meta, partial

export ApplyCacheStats

proc clearApplyCache*(grt: GraphRuntime) =
  grt.applyCacheByNode = initTable[NodeId, ApplyCacheEntry]()

proc resetApplyCacheStats*(grt: GraphRuntime) =
  grt.applyCacheStats = ApplyCacheStats()

proc resetApplyCache*(grt: GraphRuntime) =
  clearApplyCache(grt)
  resetApplyCacheStats(grt)

proc applyCacheEntry(v: Gvalue): ApplyCacheEntry =
  let grt = v.runtime
  let nodeId = v.stableNodeId
  if not grt.applyCacheByNode.hasKey(nodeId):
    grt.applyCacheByNode[nodeId] = ApplyCacheEntry()
  grt.applyCacheByNode[nodeId]

proc applyEntryMatches(entry: ApplyCacheEntry,
                       analysis: ApplyAnalysis): bool =
  if not entry.hasKey or entry.reduced == nil or analysis.fn == nil:
    return false
  entry.key.callableKey == analysis.callableKey and
    entry.key.inputKeys == analysis.inputKeys

proc setApplyReduction(entry: ApplyCacheEntry,
                       analysis: ApplyAnalysis,
                       reduced: Gvalue) =
  entry.hasKey = true
  entry.key = ApplyCacheKey(
    callableKey: analysis.callableKey,
    inputKeys: analysis.inputKeys)
  entry.reduced = reduced
  entry.partials = initTable[NodeId, Gvalue]()

proc ensureReduction(v: Gvalue): ApplyCacheEntry =
  var analysis = v.analyzeApply
  analysis.prepareApplyReduction
  if analysis.fn == nil:
    raiseUnresolvedValueError("deferred apply unresolved at eval: " & analysis.fun.nodeRepr)

  let grt = v.runtime
  result = v.applyCacheEntry
  if result.applyEntryMatches(analysis):
    grt.applyCacheStats.reduceHits.inc
    return
  grt.applyCacheStats.reduceMisses.inc

  let reduced = instantiateNormalizedBody(analysis.fn, analysis.arg)
  result.setApplyReduction(analysis, reduced)

proc requireReduction(v: Gvalue): ApplyCacheEntry =
  result = v.ensureReduction
  if result.reduced == nil:
    raiseValueError("deferred apply preparation did not produce reduced body")

proc reducedApplyExpr*(v: Gvalue): Gvalue =
  let entry = v.requireReduction
  entry.reduced

proc lookupApplyPartial(grt: GraphRuntime,
                        entry: ApplyCacheEntry,
                        target: Gvalue): tuple[found: bool, value: Gvalue] =
  let targetId = target.stableNodeId
  if entry.partials.hasKey(targetId):
    grt.applyCacheStats.partialHits.inc
    return (true, entry.partials[targetId])
  grt.applyCacheStats.partialMisses.inc
  (false, nil)

proc storeApplyPartial(entry: ApplyCacheEntry,
                       target: Gvalue,
                       value: Gvalue): Gvalue =
  entry.partials[target.stableNodeId] = value
  value

proc materializeConcreteApplyPartial(z: Gvalue,
                                     reduced: Gvalue,
                                     target: Gvalue): Gvalue =
  let grt = z.runtime
  inc grt.applyGradPrepareDepth
  defer:
    dec grt.applyGradPrepareDepth
  if grt.applyGradPrepareDepth > grt.applyGradPrepareDepthLimit:
    raiseValueError(
      "apply partial materialization exceeded depth limit " &
      $grt.applyGradPrepareDepthLimit &
      "\napply: " & z.nodeRepr &
      "\ntarget: " & target.nodeRepr)
  reduced.gradOrZero(target)

proc materializeApplyPartial(z: Gvalue,
                             reduced: Gvalue,
                             target: Gvalue): Gvalue =
  if not target.needsConcreteApplyPartial:
    return target.zeroLike
  z.materializeConcreteApplyPartial(reduced, target)

proc resolveApplyPartial(z: Gvalue,
                         entry: ApplyCacheEntry,
                         target: Gvalue): Gvalue =
  let cached = z.runtime.lookupApplyPartial(entry, target)
  if cached.found:
    return cached.value
  let partial = z.materializeApplyPartial(entry.reduced, target)
  entry.storeApplyPartial(target, partial)

proc ensurePartial*(z: Gvalue, target: Gvalue): Gvalue =
  let requiredTarget = target.requireApplyPartialTarget
  let grt = z.runtime
  let targetGrt = requiredTarget.runtime
  # Apply-partial cache keys are runtime-local stable ids.
  if grt != targetGrt:
    raiseValueError("apply partial mixes multiple graph runtimes")
  if z.isApplyPartialNode:
    let view = z.requireApplyPartialView
    return view.base.ensurePartial(view.target).gradOrZero(requiredTarget)
  let entry = z.requireReduction
  z.resolveApplyPartial(entry, requiredTarget)

proc reduceApplyPartialExpr*(v: Gvalue): Gvalue =
  let view = v.requireApplyPartialView
  view.base.ensurePartial(view.target)
