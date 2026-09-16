## Private execution graphs with reusable, typed value buffers.
## Build derivatives before planning and include every required result as a root.
import std/[tables, sets]
import core, functional, multi
import base/profile
when defined(graphPlanMemory):
  import base/alignedMem

type
  GraphPlanStats* = object
    runs*, forwards*, reuses*, sourceAudits*: int
    buffers*, arenaBytes*, peakLiveBytes*: int
    workspaces*, workspaceBytes*: int
  PlanState = enum psNew, psActive, psDone
  PlanNode = ref object
    value, shape: Gvalue
    deps, aliases, held, nested: seq[PlanNode]
    mode: GbufferMode
    external, protected, registered, dropped, shapeReady: bool
    state: PlanState
    refs, slot, pool: int
  PlanPool = object
    proto: Gvalue
    free: seq[int]
  PlanBuffer = object
    value: Gvalue
    owner: NodeKey
    fixed: bool
    bytes: int
    pool: int  # The allocation stays in its home pool across in-place transfers.
  GraphPlan* = ref object
    stats*: GraphPlanStats
    runtime: GraphRuntime
    sources, original, roots, results, copies, constants: seq[Gvalue]
    sourceKeys, overrides: NodeSet
    nodes: seq[PlanNode]
    byNode: NodeTable[PlanNode]
    instantiations: ApplyFrameCache
    workspaces: GworkCache
    buffers: seq[PlanBuffer]
    pools: seq[PlanPool]
    active: seq[PlanNode]
    feeds: seq[Gvalue]
    feedKeys: NodeSet
    versions, resultVersions, constantVersions: seq[int]
    revision, auditedRevision: uint64
    compiled, valid, running: bool
    liveBytes, opaqueDepth: int

when defined(graphPlanMemory):
  proc reportMemory*(p: GraphPlan, stage: string) =
    ## Diagnostic collections belong outside evaluation timers. Counts reuse
    ## owned metadata and do not allocate another traversal of the graph.
    var ext, deps, aliases, held, nested, ready, shapes: int
    for n in p.nodes:
      if n.external: inc ext
      deps += n.deps.len
      aliases += n.aliases.len
      held += n.held.len
      nested += n.nested.len
      # A completed shape query can return nil for unsupported storage.
      if n.shapeReady: inc ready
      if n.shape != nil: inc shapes
    let
      occupied = getOccupiedMem()
      total = getTotalMem()
      raw = getRawMemUsed()
    GC_fullCollect()
    let
      collected = getOccupiedMem()
      collectedTotal = getTotalMem()
      collectedRaw = getRawMemUsed()
    echo "graph-plan-memory stage=", stage,
      " assignedIds=", p.runtime.nextStableNodeId,
      " reachableSources=", p.sources.len, " sourceKeys=", p.sourceKeys.len,
      " planRecords=", p.nodes.len, " nodePoolBytes=", p.nodes.len * sizeof(int),
      " externalRecords=", ext,
      " dependencyEdges=", deps, " aliasEdges=", aliases,
      " heldEdges=", held, " nestedEdges=", nested,
      " initializedShapes=", ready, " shapeValues=", shapes,
      " applyVariants=", p.instantiations.len,
      " runtimeApplyEntries=", p.runtime.functional.applyCacheByNode.len,
      " runtimeRunStats=", p.runtime.runStatsByNode.len,
      " occupiedBeforeGc=", occupied, " occupiedAfterGc=", collected,
      " heapBeforeGc=", total, " heapAfterGc=", collectedTotal,
      " rawUsedBeforeGc=", raw, " rawUsedAfterGc=", collectedRaw,
      " rawAllocated=", getRawMemAllocated(),
      " arenaBytes=", p.stats.arenaBytes, " peakLiveBytes=", p.stats.peakLiveBytes,
      " workspaceBytes=", p.stats.workspaceBytes,
      " workspaces=", p.workspaces.len, " buffers=", p.buffers.len,
      " pools=", p.pools.len,
      " forwards=", p.stats.forwards, " reuses=", p.stats.reuses,
      " sourceAudits=", p.stats.sourceAudits

proc eval*(p: GraphPlan): seq[Gvalue] {.discardable.}
proc ensure(p: GraphPlan, n: PlanNode)
proc activate(p: GraphPlan, n: PlanNode)
proc drop(p: GraphPlan, n: PlanNode)

proc numeric(v: Gvalue): bool =
  if v.lambdaResultProto != nil:
    return false
  if v of Gmulti:
    let m = Gmulti(v)
    if m.isStructuralValue:
      return false
    for i in 0..<m.len:
      if not numeric(m.storedSlot(i)):
        return false
  true

proc overridden(v: Gvalue): bool =
  # Structural function values also use updated() for bindings. They must be
  # cloned, rather than mistaken for externally supplied numerical overrides.
  v.gfunc != nil and v.valueOverride and v.numeric

proc protect(p: GraphPlan, n: PlanNode) =
  if n.protected or n.external:
    return
  n.protected = true
  if n.slot >= 0:
    # A dynamically discovered opaque consumer may retain this descriptor.
    # Convert its live allocation into dedicated storage for this private node.
    p.buffers[n.slot].fixed = true
    if p.runtime.graphDebug:
      echo "[graph/plan] dedicate ", n.value.nodeRepr
  for x in n.deps:
    p.protect(x)

proc node(p: GraphPlan, v: Gvalue): PlanNode =
  let key = v.nodeKey
  if p.byNode.hasKey(key):
    return p.byNode[key]
  if v.runtime != p.runtime:
    raiseValueError("storage plan dependency mixes graph runtimes")
  let leaf = v.gfunc == nil and v.inputs.len == 0 and
    not v.staticZeroLeaf and v.restoreValue == nil and v.numeric
  # An initialized numerical leaf has no writer with which to reconstruct a
  # replacement, including leaves introduced by a nested forward.
  result = PlanNode(value: v, slot: -1, pool: -1, external: key in p.sourceKeys or leaf)
  p.byNode[key] = result
  p.nodes.add result
  if result.external:
    return
  result.mode = if v.gfunc == nil:
      (if v.staticZeroLeaf: bmZero else: bmFull)
    else: v.gfunc.bufferMode
  if v.gfunc != nil and v.gfunc.inplace.len > 0 and result.mode != bmFull:
    raiseValueError("in-place storage requires a full-write forward: " & v.gfunc.name)
  # Raw edges include fast-path operands omitted from the reachable view.
  # Most nodes have one or two inputs; allocate a set only for wide input lists.
  let n = result
  let wide = v.inputs.len > 16
  var seen: NodeSet
  if wide: seen = initHashSet[NodeKey](v.inputs.len)
  template addChild(child: Gvalue) =
    if wide:
      if seen.markSeenNode(child): n.deps.add p.node(child)
    else:
      var found = false
      for old in n.deps:
        if old.value.nodeKey == child.nodeKey:
          found = true
          break
      if not found: n.deps.add p.node(child)
  for x in v.inputs:
    addChild(x)
  if v.gfunc != nil and v.gfunc.inputView != nil:
    v.walkInputView(iwmReachable, proc(x: Gvalue) =
      addChild(x))
  if result.mode == bmAlias:
    let indices = v.gfunc.aliasInputs
    for i, x in v.inputs:
      if indices.len == 0 or i in indices:
        let a = p.node(x)
        if a notin result.aliases: result.aliases.add a
  if result.mode == bmOpaque:
    for x in result.deps: p.protect(x)
  if p.opaqueDepth > 0:
    p.protect(result)

proc activate(p: GraphPlan, n: PlanNode) =
  if n.external:
    n.registered = true
    return
  if n.registered and not (n.dropped and n.state == psNew):
    return
  n.registered = true
  n.dropped = false
  for x in n.deps:
    p.activate(x)
    inc x.refs

proc unbind(p: GraphPlan, n: PlanNode) =
  if n.slot < 0 or p.buffers[n.slot].fixed:
    return
  let i = n.slot
  n.value.bindBuffer(n.shape)
  n.value.valueReady = false
  n.slot = -1
  if p.buffers[i].owner == n.value.nodeKey:
    p.buffers[i].owner = nil
    p.liveBytes -= p.buffers[i].bytes
    p.pools[p.buffers[i].pool].free.add i

proc releaseInputs(p: GraphPlan, n: PlanNode, keepAliases: bool) =
  if n.dropped:
    return
  n.dropped = true
  for x in n.deps:
    if not keepAliases or x notin n.held:
      p.drop(x)

proc drop(p: GraphPlan, n: PlanNode) =
  dec n.refs
  if n.refs > 0 or n.external:
    return
  if n.refs < 0:
    raiseError("storage plan released an input twice")
  if n.state == psNew:
    # A branch which no longer has a consumer never executes. Cancel its
    # dependency edges as well, releasing common subexpressions when possible.
    p.releaseInputs(n, false)
  elif n.state == psDone:
    if n.mode == bmAlias:
      for x in n.held: p.drop(x)
      n.value.valueReady = false
    p.unbind(n)

proc compatible(a, b: Gvalue): bool =
  a.bufferCompatible(b) and b.bufferCompatible(a)

proc prepare(n: PlanNode) =
  n.value.ensureStorage
  if n.mode == bmZero and n.value.bufferBytes > 0:
    n.value.clearBuffer

proc acquire(p: GraphPlan, n: PlanNode) =
  let v = n.value
  if n.protected or n.mode == bmOpaque:
    n.prepare
    return
  if not n.shapeReady:
    n.shape = v.bufferProto
    n.shapeReady = true
    if n.shape != nil:
      for i in 0..<p.pools.len:
        if compatible(n.shape, p.pools[i].proto):
          n.pool = i
          break
      if n.pool < 0:
        n.pool = p.pools.len
        p.pools.add PlanPool(proto: n.shape)
  if n.pool < 0:
    n.prepare
    return
  var slot = -1
  var inplace = false
  if v.gfunc != nil:
    for i in v.gfunc.inplace:
      let x = p.byNode[v.inputs[i].nodeKey]
      if x.slot >= 0 and x.refs == 1 and not x.protected and
         x notin p.active and not p.buffers[x.slot].fixed and
         p.buffers[x.slot].owner == x.value.nodeKey and
         compatible(n.shape, p.buffers[x.slot].value):
        slot = x.slot
        inplace = true
        break
  if slot < 0:
    if p.pools[n.pool].free.len > 0:
      let i = p.pools[n.pool].free[^1]
      if compatible(n.shape, p.buffers[i].value):
        slot = p.pools[n.pool].free.pop
      else:
        # Compatibility with a prototype need not extend to concrete storage.
        # Leave the rejected slot at home and dedicate this node for this generation.
        n.pool = -1
        n.prepare
        return
  if slot < 0:
    # All members come from the same prototype; requesters need not be compatible.
    let storage = p.pools[n.pool].proto.bufferProto
    if storage != nil: storage.ensureStorage
    if storage == nil or not compatible(n.shape, storage):
      n.pool = -1
      n.prepare
      return
    slot = p.buffers.len
    let bytes = storage.bufferBytes
    p.buffers.add PlanBuffer(value: storage, bytes: bytes, pool: n.pool)
    p.stats.buffers = p.buffers.len
    p.stats.arenaBytes += bytes
  else:
    inc p.stats.reuses
  if not inplace:
    p.liveBytes += p.buffers[slot].bytes
    p.stats.peakLiveBytes = max(p.stats.peakLiveBytes, p.liveBytes)
  p.buffers[slot].owner = v.nodeKey
  n.slot = slot
  v.bindBuffer(p.buffers[slot].value)
  n.prepare

proc externalValue(p: GraphPlan, n: PlanNode) =
  let v = n.value
  # Explicitly overridden source expressions remain externally owned boundaries.
  # Their normal evaluator decides when a newer input supersedes the override.
  let hook = p.runtime.evalFrame
  p.runtime.evalFrame = nil
  try:
    discard v.eval
  finally:
    p.runtime.evalFrame = hook
  if p.feedKeys.markSeenNode(v):
    p.feeds.add v
  n.state = psDone

proc ensure(p: GraphPlan, n: PlanNode) =
  let v = n.value
  if n.external:
    if n.state != psDone: p.externalValue(n)
    return
  if n.state == psActive:
    raiseError("cycle detected in storage plan:\n" & v.nodeRepr)
  if n.state == psDone:
    if not v.valueReady or not v.hasStorage:
      raiseError("storage plan dependency requested after its last use:\n" & v.nodeRepr &
        "\nDeclare every numerical dependency and plan all required outputs together.")
    return
  n.state = psActive
  p.active.add n
  var maxep = v.epoch
  let rawAliases = n.mode == bmAlias and v.gfunc.inputView == nil
  if n.mode == bmAlias:
    n.held.setLen(0)
    if rawAliases:
      n.held.setLen(n.aliases.len)
      for i, x in n.aliases: n.held[i] = x
  var success = false
  try:
    v.walkInputView(iwmEval, proc(x: Gvalue) =
      let child = p.node(x)
      if v.gfunc != nil and v.gfunc.inputView != nil and child notin n.deps:
        raiseError("evaluation dependency is absent from raw/reachable inputs:\n" & v.nodeRepr)
      p.activate(child)
      p.ensure(child)
      if n.mode == bmAlias and not rawAliases and child in n.aliases and child notin n.held:
        n.held.add child
      maxep = max(maxep, x.epoch))
    v.valueReady = false
    if n.mode != bmAlias:
      p.acquire(n)
    elif not n.shapeReady:
      n.shape = v.bufferProto
      n.shapeReady = true
    let f = v.gfunc
    if f != nil and f.forward != nil:
      if n.mode == bmOpaque: inc p.opaqueDepth
      let t = getTics()
      try:
        f.forward(v)
      finally:
        let secs = ticDiffSecs(getTics(), t)
        inc p.stats.forwards
        p.runtime.runStatsByNode.mgetOrPut(v.stableNodeId, RunStat()).record(secs, f.name)
        if n.mode == bmOpaque: dec p.opaqueDepth
      v.valueOverride = false
    elif v.restoreValue != nil:
      v.restoreValue(v)
    elif n.deps.len > 0:
      raiseError("storage plan value has dependencies but no forward:\n" & v.nodeRepr)
    v.epoch = max(maxep, v.epoch)
    v.valueReady = true
    n.state = psDone
    success = true
  finally:
    for x in n.nested:
      p.drop(x)
    n.nested.setLen(0)
    discard p.active.pop
    if success:
      p.releaseInputs(n, n.mode == bmAlias)

proc nestedEval(p: GraphPlan, v: Gvalue) =
  if p.active.len == 0:
    raiseError("nested evaluation has no active storage-plan forward")
  let n = p.node(v)
  p.activate(n)
  inc n.refs
  p.active[^1].nested.add n
  p.ensure(n)

proc copiedOutput(p: GraphPlan, n: PlanNode, seen: var NodeSet): bool =
  if n.external: return true
  if not seen.markSeenNode(n.value): return false
  if n.mode in {bmOpaque, bmAlias}:
    let deps = if n.mode == bmAlias: n.held else: n.deps
    for x in deps:
      if p.copiedOutput(x, seen): return true
  false

proc detach(p: GraphPlan) =
  for n in p.nodes:
    p.unbind(n)
    if not n.external and n.mode == bmAlias and n.shape != nil:
      n.value.bindBuffer(n.shape)
      n.value.valueReady = false
  p.liveBytes = 0
  for pool in p.pools.mitems: pool.free.setLen(0)
  for i in 0..<p.buffers.len:
    if not p.buffers[i].fixed:
      p.buffers[i].owner = nil
      p.pools[p.buffers[i].pool].free.add i

proc clear*(p: GraphPlan) =
  ## Drop private execution caches and the arena. Published results retain their
  ## backing allocations; subsequent execution creates a new private generation.
  if p.running:
    raiseValueError("cannot clear a running storage plan")
  for n in p.nodes:
    if not n.external:
      n.value.releaseWork
  p.instantiations.clear
  p.workspaces.clear
  p.detach
  p.nodes.setLen(0)
  p.byNode.clear
  p.buffers.setLen(0)
  p.pools.setLen(0)
  p.roots.setLen(0)
  p.feeds.setLen(0)
  p.versions.setLen(0)
  p.stats.buffers = 0
  p.stats.arenaBytes = 0
  p.stats.peakLiveBytes = 0
  p.stats.workspaces = 0
  p.stats.workspaceBytes = 0
  p.compiled = false
  p.valid = false

proc rebuild(p: GraphPlan) =
  p.clear
  p.sources = graphValues(p.original)
  when defined(graphPlanMemory):
    p.reportMemory("traversal-complete")
  p.sourceKeys = initHashSet[NodeKey]()
  p.overrides = initHashSet[NodeKey]()
  p.constants.setLen(0)
  p.constantVersions.setLen(0)
  let rev = p.runtime.boundaryRevision
  inc p.stats.sourceAudits
  var preserve: seq[Gvalue]
  for v in p.sources:
    if v.gfunc == nil and (v.staticZeroLeaf or v.restoreValue != nil):
      p.constants.add v
      p.constantVersions.add v.epoch
    # Only preserved overrides and ordinary leaves can remain original objects.
    if v.overridden:
      p.sourceKeys.incl v.nodeKey
      p.overrides.incl v.nodeKey
      preserve.add v
    elif v.gfunc == nil and v.inputs.len == 0 and
        not v.staticZeroLeaf and v.restoreValue == nil:
      p.sourceKeys.incl v.nodeKey
  p.roots = cloneValues(p.original, preserve, copyConstants = true)
  when defined(graphPlanMemory):
    p.reportMemory("clone-complete")
  # Retired wrappers are never published again; callers reacquire p[i].
  # Readiness also invalidates dependent plans' cached feed versions.
  for v in p.results:
    v.stale = true
    v.valueReady = false
  p.results.setLen(0)
  p.copies.setLen(0)
  for v in p.roots:
    discard p.node(v)
    let pub = v.valueLike
    pub.stale = true
    pub.valueReady = false
    p.results.add pub
    p.copies.add v.valueLike
  p.revision = p.runtime.symbolicRevision
  # Cloning can advance the runtime; only the source scan was audited.
  p.auditedRevision = rev
  p.compiled = true
  when defined(graphPlanMemory):
    p.reportMemory("plan-records-complete")

proc plan*(roots: varargs[Gvalue]): GraphPlan =
  if roots.len == 0:
    raiseValueError("storage plan requires at least one result")
  let rt = sharedGraphRuntime(roots, "storage plan")
  for v in roots:
    if not numeric(v):
      raiseValueError("storage plans publish numerical values; apply function-valued results first")
  result = GraphPlan(runtime: rt, original: @roots,
    byNode: initTable[NodeKey, PlanNode](), sourceKeys: initHashSet[NodeKey](),
    overrides: initHashSet[NodeKey](), feedKeys: initHashSet[NodeKey](),
    instantiations: newTable[ApplyFrameKey, ApplyCacheEntry](),
    workspaces: newTable[GworkKey, Gwork]())
  result.rebuild

proc len*(p: GraphPlan): int = p.results.len
proc `[]`*(p: GraphPlan, i: int): Gvalue = p.results[i]

proc current(p: GraphPlan): bool =
  if not p.valid:
    return false
  for i, v in p.feeds:
    if v.gfunc != nil:
      discard v.eval
    if v.epoch != p.versions[i] or not v.valueReady or not v.hasStorage:
      if p.runtime.graphDebug:
        echo "[graph/plan] stale feed ", v.nodeRepr, " saved=", p.versions[i],
          " ready=", v.valueReady, " storage=", v.hasStorage
      return false
  for i, v in p.results:
    if v.epoch != p.resultVersions[i] or not v.valueReady or not v.hasStorage:
      if p.runtime.graphDebug:
        echo "[graph/plan] stale result ", v.nodeRepr, " saved=", p.resultVersions[i],
          " ready=", v.valueReady, " storage=", v.hasStorage
      return false
  true

proc refresh(p: GraphPlan): bool =
  if not p.compiled or p.revision != p.runtime.symbolicRevision:
    if p.runtime.graphDebug:
      echo "[graph/plan] rebuild compiled=", p.compiled, " revision=", p.revision,
        "/", p.runtime.symbolicRevision
    p.rebuild
    return true
  if p.auditedRevision == p.runtime.boundaryRevision:
    return false
  let rev = p.runtime.boundaryRevision
  inc p.stats.sourceAudits
  var overrides = initHashSet[NodeKey]()
  var constantsChanged = false
  var num = 0
  for v in p.sources:
    if v.overridden: overrides.incl v.nodeKey
    if v.gfunc == nil and (v.staticZeroLeaf or v.restoreValue != nil):
      if num >= p.constants.len or v.nodeKey != p.constants[num].nodeKey or
          v.epoch != p.constantVersions[num]:
        constantsChanged = true
      inc num
  if num != p.constants.len: constantsChanged = true
  if overrides != p.overrides or constantsChanged:
    if p.runtime.graphDebug:
      echo "[graph/plan] rebuild overrides=", overrides != p.overrides,
        " constants=", constantsChanged
    p.rebuild
    return true
  p.auditedRevision = rev
  false

proc run(p: GraphPlan): seq[Gvalue] =
  inc p.stats.runs
  var rebuilt = p.refresh
  let fresh = p.current
  # Validating an externally overridden expression can discover a newer input
  # or missing storage and end that override. Audit its boundary before reuse.
  if p.refresh: rebuilt = true
  if fresh and not rebuilt:
    return p.results
  p.valid = false
  p.running = true
  p.liveBytes = 0
  for buffer in p.buffers:
    if buffer.fixed: p.liveBytes += buffer.bytes
  p.stats.peakLiveBytes = max(p.stats.peakLiveBytes, p.liveBytes)
  p.opaqueDepth = 0
  p.feeds.setLen(0)
  p.feedKeys.clear
  for v in p.results:
    v.stale = true
    v.valueReady = false
  for n in p.nodes:
    n.state = psNew
    n.refs = 0
    n.registered = false
    n.dropped = false
    n.nested.setLen(0)
    n.held.setLen(0)
  for v in p.roots:
    let n = p.node(v)
    p.activate(n)
    inc n.refs
  p.runtime.evalFrame = proc(v: Gvalue) = p.nestedEval(v)
  p.runtime.applyFrame = p.instantiations
  p.runtime.workFrame = p.workspaces
  try:
    for v in p.roots: p.ensure(p.node(v))
    for i, v in p.roots:
      var seen = initHashSet[NodeKey]()
      if p.copiedOutput(p.node(v), seen):
        p.copies[i].ensureStorage
        p.copies[i].valCopy(v)
        p.results[i].valAlias(p.copies[i])
      else:
        p.results[i].valAlias(v)
      p.results[i].updated
    p.versions.setLen(0)
    for v in p.feeds: p.versions.add v.epoch
    p.resultVersions.setLen(0)
    for v in p.results: p.resultVersions.add v.epoch
    p.revision = p.runtime.symbolicRevision
    # A cold externalValue can expire an override after validation was skipped.
    # Leave its boundary revision pending for the next source audit.
    p.valid = true
  finally:
    p.runtime.evalFrame = nil
    p.runtime.applyFrame = nil
    p.runtime.workFrame = nil
    p.active.setLen(0)
    p.opaqueDepth = 0
    p.detach
    p.running = false
    p.stats.workspaces = p.workspaces.len
    p.stats.workspaceBytes = 0
    for work in p.workspaces.values: p.stats.workspaceBytes += work.bytes
  p.results

proc eval*(p: GraphPlan): seq[Gvalue] {.discardable.} =
  # Reject nested calls before changing previously published results.
  if p.running or p.runtime.evalFrame != nil:
    raiseValueError("storage plans execute serially; use graph apply for nested functions")
  try:
    result = p.run
  except:
    p.valid = false
    # Validation and partial publication can fail outside the execution loop.
    for v in p.results:
      v.stale = true
      v.valueReady = false
    raise
