## Log-Jacobian chain builder.
##
## For a chain v = w0 -> w1 -> ... -> wn = u where each step's total derivative
## factors through its declared `via` (`Gfunc.logdet`),
##
##   ln|det(du/dv)| = sum_k ln|det(d w_k / d w_{k-1})| = sum_k ld_k.
##
## Determinants do not localize across merging paths: anything other than a
## via chain, with cond choosing between chains, is rejected. Full contract in
## DESIGN.md "Log-Jacobian factorization".

import std/[sets, tables]
import base, traverse, cond

proc requireSquare(z, via: Gvalue) =
  # dim z == dim via is required for det(d z/d via); shape compatibility is the
  # checkable proxy.
  if not (z.copyCompatible(via) and via.copyCompatible(z)):
    raiseValueError(
      "logDetJ step is not square (node and via shapes are incompatible):" &
      "\nnode: " & z.nodeRepr &
      "\nvia: " & via.nodeRepr)

proc requireCompatibleContributions(left, right: Gvalue) =
  ## Contribution algebra is erased in core, while the public scalar wrapper
  ## restores Gscalar. Check the erased boundary before zero elision or
  ## addLike can hide/miscast a malformed hook result.
  if not (left.copyCompatible(right) and right.copyCompatible(left)):
    raiseValueError(
      "logDetJ contributions have incompatible types or shapes:" &
      "\nleft: " & left.nodeRepr &
      "\nright: " & right.nodeRepr)

proc requireSeparates(z, via, base: Gvalue) =
  ## Every iwmBackward path from z to base must pass through via, and via must
  ## itself be a backward dependency of z; a bypass or a fabricated via makes
  ## the declared per-step determinant silently wrong.
  if z.nodeKey == via.nodeKey or not z.reaches(via, iwmBackward):
    raiseValueError(
      "logDetJ: declared via is not a backward dependency of the node:" &
      "\nnode: " & z.nodeRepr &
      "\nvia: " & via.nodeRepr)
  if z.reaches(base, iwmBackward, stop = via):
    raiseValueError(
      "logDetJ: dependence on the base bypasses the declared factorization:" &
      "\nnode: " & z.nodeRepr &
      "\nvia: " & via.nodeRepr &
      "\nbase: " & base.nodeRepr)

proc resetLdjCache*(grt: GraphRuntime) =
  grt.ldjCacheByNode = initTable[NodeId, LdjCacheEntry]()

proc currentEntry(z: Gvalue): LdjCacheEntry =
  let grt = z.runtime
  let id = z.stableNodeId
  result = grt.ldjCacheByNode.getOrDefault(id)
  if result == nil or result.revision != grt.symbolicRevision:
    result = LdjCacheEntry(
      revision: grt.symbolicRevision,
      sums: initTable[NodeId, Gvalue]())
    grt.ldjCacheByNode[id] = result

proc logDetJImpl*(u: Gvalue, v: Gvalue): Gvalue =
  ## Erased chain sum for typed wrappers: nil means the empty chain (u is v, or
  ## every cond branch reaches v directly). Fails when a chain node declares no
  ## factorization, a step is not square, or v is reachable around a via.
  discard sharedGraphRuntime([u, v], "logDetJ")
  let baseId = v.stableNodeId
  var active = initHashSet[NodeKey]()

  proc build(z: Gvalue): Gvalue =
    if z.nodeKey == v.nodeKey:
      return nil
    if z.nodeKey in active:
      raiseError("cycle detected in logDetJ chain:\n" & z.nodeRepr)
    let entry = currentEntry(z)
    if entry.sums.hasKey(baseId):
      return entry.sums[baseId]
    active.incl z.nodeKey
    defer: active.excl z.nodeKey
    var sum: Gvalue
    if z.isCondNode:
      # Each branch is its own chain.
      sum = distributeCond(z, build)
    else:
      if z.gfunc == nil:
        raiseError(
          "logDetJ: chain reached a leaf before the base " &
          "(wrong base or broken chain):" &
          "\nleaf: " & z.nodeRepr &
          "\nbase: " & v.nodeRepr)
      let hook = z.gfunc.logdet
      if hook == nil:
        raiseError(
          "logDetJ: no factorization declared on the chain " &
          "(missing Gfunc.logdet):" &
          "\nnode: " & z.nodeRepr &
          "\nbase: " & v.nodeRepr)
      if entry.ld == nil:
        let (ld, via) = hook(z)
        discard sharedGraphRuntime([z, ld, via], "logDetJ hook result")
        requireSquare(z, via)
        entry.ld = ld
        entry.via = via
      requireSeparates(z, entry.via, v)
      let rest = build(entry.via)
      if rest == nil:
        sum = entry.ld
      else:
        requireCompatibleContributions(entry.ld, rest)
        sum =
          if entry.ld.isStaticZeroLeaf: rest
          elif rest.isStaticZeroLeaf: entry.ld
          else: entry.ld.addLike(entry.ld, rest)
    # `hasKey` distinguishes a cached empty-chain nil from a missing sum.
    entry.sums[baseId] = sum
    sum

  build(u)
