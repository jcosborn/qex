suite "functional reuse":
  test "source lambda remains reusable after returned lambda reclosure":
    let outer = grt.localScalar()
    let inner = grt.localScalar()
    let a = grt.toGvalue(10.0)
    let makeAdder = lambda(outer, lambda(inner, outer + inner + a))

    let add1 = apply(makeAdder, 1.0)
    discard add1.eval
    let z1 = apply(add1, 3.0)
    z1 :~ 14.0

    let add2 = apply(makeAdder, 2.0)
    let z2 = apply(add2, 3.0)
    let dz2da = z2.grad a
    z2 :~ 15.0
    dz2da :~ 1.0

  test "apply and grad caches invalidate when callable rebinds":
    grt.resetApplyCache()
    grt.resetGradCache()

    let f = local(grt.localScalar())
    let x = grt.toGvalue(3.0)
    let a = grt.toGvalue(2.0)
    let v = grt.localScalar()
    let z = apply(f, x)

    let linear = lambda(v, a * v + 1.0)
    let quadratic = lambda(v, a * v * v + 1.0)

    f.valCopy linear
    z :~ 7.0
    let dz1 = z.grad a
    dz1 :~ 3.0
    let p1 = cast[pointer](dz1)

    f.valCopy quadratic
    z :~ 19.0
    let dz2 = z.grad a
    dz2 :~ 9.0

    let grt = z.runtime
    check cast[pointer](dz2) != p1
    check grt.gradCacheStats.invalidations > 0
    check grt.applyCacheStats.reduceMisses >= 2

  test "evaluated higher-order apply refreshes after callable rebind":
    grt.resetApplyCache()
    grt.resetGradCache()

    let f = local(grt.localScalar())
    let u = grt.localScalar()
    let hof = lambda(f, lambda(u, apply(f, u) + 1.0))
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let b = grt.toGvalue(3.0)
    let linear = lambda(v, a * v)
    let quadratic = lambda(v, b * v * v)

    f.valCopy linear
    let h = apply(hof, f)
    discard h.eval
    let z = apply(h, 3.0)

    let dz1da = z.grad a
    z :~ 7.0
    dz1da :~ 3.0

    f.valCopy quadratic
    let dz2db = z.grad b
    let dz2da = z.grad a
    z :~ 28.0
    dz2db :~ 9.0
    dz2da :~ 0.0

    let grt = z.runtime
    check grt.gradCacheStats.invalidations > 0
    check grt.applyCacheStats.reduceMisses >= 2

  test "reuse forward values and gradient node on stable graph":
    grt.resetApplyCache()
    grt.resetGradCache()

    let x = grt.toGvalue(3.0)
    let y = grt.toGvalue(2.0)
    let mul = x * y
    let v = grt.localScalar()
    let z = apply(lambda(v, v * v + v), mul)

    discard z.eval
    let applyRuns = z.runCount

    let dzdx1 = z.grad x
    dzdx1 :~ 26.0

    discard dzdx1.eval

    let dzdx2 = z.grad x
    let grt = z.runtime
    check cast[pointer](dzdx2) == cast[pointer](dzdx1)
    check grt.gradCacheStats.directHits > 0
    check grt.applyCacheStats.reduceHits > 0

    discard z.eval
    check z.runCount == applyRuns

  test "capture value updates keep gradient cache live":
    grt.resetGradCache()

    let a = grt.toGvalue(2.0)
    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let z = apply(lambda(v, a * v + 1.0), x)
    let dzda = z.grad a
    let p0 = cast[pointer](dzda)

    z :~ 7.0
    dzda :~ 3.0

    a.update 4.0
    z :~ 13.0
    dzda :~ 3.0

    let dzda1 = z.grad a
    let grt = z.runtime
    check cast[pointer](dzda1) == p0
    check grt.gradCacheStats.directHits > 0
    check grt.gradCacheStats.invalidations == 0

  test "cache stats reset keeps cached entries alive":
    grt.resetApplyCache()
    grt.resetGradCache()

    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let z = apply(lambda(v, v * v + 1.0), x)
    let dzdx = z.grad x

    z :~ 10.0
    dzdx :~ 6.0

    grt.resetApplyCacheStats()
    grt.resetGradCacheStats()

    let dzdx1 = z.grad x
    let grt = z.runtime
    check cast[pointer](dzdx1) == cast[pointer](dzdx)
    check grt.gradCacheStats.directHits > 0

    grt.clearGradCache()
    grt.resetApplyCacheStats()

    let dzdx2 = z.grad x
    discard dzdx2.eval
    check grt.applyCacheStats.reduceHits > 0

  test "grad sees apply closure deps before eval":
    let a = grt.toGvalue(2.0)

    block:
      let f = local(grt.localScalar())
      let dep = apply(f, 3.0)
      let v = grt.localScalar()
      let g = lambda(v, a * v + 1.0)
      f.valCopy g

      let direct = dep.grad a
      dep :~ 7.0
      direct :~ 3.0

    a.update 4.0

    block:
      let f = local(grt.localScalar())
      let dep = apply(f, 3.0)
      let v = grt.localScalar()
      let g = lambda(v, a * v + 1.0)
      f.valCopy g

      let direct = dep.grad a
      dep :~ 13.0
      direct :~ 3.0

  test "dependency tree exposes hidden apply deps":
    let f = local(grt.localScalar())
    let x = grt.toGvalue(3.0)
    let a = grt.toGvalue(2.0)
    let v = grt.localScalar()
    f.valCopy lambda(v, a * v + 1.0)

    let z = apply(f, x)
    let rawTree = z.treeRepr
    let depTree = z.treeRepr(iwmDepend)

    check depTree.splitLines.len > rawTree.splitLines.len

  test "signature tree exposes apply capture deps without reducing":
    let f = local(grt.localScalar())
    let x = grt.toGvalue(3.0)
    let a = grt.toGvalue(2.0)
    let v = grt.localScalar()
    f.valCopy lambda(v, a * v + 1.0)

    let z = apply(f, x)
    let rawTree = z.treeRepr
    let sigTree = z.treeRepr(iwmGradSignature)

    check not rawTree.contains("2.0")
    check sigTree.contains("2.0")

  test "callable dep walk reaches boundaries through ordinary nodes":
    proc opaqueForward(v: Gvalue) =
      discard v

    let opaqueWrap = newGfunc(
      forward = opaqueForward,
      name = "opaqueWrap")

    let a = grt.toGvalue(2.0)
    let f = local(grt.localScalar())
    let v = grt.localScalar()
    f.valCopy lambda(v, a * v)

    let opaque = graphNode(Gvalue(runtime: f.runtime), @[f], opaqueWrap)
    let deps = collectCallableValueDeps([opaque], [opaque])

    check deps.len == 1
    if deps.len > 0:
      check sameNode(deps[0], a)

  test "deferred apply partial eval tree skips intermediate apply node":
    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let z = apply(lambda(v, v * v + 1.0), x)
    let dzdx = z.grad x
    let dependTree = dzdx.treeRepr(iwmDepend)
    let evalTree = dzdx.treeRepr(iwmEval)

    check dependTree.contains("applyDeferred")
    check not evalTree.contains("applyDeferred")

  test "lambda grad build does not eval apply":
    grt.resetApplyCache()

    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let e = exp(v)
    let z = apply(lambda(v, e + v), x)
    let applyRuns0 = z.runCount
    let expRuns0 = e.runCount

    let dzdx = z.grad x
    discard dzdx.grad x

    let grt = z.runtime
    check z.runCount == applyRuns0
    check e.runCount == expRuns0
    check grt.applyCacheStats.reduceMisses == 0
    check grt.applyCacheStats.partialMisses == 0

  test "lambda reuses reduced work across forward and gradient eval":
    grt.resetApplyCache()
    grt.resetGradCache()

    block:
      let x = grt.toGvalue(3.0)
      let v = grt.localScalar()
      let e = exp(v)
      let z = apply(lambda(v, e + v), x)
      let dzdx = z.grad x
      let applyRuns0 = z.runCount
      let expRuns0 = e.runCount

      discard z.eval
      check z.runCount == applyRuns0 + 1
      let expRuns1 = e.runCount
      check expRuns1 > expRuns0

      discard dzdx.eval
      check z.runCount == applyRuns0 + 1
      check e.runCount == expRuns1

    block:
      let x = grt.toGvalue(3.0)
      let v = grt.localScalar()
      let e = exp(v)
      let z = apply(lambda(v, e + v), x)
      let dzdx = z.grad x
      let applyRuns0 = z.runCount
      let expRuns0 = e.runCount

      discard dzdx.eval
      check z.runCount == applyRuns0
      let expRuns1 = e.runCount
      check expRuns1 > expRuns0

      discard z.eval
      check z.runCount == applyRuns0 + 1
      check e.runCount == expRuns1

  test "apply partial depth guard raises without poisoning cache":
    grt.resetApplyCache()

    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let z = apply(lambda(v, v * v + 1.0), x)
    let dzdx = z.grad x
    let grt = z.runtime
    let savedLimit = grt.applyGradPrepareDepthLimit

    try:
      grt.applyGradPrepareDepthLimit = 0
      expect(GraphValueError):
        discard dzdx.eval

      grt.applyGradPrepareDepthLimit = savedLimit
      dzdx :~ 6.0
    finally:
      grt.applyGradPrepareDepthLimit = savedLimit

  test "callable resolution depth overflow raises explicitly":
    let savedLimit = grt.lambdaResolveDepthLimit

    try:
      grt.lambdaResolveDepthLimit = 0
      let v = grt.localScalar()
      expect(GraphValueError):
        discard apply(lambda(v, v + 1.0), 2.0)
    finally:
      grt.lambdaResolveDepthLimit = savedLimit

  test "cond gradients stay live without cache invalidation":
    grt.resetGradCache()

    let x = grt.toGvalue(2.0)
    let y = grt.toGvalue(3.0)
    let k = grt.toGvalue(0)
    let z = cond(k, x * x, y * y)

    let dzdx = z.grad x
    let p0 = cast[pointer](dzdx)
    dzdx :~ 0.0

    k.update 1
    dzdx :~ 4.0

    let dzdx1 = z.grad x
    let grt = z.runtime
    check cast[pointer](dzdx1) == p0
    check grt.gradCacheStats.directHits > 0
    check grt.gradCacheStats.invalidations == 0

  test "cond keeps apply partials live across selector flips":
    grt.resetApplyCache()
    grt.resetGradCache()

    let x = grt.toGvalue(3.0)
    let k = grt.toGvalue(1)
    let v = grt.localScalar()
    let z = cond(k, apply(lambda(v, v * v + 1.0), x), 0.0)

    let dz = z.grad x
    let p0 = cast[pointer](dz)
    let grt = z.runtime
    check grt.applyCacheStats.reduceMisses == 0
    check grt.applyCacheStats.partialMisses == 0
    dz :~ 6.0
    let partialMisses = grt.applyCacheStats.partialMisses
    check partialMisses > 0

    k.update 0
    dz :~ 0.0
    let dzOff = z.grad x
    check cast[pointer](dzOff) == p0
    dzOff :~ 0.0

    k.update 1
    dz :~ 6.0
    let dzOn = z.grad x
    check cast[pointer](dzOn) == p0
    check grt.applyCacheStats.partialMisses == partialMisses
    check grt.gradCacheStats.directHits > 0
    check grt.gradCacheStats.invalidations == 0

  test "cond grad build does not reduce either apply branch":
    grt.resetApplyCache()

    let a = grt.toGvalue(2.0)
    let x = grt.toGvalue(3.0)
    let k = grt.toGvalue(1)
    let v = grt.localScalar()
    let t = apply(lambda(v, a * v + 1.0), x)
    let f = apply(lambda(v, a + v * v), x)
    let z = cond(k, t, f)

    discard z.grad a
    discard z.grad x

    let grt = z.runtime
    check grt.applyCacheStats.reduceMisses == 0
    check grt.applyCacheStats.partialMisses == 0

  test "cond gradient eval reduces only the selected apply branch":
    grt.resetApplyCache()

    let a = grt.toGvalue(2.0)
    let x = grt.toGvalue(3.0)
    let k = grt.toGvalue(1)
    let v = grt.localScalar()
    let z = cond(k,
      apply(lambda(v, a * v + 1.0), x),
      apply(lambda(v, a + v * v), x))
    let dzda = z.grad a

    let grt = z.runtime
    check grt.applyCacheStats.reduceMisses == 0
    dzda :~ 3.0
    check grt.applyCacheStats.reduceMisses == 1

    k.update 0
    dzda :~ 1.0
    check grt.applyCacheStats.reduceMisses == 2

  test "apply reductions stay cached across deep capture value updates":
    grt.resetApplyCache()
    grt.resetGradCache()

    let protoArg = grt.localScalar()
    let protoRet = grt.localScalar()
    let fnProto = lambda(protoArg, protoRet)
    let x = local(fnProto)
    let f = local(fnProto)
    let Y = lambda(f, apply(lambda(x, apply(f, apply(x, x))), lambda(x, apply(f, apply(x, x)))))

    let rf = local(grt.localScalar())
    let u = grt.localScalar()
    let v = grt.localScalar()
    let base = grt.toGvalue(1.0)
    let step = grt.toGvalue(1.0)
    let F = lambda(rf, lambda(u,
      cond(equal(u, 0.0), base,
        apply(lambda(v, apply(rf, v) + step), u - 1.0))))

    let z = apply(apply(Y, F), 3.0)
    let dzdstep = z.grad step
    let grt = z.runtime
    check grt.applyCacheStats.reduceMisses == 0
    dzdstep :~ 3.0
    let reduceMisses = grt.applyCacheStats.reduceMisses
    check reduceMisses > 0

    base.update 2.0
    let dzdbase = z.grad base
    check grt.applyCacheStats.reduceMisses == reduceMisses
    dzdbase :~ 1.0
    check grt.applyCacheStats.reduceMisses == reduceMisses
    z :~ 5.0
