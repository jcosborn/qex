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
    let linearReduceMisses = z.runtime.applyCacheStats.reduceMisses
    check linearReduceMisses >= 1
    let dz1 = z.grad a
    dz1 :~ 3.0
    let p1 = cast[pointer](dz1)

    f.valCopy quadratic
    z :~ 19.0
    check z.runtime.applyCacheStats.reduceMisses > linearReduceMisses
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
    let hof = lambda(f, lambda(u, Gscalar(apply(f, u)) + 1.0))
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

  test "apply argument gradient refreshes after callable rebind":
    grt.resetApplyCache()
    grt.resetGradCache()

    let f = local(grt.localScalar())
    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let z = apply(f, x)

    f.valCopy lambda(v, v + 1.0)
    let dzdx = z.grad x
    z :~ 4.0
    dzdx :~ 1.0

    f.valCopy lambda(v, v * v)
    z :~ 9.0
    dzdx :~ 6.0
    let dzdx2 = z.grad x
    dzdx2 :~ 6.0

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

  test "leaf callable wrapper apply follows direct binding":
    grt.resetApplyCache()

    let a = grt.toGvalue(2.0)
    let v = grt.localScalar()
    let fn = lambda(v, a * v + 1.0)
    let f = callableWrapperNode(grt.localScalar())

    f.valCopy fn
    let z = Gscalar(apply(f, 3.0))
    z :~ 7.0
    let reduceMisses = grt.applyCacheStats.reduceMisses
    check reduceMisses > 0

    a.update 4.0
    z :~ 13.0
    check grt.applyCacheStats.reduceMisses == reduceMisses

  test "apply dependency trees do not reduce callable bodies":
    grt.resetApplyCache()

    let f = local(grt.localScalar())
    let x = grt.toGvalue(3.0)
    let a = grt.toGvalue(2.0)
    let v = grt.localScalar()
    f.valCopy lambda(v, a * v + 1.0)

    let z = apply(f, x)
    let grt = z.runtime
    check grt.applyCacheStats.reduceMisses == 0
    check grt.applyCacheStats.partialMisses == 0

    let depTree = z.treeRepr(iwmDepend)
    let sigTree = z.treeRepr(iwmGradSignature)

    check depTree.contains("2.0")
    check sigTree.contains("2.0")
    check grt.applyCacheStats.reduceMisses == 0
    check grt.applyCacheStats.partialMisses == 0

  test "apply signature and grad planning do not reduce callable bodies":
    grt.resetApplyCache()
    grt.resetGradCache()

    let f = local(grt.localScalar())
    let x = grt.toGvalue(3.0)
    let a = grt.toGvalue(2.0)
    let v = grt.localScalar()
    f.valCopy lambda(v, a * v + 1.0)

    let z = apply(f, x)
    discard z.treeRepr(iwmGradSignature)
    discard z.grad a

    check grt.applyCacheStats.reduceMisses == 0
    check grt.applyCacheStats.partialMisses == 0

  test "duplicated capture updates keep apply and gradient caches live":
    grt.resetApplyCache()
    grt.resetGradCache()

    let a = grt.toGvalue(2.0)
    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let z = apply(lambda(v, a * v + a * v), x)
    let dzda = z.grad a
    let p0 = cast[pointer](dzda)

    z :~ 12.0
    dzda :~ 6.0
    let reduceMisses = grt.applyCacheStats.reduceMisses
    check reduceMisses > 0

    a.update 4.0
    z :~ 24.0
    dzda :~ 6.0
    let dzda1 = z.grad a
    check cast[pointer](dzda1) == p0
    check grt.applyCacheStats.reduceMisses == reduceMisses
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

  test "apply cache reset stays on the core runtime surface":
    grt.resetApplyCache()

    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let z = apply(lambda(v, v + 1.0), x)

    z :~ 4.0
    check grt.applyCacheStats.reduceMisses > 0

    grt.resetApplyCache()
    check grt.applyCacheStats.reduceHits == 0
    check grt.applyCacheStats.reduceMisses == 0
    x.update 5.0
    z :~ 6.0
    check grt.applyCacheStats.reduceHits == 0
    check grt.applyCacheStats.reduceMisses > 0

  test "apply reduction cache reuses unchanged structure after input updates":
    grt.resetApplyCache()

    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let z = apply(lambda(v, v + 1.0), x)

    z :~ 4.0
    let reduceMisses = grt.applyCacheStats.reduceMisses
    let reduceHits = grt.applyCacheStats.reduceHits
    check reduceMisses > 0

    x.update 4.0
    z :~ 5.0
    check grt.applyCacheStats.reduceMisses == reduceMisses
    check grt.applyCacheStats.reduceHits > reduceHits

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

  test "callable boundary through ordinary node stays observable via apply":
    proc opaqueForward(v: Gvalue) =
      v.valCopy v.inputs[0]

    proc reduceOpaque(v: Gvalue): Gvalue =
      v.inputs[0]

    let opaqueWrap = newGfunc(
      forward = opaqueForward,
      reduceCallable = reduceOpaque,
      name = "opaqueWrap")

    let a = grt.toGvalue(2.0)
    let f = local(grt.localScalar())
    let v = grt.localScalar()
    f.valCopy lambda(v, a * v)

    let opaque = graphNode(callableWrapperNode(grt.localScalar()), @[f], opaqueWrap)
    let z = Gscalar(apply(opaque, 3.0))

    check z.treeRepr(iwmDepend).contains("2.0")
    z :~ 6.0

    a.update 4.0
    z :~ 12.0

  test "deferred apply refreshes after local callable rebind":
    grt.resetApplyCache()
    let f = local(grt.localScalar())
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let b = grt.toGvalue(3.0)
    let x = grt.toGvalue(5.0)

    f.valCopy lambda(v, a * v)
    let z = Gscalar(apply(f, x))
    z :~ 10.0
    let misses = grt.applyCacheStats.reduceMisses
    check misses == 1

    f.valCopy lambda(v, b * v + 1.0)
    z :~ 16.0
    check grt.applyCacheStats.reduceMisses == misses + 1

  test "shared callable captures update through multiple applies":
    grt.resetApplyCache()

    let a = grt.toGvalue(2.0)
    let f = local(grt.localScalar())
    let g = local(grt.localScalar())
    let v = grt.localScalar()
    let u = grt.localScalar()

    f.valCopy lambda(v, a * v)
    g.valCopy lambda(u, a + u)

    let z = Gscalar(apply(f, 3.0)) + Gscalar(apply(g, 5.0))
    z :~ 13.0
    let reduceMisses = grt.applyCacheStats.reduceMisses
    check reduceMisses > 0

    a.update 4.0
    z :~ 21.0
    check grt.applyCacheStats.reduceMisses == reduceMisses

  test "callable boundary honors custom depend walks through apply":
    grt.resetApplyCache()

    let raw = grt.toGvalue(0.0)
    let f = local(grt.localScalar())
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let b = grt.toGvalue(3.0)
    f.valCopy lambda(v, a * v)

    proc copyForward(node: Gvalue) =
      node.valCopy f

    proc walkCustomDepend(node: Gvalue, visit: GnodeVisit) =
      discard node
      visit f

    let customCallable = newGfunc(
      forward = copyForward,
      depWalks = GdepWalks(
        eval: walkCustomDepend,
        gradSignature: walkCustomDepend,
        depend: walkCustomDepend),
      name = "custom callable boundary")
    let wrapper = graphNode(
      callableWrapperNode(grt.localScalar()),
      @[Gvalue(raw)],
      customCallable,
      "custom callable boundary")

    let z = Gscalar(apply(wrapper, 3.0))
    z :~ 6.0
    let reduceMisses = grt.applyCacheStats.reduceMisses
    check reduceMisses > 0

    f.valCopy lambda(v, b * v)
    z :~ 9.0
    check grt.applyCacheStats.reduceMisses == reduceMisses + 1

  test "callable boundary rejects nil custom depend deps through apply":
    let raw = grt.toGvalue(0.0)

    proc walkNilDepend(node: Gvalue, visit: GnodeVisit) =
      discard node
      let missing: Gvalue = nil
      visit missing

    let customCallable = newGfunc(
      depWalks = GdepWalks(depend: walkNilDepend),
      name = "nil callable dependency")
    let wrapper = graphNode(
      callableWrapperNode(grt.localScalar()),
      @[Gvalue(raw)],
      customCallable,
      "nil callable dependency")

    try:
      discard apply(wrapper, 3.0).treeRepr(iwmDepend)
      check false
    except GraphValueError as e:
      check e.msg.contains("nil dependency")

  test "apply reduction cache preserves capture binding roles":
    grt.resetApplyCache()

    let a = grt.toGvalue(2.0)
    let b = grt.toGvalue(5.0)
    let v = grt.localScalar()
    let fn = Glambda(lambda(v, a * v + b))
    let z = Gscalar(apply(fn, 3.0))

    z :~ 11.0

    var aIndex = -1
    var bIndex = -1
    for i in 0..<fn.env.len:
      if sameNode(fn.env[i].value, a):
        aIndex = i
      if sameNode(fn.env[i].value, b):
        bIndex = i
    check aIndex >= 0
    check bIndex >= 0

    fn.env[aIndex].value = b
    fn.env[bIndex].value = a
    fn.updated

    z :~ 17.0

  test "callable producer value updates refresh without re-reducing":
    grt.resetApplyCache()

    let producerDep = grt.toGvalue(0.0)
    let capture = grt.toGvalue(2.0)
    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let fn = lambda(v, capture * v)

    proc produceCallable(node: Gvalue) =
      node.valCopy fn

    let producerFunc = newGfunc(
      forward = produceCallable,
      name = "callable freshness direct producer")
    let producer = graphNode(
      callableWrapperNode(grt.localScalar()),
      @[Gvalue(producerDep)],
      producerFunc,
      "callable freshness direct producer")

    let z = Gscalar(apply(producer, x))
    z :~ 6.0
    let reduceMisses = grt.applyCacheStats.reduceMisses
    let reduceHits = grt.applyCacheStats.reduceHits

    capture.update 3.0
    producerDep.update 1.0
    z :~ 9.0
    check grt.applyCacheStats.reduceMisses == reduceMisses
    check grt.applyCacheStats.reduceHits > reduceHits

  test "callable producer capture updates stay visible without rebind":
    grt.resetApplyCache()
    grt.resetGradCache()

    let producerDep = grt.toGvalue(0.0)
    let capture = grt.toGvalue(2.0)
    let x = grt.toGvalue(4.0)
    let v = grt.localScalar()
    let fn = lambda(v, capture * v * v)

    proc produceCallable(node: Gvalue) =
      node.valCopy fn

    let producerFunc = newGfunc(
      forward = produceCallable,
      name = "callable freshness captured producer")
    let producer = graphNode(
      callableWrapperNode(grt.localScalar()),
      @[Gvalue(producerDep)],
      producerFunc,
      "callable freshness captured producer")

    let z = Gscalar(apply(producer, x))
    z :~ 32.0
    z.grad(x) :~ 16.0
    let producerRuns = producer.runCount
    let reduceMisses = grt.applyCacheStats.reduceMisses
    let reduceHits = grt.applyCacheStats.reduceHits

    capture.update 5.0
    check z.treeRepr(iwmDepend).contains("5.0")
    z :~ 80.0
    z.grad(x) :~ 40.0
    check producer.runCount == producerRuns
    check grt.applyCacheStats.reduceMisses == reduceMisses
    check grt.applyCacheStats.reduceHits > reduceHits

  test "callable freshness honors custom producer deps":
    grt.resetApplyCache()
    grt.resetGradCache()

    let selector = grt.toGvalue(0)
    let x = grt.localScalar()
    let identity = lambda(x, x)
    let square = lambda(x, x * x)

    proc chooseForward(v: Gvalue) =
      if selector.ival == 0:
        v.valCopy identity
      else:
        v.valCopy square

    proc walkSelector(node: Gvalue, visit: GnodeVisit) =
      discard node
      visit selector

    let chooseFunc = newGfunc(
      forward = chooseForward,
      depWalks = GdepWalks(
        eval: walkSelector,
        gradSignature: walkSelector,
        depend: walkSelector),
      name = "hidden callable chooser")
    let chooser = graphNode(
      callableWrapperNode(grt.localScalar()),
      newSeq[Gvalue](0),
      chooseFunc,
      "hidden callable chooser")
    let arg = grt.toGvalue(3.0)
    let z = apply(chooser, arg)

    z :~ 3.0
    let reduceMisses = grt.applyCacheStats.reduceMisses
    check reduceMisses > 0

    arg.update 4.0
    z :~ 4.0
    check grt.applyCacheStats.reduceMisses == reduceMisses
    check grt.applyCacheStats.reduceHits > 0
    z.grad(arg) :~ 1.0
    check chooser.resolveDirectLambda != nil

    selector.update 1
    check chooser.resolveDirectLambda == nil

    z :~ 16.0
    let reboundMisses = grt.applyCacheStats.reduceMisses
    check reboundMisses == reduceMisses + 1
    let reboundHits = grt.applyCacheStats.reduceHits

    arg.update 5.0
    z :~ 25.0
    check grt.applyCacheStats.reduceMisses == reboundMisses
    check grt.applyCacheStats.reduceHits > reboundHits
    z.grad(arg) :~ 10.0

  test "stale callable capture update does not affect rebound apply":
    grt.resetApplyCache()
    grt.resetGradCache()

    let selector = grt.toGvalue(0)
    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let b = grt.toGvalue(5.0)
    let linear = lambda(v, a * v)
    let quadratic = lambda(v, b * v * v)

    proc chooseForward(node: Gvalue) =
      if selector.ival == 0:
        node.valCopy linear
      else:
        node.valCopy quadratic

    proc walkSelector(node: Gvalue, visit: GnodeVisit) =
      discard node
      visit selector

    let chooseFunc = newGfunc(
      forward = chooseForward,
      depWalks = GdepWalks(
        eval: walkSelector,
        gradSignature: walkSelector,
        depend: walkSelector),
      name = "rebound stale capture chooser")
    let chooser = graphNode(
      callableWrapperNode(grt.localScalar()),
      newSeq[Gvalue](0),
      chooseFunc,
      "rebound stale capture chooser")
    let z = apply(chooser, x)

    z :~ 6.0
    let initialReduceMisses = grt.applyCacheStats.reduceMisses
    check initialReduceMisses == 1

    selector.update 1
    a.update 7.0

    z :~ 45.0
    let reboundReduceMisses = grt.applyCacheStats.reduceMisses
    check reboundReduceMisses == initialReduceMisses + 1

    z.grad(b) :~ 9.0
    z.grad(a) :~ 0.0
    check grt.applyCacheStats.reduceMisses == reboundReduceMisses

    x.update 4.0
    z :~ 80.0
    z.grad(b) :~ 16.0
    z.grad(a) :~ 0.0
    check grt.applyCacheStats.reduceMisses == reboundReduceMisses

  test "deferred apply partial eval tree skips intermediate apply node":
    let x = grt.toGvalue(3.0)
    let a = grt.toGvalue(7.0)
    let v = grt.localScalar()
    let z = apply(lambda(v, a * v + v * v + 1.0), x)
    let dzdx = z.grad x
    let dependTree = dzdx.treeRepr(iwmDepend)
    let evalTree = dzdx.treeRepr(iwmEval)

    check dependTree.contains("applyDeferred")
    check not evalTree.contains("applyDeferred")
    check evalTree.contains("7.0")
    dzdx :~ 13.0

    a.update 11.0
    dzdx :~ 17.0

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
    block:
      grt.resetApplyCache()
      grt.resetGradCache()

      let x = grt.toGvalue(3.0)
      let v = grt.localScalar()
      let e = exp(v)
      let z = apply(lambda(v, e + v), x)
      let dzdx = z.grad x
      let applyRuns0 = z.runCount

      discard z.eval
      check z.runCount == applyRuns0 + 1
      check e.runCount == 0
      let reduceMisses = grt.applyCacheStats.reduceMisses
      check reduceMisses > 0
      let reduceHits0 = grt.applyCacheStats.reduceHits

      discard dzdx.eval
      check z.runCount == applyRuns0 + 1
      check e.runCount == 0
      check grt.applyCacheStats.reduceMisses == reduceMisses
      check grt.applyCacheStats.reduceHits > reduceHits0

    block:
      grt.resetApplyCache()
      grt.resetGradCache()

      let x = grt.toGvalue(3.0)
      let v = grt.localScalar()
      let e = exp(v)
      let z = apply(lambda(v, e + v), x)
      let dzdx = z.grad x
      let applyRuns0 = z.runCount

      discard dzdx.eval
      check z.runCount == applyRuns0
      check e.runCount == 0
      let reduceMisses = grt.applyCacheStats.reduceMisses
      check reduceMisses > 0
      let reduceHits0 = grt.applyCacheStats.reduceHits

      discard z.eval
      check z.runCount == applyRuns0 + 1
      check e.runCount == 0
      check grt.applyCacheStats.reduceMisses == reduceMisses
      check grt.applyCacheStats.reduceHits > reduceHits0

  test "apply partial cache reuses materialized target partial":
    grt.resetApplyCache()
    grt.resetGradCache()

    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let z = apply(lambda(v, v * v + 1.0), x)
    let dzdx = z.grad x

    dzdx :~ 6.0
    let partialMisses = grt.applyCacheStats.partialMisses
    let partialHits = grt.applyCacheStats.partialHits
    check partialMisses > 0

    x.update 4.0
    dzdx :~ 8.0
    check grt.applyCacheStats.partialMisses == partialMisses
    check grt.applyCacheStats.partialHits > partialHits

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

  test "apply partial gradient treats callable-like target as zero":
    let f = local(grt.localScalar())
    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    f.valCopy lambda(v, v * v + 1.0)

    let z = apply(f, x)
    let dzdx = z.grad x
    let dzdxdf = dzdx.grad f

    dzdx :~ 6.0
    check dzdxdf of Gwrapper
    if dzdxdf of Gwrapper:
      check Gwrapper(dzdxdf).kind == wkLocal

  test "apply partial rejects callable-like upstream scaling":
    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let z = apply(lambda(v, v * v + 1.0), x)

    proc copyForward(node: Gvalue) =
      node.valCopy node.inputs[0]

    proc callableUpstreamBackward(zb: Gvalue,
                                  node: Gvalue,
                                  i: int,
                                  dep: Gvalue): Gvalue =
      discard zb
      discard node
      discard dep
      if i != 0:
        raiseValueError("callable-upstream test input index")
      local(grt.localScalar())

    let callableUpstreamFunc = newGfunc(
      forward = copyForward,
      backward = callableUpstreamBackward,
      name = "callable upstream")
    let consumer = graphNode(
      scalarNodeLike(x),
      @[z],
      callableUpstreamFunc,
      "callable upstream consumer")

    try:
      discard consumer.grad x
      check false
    except GraphValueError as e:
      check e.msg.contains("callable upstream")

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
    let z = cond(k, Gscalar(apply(lambda(v, v * v + 1.0), x)), 0.0)

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
        Gscalar(apply(lambda(v, Gscalar(apply(rf, v)) + step), u - 1.0)))))

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
