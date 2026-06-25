suite "functional reuse":
  test "source lambda remains reusable after returned lambda capture rebinding":
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

  test "local lambda refs reject lambda body rebinding":
    grt.resetApplyCache()
    grt.resetGradCache()

    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    let x = grt.toGvalue(3.0)
    let a = grt.toGvalue(2.0)
    let v = grt.localScalar()
    let z = apply(f, x)

    let linear = lambda(v, a * v + 1.0)
    let quadratic = lambda(v, a * v * v + 1.0)

    expect(GraphValueError):
      f.valCopy linear
    expect(GraphValueError):
      f.valCopy quadratic
    expect(GraphValueError):
      discard z.eval

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

  test "apply input-view trees do not instantiate lambda bodies":
    grt.resetApplyCache()

    let x = grt.toGvalue(3.0)
    let a = grt.toGvalue(2.0)
    let v = grt.localScalar()
    let f = lambda(v, a * v + 1.0)

    let z = apply(f, x)
    let grt = z.runtime
    check grt.functional.applyCacheStats.instantiationMisses == 0

    var sawCapture = false
    for dep in z.collectInputView(iwmReachable):
      if dep.nodeKey == a.nodeKey:
        sawCapture = true
    check sawCapture
    check grt.functional.applyCacheStats.instantiationMisses == 0

  test "apply backward input view exposes structural value targets without instantiation":
    grt.resetApplyCache()

    let x = grt.toGvalue(3.0)
    let scale = grt.toGvalue(2.0)
    let fParam = lambdaParam(grt.localScalar(), grt.localScalar())
    let u = grt.localScalar()
    let hof = lambda(fParam, Gscalar(apply(fParam, x)))
    let arg = lambda(u, scale * u)
    let call = apply(hof, arg)

    var sawX = false
    var sawScale = false
    for dep in call.collectInputView(iwmBackward):
      if dep.nodeKey == x.nodeKey:
        sawX = true
      if dep.nodeKey == scale.nodeKey:
        sawScale = true

    check sawX
    check sawScale
    check grt.functional.applyCacheStats.instantiationMisses == 0

  test "apply gradient build does not instantiate lambda bodies":
    grt.resetApplyCache()
    grt.resetGradCache()

    let x = grt.toGvalue(3.0)
    let a = grt.toGvalue(2.0)
    let v = grt.localScalar()
    let f = lambda(v, a * v + 1.0)

    let z = apply(f, x)
    discard z.grad a

    check grt.functional.applyCacheStats.instantiationMisses == 0

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
    let instantiationMisses = grt.functional.applyCacheStats.instantiationMisses
    check instantiationMisses > 0

    a.update 4.0
    z :~ 24.0
    dzda :~ 6.0
    let dzda1 = z.grad a
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

    grt.functional.applyCacheStats = ApplyCacheStats()
    grt.gradCacheStats = GradCacheStats()

    let dzdx1 = z.grad x
    let grt = z.runtime
    check cast[pointer](dzdx1) == cast[pointer](dzdx)
    check grt.gradCacheStats.directHits > 0

    grt.resetGradCache(stats = false)
    grt.functional.applyCacheStats = ApplyCacheStats()

    let dzdx2 = z.grad x
    discard dzdx2.eval
    dzdx2 :~ 6.0

  test "apply cache state is runtime local":
    let grt1 = initGraphRuntime()
    let grt2 = initGraphRuntime()
    grt1.resetApplyCache()
    grt2.resetApplyCache()

    let x1 = grt1.toGvalue(3.0)
    let v1 = grt1.localScalar()
    let z1 = apply(lambda(v1, v1 + 1.0), x1)
    z1 :~ 4.0

    check grt1.functional.applyCacheStats.instantiationMisses > 0
    check grt2.functional.applyCacheStats.instantiationMisses == 0

    grt1.functional.applyCacheStats = ApplyCacheStats()
    check grt1.functional.applyCacheStats.instantiationMisses == 0
    check grt2.functional.applyCacheStats.instantiationMisses == 0

  test "apply cache reset stays on the functional apply surface":
    grt.resetApplyCache()

    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let z = apply(lambda(v, v + 1.0), x)

    z :~ 4.0
    check grt.functional.applyCacheStats.instantiationMisses > 0

    grt.resetApplyCache()
    check grt.functional.applyCacheStats.instantiationHits == 0
    check grt.functional.applyCacheStats.instantiationMisses == 0
    x.update 5.0
    z :~ 6.0
    check grt.functional.applyCacheStats.instantiationHits == 0
    check grt.functional.applyCacheStats.instantiationMisses > 0

  test "apply instantiation cache reuses unchanged structure after input updates":
    grt.resetApplyCache()

    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let z = apply(lambda(v, v + 1.0), x)

    z :~ 4.0
    let instantiationMisses = grt.functional.applyCacheStats.instantiationMisses
    let instantiationHits = grt.functional.applyCacheStats.instantiationHits
    check instantiationMisses > 0

    x.update 4.0
    z :~ 5.0
    check grt.functional.applyCacheStats.instantiationMisses == instantiationMisses
    check grt.functional.applyCacheStats.instantiationHits > instantiationHits

  test "grad sees apply lambda capture deps before eval":
    let a = grt.toGvalue(2.0)

    block:
      let v = grt.localScalar()
      let dep = apply(lambda(v, a * v + 1.0), 3.0)

      let direct = dep.grad a
      dep :~ 7.0
      direct :~ 3.0

    a.update 4.0

    block:
      let v = grt.localScalar()
      let dep = apply(lambda(v, a * v + 1.0), 3.0)

      let direct = dep.grad a
      dep :~ 13.0
      direct :~ 3.0

  test "apply gradients stay stable across structural VJPs":
    let a = grt.toGvalue(2.0)
    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let f = lambda(v, a * v + 1.0)
    let z = Gscalar(apply(f, x))

    z.grad(a) :~ 3.0

    z.grad(x) :~ 2.0

    let vf = vjpOf(f)
    let dx = Gscalar(apply(apply(vf, x), 1.0))
    dx :~ 2.0
    dx.grad(a) :~ 1.0

  test "lambda body apply behavior stays stable across structural VJP builds":
    grt.resetApplyCache()
    grt.resetGradCache()

    let fParam = lambdaParam(grt.localScalar(), grt.localScalar())
    let x = grt.toGvalue(3.0)
    let sourceApply = Gscalar(apply(fParam, x))
    let hof = lambda(fParam, sourceApply)

    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let arg = lambda(v, a * v)
    let z = Gscalar(apply(hof, arg))

    z :~ 6.0
    z.grad(a) :~ 3.0
    grt.resetGradCache()
    z.grad(a) :~ 3.0

  test "reachable input view exposes apply capture deps without reducing":
    grt.resetApplyCache()

    let x = grt.toGvalue(3.0)
    let a = grt.toGvalue(2.0)
    let v = grt.localScalar()
    let f = lambda(v, a * v + 1.0)

    let z = apply(f, x)
    var sawCapture = false
    for dep in z.collectInputView(iwmReachable):
      if dep.nodeKey == a.nodeKey:
        sawCapture = true

    check sawCapture
    check grt.functional.applyCacheStats.instantiationMisses == 0

  test "shared lambda captures update through multiple applies":
    grt.resetApplyCache()

    let a = grt.toGvalue(2.0)
    let v = grt.localScalar()
    let u = grt.localScalar()
    let f = lambda(v, a * v)
    let g = lambda(u, a + u)

    let z = Gscalar(apply(f, 3.0)) + Gscalar(apply(g, 5.0))
    z :~ 13.0
    let instantiationMisses = grt.functional.applyCacheStats.instantiationMisses
    check instantiationMisses > 0

    a.update 4.0
    z :~ 21.0
    check grt.functional.applyCacheStats.instantiationMisses == instantiationMisses

  test "apply instantiation cache preserves capture value roles":
    grt.resetApplyCache()

    let a = grt.toGvalue(2.0)
    let b = grt.toGvalue(5.0)
    let v = grt.localScalar()
    let fn = lambda(v, a * v + b)
    let z = Gscalar(apply(fn, 3.0))

    z :~ 11.0
    a.update 4.0
    z :~ 17.0

  test "late-binding graph-produced lambda captures participate in apply gradients":
    grt.resetApplyCache()
    grt.resetGradCache()

    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    let x = grt.toGvalue(3.0)
    let z = Gscalar(apply(f, x))

    let p = lambdaParam(grt.localScalar(), grt.localScalar())
    let u = grt.localScalar()
    let hof = lambda(p, lambda(u, Gscalar(apply(p, u)) + 1.0))

    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let g = lambda(v, a * v)
    let h = apply(hof, g)

    f.valCopy(h)

    z :~ 7.0
    z.grad(a) :~ 3.0
    z.grad(x) :~ 2.0

    a.update 4.0
    z :~ 13.0
    z.grad(a) :~ 3.0
    z.grad(x) :~ 4.0

  test "late-bound local lambda eval follows capture updates":
    grt.resetApplyCache()
    grt.resetGradCache()

    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    let x = grt.toGvalue(3.0)
    let z = Gscalar(apply(f, x))

    let p = lambdaParam(grt.localScalar(), grt.localScalar())
    let u = grt.localScalar()
    let hof = lambda(p, lambda(u, Gscalar(apply(p, u)) + 1.0))
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let h = apply(hof, lambda(v, a * v))

    f.valCopy(h)
    z :~ 7.0

    let runs0 = z.runCount
    a.update 4.0
    z :~ 13.0
    check z.runCount > runs0

  test "late-bound lambda argument eval follows capture updates":
    grt.resetApplyCache()
    grt.resetGradCache()

    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    let x = grt.toGvalue(3.0)
    let p = lambdaParam(grt.localScalar(), grt.localScalar())
    let hof = lambda(p, Gscalar(apply(p, x)))
    let z = Gscalar(apply(hof, f))

    let q = lambdaParam(grt.localScalar(), grt.localScalar())
    let u = grt.localScalar()
    let maker = lambda(q, lambda(u, Gscalar(apply(q, u)) + 1.0))
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let h: GlambdaRef = GlambdaRef(apply(maker, lambda(v, a * v)))

    f.valCopy(h)
    z :~ 7.0

    let runs0 = z.runCount
    a.update 4.0
    z :~ 13.0
    check z.runCount > runs0

  test "late-binding lambda chooser exposes selected apply capture deps":
    grt.resetApplyCache()
    grt.resetGradCache()

    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    let x = grt.toGvalue(3.0)
    let z = Gscalar(apply(f, x))

    let k = grt.toGvalue(1)
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let b = grt.toGvalue(5.0)

    f.valCopy(cond(k, lambda(v, a * v), lambda(v, b * v * v)))

    z :~ 6.0
    z.grad(a) :~ 3.0
    z.grad(b) :~ 0.0

    k.update 0
    z :~ 45.0
    z.grad(a) :~ 0.0
    z.grad(b) :~ 9.0

  test "lambda cond freshness honors selector deps":
    grt.resetApplyCache()
    grt.resetGradCache()

    let selector = grt.toGvalue(0)
    let x = grt.localScalar()
    let identity = lambda(x, x)
    let square = lambda(x, x * x)
    let chooser = cond(equal(selector, 0), identity, square)
    let arg = grt.toGvalue(3.0)
    let z = apply(chooser, arg)

    z :~ 3.0
    let instantiationMisses = grt.functional.applyCacheStats.instantiationMisses
    check instantiationMisses > 0

    arg.update 4.0
    z :~ 4.0
    check grt.functional.applyCacheStats.instantiationMisses == instantiationMisses
    check grt.functional.applyCacheStats.instantiationHits > 0
    z.grad(arg) :~ 1.0

    selector.update 1

    z :~ 16.0
    let reboundInstantiationMisses = grt.functional.applyCacheStats.instantiationMisses
    check reboundInstantiationMisses > instantiationMisses
    let reboundInstantiationHits = grt.functional.applyCacheStats.instantiationHits

    arg.update 5.0
    z :~ 25.0
    check grt.functional.applyCacheStats.instantiationMisses >= reboundInstantiationMisses
    check grt.functional.applyCacheStats.instantiationHits >= reboundInstantiationHits
    z.grad(arg) :~ 10.0

  test "lambda cond exposes selected rebound captures to grad":
    grt.resetApplyCache()
    grt.resetGradCache()

    let selector = grt.toGvalue(0)
    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let b = grt.toGvalue(5.0)
    let linear = lambda(v, a * v)
    let quadratic = lambda(v, b * v * v)
    let chooser = cond(equal(selector, 0), linear, quadratic)
    let z = apply(chooser, x)

    z :~ 6.0
    let initialInstantiationMisses = grt.functional.applyCacheStats.instantiationMisses
    check initialInstantiationMisses == 1

    selector.update 1
    a.update 7.0

    z :~ 45.0
    let reboundInstantiationMisses = grt.functional.applyCacheStats.instantiationMisses
    check reboundInstantiationMisses == initialInstantiationMisses + 1

    z.grad(b) :~ 9.0
    z.grad(a) :~ 0.0
    check grt.functional.applyCacheStats.instantiationMisses >= reboundInstantiationMisses

    x.update 4.0
    z :~ 80.0
    z.grad(b) :~ 16.0
    z.grad(a) :~ 0.0
    check grt.functional.applyCacheStats.instantiationMisses >= reboundInstantiationMisses

  test "lambda cond apply gradients follow selected branch":
    grt.resetApplyCache()
    grt.resetGradCache()

    let k = grt.toGvalue(1)
    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let b = grt.toGvalue(5.0)
    let f = cond(k, lambda(v, a * v), lambda(v, b * v * v))
    let z = Gscalar(apply(f, x))

    discard f.eval
    check f.lambdaParamProto of Gscalar
    check f.lambdaResultProto of Gscalar

    z :~ 6.0
    z.grad(a) :~ 3.0
    z.grad(b) :~ 0.0

    k.update 0
    z :~ 45.0
    z.grad(a) :~ 0.0
    z.grad(b) :~ 9.0

  test "lambda cond apply VJP gradient node stays live across selector flips":
    grt.resetApplyCache()
    grt.resetGradCache()

    let k = grt.toGvalue(1)
    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let b = grt.toGvalue(5.0)
    let f = cond(k, lambda(v, a * v), lambda(v, b * v))
    let z = Gscalar(apply(f, x))
    let dzdx = z.grad x
    let p0 = cast[pointer](dzdx)

    dzdx :~ 2.0
    k.update 0
    dzdx :~ 5.0

    let dzdxAgain = z.grad x
    check cast[pointer](dzdxAgain) == p0
    check grt.gradCacheStats.directHits > 0
    check grt.gradCacheStats.invalidations == 0
    dzdxAgain :~ 5.0

  test "lambda apply gradients preserve deep branch body identity":
    grt.resetApplyCache()
    grt.resetGradCache()

    proc deepNeg(x: Gscalar): Gscalar =
      result = x
      for i in 0..<12:
        result = -result

    let k = grt.toGvalue(1)
    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let f = cond(k, lambda(v, deepNeg(v * v)), lambda(v, deepNeg(v + v)))
    let z = Gscalar(apply(f, x))
    let dzdx = z.grad x

    z :~ 9.0
    dzdx :~ 6.0

    k.update 0
    z :~ 6.0
    dzdx :~ 2.0

  test "lambda apply gradients distinguish revisited whole nodes from subnodes":
    grt.resetApplyCache()
    grt.resetGradCache()

    proc wholeReuse(v: Gscalar): Gscalar =
      let s = v * v
      let a = s * s
      a + a

    proc subnodeReuse(v: Gscalar): Gscalar =
      let s = v * v
      let a = s * s
      a + s

    let k = grt.toGvalue(1)
    let x = grt.toGvalue(2.0)
    let v = grt.localScalar()
    let f = cond(k, lambda(v, wholeReuse(v)), lambda(v, subnodeReuse(v)))
    let z = Gscalar(apply(f, x))
    let dzdx = z.grad x

    z :~ 32.0
    dzdx :~ 64.0

    k.update 0
    z :~ 20.0
    dzdx :~ 36.0

  test "lambda apply gradients separate same-named default-key custom funcs":
    grt.resetApplyCache()
    grt.resetGradCache()

    proc upstreamOrOne(zb: Gvalue, z: Gvalue): Gscalar =
      if zb == nil:
        return toGvalue(z.runtime, 1.0)
      Gscalar(zb)

    proc scaleTwoForward(node: Gvalue) =
      let x = Gscalar(node.inputs[0])
      Gscalar(node).sval = 2.0 * x.sval

    proc scaleTwoBackward(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
      if i != 0:
        raiseValueError("same-key scale two input index")
      upstreamOrOne(zb, z) * 2.0

    proc scaleFiveForward(node: Gvalue) =
      let x = Gscalar(node.inputs[0])
      Gscalar(node).sval = 5.0 * x.sval

    proc scaleFiveBackward(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
      if i != 0:
        raiseValueError("same-key scale five input index")
      upstreamOrOne(zb, z) * 5.0

    let scaleTwoFunc = Gfunc(
      forward: scaleTwoForward,
      backward: scaleTwoBackward,
      name: "same default-key scale")
    let scaleFiveFunc = Gfunc(
      forward: scaleFiveForward,
      backward: scaleFiveBackward,
      name: "same default-key scale")

    proc scaleTwo(x: Gscalar): Gscalar =
      graphNode(scalarNodeLike(x), @[Gvalue(x)], scaleTwoFunc, "same-key scale")

    proc scaleFive(x: Gscalar): Gscalar =
      graphNode(scalarNodeLike(x), @[Gvalue(x)], scaleFiveFunc, "same-key scale")

    let k = grt.toGvalue(1)
    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let f = cond(k, lambda(v, scaleTwo(v)), lambda(v, scaleFive(v)))
    let z = Gscalar(apply(f, x))
    let dzdx = z.grad x

    z :~ 6.0
    dzdx :~ 2.0

    k.update 0
    z :~ 15.0
    dzdx :~ 5.0

  test "apply keeps lambda-valued argument dependencies visible":
    grt.resetApplyCache()
    grt.resetGradCache()

    let fParam = lambdaParam(grt.localScalar(), grt.localScalar())
    let u = grt.localScalar()
    let hof = lambda(fParam, Gscalar(apply(fParam, u)))
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let b = grt.toGvalue(5.0)
    let k = grt.toGvalue(1)
    let arg = cond(k, lambda(v, a * v), lambda(v, b * v * v))
    let z = apply(hof, arg)
    let reachableDeps = z.collectInputView(iwmReachable)
    var seesArg = false
    var seesA = false
    var seesB = false
    for dep in reachableDeps:
      if dep.nodeKey == arg.nodeKey:
        seesArg = true
      if dep.nodeKey == a.nodeKey:
        seesA = true
      if dep.nodeKey == b.nodeKey:
        seesB = true
    check seesArg
    check seesA
    check seesB

  test "apply backward view keeps lambda argument captures, not lambda root":
    grt.resetApplyCache()
    grt.resetGradCache()

    let fParam = lambdaParam(grt.localScalar(), grt.localScalar())
    let u = grt.localScalar()
    let hof = lambda(fParam, Gscalar(apply(fParam, u)))
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let arg = lambda(v, a * v)
    let z = apply(hof, arg)
    let deps = z.collectInputView(iwmBackward)

    var seesArg = false
    var seesA = false
    for dep in deps:
      if dep.nodeKey == arg.nodeKey:
        seesArg = true
      if dep.nodeKey == a.nodeKey:
        seesA = true

    check not seesArg
    check seesA

  test "apply backward view deduplicates scalar arg captured by function":
    grt.resetApplyCache()
    grt.resetGradCache()

    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let f = lambda(v, x * v)
    let z = apply(f, x)
    let deps = z.collectInputView(iwmBackward)

    var xCount = 0
    for dep in deps:
      if dep.nodeKey == x.nodeKey:
        inc xCount

    check xCount == 1
    let dzdx = z.grad x
    let d2zdx2 = dzdx.grad x
    dzdx :~ 6.0
    d2zdx2 :~ 2.0

  test "scalar capture updates do not change symbolic revision":
    grt.resetApplyCache()
    grt.resetGradCache()

    let a = grt.toGvalue(2.0)
    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let z = apply(lambda(v, a * v + 1.0), x)
    let revision = grt.symbolicRevision

    a.update 5.0

    check grt.symbolicRevision == revision
    z :~ 16.0

  test "symbolic vjpOf construction does not mutate lambda revision":
    grt.resetApplyCache()
    grt.resetGradCache()

    let x = grt.toGvalue(3.0)
    let fParam = lambdaParam(grt.localScalar(), grt.localScalar())
    let hof = lambda(fParam, Gscalar(apply(fParam, x)))
    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    discard apply(hof, f)
    let revision = grt.symbolicRevision
    let epoch = f.epoch

    let vf = vjpOf(f)

    check vf.lambdaParamProto of Gscalar
    check f.epoch == epoch
    expect(GraphValueError):
      discard apply(f, 1.0).eval
    check grt.symbolicRevision == revision

  test "symbolic vjpOf construction keeps unrelated apply cache valid":
    grt.resetApplyCache()
    grt.resetGradCache()

    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let arg = lambda(v, a * v)
    let fParam = lambdaParam(grt.localScalar(), grt.localScalar())
    let hof = lambda(fParam, Gscalar(apply(fParam, x)))
    let z = Gscalar(apply(hof, arg))

    z :~ 6.0
    grt.functional.applyCacheStats = ApplyCacheStats()
    let revision = grt.symbolicRevision

    let gParam = lambdaParam(grt.localScalar(), grt.localScalar())
    let vg = vjpOf(gParam)
    let u = grt.localScalar()
    let metadataBuilder = lambda(gParam, Gscalar(apply(apply(vg, u), 1.0)))
    discard apply(metadataBuilder, arg)

    check grt.symbolicRevision == revision
    z :~ 6.0
    check grt.functional.applyCacheStats.instantiationMisses == 0

    x.update 4.0
    z :~ 8.0
    check grt.functional.applyCacheStats.instantiationMisses == 0

  test "symbolic vjpOf construction keeps unrelated grad cache valid":
    grt.resetApplyCache()
    grt.resetGradCache()

    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let arg = lambda(v, a * v)
    let fParam = lambdaParam(grt.localScalar(), grt.localScalar())
    let hof = lambda(fParam, Gscalar(apply(fParam, x)))
    let z = Gscalar(apply(hof, arg))

    z.grad(a) :~ 3.0
    grt.functional.applyCacheStats = ApplyCacheStats()
    grt.gradCacheStats = GradCacheStats()
    let revision = grt.symbolicRevision

    let gParam = lambdaParam(grt.localScalar(), grt.localScalar())
    let vg = vjpOf(gParam)
    let u = grt.localScalar()
    let metadataBuilder = lambda(gParam, Gscalar(apply(apply(vg, u), 1.0)))
    discard apply(metadataBuilder, arg)

    check grt.symbolicRevision == revision
    z.grad(a) :~ 3.0
    check grt.gradCacheStats.invalidations == 0

  test "apply cache distinguishes selected lambda identity":
    grt.resetApplyCache()
    grt.resetGradCache()

    let k = grt.toGvalue(1)
    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let f = cond(k, lambda(v, v + 1.0), lambda(v, v * v))
    let z = Gscalar(apply(f, x))

    z :~ 4.0
    grt.functional.applyCacheStats = ApplyCacheStats()

    k.update 0

    z :~ 9.0
    check grt.functional.applyCacheStats.instantiationMisses > 0

  test "apply VJP eval tree keeps apply boundary lazy":
    let x = grt.toGvalue(3.0)
    let a = grt.toGvalue(7.0)
    let v = grt.localScalar()
    let z = apply(lambda(v, a * v + v * v + 1.0), x)
    let dzdx = z.grad x
    let applyRuns0 = z.runCount

    dzdx :~ 13.0
    check z.runCount == applyRuns0

    a.update 11.0
    dzdx :~ 17.0
    check z.runCount == applyRuns0

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
    check grt.functional.applyCacheStats.instantiationMisses == 0

  test "lambda reuses instantiated work across forward and gradient eval":
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
      let instantiationMisses = grt.functional.applyCacheStats.instantiationMisses
      check instantiationMisses > 0
      let instantiationHits0 = grt.functional.applyCacheStats.instantiationHits

      discard dzdx.eval
      check z.runCount >= applyRuns0
      check e.runCount == 0
      check grt.functional.applyCacheStats.instantiationMisses >= instantiationMisses
      check grt.functional.applyCacheStats.instantiationHits >= instantiationHits0

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
      let instantiationMisses = grt.functional.applyCacheStats.instantiationMisses
      check instantiationMisses > 0
      let instantiationHits0 = grt.functional.applyCacheStats.instantiationHits

      discard z.eval
      check z.runCount == applyRuns0 + 1
      check e.runCount == 0
      check grt.functional.applyCacheStats.instantiationMisses >= instantiationMisses
      check grt.functional.applyCacheStats.instantiationHits >= instantiationHits0

  test "apply VJP stays live after input updates":
    grt.resetApplyCache()
    grt.resetGradCache()

    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let z = apply(lambda(v, v * v + 1.0), x)
    let dzdx = z.grad x

    dzdx :~ 6.0
    let dzdxPointer = cast[pointer](dzdx)

    x.update 4.0
    dzdx :~ 8.0
    let dzdxAgain = z.grad x
    check cast[pointer](dzdxAgain) == dzdxPointer

  test "whole lambda gradients are rejected":
    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let f = lambda(v, v * v + 1.0)
    let z = apply(f, x)
    let dzdx = z.grad x

    dzdx :~ 6.0
    expect(GraphValueError):
      discard dzdx.grad f

    let fParam = lambdaParam(grt.localScalar(), grt.localScalar())
    let u = grt.localScalar()
    let hof = lambda(fParam, lambda(u, Gscalar(apply(fParam, u)) + 1.0))
    let a = grt.toGvalue(2.0)
    let g = lambda(v, a * v)
    let h: GlambdaRef = GlambdaRef(apply(hof, g))

    expect(GraphValueError):
      discard h.grad(h)
    let p = grt.localScalar()
    let concrete = lambda(p, a * p)
    expect(GraphValueError):
      discard concrete.grad(concrete)

    # The unified guard covers grad over both Glambda and GlambdaRef and shares
    # one "not first-class" message.
    try:
      discard concrete.grad(concrete)
      fail()
    except GraphValueError as e:
      check "not first-class" in e.msg
    try:
      discard h.grad(h)
      fail()
    except GraphValueError as e:
      check "not first-class" in e.msg

  test "vjpOf builds ordinary structural VJP lambdas":
    let x = grt.toGvalue(3.0)
    let seed = grt.toGvalue(4.0)
    let v = grt.localScalar()
    let f = lambda(v, v * v)
    let vf = vjpOf(f)
    let dzdx = apply(apply(vf, x), seed)

    dzdx :~ 24.0

    x.update 5.0
    dzdx :~ 40.0

  test "vjpOf returned lambda keeps nested target deps live":
    let scale = grt.toGvalue(2.0)
    let x = grt.toGvalue(3.0)
    let y = grt.toGvalue(5.0)
    let seed = grt.toGvalue(4.0)
    let vx = grt.localScalar()
    let vy = grt.localScalar()
    let f = lambda(vx, lambda(vy, scale * vx * vy))
    let dx = apply(apply(apply(vjpOf(f), x), y), seed)

    dx :~ 40.0

    scale.update 3.0
    dx :~ 60.0

    y.update 7.0
    dx :~ 84.0

  test "independent structural VJP builds do not share active target deps":
    let fParam = lambdaParam(grt.localScalar(), grt.localScalar())
    let target = grt.toGvalue(3.0)
    let scale = grt.toGvalue(2.0)
    let hof = lambda(fParam, Gscalar(apply(fParam, target)))
    let v = grt.localScalar()
    let z = apply(hof, lambda(v, scale * v))
    let dzdtarget = z.grad target

    dzdtarget :~ 2.0

    let x = grt.toGvalue(4.0)
    let seed = grt.toGvalue(1.0)
    let w = grt.localScalar()
    let independent = apply(apply(vjpOf(lambda(w, w * w)), x), seed)

    independent :~ 8.0
    target.update 5.0
    independent :~ 8.0

  test "structural VJP gradients keep value targets reachable":
    let x = grt.toGvalue(3.0)
    let scale = grt.toGvalue(2.0)
    let fParam = lambdaParam(grt.localScalar(), grt.localScalar())
    let u = grt.localScalar()
    let hof = lambda(fParam, Gscalar(apply(fParam, x)))
    let arg = lambda(u, scale * u)
    let z = Gscalar(apply(hof, arg))
    let dzdscale = z.grad scale

    var seesScale = false
    for dep in dzdscale.collectInputView(iwmReachable):
      if dep.nodeKey == scale.nodeKey:
        seesScale = true
    check seesScale

    z :~ 6.0
    dzdscale :~ 3.0

    x.update 5.0
    z :~ 10.0
    dzdscale :~ 5.0

  test "structural VJP clone keeps active value targets visible":
    grt.resetApplyCache()
    grt.resetGradCache()

    proc graphSees(root: Gvalue,
                   target: Gvalue,
                   mode: InputWalkMode): bool =
      var seen: seq[Gvalue] = @[]
      var stack = @[root]
      while stack.len > 0:
        let node = stack[^1]
        stack.setLen(stack.len - 1)
        var alreadySeen = false
        for value in seen:
          if value.nodeKey == node.nodeKey:
            alreadySeen = true
            break
        if alreadySeen:
          continue
        seen.add node
        if node.nodeKey == target.nodeKey:
          return true
        node.walkInputView(mode, proc(child: Gvalue) =
          stack.add child)

    let a = grt.toGvalue(2.0)
    let x = grt.toGvalue(3.0)
    let scale = grt.toGvalue(5.0)
    let v = grt.localScalar()
    let capturedFun = lambda(v, a * a * v)
    let capturedApply = apply(capturedFun, x)

    let u = grt.localScalar()
    let outer = lambda(u, Gscalar(capturedApply) * u)
    let z = Gscalar(apply(outer, scale))
    let dzda = z.grad a
    let d2zda2 = dzda.grad a

    check graphSees(dzda, a, iwmReachable)
    check graphSees(dzda, a, iwmBackward)
    z :~ 60.0
    dzda :~ 60.0
    d2zda2 :~ 30.0

    block:
      let p = grt.localScalar()
      let q = grt.localScalar()
      let v = grt.localScalar()
      let source = lambda(v, p * v + q * v * v)
      let call = Gscalar(apply(source, p + q))
      let bothBars = call.grad(p) + call.grad(q)
      let d = lambda(p, lambda(q, bothBars))

      apply(apply(d, 2.0), 3.0) :~ 94.0
      apply(apply(d, 1.0), 4.0) :~ 112.0

  test "vjpOf of vjpOf differentiates as ordinary graph":
    let x = grt.toGvalue(3.0)
    let seed = grt.toGvalue(4.0)
    let upstream = grt.toGvalue(5.0)
    let v = grt.localScalar()
    let f = lambda(v, v * v)
    let second = apply(apply(apply(vjpOf(vjpOf(f)), x), seed), upstream)

    second :~ 40.0

    seed.update 7.0
    second :~ 70.0

  test "structural VJP clone keeps symbolic source unbound":
    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    let vf = vjpOf(f)
    let u = grt.localScalar()
    let hof = lambda(f, lambda(u, apply(apply(vf, u), 1.0)))
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let h = apply(hof, lambda(v, a * v * v))

    apply(h, 5.0) :~ 20.0
    expect(GraphValueError):
      discard apply(f, 1.0).eval

  test "structural VJP clone placeholders are independent per instantiation":
    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    let vf = vjpOf(f)
    let u = grt.localScalar()
    let hof = lambda(f, lambda(u, apply(apply(vf, u), 1.0)))
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let b = grt.toGvalue(3.0)
    let quad = lambda(v, a * v * v)
    let cubic = lambda(v, b * v * v * v)
    let dq = apply(apply(hof, quad), 5.0)
    let dc = apply(apply(hof, cubic), 5.0)

    dq :~ 20.0
    dc :~ 225.0

    a.update 4.0
    b.update 1.0
    dq :~ 40.0
    dc :~ 75.0

  test "structural VJP differentiates lambda argument captures":
    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    let u = grt.localScalar()
    let hof = lambda(f, lambda(u, Gscalar(apply(f, u)) + 1.0))
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let g = lambda(v, a * v)
    let h = apply(hof, g)

    let dhda = Gscalar(apply(h, 5.0)).grad a

    dhda :~ 5.0

    a.update 4.0
    dhda :~ 5.0

  test "lambda output gradSeeded rejects primal lambda seeds":
    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    let u = grt.localScalar()
    let hof = lambda(f, lambda(u, Gscalar(apply(f, u)) + 1.0))
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let g = lambda(v, a * v)
    let h = apply(hof, g)

    expect(GraphValueError):
      discard h.gradSeeded(a, h)

  test "lambda-valued apply root gradient is zero":
    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    let u = grt.localScalar()
    let hof = lambda(f, lambda(u, Gscalar(apply(f, u)) + 1.0))
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let g = lambda(v, a * v)
    let h = apply(hof, g)

    expect(GraphValueError):
      discard h.grad a

  test "graph-produced lambda refs are not zero cotangents":
    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    let u = grt.localScalar()
    let hof = lambda(f, lambda(u, Gscalar(apply(f, u)) + 1.0))
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let g = lambda(v, a * v)
    let h = apply(hof, g)

    check h of GlambdaRef
    check not h.isZero
    apply(h, 5.0) :~ 11.0

  test "lambda-valued cond root gradient through captures is rejected":
    let a = grt.toGvalue(2.0)
    let k = grt.toGvalue(1)
    let u = grt.localScalar()
    let v = grt.localScalar()

    let f1 = lambda(u, a * u + 1.0)
    let f2 = lambda(v, a + v)

    let h = cond(k, f1, f2)
    check h of GlambdaRef
    if h of GlambdaRef:
      discard h.eval
      check h.lambdaParamProto of Gscalar
      check h.lambdaResultProto of Gscalar
    apply(h, 3.0) :~ 7.0
    expect(GraphValueError):
      discard h.grad a

    k.update 0
    if h of GlambdaRef:
      discard h.eval
      check h.lambdaParamProto of Gscalar
      check h.lambdaResultProto of Gscalar
    apply(h, 3.0) :~ 5.0
    expect(GraphValueError):
      discard h.grad a

  test "apply VJP rejects lambda-like upstream scaling":
    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let z = apply(lambda(v, v * v + 1.0), x)

    proc copyForward(node: Gvalue) =
      node.valCopy node.inputs[0]

    proc lambdaUpstreamBackward(zb: Gvalue,
                                  node: Gvalue,
                                  i: int,
                                  input: Gvalue): Gvalue =
      discard zb
      discard node
      discard input
      if i != 0:
        raiseValueError("lambda-upstream test input index")
      lambdaParam(grt.localScalar(), grt.localScalar())

    let lambdaUpstreamFunc = Gfunc(
      forward: copyForward,
      backward: lambdaUpstreamBackward,
      name: "lambda upstream")
    let consumer = graphNode(
      scalarNodeLike(x),
      @[z],
      lambdaUpstreamFunc,
      "lambda upstream consumer")

    try:
      discard consumer.grad x
      check false
    except GraphValueError as e:
      check e.msg.contains("lambda upstream")

  test "cyclic lambda binding raises explicitly":
    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    f.valCopy(f)

    expect(GraphValueError):
      discard apply(f, 2.0).eval

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

  test "cond keeps apply VJPs live across selector flips":
    grt.resetApplyCache()
    grt.resetGradCache()

    let x = grt.toGvalue(3.0)
    let k = grt.toGvalue(1)
    let v = grt.localScalar()
    let z = cond(k, Gscalar(apply(lambda(v, v * v + 1.0), x)), 0.0)

    let dz = z.grad x
    let p0 = cast[pointer](dz)
    let grt = z.runtime
    check grt.functional.applyCacheStats.instantiationMisses == 0
    dz :~ 6.0

    k.update 0
    dz :~ 0.0
    let dzOff = z.grad x
    check cast[pointer](dzOff) == p0
    dzOff :~ 0.0

    k.update 1
    dz :~ 6.0
    let dzOn = z.grad x
    check cast[pointer](dzOn) == p0
    check grt.gradCacheStats.directHits > 0
    check grt.gradCacheStats.invalidations == 0

  test "inactive true cond apply VJP branch is not evaluated":
    grt.resetApplyCache()
    grt.resetGradCache()

    let x = grt.toGvalue(3.0)
    let k = grt.toGvalue(0)
    let v = grt.localScalar()
    let p = 1.0 / (x - 3.0)
    let z = cond(k, Gscalar(apply(lambda(v, v + 1.0), p)), 0.0)
    let dzdp = z.grad p

    dzdp :~ 0.0
    check p.runCount == 0

  test "inactive false cond apply VJP branch is not evaluated":
    grt.resetApplyCache()
    grt.resetGradCache()

    let x = grt.toGvalue(3.0)
    let k = grt.toGvalue(1)
    let v = grt.localScalar()
    let p = 1.0 / (x - 3.0)
    let z = cond(k, 0.0, Gscalar(apply(lambda(v, v + 1.0), p)))
    let dzdp = z.grad p

    dzdp :~ 0.0
    check p.runCount == 0

  test "inactive unresolved cond apply VJP branch resolves only when selected":
    grt.resetApplyCache()
    grt.resetGradCache()

    let x = grt.toGvalue(3.0)
    let k = grt.toGvalue(1)
    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    let z = cond(k, x * x, Gscalar(apply(f, x)))
    let dzdx = z.grad x

    dzdx :~ 6.0

    k.update 0
    expect(GraphValueError):
      discard dzdx.eval

  test "inactive unresolved lambda cond VJP branch resolves only when selected":
    grt.resetApplyCache()
    grt.resetGradCache()

    let x = grt.toGvalue(3.0)
    let k = grt.toGvalue(1)
    let v = grt.localScalar()
    let resolvedBody = lambda(v, v + 1.0)
    let unresolved = lambdaParam(grt.localScalar(), grt.localScalar())
    let f = cond(k, Gvalue(resolvedBody), Gvalue(unresolved))
    let z = Gscalar(apply(f, x))
    let dzdx = z.grad x

    dzdx :~ 1.0

    k.update 0
    expect(GraphValueError):
      discard dzdx.eval

  test "constant cond apply VJP does not evaluate unrelated selector":
    grt.resetApplyCache()
    grt.resetGradCache()

    let x = grt.toGvalue(3.0)
    let kSource = grt.toGvalue(1.0)
    let k = equal(kSource, 1.0)
    let trueBranch = grt.toGvalue(17.0)
    let falseBranch = grt.toGvalue(23.0)
    let v = grt.localScalar()
    let z = apply(lambda(v, cond(k, trueBranch, falseBranch)), x)

    let dzdx = z.grad x
    let kRuns0 = k.runCount
    dzdx :~ 0.0
    check k.runCount == kRuns0

  test "cond grad build does not evaluate either apply branch":
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
    check grt.functional.applyCacheStats.instantiationMisses == 0

  test "cond gradient eval materializes apply VJPs after branch expansion":
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
    check grt.functional.applyCacheStats.instantiationMisses == 0
    dzda :~ 3.0
    let firstInstantiationMisses = grt.functional.applyCacheStats.instantiationMisses

    k.update 0
    dzda :~ 1.0
    check grt.functional.applyCacheStats.instantiationMisses >= firstInstantiationMisses

  test "apply instantiations stay cached across deep capture value updates":
    grt.resetApplyCache()
    grt.resetGradCache()

    let protoArg = grt.localScalar()
    let protoRet = grt.localScalar()
    let fnProto = lambda(protoArg, protoRet)
    let x = lambdaParam(fnProto, fnProto)
    let f = lambdaParam(fnProto, fnProto)
    let Y = lambda(f, apply(lambda(x, apply(f, apply(x, x))), lambda(x, apply(f, apply(x, x)))))

    let rf = lambdaParam(grt.localScalar(), grt.localScalar())
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
    check grt.functional.applyCacheStats.instantiationMisses == 0
    dzdstep :~ 3.0
    let instantiationMisses = grt.functional.applyCacheStats.instantiationMisses
    check instantiationMisses > 0

    base.update 2.0
    let dzdbase = z.grad base
    dzdbase :~ 1.0
    z :~ 5.0
