suite "functional lambda":
  test "apply scalar and grad":
    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let z = apply(lambda(v, v + v), x * x)
    let dzdx = z.grad x
    z :~ 18.0
    dzdx :~ 12.0

  test "lambda capture by reference":
    let x = grt.toGvalue(2.0)
    let y = grt.toGvalue(4.0)
    let v = grt.localScalar()
    let z = apply(lambda(v, v * y + y), x)
    let dzdy = z.grad y
    z :~ 12.0
    dzdy :~ 3.0
    y.update 5.0
    z :~ 15.0
    dzdy :~ 3.0

  test "lambda normalization keeps captures paired with inputs":
    let y = grt.toGvalue(4.0)
    let v = grt.localScalar()
    let f = lambda(v, v * y + y)
    let z = apply(f, 3.0)

    check f.lambdaParamProto of Gscalar
    check f.lambdaResultProto of Gscalar
    z :~ 16.0

    y.update 5.0
    z :~ 20.0

  test "lambda normalization keeps capture substitution live":
    let a = grt.toGvalue(2.0)
    let b = grt.toGvalue(5.0)
    let v = grt.localScalar()
    let f = lambda(v, a * v + b)
    let z = apply(f, 3.0)

    z :~ 11.0

    a.update 4.0
    b.update 1.0
    z :~ 13.0

  test "returned lambda keeps multiple captures paired after instantiation":
    let a = grt.toGvalue(2.0)
    let b = grt.toGvalue(5.0)
    let outer = grt.localScalar()
    let inner = grt.localScalar()
    let f = lambda(outer, lambda(inner, a * outer + b * inner))
    let h = apply(f, 3.0)
    let z = apply(h, 4.0)
    let dzda = z.grad a
    let dzdb = z.grad b

    z :~ 26.0
    dzda :~ 3.0
    dzdb :~ 4.0

    a.update 7.0
    b.update 11.0
    z :~ 65.0
    dzda :~ 3.0
    dzdb :~ 4.0

  test "lambda normalization follows produced lambda ref deps":
    let a = grt.toGvalue(2.0)
    let outer = grt.localScalar()
    let inner = grt.localScalar()
    let maker = lambda(outer, lambda(inner, a * outer + inner))
    let produced = apply(maker, 3.0)
    let p = grt.localScalar()
    let wrapper = lambda(p, apply(produced, p))
    let z = apply(wrapper, 4.0)

    z :~ 10.0

    a.update 5.0
    z :~ 19.0

  test "multiple applies in one graph":
    let x = grt.toGvalue(2.0)
    let y = grt.toGvalue(7.0)
    let v = grt.localScalar()
    let f = lambda(v, v + y)
    let z = Gscalar(apply(f, x)) * Gscalar(apply(f, x + 1.0))
    let dzdy = z.grad y
    z :~ 90.0
    dzdy :~ 19.0

  test "higher order function argument":
    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    let u = grt.localScalar()
    let hof = lambda(f, lambda(u, Gscalar(apply(f, u)) + 1.0))
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let g = lambda(v, a * v)
    let z = apply(apply(hof, g), 3.0)
    let dzda = z.grad a
    z :~ 7.0
    dzda :~ 3.0
    a.update 4.0
    z :~ 13.0
    dzda :~ 3.0

  test "evaluated function-valued apply remains lambda":
    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    let u = grt.localScalar()
    let hof = lambda(f, lambda(u, Gscalar(apply(f, u)) + 1.0))
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let g = lambda(v, a * v)
    let h = apply(hof, g)
    check h of GlambdaRef
    check h.lambdaParamProto of Gscalar
    check h.lambdaResultProto of Gscalar
    discard h.eval

    let z = apply(h, 3.0)
    let dzda = z.grad a
    z :~ 7.0
    dzda :~ 3.0

  test "evaluated function-valued apply tracks lambda capture updates":
    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    let u = grt.localScalar()
    let hof = lambda(f, lambda(u, Gscalar(apply(f, u)) + 1.0))
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let g = lambda(v, a * v)
    let h = apply(hof, g)
    discard h.eval

    let z = apply(h, 3.0)
    z :~ 7.0

    a.update 4.0
    z :~ 13.0

  test "evaluated lambda refreshes only when captured values change":
    grt.resetApplyCache()

    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    let u = grt.localScalar()
    let hof = lambda(f, lambda(u, Gscalar(apply(f, u)) + 1.0))
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let g = lambda(v, a * v)
    let h = apply(hof, g)

    discard h.eval
    let applyRuns0 = h.runCount
    discard h.eval
    check h.runCount == applyRuns0

    a.update 4.0
    discard h.eval
    check h.runCount == applyRuns0 + 1

    let z = apply(h, 3.0)
    z :~ 13.0

  test "lambda ref valCopy rejects lambda bodies":
    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    let a = grt.toGvalue(2.0)
    let v = grt.localScalar()
    let epoch0 = f.epoch

    expect(GraphValueError):
      f.valCopy lambda(v, a * v + 1.0)
    expect(GraphValueError):
      discard apply(f, 1.0).eval
    check f.epoch == epoch0

  test "lambda placeholder rejects lambda and non-lambda binding":
    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    check not compiles(f.valCopy 3.0)
    check not compiles(f.valCopy 3)

    let v = grt.localScalar()
    let valid = lambda(v, v + 1.0)
    let epoch = f.epoch

    expect(GraphValueError):
      f.valCopy valid
    expect(GraphValueError):
      f.valCopy(grt.toGvalue(3.0))

    expect(GraphValueError):
      discard apply(f, 1.0).eval
    check f.epoch == epoch

  test "lambda parameter requires explicit result prototype":
    check not compiles(lambdaParam(grt.localScalar()))

  test "lambda parameter placeholder requires runtime-bearing prototype":
    check not compiles(lambdaParam(grt))
    let rawPrototype = Gvalue()
    let nilPrototype: Gvalue = nil
    expect(GraphValueError):
      discard lambdaParam(rawPrototype, grt.localScalar())
    expect(GraphValueError):
      discard lambdaParam(grt.localScalar(), rawPrototype)
    expect(GraphValueError):
      discard lambdaParam(rawPrototype, rawPrototype)
    expect(GraphValueError):
      discard lambdaParam(nilPrototype, grt.localScalar())
    expect(GraphValueError):
      discard lambdaParam(grt.localScalar(), nilPrototype)

  test "lambda placeholders require result prototypes before binding":
    let malformed = GlambdaRef(runtime: grt).assignStableNodeId
    let v = grt.localScalar()

    expect(GraphValueError):
      malformed.valCopy lambda(v, v + 1.0)

  test "malformed lambda refs reject public prototype probes":
    let malformed = GlambdaRef(runtime: grt).assignStableNodeId

    expect(GraphValueError):
      discard malformed.lambdaParamProto
    expect(GraphValueError):
      discard malformed.lambdaResultProto

  test "malformed lambda refs are incompatible":
    let malformed = GlambdaRef(runtime: grt).assignStableNodeId
    let valid = lambdaParam(grt.localScalar(), grt.localScalar())
    let v = grt.localScalar()
    let body = lambda(v, v + 1.0)

    check not malformed.copyCompatible(valid)
    check not valid.copyCompatible(malformed)
    check not malformed.copyCompatible(body)
    check not body.copyCompatible(malformed)

  test "malformed lambda refs fail early as lambda-shaped values":
    let malformed = GlambdaRef(runtime: grt).assignStableNodeId
    let x = grt.localScalar()

    expect(GraphValueError):
      discard malformed.zeroLike
    expect(GraphValueError):
      discard grad(grt.toGvalue(1.0), malformed)
    expect(GraphValueError):
      discard vjpOf(malformed)
    expect(GraphValueError):
      discard apply(malformed, x)

  test "apply materializes nested lambda results":
    let finalProto = grt.localScalar()
    let innerParam = grt.localScalar()
    let innerFnProto = lambda(innerParam, finalProto)
    let outerParam = grt.localScalar()
    let outerFn = lambda(outerParam, innerFnProto)

    let h = apply(outerFn, 2.0)
    check h of GlambdaRef
    check h.lambdaParamProto of Gscalar
    check h.lambdaResultProto of Gscalar

  test "apply accepts deeper nested lambda results":
    let innerParam = grt.localScalar()
    let innerFnProto = lambda(innerParam, innerParam)
    let middleParam = grt.localScalar()
    let middleFnProto = lambda(middleParam, innerFnProto)
    let outerParam = grt.localScalar()
    let outerFn = lambda(outerParam, middleFnProto)

    let h = apply(outerFn, 2.0)
    check h.lambdaParamProto of Gscalar
    check h.lambdaResultProto of GlambdaRef
    let z = apply(apply(h, 3.0), 4.0)
    z :~ 4.0

  test "cond rejects incompatible lambda branch prototypes":
    let scalarFn = lambdaParam(grt.localScalar(), grt.localScalar())
    let intFn = lambdaParam(grt.localInt(), grt.localInt())
    expect(GraphValueError):
      discard cond(grt.toGvalue(1), scalarFn, intFn)

  test "conditional lambda eval preserves branches until selected":
    grt.resetApplyCache()
    grt.resetGradCache()

    let k = grt.toGvalue(1.0)
    let selector = equal(k, 1.0)
    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let b = grt.toGvalue(5.0)
    let f = cond(selector, lambda(v, a * v), lambda(v, b * v))

    check f.isCondNode
    check selector.runCount == 0

    let z = Gscalar(apply(f, x))
    z :~ 6.0

    k.update 0.0
    z :~ 15.0

  test "conditional lambda applies through erased cond":
    let k = grt.toGvalue(1)
    let v = grt.localScalar()
    let f = lambda(v, v + 1.0)
    let g = lambda(v, v * 2.0)

    let erased = cond(k, f, g)
    apply(erased, 3.0) :~ 4.0

  test "lambda ref stays stale until reevaluated":
    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    let u = grt.localScalar()
    let hof = lambda(f, lambda(u, Gscalar(apply(f, u)) + 1.0))
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let g = lambda(v, a * v)
    let h = apply(hof, g)

    check h of GlambdaRef
    check h.lambdaParamProto of Gscalar
    check h.lambdaResultProto of Gscalar

    apply(h, 3.0) :~ 7.0

    discard h.eval
    let epoch0 = h.epoch
    let runs0 = h.runCount

    a.update 4.0
    check h.epoch == epoch0

    discard h.eval
    check h.runCount == runs0 + 1
    check h.epoch > epoch0
    apply(h, 3.0) :~ 13.0

  test "apply unresolved throws":
    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    let z = apply(f, 2.0)
    expect(GraphValueError):
      discard z.eval

  test "apply builds symbolically for unused incompatible direct argument":
    let v = grt.localScalar()
    let intArg = grt.toGvalue(1)
    let z = apply(lambda(v, grt.toGvalue(3.0)), intArg)

    z :~ 3.0

  test "apply builds symbolically for unused incompatible structural argument":
    let selector = grt.toGvalue(1)
    let i = grt.localInt()
    let f = cond(
      selector,
      lambda(i, grt.toGvalue(0.0)),
      lambda(i, grt.toGvalue(1.0)))
    let x = grt.toGvalue(2.0)

    apply(f, x) :~ 0.0

  test "apply gradient rejects incompatible structural argument result shape":
    grt.resetGradCache()

    let selector = grt.toGvalue(1)
    let i = grt.localInt()
    let f = cond(
      selector,
      lambda(i, grt.toGvalue(0.0)),
      lambda(i, grt.toGvalue(1.0)))
    let x = grt.toGvalue(2.0)
    let z = Gscalar(apply(f, x))

    expect(GraphValueError):
      discard grad(z, x)

  test "higher-order apply is symbolic but still rejects non-lambda calls":
    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    let u = grt.toGvalue(2.0)
    let hof = lambda(f, apply(f, u))
    let intParam = grt.localInt()
    let intFn = lambda(intParam, intParam)

    let intUse = apply(hof, intFn)
    let notLambda = apply(hof, grt.toGvalue(1.0))
    intUse :~ 2.0
    expect(GraphValueError):
      discard notLambda.eval

  test "constant higher-order lambda accepts unused incompatible argument":
    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    let hof = lambda(f, grt.toGvalue(5.0))
    let inner = lambdaParam(grt.localInt(), grt.localScalar())
    let badFn = lambda(inner, grt.toGvalue(0.0))

    apply(hof, badFn) :~ 5.0

  test "lambda ref erased binding rejects incompatible parameter prototype":
    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    let bad = lambdaParam(grt.localInt(), grt.localScalar())
    let epoch = f.epoch

    expect(GraphValueError):
      f.valCopy(Gvalue(bad))
    expect(GraphValueError):
      discard apply(f, 1.0).eval
    check f.epoch == epoch

  test "structural lambda int literal uses parameter prototype":
    let selector = grt.toGvalue(1)
    let v = grt.localScalar()
    let f = cond(
      selector,
      lambda(v, v + 1.0),
      lambda(v, v * 2.0))

    apply(f, 1) :~ 2.0

    selector.update 0
    apply(f, 1) :~ 2.0

  test "vjpOf unresolved lambda uses result cotangent prototype":
    block:
      let f = lambdaParam(grt.localScalar(), grt.localInt())
      let vf = vjpOf(f)
      let seedFunProto = vf.lambdaResultProto

      check vf.lambdaParamProto of Gscalar
      check seedFunProto != nil
      check seedFunProto.lambdaParamProto of Gint
      check seedFunProto.lambdaResultProto of Gscalar

    block:
      let f = lambdaParam(grt.localInt(), grt.localScalar())
      let vf = vjpOf(f)
      let seedFunProto = vf.lambdaResultProto

      check vf.lambdaParamProto of Gint
      check seedFunProto != nil
      check seedFunProto.lambdaParamProto of Gscalar
      check seedFunProto.lambdaResultProto of Gint

  test "vjpOf conditional lambda follows selected branch":
    let k = grt.toGvalue(1)
    let x = grt.toGvalue(3.0)
    let seed = grt.toGvalue(1.0)
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let b = grt.toGvalue(5.0)
    let f = cond(k, lambda(v, a * v), lambda(v, b * v * v))
    let vf = vjpOf(f)
    let dx = apply(apply(vf, x), seed)

    check vf.isCondNode
    dx :~ 2.0

    k.update 0
    dx :~ 30.0

  test "apply eval does not traverse unused lambda-valued argument":
    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    let hof = lambda(f, grt.toGvalue(1.0))
    let x = grt.toGvalue(3.0)
    let p = 1.0 / (x - 3.0)
    let v = grt.localScalar()
    let arg = lambda(v, p + v)
    let z = apply(hof, arg)

    z :~ 1.0
    check p.runCount == 0

  test "failed cross-runtime apply does not mutate lambda argument":
    let grt2 = initGraphRuntime()
    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    discard vjpOf(f)
    let u = grt.localScalar()
    let hof = lambda(f, lambda(u, apply(f, u)))

    let v = grt2.localScalar()
    let a = grt2.toGvalue(2.0)
    let arg = lambda(v, a * v)
    let epoch = arg.epoch

    expect(GraphValueError):
      discard apply(hof, arg)
    check arg.epoch == epoch

  test "ordinary actual does not mutate symbolic vjpOf source":
    grt.resetApplyCache()
    grt.resetGradCache()

    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    let epoch = f.epoch
    let vf = vjpOf(f)
    check f.epoch == epoch
    let u = grt.localScalar()
    let hof = lambda(f, Gscalar(apply(apply(vf, u), 1.0)))
    let actual = grt.toGvalue(2.0)
    let z = apply(hof, actual)

    expect(GraphValueError):
      discard z.eval
    check f.epoch == epoch

  test "apply evaluates nested lambdas symbolically":
    let innerProtoArg = grt.localScalar()
    let innerProtoRet = grt.localScalar()
    let scalarFnProto = lambda(innerProtoArg, innerProtoRet)
    let f = lambdaParam(grt.localScalar(), scalarFnProto)
    let hof = lambda(f, apply(apply(f, 1.0), 2.0))

    let makeArg = grt.localScalar()
    let makeInner = grt.localScalar()
    let maker = lambda(makeArg, lambda(makeInner, makeArg + makeInner))
    let z = apply(hof, maker)
    z :~ 3.0

  test "cyclic lambda binding raises explicitly":
    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    f.valCopy(f)

    expect(GraphValueError):
      discard apply(f, 1.0).eval

  test "Y combinator style recursion eval":
    let protoArg = grt.localScalar()
    let protoRet = grt.localScalar()
    let fnProto = lambda(protoArg, protoRet)
    let x = lambdaParam(fnProto, fnProto)
    let f = lambdaParam(fnProto, fnProto)
    let Y = lambda(f, apply(lambda(x, apply(f, apply(x, x))), lambda(x, apply(f, apply(x, x)))))

    let rf = lambdaParam(grt.localScalar(), grt.localScalar())
    let u = grt.localScalar()
    let v = grt.localScalar()
    let y = grt.toGvalue(4.0)
    let F = lambda(rf, lambda(u,
      cond(equal(u, 0.0), y,
        Gscalar(apply(lambda(v, Gscalar(apply(rf, v)) + Gscalar(apply(rf, v))), u - 1.0)))))

    let z = apply(apply(Y, F), 4.0)
    z :~ 64.0
