suite "functional lambda":
  test "apply scalar and grad":
    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let z = apply(lambda(v, v + v), x * x)
    let dzdx = z.grad x
    z :~ 18.0
    dzdx :~ 12.0

  test "closure capture by reference":
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
    let f = local(grt.localScalar())
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

  test "evaluated function-valued apply remains callable":
    let f = local(grt.localScalar())
    let u = grt.localScalar()
    let hof = lambda(f, lambda(u, Gscalar(apply(f, u)) + 1.0))
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let g = lambda(v, a * v)
    let h = apply(hof, g)
    check h of Gwrapper
    if h of Gwrapper:
      check Gwrapper(h).kind == wkCallable
    discard h.eval

    let z = apply(h, 3.0)
    let dzda = z.grad a
    z :~ 7.0
    dzda :~ 3.0

  test "evaluated function-valued apply tracks closure updates":
    let f = local(grt.localScalar())
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

  test "evaluated callable refreshes only when captured values change":
    grt.resetApplyCache()

    let f = local(grt.localScalar())
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

  test "local wrapper exposes current binding through traversal modes":
    let f = local(grt.localScalar())
    let a = grt.toGvalue(2.0)
    let v = grt.localScalar()
    let epoch0 = f.epoch
    f.valCopy lambda(v, a * v + 1.0)

    check f of Gwrapper
    if f of Gwrapper:
      check Gwrapper(f).kind == wkLocal
    check f.epoch > epoch0

    let rawTree = f.treeRepr
    let depTree = f.treeRepr(iwmDepend)
    let sigTree = f.treeRepr(iwmGradSignature)

    check depTree.splitLines.len > rawTree.splitLines.len
    check sigTree.splitLines.len > rawTree.splitLines.len

    let epoch1 = f.epoch
    f.valCopy lambda(v, a * v + 2.0)
    check f.epoch > epoch1

  test "callable placeholder rejects non-callable binding":
    let f = local(grt.localScalar())
    check not compiles(f.valCopy 3.0)
    check not compiles(f.valCopy 3)

    let v = grt.localScalar()
    let valid = lambda(v, v + 1.0)
    f.valCopy valid
    let bound = Gwrapper(f).bound
    let epoch = f.epoch

    expect(GraphValueError):
      f.valCopy(grt.toGvalue(3.0))

    check sameNode(Gwrapper(f).bound, bound)
    check f.epoch == epoch

  test "callable wrapper node requires explicit result prototype":
    check not compiles(callableWrapperNode())
    let missing: Gvalue = nil
    let rawPrototype = Gvalue()
    expect(GraphValueError):
      discard callableWrapperNode(missing)
    expect(GraphValueError):
      discard callableWrapperNode(rawPrototype)

  test "local callable placeholder requires runtime-bearing prototype":
    check not compiles(local(grt))
    let rawPrototype = Gvalue()
    expect(GraphValueError):
      discard local(rawPrototype)

  test "callable placeholders require result prototypes before binding":
    let malformed = Gwrapper(kind: wkLocal).attachRuntime(grt)
    let v = grt.localScalar()

    expect(GraphValueError):
      malformed.valCopy lambda(v, v + 1.0)

  test "apply result prototype validates graph values":
    let scalar = grt.localScalar()
    let intValue = grt.localInt()
    let param = grt.localScalar()
    let fn = lambda(param, param + 1.0)
    let wrapper = callableWrapperNode(grt.localScalar())

    check applyResultProto(fn) != nil
    check applyResultProto(wrapper) != nil
    check applyResultProto(scalar) == nil
    check applyResultProto(intValue) == nil

    let missing: Gvalue = nil
    let raw = Gvalue()
    expect(GraphValueError):
      discard applyResultProto(missing)
    expect(GraphValueError):
      discard applyResultProto(raw)

  test "apply result prototype materializes nested callable wrappers":
    let finalProto = grt.localScalar()
    let innerParam = grt.localScalar()
    let innerFnProto = lambda(innerParam, finalProto)
    let outerParam = grt.localScalar()
    let outerFn = lambda(outerParam, innerFnProto)

    let proto = applyResultProto(outerFn)

    check proto of Gwrapper
    if proto of Gwrapper:
      let wrapper = Gwrapper(proto)
      check wrapper.kind == wkCallable
      check wrapper.bound == nil
      check wrapper.retProto of Gscalar
      if wrapper.retProto != nil:
        check wrapper.retProto.copyCompatible(finalProto)

  test "apply result prototype depth limit rejects too-deep callable wrappers":
    let savedLimit = grt.lambdaResolveDepthLimit
    grt.lambdaResolveDepthLimit = 1
    try:
      let finalProto = grt.localScalar()
      let innerParam = grt.localScalar()
      let innerFnProto = lambda(innerParam, finalProto)
      let middleParam = grt.localScalar()
      let middleFnProto = lambda(middleParam, innerFnProto)
      let outerParam = grt.localScalar()
      let outerFn = lambda(outerParam, middleFnProto)

      expect(GraphValueError):
        discard applyResultProto(outerFn)
    finally:
      grt.lambdaResolveDepthLimit = savedLimit

  test "apply result prototype accepts nesting at the depth limit":
    let savedLimit = grt.lambdaResolveDepthLimit
    grt.lambdaResolveDepthLimit = 2
    try:
      let finalProto = grt.localScalar()
      let innerParam = grt.localScalar()
      let innerFnProto = lambda(innerParam, finalProto)
      let middleParam = grt.localScalar()
      let middleFnProto = lambda(middleParam, innerFnProto)
      let outerParam = grt.localScalar()
      let outerFn = lambda(outerParam, middleFnProto)

      let proto = applyResultProto(outerFn)

      check proto of Gwrapper
      if proto of Gwrapper:
        let outerWrapper = Gwrapper(proto)
        check outerWrapper.retProto of Gwrapper
        if outerWrapper.retProto of Gwrapper:
          let innerWrapper = Gwrapper(outerWrapper.retProto)
          check innerWrapper.retProto.copyCompatible(finalProto)
    finally:
      grt.lambdaResolveDepthLimit = savedLimit

  test "cond rejects incompatible callable branch prototypes":
    let scalarFn = local(grt.localScalar())
    let intFn = local(grt.localInt())
    expect(GraphValueError):
      discard cond(grt.toGvalue(1), scalarFn, intFn)

  test "callable wrapper stays stale until reevaluated":
    let f = local(grt.localScalar())
    let u = grt.localScalar()
    let hof = lambda(f, lambda(u, Gscalar(apply(f, u)) + 1.0))
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let g = lambda(v, a * v)
    let h = apply(hof, g)

    check h of Gwrapper
    if h of Gwrapper:
      check Gwrapper(h).kind == wkCallable

    let rawTree = h.treeRepr
    let depTree = h.treeRepr(iwmDepend)
    let sigTree = h.treeRepr(iwmGradSignature)
    check depTree.splitLines.len > rawTree.splitLines.len
    check sigTree.splitLines.len > rawTree.splitLines.len

    discard h.eval
    let epoch0 = h.epoch
    let runs0 = h.runCount

    a.update 4.0
    check h.epoch == epoch0

    discard h.eval
    check h.runCount == runs0 + 1
    check h.epoch > epoch0

  test "deferred apply unresolved throws":
    let f = local(grt.localScalar())
    let z = apply(f, 2.0)
    expect(GraphValueError):
      discard z.eval

  test "apply rejects direct lambda argument shape mismatch":
    let v = grt.localScalar()
    let intArg = grt.toGvalue(1)

    expect(GraphValueError):
      discard apply(lambda(v, v + 1.0), intArg)

  test "deferred apply rejects resolved lambda argument shape mismatch":
    let f = local(grt.localScalar())
    let x = grt.toGvalue(2.0)
    let z = apply(f, x)
    let i = grt.localInt()
    f.valCopy lambda(i, grt.toGvalue(0.0))

    expect(GraphValueError):
      discard z.eval

    let v = grt.localScalar()
    f.valCopy lambda(v, v + 1.0)
    z :~ 3.0

  test "deferred apply gradient rejects resolved lambda argument shape mismatch":
    grt.resetGradCache()

    let f = local(grt.localScalar())
    let x = grt.toGvalue(2.0)
    let z = apply(f, x)
    let dzdx = z.grad x
    let i = grt.localInt()
    f.valCopy lambda(i, grt.toGvalue(0.0))

    expect(GraphValueError):
      discard dzdx.eval

    let v = grt.localScalar()
    f.valCopy lambda(v, v * v)
    dzdx :~ 4.0

  test "apply rejects incompatible callable wrapper argument":
    let f = local(grt.localScalar())
    let u = grt.localScalar()
    let hof = lambda(f, apply(f, u))
    let intParam = grt.localInt()
    let intFn = lambda(intParam, intParam)

    expect(GraphValueError):
      discard apply(hof, intFn)
    expect(GraphValueError):
      discard apply(hof, grt.toGvalue(1.0))

  test "apply compares nested callable result prototypes":
    let innerProtoArg = grt.localScalar()
    let innerProtoRet = grt.localScalar()
    let scalarFnProto = lambda(innerProtoArg, innerProtoRet)
    let f = local(scalarFnProto)
    let hof = lambda(f, apply(apply(f, 1.0), 2.0))

    let makeArg = grt.localScalar()
    let makeInner = grt.localScalar()
    let maker = lambda(makeArg, lambda(makeInner, makeArg + makeInner))
    let z = apply(hof, maker)
    z :~ 3.0

    let badArg = grt.localScalar()
    let intParam = grt.localInt()
    let badMaker = lambda(badArg, lambda(intParam, intParam))
    expect(GraphValueError):
      discard apply(hof, badMaker)

  test "callable resolution rejects nil reductions":
    proc nilReduce(v: Gvalue): Gvalue =
      discard v
      nil

    let x = grt.localScalar()
    let badReduce = newGfunc(reduceCallable = nilReduce, name = "badReduce")
    let badNode = graphNode(scalarNodeLike(x), @[Gvalue(x)], badReduce, "badReduce")
    expect(GraphValueError):
      discard badNode.resolveCallableChain(crmReduced)

  test "Y combinator style recursion eval":
    let protoArg = grt.localScalar()
    let protoRet = grt.localScalar()
    let fnProto = lambda(protoArg, protoRet)
    let x = local(fnProto)
    let f = local(fnProto)
    let Y = lambda(f, apply(lambda(x, apply(f, apply(x, x))), lambda(x, apply(f, apply(x, x)))))

    let rf = local(grt.localScalar())
    let u = grt.localScalar()
    let v = grt.localScalar()
    let y = grt.toGvalue(4.0)
    let F = lambda(rf, lambda(u,
      cond(equal(u, 0.0), y,
        Gscalar(apply(lambda(v, Gscalar(apply(rf, v)) + Gscalar(apply(rf, v))), u - 1.0)))))

    let z = apply(apply(Y, F), 4.0)
    z :~ 64.0
