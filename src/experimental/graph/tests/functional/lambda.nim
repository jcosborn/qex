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
    let z = apply(f, x) * apply(f, x + 1.0)
    let dzdy = z.grad y
    z :~ 90.0
    dzdy :~ 19.0

  test "higher order function argument":
    let f = local(grt.localScalar())
    let u = grt.localScalar()
    let hof = lambda(f, lambda(u, apply(f, u) + 1.0))
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
    let hof = lambda(f, lambda(u, apply(f, u) + 1.0))
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let g = lambda(v, a * v)
    let h = apply(hof, g)
    check h.isCallableWrapper
    discard h.eval

    let z = apply(h, 3.0)
    let dzda = z.grad a
    z :~ 7.0
    dzda :~ 3.0

  test "evaluated function-valued apply tracks closure updates":
    let f = local(grt.localScalar())
    let u = grt.localScalar()
    let hof = lambda(f, lambda(u, apply(f, u) + 1.0))
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
    let hof = lambda(f, lambda(u, apply(f, u) + 1.0))
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
    let epoch0 = f.epochOf
    f.valCopy lambda(v, a * v + 1.0)

    check f.isLocalWrapper
    check not f.isCallableWrapper
    check f.epochOf > epoch0

    let rawTree = f.treeRepr
    let depTree = f.treeRepr(iwmDepend)
    let sigTree = f.treeRepr(iwmGradSignature)

    check depTree.splitLines.len > rawTree.splitLines.len
    check sigTree.splitLines.len > rawTree.splitLines.len

    let epoch1 = f.epochOf
    f.valCopy lambda(v, a * v + 2.0)
    check f.epochOf > epoch1

  test "callable placeholder rejects non-callable binding":
    let f = local(grt.localScalar())
    expect(GraphValueError):
      f.valCopy 3.0

  test "cond rejects incompatible callable branch prototypes":
    let scalarFn = local(grt.localScalar())
    let intFn = local(grt.localInt())
    expect(GraphValueError):
      discard cond(grt.toGvalue(1), scalarFn, intFn)

  test "callable wrapper stays stale until reevaluated":
    let f = local(grt.localScalar())
    let u = grt.localScalar()
    let hof = lambda(f, lambda(u, apply(f, u) + 1.0))
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let g = lambda(v, a * v)
    let h = apply(hof, g)

    check h.isCallableWrapper
    check not h.isLocalWrapper

    let rawTree = h.treeRepr
    let depTree = h.treeRepr(iwmDepend)
    let sigTree = h.treeRepr(iwmGradSignature)
    check depTree.splitLines.len > rawTree.splitLines.len
    check sigTree.splitLines.len > rawTree.splitLines.len

    discard h.eval
    let epoch0 = h.epochOf
    let runs0 = h.runCount

    a.update 4.0
    check h.epochOf == epoch0

    discard h.eval
    check h.runCount == runs0 + 1
    check h.epochOf > epoch0

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

  test "apply rejects malformed node arity early":
    let v = grt.localScalar()
    let z = apply(lambda(v, v + 1.0), 2.0)
    z.inputs.add grt.toGvalue(3.0)
    expect(GraphValueError):
      discard z.eval

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
        apply(lambda(v, apply(rf, v) + apply(rf, v)), u - 1.0))))

    let z = apply(apply(Y, F), 4.0)
    z :~ 64.0
