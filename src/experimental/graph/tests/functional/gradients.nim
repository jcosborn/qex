suite "functional lambda gradients":
  test "direct apply first and second derivative":
    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let z = apply(lambda(v, v * v + 2.0 * v), x)
    let dzdx = z.grad x
    let d2zdx2 = dzdx.grad x
    z :~ 15.0
    dzdx :~ 8.0
    d2zdx2 :~ 2.0

  test "chain rule through apply argument":
    let x = grt.toGvalue(2.0)
    let v = grt.localScalar()
    let z = apply(lambda(v, v + v * v), x * x)
    let dzdx = z.grad x
    let d2zdx2 = dzdx.grad x
    z :~ 20.0
    dzdx :~ 36.0
    d2zdx2 :~ 50.0

  test "apply argument also captured by lambda contributes both gradients":
    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let f = lambda(v, x * v)
    let z = apply(f, x)
    let dzdx = z.grad x
    let d2zdx2 = dzdx.grad x

    z :~ 9.0
    dzdx :~ 6.0
    d2zdx2 :~ 2.0

    x.update 4.0
    z :~ 16.0
    dzdx :~ 8.0
    d2zdx2 :~ 2.0

  test "apply argument expression and capture both contribute to same leaf":
    let a = grt.toGvalue(2.0)
    let x = grt.toGvalue(3.0)
    let v = grt.localScalar()
    let arg = a * x
    let f = lambda(v, a * v)
    let z = Gscalar(apply(f, arg))
    let dzda = z.grad a
    let dzdx = z.grad x

    z :~ 12.0
    dzda :~ 12.0
    dzdx :~ 4.0

    a.update 4.0
    z :~ 48.0
    dzda :~ 24.0
    dzdx :~ 16.0

  test "apply backward uses traversal index and verifies dependency identity":
    let a = grt.toGvalue(2.0)
    let x = grt.toGvalue(3.0)
    let unrelated = grt.toGvalue(5.0)
    let v = grt.localScalar()
    let z = Gscalar(apply(lambda(v, a * v + 1.0), x))

    let dx = z.gfunc.backward(nil, z, 0, x)
    let da = z.gfunc.backward(nil, z, 1, a)
    dx :~ 2.0
    da :~ 3.0

    expect(GraphValueError):
      discard z.gfunc.backward(nil, z, 0, unrelated)
    expect(GraphValueError):
      discard z.gfunc.backward(nil, z, 0, a)

  test "higher order apply gradient wrt lambda capture and input":
    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    let u = grt.localScalar()
    let hof = lambda(f, lambda(u, Gscalar(apply(f, u)) * Gscalar(apply(f, u))))
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let x = grt.toGvalue(3.0)
    let g = lambda(v, a * v + 1.0)
    let z = apply(apply(hof, g), x)
    let dzda = z.grad a
    let d2zda2 = dzda.grad a
    let dzdx = z.grad x
    z :~ 49.0
    dzda :~ 42.0
    d2zda2 :~ 18.0
    dzdx :~ 28.0
    a.update 4.0
    z :~ 169.0
    dzda :~ 78.0
    d2zda2 :~ 18.0
    dzdx :~ 104.0

  test "higher order direct apply gradient wrt lambda argument capture":
    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    let a = grt.toGvalue(3.0)
    let hof = lambda(f, apply(f, a))
    let v = grt.localScalar()
    let c = grt.toGvalue(2.0)
    let g = lambda(v, c * v * v)

    let z = apply(hof, g)
    let dzdc = z.grad c
    let dzda = z.grad a

    z :~ 18.0
    dzdc :~ 9.0
    dzda :~ 12.0

    a.update 4.0
    z :~ 32.0
    dzdc :~ 16.0
    dzda :~ 16.0

    c.update 5.0
    z :~ 80.0
    dzdc :~ 16.0
    dzda :~ 40.0

  test "function-position apply VJP propagates direct and argument targets":
    grt.resetApplyCache()
    grt.resetGradCache()

    let a = grt.toGvalue(2.0)
    let b = grt.toGvalue(3.0)
    let x = grt.toGvalue(4.0)
    let p = grt.localScalar()
    let v = grt.localScalar()
    let hof = lambda(p, lambda(v, a * v + p * v))
    let fun = apply(hof, b)
    let z = Gscalar(apply(fun, x))
    let dzda = z.grad a
    let dzdb = z.grad b
    let dzdx = z.grad x

    z :~ 20.0
    dzda :~ 4.0
    dzdb :~ 4.0
    dzdx :~ 5.0

    b.update 5.0
    z :~ 28.0
    dzda :~ 4.0
    dzdb :~ 4.0
    dzdx :~ 7.0

  test "function-position apply VJP includes applied argument value target":
    grt.resetApplyCache()
    grt.resetGradCache()

    let a = grt.toGvalue(2.0)
    let b = grt.toGvalue(3.0)
    let x = grt.toGvalue(4.0)
    let p = grt.localScalar()
    let v = grt.localScalar()
    let hof = lambda(p, lambda(v, a * v + p * v + p * p))
    let fun = apply(hof, b)
    let z = Gscalar(apply(fun, x))
    let dzda = z.grad a
    let dzdb = z.grad b
    let dzdx = z.grad x

    z :~ 29.0
    dzda :~ 4.0
    dzdb :~ 10.0
    dzdx :~ 5.0

    b.update 5.0
    z :~ 53.0
    dzda :~ 4.0
    dzdb :~ 14.0
    dzdx :~ 7.0

  test "function-position apply VJP includes applied argument expression target":
    grt.resetApplyCache()
    grt.resetGradCache()

    let a = grt.toGvalue(2.0)
    let b = grt.toGvalue(3.0)
    let x = grt.toGvalue(4.0)
    let p = grt.localScalar()
    let v = grt.localScalar()
    let hof = lambda(p, lambda(v, a * v + p * v + p * p))
    let fun = apply(hof, b * x)
    let z = Gscalar(apply(fun, x))
    let dzda = z.grad a
    let dzdb = z.grad b
    let dzdx = z.grad x

    z :~ 200.0
    dzda :~ 4.0
    dzdb :~ 112.0
    dzdx :~ 98.0

    b.update 5.0
    z :~ 488.0
    dzda :~ 4.0
    dzdb :~ 176.0
    dzdx :~ 242.0

  test "structural VJP build does not mutate captured apply inputs":
    grt.resetApplyCache()
    grt.resetGradCache()

    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let x = grt.toGvalue(3.0)
    let capturedApply = apply(lambda(v, a * v), x)

    let u = grt.localScalar()
    let scale = grt.toGvalue(5.0)
    let outer = lambda(u, Gscalar(capturedApply) * u)
    let z = apply(outer, scale)
    let dzda = z.grad a

    z :~ 30.0
    dzda :~ 15.0

  test "higher order direct apply sums function and lambda argument capture bars":
    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    let a = grt.toGvalue(3.0)
    let c = grt.toGvalue(2.0)
    let hof = lambda(f, c * Gscalar(apply(f, a)))
    let v = grt.localScalar()
    let g = lambda(v, c * v)

    let z = apply(hof, g)
    let dzdc = z.grad c

    z :~ 12.0
    dzdc :~ 12.0

    c.update 4.0
    z :~ 48.0
    dzdc :~ 24.0

  test "conditional structural lambda VJP contributes through capture":
    grt.resetApplyCache()
    grt.resetGradCache()

    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    let x = grt.toGvalue(2.0)
    let y = grt.toGvalue(5.0)
    let kSource = grt.toGvalue(1.0)
    let k = equal(kSource, 1.0)
    let hof = lambda(f,
      Gscalar(apply(f, x)) + cond(k, Gscalar(apply(f, y)), 0.0))

    let v = grt.localScalar()
    let c = grt.toGvalue(3.0)
    let z = apply(hof, lambda(v, c * v))

    let kRuns0 = k.runCount
    let dzdc = z.grad c
    check k.runCount == kRuns0

    dzdc :~ 7.0
    kSource.update 0.0
    dzdc :~ 2.0

  test "conditional lambda actual VJP contributes through selected capture":
    grt.resetApplyCache()
    grt.resetGradCache()

    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    let x = grt.toGvalue(3.0)
    let hof = lambda(f, Gscalar(apply(f, x)))
    let k = grt.toGvalue(1)
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let b = grt.toGvalue(5.0)
    let arg = cond(k, lambda(v, a * v), lambda(v, b * v * v))
    let z = apply(hof, arg)
    let dzda = z.grad a
    let dzdb = z.grad b

    z :~ 6.0
    dzda :~ 3.0
    dzdb :~ 0.0

    k.update 0
    z :~ 45.0
    dzda :~ 0.0
    dzdb :~ 9.0

  test "symbolic vjpOf distributes over conditional lambda substitution":
    grt.resetApplyCache()
    grt.resetGradCache()

    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    let x = grt.toGvalue(3.0)
    let seed = grt.toGvalue(1.0)
    let vf = vjpOf(f)
    let hof = lambda(f, apply(apply(vf, x), seed))
    let k = grt.toGvalue(1)
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let b = grt.toGvalue(5.0)
    let arg = cond(k, lambda(v, a * v), lambda(v, b * v * v))
    let dx = apply(hof, arg)

    dx :~ 2.0

    k.update 0
    dx :~ 30.0

  test "apply VJP value target is substituted when cloned into lambda":
    let p = grt.localScalar()
    let v = grt.localScalar()
    let body = apply(lambda(v, p * v), p).grad p
    let d = lambda(p, body)

    apply(d, 3.0) :~ 6.0

  test "nested returned lambda remaps structural VJP placeholders":
    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    let vf = vjpOf(f)
    let u = grt.localScalar()
    let w = grt.localScalar()
    let hof = lambda(f, lambda(u, lambda(w,
      apply(apply(vf, u * w), 1.0))))
    let v = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let h = apply(hof, lambda(v, a * v * v))
    let z = apply(apply(h, 3.0), 4.0)

    z :~ 48.0
    a.update 3.0
    z :~ 72.0

  test "nested lambda clone preserves shadowing":
    let outer = grt.localScalar()
    let inner = grt.localScalar()
    let a = grt.toGvalue(2.0)
    let f = lambda(outer, lambda(inner, outer * inner + inner))
    let h = apply(f, a)
    let z = apply(h, 5.0)
    let dzda = z.grad a

    z :~ 15.0
    dzda :~ 5.0

    a.update 3.0
    z :~ 20.0
    dzda :~ 5.0

  test "structural lambda VJP higher-order gradient uses propagated input":
    let f = lambdaParam(grt.localScalar(), grt.localScalar())
    let a = grt.toGvalue(5.0)
    let c = grt.toGvalue(2.0)
    let hof = lambda(f, c * Gscalar(apply(f, a)))
    let v = grt.localScalar()
    let z = apply(hof, lambda(v, c * v))

    let dzdc = z.grad c
    let d2zdc2 = dzdc.grad c

    dzdc :~ 20.0
    d2zdc2 :~ 10.0

  test "multiple applies accumulate gradients":
    let x = grt.toGvalue(2.0)
    let y = grt.toGvalue(5.0)
    let v = grt.localScalar()
    let f = lambda(v, v * y)
    let g = lambda(v, v + y)
    let z = Gscalar(apply(f, x)) * Gscalar(apply(g, x + 1.0))
    let dzdy = z.grad y
    let dzdx = z.grad x
    z :~ 80.0
    dzdy :~ 26.0
    dzdx :~ 50.0
    y.update 4.0
    z :~ 56.0
    dzdy :~ 22.0
    dzdx :~ 36.0

  test "Y combinator recursion grad wrt step and base":
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
    let dzdbase = z.grad base

    z :~ 4.0
    dzdstep :~ 3.0
    dzdbase :~ 1.0

    base.update 2.0
    z :~ 5.0
    dzdstep :~ 3.0
    dzdbase :~ 1.0

    let z0 = apply(apply(Y, F), 0.0)
    let dz0dstep = z0.grad step
    let dz0dbase = z0.grad base
    z0 :~ 2.0
    dz0dstep :~ 0.0
    dz0dbase :~ 1.0

  test "same-shaped recursive closures keep nominal VJP state separate":
    let protoArg = grt.localScalar()
    let protoRet = grt.localScalar()
    let fnProto = lambda(protoArg, protoRet)
    let x1 = lambdaParam(fnProto, fnProto)
    let f1 = lambdaParam(fnProto, fnProto)
    let y1 = lambda(f1, apply(lambda(x1, apply(f1, apply(x1, x1))), lambda(x1, apply(f1, apply(x1, x1)))))
    let x2 = lambdaParam(fnProto, fnProto)
    let f2 = lambdaParam(fnProto, fnProto)
    let y2 = lambda(f2, apply(lambda(x2, apply(f2, apply(x2, x2))), lambda(x2, apply(f2, apply(x2, x2)))))

    let rf1 = lambdaParam(grt.localScalar(), grt.localScalar())
    let rf2 = lambdaParam(grt.localScalar(), grt.localScalar())
    let u1 = grt.localScalar()
    let u2 = grt.localScalar()
    let v1 = grt.localScalar()
    let v2 = grt.localScalar()
    let base = grt.toGvalue(1.0)
    let step = grt.toGvalue(1.0)
    let F1 = lambda(rf1, lambda(u1,
      cond(equal(u1, 0.0), base,
        Gscalar(apply(lambda(v1, Gscalar(apply(rf1, v1)) + step), u1 - 1.0)))))
    let F2 = lambda(rf2, lambda(u2,
      cond(equal(u2, 0.0), base,
        Gscalar(apply(lambda(v2, Gscalar(apply(rf2, v2)) + step * 2.0), u2 - 1.0)))))

    let z = Gscalar(apply(apply(y1, F1), 2.0)) +
      Gscalar(apply(apply(y2, F2), 3.0))
    let dzdstep = z.grad step
    let dzdbase = z.grad base

    z :~ 10.0
    dzdstep :~ 8.0
    dzdbase :~ 2.0
