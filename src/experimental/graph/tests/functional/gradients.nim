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

  test "higher order apply gradient wrt closure and input":
    let f = local(grt.localScalar())
    let u = grt.localScalar()
    let hof = lambda(f, lambda(u, apply(f, u) * apply(f, u)))
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

  test "multiple applies accumulate gradients":
    let x = grt.toGvalue(2.0)
    let y = grt.toGvalue(5.0)
    let v = grt.localScalar()
    let f = lambda(v, v * y)
    let g = lambda(v, v + y)
    let z = apply(f, x) * apply(g, x + 1.0)
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
