import math, unittest

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))

import core, scalar, functional

template checkeq(ii: tuple[filename:string, line:int, column:int], sa: string, a: float, sb: string, b: float) =
  if not almostEqual(a, b, unitsInLastPlace = 64):
    checkpoint(ii.filename & ":" & $ii.line & ":" & $ii.column & ": Check failed: " & sa & " :~ " & sb)
    checkpoint("  " & sa & ": " & $a)
    checkpoint("  " & sb & ": " & $b)
    fail()

template `:~`(a: Gvalue, b: float) =
  checkeq(instantiationInfo(), astToStr a, a.eval.getfloat, astToStr b, b)

suite "functional lambda gradients":
  test "direct apply first and second derivative":
    let x = toGvalue(3.0)
    let v = localScalar()
    let z = apply(lambda(v, v * v + 2.0 * v), x)
    let dzdx = z.grad x
    let d2zdx2 = dzdx.grad x
    z :~ 15.0
    dzdx :~ 8.0
    d2zdx2 :~ 2.0

  test "chain rule through apply argument":
    let x = toGvalue(2.0)
    let v = localScalar()
    let z = apply(lambda(v, v + v * v), x * x)
    let dzdx = z.grad x
    let d2zdx2 = dzdx.grad x
    z :~ 20.0
    dzdx :~ 36.0
    d2zdx2 :~ 50.0

  test "higher order apply gradient wrt closure and input":
    let f = local(localScalar())
    let u = localScalar()
    let hof = lambda(f, lambda(u, apply(f, u) * apply(f, u)))
    let v = localScalar()
    let a = toGvalue(2.0)
    let x = toGvalue(3.0)
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
    let x = toGvalue(2.0)
    let y = toGvalue(5.0)
    let v = localScalar()
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
    let protoArg = localScalar()
    let protoRet = localScalar()
    let fnProto = lambda(protoArg, protoRet)
    let x = local(fnProto)
    let f = local(fnProto)
    let Y = lambda(f, apply(lambda(x, apply(f, apply(x, x))), lambda(x, apply(f, apply(x, x)))))

    let rf = local(localScalar())
    let u = localScalar()
    let v = localScalar()
    let base = toGvalue(1.0)
    let step = toGvalue(1.0)
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

suite "functional lambda":
  test "apply scalar and grad":
    let x = toGvalue(3.0)
    let v = localScalar()
    let z = apply(lambda(v, v + v), x * x)
    let dzdx = z.grad x
    z :~ 18.0
    dzdx :~ 12.0

  test "closure capture by reference":
    let x = toGvalue(2.0)
    let y = toGvalue(4.0)
    let v = localScalar()
    let z = apply(lambda(v, v * y + y), x)
    let dzdy = z.grad y
    z :~ 12.0
    dzdy :~ 3.0
    y.update 5.0
    z :~ 15.0
    dzdy :~ 3.0

  test "multiple applies in one graph":
    let x = toGvalue(2.0)
    let y = toGvalue(7.0)
    let v = localScalar()
    let f = lambda(v, v + y)
    let z = apply(f, x) * apply(f, x + 1.0)
    let dzdy = z.grad y
    z :~ 90.0
    dzdy :~ 19.0

  test "higher order function argument":
    let f = local(localScalar())
    let u = localScalar()
    let hof = lambda(f, lambda(u, apply(f, u) + 1.0))
    let v = localScalar()
    let a = toGvalue(2.0)
    let g = lambda(v, a * v)
    let z = apply(apply(hof, g), 3.0)
    let dzda = z.grad a
    z :~ 7.0
    dzda :~ 3.0
    a.update 4.0
    z :~ 13.0
    dzda :~ 3.0

  test "evaluated function-valued apply remains callable":
    let f = local(localScalar())
    let u = localScalar()
    let hof = lambda(f, lambda(u, apply(f, u) + 1.0))
    let v = localScalar()
    let a = toGvalue(2.0)
    let g = lambda(v, a * v)
    let h = apply(hof, g)
    check h of Gcallable
    discard h.eval

    let z = apply(h, 3.0)
    let dzda = z.grad a
    z :~ 7.0
    dzda :~ 3.0

  test "evaluated function-valued apply tracks closure updates":
    let f = local(localScalar())
    let u = localScalar()
    let hof = lambda(f, lambda(u, apply(f, u) + 1.0))
    let v = localScalar()
    let a = toGvalue(2.0)
    let g = lambda(v, a * v)
    let h = apply(hof, g)
    discard h.eval

    let z = apply(h, 3.0)
    z :~ 7.0

    a.update 4.0
    z :~ 13.0

  test "deferred apply unresolved throws":
    let f = local(localScalar())
    let z = apply(f, 2.0)
    expect(GraphValueError):
      discard z.eval

  test "Y combinator style recursion eval":
    let protoArg = localScalar()
    let protoRet = localScalar()
    let fnProto = lambda(protoArg, protoRet)
    let x = local(fnProto)
    let f = local(fnProto)
    let Y = lambda(f, apply(lambda(x, apply(f, apply(x, x))), lambda(x, apply(f, apply(x, x)))))

    let rf = local(localScalar())
    let u = localScalar()
    let v = localScalar()
    let y = toGvalue(4.0)
    let F = lambda(rf, lambda(u,
      cond(equal(u, 0.0), y,
        apply(lambda(v, apply(rf, v) + apply(rf, v)), u - 1.0))))

    let z = apply(apply(Y, F), 4.0)
    z :~ 64.0

suite "functional reuse":
  test "source lambda remains reusable after returned lambda reclosure":
    let outer = localScalar()
    let inner = localScalar()
    let a = toGvalue(10.0)
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
    resetApplyCacheStats()
    resetGradCacheStats()

    let f = local(localScalar())
    let x = toGvalue(3.0)
    let a = toGvalue(2.0)
    let v = localScalar()
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

    check cast[pointer](dz2) != p1
    check gradCacheStats.invalidations > 0
    check applyCacheStats.reduceMisses >= 2

  test "evaluated higher-order apply refreshes after callable rebind":
    resetApplyCacheStats()
    resetGradCacheStats()

    let f = local(localScalar())
    let u = localScalar()
    let hof = lambda(f, lambda(u, apply(f, u) + 1.0))
    let v = localScalar()
    let a = toGvalue(2.0)
    let b = toGvalue(3.0)
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

    check gradCacheStats.invalidations > 0
    check applyCacheStats.reduceMisses >= 2

  test "reuse forward values and gradient node on stable graph":
    resetApplyCacheStats()
    resetGradCacheStats()

    let x = toGvalue(3.0)
    let y = toGvalue(2.0)
    let mul = x * y
    let v = localScalar()
    let z = apply(lambda(v, v * v + v), mul)

    discard z.eval
    let applyRuns = z.gfunc.runCount

    let dzdx1 = z.grad x
    dzdx1 :~ 26.0

    discard dzdx1.eval

    let dzdx2 = z.grad x
    check cast[pointer](dzdx2) == cast[pointer](dzdx1)
    check gradCacheStats.directHits > 0
    check applyCacheStats.reduceHits > 0

    discard z.eval
    check z.gfunc.runCount == applyRuns

  test "grad sees apply closure deps before eval":
    let a = toGvalue(2.0)

    block:
      let f = local(localScalar())
      let dep = apply(f, 3.0)
      let v = localScalar()
      let g = lambda(v, a * v + 1.0)
      f.valCopy g

      let direct = dep.grad a
      dep :~ 7.0
      direct :~ 3.0

    a.update 4.0

    block:
      let f = local(localScalar())
      let dep = apply(f, 3.0)
      let v = localScalar()
      let g = lambda(v, a * v + 1.0)
      f.valCopy g

      let direct = dep.grad a
      dep :~ 13.0
      direct :~ 3.0

  test "lambda grad build does not eval apply":
    resetApplyCacheStats()

    let x = toGvalue(3.0)
    let v = localScalar()
    let e = exp(v)
    let z = apply(lambda(v, e + v), x)
    let applyRuns0 = z.gfunc.runCount
    let expRuns0 = e.gfunc.runCount

    let dzdx = z.grad x
    discard dzdx.grad x

    check z.gfunc.runCount == applyRuns0
    check e.gfunc.runCount == expRuns0
    check applyCacheStats.reduceMisses == 0
    check applyCacheStats.partialMisses == 0

  test "lambda reuses reduced work across forward and gradient eval":
    resetApplyCacheStats()
    resetGradCacheStats()

    block:
      let x = toGvalue(3.0)
      let v = localScalar()
      let e = exp(v)
      let z = apply(lambda(v, e + v), x)
      let dzdx = z.grad x
      let applyRuns0 = z.gfunc.runCount
      let expRuns0 = e.gfunc.runCount

      discard z.eval
      check z.gfunc.runCount == applyRuns0 + 1
      let expRuns1 = e.gfunc.runCount
      check expRuns1 > expRuns0

      discard dzdx.eval
      check z.gfunc.runCount == applyRuns0 + 1
      check e.gfunc.runCount == expRuns1

    block:
      let x = toGvalue(3.0)
      let v = localScalar()
      let e = exp(v)
      let z = apply(lambda(v, e + v), x)
      let dzdx = z.grad x
      let applyRuns0 = z.gfunc.runCount
      let expRuns0 = e.gfunc.runCount

      discard dzdx.eval
      check z.gfunc.runCount == applyRuns0
      let expRuns1 = e.gfunc.runCount
      check expRuns1 > expRuns0

      discard z.eval
      check z.gfunc.runCount == applyRuns0 + 1
      check e.gfunc.runCount == expRuns1

  test "apply partial depth guard raises without poisoning cache":
    resetApplyCacheStats()

    let x = toGvalue(3.0)
    let v = localScalar()
    let z = apply(lambda(v, v * v + 1.0), x)
    let dzdx = z.grad x
    let savedLimit = applyGradPrepareDepthLimit

    try:
      applyGradPrepareDepthLimit = 0
      expect(GraphValueError):
        discard dzdx.eval

      applyGradPrepareDepthLimit = savedLimit
      dzdx :~ 6.0
    finally:
      applyGradPrepareDepthLimit = savedLimit

  test "callable resolution depth overflow raises explicitly":
    let savedLimit = lambdaResolveDepthLimit

    try:
      lambdaResolveDepthLimit = 0
      let v = localScalar()
      expect(GraphValueError):
        discard apply(lambda(v, v + 1.0), 2.0)
    finally:
      lambdaResolveDepthLimit = savedLimit

  test "cond gradients stay live without cache invalidation":
    resetGradCacheStats()

    let x = toGvalue(2.0)
    let y = toGvalue(3.0)
    let k = toGvalue(0)
    let z = cond(k, x * x, y * y)

    let dzdx = z.grad x
    let p0 = cast[pointer](dzdx)
    dzdx :~ 0.0

    k.update 1
    dzdx :~ 4.0

    let dzdx1 = z.grad x
    check cast[pointer](dzdx1) == p0
    check gradCacheStats.directHits > 0
    check gradCacheStats.invalidations == 0

  test "cond keeps apply partials live across selector flips":
    resetApplyCacheStats()
    resetGradCacheStats()

    let x = toGvalue(3.0)
    let k = toGvalue(1)
    let v = localScalar()
    let z = cond(k, apply(lambda(v, v * v + 1.0), x), 0.0)

    let dz = z.grad x
    let p0 = cast[pointer](dz)
    check applyCacheStats.reduceMisses == 0
    check applyCacheStats.partialMisses == 0
    dz :~ 6.0
    let partialMisses = applyCacheStats.partialMisses
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
    check applyCacheStats.partialMisses == partialMisses
    check gradCacheStats.directHits > 0
    check gradCacheStats.invalidations == 0

  test "cond grad build does not reduce either apply branch":
    resetApplyCacheStats()

    let a = toGvalue(2.0)
    let x = toGvalue(3.0)
    let k = toGvalue(1)
    let v = localScalar()
    let t = apply(lambda(v, a * v + 1.0), x)
    let f = apply(lambda(v, a + v * v), x)
    let z = cond(k, t, f)

    discard z.grad a
    discard z.grad x

    check applyCacheStats.reduceMisses == 0
    check applyCacheStats.partialMisses == 0

  test "cond gradient eval reduces only the selected apply branch":
    resetApplyCacheStats()

    let a = toGvalue(2.0)
    let x = toGvalue(3.0)
    let k = toGvalue(1)
    let v = localScalar()
    let z = cond(k,
      apply(lambda(v, a * v + 1.0), x),
      apply(lambda(v, a + v * v), x))
    let dzda = z.grad a

    check applyCacheStats.reduceMisses == 0
    dzda :~ 3.0
    check applyCacheStats.reduceMisses == 1

    k.update 0
    dzda :~ 1.0
    check applyCacheStats.reduceMisses == 2

  test "apply partials refresh on deep capture updates before eval":
    resetApplyCacheStats()
    resetGradCacheStats()

    let protoArg = localScalar()
    let protoRet = localScalar()
    let fnProto = lambda(protoArg, protoRet)
    let x = local(fnProto)
    let f = local(fnProto)
    let Y = lambda(f, apply(lambda(x, apply(f, apply(x, x))), lambda(x, apply(f, apply(x, x)))))

    let rf = local(localScalar())
    let u = localScalar()
    let v = localScalar()
    let base = toGvalue(1.0)
    let step = toGvalue(1.0)
    let F = lambda(rf, lambda(u,
      cond(equal(u, 0.0), base,
        apply(lambda(v, apply(rf, v) + step), u - 1.0))))

    let z = apply(apply(Y, F), 3.0)
    let dzdstep = z.grad step
    check applyCacheStats.reduceMisses == 0
    dzdstep :~ 3.0
    let reduceMisses = applyCacheStats.reduceMisses
    check reduceMisses > 0

    base.update 2.0
    let dzdbase = z.grad base
    check applyCacheStats.reduceMisses == reduceMisses
    dzdbase :~ 1.0
    check applyCacheStats.reduceMisses > reduceMisses
    z :~ 5.0
