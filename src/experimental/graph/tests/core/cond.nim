type
  GcondWide = ref object of Gvalue
  GcondNarrow = ref object of GcondWide

method copyCompatible*(prototype: GcondWide, value: Gvalue): bool =
  value of GcondWide

method copyCompatible*(prototype: GcondNarrow, value: Gvalue): bool =
  value of GcondNarrow

suite "bool and cond":
  setup:
    let values = sampleScalarSweep()
    let a {.used.} = values.a
    let b {.used.} = values.b
    let c {.used.} = values.c
    let d {.used.} = values.d
    let x {.used.} = grt.toGvalue(a)
    let y {.used.} = grt.toGvalue(b)

  test "not":
    let f = grt.toGvalue(0)
    not(f) :~ 1
    not(not f) :~ 0
    let t = grt.toGvalue(1.0)
    not(t) :~ 0.0
    not(not t) :~ 1.0

  test "and":
    let fi = grt.toGvalue(0)
    let ti = grt.toGvalue(1)
    let t = grt.toGvalue(1.0)
    let f = grt.toGvalue(0.0)
    fi and t :~ 0.0
    t and fi :~ 0
    ti and t :~ 1.0
    t and ti :~ 1
    f and fi :~ 0
    fi and f :~ 0.0

  test "or":
    let fi = grt.toGvalue(0)
    let ti = grt.toGvalue(1)
    let t = grt.toGvalue(1.0)
    let f = grt.toGvalue(0.0)
    fi or t :~ 1.0
    t or fi :~ 1
    ti or t :~ 1.0
    t or ti :~ 1
    f or fi :~ 0
    fi or f :~ 0.0

  test "xor":
    let fi = grt.toGvalue(0)
    let ti = grt.toGvalue(1)
    let t = grt.toGvalue(1.0)
    let f = grt.toGvalue(0.0)
    fi xor t :~ 1.0
    t xor fi :~ 1
    ti xor t :~ 0.0
    t xor ti :~ 0
    f xor fi :~ 0
    fi xor f :~ 0.0

  test "strict scalar greater-than stays strict on equality":
    let sx = grt.toGvalue(2.0)
    let sy = grt.toGvalue(2.0)
    (sx > sy) :~ 0.0
    (sx >= sy) :~ 1.0
    (sx <= sy) :~ 1.0
    sy.update 3.0
    (sy > sx) :~ 1.0
    (sx > sy) :~ 0.0

  test "strict int greater-than stays strict on equality":
    let ix = grt.toGvalue(2)
    let iy = grt.toGvalue(2)
    (ix > iy) :~ 0
    (ix >= iy) :~ 1
    (ix <= iy) :~ 1
    iy.update 3
    (iy > ix) :~ 1
    (ix > iy) :~ 0

  test "comparison literal overloads preserve concrete result type":
    let sx = grt.toGvalue(2.0)
    let sy = grt.toGvalue(3.0)
    let ix = grt.toGvalue(2)
    let iy = grt.toGvalue(3)

    let slt: Gscalar = sx < sy
    let seqv: Gscalar = equal(sx, 2.0)
    let seqLeft: Gscalar = equal(2.0, sx)
    let sltRight: Gscalar = sx < 3
    let sltLeft: Gscalar = 1 < sx
    let sgeRight: Gscalar = sx >= 2.0
    let sleLeft: Gscalar = 2.0 <= sx

    slt :~ 1.0
    seqv :~ 1.0
    seqLeft :~ 1.0
    sltRight :~ 1.0
    sltLeft :~ 1.0
    sgeRight :~ 1.0
    sleLeft :~ 1.0
    equal(sx, sy) :~ 0.0

    let ilt: Gint = ix < iy
    let ieq: Gint = equal(ix, 2)
    let ieqLeft: Gint = equal(2, ix)
    let iltRight: Gint = ix < 3
    let iltLeft: Gint = 1 < ix
    let igeRight: Gint = ix >= 2
    let ileLeft: Gint = 2 <= ix

    ilt :~ 1
    ieq :~ 1
    ieqLeft :~ 1
    iltRight :~ 1
    iltLeft :~ 1
    igeRight :~ 1
    ileLeft :~ 1
    equal(ix, iy) :~ 0

  test "condi":
    let k = grt.toGvalue(0)
    let z = cond(k, x, y)
    let dx = z.grad x
    let dy = z.grad y
    z :~ b
    dx :~ 0.0
    dy :~ 1.0
    k.update 1
    z :~ a
    dx :~ 1.0
    dy :~ 0.0

  test "condParts exposes selector and branches":
    let k = grt.toGvalue(1)
    let trueBranch = x + 1.0
    let falseBranch = y + 1.0
    let z = cond(k, trueBranch, falseBranch)
    let parts = z.condParts

    check parts.selector.nodeKey == k.nodeKey
    check parts.whenTrue.nodeKey == trueBranch.nodeKey
    check parts.whenFalse.nodeKey == falseBranch.nodeKey

    expect(GraphValueError):
      discard x.condParts

  test "graphNode clears stale static zero markers":
    proc identityForward(v: Gvalue) =
      v.valCopy(v.inputs[0])

    proc identityBackward(zb: Gvalue,
                          z: Gvalue,
                          i: int,
                          input: Gvalue): Gvalue =
      discard z
      discard i
      discard input
      if zb == nil:
        return x.oneLike
      zb

    let identityFunc = Gfunc(
      forward: identityForward,
      backward: identityBackward,
      name: "zero-backed identity")
    let z = graphNode(Gscalar(x.zeroLike), @[Gvalue(x)], identityFunc, "zero-backed identity")

    check not z.isStaticZeroLeaf
    z :~ a
    z.grad(x) :~ 1.0

  test "conds":
    let k = grt.toGvalue(1.0)
    let z = cond(k, x, y)
    let dx = z.grad x
    let dy = z.grad y
    z :~ a
    dx :~ 1.0
    dy :~ 0.0
    k.update 0.0
    z :~ b
    dx :~ 0.0
    dy :~ 1.0

  test "condi 2":
    let k = grt.toGvalue(0)
    let z = cond(k, x, y)
    let z2 = z*z
    let dx = z2.grad x
    let dy = z2.grad y
    z2 :~ b*b
    dx :~ 0.0
    dy :~ 2.0*b
    k.update 1
    y.update c
    z2 :~ a*a
    dx :~ 2.0*a
    dy :~ 0.0
    k.update 0
    z2 :~ c*c
    dx :~ 0.0
    dy :~ 2.0*c

  test "conds 2":
    let k = grt.toGvalue(1.0)
    let z = cond(k, x, y)
    let z2 = z*z
    let dx = z2.grad x
    let dy = z2.grad y
    z2 :~ a*a
    dx :~ 2.0*a
    dy :~ 0.0
    k.update 0.0
    x.update d
    z2 :~ b*b
    dx :~ 0.0
    dy :~ 2.0*b
    k.update 1.0
    x.update c
    z2 :~ c*c
    dx :~ 2.0*c
    dy :~ 0.0

  test "cond gradients stay live across selector flips":
    let xLive = grt.toGvalue(2.0)
    let yLive = grt.toGvalue(3.0)
    let kLive = grt.toGvalue(0)
    let zLive = cond(kLive, xLive, yLive)
    let dxLive = zLive.grad xLive
    let dyLive = zLive.grad yLive

    zLive :~ 3.0
    dxLive :~ 0.0
    dyLive :~ 1.0

    kLive.update 1
    zLive :~ 2.0
    dxLive :~ 1.0
    dyLive :~ 0.0

  test "typed newCondNode preserves concrete result and gradients":
    let selector = grt.toGvalue(1)
    let trueBranch = grt.toGvalue(5.0)
    let falseBranch = grt.toGvalue(7.0)
    let z: Gscalar = newCondNode(selector, trueBranch, falseBranch)
    let dzTrue = z.grad trueBranch
    let dzFalse = z.grad falseBranch

    z :~ 5.0
    dzTrue :~ 1.0
    dzFalse :~ 0.0

    selector.update 0
    z :~ 7.0
    dzTrue :~ 0.0
    dzFalse :~ 1.0

  test "cond backward rejects malformed raw input index":
    let selector = grt.toGvalue(1)
    let trueBranch = grt.toGvalue(5.0)
    let falseBranch = grt.toGvalue(7.0)
    let z = cond(selector, trueBranch, falseBranch)

    expect(GraphValueError):
      discard z.gfunc.backward(nil, z, 3, trueBranch)

  test "literal cond overloads preserve concrete branch type":
    let ks = grt.toGvalue(1.0)
    let ki = grt.toGvalue(1)
    let sx = grt.toGvalue(2.0)
    let ix = grt.toGvalue(2)

    let scalarRight: Gscalar = cond(ks, sx, 0.0)
    let scalarLeft: Gscalar = cond(ki, 0.0, sx)
    let scalarIntLiteral: Gscalar = cond(ks, sx, 1)
    let intRight: Gint = cond(ki, ix, 0)
    let intLeft: Gint = cond(ks, 0, ix)

    scalarRight :~ 2.0
    scalarLeft :~ 0.0
    scalarIntLiteral :~ 2.0
    intRight :~ 2
    intLeft :~ 0

    scalarRight.grad(sx) :~ 1.0
    ks.update 0.0
    scalarRight :~ 0.0
    scalarRight.grad(sx) :~ 0.0

  test "cond rejects malformed inputs early":
    let scalarBranch: Gvalue = grt.toGvalue(2.0)
    let intBranch: Gvalue = grt.toGvalue(3)

    expect(GraphValueError):
      discard cond(grt.toGvalue(1), scalarBranch, intBranch)

  test "cond checks branch compatibility in both directions":
    let selector = grt.toGvalue(1)
    let wide = GcondWide(runtime: grt)
    let narrow = GcondNarrow(runtime: grt)

    expect(GraphValueError):
      discard newCondNode(selector, Gvalue(narrow), Gvalue(wide))

    expect(GraphValueError):
      discard newCondNode(selector, Gvalue(wide), Gvalue(narrow))

  test "cond eval shortcut":
    let t = grt.toGvalue(2.0)
    let f = grt.toGvalue(0.0)
    let t2 = t*t
    let t3 = t*t*t
    check t2.sval == 0.0
    check t3.sval == 0.0
    var tt = cond(t, t2, t3)
    tt :~ 4.0
    check t2.sval == 4.0
    check t3.sval == 0.0
    tt = cond(f, t3, t2)
    tt :~ 4.0
    check t2.sval == 4.0
    check t3.sval == 0.0

  test "treeRepr mode follows eval and input views":
    let k = grt.toGvalue(1)
    let t = grt.toGvalue(17.0)
    let f = grt.toGvalue(23.0)
    let z = cond(k, t + 1.0, f + 1.0)
    let evalTree = z.treeRepr(iwmEval)
    let reachableTree = z.treeRepr(iwmReachable)
    let defaultTree = z.treeRepr

    check evalTree.contains("17.0")
    check not evalTree.contains("23.0")
    check reachableTree.contains("17.0")
    check reachableTree.contains("23.0")
    check defaultTree == reachableTree

  test "eval tree follows cached cond selector without evaluating it":
    let k = grt.toGvalue(1)
    let selector = equal(k, 1)
    let t = grt.toGvalue(17.0)
    let f = grt.toGvalue(23.0)
    let z = cond(selector, t + 1.0, f + 1.0)

    let staleTree = z.treeRepr(iwmEval)
    check staleTree.contains("23.0")
    check not staleTree.contains("17.0")

    discard selector.eval
    let currentTree = z.treeRepr(iwmEval)
    check currentTree.contains("17.0")
    check not currentTree.contains("23.0")

  test "walkInputView exposes cond mode-specific inputs":
    let k = grt.toGvalue(1)
    let t = grt.toGvalue(17.0)
    let f = grt.toGvalue(23.0)
    let trueBranch = t + 1.0
    let falseBranch = f + 1.0
    let z = cond(k, trueBranch, falseBranch)

    proc walked(mode: InputWalkMode): seq[Gvalue] =
      var deps: seq[Gvalue] = @[]
      z.walkInputView(mode, proc(child: Gvalue) =
        deps.add child)
      deps

    let evalDeps = walked(iwmEval)
    check evalDeps.len == 2
    check evalDeps[0].nodeKey == k.nodeKey
    check evalDeps[1].nodeKey == trueBranch.nodeKey

    let reachableDeps = walked(iwmReachable)
    check reachableDeps.len == 3
    check reachableDeps[0].nodeKey == k.nodeKey
    check reachableDeps[1].nodeKey == trueBranch.nodeKey
    check reachableDeps[2].nodeKey == falseBranch.nodeKey

  test "reachable tree includes cond branches":
    let k = grt.toGvalue(1)
    let t = grt.toGvalue(17.0)
    let f = grt.toGvalue(23.0)
    let z = cond(k, t + 1.0, f + 1.0)
    let reachableTree = z.treeRepr(iwmReachable)

    check reachableTree.contains("17.0")
    check reachableTree.contains("23.0")

  test "cond grad build does not eval value or derivative graphs":
    let x2 = grt.toGvalue(2.0)
    let cnd = equal(x2 - 1.0, 0.0)
    let t2 = exp(x2)
    let f2 = x2 * x2 * x2
    let z = cond(cnd, t2, f2)

    let condRuns0 = cnd.runCount
    let subRuns0 = cnd.inputs[0].runCount
    let expRuns0 = t2.runCount
    let mulRuns0 = f2.runCount

    let dzdx = z.grad x2
    discard dzdx.grad x2

    check cnd.runCount == condRuns0
    check cnd.inputs[0].runCount == subRuns0
    check t2.runCount == expRuns0
    check f2.runCount == mulRuns0

  test "inactive true cond branch derivative is not evaluated":
    let x2 = grt.toGvalue(3.0)
    let k = grt.toGvalue(0)
    let z = cond(k, 1.0 / (x2 - 3.0), 0.0)

    z :~ 0.0
    z.grad(x2) :~ 0.0

  test "inactive false cond branch derivative is not evaluated":
    let x2 = grt.toGvalue(3.0)
    let k = grt.toGvalue(1)
    let z = cond(k, 0.0, 1.0 / (x2 - 3.0))

    z :~ 0.0
    z.grad(x2) :~ 0.0

  test "inactive cond branch derivative is not evaluated with upstream":
    let x2 = grt.toGvalue(3.0)
    let k = grt.toGvalue(0)
    let z = cond(k, 1.0 / (x2 - 3.0), 0.0)
    let dep = 2.0 * z

    dep :~ 0.0
    dep.grad(x2) :~ 0.0

  test "seeded cond gradient is lazy zero when target is unreachable":
    let x2 = grt.toGvalue(2.0)
    let kSource = grt.toGvalue(1.0)
    let cnd = equal(kSource, 1.0)
    let z = cond(cnd, grt.toGvalue(17.0), grt.toGvalue(23.0))
    let seed = grt.toGvalue(3.0)

    let cndRuns0 = cnd.runCount
    let dzdx = z.gradSeeded(x2, seed)

    check not dzdx.isCondNode
    dzdx :~ 0.0
    check cnd.runCount == cndRuns0

  test "seeded cond gradient is lazy zero for selector-only target":
    let x2 = grt.toGvalue(2.0)
    let cnd = equal(x2, 0.0)
    let z = cond(cnd, grt.toGvalue(17.0), grt.toGvalue(23.0))
    let seed = grt.toGvalue(3.0)

    let cndRuns0 = cnd.runCount
    let dzdx = z.gradSeeded(x2, seed)

    check not dzdx.isCondNode
    dzdx :~ 0.0
    check cnd.runCount == cndRuns0

  test "seeded cond gradient scales selected branch and follows selector flips":
    let x2 = grt.toGvalue(2.0)
    let y2 = grt.toGvalue(5.0)
    let k = grt.toGvalue(1)
    let seed = grt.toGvalue(3.0)
    let z = cond(k, x2 * x2, y2 * y2)
    let dzdx = z.gradSeeded(x2, seed)
    let dzdy = z.gradSeeded(y2, seed)

    dzdx :~ 12.0
    dzdy :~ 0.0

    k.update 0
    dzdx :~ 0.0
    dzdy :~ 30.0

    seed.update 4.0
    dzdy :~ 40.0

  test "cond backward reuses cached branch adjoint as branch upstream":
    let x2 = grt.toGvalue(2.0)
    let k = grt.toGvalue(1)
    let p = x2 * x2
    let z: Gscalar = cond(k, p, 0.0)

    z.grad(p) :~ 1.0
    z.grad(x2) :~ 4.0

  test "cond backward does not double count cached intermediate adjoint":
    let y2 = grt.toGvalue(2.0)
    let c2 = grt.toGvalue(0.5)
    let k = grt.toGvalue(1)
    let v = y2 * y2
    let dep = (y2 - c2) * (v - cond(k, v + v, v))

    dep.grad(v) :~ -(2.0 - 0.5)
    dep.grad(y2) :~ -(2.0 * 2.0) - (2.0 - 0.5) * 2.0 * 2.0

  test "cond backward caches complete intermediate adjoints":
    let y2 = grt.toGvalue(0.5)
    let c2 = grt.toGvalue(3.0)
    let k = grt.toGvalue(0)
    let e = y2 * y2
    let denom = e * e + 1.0
    let dep = y2 / denom - exp(cond(k, c2, e))
    let ev = 0.5 * 0.5
    let expected = -2.0 * 0.5 * ev / ((ev * ev + 1.0) * (ev * ev + 1.0)) -
      exp(ev)

    discard dep.grad(y2).eval
    dep.grad(e) :~ expected

  test "cond eval short-circuits value and gradient graphs":
    block:
      let x2 = grt.toGvalue(2.0)
      let cnd = equal(x2 - 1.0, 0.0)
      let t2 = exp(x2)
      let f2 = x2 * x2 * x2
      let z = cond(cnd, t2, f2)
      let condRuns0 = cnd.runCount
      let subRuns0 = cnd.inputs[0].runCount
      let expRuns0 = t2.runCount
      let mulRuns0 = f2.runCount

      z :~ 8.0

      check cnd.runCount > condRuns0
      check cnd.inputs[0].runCount > subRuns0
      check t2.runCount == expRuns0
      check f2.runCount > mulRuns0

    block:
      let x2 = grt.toGvalue(2.0)
      let cnd = equal(x2 - 1.0, 0.0)
      let t2 = exp(x2)
      let f2 = x2 * x2 * x2
      let z = cond(cnd, t2, f2)
      let dzdx = z.grad x2
      let condRuns0 = cnd.runCount
      let subRuns0 = cnd.inputs[0].runCount
      let expRuns0 = t2.runCount
      let mulRuns0 = f2.runCount

      dzdx :~ 12.0

      check cnd.runCount > condRuns0
      check cnd.inputs[0].runCount > subRuns0
      check t2.runCount == expRuns0
      check f2.runCount == mulRuns0

    block:
      let x2 = grt.toGvalue(1.0)
      let cnd = equal(x2 - 1.0, 0.0)
      let t2 = exp(x2)
      let f2 = x2 * x2 * x2
      let z = cond(cnd, t2, f2)
      let condRuns0 = cnd.runCount
      let subRuns0 = cnd.inputs[0].runCount
      let expRuns0 = t2.runCount
      let mulRuns0 = f2.runCount

      z :~ exp(1.0)

      check cnd.runCount > condRuns0
      check cnd.inputs[0].runCount > subRuns0
      check t2.runCount > expRuns0
      check f2.runCount == mulRuns0

    block:
      let x2 = grt.toGvalue(1.0)
      let cnd = equal(x2 - 1.0, 0.0)
      let t2 = exp(x2)
      let f2 = x2 * x2 * x2
      let z = cond(cnd, t2, f2)
      let dzdx = z.grad x2
      let condRuns0 = cnd.runCount
      let subRuns0 = cnd.inputs[0].runCount
      let expRuns0 = t2.runCount
      let mulRuns0 = f2.runCount

      dzdx :~ exp(1.0)

      check cnd.runCount > condRuns0
      check cnd.inputs[0].runCount > subRuns0
      check t2.runCount > expRuns0
      check f2.runCount == mulRuns0

  test "higher-order grad through cond":
    let x3 = grt.toGvalue(2.0)
    let c3 = equal(x3 - 1.0, 0.0)
    let z3 = cond(c3, x3 * x3, x3 * x3 * x3)
    let dzdx3 = z3.grad x3
    let d2zdx23 = dzdx3.grad x3

    dzdx3 :~ 12.0
    d2zdx23 :~ 12.0

    x3.update 1.0
    dzdx3 :~ 2.0
    d2zdx23 :~ 2.0
