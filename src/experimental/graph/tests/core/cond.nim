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

    expect(GraphValueError):
      discard newCondNode(nil, scalarBranch, scalarBranch)
    expect(GraphValueError):
      discard newCondNode(grt.toGvalue(1), nil, scalarBranch)

  test "walkedInputs rejects nil callbacks":
    expect(GraphValueError):
      discard walkedInputs(nil)

  test "cond eval shortcut":
    let t = grt.toGvalue(2.0)
    let f = grt.toGvalue(0.0)
    let t2 = t*t
    let t3 = t*t*t
    check t2.getfloat == 0.0
    check t3.getfloat == 0.0
    var tt = cond(t, t2, t3)
    tt :~ 4.0
    check t2.getfloat == 4.0
    check t3.getfloat == 0.0
    tt = cond(f, t3, t2)
    tt :~ 4.0
    check t2.getfloat == 4.0
    check t3.getfloat == 0.0

  test "treeRepr mode follows eval and dependency walks":
    let k = grt.toGvalue(1)
    let t = grt.toGvalue(17.0)
    let f = grt.toGvalue(23.0)
    let z = cond(k, t + 1.0, f + 1.0)
    let evalTree = z.treeRepr(iwmEval)
    let dependTree = z.treeRepr(iwmDepend)

    check evalTree.contains("17.0")
    check not evalTree.contains("23.0")
    check dependTree.contains("17.0")
    check dependTree.contains("23.0")

  test "signature tree omits cond branches":
    let k = grt.toGvalue(1)
    let t = grt.toGvalue(17.0)
    let f = grt.toGvalue(23.0)
    let z = cond(k, t + 1.0, f + 1.0)
    let sigTree = z.treeRepr(iwmGradSignature)

    check not sigTree.contains("17.0")
    check not sigTree.contains("23.0")

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
      check f2.runCount > mulRuns0

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
