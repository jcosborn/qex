proc mixedScalarGaugeForward(v: Gvalue) =
  let z = v.requireMultiValue("mixed scalar/gauge forward")
  for i in 0..<z.inputs.len:
    z.storedSlot(i).valCopy(z.inputs[i])

proc mixedScalarGaugeBackward(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard dep
  let upstream = requireMultiUpstream(zb, "mixed scalar/gauge backward")
  if i < 0 or i >= z.inputs.len:
    raiseValueError("mixed scalar/gauge backward input index out of range: " & $i)
  upstream[i]

let mixedScalarGaugeFunc = newGfunc(
  forward = mixedScalarGaugeForward,
  backward = mixedScalarGaugeBackward,
  name = "mixedScalarGauge")

proc newScalarGaugeMulti(s: Gscalar, g: Ggauge, label: string): Gmulti =
  newMultiOutputNode([Gvalue(s), Gvalue(g)], [Gvalue(s), Gvalue(g)], mixedScalarGaugeFunc, label)

suite "gauge basic":
  setup:
    let gg {.used.} = grt.toGvalue(g)
    let gu {.used.} = grt.toGvalue(u)
    let gp {.used.} = grt.toGvalue(p)
    let gq {.used.} = grt.toGvalue(q)
    let gm {.used.} = grt.toGvalue(m)
    let x {.used.} = grt.toGvalue(a)
    let y {.used.} = grt.toGvalue(b)

  test "norm2":
    let n2 = gg.norm2
    let p2 = gp.norm2
    let dp = grad(0.5 * p2, gp)
    n2 :~ 4.0*float(nc*vol)
    dp.norm2 :~ p2
    norm2(dp-gp) :~ 0
    ckgrad(norm2, gm, gq)

  test "redot":
    let n2 = gg.redot gg
    let p2 = gp.redot gp
    let dp = grad(0.5 * p2, gp)
    n2.eval :~ 4.0*float(nc*vol)
    dp.norm2 :~ p2
    norm2(dp-gp).eval :~ 0
    let pq = gp.redot gq
    norm2(grad(pq, gp) - gq) :< 1e-26
    norm2(grad(pq, gq) - gp) :< 1e-26
    ckgrad2(redot, gp, gq, gg, gu)

  test "g+g propagates gauge upstream to both inputs":
    let z = redot(gp + gq, gu)

    norm2(grad(z, gp) - gu) :< 1e-26
    norm2(grad(z, gq) - gu) :< 1e-26

  test "gauge operators require explicit erased operand casts":
    let erasedGauge: Gvalue = gq
    let erasedScalar: Gvalue = x
    let castGauge = erasedGauge.requireGauge("erased gauge test right")
    let castScalar = erasedScalar.requireScalar("erased gauge test scalar")

    check not compiles(gp + erasedGauge)
    check not compiles(gp * erasedGauge)
    check not compiles(gp - erasedGauge)
    check not compiles(gp - erasedScalar)
    check not compiles(redot(gp, erasedGauge))
    check not compiles(expDeriv(gp, erasedGauge))

    let addResult: Ggauge = gp + castGauge
    let mulResult: Ggauge = gp * castGauge
    let subGaugeResult: Ggauge = gp - castGauge
    let subScalarResult: Ggauge = gp - castScalar
    let dotResult: Gscalar = redot(gp, castGauge)
    let derivResult: Ggauge = expDeriv(gp, castGauge)

    norm2(addResult - (gp + gq)) :< 1e-26
    norm2(mulResult - (gp * gq)) :< 1e-26
    norm2(subGaugeResult - (gp - gq)) :< 1e-26
    norm2(subScalarResult - (gp - x)) :< 1e-26
    dotResult :~ redot(gp, gq)
    norm2(derivResult - expDeriv(gp, gq)) :< 1e-26

  test "gauge validator rejects wrong value type":
    let scalarValue: Gvalue = x
    expect(GraphValueError):
      discard scalarValue.requireGauge("gaugeSnapshot")

  test "gauge gradient add labels wrong operand side":
    try:
      discard gg.addLike(Gvalue(x), Gvalue(gg))
      fail()
    except GraphValueError as e:
      check e.msg.contains("gauge gradient add left")

    try:
      discard gg.addLike(Gvalue(gg), Gvalue(x))
      fail()
    except GraphValueError as e:
      check e.msg.contains("gauge gradient add right")

  test "explicit gauge forward hooks reject wrong result node type":
    proc wrongGaugeResultForward(v: Gvalue) =
      discard v.requireUnaryNodeView(Ggauge, "wrong gauge result forward")
      discard v.requireGauge("wrong gauge result forward result")

    let wrongGaugeResultFunc = newGfunc(
      forward = wrongGaugeResultForward,
      name = "wrongGaugeResult")
    let malformed = graphNode(
      scalarNodeLike(x),
      @[Gvalue(gg)],
      wrongGaugeResultFunc,
      "wrong gauge result")

    try:
      discard malformed.eval
      fail()
    except GraphValueError as e:
      check e.msg.contains("wrong gauge result forward")
      check e.msg.contains("expects gauge value")

  test "toGvalue owns a snapshot of caller gauge storage":
    var localGauge = lo.newgauge
    threads:
      for mu in 0..<localGauge.len:
        localGauge[mu] := g[mu]

    let wrapped = grt.toGvalue(localGauge)
    let snapshot = wrapped.gaugeNodeLike
    snapshot.valCopy(wrapped)
    snapshot.updated

    threads:
      for mu in 0..<localGauge.len:
        localGauge[mu] := 0.0

    norm2(wrapped - snapshot) :< 1e-26

  test "gaugeSnapshot returns a snapshot of graph gauge storage":
    let wrapped = grt.toGvalue(g)
    let before = wrapped.norm2.eval.sval
    var snapshot = wrapped.gaugeSnapshot
    snapshot.zeroGaugeStorage

    wrapped.norm2 :~ before

  test "mutateGauge marks graph-owned gauge storage fresh":
    let wrapped = grt.toGvalue(g)
    let n2 = wrapped.norm2
    n2 :~ 4.0*float(nc*vol)

    mutateGauge(wrapped, storage):
      storage.zeroGaugeStorage

    n2 :~ 0.0

  test "gauge copy paths reject incompatible shapes before copying":
    let lo2 = @[4,4,4,8].newLayout
    let g2 = lo2.newgauge
    let left = grt.toGvalue(zeroGaugeLike(g))
    let right = grt.toGvalue(zeroGaugeLike(g2))

    expect(GraphValueError):
      left.update(g2)
    expect(GraphValueError):
      left.valCopy(right)

  test "binary gauge ops reject incompatible layouts at construction":
    let lo2 = @[4,4,4,8].newLayout
    let g2 = lo2.newgauge
    let other = grt.toGvalue(zeroGaugeLike(g2))
    let erasedOther: Gvalue = other

    expect(GraphValueError):
      discard gg + other
    expect(GraphValueError):
      discard gg * other
    expect(GraphValueError):
      discard gg - other
    expect(GraphValueError):
      discard redot(gg, other)
    expect(GraphValueError):
      discard expDeriv(gg, other)
    expect(GraphValueError):
      discard gg + erasedOther.requireGauge("erased shape right")

  test "gauge numeric literal overloads stay explicit":
    let one = grt.toGvalue(1.0)
    let two = grt.toGvalue(2.0)
    let shiftedFloat: Ggauge = gp - 1.0
    let shiftedInt: Ggauge = gp - 1
    let addedFloat: Ggauge = 1.0 + gp
    let addedInt: Ggauge = 1 + gp
    let scaledFloat: Ggauge = 2.0 * gp
    let scaledInt: Ggauge = 2 * gp

    norm2(shiftedFloat - (gp - one)) :< 1e-26
    norm2(shiftedInt - (gp - one)) :< 1e-26
    norm2(addedFloat - (one + gp)) :< 1e-26
    norm2(addedInt - (one + gp)) :< 1e-26
    norm2(scaledFloat - (two * gp)) :< 1e-26
    norm2(scaledInt - (two * gp)) :< 1e-26

  test "retr":
    let rtp = gp.retr
    let n2 = retr(gg * gg.adj)
    rtp*rtp :< 1e-20
    n2.eval :~ 4.0*float(nc*vol)
    let p2 = retr(gp * gq.adj)
    p2 :~ redot(gp, gq)
    norm2(grad(p2, gp) - gq) :< 1e-26
    norm2(grad(p2, gq) - gp) :< 1e-26
    ckgrad(retr, gp, gq)

  test "adj":
    norm2(gg.adj*gg - 1.0)/float(4*nc*vol) :< 1e-22
    norm2(gg*gg.adj - 1.0)/float(4*nc*vol) :< 1e-22
    norm2(gp.adj + gp) :< 1e-26
    norm2(grad(gp.adj.norm2, gp) - 2.0*gp) :< 1e-26
    ckgradm(adj, gg, gp, gq)

  test "neg":
    norm2(gp.adj - (-gp)) :< 1e-26
    norm2(-gp) :~ gp.norm2
    ckgradm(`-`, gg, gp, gq)

  test "addsg":
    let p2 = norm2(x+gp)
    grad(p2, x) :~ retr(2.0*(a+gp))
    norm2(grad(p2, gp) - 2.0*(a+gp)) :< 1e-26
    ckgradm2(`+`, x, gp, y, gq, gg)

  test "addgg":
    let pq = norm2(gp+gq)
    norm2(grad(pq, gp) - 2.0*(gp+gq)) :< 1e-26
    norm2(grad(pq, gq) - 2.0*(gp+gq)) :< 1e-26
    ckgradm2(`+`, gq, gp, gu, gg, gm)

  test "mulsg":
    let p2 = norm2(x*gp)
    grad(p2, x) :~ 2.0*a*gp.norm2
    norm2(grad(p2, gp) - 2.0*a*a*gp) :< 1e-26
    ckgradm2(`*`, x, gp, y, gq, gg)

  test "mulgg":
    let pq = norm2(gp*gq)
    norm2(grad(pq, gp) - 2.0*gp*gq*gq.adj) :< 1e-24
    norm2(grad(pq, gq) - 2.0*gp.adj*gp*gq) :< 1e-24
    ckgradm2(`*`, gq, gp, gu, gg, gm)

  test "subgs":
    let p2 = norm2(gp-x)
    grad(p2, x) :~ retr(-2.0*(gp-a))
    norm2(grad(p2, gp) - 2.0*(gp-x)) :< 1e-26
    ckgradm2(`-`, gp, x, gq, y, gg)

  test "subgg":
    let pq = norm2(gp-gq)
    norm2(grad(pq, gp) - 2.0*(gp-gq)) :< 1e-26
    norm2(grad(pq, gq) - 2.0*(gq-gp)) :< 1e-26
    ckgradm2(`-`, gq, gp, gu, gg, gm)

  test "exp":
    let egp = exp(gp)
    norm2(egp.adj*egp - 1.0) :< 1e-20
    norm2(egp*egp.adj - 1.0) :< 1e-20
    ckgradm(exp, gm, 0.1*gp, gg)

  test "expDeriv backward reports unsupported paths for both inputs":
    let loss = expDeriv(gp, gm).norm2

    try:
      discard loss.grad gp
      fail()
    except GraphValueError as e:
      check e.msg.contains("force-direction input")

    try:
      discard loss.grad gm
      fail()
    except GraphValueError as e:
      check e.msg.contains("matrix exponential")

  test "projTAH":
    let gt = gg.projTAH
    let tgt = gt.retr
    tgt*tgt :< 1e-26
    ckgradm(projTAH, gg, gp, gu)

  test "zero-valued gauges carry zero buffers through copies and arithmetic":
    let zero = grt.toGvalue(zeroGaugeLike(g))
    let zeroCopy = Ggauge(zero.newOneOf)
    zeroCopy.valCopy(zero)

    zero.norm2 :~ 0.0
    zeroCopy.norm2 :~ 0.0
    norm2((zero + gg) - gg) :< 1e-26

  test "cond rejects incompatible gauge layouts":
    let lo2 = @[4,4,4,8].newLayout
    let g2 = lo2.newgauge

    expect(GraphValueError):
      discard cond(
        grt.toGvalue(1),
        grt.toGvalue(zeroGaugeLike(g)),
        grt.toGvalue(zeroGaugeLike(g2)))

  test "multi output preserves heterogeneous gauge layouts":
    let lo2 = @[4,4,4,8].newLayout
    let g2 = lo2.newgauge
    let left = grt.toGvalue(zeroGaugeLike(g))
    let right = grt.toGvalue(zeroGaugeLike(g2))
    let slots = [Gvalue(left), Gvalue(right)]
    let mixed = newMultiOutputNode(slots, newSeq[Gvalue](0), nil, "mixed gauge")
    let first = Ggauge(mixed.storedSlot(0))
    let second = Ggauge(mixed.storedSlot(1))

    check first.copyCompatible(left)
    check second.copyCompatible(right)
    check not first.copyCompatible(right)
    check not second.copyCompatible(left)

  test "multi add supports heterogeneous scalar and gauge slots":
    let left = newScalarGaugeMulti(x, gp, "mixed left")
    let right = newScalarGaugeMulti(y, gq, "mixed right")
    let added: Gmulti = left + right
    let addedGauge: Ggauge = Ggauge(added[1])

    added[0] :~ a + b
    norm2(addedGauge - (gp + gq)) :< 1e-26
    grad(added[0], x) :~ 1.0
    grad(added[0], y) :~ 1.0
