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

  test "projTAH":
    let gt = gg.projTAH
    let tgt = gt.retr
    tgt*tgt :< 1e-26
    ckgradm(projTAH, gg, gp, gu)

  test "zero-valued gauges carry zero buffers through copies and arithmetic":
    let zero = grt.toGvalue(g, isZero = true)
    let zeroCopy = Ggauge(zero.newOneOf)
    zeroCopy.valCopy(zero)

    zero.norm2 :~ 0.0
    zeroCopy.norm2 :~ 0.0
    norm2((zero + gg) - gg) :< 1e-26

  test "cond rejects incompatible gauge layouts":
    let lo2 = @[4,4,4,8].newLayout
    let g2 = lo2.newgauge

    expect(GraphValueError):
      discard cond(grt.toGvalue(1), grt.toGvalue(g, isZero = true), grt.toGvalue(g2, isZero = true))

  test "multi output preserves heterogeneous gauge layouts":
    let lo2 = @[4,4,4,8].newLayout
    let g2 = lo2.newgauge
    let left = grt.toGvalue(g, isZero = true)
    let right = grt.toGvalue(g2, isZero = true)
    let slots = [Gvalue(left), Gvalue(right)]
    let mixed = newMultiOutputNode(slots, slots, nil, "mixed gauge")
    let first = Ggauge(mixed.slotValue(0))
    let second = Ggauge(mixed.slotValue(1))

    check first.copyCompatible(left)
    check second.copyCompatible(right)
    check not first.copyCompatible(right)
    check not second.copyCompatible(left)
