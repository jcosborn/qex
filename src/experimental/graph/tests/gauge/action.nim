suite "gauge action":
  let gplaq = block:
    var pl = 0.0
    for t in g.plaq:
      pl += t
    pl

  setup:
    let gg {.used.} = grt.toGvalue(g)
    let gu {.used.} = grt.toGvalue(u)
    let gm {.used.} = grt.toGvalue(m)

  test "wilson action":
    let beta = 5.4
    let c = grt.actWilson(beta)
    let s = gaugeAction(c, gg)
    s :~ -gplaq*float(6*vol*beta)
    proc act(x: Gvalue): Gvalue = gaugeAction(c, x)
    ckgrad(act, gg, gu)

  test "wilson force":
    let beta = 5.4
    let c = grt.actWilson(beta)
    proc act(x: Gvalue): Gvalue = gaugeAction(c, x)
    proc force(x: Gvalue): Gvalue = gaugeForce(c, x)
    ckforce(act, force, gg, 10.0*gm)

  test "wilson force gradient":
    let beta = 5.4
    let c = grt.actWilson(beta)
    proc force(x: Gvalue): Gvalue = gaugeForce(c, x)
    ckgradm(force, gg, gu, gm)

  test "wilson force gradient recomp":
    let beta = 5.4
    let c = grt.actWilson(beta)
    let a = gaugeAction(c, gg)
    let f2 = gaugeForce(c, gg).norm2
    let df2 = grad(f2, gg).norm2
    let rs1 = [a.eval.getfloat, f2.eval.getfloat, df2.eval.getfloat]
    c.updated
    gg.updated
    let rs2 = [a.eval.getfloat, f2.eval.getfloat, df2.eval.getfloat]
    c.updated
    gg.updated
    let rs3 = [a.eval.getfloat, f2.eval.getfloat, df2.eval.getfloat]
    check rs1 == rs2
    check rs1 == rs3

  test "gaugeAction rejects coefficient gradients through the action layer":
    let beta = grt.toGvalue(5.4)
    let c = actWilson(beta)
    expect(GraphValueError):
      discard gaugeAction(c, gg).grad beta

  test "gaugeActionDeriv rejects coefficient gradients through the action layer":
    let beta = grt.toGvalue(5.4)
    let c = actWilson(beta)
    expect(GraphValueError):
      discard gaugeActionDeriv(c, gg).grad beta
