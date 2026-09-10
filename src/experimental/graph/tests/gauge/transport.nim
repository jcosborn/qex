suite "gauge transport":
  setup:
    let gg {.used.} = grt.toGvalue(g)
    let gu {.used.} = grt.toGvalue(u)
    let gm {.used.} = grt.toGvalue(m)

  test "hop transport is unitary":
    let f = linkField(gu, 0)
    for mu in 0..<g.len:
      norm2(hop(gg, hop(gg, f, mu, 1), mu, -1) - f) :< 1e-18
      norm2(hop(gg, hop(gg, f, mu, -1), mu, 1) - f) :< 1e-18

  test "hop gradients":
    proc t1(x: Ggauge): Gscalar = retr(hop(x, linkField(x, 1), 0, 1))
    ckgrad(t1, gg, gu)
    proc t2(x: Ggauge): Gscalar = redot(linkField(gm, 2), hop(x, linkField(x, 1), 2, -1))
    ckgrad(t2, gg, gu)

  test "transport inputs fail before allocating path nodes":
    let
      f = linkField(gu, 0)
      nextNode = grt.nextStableNodeId
    expect(GraphValueError):
      discard hop(gg, f, 0, 0)
    expect(GraphValueError):
      discard hop(gg, f, g.len, 1)
    expect(GraphValueError):
      discard transport(gg, f, [0])
    expect(GraphValueError):
      discard transport(gg, f, [g.len + 1])
    expect(GraphValueError):
      discard transport(gg, f, [-g.len - 1])
    expect(GraphValueError):
      discard plaqPath(-1, 0)
    expect(GraphValueError):
      discard plaqPath(1, 1)
    check grt.nextStableNodeId == nextNode

  test "functional clones own shift and hop workspaces":
    proc shiftedHop(x: Ggauge): Gscalar =
      redot(
        shift(linkField(x, 0), 1, 1),
        hop(x, linkField(x, 2), 3, -1))
    let
      xg = grt.toGvalue(g)
      xu = grt.toGvalue(u)
      param = Ggauge(xg.newOneOf)
      fn = lambda(param, shiftedHop(param))
      gotG = Gscalar(apply(fn, xg))
      gotU = Gscalar(apply(fn, xu))
      refG = shiftedHop(xg)
      refU = shiftedHop(xu)
      gotGradG = grad(gotG, xg)
      refGradG = grad(refG, xg)
    (gotG - refG) :< 1e-8
    (gotU - refU) :< 1e-8
    norm2(gotGradG - refGradG) :< 1e-20
    xg.update u
    xu.update g
    (gotU - refU) :< 1e-8
    (gotG - refG) :< 1e-8
    norm2(gotGradG - refGradG) :< 1e-20

  test "wilson line gradients":
    proc w1(x: Ggauge): Gscalar = retr(wilsonLine(x, plaqPath(0, 1)))
    ckgrad(w1, gg, gu)
    proc w2(x: Ggauge): Gscalar = retr(wilsonLine(x, [1, 2, 3, -1, -2, -3]))
    ckgrad(w2, gg, gu)
    proc w3(x: Ggauge): Gscalar = redot(linkField(gm, 0), wilsonLine(x, [4, 1, -4]))
    ckgrad(w3, gg, gu)
