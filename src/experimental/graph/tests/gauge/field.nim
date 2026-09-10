suite "gauge field ops":
  setup:
    let gg {.used.} = grt.toGvalue(g)
    let gu {.used.} = grt.toGvalue(u)
    let gm {.used.} = grt.toGvalue(m)
    let gq {.used.} = grt.toGvalue(q)

  test "direction extraction and injection":
    var s = retr(linkField(gg, 0))
    for mu in 1..<g.len:
      s = s + retr(linkField(gg, mu))
    (s - retr(gg)) :< 1e-6
    let f = linkField(gg, 1)
    let inj = injectLink(f, 1, gg)
    (retr(inj) - retr(f)) :< 1e-8
    (norm2(inj) - norm2(f)) :< 1e-8

  test "field directions are validated before storage access":
    let f = linkField(gg, 0)
    expect(GraphValueError):
      discard linkField(gg, -1)
    expect(GraphValueError):
      discard linkField(gg, g.len)
    expect(GraphValueError):
      discard shift(f, -1, 1)
    expect(GraphValueError):
      discard shift(f, g[0].l.nDim, 1)

  test "field algebra matches gauge algebra per direction":
    let half = scalar.toGvalue(grt, 0.5)
    proc terms(mu: int): array[5, Gscalar] =
      let a = linkField(gg, mu)
      let b = linkField(gu, mu)
      [retr(a * b), redot(a, b), retr(adj(a)), norm2(projTAH(a)),
       retr(half + a)]
    var s = terms(0)
    for mu in 1..<g.len:
      let t = terms(mu)
      for k in 0..<s.len:
        s[k] = s[k] + t[k]
    (s[0] - retr(gg * gu)) :< 1e-6
    (s[1] - redot(gg, gu)) :< 1e-6
    (s[2] - retr(adj(gg))) :< 1e-6
    (s[3] - norm2(projTAH(gg))) :< 1e-6
    (s[4] - retr(half + gg)) :< 1e-6

  test "shift roundtrip, sum invariance, and adjoint":
    let f = linkField(gg, 2)
    for dir in 0..<g.len:
      (retr(shift(f, dir, 1)) - retr(f)) :< 1e-8
      norm2(shift(shift(f, dir, 1), dir, -1) - f) :< 1e-24
    let a = linkField(gg, 0)
    let b = linkField(gu, 1)
    (redot(shift(a, 3, 1), b) - redot(a, shift(b, 3, -1))) :< 1e-8

  test "field op gradients":
    proc fs1(x: Ggauge): Gscalar = retr(linkField(x, 0) * shift(linkField(x, 1), 0, 1))
    ckgrad(fs1, gg, gu)
    proc fs2(x: Ggauge): Gscalar = norm2(projTAH(linkField(x, 2) * adj(shift(linkField(x, 3), 2, 1))))
    ckgrad(fs2, gg, gu)
    proc fs3(x: Ggauge): Gscalar = redot(linkField(x, 0), shift(linkField(x, 1), 2, -1))
    ckgrad(fs3, gg, gu)
    proc fs4(x: Ggauge): Gscalar = norm2(0.5 + linkField(x, 1))
    ckgrad(fs4, gg, gu)

  test "field ops differentiate to third order":
    proc fs(x: Ggauge): Gscalar = norm2(projTAH(linkField(x, 2) * adj(shift(linkField(x, 3), 2, 1))))
    proc f2(x: Ggauge): Gscalar = redot(gm, grad(fs(x), x))
    ckgrad(f2, gg, gu)
    proc f3(x: Ggauge): Gscalar = redot(gq, grad(f2(x), x))
    ckgrad(f3, gg, gu)
