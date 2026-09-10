# Equivalences in the public definition order of gauge/fused_ops.nim.
# Reference expressions use gauge/basic_ops.nim.

suite "gauge fused vs basic":
  setup:
    let gg {.used.} = grt.toGvalue(g)
    let gu {.used.} = grt.toGvalue(u)
    let gp {.used.} = grt.toGvalue(p)
    let gq {.used.} = grt.toGvalue(q)
    let gm {.used.} = grt.toGvalue(m)
    let x {.used.} = grt.toGvalue(a)

  test "axpy(a, x, y) = a * x + y":
    let rf = axpy(x, gm, gg)
    let rg = x * gm + gg
    norm2(rf - rg) :< 1e-26
    let srf = retr(rf * gu)
    let srg = retr(rg * gu)
    grad(srf, x) :~ grad(srg, x)
    norm2(grad(srf, gm) - grad(srg, gm)) :< 1e-26
    norm2(grad(srf, gg) - grad(srg, gg)) :< 1e-26

  test "adjmul(x, y) = adj(x) * y":
    ckbinarynorm2grad(gg.adjmul gu, gg.adj * gu, gg, gu, 1e-25)

  test "muladj(x, y) = x * adj(y)":
    ckbinarynorm2grad(gg.muladj gu, gg * gu.adj, gg, gu, 1e-25)

  test "contractProjTAH(x, y) = projTAH(x * adj(y))":
    ckbinarynorm2grad(contractProjTAH(gg, gu), projTAH(gg * gu.adj), gg, gu, 1e-26)

  test "contractProjTAH subset matches masked basic ops":
    const
      parity = 1
      dir = 2
    let
      zero = grt.toGvalue(zeroGaugeLike(g))
      rf = contractProjTAH(gg, gu, parity, dir)
      rr = blendSubset(parity, dir, projTAH(gg * gu.adj), zero)
      sf = redot(rf, gm)
      sr = redot(rr, gm)
    norm2(rf - rr) :< 1e-26
    norm2(grad(sf, gg) - grad(sr, gg)) :< 1e-26
    norm2(grad(sf, gu) - grad(sr, gu)) :< 1e-26
    gg.update(u)
    norm2(rf - rr) :< 1e-26

  test "contractProjTAH consumes gauge sums without materializing them":
    let
      summed = (gg + gp) + (gq + gp)
      rf = contractProjTAH(summed, gu)
    discard rf.eval
    check summed.runCount == 0

    let
      rr = projTAH(summed * gu.adj)
      sf = redot(rf, gm)
      sr = redot(rr, gm)
    norm2(rf - rr) :< 5e-26
    norm2(grad(sf, summed) - grad(sr, summed)) :< 1e-26
    norm2(grad(sf, gg) - grad(sr, gg)) :< 1e-26
    norm2(grad(sf, gp) - grad(sr, gp)) :< 1e-26
    norm2(grad(sf, gq) - grad(sr, gq)) :< 1e-26
    norm2(grad(sf, gu) - grad(sr, gu)) :< 1e-26
    gg.update(u)
    norm2(rf - rr) :< 5e-26

  test "contractProjTAH subset consumes gauge sums without materializing them":
    const
      parity = 1
      dir = 2
    let
      zero = grt.toGvalue(zeroGaugeLike(g))
      summed = (gg + gp) + gq
      rf = contractProjTAH(summed, gu, parity, dir)
    discard rf.eval
    check summed.runCount == 0

    let
      rr = blendSubset(parity, dir, projTAH(summed * gu.adj), zero)
      sf = redot(rf, gm)
      sr = redot(rr, gm)
    norm2(rf - rr) :< 1e-26
    norm2(grad(sf, summed) - grad(sr, summed)) :< 1e-26
    norm2(grad(sf, gg) - grad(sr, gg)) :< 1e-26
    norm2(grad(sf, gp) - grad(sr, gp)) :< 1e-26
    norm2(grad(sf, gq) - grad(sr, gq)) :< 1e-26
    norm2(grad(sf, gu) - grad(sr, gu)) :< 1e-26

  test "contractProjTAH packed slots stay stable across gauge updates":
    let rf = contractProjTAH(gg, gu)
    let rg = projTAH(gg * gu.adj)

    norm2(rf - rg) :< 1e-26
    gu.update q
    norm2(rf - rg) :< 1e-26

  test "contractProjTAH does not read packed storage before first eval":
    let rf = contractProjTAH(gg, gu)
    let rg = projTAH(gg * gu.adj)

    gu.update q
    norm2(rf - rg) :< 1e-26

  test "contractProjTAH shared backward helper stays correct across outputs":
    let rf = contractProjTAH(gg, gu)
    let rg = projTAH(gg * gu.adj)

    let srf = retr(rf * gp)
    let srg = retr(rg * gp)
    let trf = retr(rf * gq)
    let trg = retr(rg * gq)

    let dsrfgg = grad(srf, gg)
    let dsrggg = grad(srg, gg)
    let dtrfgg = grad(trf, gg)
    let dtrggg = grad(trg, gg)
    let dsrfgu = grad(srf, gu)
    let dsrggu = grad(srg, gu)
    let dtrfgu = grad(trf, gu)
    let dtrggu = grad(trg, gu)

    norm2(dsrfgg - dsrggg) :< 1e-26
    norm2(dtrfgg - dtrggg) :< 1e-26
    norm2(dsrfgu - dsrggu) :< 1e-26
    norm2(dtrfgu - dtrggu) :< 1e-26

  test "axexp(a, x) = exp(a * x)":
    let rf = axexp(x, gm)
    let rg = exp(x * gm)
    norm2(rf - rg) :< 1e-26
    let srf = retr(rf * gu)
    let srg = retr(rg * gu)
    grad(srf, x) :~ grad(srg, x)
    norm2(grad(srf, gm) - grad(srg, gm)) :< 1e-26

  test "axexp packed slots stay stable across scalar updates":
    let rf = axexp(x, gm)
    let rg = exp(x * gm)

    norm2(rf - rg) :< 1e-26
    x.update 0.25
    norm2(rf - rg) :< 1e-26

  test "axexp does not read packed storage before first eval":
    let rf = axexp(x, gm)
    let rg = exp(x * gm)

    x.update 0.25
    norm2(rf - rg) :< 1e-26

  test "axexp shared backward helper stays correct across outputs":
    let rf = axexp(x, gm)
    let rg = exp(x * gm)

    let srf = retr(rf * gu)
    let srg = retr(rg * gu)
    let trf = retr(rf * gp)
    let trg = retr(rg * gp)

    let dsrfx = grad(srf, x)
    let dsrgx = grad(srg, x)
    let dtrfx = grad(trf, x)
    let dtrgx = grad(trg, x)
    let dsrfgm = grad(srf, gm)
    let dsrggm = grad(srg, gm)
    let dtrfgm = grad(trf, gm)
    let dtrggm = grad(trg, gm)

    dsrfx :~ dsrgx
    dtrfx :~ dtrgx
    norm2(dsrfgm - dsrggm) :< 1e-26
    norm2(dtrfgm - dtrggm) :< 1e-26

  test "axexpmuly(a, x, y) = exp(a * x) * y":
    # gm is anti-Hermitian and traceless for SU(N), as axexpmuly requires.
    let rf: Ggauge = axexpmuly(x, gm, gg)
    let rg = exp(x * gm) * gg
    check rf.runtime == grt
    norm2(rf - rg) :< 1e-26
    let srf = retr(rf * gu)
    let srg = retr(rg * gu)
    grad(srf, x) :~ grad(srg, x)
    norm2(grad(srf, gm) - grad(srg, gm)) :< 1e-26
    norm2(grad(srf, gg) - grad(srg, gg)) :< 1e-26

  test "axexpmuly subset matches masked basic ops":
    const
      parity = 1
      dir = 2
    let
      zero = grt.toGvalue(zeroGaugeLike(g))
      rf = axexpmuly(x, gm, gg, parity, dir)
      rr = blendSubset(parity, dir, exp(x * gm) * gg, zero)
      sf = redot(rf, gu)
      sr = redot(rr, gu)
    norm2(rf - rr) :< 1e-26
    grad(sf, x) :~ grad(sr, x)
    norm2(grad(sf, gm) - grad(sr, gm)) :< 1e-26
    norm2(grad(sf, gg) - grad(sr, gg)) :< 1e-26

  test "axexpmuly shared result stays correct across outputs and updates":
    let rf = axexpmuly(x, gm, gg)
    let rg = exp(x * gm) * gg

    let srf = retr(rf * gu)
    let srg = retr(rg * gu)
    let trf = retr(rf * gp)
    let trg = retr(rg * gp)

    let dsrfgg = grad(srf, gg)
    let dsrggg = grad(srg, gg)
    let dtrfgg = grad(trf, gg)
    let dtrggg = grad(trg, gg)
    let dsrfx = grad(srf, x)
    let dsrgx = grad(srg, x)
    let dtrfx = grad(trf, x)
    let dtrgx = grad(trg, x)
    let dsrfgm = grad(srf, gm)
    let dsrggm = grad(srg, gm)
    let dtrfgm = grad(trf, gm)
    let dtrggm = grad(trg, gm)

    norm2(rf - rg) :< 1e-26
    norm2(dsrfgg - dsrggg) :< 1e-26
    norm2(dtrfgg - dtrggg) :< 1e-26
    dsrfx :~ dsrgx
    dtrfx :~ dtrgx
    norm2(dsrfgm - dsrggm) :< 1e-26
    norm2(dtrfgm - dtrggm) :< 1e-26

    x.update 0.25

    norm2(rf - rg) :< 1e-26
    norm2(dsrfgg - dsrggg) :< 1e-26
    norm2(dtrfgg - dtrggg) :< 1e-26
    dsrfx :~ dsrgx
    dtrfx :~ dtrgx
    norm2(dsrfgm - dsrggm) :< 1e-26
    norm2(dtrfgm - dtrggm) :< 1e-26

  test "axexpmuly does not read packed storage before first eval":
    let rf = axexpmuly(x, gm, gg)
    let rg = exp(x * gm) * gg

    x.update 0.25
    norm2(rf - rg) :< 1e-26
