suite "gauge fused":
  setup:
    let gg {.used.} = grt.toGvalue(g)
    let gu {.used.} = grt.toGvalue(u)
    let gp {.used.} = grt.toGvalue(p)
    let gq {.used.} = grt.toGvalue(q)
    let gm {.used.} = grt.toGvalue(m)
    let x {.used.} = grt.toGvalue(a)
    let y {.used.} = grt.toGvalue(b)

  test "adjmul":
    ckbinarynorm2grad(gg.adjmul gu, gg.adj * gu, gg, gu, 1e-25)
    ckgradm2(adjmul, gg, gu, gp, gq, gm)

  test "muladj":
    ckbinarynorm2grad(gg.muladj gu, gg * gu.adj, gg, gu, 1e-25)
    ckgradm2(muladj, gg, gu, gp, gq, gm)

  test "contractProjTAH":
    ckbinarynorm2grad(contractProjTAH(gg, gu), projTAH(gg * gu.adj), gg, gu, 1e-26)
    ckgradm2(contractProjTAH, gg, gu, gp, gq, gm)

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

  test "axexp":
    let rf = axexp(x, gm)
    let rg = exp(x*gm)
    norm2(rf - rg) :< 1e-26
    let srf = retr(rf * gu)
    let srg = retr(rg * gu)
    grad(srf, x) :~ grad(srg, x)
    norm2(grad(srf, gm) - grad(srg, gm)) :< 1e-26
    ckgradm2(axexp, x, gm, y, 0.05*gq, gp)

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

  test "axexpmuly pack shares helper and result without breaking derivatives":
    let packed = axexpmulyPack(x, gm, gg)
    let packedExp = packed[0]
    let packedResult = packed[1]
    let refExp = axexp(x, gm)
    let refResult = exp(x * gm) * gg
    let packedScalar = retr(packedExp * gu) + retr(packedResult * gp)
    let refScalar = retr(refExp * gu) + retr(refResult * gp)

    norm2(packedExp - refExp) :< 1e-26
    norm2(packedResult - refResult) :< 1e-26
    norm2(packedResult - axexpmuly(x, gm, gg)) :< 1e-26
    packedScalar :~ refScalar
    grad(packedScalar, x) :~ grad(refScalar, x)
    norm2(grad(packedScalar, gm) - grad(refScalar, gm)) :< 1e-25
    norm2(grad(packedScalar, gg) - grad(refScalar, gg)) :< 1e-26

    x.update 0.25

    norm2(packedExp - refExp) :< 1e-26
    norm2(packedResult - refResult) :< 1e-26
    norm2(packedResult - axexpmuly(x, gm, gg)) :< 1e-26
    packedScalar :~ refScalar
    grad(packedScalar, x) :~ grad(refScalar, x)
    norm2(grad(packedScalar, gm) - grad(refScalar, gm)) :< 1e-25
    norm2(grad(packedScalar, gg) - grad(refScalar, gg)) :< 1e-26

  test "axexpmuly":
    let rf = axexpmuly(x, gm, gg)
    let rg = exp(x*gm)*gg
    norm2(rf - rg) :< 1e-26
    let srf = retr(rf * gu)
    let srg = retr(rg * gu)
    grad(srf, x) :~ grad(srg, x)
    norm2(grad(srf, gm) - grad(srg, gm)) :< 1e-26
    norm2(grad(srf, gg) - grad(srg, gg)) :< 1e-26
    ckgradm3(axexpmuly, x, gm, gu, y, 0.05*gq, gg, gp)

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
