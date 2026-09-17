import qex, gauge/plaquette
import testutils

proc run(lo: auto, nc: static int) =
  let
    g = lo.newGauge(nc)
    gp = lo.newGauge(nc)
    gm = lo.newGauge(nc)
    f = lo.newGauge(nc)
    e = lo.newGauge(nc)
    fp = lo.newGauge(nc)
    fm = lo.newGauge(nc)
    ds = @[lo.newGauge(nc), lo.newGauge(nc), lo.newGauge(nc)]
    eps = 1.0e-4
  var r = lo.newRNGField(Philox4x64, 91820126'u64)
  threads:
    g.random r
    for k in 0..<ds.len:
      ds[k].gaussian r
      for mu in 0..<g.len:
        ds[k][mu] *= 0.1
    for mu in 0..<g.len:
      g[mu] += ds[0][mu]
  let
    wa = newPlaqWork(g[0], action=true)
    ws = newPlaqWork(g[0])
  suite "Plaquette stencil Nc=" & $nc:
    test "plaquette and staple sums match production kernels":
      let
        a = wa.plaqSum(g)
        want = -float(nc)*GaugeActionCoeffs(plaq: 1.0).gaugeAction1(g)
      check abs(a-want) < 1.0e-11*(1.0+abs(want))
      let empty: seq[type(g)] = @[]
      ws.stapleSum(g, empty, f)
      GaugeActionCoeffs(plaq: float(nc)).gaugeActionDeriv(g, e)
      check relativeDiff(f, e) < 2.0e-13

    test "directional jets and seed symmetry through cubic degree":
      for order in 1..3:
        let
          w = newPlaqWork(g[0], order)
          prev = newPlaqWork(g[0], order-1)
          seeds = ds[0..<order]
          before = ds[0..<order-1]
        w.stapleSum(g, seeds, f)
        threads:
          for mu in 0..<g.len:
            gp[mu] := g[mu]+eps*ds[order-1][mu]
            gm[mu] := g[mu]-eps*ds[order-1][mu]
        prev.stapleSum(gp, before, fp)
        prev.stapleSum(gm, before, fm)
        threads:
          for mu in 0..<g.len:
            e[mu] := (0.5/eps)*(fp[mu]-fm[mu])
        check relativeDiff(f, e) < 2.0e-9
        var rev = seeds
        swap(rev[0], rev[^1])
        w.stapleSum(g, rev, e)
        check relativeDiff(f, e) < 2.0e-13
        if order == 1:
          threads:
            for mu in 0..<g.len:
              e[mu] := 0
          GaugeActionCoeffs(plaq: -float(nc)).gaugeDerivDeriv2(g, ds[0], e)
          check relativeDiff(f, e) < 2.0e-13
        if order == 3:
          w.stapleSum(gp, seeds, e)
          check relativeDiff(f, e) == 0

    test "repeated seeds, workspace clones, and zero fourth jet":
      let
        w = newPlaqWork(g[0], 2)
        clone = w.newOneOf
        same = @[ds[0], ds[0]]
      w.stapleSum(g, same, f)
      clone.stapleSum(gp, same, fp)
      w.stapleSum(g, same, e)
      check relativeDiff(f, e) == 0
      threads:
        for mu in 0..<g.len:
          gp[mu] := g[mu]+eps*ds[0][mu]
          gm[mu] := g[mu]-eps*ds[0][mu]
      stapleSum(gp, @[ds[0]], fp)
      stapleSum(gm, @[ds[0]], fm)
      threads:
        for mu in 0..<g.len:
          e[mu] := (0.5/eps)*(fp[mu]-fm[mu])
      check relativeDiff(f, e) < 2.0e-9
      stapleSum(g, @[ds[0],ds[1],ds[2],ds[0]], f)
      threads:
        for mu in 0..<g.len:
          e[mu] := 0
      check relativeDiff(f, e) == 0

    test "workspaces reject incompatible bundles and aliased outputs":
      var lat = newSeq[int](lo.nDim)
      for mu, n in lo.physGeom: lat[mu] = int(n)
      lat[0] *= 2
      let wrong = lat.newLayout.newGauge(nc)
      expect(ValueError): discard wa.plaqSum(wrong)
      expect(ValueError): stapleSum(g, wrong)
      expect(ValueError): stapleSum(g, @[wrong], f)
      threads:
        for mu in 0..<g.len: fp[mu] := g[mu]
      expect(ValueError): stapleSum(g, g)
      for order in 1..3:
        let w = newPlaqWork(g[0], order)
        expect(ValueError): w.stapleSum(g, ds[0..<order], ds[0])
      check relativeDiff(g, fp) == 0.0
      let w = newPlaqWork(g[0], 3)
      w.stapleSum(g, ds, e)
      w.stapleSum(gp, ds, gp)  # The third jet reads seeds only.
      check relativeDiff(gp, e) == 0.0

qexInit()
letParam:
  lat = latticeFromLocalLattice(@[4,4,4,4], nRanks)
let lo = lat.newLayout
run(lo, 3)
run(lo, 1)

qexFinalize()
