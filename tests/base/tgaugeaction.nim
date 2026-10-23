import qex
import algorithms/numdiff
import testutils

qexInit()

letParam:
  lat = latticeFromLocalLattice(@[4,4,4,4], nRanks)
let
  (lo, g, r0) = setupLattice(lat)
  h = lo.newGauge
  base = lo.newGauge
var r = r0
threads:
  h.randomTAH r
  base.random r

suite "Gauge action subset kernels":
  test "Wilson fast path and rectangle fallback match gaugeAction2":
    let
      cw = GaugeActionCoeffs(plaq: 5.4)
      cr = Symanzik(5.4)
    check(cw.gaugeAction1(g) ~ cw.gaugeAction2(g))
    check(cr.gaugeAction1(g) ~ cr.gaugeAction2(g))

  test "subset derivative overwrites the masked output":
    let c = GaugeActionCoeffs(plaq: 5.4)
    var full = lo.newGauge
    threads:
      for mu in 0..<g.len:
        full[mu] := 0
    c.gaugeDeriv2(g, full)

    for p in 0..1:
      let sub = lo.getSubset(if p == 0: "even" else: "odd")
      for d in 0..<g.len:
        var
          got = lo.newGauge
          expected = lo.newGauge
        threads:
          for mu in 0..<g.len:
            got[mu] := base[mu]
            expected[mu] := base[mu]
          threadBarrier()
          expected[d] := 0
          threadBarrier()
          for x in sub:
            expected[d][x] := full[d][x]
        c.gaugeDeriv2Subset(g, got, p, d)
        check(got ~ expected)

  test "subset Hessian matches a masked full Hessian":
    let c = GaugeActionCoeffs(plaq: 5.4)
    for p in 0..1:
      let sub = lo.getSubset(if p == 0: "even" else: "odd")
      for d in 0..<g.len:
        var
          hm = lo.newGauge
          got = lo.newGauge
          expected = lo.newGauge
        threads:
          for mu in 0..<g.len:
            hm[mu] := 0
            got[mu] := 0
            expected[mu] := 0
          threadBarrier()
          for x in sub:
            hm[d][x] := h[d][x]
        c.gaugeDerivDeriv2Subset(g, h, got, p, d)
        c.gaugeDerivDeriv2(g, hm, expected)
        check(got ~ expected)

  test "subset derivative and Hessian match finite differences":
    let
      c = GaugeActionCoeffs(plaq: 5.4)
      p = 1
      d = min(2, g.len - 1)
      sub = lo.getSubset("odd")
    var
      gt = lo.newGauge
      ds = lo.newGauge
      hs = lo.newGauge
    threads:
      for mu in 0..<g.len:
        ds[mu] := 0
        hs[mu] := 0
    c.gaugeDeriv2Subset(g, ds, p, d)

    proc actionAt(t: float): float =
      threads:
        for mu in 0..<g.len:
          gt[mu] := g[mu]
        threadBarrier()
        for x in sub:
          gt[d][x] := g[d][x] + t*h[d][x]
      c.gaugeAction1(gt)

    var num, err: float
    ndiff(num, err, actionAt, 0.0, 0.1, ordMax = 4)
    let ana = redot(h[d], ds[d])
    check abs(num-ana) <= max(1.0e-8, 32.0*err)

    c.gaugeDerivDeriv2Subset(g, h, hs, p, d)
    proc derivativeAt(t: float): float =
      threads:
        for mu in 0..<g.len:
          gt[mu] := g[mu] + t*base[mu]
      c.gaugeDeriv2Subset(gt, ds, p, d)
      redot(h[d], ds[d])

    ndiff(num, err, derivativeAt, 0.0, 0.1, ordMax = 4)
    var hana = 0.0
    for mu in 0..<g.len:
      hana += redot(base[mu], hs[mu])
    check abs(num-hana) <= max(1.0e-8, 32.0*err)

  test "subset Hessian accumulation variants":
    let
      c = GaugeActionCoeffs(plaq: 5.4)
      p = 1
      d = min(2, g.len - 1)
      sub = lo.getSubset("odd")
    var delta = lo.newGauge
    threads:
      for mu in 0..<g.len:
        delta[mu] := 0
    c.gaugeDerivDeriv2Subset(g, h, delta, p, d)

    block:
      var
        got = lo.newGauge
        expected = lo.newGauge
      threads:
        for mu in 0..<g.len:
          got[mu] := base[mu]
          expected[mu] := base[mu] + delta[mu]
      c.gaugeDerivDeriv2SubsetAdd(g, h[d], got, p, d)
      check(got ~ expected)

    block:
      var
        got = lo.newGauge
        expected = lo.newGauge
      threads:
        for mu in 0..<g.len:
          got[mu] := 0
          expected[mu] := base[mu]
        threadBarrier()
        for x in sub:
          got[d][x] := h[d][x]
          expected[d][x] := h[d][x]
        threadBarrier()
        for mu in 0..<g.len:
          expected[mu] += delta[mu]
      c.gaugeDerivDeriv2SubsetAddBase(g, h[d], base, got, p, d)
      check(got ~ expected)
      # f = base, as in the stout gradient over a donated upstream:
      # f[S] += H[S] and f[~S] = base[~S] + H[~S] give f = base + H.
      var inplace = lo.newGauge
      threads:
        for mu in 0..<g.len:
          inplace[mu] := base[mu]
          expected[mu] := base[mu] + delta[mu]
      c.gaugeDerivDeriv2SubsetAddBase(g, h[d], inplace, inplace, p, d)
      check(inplace ~ expected)

    block:
      let hs = @[h[d], base[d], h[d]]
      var
        hm = lo.newGauge
        got = lo.newGauge
        expected = lo.newGauge
        w = g[d].newOneOf
      threads:
        for mu in 0..<g.len:
          hm[mu] := 0
          got[mu] := 0
          expected[mu] := 0
        threadBarrier()
        for x in sub:
          hm[d][x] := h[d][x] + base[d][x] + h[d][x]
      c.gaugeDerivDeriv2SubsetSum(g, hs, w, got, p, d)
      c.gaugeDerivDeriv2Subset(g, hm, expected, p, d)
      check(got ~ expected)

  test "subset sums preserve operand order and scratch outside the parity":
    let d = min(1, g.len-1)
    let w = g[d].newOneOf
    let want = g[d].newOneOf
    let got = lo.newGauge
    let expected = lo.newGauge
    let work = newLoopWork(g[0])
    for c in [GaugeActionCoeffs(), GaugeActionCoeffs(plaq: 0.8),
              GaugeActionCoeffs(plaq: 0.8, rect: -0.11),
              GaugeActionCoeffs(plaq: 0.8, rect: -0.11, pgm: 0.07)]:
      for parity in 0..1:
        let sub = lo.getSubset(if parity == 0: "even" else: "odd")
        let sw = newSubsetWork(g[0], parity, d)
        for seeds in [@[h[d]], @[h[d], base[d], g[d]]]:
          threads:
            w := base[d]
            want := base[d]
            for mu in 0..<g.len: got[mu] := base[mu]
            threadBarrier()
            for x in sub:
              want[x] := seeds[0][x]
              for k in 1..<seeds.len: want[x] += seeds[k][x]
          c.gaugeDerivDeriv2SubsetSum(g, seeds, w, got, parity, d, sw, (23,0), work=work)
          c.gaugeDerivDeriv2Subset(g, want, expected, parity, d, sw, (23,0), work=work)
          check norm2(w-want) == 0
          # Two instantiations of one kernel, which compilers may round differently.
          for mu in 0..<g.len:
            check norm2(got[mu]-expected[mu]) <= (4*epsilon(float))^2*(1+norm2(expected[mu]))

  test "subset works sharing their call buffers keep their link boundaries":
    # Derivatives at A and B, then the Hessians at B and A under the tags of
    # their derivatives: A's reads the link boundaries it kept while B's calls
    # ran through the shared buffers, B's Hessian through A's derivative ones.
    let
      c = GaugeActionCoeffs(plaq: 5.4)
      sh = newSubsetShifts(g[0])
      d1 = min(1, g.len - 1)
      ta: LoopTag = (1, 1)
      tb: LoopTag = (2, 1)
    for (a, b) in [((0, 0), (1, d1)), ((1, d1), (1, 0))]:
      let
        wa = newSubsetWork(g[0], a[0], a[1], sh)
        wb = newSubsetWork(g[0], b[0], b[1], sh)
      var
        fa = g[0].newOneOf
        fb = g[0].newOneOf
        ga = lo.newGauge
        gb = lo.newGauge
        ea = lo.newGauge
        eb = lo.newGauge
      c.gaugeDeriv2SubsetWork(g, fa, wa, a[0], a[1], true, ta)
      c.gaugeDeriv2SubsetWork(g, fb, wb, b[0], b[1], true, tb)
      c.gaugeDerivDeriv2Subset(g, h, gb, b[0], b[1], wb, tb)
      c.gaugeDerivDeriv2Subset(g, h, ga, a[0], a[1], wa, ta)
      c.gaugeDerivDeriv2Subset(g, h, ea, a[0], a[1])
      c.gaugeDerivDeriv2Subset(g, h, eb, b[0], b[1])
      check(ga ~ ea)
      check(gb ~ eb)

  test "subset kernels reject adjoint coefficients":
    let c = GaugeActionCoeffs(adjplaq: 1.0)
    var
      outg = lo.newGauge
      w = g[0].newOneOf
    expect(ValueError):
      c.gaugeDeriv2Subset(g, outg, 0, 0)
    expect(ValueError):
      c.gaugeDerivDeriv2Subset(g, h, outg, 0, 0)
    expect(ValueError):
      c.gaugeDerivDeriv2SubsetAdd(g, h[0], outg, 0, 0)
    expect(ValueError):
      c.gaugeDerivDeriv2SubsetAddBase(g, h[0], base, outg, 0, 0)
    expect(ValueError):
      c.gaugeDerivDeriv2SubsetSum(g, @[h[0]], w, outg, 0, 0)

qexFinalize()
