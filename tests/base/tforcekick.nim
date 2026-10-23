import qex
import testutils

qexInit()
# Four outer sites at VLEN=4: six threads include empty partitions.
let small = newLayout(latticeFromLocalLattice([4,4], nRanks), 4)
let full = newLayout(latticeFromLocalLattice([4,4,4,4], nRanks))

proc run(lo: auto, nc: static int) =
  let g = lo.newGauge(nc)
  let p = lo.newGauge(nc)
  let f = lo.newGauge(nc)
  let q = lo.newGauge(nc)
  let want = lo.newGauge(nc)
  let work = newLoopWork(g[0])
  var rng = lo.newRNGField(Philox4x64, 20261023'u64)
  threads:
    g.random rng
    p.randomTAH rng
  let cases = [GaugeActionCoeffs(), GaugeActionCoeffs(plaq: 5.4),
    Symanzik(5.4), GaugeActionCoeffs(plaq: 0.8, rect: -0.11, pgm: 0.07)]
  suite "Force kick Nc=" & $nc & " lattice=" & $lo.physGeom:
    test "momentum and all statistics match the stored force":
      for c in cases:
        c.force(g, f)
        var sum = 0.0
        var low = norm2(f[0]{0})
        var high = low
        for mu in 0..<f.len:
          for i in 0..<lo.nSites:
            let n = norm2(f[mu]{i})
            sum += n
            low = min(low, n)
            high = max(high, n)
        rankSum(sum)
        rankMax(high)
        var neg = -low
        rankMax(neg)
        for val in [0.0, -0.125]:
          let a: float = val
          threads:
            for mu in 0..<f.len: want[mu] := a*f[mu]+p[mu]
          let got = c.forceKick(g, p, q, a, work=work, tag=(1,1))
          check relativeDiff(q,want) < 2e-12
          withCT(2e-12):
            check got.s2 ~ sum
            check got.n2min ~ -neg
            check got.n2max ~ high

    test "momentum may own the output and links remain disjoint":
      for c in cases:
        let expected = c.forceKick(g, p, want, -0.125, work=work, tag=(1,1))
        threads:
          for mu in 0..<p.len: q[mu] := p[mu]
        let got = c.forceKick(g, q, q, -0.125, work=work, tag=(1,1))
        check relativeDiff(q,want) == 0
        check got == expected
        expect ValueError: discard c.forceKick(g, p, g, -0.125, work=work)

run(small,1)
run(small,3)
run(full,1)
run(full,3)
qexFinalize()
