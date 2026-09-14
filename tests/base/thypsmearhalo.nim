import qex, gauge/[hypsmear, hypsmear2], comms/halo
import std/importutils
import testutils

qexInit()
letParam:
  expectRanks = nRanks
check nRanks == expectRanks
letParam:
  lat = latticeFromLocalLattice(@[4,4,4,4], nRanks)
let lo = lat.newLayout
var rng = lo.newRNGField(MRG32k3a, 987654321'u)
let
  g = lo.newGauge
  chain = lo.newGauge
  reference = lo.newGauge
  candidate = lo.newGauge
  fReference = lo.newGauge
  fCandidate = lo.newGauge
  delta = lo.newGauge
  coef = HypCoefs(alpha1: 0.4, alpha2: 0.5, alpha3: 0.5)
var info: PerfInfo
threads:
  g.gaussian rng
  chain.gaussian rng

proc finiteBits(bits: uint64): bool =
  (bits and 0x7ff0000000000000'u64) != 0x7ff0000000000000'u64

proc compare(x, y: auto, tolerance: float) =
  var invalid = 0
  for mu in 0..<x.len:
    # Check stored components before arithmetic: fast-math can assume a
    # floating-point function's argument/result is finite and erase its guard.
    for e in 0..<lo.nSitesOuter:
      for i in 0..<x[mu][0].nrows:
        for j in 0..<x[mu][0].ncols:
          for k in 0..<lo.V:
            if not finiteBits(cast[uint64](x[mu][e][i,j].re[k])): inc invalid
            if not finiteBits(cast[uint64](x[mu][e][i,j].im[k])): inc invalid
    threads:
      delta[mu] := x[mu] - y[mu]
    let error = delta[mu].norm2 / max(1.0, y[mu].norm2)
    check error < tolerance
  rankSum(invalid)
  check invalid == 0

suite "HYP halo derivative domain":
  test "undefined halo intermediates never contribute to physical forces":
    let pullback = coef.smearGetForce(g, reference, info)
    pullback(fReference, chain)
    let ht = newHypTemps(g)
    privateAccess(ht.type)
    let undefined = cast[float](0x7ff8000000000001'u64)
    for repeat in 0..1:
      # Forward only constructs a subset of these padded intermediates.
      # Poisoning the rest pins the reverse pass's dependency domain.
      for mu in 0..<4:
        for nu in 0..<4:
          if mu == nu: continue
          ht.h1x[mu][nu].halo := undefined
          ht.h1[mu][nu].halo := undefined
          ht.h2x[mu][nu].halo := undefined
          ht.h2[mu][nu].halo := undefined
      ht.smear(coef, candidate)
      ht.force(coef, fCandidate, chain)
      compare(candidate, reference, 1e-18)
      compare(fCandidate, fReference, 1e-16)

qexFinalize()
