import testutils
import qex, simd

qexInit()

template testPrec(S, D: untyped, n: static int) =
  test "S" & $n & " <-> D" & $n:
    var a: array[n, float64]
    for i in 0..<n: a[i] = float64(i+1) / 8.0
    var d {.noInit.}: D
    d := a
    var s {.noInit.}: S
    s := 0
    s := d
    for i in 0..<n:
      check s[i] == float32(a[i])
    var d2 {.noInit.}: D
    d2 := 0
    d2 := s
    for i in 0..<n:
      check d2[i] == a[i]
    let s3 = toSingle(d)
    let d3 = toDouble(s)
    for i in 0..<n:
      check s3[i] == float32(a[i])
      check d3[i] == a[i]

suite "mixed precision SIMD assignment":
  when declared(SimdS1) and declared(SimdD1): testPrec(SimdS1, SimdD1, 1)
  when declared(SimdS2) and declared(SimdD2): testPrec(SimdS2, SimdD2, 2)
  when declared(SimdS4) and declared(SimdD4): testPrec(SimdS4, SimdD4, 4)
  when declared(SimdS8) and declared(SimdD8): testPrec(SimdS8, SimdD8, 8)
  when declared(SimdS16) and declared(SimdD16): testPrec(SimdS16, SimdD16, 16)

let lat = latticeFromLocalLattice([8,8,8,8], nRanks)
var
  l = newLayout(lat)
  g = l.newGauge
  r = newRNGField(RngMilc6, l, 987654321)
threads: g.random r

suite "mixed precision gauge fields":
  test "double to single to double":
    let gd = g.newGaugeS.newGauge
    let gd2 = gd.newGaugeS.newGauge
    var d = l.newGauge
    for mu in 0..<g.len:
      d[mu] := gd[mu]
      d[mu] -= g[mu]
      check(d[mu].norm2 <= 1e-14 * g[mu].norm2)  # single precision rounding only
      d[mu] := gd2[mu]
      d[mu] -= gd[mu]
      check(d[mu].norm2 == 0)  # a second round trip changes nothing

qexFinalize()
