import unittest
import qex
import simd/simdArray

proc test2(Smd: typedesc, lat: array): float =
  const vl = Smd.numNumbers
  let nd = lat.len
  var lo = newLayout(lat, vl)
  type LatReal = Field[vl, Smd]
  var x,y,z: LatReal
  x.new(lo)
  y.new(lo)
  z.new(lo)
  let f = newShifters(x, 1)
  let b = newShifters(x, -1)
  for e in x:
    var t: array[vl,int]
    for i in 0..<nd:
      t = 10*t + lo.vcoords(i,e)
    x[e] := t
  for mu in 0..<nd:
    y := f[mu] ^* x
    z := b[mu] ^* y
    let r = norm2(z-x)
    if r != 0.0:
      echo "mu: ", mu, "  r: ", r
    threadMaster:
      result += r

template test1(Smd: typedesc, lat: array) =
  test "lattice: " & $lat:
    let r = test2(Smd, lat)
    check(r == 0.0)

qexInit()

makeSimdArray(SimdD1, 1, float64)
template isWrapper*(x: SimdD1): untyped = false
template toDoubleImpl*(x: SimdD1): untyped = x
suite "SimdD1":
  test1(SimdD1, [8])
  test1(SimdD1, [8,8])
  test1(SimdD1, [8,8,8])
  test1(SimdD1, [8,8,8,8])
  test1(SimdD1, [8,8,8,8,8])

makeSimdArray(SimdD2, 2, float64)
template isWrapper*(x: SimdD2): untyped = false
template toDoubleImpl*(x: SimdD2): untyped = x
suite "SimdD2":
  test1(SimdD2, [16])
  test1(SimdD2, [8,8])
  test1(SimdD2, [8,8,8])
  test1(SimdD2, [8,8,8,8])
  test1(SimdD2, [8,8,8,8,8])

#makeSimdArray(SimdD4, 4, float64)
#template isWrapper*(x: SimdD4): untyped = false
#template toDoubleImpl*(x: SimdD4): untyped = x
#test(SimdD4, [8,8])

#makeSimdArray(SimdD8, 8, float64)
#template isWrapper*(x: SimdD8): untyped = false
#template toDoubleImpl*(x: SimdD8): untyped = x
#test(SimdD8, [8,8,8])
#test(SimdD8, [8,8,8,8])

qexFinalize()
