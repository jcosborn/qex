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
  for e in x:
    var t: array[vl,int]
    for i in 0..<nd:
      t = 1 + 10*t + lo.vcoords(i,e)
    x[e] := t
  var res = 0.0
  for mu in 0..<nd:
    let dmax = lat[mu] - 1
    for d in 1..<dmax:
      let f = newShifter(x, mu, d)
      let b = newShifter(x, mu, -d)
      threads:
        y := f ^* x
        z := b ^* y
        let r = norm2(z-x)
        if r != 0.0:
          echo "mu: ", mu, "  d: ", d, "  r: ", r
          for i in x.sites:
            var xi,yi,zi: float
            xi := x{i}
            yi := y{i}
            zi := z{i}
            echoAll xi, " ", yi, " ", zi
        threadMaster:
          res += r
  result = res

template test1(Smd: typedesc, lat: array) =
  test "lattice: " & $lat:
    let r = test2(Smd, lat)
    check(r == 0.0)

qexInit()

template makeSimdArrayX(T,N,B: untyped) {.dirty.} =
  makeSimdArray(`T X`, N, B)
  type T = Simd[`T X`]

makeSimdArrayX(SimdD1, 1, float64)
suite "SimdD1":
  test1(SimdD1, [8])
  test1(SimdD1, [9,8])
  test1(SimdD1, [8,9])
  test1(SimdD1, [7,8,9])
  test1(SimdD1, [7,8,9,10])
  test1(SimdD1, [7,8,9,10,11])

makeSimdArrayX(SimdD2, 2, float64)
makeSimdArrayX(SimdD4, 4, float64)
makeSimdArrayX(SimdD8, 8, float64)
if nRanks == 1:
  suite "SimdD2":
    test1(SimdD2, [16])
    test1(SimdD2, [8,8])
    test1(SimdD2, [8,8,8])
    test1(SimdD2, [8,8,8,8])
    test1(SimdD2, [8,8,8,8,8])

  #makeSimdArray(SimdD4, 4, float64)
  #template isWrapper*(x: SimdD4): untyped = false
  #template toDoubleImpl*(x: SimdD4): untyped = x
  suite "SimdD4":
    test1(SimdD4, [8,8])

  #makeSimdArray(SimdD8, 8, float64)
  #template isWrapper*(x: SimdD8): untyped = false
  #template toDoubleImpl*(x: SimdD8): untyped = x
  suite "SimdD8":
    test1(SimdD8, [8,8,8])
    test1(SimdD8, [8,8,8,8])

qexFinalize()
