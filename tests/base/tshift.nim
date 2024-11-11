import testutils
import qex, simd, simd/simdArray

# tests
#  b f = 1
#  f^n = fn
#  f^L = 1
#  bb ba fb fa = 1

proc set(x: Field, offset: seq[int]) =
  let lo = x.l
  const vl = lo.V
  let lat = lo.physGeom
  let nd = lat.len
  threads:
    for e in x:
      var t: array[vl,int]
      for i in 0..<nd:
        t = 1 + 10*t + ((offset[i] + lat[i] + lo.vcoords(i,e)) mod lat[i])
      x[e] := t

proc testf(x,y,z: auto, mu,d: int): float =
  let f = newShifter(x, mu, d)
  var offs = newSeq[int](x.l.nDim)
  offs[mu] = d
  z.set(offs)
  var res = 0.0
  threads:
    y := f ^* x
    let r = norm2(z-y)
    threadMaster:
      res = r
  if res != 0.0:
    echo "mu: ", mu, "  d: ", d, "  r: ", res
    checkeq(y, z)
  result = res

proc testfb(x,y,z: auto, mu,d: int): float =
  let f = newShifter(x, mu, d)
  let b = newShifter(x, mu, -d)
  var res = 0.0
  threads:
    y := f ^* x
    z := b ^* y
    let r = norm2(z-x)
    threadMaster:
      res = r
  if res != 0.0:
    echo "mu: ", mu, "  d: ", d, "  r: ", res
    checkeq(x, z)
    #for i in x.sites:
    #  var xi,yi,zi: float
    #  xi := x{i}
    #  yi := y{i}
    #  zi := z{i}
    #  echoAll xi, " ", yi, " ", zi
  result = res

proc test2[N,T](Smd: typedesc, lat: array[N,T]): float =
  const vl = int Smd.numNumbers
  const nd = lat.len
  var lo = newLayout(lat, vl)
  type LatReal = Field[vl, Smd]
  var x,y,z: LatReal
  x.new(lo)
  y.new(lo)
  z.new(lo)
  var offs = newSeq[int](nd)
  x.set(offs)
  for mu in 0..<nd:
    var dmax = lat[mu]
    if nRanks>1: dmax = lo.outerGeom[mu]
    for d in 1..dmax:
      result += testf(x,y,z, mu, d)
      result += testfb(x,y,z, mu, d)

template test1(Smd: typedesc, lat0: array) =
  let lat = latticeFromLocalLattice(lat0, nRanks)
  test $Smd & " lattice: " & $lat:
    let r = test2(Smd, lat)
    check(r == 0.0)

template testS1(Smd: typedesc) =
  suite $Smd:
    test1(Smd, [8])
    test1(Smd, [9,8])
    test1(Smd, [8,9])
    test1(Smd, [7,8,9])
    test1(Smd, [7,8,9,10])
    test1(Smd, [7,8,9,10,11])

template testS2(Smd: typedesc) =
  suite $Smd:
    test1(Smd, [16])
    test1(Smd, [8,8])
    test1(Smd, [9,16])
    test1(Smd, [16,9])
    test1(Smd, [8,8,8])
    test1(Smd, [8,8,8,8])
    test1(Smd, [8,8,8,8,8])

template testS4(Smd: typedesc) =
  suite $Smd:
    test1(Smd, [8,16])
    test1(Smd, [16,8])
    test1(Smd, [8,8,8])
    test1(Smd, [8,8,8,8])
    test1(Smd, [8,8,8,8,8])

template testS8(Smd: typedesc) =
  suite $Smd:
    test1(Smd, [8,8,16])
    test1(Smd, [8,16,8])
    test1(Smd, [16,8,8])
    test1(Smd, [8,8,8,8])
    test1(Smd, [8,8,8,8,8])

template testS16(Smd: typedesc) =
  suite $Smd:
    test1(Smd, [8,8,8,16])
    test1(Smd, [8,8,16,8])
    test1(Smd, [8,16,8,8])
    test1(Smd, [16,8,8,8])
    test1(Smd, [8,8,8,8,8])

qexInit()

template makeSimdArrayX(T,N,B: untyped) {.dirty.} =
  #makeSimdArray(`T X`, N, B)
  #type T = Simd[`T X`]
  #template toDoubleImpl(x: `T X`): untyped = x  # always already double
  type T = Simd[SimdArrayObj[N,B]]

#testS1(float)

makeSimdArrayX(SD1, 1, float)
testS1(SD1)
when declared(SimdD1):
  testS1(SimdD1)

makeSimdArrayX(SS1, 1, float32)
testS1(SS1)
when declared(SimdS1):
  testS1(SimdS1)

makeSimdArrayX(SD2, 2, float)
testS2(SD2)
when declared(SimdD2):
  testS2(SimdD2)

makeSimdArrayX(SS2, 2, float32)
testS2(SS2)
when declared(SimdS2):
  testS2(SimdS2)

makeSimdArrayX(SD4, 4, float)
testS4(SD4)
when declared(SimdD4):
  testS4(SimdD4)

makeSimdArrayX(SS4, 4, float32)
testS4(SS4)
when declared(SimdS4):
  testS4(SimdS4)

makeSimdArrayX(SD8, 8, float)
testS8(SD8)
when declared(SimdD8):
  testS8(SimdD8)

makeSimdArrayX(SS8, 8, float32)
testS8(SS8)
when declared(SimdS8):
  testS8(SimdS8)

makeSimdArrayX(SD16, 16, float)
testS16(SD16)
when declared(SimdD16):
  testS16(SimdD16)

makeSimdArrayX(SS16, 16, float32)
testS16(SS16)
when declared(SimdS16):
  testS16(SimdS16)

qexFinalize()
