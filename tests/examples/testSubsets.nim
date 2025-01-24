import qex
import physics/qcdTypes

proc test =
  qexInit()

  var lat = [8,8,8,8]
  var lo = newLayout(lat)
  let vnc = lo.physVol * getDefaultNc()

  var v1 = lo.ColorVector()
  var v2 = lo.ColorVector()
  var v3 = lo.ColorVector()
  var v1e = v1.even

  # start thread block ("omp parallel")
  threads:
    v1 := 1
    threadBarrier()  # required when changing subsets if potential race condition
    v1.even := 2
    threadBarrier()

    # loop over all (vectorized) site indices in v1
    for e in v1.all:
      v2[e] := 2 * v1[e]
      v3[e] += v2[e]

    echo v1.norm2, "\t", v2.norm2, "\t", v3.norm2
    doAssert v1.norm2 ==  2.5*vnc   # even: 2, odd: 1
    doAssert v2.norm2 == 10.0*vnc   # even: 4, odd: 2
    doAssert v3.norm2 == 10.0*vnc   # even: 4, odd: 2

    threadBarrier()
    v1e := 3
    threadBarrier()
    for e in v1.even:
      v2[e] := 2 * v1[e]
      v3[e] += v2[e]
    threadBarrier()

    echo v1.norm2, "\t", v2.norm2, "\t", v3.norm2
    doAssert v1.norm2 ==  5.0*vnc   # even:  3, odd: 1
    doAssert v2.norm2 == 20.0*vnc   # even:  6, odd: 2
    doAssert v3.norm2 == 52.0*vnc   # even: 10, odd: 2

  qexFinalize()

test()
