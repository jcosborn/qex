import qex
import physics/qcdTypes

proc test =
  qexInit()

  var lat = [4,4,4,4]
  var lo = newLayout(lat)

  var v1 = lo.ColorVector()
  var v2 = lo.ColorVector()
  var v3 = lo.ColorVector()

  # single-site ColorVector for sum of v1 over sites
  var s: type(simdSum(v1[0]))

  # start thread block ("omp parallel")
  threads:
    # initialize some fields
    v1 := 1
    v3 := 1

    # loop over all (vectorized) site indices in v1
    for e in v1.all:
      v2[e] := 2 * v1[e]
      v3[e] += v2[e]

    echo v1.norm2, "\t", v2.norm2, "\t", v3.norm2

    # thread local temporary for summing v1
    # Nim automatically initializes to binary '0', unless told not to
    var t: type(v1[0])

    # loop over even (vectorized) site indices in v1
    # barrier is for safety since we are switching subsets
    # the barrier may be automatically inserted in the future
    threadBarrier()
    for e in v1.even:
      v2[e] := 2 * v1[e]
      v3[e] += v2[e]
      t += v1[e]
    threadBarrier()

    # sum over SIMD vector lanes
    var tt = simdSum(t)
    # sum over threads and ranks, result overwrites tt
    threadRankSum(tt)
    # only master thread assigns result to s
    threadMaster:
      s := tt

    echo v1.norm2, "\t", v2.norm2, "\t", v3.norm2

  echo s.toString

  qexFinalize()

test()
