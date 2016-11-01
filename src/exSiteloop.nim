import qex
import qcdTypes

proc test =
  qexInit()

  var lat = [4,4,4,4]
  var lo = newLayout(lat)

  var v1 = lo.ColorVector()
  var v2 = lo.ColorVector()
  var v3 = lo.ColorVector()

  var s: type(simdSum(v1[0]))

  threads:
    v1 := 1
    for e in v1.all:
      v2[e] := v1[e]
      v3[e] := 2 * v1[e]
    echo v1.norm2, "\t", v2.norm2, "\t", v3.norm2

    var t: type(v1[0])
    threadBarrier()
    for e in v1.even:
      v2[e] := 2 * v1[e]
      v3[e] := v1[e]
      t += v1[e]
    var tt = simdSum(t)
    threadRankSum(tt)
    threadMaster:
      s := tt
    threadBarrier()
    echo v1.norm2, "\t", v2.norm2, "\t", v3.norm2

  echo s.toString

  qexFinalize()

test()
