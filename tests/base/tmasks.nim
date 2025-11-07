import qex
import testutils

qexInit()

# Accept command line arguments and set up parameters
let
  #source_type = strParam("source_type", "Z4") # Z4, Z2, U1, Gauss
  dilute_eo = parseDilution "EO"
  dilute_corner = parseDilution "CORNER"

let
  lat = @[8,8,8,8]
  nt = lat[^1]
var
  lo = lat.newLayout
  g = lo.newGauge
  r = newRNGField(RngMilc6, lo)

threads: g.random r

var
  eta = lo.ColorVector        # Noise source
  tmps = lo.ColorVector

suite "test random source and masked assignment":
  test "random z4 source norm2":
    var n = 0.0
    threads:
      eta.z4 r
      var n2 = eta.norm2
      threadMaster:
        n = n2
    check(n ~ float(3*lat[0]*lat[1]*lat[2]*lat[3]))

  test "dilution EO and masked assignment":
    for t in 0..<nt:
      for dl in dilution(dilute_eo):
        var n = 0.0
        threads:
          tmps := 0
          threadBarrier()
          for i in tmps.sites(dl):
            if lo.coords[^1][i] == t:
              tmps{i} := eta{i}
          threadBarrier()
          var n2 = tmps.norm2
          threadMaster:
            n = n2
        check(n ~ float(3*lat[0]*lat[1]*lat[2])/2)

  test "dilution CORNER and masked assignment":
    for t in 0..<nt:
      for dl in dilution(dilute_corner):
        var n = 0.0
        threads:
          tmps := 0
          threadBarrier()
          for i in tmps.sites(dl):
            if lo.coords[^1][i] == t:
              tmps{i} := eta{i}
          threadBarrier()
          var n2 = tmps.norm2
          threadMaster:
            n = n2
        check(n ~ float(3*lat[0]*lat[1]*lat[2])/8)

qexFinalize()
