import qex
import testutils

qexInit()

# Accept command line arguments and set up parameters
let
  source_type = strParam("source_type", "Z4") # Z4, Z2, U1, Gauss
  dilute_type = strParam("dilute_type", "EO").parseDilution # EO, CORNER

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
    threads:
      eta.z4 r
      var n = eta.norm2
      threadSingle:
        check(n ~ float(3*lat[0]*lat[1]*lat[2]*lat[3]))

  test "dilution and masked assignment":
    for t in 0..<nt:
      for dl in dilution(dilute_type):
        threads:
          tmps := 0
          threadBarrier()
          for i in tmps.sites(dl):
            if lo.coords[^1][i] == t:
              # tmps{i} := eta{i}    # Doesn't work
              forO c, 0, tmps{0}.len-1:
                tmps{i}[c].re := eta{i}[c].re
                tmps{i}[c].im := eta{i}[c].im
          threadBarrier()
          var n = tmps.norm2
          threadSingle:
            check(n ~ float(3*lat[0]*lat[1]*lat[2])/2)

qexFinalize()
