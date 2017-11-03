import os
import qex/base
export base
import qex/layout
export layout
#import qex/shifts
#export shifts
import qex/field
export field
import qex/io
export io
import qex/gauge
export gauge
import qex/physics/qcdTypes, qex/physics/stagD, qex/physics/hisqLinks
export qcdTypes, stagD, hisqLinks
import qex/rng
export rng
import qex/eigens
export eigens
import qex/algorithms/dilution
export dilution

when isMainModule:
  qexInit()
  echo "rank ", myRank, "/", nRanks
  var lat = [4,4,4,4]
  var lo = newLayout(lat)
  lo.makeShift(0,1)
  lo.makeShift(3,-2,"even")
  qexFinalize()
