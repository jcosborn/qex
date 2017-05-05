import os
import base
export base
import layout
export layout
#import shifts
#export shifts
import field
export field
import io
export io

when isMainModule:
  qexInit()
  echo "rank ", myRank, "/", nRanks
  var lat = [4,4,4,4]
  var lo = newLayout(lat)
  lo.makeShift(0,1)
  lo.makeShift(3,-2,"even")
  qexFinalize()
