import os
import globals
export globals

import threading
export threading
import comms
export comms
import layout
export layout
import shifts
export shifts
import field
export field
import reader
export reader
import qexInternal
export qexInternal

when isMainModule:
  qexInit()
  echo "rank ", myRank, "/", nRanks
  var lat = [4,4,4,4]
  var lo = newLayout(lat)
  lo.makeShift(0,1)
  lo.makeShift(3,-2,"even")
  qexFinalize()
