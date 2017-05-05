#import base/globals
#import base/threading
#export threading
#import comms/comms
#export comms
import base/qexInternal
export qexInternal
import base/stdUtils
export stdUtils
import base/metaUtils
export metaUtils
import base/basicOps
export basicOps
import base/alignedMem
export alignedMem
#import base/profile
#export profile

when isMainModule:
  qexInit()
  echo "rank ", myRank, "/", nRanks
  #var lat = [4,4,4,4]
  #var lo = newLayout(lat)
  #lo.makeShift(0,1)
  #lo.makeShift(3,-2,"even")
  qexFinalize()
