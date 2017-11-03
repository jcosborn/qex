import qex
import qex/physics/qcdTypes

qexInit()

var lat = [16,16,16,16]
var lo = newLayout(lat)
var v = lo.ColorVector()
var v1x = sqrt(3*lo.physVol.float)
var n = 1000

#threads:
#  v := 1
#  let v2 = v.norm2
#  let v1 = sqrt(v2)
#  threadMaster:
#    v1x = v1

threads:
  for i in 1..n:
    v := i
    let v2 = v.norm2
    let v1 = sqrt(v2)
    let v1g = i.float*v1x
    if abs(v1-v1g)>1e-10:
      echo "error: ", v1, " : ", v1g

qexFinalize()
