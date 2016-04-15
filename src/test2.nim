import qex
import qcdTypes
import gaugeUtils
import os

proc test() =
  #let defaultGaugeFile = "l88.scidac"
  let defaultLat = @[8,8,8,8]
  defaultSetup()
  threads:
    for i in 0..<4:
      #g[i] := 1
      #var tr:type(g[0][0][0,0])
      #for x in g[i].all:
      #  for c in 0..<nc:
      #    #echo g[i][x][c,c]
      #    tr += g[i][x][c,c]
      #echo i, ": ", vecSum(tr)
      echo i, ": ", trace(g[i])/lo.physVol
  var pl = plaq(g)
  echo pl
  #for i in 0..2:
  #  echo g[0][0][i,i]
  #echo g[0][0][0,0]
  #echo g[0][256][0,0]

qexInit()
test()
qexFinalize()
