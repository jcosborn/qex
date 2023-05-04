import qex
import testutils
import strutils

proc echoRankS(c: Comm, xs: varargs[string,`$`]) =
  for i in 0..<nRanks:
    c.barrier
    if myRank == i:
      echoRank xs.join
  c.barrier

qexInit()
var c = getDefaultComm()
echo "rank ",myRank," / ",nRanks
echo "thread ",threadNum," / ",numThreads
let seed:uint64 = 7_005_003_002_001_000_000u64 + myRank.uint64
c.echoRankS "seed: ",seed

var
  lat = [8,8,8,8]
  #lat = latticeFromLocalLattice([8,8,8,8], nRanks)
  lo = newLayout(lat)
  p = lo.newGauge
  r = lo.newRNGField(RngMilc6, seed)
  R:RngMilc6
if p[0].numberType is float32: CT = 1e-9
R.seed(seed, 987654321)

p.randomTAH r
var p2:float
threads:
  var p2t = 0.0
  for i in 0..<p.len:
    p2t += p[i].norm2
  threadmaster: p2 = p2t
echo "p2: ",p2

var u = R.uniform
c.echoRankS "urand: ",u

var us = newseq[typeof u](nRanks)

if myRank == 0:
  for i in 1..<nRanks:
    c.pushRecv i, us[i]
else:
  c.pushSend 0, u
c.waitAll

if myRank == 0:
  us[0] = u
  suite "RNG seeding always from rank 0":
    test "random field":
      when defined(FUELCompat):
        let p2Expected = 130731.9294737126
      else:
        let p2Expected = 131563.7475902051
      check(p2 ~ p2Expected)
    test "global random number":
      let uExpected = 0.7708062529563904
      for i in 0..<nRanks:
        echo "u[",i,"]: ",us[i]
        check(us[i] ~ uExpected)

qexFinalize()
