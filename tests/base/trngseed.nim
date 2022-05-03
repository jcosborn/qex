import qex, comms/qmp
import testutils
import strutils

proc echoRankS(xs: varargs[string,`$`]) =
  for i in 0..<nRanks:
    QMP_barrier()
    if myRank == i:
      echoRank xs.join
  QMP_barrier()

qexInit()
echo "rank ",myRank," / ",nRanks
echo "thread ",threadNum," / ",numThreads
let seed:uint64 = 7_005_003_002_001_000_000u64 + myRank.uint64
echoRankS "seed: ",seed

var
  lo = newLayout([8,8,8,8])
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
echoRankS "urand: ",u

var us = newseq[typeof u](nRanks)

if nRanks > 1:
  var
    uh:QMP_msghandle_t
    um:QMP_msgmem_t
    mm = newseq[QMP_msgmem_t](nRanks-1)
  if myRank == 0:
    var mh = newseq[QMP_msghandle_t](nRanks-1)
    for i in 1..<nRanks:
      mm[i-1] = QMP_declare_msgmem(us[i].addr, sizeof(u).csize_t)
      mh[i-1] = QMP_declare_receive_from(mm[i-1], i.cint, 0)
    uh = QMP_declare_multiple(mh[0].addr, cint(nRanks-1))
  else:
    um = QMP_declare_msgmem(u.addr, sizeof(u).csize_t)
    uh = QMP_declare_send_to(um, 0, 0)
  discard QMP_start(uh)
  discard QMP_wait(uh)
  QMP_free_msghandle(uh)
  if myRank == 0:
    for i in 0..<nRanks-1:
      QMP_free_msgmem(mm[i])
  else:
    QMP_free_msgmem(um)

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
