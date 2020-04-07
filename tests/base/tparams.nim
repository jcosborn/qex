import testutils

import qex
qexInit()

suite "Test of params":
  test "default int":
    let t0 = myRank + 10
    let t1 = intParam("_", t0)
    echoAll myRank, "/", nRanks, ": ", t0, "  ", t1
    check(t1==10)

  test "default float":
    let t0 = myRank + 10.0
    let t1 = floatParam("_", t0)
    echoAll myRank, "/", nRanks, ": ", t0, "  ", t1
    check(t1==10.0)

  test "default string":
    let t0 = $(myRank + 10)
    let t1 = stringParam("_", t0)
    echoAll myRank, "/", nRanks, ": ", t0, "  ", t1
    check(t1=="10")

qexFinalize()
