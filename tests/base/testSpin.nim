import times
import qex

proc test =
  #destructors.newSeq(0)
  #var lat = [4,4,2,2]
  #var lat = [4,4,4,4]
  #var lat = [8,8,4,4]
  var lat = [8,8,8,8]
  #var lat = [16,16,8,8]
  #var lat = [16,16,16,16]
  #var lat = [32,32,16,16]
  var lo = newLayout(lat)
  #layout.makeShift(0,1)
  #layout.makeShift(3,-2,"even")
  var m1 = lo.ColorMatrix()
  var d1 = lo.DiracFermion()
  var d2 = lo.DiracFermion()

  threads:
    m1 := 0
    d1 := 0
    d2 := 0
  echo "done init"

  threads:
    m1["odd"] := 1
    d1["odd"] := 1
    d2 := 2 * d1
    threadBarrier()
    echo m1.norm2/lo.nSites.float
    echo d1.norm2/lo.nSites.float
    echo d2.norm2/lo.nSites.float

  let nrep = int(1e9/lo.nSites.float)
  #let nrep = 1
  let t0 = epochTime()
  threads:
    for i in 1..nrep:
      d2 := m1 * d1
  let t1 = epochTime()
  echo "time: ", (t1-t0)
  echo d1[0][0][0]
  let ns = lo.nSitesOuter-1
  echo d1[ns][0][0]
  echo d2[ns][0][0]
  echo "mflops: ", (66e-6*lo.nSites.float*nrep.float)/(t1-t0)

  echo("mem: (used+free)/total: (", getOccupiedMem(), "+", getFreeMem(), ")/",
       getTotalMem())

  echo GC_getStatistics()
  GC_fullCollect()
  echo GC_getStatistics()
  #echo "destructors: ", destructors.len
  #for f in destructors: f()
  #echo GC_getStatistics()
  #GC_fullCollect()
  #echo GC_getStatistics()
  d1 = nil
  d2 = nil
  m1 = nil
  echo GC_getStatistics()
  GC_fullCollect()
  echo GC_getStatistics()
  echo("mem: (used+free)/total: (", getOccupiedMem(), "+", getFreeMem(), ")/",
       getTotalMem())


qexInit()
echo "rank ", myRank, "/", nRanks
echo threadNum, "/", numThreads
test()
qexFinalize()
