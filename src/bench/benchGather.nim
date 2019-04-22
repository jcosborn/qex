import qex, comms/gather

qexInit()
echo "rank ", myRank, "/", nRanks
var dlat = intSeqParam("llat", @[12,12,12,12])
var trank = nRanks
var mu = 0
while trank>1:
  trank = trank div 2
  dlat[mu] *= 2
  mu = (mu+1) mod dlat.len
var lat = intSeqParam("lat", dlat)

var t0 = intParam("t0", lat[^1] div 4)
var t1 = intParam("t1", (3*lat[^1]) div 4)

var lat2 = lat
lat2[^1] = t1-t0+1

var lo = newLayout(lat, 1)
var lo2 = newLayout(lat2, 1)

var f = lo.DiracFermion1()
var f2 = lo2.DiracFermion1()

tic()
var rl = newSeq[RecvList](0)
var x = newSeq[int32](lat.len)
for i in 0..<lo2.nSites:
  lo2.coord(x, i)
  let ri = lo.rankIndex(x)
  rl.add RecvList(didx: i.int32, srank: ri.rank.int32, sidx: ri.index.int32)
toc("rl setup")

let c = getComm()
let gm = c.makeGatherMap(rl)
toc("makeGatherMap")

template `&&`(x: Field): untyped = cast[pointer](unsafeAddr(x[0]))

f := 1
f2 := 2
let fn1 = f.norm2
let f2n1 = f2.norm2
echo "f: ", fn1
echo "f2: ", f2n1
toc("begin gather")
c.gather(gm, sizeof(f[0]), &&f, &&f2)
toc("end gather")
let fn2 = f.norm2
let f2n2 = f2.norm2
echo "f: ", fn2
echo "f2: ", f2n2
let f2n = (fn1*lat2[^1].float)/lat[^1].float
echo "f2: ", f2n

f2 := 2
toc("begin gatherReversed")
c.gatherReversed(gm, sizeof(f[0]), &&f2, &&f)
toc("end gatherReversed")
let fn3 = f.norm2
let f2n3 = f2.norm2
echo "f: ", fn3
echo "f2: ", f2n3
let fn = fn1 - f2n + f2n3
echo "f: ", fn

echoTimers()
qexFinalize()
