import qex
import gather

qexInit()
echo "rank ", myRank, "/", nRanks
var lat = intSeqParam("lat", @[8,8,8,16])
var t0 = intParam("t0", lat[^1] div 4)
var t1 = intParam("t1", (3*lat[^1]) div 4)

var lat2 = lat
lat2[^1] = t1-t0+1

var lo = newLayout(lat, 1)
var lo2 = newLayout(lat2, 1)

var f = lo.ColorVector1()
var f2 = lo2.ColorVector1()

var rl = newSeq[RecvList](0)
var x = newSeq[int32](lat.len)
for i in 0..<lo2.nSites:
  lo2.coord(x, i)
  let ri = lo.rankIndex(x)
  rl.add RecvList(didx: i.int32, srank: ri.rank.int32, sidx: ri.index.int32)

let c = getComm()
let gm = c.makeGatherMap(rl)

template `&&`(x: Field): untyped = cast[pointer](unsafeAddr(x[0]))

f := 1
f2 := 2
let fn1 = f.norm2
let f2n1 = f2.norm2
echo "f: ", fn1
echo "f2: ", f2n1
c.gather(gm, sizeof(f[0]), &&f, &&f2)
let fn2 = f.norm2
let f2n2 = f2.norm2
echo "f: ", fn2
echo "f2: ", f2n2
let f2n = (fn1*lat2[^1].float)/lat[^1].float
echo "f2: ", f2n

f2 := 2
c.gatherReversed(gm, sizeof(f[0]), &&f2, &&f)
let fn3 = f.norm2
let f2n3 = f2.norm2
echo "f: ", fn3
echo "f2: ", f2n3
let fn = fn1 - f2n + f2n3
echo "f: ", fn


qexFinalize()
