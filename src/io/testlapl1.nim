import qex
import parallelIo
import modfile

qexInit()
echo "rank ", myRank, "/", nRanks
defaultSetup()

var cv = lo.ColorVector()

var lo1 = newLayout(lat, 1, lo.rankGeom, @[1,1,1,1])
var cv1 = lo1.ColorVector1()

var srcfn = stringParam("srcfn", "colorvec.mod")
var pr = openRead(srcfn)
var srchdr = pr.modReadHeader()





#[
var f = lo.DiracFermion1()
f := 1
echo f.norm2
var f2 = lo.DiracFermion1()
echo f2.norm2

var t0 = lat[^1] - 1
var size = lat
var offset = newSeq[int](lat.len)
var ioranks = intSeqParam("ior", @[0])
#size[^1] = 1
#offset[^1] = t0

echo "lat: ", lat
echo "size: ", size
echo "offset: ", offset
echo "ioranks: ", ioranks
var wm = lo.setupWrite(size, offset, ioranks)

var pw = openCreate("pwtest.bin")
pw.writeSingle("header::")
pw.write(f, wm)
pw.writeSingle("::midd::")
f := 2
pw.write(f, wm)
pw.writeSingle("::footer")
pw.close()

var pr = openRead("pwtest.bin")
var header = newString(8)
pr.readSingle(header)
echo "header: '", header, "'"
pr.read(f2, wm)
echo f2.norm2
pr.readSingle(header)
echo "header: '", header, "'"
pr.read(f2, wm)
echo f2.norm2
pr.readSingle(header)
echo "header: '", header, "'"
pr.close()
]#

qexFinalize()
