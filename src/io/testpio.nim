import qex
import parallelIo

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
var lo = newLayout(lat, 1)

var f = lo.DiracFermion1()
f := 1
echo f.norm2
var f2 = lo.DiracFermion1()
echo f2.norm2

var t0 = lat[^1] - 1
var size = lat
var offset = newSeq[int](lat.len)
var nio = intParam("nio", 1)
var dior = @[0]
for i in 1..<nio:
  dior.add( (i*(nRanks-1)) div (nio-1) )
var ioranks = intSeqParam("ior", dior)
#size[^1] = 1
#offset[^1] = t0

var fn = stringParam("fn", "pwtest.bin")

echo "lat: ", lat
echo "size: ", size
echo "offset: ", offset
echo "ioranks: ", ioranks
tic()
var wm = lo.setupWrite(size, offset, ioranks)
toc("setupWrite")

var pw = openCreate(fn)
toc("openCreate")
pw.writeSingle("header::")
toc("write header")
pw.write(f, wm)
toc("write data")
pw.writeSingle("::midd::")
toc("write mid")
f := 2
toc("set field")
pw.write(f, wm)
toc("write data2")
pw.writeSingle("::footer")
toc("write footer")
pw.close()
toc("write close")

var pr = openRead(fn)
toc("openRead")
var header = newString(8)
pr.readSingle(header)
toc("read header")
echo "header: '", header, "'"
pr.read(f2, wm)
toc("read data")
echo f2.norm2
toc("echo norm2")
pr.readSingle(header)
toc("read mid")
echo "header: '", header, "'"
pr.read(f2, wm)
toc("read data2")
echo f2.norm2
toc("echo norm2")
pr.readSingle(header)
toc("read footer")
echo "header: '", header, "'"
pr.close()
toc("read close")

echoTimers()
qexFinalize()
