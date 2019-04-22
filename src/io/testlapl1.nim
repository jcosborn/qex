import qex
import parallelIo
import modfile
import physics/wilsonD
import physics/wilsonSolve
import xmlparser, xmltree, strutils, sequtils, endians
template `&&`(x: int32): untyped = cast[pointer](unsafeAddr(x))

qexInit()
echo "rank ", myRank, "/", nRanks
defaultSetup()

var src = lo.DiracFermion()
var dest = lo.DiracFermion()

var lo1 = newLayout(lat, 1, lo.rankGeom, @[1,1,1,1])
var cv1 = lo1.ColorVectorS1()

var srcfn = stringParam("srcfn", "colorvec.mod")
var pr = openRead(srcfn)
var srchdr = pr.modReadHeader()
var srcmap = pr.modReadMap(srchdr.mapstart)

var xml = parseXml(srchdr.userdata)
let ls = xml.findAll("lattSize")[0][0].text.split.mapit(parseInt(it))
let nv = xml.findAll("num_vecs")[0][0].text.parseInt
let dd = xml.findAll("decay_dir")[0][0].text.parseInt
let wght = xml.findAll("Weights")[0]
var weights = newSeq[seq[float]](nv)
for i in 0..<nv:
  weights[i] = wght[i][0].text.split.mapit(parseFloat(it))
echo "lattice size: ", ls
echo "num vecs: ", nv
echo "decay dir: ", dd
echo "weights: ", weights
let nt = ls[^1]

var nio = intParam("nio", sqrt(nRanks.float).int)
var ioranks = @[0]
for i in 1..<nio:
  ioranks.add( (i*(nRanks-1)) div (nio-1) )
var size = ls
size[^1] = 1
var offset = @[0,0,0,0]

var w = newWilson(g, src)
var m = floatParam("mass", 0.1)
var sp = initSolverParams()

proc setSpin(src: Field, v1: Field2, sp: int) =
  for i in src.sites:
    var cv: array[4,cint]
    src.l.coord(cv,(src.l.myRank,i))
    let j = v1.l.rankIndex(cv).index
    src{i}[sp] := v1[j]

pr.setBig32
#for t in 0..<nt:
for t in 0..<1:
  offset[^1] = t
  var wm = lo1.setupWrite(size, offset, ioranks)
  #var gm = setupGatherHyper(lo2, lo, offset)
  #for i in 0..<nv:
  for i in 0..<1:
    let p = srcmap.getPos(t,i)
    #echo i, " ", t, " ", p
    pr.seekSet(p)
    cv1 := 0
    pr.read(cv1, wm)
    echo cv1.norm2
    for s in 0..<4:
      src := 0
      src.setSpin(cv1, s)
      echo src.norm2
      #w.solve(dest, src, m, sp)
      #gm.gather(prop[n], dest)
pr.setSwap(0)



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
