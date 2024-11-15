import qex
import io/parallelIo, comms/gather
import io/modfile
import io/timesliceIo
import physics/wilsonD
import physics/wilsonSolve
import contract
import xmlparser, xmltree, strutils, sequtils, times
template `&&`(x: int32): untyped = cast[pointer](unsafeAddr(x))
template `&&`(x: Field): untyped = cast[pointer](unsafeAddr(x[0]))

template TAG(x: untyped): untyped = xmltree.`<>`(x)
template TXT(x: untyped): untyped = xmltree.newText(x)

proc makeTestEigFile(fn: string, f: Field, n: int, tsio: TimesliceIo) =
  let lattSize = f.l.physGeom.join(" ")
  let rankGeom = f.l.rankGeom.join(" ")
  let localGeom = f.l.localGeom.join(" ")
  let numVecs = $n
  #let runDate = "05 Jul 13 16:45:15 EDT"
  let runDate = $now()
  let totalVolume = $f.l.physVol
  var ud =
    TAG MODMetaData(
      TAG id(TXT "eigenVecsTimeSlice"),
      TAG lattSize(TXT lattSize),
      TAG decay_dir(TXT "3"),
      TAG num_vecs(TXT numVecs),
      TAG ProgramInfo(
        TAG code_version(
          TAG basePrecision(TXT "32")
        ),
        TAG run_date(TXT runDate),
        TAG Setgeom(
          TAG latt_size(TXT lattSize),
          TAG logical_size(TXT rankGeom),
          TAG subgrid_size(TXT localGeom),
          TAG total_volume(TXT totalVolume),
          TAG subgrid_volume(TXT totalVolume)
        )
      )
    )
  var weights = TAG Weights()
  let nt = f.l[^1]
  var ws = newSeq[float](nt)
  for i in 0..<n:
    for j in 0..<nt:
      ws[j] = float(i+j+1)
    let w = ws.join(" ")
    weights.add TAG elem(TXT w)
  ud.add weights
  #echo ud

  var mw = newModFileWriter(fn, $ud)
  var r = newRngField(f.l, RngMilc6)
  for i in 0..<n:
    f.gaussian(r)
    for t in 0..<nt:
      #echo "begin write: ", mw.w.pos
      mw.beginWrite(packKey(@[t,i]))
      tsio.write(mw.w, f, t)
      mw.endWrite()
      #echo "end write: ", mw.w.pos, "  cksum: ", mw.w.crc32
  mw.close()

qexInit()
#echo "rank ", myRank, "/", nRanks
defaultSetup()
var comm = getComm()

var tsrc = intParam("tsrc", 0)
var dt = intParam("dt", (lat[^1] div 2))
var nv0 = intParam("nv", 10)
var srcfn = stringParam("srcfn", "testevecs.mod")

var lat2 = lat
lat2[^1] = dt

var lo1 = newLayout(lat, 1, lo.rankGeom, @[1,1,1,1])
var rg2 = intSeqParam("rg2", newSeq[int](0))
var lo2 = newLayout(lat2, 1, rg2)

var cv1 = lo1.ColorVectorS1()
var src = lo.DiracFermion()
var dest = lo.DiracFermion()
var dest1 = lo1.DiracFermion1()
var dest2 = lo2.DiracFermion1()

tic()
let tsio = newTimesliceIo(lo1)
toc("newTimesliceIo")
echo "ioRanks: ", tsio.ioRanks

makeTestEigFile(srcfn, cv1, nv0, tsio)
toc("makeTestEigFile")

var mr = newModFileReader(srcfn)
toc("newModFileReader")

var xml = parseXml(mr.hdr.userdata)
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
echoParams()

let nt = ls[^1]

#var w = newWilson(g, src)
#var m = floatParam("mass", 0.1)
#var sp = initSolverParams()

proc getLapl(cv: Field; mr: var ModFileReader; tsio: TimesliceIo; t,n: int) =
  let p = mr.map.getPos(t,n)
  #echo i, " ", t, " ", p
  tic()
  mr.r.seekSet(p)
  tsio.read(mr.r, cv1, t)
  toc("tsio.read")
  #echo t, ",", n, ":", p, ": ", cv1.norm2

var lomap1 = newSeq[int32](lo.nsites)
for i in lo.sites:
  var cv: array[4,cint]
  lo.coord(cv,(lo.myRank,i))
  let j = lo1.rankIndex(cv).index
  lomap1[i] = j.int32

proc setSpin(src: Field, v1: Field2, sp: int) =
  tic()
  threads:
    for i in src.sites:
      let j = lomap1[i]
      src{i}[sp] = v1[j]
  toc("setSpin")

proc getProp(dest: Field, src: Field2, cv: Field3, s: int) =
  tic()
  threads:
    src := 0
  src.setSpin(cv, s)
  #echo src.norm2
  #w.solve(dest, src, m, sp)
  threads:
    dest := src
  toc("getProp")

toc("begin gm")
var rl = newSeq[RecvList](0)
var x = newSeq[int32](lat.len)
for i in 0..<lo2.nSites:
  lo2.coord(x, i)
  x[^1] = ((x[^1] + tsrc) mod lo1[^1]).int32
  let ri = lo1.rankIndex(x)
  rl.add RecvList(didx: i.int32, srank: ri.rank.int32, sidx: ri.index.int32)
toc("RecvList")
let gm = comm.makeGatherMap(rl)
toc("makeGatherMap")

type CvType = type(cv1[0])
type CvArray = ptr UncheckedArray[CvType]
type Cmplx = type(cv1[0][0])
type CmplxArray = ptr UncheckedArray[Cmplx]
let nmomsh = 1
let localSites2 = lo2.nSites
let ncv0 = localSites2 * 16
let ncv1 = localSites2 * 16 * nv
let ncv2 = ncv1 * nmomsh
let prp1 = cast[CvArray](alloc(ncv1*sizeof(CvType)))
let prp2 = cast[CvArray](alloc(ncv2*sizeof(CvType)))

proc chopLat(dest: Field, s,n: int) =
  tic()
  threads:
    tfor i, 0..<lo1.nsites:
      dest1[i] := dest{i}
  toc("chopLat copy")
  comm.gather(gm, sizeof(dest2[0]), &&dest2, &&dest1)
  toc("chopLat gather")
  #echo "chop: ", dest.norm2, " -> ", dest2.norm2
# [vn][srcspin][destspin][site][cv]
# [mom,shift][vn][srcspin][destspin][site][cv]
  let off1 = (n*16 + s*4)*localSites2
  let p1 = cast[CvArray](addr prp1[off1])
  let p2 = cast[CvArray](addr prp2[off1])
  threads:
    tfor i, 0..<localSites2:
      for j in 0..<4:
        p1[j*localSites2+i] := dest2[i][j]
        p2[j*localSites2+i] := dest2[i][j]
  toc("chopLat pack")

let nccv1 = 16 * nv
let nccv2 = nccv1 * nmomsh
let gp = cast[CmplxArray](alloc(nccv1*nccv2*sizeof(Cmplx)))
proc naiveContract() =
  tic()
  #[
  #for i in 0..<nccv1:
  tfor i, 0, nccv1-1:
    let p1 = cast[CvArray](addr prp1[i*localSites2])
    let g1 = cast[CmplxArray](addr gp[i*nccv2])
    for j in 0..<nccv2:
      let p2 = cast[CvArray](addr prp2[j*localSites2])
      var t: type(dot(p1[i], p2[j]))
      for k in 0..<localSites2:
        t += dot(p1[i], p2[j])
      g1[j] = t
  ]#
  #template numberType(x: ComplexType): untyped = numberType(x.re)
  #template numberType(x: Color): untyped = numberType(x[])
  template `&|`(x: typed): untyped = cast[ptr CmplxA[numberType(x[0])]](x)
  cmatmul(&|gp, &|prp1, &|prp2, nccv1, nccv2, 3*localSites2)
  toc("naiveContract")
  let n = nccv1*nccv2*2
  let p = cast[ptr type(gp[0].re)](gp)
  comm.allReduce(p, n)
  toc("reduce")

toc("begin loop")

for i in 0..<nv:
  threads:
    cv1 := 0
  getLapl(cv1, mr, tsio, tsrc, i)
  #echo cv1.norm2
  for s in 0..<4:
    getProp(dest, src, cv1, s)
    chopLat(dest, s, i)

toc("end loop")

naiveContract()

#[
toc("begin loop")
for t in 0..<nt:
#for t in 0..<1:
  for i in 0..<nv:
  #for i in 0..<1:
    cv1 := 0
    getLapl(cv1, mr, tsio, t, i)
    #echo cv1.norm2
    for s in 0..<4:
      getProp(dest, src, cv1, s)
      chopLat(dest)
toc("end loop")
]#

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

toc("end")
echoTimers()
qexFinalize()
