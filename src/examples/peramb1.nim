import qex
import comms/gather
import io / [parallelIo, timesliceIo, modfile]
import physics / [wilsonD, wilsonSolve]
import contract, modeigs1
import xmlparser, xmltree, strutils, sequtils, endians, times, strformat
template `&&`(x: int32): untyped = cast[pointer](unsafeAddr(x))
template `&&`(x: Field): untyped = cast[pointer](unsafeAddr(x[0]))

qexInit()
#echo "rank ", myRank, "/", nRanks
defaultSetup()

var tsrc = intParam("tsrc", 0)
var dt = intParam("dt", (lat[^1] div 2))
var nv0 = intParam("nv", 10)
var srcfn = stringParam("srcfn", "testevecs.mod")

echo "lo1:"
var lo1 = newLayout(lat, 1, lo.rankGeom, @[1,1,1,1])
var cv1 = lo1.ColorVectorS1()

tic()
let tsio1 = newTimesliceIo(lo1)
toc("newTimesliceIo")
echo "ioRanks: ", tsio1.ioRanks

makeTestEigFile(srcfn, cv1, nv0, tsio1)
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

proc getLapl(cv: Field; mr: var ModFileReader; tsio: TimesliceIo;
             t,n,td: int) =
  let p = mr.map.getPos(t,n)
  #echo i, " ", t, " ", p
  tic()
  mr.r.seekSet(p)
  tsio.read(mr.r, cv, td)
  toc("tsio.read")
  #echo t, ",", n, ":", p, ": ", cv1.norm2

var comm = getComm()
var lat2 = lat
lat2[^1] = dt
var rg2 = intSeqParam("rg2", newSeq[int](0))
echo "lo2:"
var lo2 = newLayout(lat2, 1, rg2)
let tsio2 = newTimesliceIo(lo2)

type CvType = type(cv1[0])
type CvArray = ptr UncheckedArray[CvType]

let localSites2 = lo2.nSites
let nx2 = localSites2 div dt
let ncvsink = localSites2 * nv
let sink = cast[CvArray](alloc(ncvsink*sizeof(CvType)))  # [dt][nv][nx2][cv]

var map2 = newSeq[int32](localSites2)
var map2s = newSeq[int32](localSites2)
var xn = newSeq[int32](dt)
let nvx2 = nv * nx2
let nvx2s = nv * nx2 * 16
for i in 0..<localSites2:
  let t = lo2.coords[^1][i]
  map2[i] = int32(t*nvx2 + xn[t])
  map2s[i] = int32(t*nvx2s + xn[t])
  inc xn[t]
toc("make map2")

# read in the eigenvectors for the sink
var cv2 = lo2.ColorVectorS1()
for v in 0..<nv:
  let kv0 = v*nx2
  for t in 0..<dt:
    let ts = (tsrc + t) mod lo[^1]
    cv2.getLapl(mr, tsio2, ts, v, t)
    #echo &"{t} {v} {cv2.norm2}"
  threads:
    tfor i, 0..<localSites2:
      let k = kv0 + map2[i]
      sink[k] := cv2[i]
  #echo "cv201:"
  #echo cv2[0]
  #echo cv2[1]
#echo "sink01:"
#echo sink[0]
#echo sink[1]
#echo "----"
toc("read sink")

const nc = cv1[0].len
type Cmplx = type(cv1[0][0])
type CmplxArray = ptr UncheckedArray[Cmplx]

#[
let r = cast[CmplxArray](alloc(nv*nv*sizeof(Cmplx)))
for t in 0..<dt:
  let p = cast[CmplxArray](addr sink[t*nvx2])
  cmatmul(r, p, p, nv, nv, nx2*nc)
  for i in 0..<nv:
    for j in 0..<nv:
      echo &"{i} {j} {r[i*nv+j]}"
]#

var lomap1 = newSeq[int32](lo.nsites)
for i in lo.sites:
  var cv: array[4,cint]
  lo.coord(cv,(lo.myRank,i))
  let j = lo1.rankIndex(cv).index
  lomap1[i] = j.int32
toc("make map1")

let ncv0 = localSites2 * 16
let ncv1 = localSites2 * 16 * nv
let prp1 = cast[CvArray](alloc(ncv1*sizeof(CvType)))
# [dt][spinsrc][spindest][nv][nx2][cv]

var dest1 = lo1.DiracFermion1()
var dest2 = lo2.DiracFermion1()
proc chopLat(gm: GatherMap, dest: Field, s,n: int) =
  tic()
  threads:
    tfor i, 0..<lo.nsites:
      let k = lomap1[i]
      dest1[k] := dest{i}
  toc("chopLat copy")
  comm.gather(gm, sizeof(dest2[0]), &&dest2, &&dest1)
  toc("chopLat gather")
  echo "chop: ", dest.norm2, " -> ", dest2.norm2
  #echo dest1[0]
  #echo dest1[1]
  #echo dest2[0]
  #echo dest2[1]
  let kv0 = (s*4*nv + n)*nx2
  threads:
    tfor i, 0..<localSites2:
      let k = kv0 + map2s[i]
      for j in 0..<4:
        let l = k + j*nvx2
        prp1[l] := dest2[i][j]
  #[
  let kv0 = n*nx2
  threads:
    tfor i, 0..<localSites2:
      let k = kv0 + map2[i]
      prp1[k] := dest2[i][s]
  ]#
  #echo dest1[0][0]
  #echo dest2[0][0]
  #echo "prp101:"
  #echo prp1[0]
  #echo prp1[1]
  toc("chopLat pack")

let nccv1 = 16 * nv * nv
let prmb = cast[CmplxArray](alloc(nccv1*sizeof(Cmplx)))
proc naiveContract() =
  let n = nccv1*2
  let pr = cast[ptr type(prmb[0].re)](prmb)
  tic()
  for t in 0..<dt:
    tic()
    let s = cast[CmplxArray](addr sink[t*nvx2])
    let p = cast[CmplxArray](addr prp1[t*nvx2s])
    cmatmul(prmb, p, s, 16*nv, nv, nx2*nc)
    toc("matmul")
    comm.allReduce(pr, n)
    toc("reduce")
    #for i in 0..<16*nv:
    #  for j in 0..<nv:
    #    echo &"{i} {j} {prmb[i*nv+j]}"
  toc("naiveContract")


var src = lo.DiracFermion()
var dest = lo.DiracFermion()
let nt = ls[^1]

#var w = newWilson(g, src)
#var m = floatParam("mass", 0.1)
#var sp = initSolverParams()

proc setSpin(src: Field, v1: Field2, sp: int) =
  tic()
  threads:
    for i in src.sites:
      let j = lomap1[i]
      src{i}[sp] := v1[j]
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

for i in 0..<nv:
  threads:
    cv1 := 0
  getLapl(cv1, mr, tsio1, tsrc, i, tsrc)
  #echo "getLapl:"
  #echo cv1[0]
  #echo cv1[1]
  #echo cv1.norm2
  for s in 0..<4:
    getProp(dest, src, cv1, s)
    #echo dest[0][s]
    #echo dest[1][s]
    chopLat(gm, dest, s, i)
toc("end loop")

for i in 0..<nx2:
  let t = sink[i] - prp1[i]
  let s = t.norm2
  if s>1e-5: echo &"{i} {s}"
  #echo sink[i]
  #echo prp1[i]

toc("test1")

naiveContract()

toc("end")
echoTimers()
qexFinalize()
