import qex
import times, strformat, strutils, sequtils

type Dat = object
  vol: int
  bytes: int
  nio: int
  wsecs: float
  rsecs: float
var dat = newSeq[Dat]()
var byts = newSeq[int]()

proc test(lat: seq[int]) =
  let v = lat.foldl(a * b)
  let vpr = float(v) / float(nRanks)
  #echo "volume per rank: ", vpr
  if vpr < 1024 or vpr > 16*1024*1024: return
  var lo = newLayout(lat)
  var rs = newRNGField(RngMilc6, lo, intParam("seed", 987654321).uint64)
  var r0 = lo.RealS()
  var r1 = lo.RealS()
  var fn = "testqio.lat"
  let bytes = v * sizeof(r0.numberType)
  byts.add bytes
  echo "data size: ", bytes.formatSize

  for nio in 1..nRanks:
    echo "nio: ", nio
    var d = Dat(vol:v, bytes:bytes, nio:nio)
    threads:
      r0.gaussian rs
      r1 := 0

    block:
      setNumWriteRanks(nio)
      var wr = lo.newWriter(fn, "test")
      var t0 = epochTime()
      wr.write(r0)
      var t1 = epochTime()
      if wr.status!=0:
        echo "*** write status: ", wr.status
      echo "write seconds: ", $(t1-t0)
      wr.close
      if wr.status!=0:
        echo "*** write status: ", wr.status
      d.wsecs = t1 - t0

    block:
      setNumReadRanks(nio)
      var rd = lo.newReader(fn)
      var t0 = epochTime()
      rd.read(r1)
      var t1 = epochTime()
      if rd.status!=0:
        echo "*** read status: ", rd.status
      echo "read seconds: ", $(t1-t0)
      rd.close
      if rd.status!=0:
        echo "*** read status: ", rd.status
      d.rsecs = t1 - t0

    dat.add d
    let n0 = r0.norm2
    r1 -= r0
    let d2 = r1.norm2
    echo "norm2: ", n0
    echo "diff2: ", d2
    doAssert(d2 == 0.0)

proc run =
  qexInit()
  test(@[8,8,8,8])
  test(@[16,16,16,16])
  test(@[24,24,24,24])
  #test(@[32,32,32,32])
  #test(@[48,48,48,48])
  #test(@[64,64,64,64])
  test(@[96,96,96,96])
  test(@[128,128,128,128])
  test(@[192,192,192,192])
  test(@[256,256,256,256])
  var best = newSeq[string]()
  for b in byts:
    var wbmax,rbmax = 0.0
    var wnmax,rnmax = -1
    let bsi = b.formatSize
    echo "data size: ", bsi
    echo "IoRanks Write MB/s  Read MB/s"
    for i,d in dat:
      if d.bytes == b:
        let n = d.nio
        let wb = (1e-6*b)/d.wsecs
        let rb = (1e-6*b)/d.rsecs
        if wb>wbmax:
          wbmax = wb
          wnmax = n
        if rb>rbmax:
          rbmax = rb
          rnmax = n
        echo &"{n:7d} {wb:10.3f} {rb:10.3f}"
    best.add &"{bsi:>10s} {wnmax:7d} {wbmax:10.3f} {rnmax:7d} {rbmax:10.3f}"
  echo "Best:"
  echo "     Bytes IoRanks Write MB/s IoRanks Read MB/s"
  for s in best:
    echo s

  qexFinalize()

run()
