import base
import layout
import qio
#import strutils
import field
import os, times
import iocommon
import qioInternal

type Writer*[V: static[int]] = ref object
  layout*: Layout[V]
  latsize*: seq[cint]
  status*: int
  fileName*: string
  nioranklist*: seq[int32]
  iogeom*: seq[int32]
  iorank*: seq[int32]
  qw: ptr QIO_Writer

proc ioLayout(x: Layout): Layout =
  var ioLayoutV{.global.}: Layout
  result = ioLayoutV
  if not x.isNil: ioLayoutV = x
template setLayout(x: Writer): untyped =
  discard ioLayout(x.layout)
template getLayout(x: static[int]): untyped =
  ioLayout(Layout[x](nil))

var writenodes = -1
var wiorank: ptr seq[int32]
proc getNumWriteRanks*(): int = writenodes
proc setNumWriteRanks*(n: int) =
  writenodes = n
proc ioWriteRank(node: cint): cint =
  wiorank[][node]
proc ioMasterRank(): cint = 0

proc open(wr: var Writer; ql: var QIO_Layout, md: string) =
  let nd = wr.layout.nDim.cint
  ql.latdim = nd
  wr.latsize.newSeq(nd)
  for i in 0..<nd:
    wr.latsize[i] = wr.layout.physGeom[i].cint
  ql.latsize = wr.latsize[0].addr
  ql.volume = wr.layout.physVol.csize_t
  ql.sites_on_node = wr.layout.nSites.csize_t
  ql.this_node = wr.layout.myRank.cint
  ql.number_of_nodes = wr.layout.nRanks.cint

  wr.nioranklist = listNumIoRanks(wr.layout.rankGeom)
  if writenodes<=0:
    let n = wr.nioranklist.len
    let k = min(n-1, (n div 2)+1)
    writenodes = wr.niorankList[k]
  var wrnodes = getClosestNumRanks(wr.nioranklist, int32 writenodes)
  echo "Write num nodes: ", wrnodes
  wr.iogeom = getIoGeom(wr.layout.rankGeom, wrnodes)
  echo "Write geom: ", wr.iogeom
  wr.iorank = getIoRanks(wr.layout, wr.iogeom)
  echo "Write IO ranks: ", wr.iorank
  wiorank = addr wr.iorank

  var fs: QIO_Filesystem
  fs.my_io_node = ioWriteRank
  fs.master_io_node = ioMasterRank

  let volfmt = QIO_SINGLEFILE
  var oflag: QIO_Oflag
  #oflag.serpar = QIO_SERIAL;
  oflag.serpar = QIO_PARALLEL
  oflag.mode = QIO_TRUNC
  #oflag.ildgstyle = QIO_ILDGLAT
  oflag.ildgstyle = QIO_ILDGNO
  #oflag.ildgLFN = NULL

  var qioMd = QIO_string_create()
  QIO_string_set(qioMd, md)
  wr.setLayout
  wr.qw = QIO_open_write(qioMd, cstring wr.filename, volfmt,ql.addr,fs.addr,oflag.addr)
  QIO_string_destroy(qioMd)

  if wr.qw==nil: wr.status = -1

template newWriter*[V: static[int]](l: Layout[V]; fn,md: string): untyped =
  proc wioNodeNumber2(x: ptr ConstInt):cint =
    rankIndex(getLayout(V), x.constcast).rank.cint
  proc wioNodeIndex2(x: ptr ConstInt):cint =
    rankIndex(getLayout(V), x.constcast).index.cint
  proc wioGetCoords2(x: ptr cint; node: cint; index: cint) =
    getLayout(V).coord(x, (node,index))
  proc wioNumSites2(node: cint):cint =
    return getLayout(V).nSites.cint
  var ql: QIO_Layout
  ql.node_number = wioNodeNumber2
  ql.node_index = wioNodeIndex2
  ql.get_coords = wioGetCoords2
  ql.num_sites = wioNumSites2

  var wr: Writer[V]
  wr.new
  wr.layout = l
  wr.fileName = fn
  setLayout(wr)
  wr.open ql, md
  wr

proc close*(wr: var Writer) =
  wr.setLayout
  wiorank = addr wr.iorank
  wr.status = QIO_close_write(wr.qw)

proc get[T](buf: cstring; index: csize_t; count: cint; arg: pointer) =
  type destT = cArray[IOtype(T)]
  type srcT1 = cArray[T]
  type srcT = cArray[ptr srcT1]
  let dest = cast[ptr destT](buf)
  let src = cast[ptr srcT](arg)
  let vi = index div simdLength(T).csize_t
  let vl = int(index mod simdLength(T).csize_t)
  #let vlm = 1 shl vl
  #var s: ptr destT
  for i in 0..<count:
    #vcopy(dest[i], vl, src[i][vi])
    #let t = masked(src[i][vi], vlm)
    #dest[i] := t
    dest[i] := indexed(src[i][vi], asSimd(vl))

proc getP[T](buf: cstring; index: csize_t; count: cint; arg: pointer) =
  type destT = cArray[IOtypeP(T)]
  type srcT1 = cArray[T]
  type srcT = cArray[ptr srcT1]
  let dest = cast[ptr destT](buf)
  let src = cast[ptr srcT](arg)
  let vi = index div simdLength(T).csize_t
  let vl = int(index mod simdLength(T).csize_t)
  #let vlm = 1 shl vl
  #var s: ptr destT
  for i in 0..<count:
    #vcopy(dest[i], vl, src[i][vi])
    #dest[i] := masked(src[i][vi], vlm)
    dest[i] := indexed(src[i][vi], asSimd(vl))

proc write[T](wr: var Writer, v: var openArray[ptr T], lat: openArray[int],
              md="", ps="") =
  wr.setLayout
  wiorank = addr wr.iorank

  var qioMd = QIO_string_create()
  QIO_string_set(qioMd, md)

  var size = sizeOf(v[0][]) div wr.layout.nSitesInner
  var vsize = size*v.len
  var wordSize = sizeOf(numberType(v[0][]))
  #echo "size: ", size, "  vsize: ", vsize, "  wordsize: ", wordSize

  var precs0 = case wordSize:
    of 4: "F"
    of 8: "D"
    else: ""
  var precs = precs0
  if ps != "": precs = ps
  if (precs != "F" and precs != "D") or (precs0 != "F" and precs0 != "D"):
    echo "ERROR: Unsupported precision precs=", precs,
     " precs0=", precs0, " with wordSize=", wordSize
    qexExit 1

  var nc = v[0][].getNc.cint
  var ns = v[0][].getNs.cint
  var nv = v.len.cint
  let nd = lat.len
  var
    lower = newseq[cint](nd)
    upper = newseq[cint](nd)
  for i in 0..<nd:
    lower[i] = 0.cint
    upper[i] = lat[i].cint

  if precs == precs0:
    #var datatype = "QEX_" & type(v[0][]).IOtype.name
    var datatype = type(v[0][]).IOtype.IOname
    var recInfo = QIO_create_record_info(QIO_FIELD, lower[0].addr,
                                         upper[0].addr, nd.cint, cstring datatype,
                                         cstring precs, nc, ns, size.cint, nv)
    wr.status = QIO_write(wr.qw, recInfo, qioMd, get[T], vsize.csize_t,
                          wordSize.cint, v[0].addr)
    QIO_destroy_record_info(recInfo)
  else:
    #var datatype = "QEX_" & type(v[0][]).IOtypeP.name
    var datatype = type(v[0][]).IOtypeP.IOname
    var recWordSize = case precs:
      of "F": 4
      of "D": 8
      else: 0
    vsize = (recWordSize*vsize) div wordSize
    size = vsize div v.len
    wordSize = recWordSize
    var recInfo = QIO_create_record_info(QIO_FIELD, lower[0].addr,
                                         upper[0].addr, nd.cint, cstring datatype,
                                         cstring precs, nc, ns, size.cint, nv)
    wr.status = QIO_write(wr.qw, recInfo, qioMd, getP[T], vsize.csize_t,
                          wordSize.cint, v[0].addr)
    QIO_destroy_record_info(recInfo);

  QIO_string_destroy(qioMd);

proc write*[V:static[int],T](wr:var Writer[V]; v:openArray[Field[V,T]];
                             md = ""; prec = "") =
  let nv = v.len
  var f = newSeq[type(v[0][0].addr)](nv)
  for i in 0..<nv: f[i] = v[i][0].addr
  wr.write(f,v[0].l.physGeom,md,prec)

proc write*(wr: var Writer; v: Field; md=""; prec="") =
  var f = @[ v[0].addr ]
  wr.write(f, v.l.physGeom, md, prec)

when isMainModule:
  import gauge, rng, reader
  import physics/qcdTypes
  proc tr(f: any) =
    when f is array or f is seq:
      for i in 0..<f.len:
        echo "tr", i, ": ", trace(f[i])
    else:
      echo "tr: ", trace(f)
  const fn = "testlat0.bin"
  proc readTest(f: any) =
    when f is array or f is seq:
      let lo = f[0].l
    else:
      let lo = f.l
    var rd = lo.newReader(fn)
    echo rd.fileMetadata
    echo rd.recordMetadata
    echo rd.recordDate
    echo rd.datatype
    echo rd.precision
    echo rd.colors
    echo rd.spins
    echo rd.typesize
    echo rd.datacount
    rd.read(f)
    #rd.nextRecord()
    #echo rd.status
    #echo rd.recordMetadata
    rd.close()
    f.tr
  proc writeTest(f,prec: any) =
    when f is array or f is seq:
      let lo = f[0].l
    else:
      let lo = f.l
    var wr = lo.newWriter(fn, "filemd")
    wr.write(f, "recordmd", prec)
    wr.close

  #if not fileExists(fn):
  #  echo "gauge file not found: ", fn
  #  qexExit()
  qexInit()
  echo "rank ", myRank, "/", nRanks
  var lat = [8,8,8,8]
  var lo = newLayout(lat)
  var lo1 = newLayout(lat, 1)
  var rnd = newRNGField(RngMilc6, lo1)

  proc test(f: any) =
    f.gaussian rnd
    f.tr
    echo "#### Write double precision"
    f.writeTest "D"
    f.gaussian rnd
    f.readTest
    echo "#### Write single precision"
    f.writeTest "F"
    f.gaussian rnd
    f.readTest

  block:
    echo "testing array[4, ColorMatrix]"
    var g: array[4,type(lo.ColorMatrix())]
    for i in 0..<4:
      g[i] = lo.ColorMatrix()
    test(g)

  block:
    echo "testing Complex"
    var c = lo.Complex()
    test(c)

  qexFinalize()
