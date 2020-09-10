import base
import layout
import qio
import strutils
import macros
import field
import os

type Writer*[V: static[int]] = ref object
  layout*: Layout[V]
  latsize*: seq[cint]
  status*: int
  fileName*: string
  ioranks*: cint
  iorank*: seq[cint]
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
var wiorank: ptr seq[cint]
proc getNumWriteRanks*(): int = writenodes
proc setNumWriteRanks*(n: int) =
  writenodes = n
proc ioWriteRank(node: cint): cint =
  #cint( writenodes * (node div writenodes) )
  #let k = (writenodes*node) div nRanks
  #cint( (k*nRanks+writenodes-1) div writenodes )
  wiorank[][node]
proc ioMasterRank(): cint = 0.cint

proc open(wr: var Writer; ql: var QIO_Layout, md: string) =
  let nd = wr.layout.nDim.cint
  ql.latdim = nd
  wr.latsize.newSeq(nd)
  var iogeom = newSeq[cint](nd)
  for i in 0..<nd:
    wr.latsize[i] = wr.layout.physGeom[i].cint
    iogeom[i] = if 2*i<nd: 1.cint else: wr.layout.rankGeom[i].cint
  ql.latsize = wr.latsize[0].addr
  ql.volume = wr.layout.physVol.csize_t
  ql.sites_on_node = wr.layout.nSites.csize_t
  ql.this_node = wr.layout.myRank.cint
  ql.number_of_nodes = wr.layout.nRanks.cint
  wr.iorank = newSeq[cint](wr.layout.nranks)
  var coords = newSeq[cint](nd)
  for r in 0..<wr.layout.nranks:
    wr.layout.coord(coords, r, 0)
    for i in 0..<nd:
      let t = (coords[i]*iogeom[i]) div wr.layout.physGeom[i]
      coords[i] = (t*wr.layout.physGeom[i]).int32 div iogeom[i]
    let ri = wr.layout.rankindex(coords)
    wr.iorank[r] = ri.rank.int32
  echo wr.iorank
  wiorank = addr wr.iorank
  echo "Write geom: ", iogeom
  #if writenodes<=0:
    #writenodes = int( sqrt(8*ql.number_of_nodes.float) )
    #writenodes = max(1,min(ql.number_of_nodes,writenodes))
    #let rg = wr.layout.rankGeom
    #writenodes = rg[^1]

  var fs: QIO_Filesystem
  fs.my_io_node = ioWriteRank
  fs.master_io_node = ioMasterRank

  let volfmt = QIO_SINGLEFILE
  var oflag: QIO_Oflag
  #oflag.serpar = QIO_SERIAL;
  oflag.serpar = QIO_PARALLEL
  oflag.mode = QIO_TRUNC
  oflag.ildgstyle = QIO_ILDGLAT
  #oflag.ildgLFN = NULL

  var qioMd = QIO_string_create()
  QIO_string_set(qioMd, md)
  wr.setLayout
  wr.qw = QIO_open_write(qioMd, wr.filename, volfmt,ql.addr,fs.addr,oflag.addr)
  QIO_string_destroy(qioMd)

  if wr.qw==nil: wr.status = -1

template newWriter*[V: static[int]](l: Layout[V]; fn,md: string): untyped =
  proc wioNodeNumber2(x: ptr cint):cint =
    rankIndex(getLayout(V), x).rank.cint
  proc wioNodeIndex2(x: ptr cint):cint =
    rankIndex(getLayout(V), x).index.cint
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
  wr.status = QIO_close_write(wr.qw)

import typetraits
import qioInternal

proc get[T](buf: cstring; index: csize_t; count: cint; arg: pointer) =
  type destT = cArray[IOtype(T)]
  type srcT1 = cArray[T]
  type srcT = cArray[ptr srcT1]
  let dest = cast[ptr destT](buf)
  let src = cast[ptr srcT](arg)
  let vi = index div simdLength(T).csize_t
  let vl = int(index mod simdLength(T).csize_t)
  let vlm = 1 shl vl
  #var s: ptr destT
  for i in 0..<count:
    #vcopy(dest[i], vl, src[i][vi])
    let t = masked(src[i][vi], vlm)
    dest[i] := t

proc getP[T](buf: cstring; index: csize_t; count: cint; arg: pointer) =
  type destT = cArray[IOtypeP(T)]
  type srcT1 = cArray[T]
  type srcT = cArray[ptr srcT1]
  let dest = cast[ptr destT](buf)
  let src = cast[ptr srcT](arg)
  let vi = index div simdLength(T).csize_t
  let vl = int(index mod simdLength(T).csize_t)
  let vlm = 1 shl vl
  #var s: ptr destT
  for i in 0..<count:
    #vcopy(dest[i], vl, src[i][vi])
    dest[i] := masked(src[i][vi], vlm)

proc write[T](wr: var Writer, v: var openArray[ptr T], lat: openArray[int],
              md="", ps="") =
  wr.setLayout
  wiorank = addr wr.iorank

  var qioMd = QIO_string_create()
  QIO_string_set(qioMd, md)

  #wr.writeRecordInfo()
  #let recWordSize = getWordSize(r)
  var size = sizeOf(v[0][]) div wr.layout.nSitesInner
  var vsize = size*v.len
  var wordSize = sizeOf(numberType(v[0][]))
  #echo "ws: ", wordSize, "  rws: ", recWordSize

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
    var datatype = "QEX_" & type(v[0][]).IOtype.name
    var recInfo = QIO_create_record_info(QIO_FIELD, lower[0].addr,
                                         upper[0].addr, nd.cint, datatype,
                                         precs, nc, ns, size.cint, nv)
    wr.status = QIO_write(wr.qw, recInfo, qioMd, get[T], vsize.csize_t,
                          wordSize.cint, v[0].addr)
    QIO_destroy_record_info(recInfo)
  else:
    var datatype = "QEX_" & type(v[0][]).IOtypeP.name
    var recWordSize = case precs:
      of "F": 4
      of "D": 8
      else: 0
    vsize = (recWordSize*vsize) div wordSize
    size = vsize div v.len
    wordSize = recWordSize
    var recInfo = QIO_create_record_info(QIO_FIELD, lower[0].addr,
                                         upper[0].addr, nd.cint, datatype,
                                         precs, nc, ns, size.cint, nv)
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
