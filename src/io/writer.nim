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
proc ioWriteRank(node: cint): cint =
  cint( writenodes * (node div writenodes) )
proc ioMasterRank(): cint = 0.cint

proc open(wr: var Writer; ql: var QIO_Layout, md: string) =
  let nd = wr.layout.nDim.cint
  ql.latdim = nd
  wr.latsize.newSeq(nd)
  for i in 0..<nd:
    wr.latsize[i] = wr.layout.physGeom[i].cint
  ql.latsize = wr.latsize[0].addr
  ql.volume = wr.layout.physVol
  ql.sites_on_node = wr.layout.nSites
  ql.this_node = wr.layout.myRank.cint
  ql.number_of_nodes = wr.layout.nRanks.cint
  if writenodes<=0:
    writenodes = 1 + int( sqrt(ql.number_of_nodes.float) )

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
    var li: LayoutIndexQ
    layoutIndexQ(getLayout(V).lq.addr, li.addr, x)
    result = li.rank
  proc wioNodeIndex2(x: ptr cint):cint =
    var li: LayoutIndexQ
    layoutIndexQ(getLayout(V).lq.addr, li.addr, x)
    return li.index
  proc wioGetCoords2(x: ptr cint; node: cint; index: cint) =
    var li = LayoutIndexQ((node, index))
    layoutCoordQ(getLayout(V).lq.addr, x, li.addr)
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
  wr.setLayout
  wr.open ql, md
  wr

proc close*(wr: var Writer) =
  wr.setLayout
  wr.status = QIO_close_write(wr.qw)

proc write[T](wr: var Writer, v: var openArray[ptr T], md="", ps="") =
  wr.setLayout

  var qioMd = QIO_string_create()
  QIO_string_set(qioMd, md)

  #wr.writeRecordInfo()
  #let recWordSize = getWordSize(r)
  var size = sizeOf(v[0][]) div wr.layout.nSitesInner
  var vsize = size*v.len
  var wordSize = sizeOf(numberType(v[0][]))
  #echo "ws: ", wordSize, "  rws: ", recWordSize

  var precs0 = precString(wordSize)
  var precs = precs0
  if ps != "": precs = ps

  var nc = v[0][].nc
  var ns = v[0][].ns
  var qpct = qdpPCType(v[0][])
  var nv = v.len

  if precs == precs0:
    proc get(buf: cstring; index: csize; count: cint; arg: pointer) =
      type destT{.unchecked.} = array[0,IOtype(T)]
      type srcT1{.unchecked.} = array[0,T]
      type srcT{.unchecked.} = array[0,ptr srcT1]
      let dest = cast[ptr destT](buf)
      let src = cast[ptr srcT](arg)
      let vi = index div simdLength(T)
      let vl = index mod simdLength(T)
      for i in 0..<count:
        vcopy(dest[i], vl, src[i][vi])
    var recInfo = QIO_create_record_info(QIO_FIELD, 0, 0, 0, qpct, precs,
                                         nc, ns, size, nv)
  status = QIO_write(qdpw->qiow, rec_info, qio_md,
                     get, qf->size*count, qf->word_size, (void *)qf);


var status = QDP_write_check(qw, qioMd, QIO_FIELD, get, &qf, nv,
                                                       recInfo);

    QIO_destroy_record_info(rec_info);

r.status = QIO_write(r.qr, r.recordInfo.addr, qioMd, put, vsize.csize,
                        wordSize.cint, v[0].addr)
  else:
    vsize = (recWordSize*vsize) div wordSize
    wordSize = recWordSize
    #echo "ws: ", wordSize, "  rws: ", recWordSize, "  vs: ", vsize
    proc put(buf:cstring; index:csize; count:cint; arg:pointer) =
      type srcT{.unchecked.} = array[0,IOtypeP(T)]
      type destT1{.unchecked.} = array[0,T]
      type destT{.unchecked.} = array[0,ptr destT1]
      let src = cast[ptr srcT](buf)
      let dest = cast[ptr destT](arg)
      let vi = index div simdLength(T)
      let vl = index mod simdLength(T)
      for i in 0..<count:
        vcopy(dest[i][vi], vl, src[i])
    r.status = QIO_write(r.qr, r.recordInfo.addr, qioMd, put, vsize.csize,
                        wordSize.cint, v[0].addr)

  r.recordMd = toString(qioMd)
  QIO_string_destroy(qioMd);
  r.recordInfoValid = true

when isMainModule:
  import qex
  qexInit()
  echo "rank ", myRank, "/", nRanks
  var lat = [8,8,8,8]
  var lo = newLayout(lat)
  var g: array[4,type(lo.ColorMatrix())]
  for i in 0..<4:
    g[i] = lo.ColorMatrix()
  g.random
  var fn = "testlat0.bin"
  var wr = lo.newWriter(fn, "filemd")

  wr.write(g, "recordmd")

  wr.close

  if not fileExists(fn):
    echo "gauge file not found: ", fn
    qexExit()
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
  rd.read(g)
  #rd.nextRecord()
  #echo rd.status
  #echo rd.recordMetadata
  rd.close()
  var tr:type(g[0][0][0,0])
  for i in 0..<4:
    for x in g[i].all:
      for c in 0..<g[i][0].ncols:
        #echo g[i][x][c,c]
        tr += g[i][x][c,c]
  echo tr
  qexFinalize()
