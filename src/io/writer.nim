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
  setLayout(wr)
  wr.open ql, md
  wr

proc close*(wr: var Writer) =
  wr.setLayout
  wr.status = QIO_close_write(wr.qw)

import typetraits
import qioInternal

template vcopyImpl(dest:typed; l:int; src:typed):untyped =
  for i in 0..<dest.nrows:
    for j in 0..<dest.ncols:
      dest[i,j].re = src[i,j].re[l]
      dest[i,j].im = src[i,j].im[l]
proc vcopy(dest:var SColorMatrix; l:int; src:SColorMatrixV) =
  vcopyImpl(dest,l,src)
proc vcopy(dest:var DColorMatrix; l:int; src:DColorMatrixV) =
  vcopyImpl(dest,l,src)
proc vcopy(dest:var SColorMatrix; l:int; src:DColorMatrixV) =
  vcopyImpl(dest,l,src)
proc vcopy(dest:var DColorMatrix; l:int; src:SColorMatrixV) =
  vcopyImpl(dest,l,src)

proc write[T](wr: var Writer, v: var openArray[ptr T], lat: openArray[int], md="", ps="") =
  wr.setLayout

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
    echo "ERROR: Unsupported precision precs=",precs," precs0=",precs0," with wordSize=",wordSize
    qexExit 1

  var nc = v[0][].ncols.cint
  var ns = 0.cint    # Doesn't matter for gauge fields.
  var datatype = "QEX" & type(v[0][]).name
  var nv = v.len.cint
  let nd = lat.len
  var
    lower = newseq[cint](nd)
    upper = newseq[cint](nd)
  for i in 0..<nd:
    lower[i] = 0.cint
    upper[i] = lat[i].cint

  if precs == precs0:
    proc get(buf: cstring; index: csize; count: cint; arg: pointer) =
      type destT{.unchecked.} = array[0,IOtype(T)]
      type srcT1{.unchecked.} = array[0,T]
      type srcT{.unchecked.} = array[0,ptr srcT1]
      let dest = cast[ptr destT](buf)
      let src = cast[ptr srcT](arg)
      let vi = index div simdLength(T)
      let vl = int(index mod simdLength(T))
      for i in 0..<count:
        vcopy(dest[i], vl, src[i][vi])
    var recInfo = QIO_create_record_info(QIO_FIELD, lower[0].addr, upper[0].addr, nd.cint,
      datatype, precs, nc, ns, size.cint, nv)
    wr.status = QIO_write(wr.qw, recInfo, qioMd, get, vsize.csize, wordSize.cint, v[0].addr)
    QIO_destroy_record_info(recInfo);
  else:
    var recWordSize = case precs:
      of "F": 4
      of "D": 8
      else: 0
    vsize = (recWordSize*vsize) div wordSize
    size = vsize div v.len
    wordSize = recWordSize
    # echo "ws: ", wordSize, "  rws: ", recWordSize, "  vs: ", vsize
    proc get(buf: cstring; index: csize; count: cint; arg: pointer) =
      type destT{.unchecked.} = array[0,IOtypeP(T)]
      type srcT1{.unchecked.} = array[0,T]
      type srcT{.unchecked.} = array[0,ptr srcT1]
      let dest = cast[ptr destT](buf)
      let src = cast[ptr srcT](arg)
      let vi = index div simdLength(T)
      let vl = int(index mod simdLength(T))
      for i in 0..<count:
        vcopy(dest[i], vl, src[i][vi])
    var recInfo = QIO_create_record_info(QIO_FIELD, lower[0].addr, upper[0].addr, nd.cint,
      datatype, precs, nc, ns, size.cint, nv)
    wr.status = QIO_write(wr.qw, recInfo, qioMd, get, vsize.csize, wordSize.cint, v[0].addr)
    QIO_destroy_record_info(recInfo);

  QIO_string_destroy(qioMd);

proc write*[V:static[int],T](wr:var Writer[V]; v:openArray[Field[V,T]]; md = ""; prec = "") =
  let nv = v.len
  var f = newSeq[type(v[0][0].addr)](nv)
  for i in 0..<nv: f[i] = v[i][0].addr
  wr.write(f,v[0].l.physGeom,md,prec)

when isMainModule:
  import gauge, rng, reader
  proc t(g:any) =
    var tr:type(g[0][0][0,0])
    for i in 0..<4:
      for x in g[i].all:
        for c in 0..<g[i][0].ncols:
          #echo g[i][x][c,c]
          tr += g[i][x][c,c]
    echo tr
  const fn = "testlat0.bin"
  proc readTest(lo:any) =
    var g: array[4,type(lo.ColorMatrix())]
    for i in 0..<4:
      g[i] = lo.ColorMatrix()
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
    g.t
  proc writeTest(g,prec:any) =
    var wr = g[0].l.newWriter(fn, "filemd")
    wr.write(g, "recordmd", prec)
    wr.close

  if not fileExists(fn):
    echo "gauge file not found: ", fn
    qexExit()
  qexInit()
  echo "rank ", myRank, "/", nRanks
  var lat = [8,8,8,8]
  var lo = newLayout(lat)
  var g: array[4,type(lo.ColorMatrix())]
  for i in 0..<4:
    g[i] = lo.ColorMatrix()
  g.random
  g.t

  echo "#### Write double precision"
  g.writeTest "D"
  lo.readTest

  echo "#### Write single precision"
  g.writeTest "F"
  lo.readTest

  qexFinalize()
