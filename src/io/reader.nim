import base
import layout
import qio
import strutils
import macros
import field
#import comms
import os
#import stdUtils
import iocommon

proc getFileLattice*(s:string):seq[int] =
  if not existsFile(s):
    warn "file not found: ", s
    return
  var ql:QIO_Layout
  ql.this_node = myRank.cint
  ql.number_of_nodes = nRanks.cint
  var qs = QIO_string_create()
  var qr = QIO_open_read(qs, s, ql.addr, nil, nil)
  let nd = QIO_get_reader_latdim(qr)
  #echo nd
  let ls = QIO_get_reader_latsize(qr)
  result.newSeq(nd)
  for i in 0..<nd:
    result[i] = ls[i]
  discard QIO_close_read(qr)

type Reader*[V:static[int]] = ref object
  layout*:Layout[V]
  latsize*:seq[cint]
  status*:int
  fileName*:string
  fileMetadata*:string
  recordInfoValid:bool
  recordMd:string
  recordInfo:QIO_RecordInfo
  nioranklist*: seq[int32]
  iogeom*: seq[int32]
  iorank*: seq[int32]
  qr:ptr QIO_Reader
  verb*: cint

proc IOverb*(level: int) =
  discard QIO_verbose(level.cint)

proc ioLayout(x: Layout): Layout =
  var ioLayoutV{.global.}: Layout
  result = ioLayoutV
  if not x.isNil: ioLayoutV = x
template setLayout(x: Reader): untyped =
  discard ioLayout(x.layout)
template getLayout(x: static[int]): untyped =
  ioLayout(Layout[x](nil))

var readnodes = -1
var riorank: ptr seq[int32]
proc getNumReadRanks*(): int = readnodes
proc setNumReadRanks*(n: int) =
  readnodes = n
proc ioReadRank(node: cint): cint =
  riorank[][node]
proc ioMasterRank(): cint = 0

proc toString(qs:ptr QIO_String):string =
  let n = qs.length.int - 2  # seems to have 2 byte padding
  result = spaces(n)
  for i in 0..<n:
    result[i] = qs.string[i]

proc open(rd: var Reader; ql: var QIO_Layout) =
  let nd = rd.layout.nDim.cint
  ql.latdim = nd
  rd.latsize.newSeq(nd)
  for i in 0..<nd:
    rd.latsize[i] = rd.layout.physGeom[i].cint
  ql.latsize = rd.latsize[0].addr
  ql.volume = rd.layout.physVol.csize_t
  ql.sites_on_node = rd.layout.nSites.csize_t
  ql.this_node = rd.layout.myRank.cint
  ql.number_of_nodes = rd.layout.nRanks.cint

  rd.nioranklist = listNumIoRanks(rd.layout.rankGeom)
  if readnodes<=0:
    let n = rd.nioranklist.len
    readnodes = rd.niorankList[n div 2]
  var rdnodes = getClosestNumRanks(rd.nioranklist, int32 readnodes)
  echo "Read num nodes: ", rdnodes
  rd.iogeom = getIoGeom(rd.layout.rankGeom, rdnodes)
  echo "Read geom: ", rd.iogeom
  rd.iorank = getIoRanks(rd.layout, rd.iogeom)
  echo "Read IO ranks: ", rd.iorank
  riorank = addr rd.iorank

  var fs: QIO_Filesystem
  fs.my_io_node = ioReadRank
  fs.master_io_node = ioMasterRank

  var iflag: QIO_Iflag
  #iflag.serpar = QIO_SERIAL;
  iflag.serpar = QIO_PARALLEL
  #//iflag.volfmt = QIO_UNKNOWN;
  iflag.volfmt = QIO_SINGLEFILE

  var qioMd = QIO_string_create()
  rd.qr = QIO_open_read(qioMd, rd.fileName, ql.addr, fs.addr, iflag.addr)
  rd.fileMetadata = toString(qioMd)
  QIO_string_destroy(qioMd)
  rd.recordInfoValid = false
  if rd.qr==nil: rd.status = -1

template newReader*[V: static[int]](l: Layout[V]; fn: string): untyped =
  var rd: Reader[V]
  if not fileExists(fn):
    warn "file not found: ", fn
  else:
    proc ioNodeNumber2(x:ptr cint):cint =
      rankIndex(getLayout(V), x).rank.cint
    proc ioNodeIndex2(x:ptr cint):cint =
      rankIndex(getLayout(V), x).index.cint
    proc ioGetCoords2(x:ptr cint; node: cint; index: cint) =
      getLayout(V).coord(x, (node,index))
    proc ioNumSites2(node: cint):cint =
      return getLayout(V).nSites.cint
    var ql: QIO_Layout
    ql.node_number = ioNodeNumber2
    ql.node_index = ioNodeIndex2
    ql.get_coords = ioGetCoords2
    ql.num_sites = ioNumSites2
    rd.new()
    rd.layout = l
    rd.fileName = fn
    setLayout(rd)
    open(rd, ql)
  rd

proc close*(rd: var Reader) =
  rd.setLayout()
  riorank = addr rd.iorank
  rd.status = QIO_close_read(rd.qr)

proc readRecordInfo(rd: var Reader) =
  rd.setLayout()
  riorank = addr rd.iorank
  var qioMd = QIO_string_create()
  rd.status = QIO_read_record_info(rd.qr, rd.recordInfo.addr, qioMd)
  rd.recordMd = toString(qioMd)
  QIO_string_destroy(qioMd);
  rd.recordInfoValid = true

proc nextRecord*(rd: var Reader) =
  rd.setLayout()
  riorank = addr rd.iorank
  rd.status = QIO_next_record(rd.qr)
  rd.recordInfoValid = false

proc recordMetadata*(rd: var Reader):string =
  rd.setLayout()
  riorank = addr rd.iorank
  if not rd.recordInfoValid: rd.readRecordInfo()
  return rd.recordMd

template recGet(f,t:untyped):untyped =
  proc f*(r:var Reader):auto =
    r.setLayout()
    if not r.recordInfoValid: r.readRecordInfo()
    let x = `"QIO_get_" f`(r.recordInfo.addr)
    return t(x)
recGet(record_date, `$`)
recGet(datatype, `$`)
recGet(precision, `$`)
recGet(colors, int)
recGet(spins, int)
recGet(typesize, int)
recGet(datacount, int)

proc getWordSize(rd: var Reader):int =
  #var result = 0
  let c = precision(rd)
  case c:
    of "F": result = 4
    of "D": result = 8
    else: result = 0
  #result

import qioInternal

proc put[T](buf: cstring; index: csize_t, count: cint; arg: pointer) =
  type srcT = cArray[IOtype(T)]
  type destT1 = cArray[T]
  type destT = cArray[ptr destT1]
  let src = cast[ptr srcT](buf)
  let dest = cast[ptr destT](arg)
  let vi = index div simdLength(T).csize_t
  let vl = int(index mod simdLength(T).csize_t)
  let vlm = 1 shl vl
  for i in 0..<count:
    masked(dest[i][vi], vlm) := src[i]

proc putP[T](buf: cstring; index: csize_t, count: cint; arg: pointer) =
  type srcT = cArray[IOtypeP(T)]
  type destT1 = cArray[T]
  type destT = cArray[ptr destT1]
  let src = cast[ptr srcT](buf)
  let dest = cast[ptr destT](arg)
  let vi = index div simdLength(T).csize_t
  let vl = int(index mod simdLength(T).csize_t)
  let vlm = 1 shl vl
  for i in 0..<count:
    masked(dest[i][vi], vlm) := src[i]

proc read[T](rd: var Reader, v: var openArray[ptr T]) =
  rd.setLayout()
  riorank = addr rd.iorank
  var qioMd = QIO_string_create()

  rd.readRecordInfo()
  let recWordSize = getWordSize(rd)
  var size = sizeOf(v[0][]) div rd.layout.nSitesInner
  var vsize = size*v.len
  var wordSize = sizeOf(numberType(v[0][]))
  if wordSize==recWordSize:
    rd.status = QIO_read(rd.qr, rd.recordInfo.addr, qioMd, put[T],
                         vsize.csize_t, wordSize.cint, v[0].addr)
  else:
    vsize = (recWordSize*vsize) div wordSize
    wordSize = recWordSize
    rd.status = QIO_read(rd.qr, rd.recordInfo.addr, qioMd, putP[T],
                         vsize.csize_t, wordSize.cint, v[0].addr)

  rd.recordMd = toString(qioMd)
  QIO_string_destroy(qioMd);
  rd.recordInfoValid = true

proc read*[V:static[int],T](rd: var Reader[V]; v: openArray[Field[V,T]]) =
  let nv = v.len
  var f = newSeq[type(v[0][0].addr)](nv)
  for i in 0..<nv: f[i] = v[i][0].addr
  rd.read(f)

proc read*(rd: var Reader; v: Field) =
  var f = @[ v[0].addr ]
  rd.read(f)

when isMainModule:
  qexInit()
  echo "rank ", myRank, "/", nRanks
  var lat = [8,8,8,8]
  var lo = newLayout(lat)
  var g:array[4,type(lo.ColorMatrix())]
  for i in 0..<4: g[i] = lo.ColorMatrix()
  #var fn = "l88.scidac"
  const fn = "testlat0.bin"
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
