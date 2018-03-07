import base
import layout
import qio
import strutils
import macros
import field
#import comms
import os
#import stdUtils

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
  qr:ptr QIO_Reader

proc ioLayout(x: Layout): Layout =
  var ioLayoutV{.global.}: Layout
  result = ioLayoutV
  if not x.isNil: ioLayoutV = x
template setLayout(x: Reader): untyped =
  discard ioLayout(x.layout)
template getLayout(x: static[int]): untyped =
  ioLayout(Layout[x](nil))

var readnodes{.global.} = -1
proc ioReadRank(node: cint): cint =
  cint( readnodes * (node div readnodes) )
proc ioMasterRank(): cint = 0.cint

proc toString(qs:ptr QIO_String):string =
  let n = qs.length.int
  result = spaces(n)
  for i in 0..<n:
    result[i] = qs.string[i]

proc open(r:var Reader; ql:var QIO_Layout) =
  let nd = r.layout.nDim.cint
  ql.latdim = nd
  r.latsize.newSeq(nd)
  for i in 0..<nd:
    r.latsize[i] = r.layout.physGeom[i].cint
  ql.latsize = r.latsize[0].addr
  ql.volume = r.layout.physVol
  ql.sites_on_node = r.layout.nSites
  ql.this_node = r.layout.myRank.cint
  ql.number_of_nodes = r.layout.nRanks.cint
  if readnodes<=0:
    readnodes = 1 + int( sqrt(ql.number_of_nodes.float) )

  var fs:QIO_Filesystem
  fs.my_io_node = ioReadRank
  fs.master_io_node = ioMasterRank

  var iflag:QIO_Iflag
  #iflag.serpar = QIO_SERIAL;
  iflag.serpar = QIO_PARALLEL
  #//iflag.volfmt = QIO_UNKNOWN;
  iflag.volfmt = QIO_SINGLEFILE

  var qioMd = QIO_string_create()
  r.qr = QIO_open_read(qioMd, r.fileName, ql.addr, fs.addr, iflag.addr)
  r.fileMetadata = toString(qioMd)
  QIO_string_destroy(qioMd)
  r.recordInfoValid = false
  if r.qr==nil: r.status = -1

template newReader*[V: static[int]](l: Layout[V]; fn: string): untyped =
  var rd: Reader[V]
  if not existsFile(fn):
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

proc close*(r:var Reader) =
  r.setLayout()
  r.status = QIO_close_read(r.qr)

proc readRecordInfo(r:var Reader) =
  r.setLayout()
  var qioMd = QIO_string_create()
  r.status = QIO_read_record_info(r.qr, r.recordInfo.addr, qioMd)
  r.recordMd = toString(qioMd)
  QIO_string_destroy(qioMd);
  r.recordInfoValid = true

proc nextRecord*(r:var Reader) =
  r.setLayout()
  r.status = QIO_next_record(r.qr)
  r.recordInfoValid = false

proc recordMetadata*(r:var Reader):string =
  r.setLayout()
  if not r.recordInfoValid: r.readRecordInfo()
  return r.recordMd

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

proc getWordSize(r:var Reader):int =
  #var result = 0
  let c = precision(r)
  case c:
    of "F": result = 4
    of "D": result = 8
    else: result = 0
  #result

import qioInternal

template vcopyImpl(dest:typed; l:int; src:typed):untyped =
  for i in 0..<dest.nrows:
    for j in 0..<dest.ncols:
      dest[i,j].re[l] = src[i,j].re
      dest[i,j].im[l] = src[i,j].im
proc vcopy(dest:var SColorMatrixV; l:int; src:SColorMatrix) =
  vcopyImpl(dest,l,src)
proc vcopy(dest:var DColorMatrixV; l:int; src:DColorMatrix) =
  vcopyImpl(dest,l,src)
proc vcopy(dest:var SColorMatrixV; l:int; src:DColorMatrix) =
  vcopyImpl(dest,l,src)
proc vcopy(dest:var DColorMatrixV; l:int; src:SColorMatrix) =
  vcopyImpl(dest,l,src)

proc read[T](r:var Reader, v:var openArray[ptr T]) =
  r.setLayout()
  var qioMd = QIO_string_create()

  r.readRecordInfo()
  let recWordSize = getWordSize(r)
  var size = sizeOf(v[0][]) div r.layout.nSitesInner
  var vsize = size*v.len
  var wordSize = sizeOf(numberType(v[0][]))
  #echo "ws: ", wordSize, "  rws: ", recWordSize
  if wordSize==recWordSize:
    proc put(buf:cstring; index:csize; count:cint; arg:pointer) =
      type srcT{.unchecked.} = array[0,IOtype(T)]
      type destT1{.unchecked.} = array[0,T]
      type destT{.unchecked.} = array[0,ptr destT1]
      let src = cast[ptr srcT](buf)
      let dest = cast[ptr destT](arg)
      let vi = index div simdLength(T)
      let vl = int(index mod simdLength(T))
      for i in 0..<count:
        vcopy(dest[i][vi], vl, src[i])
    r.status = QIO_read(r.qr, r.recordInfo.addr, qioMd, put, vsize.csize,
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
      let vl = int(index mod simdLength(T))
      for i in 0..<count:
        vcopy(dest[i][vi], vl, src[i])
    r.status = QIO_read(r.qr, r.recordInfo.addr, qioMd, put, vsize.csize,
                        wordSize.cint, v[0].addr)

  r.recordMd = toString(qioMd)
  QIO_string_destroy(qioMd);
  r.recordInfoValid = true

proc read*[V:static[int],T](r:var Reader[V]; v:openArray[Field[V,T]]) =
  let nv = v.len
  var f = newSeq[type(v[0][0].addr)](nv)
  for i in 0..<nv: f[i] = v[i][0].addr
  r.read(f)

when isMainModule:
  qexInit()
  echo "rank ", myRank, "/", nRanks
  var lat = [8,8,8,8]
  var lo = newLayout(lat)
  var g:array[4,type(lo.ColorMatrix())]
  for i in 0..<4: g[i] = lo.ColorMatrix()
  var fn = "l88.scidac"
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
