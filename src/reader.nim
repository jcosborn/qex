import layout
import qio
import strutils
import macros
import field
import comms
import os
import stdUtils

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

proc ioLayout[V:static[int]](x:Layout[V]):Layout[V] =
  var ioLayoutV{.global.}:Layout[x.V]
  result = ioLayoutV
  if not x.isNil: ioLayoutV = x
template setLayout(x:Reader):untyped =
  discard ioLayout(x.layout)
proc ioNodeNumber[V:static[int]](x:ptr cint):cint =
  result = 0
  #var li:LayoutIndexQ
  #layoutIndexQ(ioLayout[V](nil).lq.addr, li.addr, x)
  #result = li.rank
  discard
proc ioNodeIndex[V:static[int]](x:ptr cint):cint =
  var li:LayoutIndexQ
  layoutIndexQ(ioLayout[V]().lq.addr, li.addr, x)
  return li.index
proc ioGetCoords[V:static[int]](x:ptr cint; node: cint; index: cint) =
  var li = LayoutIndexQ((node, index))
  layoutCoordQ(ioLayout[V]().lq.addr, x, li.addr)
proc ioNumSites[V:static[int]](node: cint):cint =
  return ioLayout[V]().nSites.cint
proc ioReadRank*(node: cint): cint = 0.cint
proc ioMasterRank*(): cint = 0.cint

proc toString(qs:ptr QIO_String):string =
  let n = qs.length.int
  result = spaces(n)
  for i in 0..<n:
    result[i] = qs.string[i]

proc open(r:var Reader; ql:var QIO_Layout) =
  #r.setLayout()
  #dumpTree: "test"
  #var layout:QIO_Layout
  #layout.node_number = ioNodeNumber2
  #layout.node_index = ioNodeIndex
  #layout.get_coords = ioGetCoords
  #layout.num_sites = ioNumSites
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
  #echo ql.volume
  #echo ql.sites_on_node
  #echo ql.this_node
  #echo ql.number_of_nodes
  
  var fs:QIO_Filesystem
  fs.my_io_node = ioReadRank
  fs.master_io_node = ioMasterRank
  
  var iflag:QIO_Iflag
  iflag.serpar = QIO_SERIAL;
  #iflag.serpar = QIO_PARALLEL
  #//iflag.volfmt = QIO_UNKNOWN;
  iflag.volfmt = QIO_SINGLEFILE

  var qioMd = QIO_string_create()
  r.qr = QIO_open_read(qioMd, r.fileName, ql.addr, fs.addr, iflag.addr)
  r.fileMetadata = toString(qioMd)
  QIO_string_destroy(qioMd)
  r.recordInfoValid = false
  if r.qr==nil: r.status = -1

template newReader*[V:static[int]](l:Layout[V]; fn:string):expr =
  var iol{.global.}:type(l) = nil
  var rd:Reader[V]
  if not existsFile(fn):
    warn "file not found: ", fn
  else:
    iol = l
    proc ioNodeNumber2(x:ptr cint):cint =
      var li:LayoutIndexQ
      layoutIndexQ(iol.lq.addr, li.addr, x)
      result = li.rank
    proc ioNodeIndex2(x:ptr cint):cint =
      var li:LayoutIndexQ
      layoutIndexQ(iol.lq.addr, li.addr, x)
      return li.index
    proc ioGetCoords2(x:ptr cint; node: cint; index: cint) =
      var li = LayoutIndexQ((node, index))
      layoutCoordQ(iol.lq.addr, x, li.addr)
    proc ioNumSites2(node: cint):cint =
      return iol.nSites.cint
    var ql:QIO_Layout
    ql.node_number = ioNodeNumber2
    ql.node_index = ioNodeIndex2
    ql.get_coords = ioGetCoords2
    ql.num_sites = ioNumSites2
    rd.new()
    rd.layout = l
    rd.fileName = fn
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
  if not r.recordInfoValid: r.readRecordInfo()
  return r.recordMd

template recGet(f,t:untyped):untyped =
  proc f*(r:var Reader):auto =
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

import qex
import qcdTypes

template IOtype(x:typedesc[SColorMatrixV]):typedesc = SColorMatrix
template IOtype(x:typedesc[DColorMatrixV]):typedesc = DColorMatrix
template IOtypeP(x:typedesc[SColorMatrixV]):typedesc = DColorMatrix
template IOtypeP(x:typedesc[DColorMatrixV]):typedesc = SColorMatrix

template vcopyImpl(dest:typed; l:int; src:typed):untyped =
  for i in 0..<nc:
    for j in 0..<nc:
      dest[i][j].re[l] = src[i][j].re
      dest[i][j].im[l] = src[i][j].im
proc vcopy(dest:var SColorMatrixV; l:int; src:SColorMatrix) =
  vcopyImpl(dest,l,src)
proc vcopy(dest:var DColorMatrixV; l:int; src:DColorMatrix) =
  vcopyImpl(dest,l,src)
proc vcopy(dest:var SColorMatrixV; l:int; src:DColorMatrix) =
  vcopyImpl(dest,l,src)
proc vcopy(dest:var DColorMatrixV; l:int; src:SColorMatrix) =
  vcopyImpl(dest,l,src)

proc read[T](r:var Reader, v:var openArray[ptr T]) =
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
      let vl = index mod simdLength(T)
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
      let vl = index mod simdLength(T)
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
  var rd = lo.newReader("l88.scidac")
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
      for c in 0..<nc:
        #echo g[i][x][c,c]
        tr += g[i][x][c,c]
  echo tr
  qexFinalize()


discard """
int
QDP_read_check(QDP_Reader *qdpr, QDP_String *md, int globaldata,
               void (* put)(char *buf, size_t index, int count, void *qfin),
	       struct QDP_IO_field *qf, int count, QIO_RecordInfo *cmp_info)
{
  QIO_RecordInfo *rec_info;
  QIO_String *qio_md = QIO_string_create();
  int status;
  
  iolat = qdpr->lat;
  rec_info = QIO_create_record_info(0, 0, 0, 0, "", "", 0, 0, 0, 0);

  status = QIO_read(qdpr->qior, rec_info, qio_md, put, qf->size*count,
                    qf->word_size, (void *)qf);
  QDP_string_set(md, QIO_string_ptr(qio_md));
  QIO_string_destroy(qio_md);
  
  /* Check for consistency */
  if(QIO_compare_record_info(rec_info, cmp_info)) {
      status = 1;
  }

  QIO_destroy_record_info(rec_info);

  return status;
}

static void
QDP$PC_vput_$ABBR(char *buf, size_t index, int count, void *qfin)
{
  struct QDP_IO_field *qf = qfin;
  $QDPPCTYPE **field = ($QDPPCTYPE **)(qf->data);
  #define nc qf->nc
  $QLAPCTYPE($NCVAR(*dest));
  $QLAPCTYPE($NCVAR(*src)) = (void *)buf;
  #undef nc
  int i;

  /* For the site specified by "index", move an array of "count" data
  from the read buffer to an array of fields */

  for(i=0; i<count; i++) {
      //dest = QDP$PC_expose_$ABBR(NC field[i] ) + index;
      dest = QDP_offset_data(field[i],index);
      QDPIO_put_$ABBR(NC, $P, $PC, dest, src+i, qf->nc, qf->ns);
      //QDP$PC_reset_$ABBR(NC field[i] );
  }
}

/* Internal factory function for global data */
static void
QDP$PC_vput_$QLAABBR(char *buf, size_t index, int count, void *qfin)
{
  struct QDP_IO_field *qf = qfin;
  #define nc qf->nc
  $QLAPCTYPE($NCVAR(*dest)) = (void *)(qf->data);
  $QLAPCTYPE($NCVAR(*src)) = (void *)buf;
  #undef nc
  int i;

  for(i=0; i<count; i++) {
      QDPIO_put_$ABBR(NC, $P, $PC, (dest+i), (src+i), qf->nc, qf->ns);
  }
}

int
QDP$PC_vread_$ABBR(QDP_Reader *qdpr, QDP_String *md, $QDPPCTYPE *field[],
                   int nv)
{
  int status;
  TGET;
  ONE {
    struct QDP_IO_field qf;
    QIO_RecordInfo *cmp_info;

    qf.data = (char *) field;
    qf.size = QDPIO_size_$ABBR($P, QDP_get_nc(field[0]), QLA_Ns);
    qf.nc = QDPIO_nc_$ABBR(QDP_get_nc(field[0]));
    qf.ns = QDPIO_ns_$ABBR(QLA_Ns);
    qf.word_size = WS;

    QDP_set_iolat(qdpr->lat);
    cmp_info = QIO_create_record_info(QIO_FIELD, 0, 0, 0,
	                              "$QDPPCTYPE", "$P", qf.nc,
				      qf.ns, qf.size, nv);
    
    for(int i=0; i<nv; i++) QDP_prepare_dest( &field[i]->dc );
    
    status = QDP_read_check(qdpr, md, QIO_FIELD, QDP$PC_vput_$ABBR, &qf,
                            nv, cmp_info);
    
    QIO_destroy_record_info(cmp_info);

    SHARE_SET(&status);
    TBARRIER;
  } else {
    int *p;
    TBARRIER;
    SHARE_GET(p);
    status = *p;
  }
  TBARRIER;

  return status;
}

int
QDP$PC_read_$ABBR(QDP_Reader *qdpr, QDP_String *md, $QDPPCTYPE *field)
{
  $QDPPCTYPE *temp[1];
  temp[0] = field;
  return QDP$PC_vread_$ABBR(qdpr, md, temp, 1);
}

int
QDP$PC_vread_$QLAABBR($NC QDP_Reader *qdpr, QDP_String *md,
                      $QLAPCTYPE($NCVAR(*array)), int n)
{
  int status;
  TGET;
  ONE {
    struct QDP_IO_field qf;
    QIO_RecordInfo *cmp_info;

    qf.data = (char *) array;
    qf.size = QDPIO_size_$ABBR($P, $QDP_NC, QLA_Ns);
    qf.nc = QDPIO_nc_$ABBR($QDP_NC);
    qf.ns = QDPIO_ns_$ABBR(QLA_Ns);
    qf.word_size = WS;

    QDP_set_iolat(qdpr->lat);
    cmp_info = QIO_create_record_info(QIO_GLOBAL, 0, 0, 0,
	                              "$QDPPCTYPE", "$P", qf.nc,
				      qf.ns, qf.size, n);
    
    status = QDP_read_check(qdpr, md, QIO_GLOBAL, QDP$PC_vput_$QLAABBR, &qf,
                            n, cmp_info);

    QIO_destroy_record_info(cmp_info);
    
    SHARE_SET(&status);
    TBARRIER;
  } else {
    int *p;
    TBARRIER;
    SHARE_GET(p);
    status = *p;
  }
  TBARRIER;

  return status;
}
"""
