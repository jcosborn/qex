import base
import layout
import strutils, strformat
import field
import os
import scidacio
import qioInternal

var verb = 0

type Writer*[V:static[int]] = ref object
  layout*: Layout[V]
  latsize*: seq[int]
  status*: int
  fileName*: string
  fileMetadata*: string
  #recordInfoValid:bool
  #recordMd:string
  #nioranklist*: seq[int32]
  #iogeom*: seq[int32]
  #iorank*: seq[int32]
  verb*: int
  sw*: ScidacWriter
  atNextRecord*: bool

proc writeVerb*(level: int) =
  verb = level

proc newWriter*[V: static[int]](l: Layout[V]; fn,md: string): Writer[V] =
  #if fileExists(fn):
  #  warn "overwriting existing file: ", fn
  var sw = newScidacWriter(fn, l.physGeom, md, verb=verb)
  result.new
  result.layout = l
  result.latsize = l.physGeom
  result.status = 0
  result.fileName = fn
  result.fileMetadata = md
  result.verb = verb
  result.sw = sw
  result.atNextRecord = true

proc close*(wr: var Writer) =
  wr.sw.close

#template recGet(f,t:untyped):untyped =
#  proc f*(r: var Writer): auto =
    #if not r.recordInfoValid: r.readRecordInfo()
    #let x = `"QIO_get_" f`(r.recordInfo.addr)
    #return t(x)
#recGet(record_date, `$`)
#recGet(datatype, `$`)
#recGet(precision, `$`)
#recGet(colors, int)
#recGet(spins, int)
#recGet(typesize, int)
#recGet(datacount, int)
#proc recordDate*(rd: var Writer): string =
#  rd.sr.record.date

proc hyperindex(x: seq[cint]; subl,offs: seq[int]): int =
  let n = x.len
  for i in countdown(n-1,0):
    result = result*subl[i] + (x[i]-offs[i])

proc get[T](buf: cstring; index: int; count: int; arg: openArray[ptr T]) =
  type destT = cArray[IOtype(T)]
  type srcT1 = cArray[T]
  type srcT = cArray[ptr srcT1]
  let dest = cast[ptr destT](buf)
  let src = cast[ptr srcT](arg)
  let vi = index div simdLength(T)
  let vl = index mod simdLength(T)
  for i in 0..<count:
    dest[][i] := src[][i][][vi][asSimd(vl)]

proc getP[T](buf: cstring; index: int; count: int; arg: openArray[ptr T]) =
  type destT = cArray[IOtypeP(T)]
  type srcT1 = cArray[T]
  type srcT = cArray[ptr srcT1]
  let dest = cast[ptr destT](buf)
  let src = cast[ptr srcT](arg)
  let vi = index div simdLength(T)
  let vl = index mod simdLength(T)
  for i in 0..<count:
    dest[i] := src[][i][][vi][asSimd(vl)]

template `&`[T](x: seq[T]): untyped = cast[ptr T](unsafeaddr x[0])
template `+`(x: ptr char, i: SomeInteger): untyped =
  cast[ptr char](cast[ByteAddress](x) + ByteAddress(i))

proc write[T](wr: var Writer, v: var openArray[ptr T], lat: openArray[int],
              md="", ps="") =
  let t0 = getMonoTime()
  let tprec = if sizeof(numberType(T)) == 4: 1 else: 2
  let ioprec = case ps:
                 of "": tprec
                 of "F": 1
                 of "D": 2
                 else: 0
  if ioprec == 0:
    echo "ERROR: Unknown precision ps: ", ps
    qexExit 1
  let objcount = v.len
  let nsites = wr.layout.nsites
  let iotypebytes = (sizeof(IOtype(T)) div tprec) * ioprec
  let iositebytes = objcount * iotypebytes
  let iobytes = nsites * iositebytes
  var nc = v[0][].getNc.cint
  var ns = v[0][].getNs.cint
  let nd = lat.len
  var x = newSeq[cint](nd)
  let this_node = wr.layout.myrank
  var sublattice = newSeq[int](nd)
  var hypermin = @lat
  var hypermax = newSeq[int](nd)
  for i in 0 ..< nsites:
    wr.layout.coord(x, this_node, i)
    #echo i, " ", x
    for j in 0..<nd:
      hypermin[j] = min(hypermin[j], x[j])
      hypermax[j] = max(hypermax[j], x[j])
  #echo hypermin, " ", hypermax
  for j in 0..<nd:
    sublattice[j] = hypermax[j] - hypermin[j] + 1
  wr.sw.setRecord()
  wr.sw.record.precision = (if ioprec == 1: "F" else: "D")
  wr.sw.record.colors = nc
  wr.sw.record.spins = ns  # FIXME: no spin == -1?
  wr.sw.record.datacount = objcount
  wr.sw.record.typesize = iotypebytes
  if ioprec==tprec:
    #var datatype = "QDP_" & type(v[0][]).IOtype.name
    var datatype = type(v[0][]).IOtype.IOname
    wr.sw.record.datatype = datatype
  else:
    #var datatype = "QDP_" & type(v[0][]).IOtypeP.name
    var datatype = type(v[0][]).IOtypeP.IOname
    wr.sw.record.datatype = datatype
  var buf = create(char, iobytes)

  let t1 = getMonoTime()
  for i in countup(0'i32, nsites.int32-1):
    wr.layout.coord(x, this_node, i)
    let j = hyperindex(x, sublattice, hypermin)
    let tbuf = cast[cstring](buf + j*iositebytes)
    if ioprec==tprec:
      get(tbuf, i, objcount, v)
    else:
      getP(tbuf, i, objcount, v)

  let t2 = getMonoTime()
  wr.sw.initWriteBinary(md)
  wr.sw.writeBinary(buf, sublattice, hypermin)
  dealloc(buf)

  let t3 = getMonoTime()
  wr.sw.finishWriteBinary()
  let t4 = getMonoTime()
  if wr.verb > 0:
    echo &"write seconds layout: {t1-t0} get: {t2-t1} write: {t3-t2} finish: {t4-t3}"

proc write*[V:static[int],T](wr: var Writer[V]; v: openArray[Field[V,T]];
                             md = ""; prec = "") =
  let nv = v.len
  var f = newSeq[type(v[0][0].addr)](nv)
  for i in 0..<nv: f[i] = v[i][0].addr
  wr.write(f, v[0].l.physGeom, md, prec)

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
    #echo rd.datatype
    #echo rd.precision
    #echo rd.colors
    #echo rd.spins
    #echo rd.typesize
    #echo rd.datacount
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
