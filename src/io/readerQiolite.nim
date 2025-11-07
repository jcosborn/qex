import base
import layout
import strutils, strformat
import macros
import field
import os
import scidacio

var verb = 0

proc getFileLattice*(fn: string): seq[int] =
  if not fileExists(fn):
    warn "file not found: ", fn
    return
  var sr = newScidacReader(fn)
  result = sr.lattice
  sr.close

type Reader*[V:static[int]] = ref object
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
  sr*: ScidacReader
  atNextRecord*: bool

proc IOverb*(level: int) =
  verb = level

proc newReader*[V: static[int]](l: Layout[V]; fn: string): Reader[V] =
  if not fileExists(fn):
    warn "file not found: ", fn
    return
  var sr = newScidacReader(fn, verb)
  result.new
  result.layout = l
  result.latsize = sr.lattice
  result.status = 0
  result.fileName = fn
  result.fileMetadata = sr.fileMd
  result.verb = verb
  result.sr = sr
  result.atNextRecord = true

proc close*(rd: var Reader) =
  rd.sr.close

proc nextRecord*(rd: var Reader) =
  rd.sr.nextRecord
  rd.atNextRecord = true

proc recordMetadata*(rd: var Reader): string =
  rd.sr.recordMd

#template recGet(f,t:untyped):untyped =
#  proc f*(r: var Reader): auto =
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
proc recordDate*(rd: var Reader): string =
  rd.sr.record.date

proc hyperindex(x: seq[cint]; subl,offs: seq[int]): int =
  let n = x.len
  for i in countdown(n-1,0):
    result = result*subl[i] + (x[i]-offs[i])

import qioInternal

proc put[T](buf: cstring; index: int, count: int; arg: openArray[ptr T]) =
  type srcT = cArray[IOtype(T)]
  type destT1 = cArray[T]
  type destT = cArray[ptr destT1]
  let src = cast[ptr srcT](buf)
  let dest = cast[ptr destT](arg)
  let vi = index div simdLength(T)
  let vl = index mod simdLength(T)
  #let vlm = 1 shl vl
  for i in 0..<count:
    #masked(dest[i][vi], vlm) := src[i]
    #mindexed(dest[i][vi], asSimd(vl)) := src[i]
    dest[][i][][vi][asSimd(vl)] = src[][i]

proc putP[T](buf: cstring; index: int, count: int; arg: openArray[ptr T]) =
  type srcT = cArray[IOtypeP(T)]
  type destT1 = cArray[T]
  type destT = cArray[ptr destT1]
  let src = cast[ptr srcT](buf)
  let dest = cast[ptr destT](arg)
  let vi = index div simdLength(T)
  let vl = index mod simdLength(T)
  #let vlm = 1 shl vl
  for i in 0..<count:
    #masked(dest[i][vi], vlm) := src[i]
    #mindexed(dest[i][vi], asSimd(vl)) := src[i]
    dest[i][][vi][asSimd(vl)] = src[i]

template `&`[T](x: seq[T]): untyped = cast[ptr T](unsafeaddr x[0])
template `+`(x: ptr char, i: SomeInteger): untyped =
  cast[ptr char](cast[uint](x) + uint(i))

proc read[T](rd: var Reader, v: openArray[ptr T]) =
  let t0 = getMonoTime()
  if not rd.atNextRecord: rd.nextRecord
  template r: untyped = rd.sr.record
  let nsites = rd.layout.nsites
  let objcount = r.datacount
  let datum_size = objcount * r.typesize
  #echo "datum_size: ", datum_size, "  objcount: ", objcount
  let nbytes = nsites * datum_size
  let ioprec = if rd.sr.record.precision=="F": 1 else: 2
  let tprec = if sizeof(numberType(T)) == 4: 1 else: 2
  let ioelem = r.typesize / (4*ioprec)
  let telem = sizeof(IOtype(T)) / (4*tprec)
  #let telem = numNumbers IOtype(T)
  if ioelem != telem:
    qexError &"field data elements {telem} != reader data elements {ioelem}"
  let nd = rd.layout.nDim
  var x = newSeq[cint](nd)
  let this_node = rd.layout.myrank
  var sublattice = newSeq[int](nd)
  var hypermin = rd.sr.lattice
  var hypermax = newSeq[int](nd)
  for i in 0 ..< nsites:
    rd.layout.coord(x, this_node, i)
    #echo i, " ", x
    for j in 0..<nd:
      hypermin[j] = min(hypermin[j], x[j])
      hypermax[j] = max(hypermax[j], x[j])
  #echo hypermin, " ", hypermax
  for j in 0..<nd:
    sublattice[j] = hypermax[j] - hypermin[j] + 1
  var buf = create(char, nbytes)
  let t1 = getMonoTime()
  rd.sr.readBinary(buf, sublattice, hypermin)
  let t2 = getMonoTime()
  for i in countup(0'i32, nsites.int32-1):
    rd.layout.coord(x, this_node, i)
    let j = hyperindex(x, sublattice, hypermin)
    let tbuf = cast[cstring](buf + j*datum_size.int)
    if ioprec==tprec:
      put(tbuf, i, objcount, v)
    else:
      putP(tbuf, i, objcount, v)
  dealloc(buf)
  let t3 = getMonoTime()
  rd.sr.finishReadBinary
  rd.atNextRecord = false
  let t4 = getMonoTime()
  if verb > 0:
    echo &"read seconds layout: {t1-t0} read: {t2-t1} put: {t3-t2} finish: {t4-t3}"

proc check(rd: Reader, l: Layout) =
  var err = false
  let nd = l.nDim
  if rd.sr.lattice.len != nd:
    err = true
  else:
    for i in 0..<nd:
      if rd.sr.lattice[i] != l.physGeom[i]:
        err = true
  if err:
    qexError &"field lattice {l.physGeom} != reader lattice {rd.sr.lattice}"

proc read*[V:static[int],T](rd: var Reader[V]; v: openArray[Field[V,T]]) =
  let nv = v.len
  if rd.sr.record.datacount != nv:
    qexError &"field count {nv} != reader count {rd.sr.record.datacount}"
  for i in 0..<nv:
    check(rd, v[i].l)
  var f = newSeq[type(v[0][0].addr)](nv)
  for i in 0..<nv: f[i] = v[i][0].addr
  rd.read(f)

proc read*(rd: var Reader; v: Field) =
  if rd.sr.record.datacount != 1:
    qexError &"field count 1 != reader count {rd.sr.record.datacount}"
  check(rd, v.l)
  var f = @[ v[0].addr ]
  rd.read(f)

when isMainModule:
  qexInit()
  echo "rank ", myRank, "/", nRanks
  var lat = [8,8,8,8]
  var lo = newLayout(lat)
  var g:array[4,type(lo.ColorMatrix())]
  for i in 0..<4: g[i] = lo.ColorMatrix()
  var fn = stringParam("fn", "testlat0.bin")
  if not fileExists(fn):
    echo "gauge file not found: ", fn
    qexExit()
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
  rd.read(g)
  #rd.nextRecord()
  #echo rd.status
  #echo rd.recordMetadata
  rd.close()
  var unit0: evalType(g[0][0].norm2)
  var tr0: evalType(g[0][0][0,0])
  let nc = g[0][0].ncols
  for i in 0..<4:
    for x in g[i]:
      let t = g[i][x].adj * g[i][x] - 1
      unit0 += t.norm2
      for c in 0..<nc:
        #echo g[i][x][c,c]
        tr0 += g[i][x][c,c]
  let s = 1.0/(nc * lo.nDim * lo.physVol)
  var unit = simdReduce unit0
  rankSum unit
  echo "unitary error: ", s*unit
  var tr = simdReduce tr0
  rankSum tr
  echo "normalized trace: ", s*tr
  qexFinalize()
