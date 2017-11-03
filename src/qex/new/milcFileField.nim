import os
import streams
import strUtils
import qex/base
import fieldProxy
import memfiles
import qex/maths
import qex/maths/types
import layoutUtils
import milcIO

type
  MilcFileFieldObj[T] = object
    hdr*: MilcHeader
    nsites*: int
    m*: Memfile
    needsCksum*: bool
    p*: ptr cArray[T]
  MilcFileFieldObjRef[T] = ref MilcFileFieldObj[T]
  MilcFileField[T] = FieldProxy[MilcFileFieldObjRef[T]]
  MilcG* = object
    dir: array[4,MatrixArray[3,3,ComplexType[float32]]]

proc newFieldImpl(x: MilcFileFieldObj, y: typedesc): MilcFileField[y] =
  discard
proc closeMilcFileFieldObjRef[T](x: MilcFileFieldObjRef[T]) =
  if x.needsCksum:
    let n = 4 * 18 * x.nsites
    let c = getMilcChecksum(cast[ptr carray[uint32]](x.p), n)
    x.hdr.checksum = c
    (cast[ptr type(x.hdr)](x.m.mem))[] = x.hdr
  memfiles.close(x[].m)
proc close[T](x: var MilcFileField[T]) = closeMilcFileFieldObjRef(x[])
proc openMilcFileGaugeField*(fn: string, mode=fmRead): MilcFileField[MilcG] =
  if mode notin {fmRead,fmReadWrite}:
    echo "openMilcFileGaugeField mode must be fmRead or fmReadWrite"
    quit -1
  result[].new(closeMilcFileFieldObjRef[MilcG])
  result[].m = memfiles.open(fn, mode)
  result[].hdr = (cast[ptr type(result[].hdr)](result[].m.mem))[]
  result[].nsites = volume(result[].hdr.dims)
  result[].needsCksum = false
  echo result[].hdr
  let p0 = cast[ptr cArray[char]](result[].m.mem)
  result[].p = cast[ptr cArray[MilcG]](addr p0[sizeof(result[].hdr)])
proc newMilcFileGaugeField*(fn: string, dims: openarray[SomeInteger],
                            timestamp: string = nil): MilcFileField[MilcG] =
  result[].new(closeMilcFileFieldObjRef[MilcG])
  let nsites = volume(dims)
  result[].nsites = nsites
  let size = nsites * sizeof(MilcG) + sizeof(result[].hdr)
  result[].m = memfiles.open(fn, fmWrite, newFileSize=size)
  result[].needsCksum = true
  result[].hdr.initMilcHeader(dims, timestamp)
  (cast[ptr type(result[].hdr)](result[].m.mem))[] = result[].hdr
  let p0 = cast[ptr cArray[char]](result[].m.mem)
  result[].p = cast[ptr cArray[MilcG]](addr p0[sizeof(result[].hdr)])

template `[]`(x: MilcFileFieldObj, y: typed): untyped =
  x.p[y]
template `[]=`(x: MilcFileFieldObj, y: typed, z: untyped): untyped =
  x.p[y] = z
template `len`*(x: MilcFileFieldObj): untyped = x.nsites
template `[]`(x: MilcFileFieldObjRef, y: typed): untyped =
  x.p[y]
template `[]=`(x: MilcFileFieldObjRef, y: typed, z: untyped): untyped =
  x.p[y] = z
template `len`*(x: MilcFileFieldObjRef): untyped = x.nsites
#iterator indices(x: FieldObj): int =
#  countup(0,99)
#iterator indices(x: FieldObj): int =
#  countup(0,99)
template indices*(x: MilcFileFieldObj): untyped = 0..<x.nsites
template indices*(x: MilcFileFieldObjRef): untyped = 0..<x.nsites
#fieldScalarOverloads(SomeNumber)

template `:=`*(x,y: MilcG) =
  x.dir[0] := y.dir[0]
  x.dir[1] := y.dir[1]
  x.dir[2] := y.dir[2]
  x.dir[3] := y.dir[3]

when isMainModule:
  var home = getHomeDir()
  #var fn = home / "lqcd/milc/milc_qcd-7.8.0/binary_samples/lat.sample.l4444"
  var fn = home / "lqcd/milc/milc_qcd-git-201510211317/binary_samples/lat.sample.l4444"
  var f = openMilcFileGaugeField(fn)
  echo f.len
  echo f[0]
  echo f[0].dir[0]
  var s = [0.0,0,0,0]
  for i in f[].indices:
    for j in 0..3:
      let m = f[i].dir[j]
      let t = (m.adj * m - 1).norm2
      #echo t
      s[j] += t
  echo s

  let ts = f[].hdr.timestamp
  var f2 = newMilcFileGaugeField("test.lat", f[].hdr.dims, ts)
  f2 := f

  f2.close

  #[
  var hdr: MilcHeader
  var n = f.readData(addr hdr, sizeof hdr)
  echo hdr.magicNumber
  for d in hdr.dims:
    echo d
  echo hdr.timestamp
  echo hdr.order
  echo hdr.checksum
  var nsites = product(hdr.dims)
  echo nsites
  var nreals = 4*18*nsites
  echo nreals
  var dat = newSeq[float32](nreals)
  n = f.readData(addr dat[0], dat.len*sizeof dat[0])
  #let ck = getMilcChecksum(cast[ptr carray[uint32]](addr dat[0]), nreals, 0)
  let ck = getMilcChecksum(dat)
  echo ck
  #for d in dat:
  #  echo d
  ]#
