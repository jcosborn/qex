import os
import streams
import strUtils
import stdUtils
import fieldProxy
import memfiles

type
  MmapFieldObj[T] = object
    n: int
    m: Memfile
    p: ptr cArray[T]
  MmapField[T] = FieldProxy[MmapFieldObj[T]]

proc newFieldImpl(x: MmapFieldObj, y: typedesc): MmapField[y] =
  discard
proc newMmapField[T](fn: string, n: int): MmapField[T] =
  result[].n = n
  let bytes = n*sizeof(T)
  result[].m = memfiles.open(fn, mode=fmReadWrite, newFileSize=bytes)
  result[].p = cast[ptr cArray[T]](result[].m.mem)
proc close[T](x: var MmapField[T]) =
  memfiles.close(x[].m)

template `[]`(x: MmapFieldObj, y: typed): untyped =
  x.p[y]
template `[]=`(x: MmapFieldObj, y: typed, z: untyped): untyped =
  x.p[y] = z
template `len`*(x: MmapFieldObj): untyped = x.n
#iterator indices(x: FieldObj): int =
#  countup(0,99)
#iterator indices(x: FieldObj): int =
#  countup(0,99)
template indices(x: MmapFieldObj): untyped = 0..<x.n
fieldScalarOverloads(SomeNumber)


when isMainModule:
  proc product[T](x: openarray[T]): auto =
    result = x[0]
    for i in 1..<x.len: result *= x[i]
  var home = getHomeDir()
  var fn = home / "lqcd/milc/milc_qcd-7.8.0/binary_samples/lat.sample.l4444"
  var f = newFileStream(fn)
  var hdr: MilcHeader
  var n = f.readData(addr hdr, sizeof hdr)
  echo hdr.magicNumber
  for d in hdr.dims:
    echo d
  echo hdr.timeStamp
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
