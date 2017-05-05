import os
import streams
import strUtils
import base
import times

const MILC_MAGIC_NUMBER = 20103

type
  MilcChecksum* = object
    sum29: uint32
    sum31: uint32
  MilcHeader* = object
    magicNumber*: int32         # Identifies file format
    dims*: array[4,int32]       # Full lattice dimensions
    timeStampChars*: array[64,char]  # Date and time stamp - used to
                                # check consistency between the
                                # ASCII header file and the
                                # lattice file
    order*: int32               # 0 means no coordinate list is
                                # attached and the values are in
                                # coordinate serial order.
                                # Nonzero means that a
                                # coordinate list is attached,
                                # specifying the order of values
    checksum*: MilcChecksum

proc initMilcHeader*(x: var MilcHeader, dims: openarray[SomeInteger],
                     timestamp: string = nil) =
  var ts = timestamp
  if ts.isNil:
    let t = getLocalTime getTime()
    ts = format(t , "ddd MMM dd HH:mm:ss yyyy")
  x.magicNumber = MILC_MAGIC_NUMBER
  for i in 0..<min(4,dims.len):
    x.dims[i] = dims[i]
  for i in 0..<min(64,ts.len):
    x.timestampChars[i] = ts[i]
  x.order = 0
proc initMilcHeader*(x: var MilcHeader, dims: openarray[SomeInteger]) =
  initMilcHeader(x, dims, nil)

proc timestamp*(x: MilcHeader): string =
  result = join(x.timeStampChars)

proc getMilcChecksum*(x: ptr carray[uint32], n: int, offset=0): MilcChecksum =
  template `^=`(x,y: untyped) = x = x xor y
  for i in 0..<n:
    let t = x[i]
    let k = offset + i
    let k29 = k mod 29
    let k31 = k mod 31
    result.sum29 ^= (t shl k29) + (t shr (32-k29))
    result.sum31 ^= (t shl k31) + (t shr (32-k31))
template getMilcChecksum*[T](x: openArray[T], offset=0): MilcChecksum =
  let n = (x.len*sizeof(x[0])) div sizeof(uint32)
  getMilcChecksum(cast[ptr carray[uint32]](addr x[0]), n, offset)

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
