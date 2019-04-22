import parallelIo
import endians, tables, strutils

type
  ModFileHeader* = object
    magic*: string
    version*: int32
    userdata*: string
    mapstart*: int
  ModFileMap* = TableRef[string,int]

proc getPos*(m: ModFileMap, x: any): int =
  let l = sizeof(x)
  var s = newString(l)
  copyMem(s[0].addr, unsafeAddr(x), l)
  #echo s.toHex
  m[s]

proc getPos*(m: ModFileMap, a,b: int): int =
  let aa = toBigEndian(a.int32)
  let bb = toBigEndian(b.int32)
  m.getPos((aa,bb))

proc readString(pr: var ParallelReader): string =
  var bytes = pr.readBigInt32()
  result = newString(bytes)
  pr.readAll(result)

proc modReadHeader*(pr: var ParallelReader): ModFileHeader =
  result.magic = pr.readString()
  result.version = pr.readBigInt32()
  result.userdata = pr.readString()
  var dum = pr.readBigInt64()
  result.mapstart = pr.readBigInt64().int
  #echo result

proc modReadMap*(pr: var ParallelReader, mdstart: int): ModFileMap =
  result = newTable[string,int]()
  pr.seekSet(mdstart)
  let num = pr.readBigInt32()
  #echo num
  for i in 0..<num:
    let k = pr.readString()
    var dum = pr.readBigInt64()
    let v = pr.readBigInt64().int
    #echo k.toHex()
    #echo i, ": ", v
    result.add(k, v)
