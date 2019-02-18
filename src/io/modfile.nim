import parallelIo
import endians

type
  ModFileHeader* = object
    magic*: string
    metadata*: string
    datastart*: int
    dataend*: int

proc modReadHeader*(pr: var ParallelReader): ModFileHeader =
  var bytes1 = pr.readBigInt32()
  echo bytes1
  var magic = newString(bytes1)
  pr.readAll(magic)
  #echo magic
  var dum1 = pr.readBigInt32()
  var bytes2 = pr.readBigInt32()
  echo dum1, " ", bytes2
  var md = newString(bytes2)
  pr.readAll(md)
  #echo md
  var dum2 = pr.readBigInt32()
  var dum3 = pr.readBigInt32()
  var dum4 = pr.readBigInt32()
  var bytes3 = pr.readBigInt32()
  echo dum2, " ", dum3, " ", dum4, " ", bytes3
  result.magic = magic
  result.metadata = md
  result.datastart = 4 + bytes1 + 8 + bytes2 + 16
  result.dataend = bytes3
  echo result
