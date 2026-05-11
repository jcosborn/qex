import tables, strformat, strutils, algorithm
import backend/accelbase

# kinds
# temp: can be freed any time
# locked: can't be freed
# hostbacked: can copy back to host (need to check dirty?)
# otherwise: needs to be saved to temp

# keys
# host ptr
# string

type
  gmFlags* = enum
    gmValidFlags, gmValidBytes,
    gmCpuRead, gmCpuWrite, gmCpuWriteOnce,
    gmGpuRead, gmGpuWrite, gmGpuWriteOnce
const defaultFlags = {gmValidFlags, gmCpuRead, gmCpuWrite, gmGpuRead, gmGpuWrite}

type
  GpuMem* = object
    p*: pointer
    bytes*: int
    lastGet*: int  # kernelCallCount of last access
    getCount*: int  # counts times table entry was accessed
    kernelCount*: int  # counts number of unique kernels counts table entry was accessed
    lastCopyIn*: int  # kernelCallCount of last copy to gpu
    lastCopyOut*: int  # kernelCallCount of last copy from gpu
    copyInCount*: int
    copyOutCount*: int
    flags*: set[gmFlags]
    #noReadCpu* {.bitsize:1.}: bool # won't read from CPU field after next kernel
    #noReadGpu* {.bitsize:1.}: bool # won't read from GPU field in next kernel
    #noWriteCpu* {.bitsize:1.}: bool # haven't written to CPU field since last sync
    #noWriteGpu* {.bitsize:1.}: bool # won't write to GPU field in next kernel
    tag*: seq[string]
proc summary*(x: GpuMem): string =
  #result = &"""bytes {x.bytes:10d}  lastGet {x.lastGet:5d}  gets {x.getCount:5d}  tag {x.tag.join(":")}"""
  #result = &"""bytes {x.bytes:10d}  kernels {x.kernelCount:5d}  cpyIn {x.copyInCount:5d}  cpyOut {x.copyOutCount:5d}  tag {x.tag.join(":")}"""
  result = &"""{x.bytes:10d} bytes {x.kernelCount:5d} kernels {x.copyInCount:5d} cpyIn {x.copyOutCount:5d} cpyOut {x.tag.join(":")}"""

const statsHeader = "     bytes  kerns  cpyIn cpyOut   tag"
proc stats(x: GpuMem): string =
  result = &"""{x.bytes:10d} {x.kernelCount:6d} {x.copyInCount:6d} {x.copyOutCount:6d} {x.tag.join(":")}"""

#[
template cpuUnused*(x: GpuMem | ptr GpuMem) =
  x.noReadCpu = true
  x.noWriteCpu = true
template gpuUnused*(x: GpuMem | ptr GpuMem) =
  x.noReadGpu = true
  x.noWriteGpu = true
template cpuReadOnly*(x: GpuMem | ptr GpuMem) =
  x.noReadCpu = false
  x.noWriteCpu = true
template gpuReadOnly*(x: GpuMem | ptr GpuMem) =
  x.noReadGpu = false
  x.noWriteGpu = true
template cpuWriteOnly*(x: GpuMem | ptr GpuMem) =
  x.noReadCpu = true
  x.noWriteCpu = false
template gpuWriteOnly*(x: GpuMem | ptr GpuMem) =
  x.noReadGpu = true
  x.noWriteGpu = false
template cpuReadWrite*(x: GpuMem | ptr GpuMem) =
  x.noReadCpu = false
  x.noWriteCpu = false
template gpuReadWrite*(x: GpuMem | ptr GpuMem) =
  x.noReadGpu = false
  x.noWriteGpu = false
]#

template isNew*(x: GpuMem | ptr GpuMem): bool =
  x.getCount == 1
  #x.p == nil

proc needsCopyIn*(x: GpuMem | ptr GpuMem): bool =
  if x.lastCopyIn < kernelCallCount:
    #result= (not x.noReadGpu) and (not x.noWriteCpu)
    result = (gmGpuRead in x.flags) and
      ((gmCpuWrite in x.flags) or (gmCpuWriteOnce in x.flags))

proc needsCopyOut*(x: GpuMem | ptr GpuMem): bool =
  if x.lastCopyOut < kernelCallCount:
    #result = (not x.noReadCpu) and (not x.noWriteGpu)
    result = {gmCpuRead, gmGpuWrite} <= x.flags

proc touch*(x: ptr GpuMem) =
  inc x.getCount
  if x.lastGet < kernelCallCount:
    x.lastGet = kernelCallCount
    inc x.kernelCount

template wasCopiedIn*(x: ptr GpuMem) =
  x.lastCopyIn = kernelCallCount
  inc x.copyInCount
  if gmCpuWriteOnce in x.flags:
    x.flags.excl {gmCpuWriteOnce,gmCpuWrite}

template wasCopiedOut*(x: ptr GpuMem) =
  x.lastCopyOut = kernelCallCount
  inc x.copyOutCount

proc copyIn*(x: GpuMem | ptr GpuMem, p: pointer) =
  if x.needsCopyIn:
    gpuMemCpyToGpu(x.p, p, x.bytes)
    x.wasCopiedIn

proc copyOut*(x: GpuMem | ptr GpuMem, p: pointer) =
  if x.needsCopyOut:
    gpuMemCpyToCpu(p, x.p, x.bytes)
    x.wasCopiedOut

#proc newGpuMem*(gpuBytes: int): GpuMem =
#  result.p = gpuMalloc(gpuBytes)
#  result.bytes = gpuBytes
#  result.count

var gpuMemTable = newTable[pointer, GpuMem]()
#var gpuMemNils = newSeq[GpuMem]()
var gpuMemTag = newSeq[string](0)
#proc setGpuMemTag*(s: string) =
#  gpuMemTag = s
proc pushGpuMemTag*(s: string) =
  gpuMemTag.add s
proc popGpuMemTag*() =
  gpuMemTag.setLen(gpuMemTag.len - 1)

proc freeGpuMem*(cpuPtr: pointer) =
  mixin gpuFree
  var x: GpuMem
  if gpuMemTable.pop(cpuPtr, x):
    gpuFree(x.p)

proc getGpuMem*(cpuPtr: pointer): ptr GpuMem =
  result = addr gpuMemTable[cpuPtr]
  result.touch
  #echo result[]

proc getGpuMemDef(cpuPtr: pointer): ptr GpuMem =
  result = addr gpuMemTable.mgetOrPut(cpuPtr)
  if gmValidFlags notin result.flags:
    result.flags = defaultFlags

proc gpuMemFlagsExcl*(cpuPtr: pointer, flags: set[gmFlags]) =
  var pgm = getGpuMemDef(cpuPtr)
  pgm.flags.excl flags

proc gpuMemFlagsIncl*(cpuPtr: pointer, flags: set[gmFlags]) =
  var pgm = getGpuMemDef(cpuPtr)
  pgm.flags.incl flags

var gb = 0
proc getGpuMemImpl*(cpuPtr: pointer, gpuBytes: int): ptr GpuMem =
  result = getGpuMemDef(cpuPtr)
  #echo gpuMemTag
  if gmValidBytes notin result.flags:
    result.flags.incl gmValidBytes
    result.bytes = gpuBytes
    result.tag = gpuMemTag
  else:
    if result.bytes != gpuBytes:
      gb = result.bytes
      return nil
  if result.p == nil:
    result.p = gpuMalloc(gpuBytes)
  result.touch
  #echo result[]
template getGpuMem*(cpuPtr: pointer, gpuBytes: int): ptr GpuMem =
  let p = getGpuMemImpl(cpuPtr, gpuBytes)
  if p == nil:
    echo "getGpuMem bytes mismatch result.bytes: ", gb, "  gpuBytes: ", gpuBytes
    echo instantiationInfo()
    doAssert(gb == gpuBytes)
  p

proc dumpGpuMem*(): string =
  var gms = newSeq[GpuMem](0)
  var mem = 0
  for v in gpuMemTable.values:
    mem += v.bytes
    gms.add v
  gms.sort do (x, y: GpuMem) -> int:
    result = cmp(x.tag.join, y.tag.join)
    if result == 0:
      result = cmp(x.bytes, y.bytes)
  result = "GpuMem items: " & $gpuMemTable.len & "  bytes: " & ($mem).insertSep(',')
  result &= "\n " & statsHeader
  for v in gms:
    result &= "\n " & v.stats

proc freeAllGpuMem*() =
  mixin gpuFree
  for v in gpuMemTable.values:
    gpuFree(v.p)
  gpuMemTable.clear

#[
proc gpuMem*(gpuBytes: int, cpuPtr: pointer = nil): pointer =
  if cpuPtr == nil:
    var t: GpuMem
    t.p = gpuMalloc(gpuBytes)
    t.bytes = gpuBytes
    t.lastGet = 
    gpuMemNils.add 
    let t = gpuMemGet(cpuPtr, gpuBytes)
    result = t.p
]#

# setGpuMemTag(string)
# store string in gpumem
