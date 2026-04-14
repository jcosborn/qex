import tables, strformat, strutils, algorithm

# kinds
# temp: can be freed any time
# locked: can't be freed
# hostbacked: can copy back to host (need to check dirty?)
# otherwise: needs to be saved to temp

# keys
# host ptr
# string

type
  GpuMem* = object
    p*: pointer
    bytes*: int
    lastGet*: int
    useCount*: int  # counts times used in gpu kernel
    lastCopyIn*: int  # useCount of last copy to gpu
    lastCopyOut*: int  # useCount of last copy from gpu
    noReadCpu* {.bitsize:1.}: bool # won't read from CPU field after next kernel
    noReadGpu* {.bitsize:1.}: bool # won't read from GPU field in next kernel
    noWriteCpu* {.bitsize:1.}: bool # haven't written to CPU field since last sync
    noWriteGpu* {.bitsize:1.}: bool # won't write to GPU field in next kernel
    tag*: seq[string]
proc summary*(x: GpuMem): string =
  result = &"""bytes {x.bytes:10d}  lastGet {x.lastGet:5d}  useCount {x.useCount:5d}  tag {x.tag.join(":")}"""

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

template needsCopyIn*(x: GpuMem | ptr GpuMem): bool =
  let cpy = (not x.noReadGpu) and (not x.noWriteCpu)
  cpy

template needsCopyOut*(x: GpuMem | ptr GpuMem): bool =
  let cpy = (not x.noReadCpu) and (not x.noWriteGpu)
  cpy

template touch*(x: ptr GpuMem) =
  x.lastGet = kernelCallCount
  inc x.useCount

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

import backend/accelbase

proc freeGpuMem*(cpuPtr: pointer) =
  mixin gpuFree
  var x: GpuMem
  if gpuMemTable.pop(cpuPtr, x):
    gpuFree(x.p)

proc getGpuMem*(cpuPtr: pointer, gpuBytes: int): ptr GpuMem =
  result = addr gpuMemTable.mgetOrPut(cpuPtr)
  if result.p == nil:
    result.bytes = gpuBytes
    result.p = gpuMalloc(gpuBytes)
    result.tag = gpuMemTag
  else:
    if result.bytes != gpuBytes:
      echo "getGpuMem bytes mismatch result.bytes: ", result.bytes, "  gpuBytes: ", gpuBytes
      doAssert(result.bytes == gpuBytes)
  result.touch
  #echo result[]

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
  for v in gms:
    result &= "\n " & v.summary

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
