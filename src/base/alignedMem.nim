import strUtils
import stdUtils
import system/ansi_c

type
  RawMemStats* = object
    allocated*: int
    used*: int
    maxUsed*: int
  RawMem* = object
    size*: int
    p*: ptr UncheckedArray[char]
  RawMemRef* = ref RawMem
  alignedMem*[T] = object
    len*: int
    align*: int
    stride*: int
    bytes*: int
    #mem*: ref cArray[char]
    mem*: RawMemRef
    data*: ptr cArray[T]

var rms: RawMemStats
var rmGcThreshold = 1024*1024*1024   # 1 GB

proc newRawMem*(size: int, zero=true): RawMem =
  if size>0:
    if rms.used>rmGcThreshold:
      #echo GC_getStatistics()
      GC_fullCollect()
      #echo GC_getStatistics()
    let p = if zero: c_calloc(1,csize_t(size)) else: c_malloc(csize_t(size))
    if not p.isNil:
      result.p = cast[type result.p](p)
      result.size = size
      rms.allocated += size
      rms.used += size
      rms.maxUsed = max(rms.maxUsed, rms.used)
proc free*(rm: var RawMem) =
  if rm.size>0:
    c_free(rm.p)
    rms.used -= rm.size
    rm.size = 0
    rm.p = nil
template data*(rm: RawMem): untyped = rm.p

proc getRawMemAllocated*(): int = rms.allocated
proc getRawMemUsed*(): int = rms.used
proc getRawMemMaxUsed*(): int = rms.maxUsed
proc getRawMemGcThreshold*(): int = rmGcThreshold
proc setRawMemGcThreshold*(t: int) = rmGcThreshold = t

proc freeRawMemRef*(rm: RawMemRef) =
  #echo "freeRawmemref"
  free(rm[])
proc newRawMemRef*(size: int): RawMemRef =
  #echo "newRawMemRef"
  result.new(freeRawMemRef)
  result[] = newRawMem(size)
proc newRawMemRefU*(size: int): RawMemRef =
  #echo "newRawMemRef"
  result.new(freeRawMemRef)
  result[] = newRawMem(size, zero=false)
template data*(rm: RawMemRef): untyped = rm[].p


proc unsafeNewU*[T](a: var ref T, size: Natural) =
  {.emit: "N_NIMCALL(void*, newObjNoInit)(TNimType* typ0, NI size0);".}
  {.emit: "#define newObj newObjNoInit".}
  unsafeNew(a, size)
  {.emit: "#undef newObj".}

proc ptrAlign[T](p:ptr T; a:int):ptr T =
  let x = cast[int](p)
  let a1 = a - 1
  let y = x + (a1-((x+a1) mod a))
  #echo x, ":", y
  result = cast[type(result)](y)

proc new*[T](t:var alignedMem[T], n:int, align:int=64) =
  t.len = n
  t.align = align
  t.stride = sizeof(T)
  t.bytes = t.len * t.stride + t.align
  #unsafeNew(t.mem, t.bytes)
  #t.data = ptrAlign(cast[ptr cArray[T]](t.mem[0].addr), align)
  t.mem = newRawMemRef(t.bytes)
  t.data = ptrAlign(cast[ptr UncheckedArray[T]](t.mem.data), align)
proc newU*[T](t:var alignedMem[T], n:int, align:int=64) =
  t.len = n
  t.align = align
  t.stride = sizeof(T)
  t.bytes = t.len * t.stride + t.align
  #unsafeNewU(t.mem, t.bytes)
  #t.data = ptrAlign(cast[ptr cArray[T]](t.mem[0].addr), align)
  t.mem = newRawMemRefU(t.bytes)
  t.data = ptrAlign(cast[ptr UncheckedArray[T]](t.mem.data), align)
proc newAlignedMem*[T](t:var alignedMem[T], n:int, align:int=64) =
  new(t, n, align)
proc newAlignedMem*[T](n:int, align:int=64): alignedMem[T] =
  newAlignedMem[T](result, n, align)
proc newAlignedMemU*[T](t:var alignedMem[T], n:int, align:int=64) =
  newU(t, n, align)
proc newAlignedMemU*[T](n:int, align:int=64): alignedMem[T] =
  newAlignedMemU[T](result, n, align)

template low*(s:alignedMem):untyped = 0
template high*(s:alignedMem):untyped = s.len-1
#proc `[]`*[T](s:alignedMem[T], i:SomeInteger):var T =
#  result = s.data[i]
#template `[]`*[T](s:alignedMem[T], i:SomeInteger):untyped = s.data[i]
template `[]`*[T](s: alignedMem[T], i:SomeInteger):untyped = s.data[i]
template `[]=`*[T](s:var alignedMem[T], i:SomeInteger, v:typed) =
  s.data[i] = v

when isMainModule:
  template stats =
    echo "allocated: ", getRawMemAllocated()
    echo "used: ", getRawMemUsed()
    echo "maxUsed: ", getRawMemMaxUsed()

  proc test =
    stats()
    var x: alignedMem[float]
    newAlignedMem(x, 10)
    let c0 = cast[int](x.mem.data)
    echo c0, " ", toHex(c0,8)
    let x0 = cast[int](x[0].addr)
    echo x0, " ", toHex(x0,8)
    for i in x.low..x.high:
      x[i] = float(i)
    for i in 0..<x.len:
      assert(x[i] == float(i))
    stats()

  test()
  stats()
  GC_fullCollect()
  stats()
