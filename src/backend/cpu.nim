import base/threading
import backend/expr
import macros

#type
#  Buffer*[T] = object
#    data: ptr UncheckedArray[T]
#    n: int

#proc init[T](b: Buffer[T], n: int) =
#  let p = allocShared(n*sizeof(T))
#  b.data = cast[type b.data](p)
#  b.n = n

#proc free[T](b: Buffer[T]) =
#  free b.data

#template onGpu*(x:untyped) =
#  threads:
#    template getThreadNum():auto = threadNum
#    template getNumThreads():auto = numThreads
#    x

template gpuMalloc*(size:SomeInteger):pointer = alloc(size)
template gpuFree*(device_ptr:pointer) = dealloc(device_ptr)
proc gpuMemCpyToCPU*(dst: pointer, src: pointer; length: SomeInteger) =
  copymem(dst, src, length)
proc gpuMemCpyToGPU*(dst: pointer, src: pointer; length: SomeInteger) =
  copymem(dst, src, length)

template getNumThreads*:auto = numThreads
template getThreadNum*:auto = threadNum


proc genCpuPrepare(n:seq[NimNode]):NimNode =
  mixin toGpu
  template r(x,v:untyped):untyped =
    var v = toGpu(x)
  result = newstmtlist()
  for c in n:
    result.add getast r(c[0],c[1])

proc genCpuFinalize(n:seq[NimNode]):NimNode =
  template r(x,v:untyped):untyped =
    fromGpu(x,v)
  result = newstmtlist()
  for c in n:
    result.add getast r(c[0],c[1])

macro onGpu*(body: untyped): auto =
  proc deref(x,g,i:NimNode):auto = newCall("getGpu",x,g)
  template target(cpuPrepare, cpuFinalize, body: untyped) =
    mixin toGpu, getGpu, fromGpu
    proc gpuProc {.gensym.} =
      cpuPrepare  # a let section declare and save device pointers
      threads:
        body
      cpuFinalize
    gpuProc()
  let
    v = prepareVars(body, deref)  # gather gpu pointers in symbols, body is changed accordingly
    cpuPrepare = genCpuPrepare v
    cpuFinalize = genCpuFinalize v
  result = getast(target(cpuPrepare, cpuFinalize, body))
  echo result.repr

template onGpu*(n,x:untyped) = onGpu(x)
template onGpu*(n,t,x:untyped) = onGpu(x)

