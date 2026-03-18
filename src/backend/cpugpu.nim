import base/threading

type
  CpuGpu*[C,G] = object
    cpu*: C
    gpu*: G
    noReadCpu* {.bitsize:1.}: bool # won't read from CPU field after next kernel
    noReadGpu* {.bitsize:1.}: bool # won't read from GPU field in next kernel
    noWriteCpu* {.bitsize:1.}: bool # haven't written to CPU field since last sync
    noWriteGpu* {.bitsize:1.}: bool # won't write to GPU field in next kernel
    useCount*: int  # counts times used in gpu kernel
    lastCopyIn*: int  # useCount of last copy to gpu
    lastCopyOut*: int  # useCount of last copy from gpu
  # toGpu is called before kernel if it will be read on the gpu and was written on cpu since last sync
  # fromGpu is called after kernel if it was written on gpu and will be read on cpu later

template cpuUnused*(x: CpuGpu) =
  x.noReadCpu = true
  x.noWriteCpu = true
template gpuUnused*(x: CpuGpu) =
  x.noReadGpu = true
  x.noWriteGpu = true
template cpuReadOnly*(x: CpuGpu) =
  x.noReadCpu = false
  x.noWriteCpu = true
template gpuReadOnly*(x: CpuGpu) =
  x.noReadGpu = false
  x.noWriteGpu = true
template cpuWriteOnly*(x: CpuGpu) =
  x.noReadCpu = true
  x.noWriteCpu = false
template gpuWriteOnly*(x: CpuGpu) =
  x.noReadGpu = true
  x.noWriteGpu = false
template cpuReadWrite*(x: CpuGpu) =
  x.noReadCpu = false
  x.noWriteCpu = false
template gpuReadWrite*(x: CpuGpu) =
  x.noReadGpu = false
  x.noWriteGpu = false

proc toGpu*(x: var CpuGpu): auto {.discardable.} =
  mixin toGpu
  # copy if gpuR and cpuW
  let cpy = (not x.noReadGpu) and (not x.noWriteCpu)
  threadSingle:
    inc x.useCount
    if cpy: x.lastCopyIn = x.useCount
  x.gpu.toGpu(x, cpy)
  x.gpu

template getGpu*[C,G](x: CpuGpu[C,G], g: G): auto =
  getGpu(x, g)

proc fromGpu*[C,G](x: var CpuGpu[C,G], g: G) =
  mixin fromGpu
  # copy if cpuR and gpuW
  let cpy = (not x.noReadCpu) and (not x.noWriteGpu)
  threadSingle:
    if cpy: x.lastCopyOut = x.useCount
  x.fromGpu(g, cpy)

proc fromGpu*[C,G](x: var CpuGpu[C,G]) =
  # copy if cpuR and gpuW
  let cpy = (not x.noReadCpu) and (not x.noWriteGpu)
  threadSingle:
    if cpy: x.lastCopyOut = x.useCount
  x.fromGpu(x.gpu, cpy)

template wasCopiedIn*(x: CpuGpu): bool =
  x.lastCopyIn == x.useCount
template wasCopiedOut*(x: CpuGpu): bool =
  x.lastCopyOut == x.useCount

# cpuSyncRead
# cpuWasWritten
# gpuSyncRead
# gpuWasWritten

when isMainModule:
  import backend/accel

  type
    GpuArray[N:static int, T] = distinct array[N,T]
  template `[]`*[N:static int, T](x: GpuArray[N,T]): auto = (array[N,T])x
  template `[]`*(x: GpuArray, i: typed): auto = x[][i]
  template low(x: GpuArray): auto = x[].low
  template high(x: GpuArray): auto = x[].high

  type
    cgarray[N:static int, T] = CpuGpu[array[N,T], ptr GpuArray[N,T]]
  template `[]`*(x: cgarray, i: SomeNumber): auto = x.cpu[i]
  template `[]=`*(x: cgarray, i: SomeNumber, y: auto): auto =
    x.cpu[i] = y
  template toGpu*(g: var ptr GpuArray, x: CpuGpu, cpy: bool) =
    if g.isNil:
      when useGPU:
        g.gpuMalloc
      else:
        g = (typeof g)(addr x.cpu)
    when useGPU:
      if cpy:
        gpuMemCpyToGPU(g, addr x.cpu, sizeof(g[]))
  template getGpu*(x: cgarray, g: ptr GpuArray): auto = g[]
  template fromGpu*(x: var cgarray, g: ptr GpuArray, cpy: bool) =
    when useGPU:
      if cpy: gpuMemCpyToCPU(addr x.cpu[0], g, sizeof(x.cpu))

  const N = 1024
  var x: cgarray[N, float]
  var y: cgarray[N, float]
  var z: cgarray[N, float]
  for i in 0..<N:
    x[i] = 1.0 + i.float
    y[i] = 2.0 + i.float
  #x.gpuReadOnly
  #y.gpuReadOnly
  #z.gpuWriteOnly
  onGpu:
    let s = getNumThreads()
    var i = getThreadNum()
    while i < N:
      z[i] = x[i] * y[i]
      i += s
  echo "z[N-1]: ", z[N-1]
