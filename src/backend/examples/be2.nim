import backend/accel

type
  CpuGpu*[C,G] = object
    cpu*: C
    gpu*: G
  cgarray[N:static int, T] = CpuGpu[array[N,T], ptr array[N,T]]
template `[]`*(x: cgarray, i: SomeNumber): auto = x.cpu[i]
template `[]=`*(x: cgarray, i: SomeNumber, y: auto): auto =
  x.cpu[i] = y


template toGpu*(x: cgarray): auto =
  if x.gpu.isNil:
    x.gpu.gpuMalloc
  gpuMemCpyToGPU(x.gpu, addr x.cpu[0], sizeof(x.cpu))
  x.gpu
template getGpu*(x: cgarray, g: ptr array): auto = g[]
template fromGpu*(x: cgarray, g: ptr array) =
  gpuMemCpyToCPU(addr x.cpu[0], x.gpu, sizeof(x.cpu))


block:
  const N = 1024
  var x: cgarray[N,float]
  var y: cgarray[N,float]
  var z: cgarray[N,float]
  for i in 0..<N:
    x[i] = 1.0 + i.float
    y[i] = 2.0 + i.float
  onGpu:
    let s = getNumThreads()
    var i = getThreadNum()
    while i < N:
      z[i] = x[i] * y[i]
      i += s
  echo "z[N-1]: ", z[N-1]
