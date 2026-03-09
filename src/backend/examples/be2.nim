import backend/[accel,cpugpu]
import bench/commonBench
import base

type
  myarray[T] = object
    n: int
    p: ptr UncheckedArray[T]
template bytes*[T](x: myarray[T]): auto = x.n * sizeof(T)
template `[]`*(x: myarray, i: SomeNumber): auto = x.p[][i]
template `[]=`*(x: myarray, i: SomeNumber, y: auto): auto =
  x.p[][i] = y

type
  cgarray[T] = CpuGpu[myarray[T], myarray[T]]
proc newCgarray*[T](n: int): cgarray[T] =
  result.cpu.n = n
  result.cpu.p = cast[ptr UncheckedArray[T]](T.createU(n))
  result.gpu.n = n
template `[]`*(x: cgarray, i: SomeNumber): auto = x.cpu[i]
template `[]=`*(x: cgarray, i: SomeNumber, y: auto): auto =
  x.cpu[i] = y

template toGpu*(g: var myarray, x: cgarray, cpy: bool) =
  if g.p.isNil:
    when backendIsGpu:
      g.p = cast[type g.p](gpuMalloc(g.bytes))
    else:
      g.p = x.cpu.p
  when backendIsGpu:
    if cpy:
      gpuMemCpyToGPU(g.p, x.cpu.p, g.bytes)

template getGpu*(x: cgarray, g: myarray): auto = g

template fromGpu*(x: var cgarray, g: myarray, cpy: bool) =
  when backendIsGpu:
    if cpy:
      gpuMemCpyToCPU(x.cpu.p, g.p, x.cpu.bytes)

block:
  commsInit()
  let N = intParam("n", 1024*1024)
  var x = newCgarray[float](N)
  var y = newCgarray[float](N)
  var z = newCgarray[float](N)
  template initX(i: int): float = 1.0 + i.float
  template initY(i: int): float = 2.0 + i.float
  for i in 0..<N:
    x[i] = initX(i)
    y[i] = initY(i)
  x.gpuReadOnly
  y.gpuReadOnly
  z.gpuWriteOnly
  var b = newBench()
  benchSingle(b):
    let nrep = b.nrep
    #echo nrep
    onGpu:
      for rep in 1..nrep:
        for i in gpuRange(N):
          z[i] = x[i] * y[i]
  echo "z[N-1]: ", z[N-1]
  var errcnt = 0
  for i in 0..<N:
    let r = initX(i) * initY(i)
    if z[i] != r:
      if errcnt < 1:
        echo "Error: ", i, " ", z[i], " ", r, " ", z[i]-r
      inc errcnt
  if errcnt > 0:
    echo "Error count: ", errcnt
  let flops = N
  let bytes = 3*N*sizeof(x[0])
  echo "Gflops: ", flops * b.perNs
  echo "Gbytes: ", bytes * b.perNs
  doAssert(errcnt==0)
  commsFinalize()
