import backend/[accel,cpugpu]
import bench/commonBench
import base
import strformat

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
  let N = intParam("n", -1)
  let Nmin = intParam("nmin", if N>0: N else: 1024)
  let Nmax = intParam("nmax", if N>0: N else: 128*1024*1024)
  type Real = float32
  var x = newCgarray[Real](Nmax)
  var y = newCgarray[Real](Nmax)
  var z = newCgarray[Real](Nmax)
  template initX(i: int): Real = Real(1 + i)
  template initY(i: int): Real = Real(2 + i)
  threads:
    for i in threadRange(Nmax):
      x[i] = initX(i)
      y[i] = initY(i)
  x.gpuReadOnly
  y.gpuReadOnly
  z.gpuWriteOnly

  proc test(n: int) =
    var b = newBench()
    benchSingle(b):
      let nrep = b.nrep
      #echo nrep
      onGpu:
        for rep in 1..nrep:
          for i in gpuRange(n):
            z[i] = x[i] * y[i]
    #echo "z[n-1]: ", z[n-1]
    var errcnt = 0
    for i in 0..<n:
      let r = initX(i) * initY(i)
      if z[i] != r:
        if errcnt < 1:
          echo "Error: ", i, " ", z[i], " ", r, " ", z[i]-r
        inc errcnt
    if errcnt > 0:
      echo "Error count: ", errcnt
    let bytes = 3*n*sizeof(x[0])
    let flops = n
    let memMB = 1e-6 * bytes
    let bwGB = bytes * b.perNs
    let flopsGB = flops * b.perNs
    echo &"{memMB:8.3f} {bwGB:8.3f} {flopsGB:8.3f}"
    doAssert(errcnt==0)

  echo "gpuDefaultNumThreads: ", gpuDefaultNumThreads()
  echo "      MB     GB/s  GFlop/s    (z=x*y ", (if sizeof(Real)==8:"double)" else:"single)")
  var n = Nmin
  while n <= Nmax:
    test(n)
    n *= 2

  commsFinalize()
