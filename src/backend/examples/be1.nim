import backend/[accel,cpugpu]
import bench/commonBench
import base
import strformat

type
  Myarray[T] = object
    n: int
    p: ptr UncheckedArray[T]
template bytes*[T](x: Myarray[T]): auto = x.n * sizeof(T)
template `[]`*(x: Myarray, i: SomeInteger): auto = x.p[][i]
template `[]=`*(x: Myarray, i: SomeInteger, y: auto): auto =
  x.p[][i] = y

type
  Cgarray[T] = CpuGpu[Myarray[T], Myarray[T]]
proc newCgarray*[T](n: int): Cgarray[T] =
  result.cpu.n = n
  result.cpu.p = cast[ptr UncheckedArray[T]](T.createU(n))
  result.gpu.n = n
template `[]`*(x: Cgarray, i: SomeInteger): auto = x.cpu[i]
template `[]=`*(x: Cgarray, i: SomeInteger, y: auto): auto =
  x.cpu[i] = y

template toGpu*(g: var Myarray, x: Cgarray, cpy: bool) =
  threadSingle:
    if g.p.isNil:
      when backendIsGpu:
        g.p = cast[type g.p](gpuMalloc(g.bytes))
      else:
        g.p = x.cpu.p
    when backendIsGpu:
      if cpy:
        #echo "Copy in"
        gpuMemCpyToGPU(g.p, x.cpu.p, g.bytes)

template getGpu*(x: Cgarray, g: Myarray): auto = g

template fromGpu*(x: var Cgarray, g: Myarray, cpy: bool) =
  when backendIsGpu:
    threadSingle:
      if cpy:
        #echo "Copy out"
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
  template initX(i: int): Real = Real(1 + i mod 99)
  template initY(i: int): Real = Real(2 + i mod 99)
  threads:
    for i in threadRange(Nmax):
      x[i] = initX(i)
      y[i] = initY(i)

  proc check(n: int) =
    var errcnt = 0
    for i in 0..<n:
      #let r = initX(i) * initY(i)
      let r = x[i] * y[i]
      if z[i] != r:
        if errcnt < 10:
          echo "Error: ", i, " ", z[i], " ", r, " ", z[i]-r
        inc errcnt
    if errcnt > 0:
      echo "Error count: ", errcnt
    doAssert(errcnt==0)

  proc perf(n: int, b: Bench) =
    let bytes = 3*n*sizeof(x[0])
    let flops = n
    let memMB = 1e-6 * bytes
    let bwGB = bytes * b.perNs
    let flopsGB = flops * b.perNs
    echo &"{memMB:8.3f} {bwGB:9.3f} {flopsGB:8.3f}"

  proc testGpu(n: int) =
    #echo "testGpu: ", n
    x.cpuWriteOnly; y.cpuWriteOnly; z.cpuReadOnly
    x.gpuReadOnly; y.gpuReadOnly; z.gpuWriteOnly
    template set(a,b: SomeNumber, c = 0) =
      threads:
        for i in threadRange(n):
          x[i] += a
          y[i] += b
          z[i] += c
    template run(body: untyped) =
      #echo "run"
      onGpu:
        for i in gpuRange(n):
          z[i] = x[i] * y[i]
      body
      check(n)
    template runnw(body1,body2: untyped) =
      let fin = onGpuNowait:
        for i in gpuRange(n):
          z[i] = x[i] * y[i]
      body1  # runs on CPU while kernel running on GPU
      fin()  # wait for GPU and sync data
      body2
      check(n)
    run: discard
    doAssert(x.wasCopiedIn and not x.wasCopiedOut)
    doAssert(y.wasCopiedIn and not y.wasCopiedOut)
    doAssert(z.wasCopiedOut and not z.wasCopiedIn)
    set(1,2)
    run: discard
    doAssert(x.wasCopiedIn and not x.wasCopiedOut)
    doAssert(y.wasCopiedIn and not y.wasCopiedOut)
    doAssert(z.wasCopiedOut and not z.wasCopiedIn)
    z.cpuWriteOnly  # disable copy out to z
    run: discard
    doAssert(x.wasCopiedIn and not x.wasCopiedOut)
    doAssert(y.wasCopiedIn and not y.wasCopiedOut)
    doAssert(not z.wasCopiedOut and not z.wasCopiedIn)
    z.cpuReadOnly  # enable copy out to z
    set(-1,-1)
    run: discard
    doAssert(x.wasCopiedIn and not x.wasCopiedOut)
    doAssert(y.wasCopiedIn and not y.wasCopiedOut)
    doAssert(z.wasCopiedOut and not z.wasCopiedIn)
    x.cpuReadOnly; y.cpuReadOnly  # disable copy in for x and y
    run: discard
    doAssert(not x.wasCopiedIn and not x.wasCopiedOut)
    doAssert(not y.wasCopiedIn and not y.wasCopiedOut)
    doAssert(z.wasCopiedOut and not z.wasCopiedIn)
    x.cpuWriteOnly; y.cpuWriteOnly  # enable copy in for x and y
    set(1,1)
    run: discard
    doAssert(x.wasCopiedIn and not x.wasCopiedOut)
    doAssert(y.wasCopiedIn and not y.wasCopiedOut)
    doAssert(z.wasCopiedOut and not z.wasCopiedIn)
    set(-1,-1)
    runnw: discard
    do: discard
    doAssert(x.wasCopiedIn and not x.wasCopiedOut)
    doAssert(y.wasCopiedIn and not y.wasCopiedOut)
    doAssert(z.wasCopiedOut and not z.wasCopiedIn)

  proc testNested(n: int) =
    #echo "testGpu: ", n
    x.cpuWriteOnly; y.cpuWriteOnly; z.cpuReadOnly
    x.gpuReadOnly; y.gpuReadOnly; z.gpuWriteOnly
    template set(a,b: SomeNumber, c = 0) =
      threads:
        for i in threadRange(n):
          x[i] += a
          y[i] += b
          z[i] += c
    template setr(a,b: SomeNumber, c = 0) =
      for i in threadRange(n):
        x[i] += a
        y[i] += b
        z[i] += c
    template run(body: untyped) =
      threads:
        onGpu:
          for i in gpuRange(n):
            z[i] = x[i] * y[i]
        body
      check(n)
    template runnw(body1,body2: untyped) =
      threads:
        let fin = onGpuNowait:
          for i in gpuRange(n):
            z[i] = x[i] * y[i]
        body1  # runs on CPU while kernel running on GPU
        fin()  # wait for GPU and sync data
        body2
      check(n)
    run: discard
    doAssert(x.wasCopiedIn and not x.wasCopiedOut)
    doAssert(y.wasCopiedIn and not y.wasCopiedOut)
    doAssert(z.wasCopiedOut and not z.wasCopiedIn)
    set(1,2)
    run: discard
    doAssert(x.wasCopiedIn and not x.wasCopiedOut)
    doAssert(y.wasCopiedIn and not y.wasCopiedOut)
    doAssert(z.wasCopiedOut and not z.wasCopiedIn)
    z.cpuWriteOnly  # disable copy out to z
    run: discard
    doAssert(x.wasCopiedIn and not x.wasCopiedOut)
    doAssert(y.wasCopiedIn and not y.wasCopiedOut)
    doAssert(not z.wasCopiedOut and not z.wasCopiedIn)
    z.cpuReadOnly  # enable copy out to z
    set(-1,-1)
    run: discard
    doAssert(x.wasCopiedIn and not x.wasCopiedOut)
    doAssert(y.wasCopiedIn and not y.wasCopiedOut)
    doAssert(z.wasCopiedOut and not z.wasCopiedIn)
    x.cpuReadOnly; y.cpuReadOnly  # disable copy in for x and y
    run: discard
    doAssert(not x.wasCopiedIn and not x.wasCopiedOut)
    doAssert(not y.wasCopiedIn and not y.wasCopiedOut)
    doAssert(z.wasCopiedOut and not z.wasCopiedIn)
    x.cpuWriteOnly; y.cpuWriteOnly  # enable copy in for x and y
    set(1,1)
    run: discard
    doAssert(x.wasCopiedIn and not x.wasCopiedOut)
    doAssert(y.wasCopiedIn and not y.wasCopiedOut)
    doAssert(z.wasCopiedOut and not z.wasCopiedIn)
    set(-1,-1)
    runnw: discard
    do: discard
    doAssert(x.wasCopiedIn and not x.wasCopiedOut)
    doAssert(y.wasCopiedIn and not y.wasCopiedOut)
    doAssert(z.wasCopiedOut and not z.wasCopiedIn)

  proc benchCpu(n: int) =
    var b = newBench()
    benchSingle(b):
      let nrep = b.nrep
      threads:
        for rep in 1..nrep:
          for i in threadRange(n):
            z[i] = x[i] * y[i]
    perf(n, b)
    check(n)

  proc benchGpu(n: int) =
    var b = newBench()
    benchSingle(b):
      let nrep = b.nrep
      onGpu(n):
        for rep in 1..nrep:
          for i in gpuRange(n):
            z[i] = x[i] * y[i]
    perf(n, b)
    check(n)

  proc benchNested(n: int) =
    var b = newBench()
    benchSingle(b):
      let nrep = b.nrep
      threads:
        onGpu(n):
          for rep in 1..nrep:
            for i in gpuRange(n):
              z[i] = x[i] * y[i]
    perf(n, b)
    check(n)

  var ns = newSeq[int](0)
  block:
    var n = Nmin
    while n <= Nmax:
      ns.add n
      n *= 2

  echo "Testing GPU"
  for n in ns: testGpu(n)
  #testGpu(Nmin)
  echo "Testing Nested"
  for n in ns: testNested(n)

  echo "CPU   MB      GB/s  GFlop/s    (z=x*y ", (if sizeof(Real)==8:"double)" else:"single)")
  for n in ns: benchCpu(n)
  echo "GPU   MB      GB/s  GFlop/s    (z=x*y ", (if sizeof(Real)==8:"double)" else:"single)")
  for n in ns: benchGpu(n)
  echo "Nest  MB      GB/s  GFlop/s    (z=x*y ", (if sizeof(Real)==8:"double)" else:"single)")
  for n in ns: benchNested(n)

  commsFinalize()




#[
import backend/accel

block:
  var x = 1.0'f32
  let y = 2.0'f32
  const z = 3.0'f32
  echo "x: ", x
  onGpu:
    x = y * z
  echo "x: ", x
  doAssert x == 6.0
]#
