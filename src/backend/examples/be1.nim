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
proc destroy*(x: var Cgarray) =
  when backendIsGpu:
    if x.gpu.p != nil:
      gpuFree(x.gpu.p)
  if x.cpu.p != nil:
    dealloc(x.cpu.p)
template `[]`*(x: Cgarray, i: SomeInteger): auto = x.cpu[i]
template `[]=`*(x: Cgarray, i: SomeInteger, y: auto): auto =
  x.cpu[i] = y

template toGpu*(g: var Myarray, x: Cgarray, cpy: bool) =
  threadSingle:
    if g.p.isNil:
      when backendIsGpu:
        #echo "gpuMalloc"
        g.p = cast[type g.p](gpuMalloc(g.bytes))
        #echo "... done"
      else:
        g.p = x.cpu.p
    when backendIsGpu:
      if cpy:
        #echo "gpuMemCpyToGPU"
        gpuMemCpyToGPU(g.p, x.cpu.p, g.bytes)
        #echo "... done"

template getGpu*(x: Cgarray, g: Myarray): auto = g

template fromGpu*(x: var Cgarray, g: Myarray, cpy: bool) =
  when backendIsGpu:
    threadSingle:
      if cpy:
        #echo "Copy out"
        gpuMemCpyToCPU(x.cpu.p, g.p, x.cpu.bytes)

block:
  threadsInit()
  commsInit()
  let N = intParam("n", -1)
  let Nmin = intParam("nmin", if N>0: N else: 1024)
  let Nmax = intParam("nmax", if N>0: N else: 128*1024*1024)
  var ns = newSeq[int](0)
  block:
    var n = Nmin
    while n < Nmax:
      ns.add n
      ns.add int(round(n*sqrt(2.0)))
      n *= 2
    ns.add n
  echo ns
  type
    BenchResult = object
      name: string
      n: int
      threadLoc: string
      b: Bench
  var br: BenchResult
  var brs = newSeq[BenchResult](0)

  proc runtest(Real: typedesc) =
    var x = newCgarray[Real](Nmax)
    var y = newCgarray[Real](Nmax)
    var z = newCgarray[Real](Nmax)
    template initX(i: int): Real = Real(1 + i mod 99)
    template initY(i: int): Real = Real(2 + i mod 99)
    threads:
      for i in threadRange(Nmax):
        x[i] = initX(i)
        y[i] = initY(i)

    proc check2(n: int, ii: auto) =
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
        echo "n: ", n, "  ", ii
      doAssert(errcnt==0)
    template check(n: int) =
      check2(n, instantiationInfo())

    proc mem(br: BenchResult): float=
      let n = br.n
      let bytes = 3*n*sizeof(x[0])
      let memMB = 1e-6 * bytes
      memMB

    proc bw(br: BenchResult): float =
      let n = br.n
      let bytes = 3*n*sizeof(x[0])
      let bwGB = bytes * br.b.perNs
      bwGB

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
          for i in threadRangeV(n):
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
      br.name = "CPU"
      br.n = n
      block:
        br.threadLoc = "out"
        threads: (for i in threadRangeAligned(n,16): z[i] = 0)
        var b = newBench()
        benchSingle(b):
          let nrep = b.nrep
          threads:
            for rep in 0..<nrep:
              for i in threadRangeAligned(n,16):
                z[i] = x[i] * y[i]
        #perf(n, b)
        br.b = b
        brs.add br
        check(n)
      block:
        br.threadLoc = "in"
        threads: (for i in threadRangeAligned(n,16): z[i] = 0)
        var b = newBench()
        benchSingle(b):
          let nrep = b.nrep
          for rep in 0..<nrep:
            threads:
              for i in threadRangeAligned(n,16):
                z[i] = x[i] * y[i]
        #perf(n, b)
        br.b = b
        brs.add br
        check(n)

    proc benchGpu(n: int) =
      br.name = "GPU"
      br.n = n
      block:
        br.threadLoc = "out"
        threads: (for i in threadRangeAligned(n,16): z[i] = 0)
        z.copyToGpu
        var b = newBench()
        benchSingle(b):
          let nrep = b.nrep
          onGpu(n):
            for rep in 0..<nrep:
              for i in gpuRange(n):
                z[i] = x[i] * y[i]
        #perf(n, b)
        br.b = b
        brs.add br
        check(n)
      when false: #block:
        br.threadLoc = "in"
        threads: (for i in threadRangeAligned(n,16): z[i] = 0)
        x.copyToGpu; x.cpuUnused
        y.copyToGpu; y.cpuUnused
        z.copyToGpu; z.cpuUnused
        var cb = newSeq[proc()](0)
        var cbr = 0
        var b = newBench()
        let maxrun = 1
        benchSingle(b):
          let nrep = b.nrep
          cb.setLen(nrep)
          for rep in 0 ..< nrep:
            cb[rep] = onGpuNowait(n):
              for i in gpuRange(n):
                z[i] = x[i] * y[i]
            if (rep+1) mod maxrun == 0:
              for i in cbr .. rep : cb[i]()
              cbr = rep + 1
          for i in cbr ..< nrep : cb[i]()
          #for rep in 0 ..< nrep: cb[rep]()
        #perf(n, b)
        br.b = b
        brs.add br
        x.cpuWriteOnly; y.cpuWriteOnly
        z.copyFromGpu; z.cpuReadOnly
        check(n)
      block:
        br.threadLoc = "ins"
        threads: (for i in threadRangeAligned(n,16): z[i] = 0)
        x.copyToGpu; x.cpuUnused
        y.copyToGpu; y.cpuUnused
        z.copyToGpu; z.cpuUnused
        var b = newBench()
        benchSingle(b):
          let nrep = b.nrep
          for rep in 0 ..< nrep:
            onGpu(n):
              for i in gpuRange(n):
                z[i] = x[i] * y[i]
        #perf(n, b)
        br.b = b
        brs.add br
        x.cpuWriteOnly; y.cpuWriteOnly
        z.copyFromGpu; z.cpuReadOnly
        check(n)

    proc benchNested(n: int) =
      br.name = "Nest"
      br.n = n
      block:
      #when false:
        br.threadLoc = "out"
        threads: (for i in threadRangeAligned(n,16): z[i] = 0)
        z.copyToGpu
        var b = newBench()
        benchSingle(b):
          let nrep = b.nrep
          threads:
            onGpu(n):
              for rep in 0 ..< nrep:
                for i in gpuRange(n):
                  z[i] = x[i] * y[i]
        #perf(n, b)
        br.b = b
        brs.add br
        check(n)
      block:
      #when false:
        br.threadLoc = "in"
        threads: (for i in threadRangeAligned(n,16): z[i] = 0)
        x.copyToGpu; x.cpuUnused
        y.copyToGpu; y.cpuUnused
        z.copyToGpu; z.cpuUnused
        var b = newBench()
        benchSingle(b):
          let nrep = b.nrep
          threads:
            var cb = newSeq[proc()](nrep)
            for rep in 0 ..< nrep:
              cb[rep] = onGpuNowait(n):
                for i in gpuRange(n):
                  z[i] = x[i] * y[i]
            for rep in 0 ..< nrep:
              cb[rep]()
        #perf(n, b)
        br.b = b
        brs.add br
        x.cpuWriteOnly; y.cpuWriteOnly
        z.copyFromGpu; z.cpuReadOnly
        check(n)
      block:
      #when false:
        br.threadLoc = "ins"
        threads: (for i in threadRangeAligned(n,16): z[i] = 0)
        x.copyToGpu; x.cpuUnused
        y.copyToGpu; y.cpuUnused
        z.copyToGpu; z.cpuUnused
        var b = newBench()
        benchSingle(b):
          let nrep = b.nrep
          threads:
            for rep in 0 ..< nrep:
              onGpu(n):
                for i in gpuRange(n):
                  z[i] = x[i] * y[i]
        #perf(n, b)
        br.b = b
        brs.add br
        x.cpuWriteOnly; y.cpuWriteOnly
        z.copyFromGpu; z.cpuReadOnly
        check(n)

    echo "Testing GPU ", $Real
    for n in ns: testGpu(n)
    #testGpu(Nmin)
    echo "Testing Nested ", $Real
    for n in ns: testNested(n)

    proc echobr =
      var lastn = 0
      var s = ""
      for br in brs:
        if br.n != lastn:
          if lastn > 0: echo s
          s = &"{br.mem:8.3f}"
          lastn = br.n
        s &= &" {br.bw:9.3f}"
      echo s

    brs.setLen(0)
    for n in ns: benchCpu(n)
    #echo "CPU   MB      GB/s  GFlop/s    (z=x*y ", $Real, ")"
    echo "CPU   MB      GB/s      (z=x*y ", $Real, ")"
    echobr()

    brs.setLen(0)
    for n in ns: benchGpu(n)
    #echo "GPU   MB      GB/s  GFlop/s    (z=x*y ", $Real, ")"
    echo "GPU   MB      GB/s      (z=x*y ", $Real, ")"
    echobr()

    brs.setLen(0)
    for n in ns: benchNested(n)
    #echo "Nest  MB      GB/s  GFlop/s    (z=x*y ", $Real, ")"
    echo "Nest  MB      GB/s      (z=x*y ", $Real, ")"
    echobr()


  runtest(float64)
  runtest(float32)

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
