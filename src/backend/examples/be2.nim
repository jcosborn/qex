import backend/[accel,cpugpu]
import bench/commonBench
import base
import strformat

type
  SiteV[V:static int] = distinct int
template `[]`*(x: SiteV): int = int(x)
#converter toInt(x: SiteV): int {.inline.} = int(x)
template `[]`[T;V:static int](x: ptr UncheckedArray[T], i: SiteV[V]): ptr array[V,T] =
  cast[ptr array[V,T]](addr x[][V*i[]])
template `[]=`[T;V:static int](r: ptr UncheckedArray[T], i: SiteV[V], x: array[V,T]) =
  let i0 = V*i[]
  for k in 0..<V:
    r[i0+k] = x[k]

when backendIsGpu:
  iterator gpuSites(n:int, V:static int): SiteV[1] =
    for s in gpuRange(n*V):
      yield SiteV[1](s)
else:
  iterator gpuSites(n:int, V:static int): SiteV[V] =
    for s in gpuRange(n):
      yield SiteV[V](s)


proc `+=`*[N:static int;T](r: var array[N,T], x: array[N,T]) =
  for i in 0..<N:
    r[i] += x[i]
proc `*`*[N:static int;T](x: ptr array[N,T], y: ptr array[N,T]): array[N,T] =
  for i in 0..<N:
    result[i] = x[i] * y[i]

const V = 16
proc eq[T](r,x: ptr UncheckedArray[T], s: int) =
  for i in 0..<V:
    let k = s+i
    r[k] = x[i]
proc eq[N:static int; T](r: ptr UncheckedArray[T], x: array[N,T], s: int) =
  for i in 0..<V:
    let k = s+i
    r[k] = x[i]
proc peq[T](r,x: ptr UncheckedArray[T], s: int) =
  for i in 0..<V:
    let k = s+i
    r[k] += x[i]
proc peq[N:static int; T](r: ptr UncheckedArray[T], x: array[N,T], s: int) =
  for i in 0..<V:
    let k = s+i
    r[k] += x[i]
proc eqmul[T](r,x,y: ptr UncheckedArray[T], s: int) =
  for i in 0..<V:
    let k = s+i
    r[i] = x[k] * y[k]
proc eqmul[N:static int; T](r: var array[N,T], x,y: ptr UncheckedArray[T], s: int) =
  for i in 0..<V:
    let k = s+i
    r[i] = x[k] * y[k]
proc peqmul[T](r,x,y: ptr UncheckedArray[T], s: int) =
  for i in 0..<V:
    let k = s+i
    r[i] += x[k] * y[k]
proc peqmul[N:static int; T](r: var array[N,T], x,y: ptr UncheckedArray[T], s: int) =
  for i in 0..<V:
    let k = s+i
    r[i] += x[k] * y[k]

type
  Mat[N,M:static int, T] = object
    count: int
    dat: ptr UncheckedArray[T]
    mat: array[N*M, ptr UncheckedArray[T]]
proc newMat[N,M:static int, T](c: int, p: ptr UncheckedArray[T]): Mat[N,M,T] =
  result.count = c
  result.dat = p
  var t = p
  for i in 0..<N:
    for j in 0..<M:
      result.mat[i*M+j] = t
      t = cast[ptr UncheckedArray[T]](addr t[][c])
proc newMat[N,M:static int, T](c: int): Mat[N,M,T] =
  let p = cast[ptr UncheckedArray[T]](T.createU(N*M*c))
  newMat[N,M,T](c, p)
template bytes[N,M:static int, T](x: Mat[N,M,T]): int = (N*M*sizeof(T))*x.count
template `[]`*[N,M:static int, T](x: Mat[N,M,T], i,j: SomeNumber): auto = x.mat[i*M+j]

type
  cmatarr[N:static int; T] = Mat[N,N,T]
  cgmatarrObj[N:static int; T] = CpuGpu[cmatarr[N,T], cmatarr[N,T]]
  cgmatarr[N:static int; T] = ref CpuGpu[cmatarr[N,T], cmatarr[N,T]]
#proc destroy*(x: var cgmatarr) =
#  when backendIsGpu:
#    if x.gpu.dat != nil:
#      gpuFree(x.gpu.dat)
#  if x.cpu.dat != nil:
#    dealloc(x.cpu.dat)
proc finalize[T:cgmatarr](x: T) =
  #echo "finalize cgmatarr"
  when backendIsGpu:
    if x.gpu.dat != nil:
      gpuFree(x.gpu.dat)
  if x.cpu.dat != nil:
    dealloc(x.cpu.dat)
proc newCgmatarr*[N:static int; T](n: int): cgmatarr[N,T] =
  result.new(finalize[typeof result])
  result.cpu = newMat[N,N,T](n)
  when backendIsGpu:
    let gp = cast[ptr UncheckedArray[T]](gpuMalloc(n*N*N*sizeof(T)))
  else:
    let gp = result.cpu.dat
  result.gpu = newMat[N,N,T](n, gp)

template `[]`*(x: cgmatarr, i,j: SomeNumber): auto = x.cpu[i,j]

template toGpu*(g: var Mat, x: cgmatarr, cpy: bool) =
  when backendIsGpu:
    if cpy:
      gpuMemCpyToGPU(g.dat, x.cpu.dat, g.bytes)

template getGpu*(x: cgmatarr, g: Mat): auto = g

template fromGpu*(x: cgmatarr, g: Mat, cpy: bool) =
  when backendIsGpu:
    if cpy:
      gpuMemCpyToCPU(x.cpu.dat, g.dat, x.cpu.bytes)

block:
  threadsInit()
  commsInit()
  let N0 = intParam("n", -1)
  let Nmin = intParam("nmin", if N0>0: N0 else: 1024)
  let Nmax = intParam("nmax", if N0>0: N0 else: 2*1024*1024)
  var ns = newSeq[int](0)
  block:
    var n = Nmin
    while n < Nmax:
      ns.add n
      ns.add int(round(n*sqrt(2.0)))
      n *= 2
    ns.add n
  echo ns

  proc runtest(Nc: static int; Real: typedesc; n: int) =
    let np = n + 3  # padded
    var x = newCgmatarr[Nc,Real](np*V)
    var y = newCgmatarr[Nc,Real](np*V)
    var z = newCgmatarr[Nc,Real](np*V)
    template initX(s,ic,jc: int): Real = 1.0.Real + s.Real + ic.Real + jc.Real
    template initY(s,ic,jc: int): Real = 2.0.Real + s.Real + ic.Real + jc.Real
    threads:
      for ic in 0..<Nc:
        for jc in 0..<Nc:
          for s in threadRange(n*V):
            x[ic,jc][s] := initX(s,ic,jc)
            y[ic,jc][s] := initY(s,ic,jc)
    x.gpuReadOnly
    y.gpuReadOnly
    z.gpuWriteOnly
    var b = newBench()
    benchSingle(b):
      let nrep = b.nrep
      onGpu(n*V):
        for rep in 0..<nrep:
          when false:
            for s in gpuSites(n,V):
              for ic in 0..<Nc:
                for jc in 0..<Nc:
                  var r = x[ic,0][s] * y[0,jc][s]
                  for kc in 1..<Nc:
                    r += x[ic,kc][s] * y[kc,jc][s]
                  z[ic,jc][s] = r
          else:
            when backendIsGpu:
              for s in gpuRange(n*V):
                for ic in 0..<Nc:
                  for jc in 0..<Nc:
                    var r = x[ic,0][s] * y[0,jc][s]
                    for kc in 1..<Nc:
                      r += x[ic,kc][s] * y[kc,jc][s]
                    z[ic,jc][s] = r
            else:
              for so in gpuRange(n):
                for ic in 0..<Nc:
                  for jc in 0..<Nc:
                    var r {.noInit.}: array[V,Real]
                    var p = cast[ptr UncheckedArray[Real]](addr r[0])
                    eqmul(p, x[ic,0], y[0,jc], V*so)
                    for kc in 1..<Nc:
                      peqmul(p, x[ic,kc], y[kc,jc], V*so)
                    eq(z[ic,jc], p, V*so)
    var errcnt = 0
    for s in 0..<n*V:
      for ic in 0..<Nc:
        for jc in 0..<Nc:
          var r = x[ic,0][s] * y[0,jc][s]
          for kc in 1..<Nc:
            r += x[ic,kc][s] * y[kc,jc][s]
          let t = z[ic,jc][s]
          let d = t - r
          if abs(d) > 10.0 * Real.epsilon * abs(r):
            if errcnt < 10:
              echo &"Error: {s},{ic},{jc} {t} {r} {t-r}"
            inc errcnt
    if errcnt > 0:
      echo "Error count: ", errcnt
    let bytes = 3*n*V*Nc*Nc*sizeof(Real)
    let flops = n*V*Nc*Nc*(2*Nc-1)
    let memMB = 1e-6 * bytes
    let gb = bytes * b.perNs
    let gf = flops * b.perNs
    echo &"{memMB:8.3f} {gb:8.3f} {gf:8.3f}"
    doAssert(errcnt==0)

  #runtest(float64)
  for n in ns:
    runtest(3, float32, n)


  commsFinalize()
