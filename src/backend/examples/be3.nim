import backend/[accel,cpugpu]
import bench/commonBench
import base
import strformat

const V = 16
proc eq[T](r,x: ptr UncheckedArray[T], s: int) =
  for i in 0..<V:
    let k = s+i
    r[k] = x[i]
proc eqmul[T](r,x,y: ptr UncheckedArray[T], s: int) =
  for i in 0..<V:
    let k = s+i
    r[i] = x[k] * y[k]
proc peqmul[T](r,x,y: ptr UncheckedArray[T], s: int) =
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

const Nc = 3
type
  Real = float32
  cmatarr = Mat[Nc,Nc,Real]
  cgmatarr = CpuGpu[cmatarr, cmatarr]
proc newCgmatarr*(n: int): cgmatarr =
  result.cpu = newMat[Nc,Nc,Real](n)
  when backendIsGpu:
    let gp = cast[ptr UncheckedArray[Real]](gpuMalloc(n*Nc*Nc*sizeof(Real)))
  else:
    let gp = result.cpu.dat
  result.gpu = newMat[Nc,Nc,Real](n, gp)
template `[]`*(x: cgmatarr, i,j: SomeNumber): auto = x.cpu[i,j]

template toGpu*(g: var Mat, x: cgmatarr, cpy: bool) =
  when backendIsGpu:
    if cpy:
      gpuMemCpyToGPU(g.dat, x.cpu.dat, g.bytes)

template getGpu*(x: cgmatarr, g: Mat): auto = g

template fromGpu*(x: var cgmatarr, g: Mat, cpy: bool) =
  when backendIsGpu:
    if cpy:
      gpuMemCpyToCPU(x.cpu.dat, g.dat, x.cpu.bytes)

block:
  commsInit()
  let N = intParam("n", 1024)
  var x = newCgmatarr(N*V)
  var y = newCgmatarr(N*V)
  var z = newCgmatarr(N*V)
  #template initX(s,ic,jc: int): Real = 1.0.Real + s.Real + ic.Real + jc.Real
  #template initY(s,ic,jc: int): Real = 2.0.Real + s.Real + ic.Real + jc.Real
  template initX(s,ic,jc: int): Real = 3.0.Real
  template initY(s,ic,jc: int): Real = 3.0.Real
  for ic in 0..<Nc:
    for jc in 0..<Nc:
      for s in 0..<N*V:
        x[ic,jc][s] := initX(s,ic,jc)
        y[ic,jc][s] := initY(s,ic,jc)
  #for ic in 0..<Nc:
  #  for jc in 0..<Nc:
  #    for s in 0..<N*V:
  #      var r = x[ic,0][s] * y[0,jc][s]
  #      for kc in 1..<Nc:
  #        r += x[ic,kc][s] * y[kc,jc][s]
  #      z[ic,jc][s] = r
  x.gpuReadOnly
  y.gpuReadOnly
  z.gpuWriteOnly
  var b = newBench()
  benchSingle(b):
    let nrep = b.nrep
    #echo nrep
    onGpu:
      for rep in 1..nrep:
        when backendIsGpu:
          for s in gpuRange(N*V):
            for ic in 0..<Nc:
              for jc in 0..<Nc:
                var r = x[ic,0][s] * y[0,jc][s]
                for kc in 1..<Nc:
                  r += x[ic,kc][s] * y[kc,jc][s]
                z[ic,jc][s] = r
        else:
          for so in gpuRange(N):
            for ic in 0..<Nc:
              for jc in 0..<Nc:
                var r {.noInit.}: array[V,Real]
                var p = cast[ptr UncheckedArray[Real]](addr r[0])
                eqmul(p, x[ic,0], y[0,jc], V*so)
                for kc in 1..<Nc:
                  peqmul(p, x[ic,kc], y[kc,jc], V*so)
                eq(z[ic,jc], p, V*so)
  #echo "z[N-1]: ", z[N-1]
  var errcnt = 0
  for s in 0..<N*V:
    for ic in 0..<Nc:
      for jc in 0..<Nc:
        var r = initX(s,ic,0) * initY(s,0,jc)
        for kc in 1..<Nc:
          r += initX(s,ic,kc) * initY(s,kc,jc)
        let t = z[ic,jc][s]
        let d = t - r
        if abs(d) > 10.0 * Real.epsilon * abs(r):
          if errcnt < 10:
            echo &"Error: {s},{ic},{jc} {t} {r} {t-r}"
          inc errcnt
  if errcnt > 0:
    echo "Error count: ", errcnt
  let flops = N*V*Nc*Nc*(2*Nc-1)
  let bytes = 3*N*V*Nc*Nc*sizeof(Real)
  echo "Gflops: ", flops * b.perNs
  echo "Gbytes: ", bytes * b.perNs
  doAssert(errcnt==0)
  commsFinalize()
