import qex
import backend/[accel,cpugpu]

type
  #FieldProxy[T] = object
  #  f: T
  GpuFieldObj*[V:static int, T] = object
    n*: int
    p*: ptr UncheckedArray[T]
  GpuFieldExpr*[OP:static string, A] = object
    args*: A
  GpuFieldExpr2*[OP:static string, A] = GpuFieldExpr[OP,A]
  GpuField*[V:static int, T] = GpuFieldObj[V,T]
  GpuField2*[V:static int, T] = GpuFieldObj[V,T]
  SomeGpuField* = GpuField | GpuFieldExpr
  SomeGpuField2* = GpuField2 | GpuFieldExpr2
template gpuSites*(x: GpuField): auto = gpuSites(x.n, x.V)

proc bytes*[T:GpuField](x: T): int = x.n * sizeof(T.T)

template index*[T;V,L:static int](X: typedesc[GpuField[V,T]], I: typedesc[SiteV[L]]): typedesc =
  when V == L:
    ptr T
  else:
    index(T, Simd[int])
template `[]`*[T;V,L:static int](x: GpuField[V,T], i: SiteV[L]): auto =
  when V == L:
    x.p[][i[]]
  else:
    let s = i[] div V
    let v = i[] mod V
    x.p[][s][asSimd(v)]

proc `[]=`*[T;V:static int](r: GpuField[V,T], i: SiteV[1], x: auto) =
  let s = i[] div V
  let v = i[] mod V
  r.p[][s][asSimd(v)] = x

proc `[]=`*[T;V:static int](r: GpuField[V,T], i: SiteV[V], x: auto) =
  r.p[][i[]] := x

proc `:=`*(r: GpuField, x: SomeGpuField) =
  for s in gpuSites(r.n, r.V):
    r[s] = x[s]

proc `*`*[X:SomeNumber,Y:SomeGpuField](x: X, y: Y): GpuFieldExpr["*",(X,Y)] =
  result.args[0] = x
  result.args[1] = y
template index*[X:SomeNumber,Y:SomeGpuField](E: typedesc[GpuFieldExpr["*",(X,Y)]], I: typedesc[SiteV]):
  typedesc = X * index(Y,I)
proc `[]`*[X:SomeNumber,Y:SomeGpuField,I:SiteV](e: GpuFieldExpr["*",(X,Y)], i: I): auto =
  type R = index(typeof e, I)
  var r {.noInit.}: R
  r.mul(e.args[0], e.args[1][i])
  r

proc `+`*[X:SomeGpuField,Y:SomeGpuField](x: X, y: Y): GpuFieldExpr["+",(X,Y)] =
  result.args[0] = x
  result.args[1] = y
template index*[X:SomeGpuField,Y:SomeGpuField](E: typedesc[GpuFieldExpr["+",(X,Y)]], I: typedesc[SiteV]):
  typedesc = index(X,I) + index(Y,I)
proc `[]`*[X:SomeGpuField,Y:SomeGpuField,I:SiteV](e: GpuFieldExpr["+",(X,Y)], i: I): auto =
  type R = index(typeof e, I)
  var r {.noInit.}: R
  r.add(e.args[0][i], e.args[1][i])
  r

proc `*`*[X:SomeGpuField,Y:SomeGpuField2](x: X, y: Y): GpuFieldExpr["*",(X,Y)] =
  result.args[0] = x
  result.args[1] = y
template index*[X:SomeGpuField,Y:SomeGpuField2](E: typedesc[GpuFieldExpr["*",(X,Y)]], I: typedesc[SiteV]):
  typedesc = index(X,I) * index(Y,I)
proc `[]`*[X:SomeGpuField,Y:SomeGpuField2,I:SiteV](e: GpuFieldExpr["*",(X,Y)], i: I): auto =
  type R = index(typeof e, I)
  var r {.noInit.}: R
  r.mul(e.args[0][i], e.args[1][i])
  r

type
  CgField*[T] = CpuGpu[T, gpuType(T)]
  CgFld*[C:Field,G:GpuField] = CpuGpu[C,G]
  CgFld2*[C:Field,G:GpuField] = CpuGpu[C,G]
template `[]`*(x: CgFld, i: int): auto = x.cpu[i]
template numberType*[T:CgFld](x: T): typedesc = numberType(T.C)
template `:=`*(r: CgFld, x: SomeNumber) =
  r.cpu := x
template gpuSites*(x: CgFld): auto = gpuSites(x.gpu.n, x.gpu.V)

proc destroy*[C:Field,G:GpuField](x: var CpuGpu[C,G]) =
  when backendIsGpu:
    if x.gpu.p != nil:
      gpuFree(x.gpu.p)
proc newCgField*(c: Field): auto =
  type T = typeof c
  var r: CgField[T]
  r.cpu = c
  r.gpu.n = c.l.nSitesOuter
  when backendIsGpu:
    r.gpu.p = cast[type r.gpu.p](gpuMalloc(r.gpu.bytes))
  else:
    r.gpu.p = cast[type r.gpu.p](addr r.cpu[0])
  r
proc CgColorVector*(lo: Layout): auto = newCgField(lo.ColorVector())
proc CgColorMatrix*(lo: Layout): auto = newCgField(lo.ColorMatrix())
proc CgColorVectorS*(lo: Layout): auto = newCgField(lo.ColorVectorS())
proc CgColorMatrixS*(lo: Layout): auto = newCgField(lo.ColorMatrixS())
proc CgColorVectorD*(lo: Layout): auto = newCgField(lo.ColorVectorD())
proc CgColorMatrixD*(lo: Layout): auto = newCgField(lo.ColorMatrixD())

template toGpu*(g: var GpuField, x: CgFld, cpy: bool) =
  when backendIsGpu:
    if cpy:
      gpuMemCpyToGPU(g.p, addr x.cpu[0], g.bytes)

template getGpu*(x: CgFld, g: GpuField): auto = g

template fromGpu*(x: CgFld, g: GpuField, cpy: bool) =
  when backendIsGpu:
    if cpy:
      gpuMemCpyToCPU(addr x.cpu[0], g.p, g.bytes)

template `*`*(x: SomeNumber, y: CgFld): auto = x * y.cpu
template `*`*(x: CgFld, y: CgFld2): auto = x.cpu * y.cpu
template `+`*(x: SomeField, y: CgFld): auto = x + y.cpu
template `:=`*(x: var CgFld, y: SomeField) =
  x.cpu := y
template `:=`*(x: var CgFld, y: CgFld) =
  x.cpu := y.cpu


proc toGpu*(x: Field): auto =
  var g: gpuType(typeof x)
  g.n = x.l.nSitesOuter
  when backendIsGpu:
    g.p = cast[type g.p](gpuMalloc(g.bytes))
    gpuMemCpyToGPU(g.p, addr x[0], g.bytes)
  else:
    g.p = cast[type g.p](addr x[0])
  g

template getGpu*(x: Field, g: GpuField): auto = g

proc fromGpu*(x: Field, g: GpuField) =
  when backendIsGpu:
    gpuMemCpyToCPU(addr x[0], g.p, g.bytes)

proc toGpu*(g: gpuSeq[GpuField], x: seq[Field]) =
  tic("toGpuSeqField")
  for i in 0..<g.n:
    #g[i] = toGpu(x[i])
    let t = toGpu(x[i])
    gpuMemCpyToGpu(addr g[i], addr t, sizeof(t))
  toc("end")

proc copyFromGpu*(x: seq[Field], g: gpuSeq[GpuField]) =
  tic("copyFromGpuSeqField")
  #for i in 0..<x.len:
  #  var t {.noInit.}: typeof g[i]
  #  gpuMemCpyToCpu(addr t, addr g[i], sizeof(t))
  #  x[i].fromGpu(t)
  var t = newSeq[GpuField](g.n)
  gpuMemCpyToCpu(addr t[0], addr g[0], g.bytes)
  for i in 0..<x.len:
    x[i].fromGpu(t[i])
  toc("end")
