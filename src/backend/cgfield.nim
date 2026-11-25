import qex
import backend/[accel,cpugpu]
import strutils

type
  GpuFieldData*[V:static int, T] = object
    n*: int   # outer sites
    ne*: int  # even outer sites, which come first
    p*: ptr UncheckedArray[T]
  GpuFieldObj*[V:static int, T] = object
    d*: GpuFieldData[V,T]
    pgm*: ptr GpuMem
  GpuField*[V:static int, T] = GpuFieldObj[V,T]
  GpuSubField*[V:static int, T] = object
    ## the outer sites lo..<hi of a field: x.even, x.odd, x.all
    f*: GpuField[V,T]
    lo*, hi*: int
  GpuFieldExpr*[OP:static string, A] = object
    ## an expression of fields and numbers in a kernel, evaluated at a site
    args*: A
  SomeGpuField* = GpuField | GpuSubField | GpuFieldExpr
template n*(x: GpuField): auto = x.d.n
template p*(x: GpuField): auto = x.d.p
proc gpuFlagsExcl*(x: Field, f: set[gmFlags]) =
  gpuMemFlagsExcl(addr x[0], f)
proc gpuFlagsIncl*(x: Field, f: set[gmFlags]) =
  gpuMemFlagsIncl(addr x[0], f)

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

proc `[]=`*[T;V,L:static int](r: GpuField[V,T], i: SiteV[L], x: auto) =
  when V == L:  # a whole site, also the one lane of V = 1
    r.p[][i[]] := x
  else:
    let s = i[] div V
    let v = i[] mod V
    r.p[][s][asSimd(v)] = x

# In a kernel each thread owns a fixed share of the sites of a field, and the
# loops over the field and over its subsets visit only those sites, so a
# kernel may run statements one after the other that read and write the
# same sites.  On a device the share is every nt-th lane from lane t; with
# the CPU backend, whole SIMD sites, a contiguous block of the even outer
# sites and one of the odd ones.

iterator ownedLanes(lo, hi, t, nt: int): int =
  ## lanes t + k nt in lo..<hi
  var i = t
  if i < lo: i += ((lo - i + nt - 1) div nt)*nt
  while i < hi:
    yield i
    i += nt

iterator ownedOuter(n, ne, lo, hi, t, nt: int): int =
  ## outer sites in lo..<hi of the block of thread t of nt in the even sites
  ## 0..<ne and of its block in the odd sites ne..<n
  for i in max(lo, (ne*t) div nt) ..< min(hi, (ne*(t+1)) div nt): yield i
  let no = n - ne
  for i in max(lo, ne + (no*t) div nt) ..< min(hi, ne + (no*(t+1)) div nt): yield i

when backendIsGpu:
  iterator items*[V:static int, T](x: GpuField[V,T]): SiteV[1] =
    ## in a kernel: the sites of x that this thread owns, one lane per index
    for i in ownedLanes(0, x.n*V, gpuThreadNum, gpuNumThreads): yield SiteV[1](i)
  iterator items*[V:static int, T](x: GpuSubField[V,T]): SiteV[1] =
    for i in ownedLanes(x.lo*V, x.hi*V, gpuThreadNum, gpuNumThreads): yield SiteV[1](i)
elif gpuSimd:
  iterator items*[V:static int, T](x: GpuField[V,T]): SiteV[V] =
    ## in a kernel: the sites of x that this thread owns, whole SIMD sites
    for i in ownedOuter(x.n, x.d.ne, 0, x.n, gpuThreadNum, gpuNumThreads): yield SiteV[V](i)
  iterator items*[V:static int, T](x: GpuSubField[V,T]): SiteV[V] =
    for i in ownedOuter(x.f.n, x.f.d.ne, x.lo, x.hi, gpuThreadNum, gpuNumThreads): yield SiteV[V](i)
else:  # lane by lane, through the views of lanes as on a device
  iterator items*[V:static int, T](x: GpuField[V,T]): SiteV[1] =
    for i in ownedOuter(x.n, x.d.ne, 0, x.n, gpuThreadNum, gpuNumThreads):
      for l in 0..<V: yield SiteV[1](V*i + l)
  iterator items*[V:static int, T](x: GpuSubField[V,T]): SiteV[1] =
    for i in ownedOuter(x.f.n, x.f.d.ne, x.lo, x.hi, gpuThreadNum, gpuNumThreads):
      for l in 0..<V: yield SiteV[1](V*i + l)

template gpuSites*(x: GpuField): untyped =
  items((let gpuF = x; gpuF))  # a copy: see gpuRange of backend/accel

# In a generic proc Nim looks up the implicit items of `for e in f` where the
# proc is instantiated, so the generic code of a library writes items(f).
template even*[V:static int, T](x: GpuField[V,T]): GpuSubField[V,T] =
  GpuSubField[V,T](f: x, lo: 0, hi: x.d.ne)
template odd*[V:static int, T](x: GpuField[V,T]): GpuSubField[V,T] =
  GpuSubField[V,T](f: x, lo: x.d.ne, hi: x.n)
template all*[V:static int, T](x: GpuField[V,T]): GpuSubField[V,T] =
  GpuSubField[V,T](f: x, lo: 0, hi: x.n)
template `[]`*(x: GpuSubField, i: SiteV): untyped = x.f[i]

# Expressions of fields and numbers in a kernel: x[i] of an expression is the
# value at site index i, computed with the operators of the site types.

template atSite*(x: SomeNumber; i: SiteV): untyped = x
template atSite*(x: GpuField|GpuSubField|GpuFieldExpr; i: SiteV): untyped = x[i]

template `[]`*[OP:static string, A](e: GpuFieldExpr[OP,A]; i: SiteV): untyped =
  when OP == "+": atSite(e.args[0], i) + atSite(e.args[1], i)
  elif OP == "-": atSite(e.args[0], i) - atSite(e.args[1], i)
  elif OP == "*": atSite(e.args[0], i) * atSite(e.args[1], i)
  elif OP == "neg": -atSite(e.args[0], i)
  else: atSite(e.args[0], i).adj

# Nim binds a type class to one type per signature: the second operand of a
# binary operator has classes of its own names, so a matrix field times a
# vector field matches.
type
  GpuField2*[V:static int, T] = GpuFieldObj[V,T]
  GpuSubField2*[V:static int, T] = GpuSubField[V,T]
  GpuFieldExpr2*[OP:static string, A] = GpuFieldExpr[OP,A]
  SomeGpuField2* = GpuField2 | GpuSubField2 | GpuFieldExpr2

proc `+`*[X:SomeGpuField, Y:SomeGpuField2](x: X; y: Y): GpuFieldExpr["+",(X,Y)] = result.args = (x, y)
proc `+`*[X:SomeGpuField, Y:SomeNumber](x: X; y: Y): GpuFieldExpr["+",(X,Y)] = result.args = (x, y)
proc `+`*[X:SomeNumber, Y:SomeGpuField](x: X; y: Y): GpuFieldExpr["+",(X,Y)] = result.args = (x, y)
proc `-`*[X:SomeGpuField, Y:SomeGpuField2](x: X; y: Y): GpuFieldExpr["-",(X,Y)] = result.args = (x, y)
proc `-`*[X:SomeGpuField, Y:SomeNumber](x: X; y: Y): GpuFieldExpr["-",(X,Y)] = result.args = (x, y)
proc `-`*[X:SomeNumber, Y:SomeGpuField](x: X; y: Y): GpuFieldExpr["-",(X,Y)] = result.args = (x, y)
proc `*`*[X:SomeGpuField, Y:SomeGpuField2](x: X; y: Y): GpuFieldExpr["*",(X,Y)] = result.args = (x, y)
proc `*`*[X:SomeGpuField, Y:SomeNumber](x: X; y: Y): GpuFieldExpr["*",(X,Y)] = result.args = (x, y)
proc `*`*[X:SomeNumber, Y:SomeGpuField](x: X; y: Y): GpuFieldExpr["*",(X,Y)] = result.args = (x, y)
proc `-`*[X:SomeGpuField](x: X): GpuFieldExpr["neg",(X,)] = result.args = (x,)
proc adj*[X:SomeGpuField](x: X): GpuFieldExpr["adj",(X,)] = result.args = (x,)

template lanes(x: GpuField): int = x.n*x.V
template lanes(x: GpuSubField): int = (x.hi - x.lo)*x.f.V

var viewSum: GpuSum[float]  # of the sums of field views
var viewSumN = 0  # the threads it has room for, the largest launch so far

template siteSum(x: GpuField|GpuSubField; e, body: untyped): float =
  ## outside kernels and threads blocks: the sum of body over the sites e of
  ## x and the ranks, in an order fixed by the thread count of the kernel
  ## (GpuSum)
  block:
    let gpuX = x
    let gpuL = lanes(gpuX)
    if viewSumN < gpuL:
      if viewSumN > 0: viewSum.free
      viewSum = newGpuSum[float](gpuL)
      viewSumN = gpuL
    let gpuS = viewSum
    onGpu(gpuL):
      var s = 0.0
      for e in items(gpuX): s += simdSum(body)
      gpuS.reduce s
    var r = gpuS.value
    getDefaultComm().allReduce(r)
    r

proc norm2*(x: GpuField|GpuSubField): float =
  ## outside kernels and threads blocks: |x|^2 over the sites and the ranks
  siteSum(x, e, norm2(x[e]))

proc redot*(x: GpuField|GpuSubField; y: GpuField2|GpuSubField2): float =
  ## outside kernels and threads blocks: Re (x, y) over the sites of x and
  ## the ranks
  siteSum(x, e, redot(x[e], y[e]))

template assignOp(op: untyped) =
  template op*(r: GpuField|GpuSubField; x: SomeGpuField|SomeNumber) =
    ## in a kernel: the statement at the sites of r that this thread owns
    block:
      let gpuR = r  # a copy: see gpuRange of backend/accel
      let gpuX = x
      for gpuE in items(gpuR):
        op(gpuR[gpuE], atSite(gpuX, gpuE))
assignOp(`:=`)
assignOp(`+=`)
assignOp(`-=`)
assignOp(`*=`)


proc toGpu*(g: var GpuField, c: Field) =
  g.d.n = c.l.nSitesOuter
  g.d.ne = c.l.nEvenOuter
  when not beSharedMem:
    tic("toGpuField")
    pushGpuMemTag("Field"&capitalizeAscii(displayName(typeof(g.T))))
    let pgm = getGpuMem(addr c[0], g.bytes)
    popGpuMemTag()
    g.d.p = cast[type g.p](pgm.p)
    g.pgm = pgm
    toc("getGpuMem")
    pgm.copyIn(addr c[0])
    toc("copyIn",flops=pgm.bytes)
  else:
    g.d.p = cast[type g.p](addr c[0])

proc gpuFieldView*[V:static int](lo: Layout[V]; T: typedesc; p: pointer): auto =
  ## a field view of device memory p that holds a field of lo with element
  ## type T, as gpuMalloc gives it
  var g: GpuField[V, gpuType(T)]
  g.d.n = lo.nSitesOuter
  g.d.ne = lo.nEvenOuter
  g.d.p = cast[typeof(g.d.p)](p)
  g

proc toGpu*(c: Field): auto =
  var g: gpuType(typeof c)
  g.toGpu(c)
  g

template getGpu*(x: Field, g: GpuField): auto = g

proc fromGpu*(x: Field) =
  when not beSharedMem:
    tic("fromGpuField")
    let pgm = getGpuMem(addr x[0])
    pgm.copyOut(addr x[0])
    toc("copyOut",flops=pgm.bytes)

template fromGpu*(x: Field, g: GpuField) = fromGpu(x)

proc toGpu*(g: GpuSeq[GpuField], c: seq[Field], pgm: ptr GpuMem) =
  if pgm.needsCopyIn:
    tic("toGpuSeqFieldNew")
    var t = newSeq[typeof(g[0])](g.n)
    for i in 0..<g.n:
      t[i].toGpu(c[i])
    toc("toGpu")
    pgm.copyIn(addr t[0])
    toc("copyIn")
  else:
    tic("toGpuSeqFieldOld")
    for i in 0..<g.n:
      discard toGpu(c[i])
    toc("toGpu")

proc fromGpu*(c: var seq[Field], g: GpuSeq[GpuField], pgm: ptr GpuMem) =
  if pgm.needsCopyOut:
    tic("fromGpuSeqFieldNew")
    var t = newSeq[typeof(g[0])](g.n)
    pgm.copyOut(addr t[0])
    toc("copyOut")
    for i in 0..<g.n:
      c[i].fromGpu(t[i])
    toc("fromGpu")
  else:
    tic("fromGpuSeqFieldOld")
    for i in 0..<g.n:
      fromGpu(c[i])
    toc("fromGpu")

#[
proc copyFromGpu*(x: seq[Field], g: GpuSeq[GpuField]) =
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
]#






type
  CgFieldObj*[T] = CpuGpu[T, gpuType(T)]
  CgField*[T] = ref CpuGpu[T, gpuType(T)]
  CgFld*[C:Field,G:GpuField] = ref CpuGpu[C, G]
  CgFld2*[C:Field,G:GpuField] = ref CpuGpu[C, G]
template `[]`*(x: CgFld, i: int): auto = x.cpu[i]
template numberType*[T:CgFld](x: typedesc T): typedesc = numberType(T.C)
#template numberType*[T:CgFld](x: T): typedesc = numberType(T.C)
template `:=`*(r: CgFld, x: SomeNumber) =
  r.cpu := x
template gpuSites*(x: CgFld): untyped = items(x.gpu)

template displayName*[T](t: typedesc[Color[T]]): string =
  "Color" & displayName(T)
template displayName*[N:static int,T](t: typedesc[VectorArray[N,T]]): string =
  "Vector" & $N
template displayName*[N,M:static int,T](t: typedesc[MatrixArray[N,M,T]]): string =
  "Matrix" & $N & "," & $M
template displayName*[V:static int,T](t: typedesc[GpuField[V,T]]): string =
  "Field" & displayName(T)

proc newCgField*(c: Field): auto =
  ## The device copy of c goes with c's storage, which other CgFields and
  ## kernels may share.
  type C = typeof c
  type CG = CgFieldObj[C]
  var r: ref CG
  r.new
  r.cpu = c
  r.gpu.d.n = c.l.nSitesOuter
  r.gpu.d.ne = c.l.nEvenOuter
  when not beSharedMem:
    #r.gpu.p = cast[type r.gpu.p](gpuMalloc(r.gpu.bytes))
    pushGpuMemTag("Field"&capitalizeAscii(displayName(typeof(r.gpu.T))))
    let pgm = getGpuMem(addr c[0], r.gpu.bytes)
    popGpuMemTag()
    r.gpu.d.p = cast[type r.gpu.p](pgm.p)
    r.gpu.pgm = pgm
    #echo "malloc: ", cast[int](r.gpu.p)
  else:
    r.gpu.d.p = cast[type r.gpu.p](addr r.cpu[0])
  r
proc CgColorVector*(lo: Layout): auto = newCgField(lo.ColorVector())
proc CgColorMatrix*(lo: Layout): auto = newCgField(lo.ColorMatrix())
proc CgColorVectorS*(lo: Layout): auto = newCgField(lo.ColorVectorS())
proc CgColorMatrixS*(lo: Layout): auto = newCgField(lo.ColorMatrixS())
proc CgColorVectorD*(lo: Layout): auto = newCgField(lo.ColorVectorD())
proc CgColorMatrixD*(lo: Layout): auto = newCgField(lo.ColorMatrixD())

template toGpu*(g: GpuField, x: CgFld, cpy: bool) =
  when not beSharedMem:
    g.pgm.touch
    if cpy:
      tic("toGpuCgField")
      gpuMemCpyToGPU(g.p, addr x.cpu[0], g.bytes)
      g.pgm.wasCopiedIn
      toc("gpuMemCpyToGPU")

template getGpu*(x: CgFld, g: GpuField): auto = g

template fromGpu*(x: CgFld, g: GpuField, cpy: bool) =
  when not beSharedMem:
    if cpy:
      tic("fromGpuCgField")
      gpuMemCpyToCPU(addr x.cpu[0], g.p, g.bytes)
      g.pgm.wasCopiedOut
      toc("gpuMemCpyToCPU")

template `*`*(x: SomeNumber, y: CgFld): auto = x * y.cpu
template `*`*(x: CgFld, y: CgFld2): auto = x.cpu * y.cpu
template `+`*(x: SomeField, y: CgFld): auto = x + y.cpu
template `:=`*(x: var CgFld, y: SomeField) =
  x.cpu := y
template `:=`*(x: var CgFld, y: CgFld) =
  x.cpu := y.cpu
