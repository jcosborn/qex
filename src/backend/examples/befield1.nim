import qex
import physics/qcdTypes
import backend/[accel,cpugpu]
import bench/commonBench
#import strformat
import macros
import base/metaUtils
import parseUtils
import sequtils, strutils

template isWrapper*(x: typedesc[array]): bool = false

template `*`*[T:SomeNumber](x: typedesc[T], y: typedesc[T]): typedesc =
  T
template `*`*[N:static int,X,Y](x: typedesc[array[N,X]], y: typedesc[array[N,Y]]): typedesc =
  array[N, X * Y]
template `*`*[X:SomeNumber,N:static int,T](x: typedesc[X], y: typedesc[array[N,T]]): typedesc =
  array[N, X * T]
template `*`*[X:SomeNumber,Y:Simd](x: typedesc[X], y: typedesc[Y]): typedesc =
  asSimd(X * Y[])
template `*`*[X:SomeNumber,Y:ComplexObj](x: typedesc[X], y: typedesc[Y]): typedesc =
  ComplexObj[X*Y.TR,X*Y.TI]
template `*`*[X:SomeNumber,Y:AsComplex](x: typedesc[X], y: typedesc[Y]): typedesc =
  asComplex(X * Y[])
template `*`*[X:SomeNumber,Y:VectorArrayObj](x: typedesc[X], y: typedesc[Y]): typedesc =
  VectorArrayObj[Y.I, X * Y.T]
template `*`*[X:SomeNumber,Y:AsVector](x: typedesc[X], y: typedesc[Y]): typedesc =
  asVector(X * Y[])
template `*`*[X:SomeNumber,Y:Color](x: typedesc[X], y: typedesc[Y]): typedesc =
  asColor(X * Y[])
template `*`*[X:SomeNumber,Y:Color](x: typedesc[X], y: typedesc[ptr Y]): typedesc =
  #static: echo $Y, "  ", $(Y[])
  asColor(X * Y[])

template `*`*[X:ComplexObj,Y:ComplexObj2](x: typedesc[X], y: typedesc[Y]): typedesc =
  ComplexObj[X.TR*Y.TR-X.TI*Y.TI,X.TR*Y.TI+X.TI*Y.TR]
template `*`*[X:AsComplex,Y:AsComplex2](x: typedesc[X], y: typedesc[Y]): typedesc =
  asComplex(X[] * Y[])
template `*`*[N,M:static int;X,Y](x: typedesc[MatrixArrayObj[N,M,X]], y: typedesc[VectorArrayObj[M,Y]]): typedesc =
  VectorArrayObj[N, X * Y]
template `*`*[X:AsMatrix,Y:AsVector](x: typedesc[X], y: typedesc[Y]): typedesc =
  asVector(X[] * Y[])

template `+`*[T:SomeNumber](x: typedesc[T], y: typedesc[T]): typedesc =
  T
template `-`*[T:SomeNumber](x: typedesc[T], y: typedesc[T]): typedesc =
  T

template `+`*[N:static int,X,Y](x: typedesc[array[N,X]], y: typedesc[array[N,Y]]): typedesc =
  array[N, X + Y]
template `-`*[N:static int,X,Y](x: typedesc[array[N,X]], y: typedesc[array[N,Y]]): typedesc =
  array[N, X - Y]
template `+`*[X:Simd,Y:Simd](x: typedesc[X], y: typedesc[Y]): typedesc =
  asSimd(X[] + Y[])
template `-`*[X:Simd,Y:Simd](x: typedesc[X], y: typedesc[Y]): typedesc =
  asSimd(X[] - Y[])
template `-`*[X:ComplexObj,Y:ComplexObj](x: typedesc[X], y: typedesc[Y]): typedesc =
  ComplexObj[X.TR-Y.TR,X.TI-Y.TI]
template `+`*[X:ComplexObj,Y:ComplexObj](x: typedesc[X], y: typedesc[Y]): typedesc =
  ComplexObj[X.TR+Y.TR,X.TI+Y.TI]
template `+`*[X:AsComplex,Y:AsComplex](x: typedesc[X], y: typedesc[Y]): typedesc =
  asComplex(X[] + Y[])
template `+`*[I:static int,X,Y](x: typedesc[VectorArrayObj[I,X]], y: typedesc[VectorArrayObj[I,Y]]):
  typedesc = VectorArrayObj[I, X + Y]
template `+`*[X:AsVector,Y:AsVector](x: typedesc[X], y: typedesc[Y]): typedesc =
  asVector(X[] + Y[])
template `+`*[X:Color,Y:Color2](x: typedesc[X], y: typedesc[Y]): typedesc =
  asColor(X[] + Y[])
template `+`*[X:Color,Y:Color2](x: typedesc[X], y: typedesc[ptr Y]): typedesc =
  asColor(X[] + Y[])

template `*`*[X:Simd,Y:Simd](x: typedesc[X], y: typedesc[Y]): typedesc =
  asSimd(X[] * Y[])
template `*`*[X:Color,Y:Color2](x: typedesc[X], y: typedesc[Y]): typedesc =
  asColor(X[] * Y[])
template `*`*[X:Color,Y:Color2](x: typedesc[X], y: typedesc[ptr Y]): typedesc =
  asColor(X[] * Y[])
template `*`*[X:Color,Y:Color2](x: typedesc[ptr X], y: typedesc[ptr Y]): typedesc =
  asColor(X[] * Y[])

template toSingleImpl*[N:static int](x: array[N,float32]): auto = x
template toDoubleImpl*[N:static int](x: array[N,float]): auto = x
proc `:=`*(r: var Color, x: ptr Color2) =
  r := x[]

proc mul*[R:array,Y:array](r: var R, x: SomeNumber, y: Y) =
  for i in 0..<r.len:
    r[i] = x * y[i]
proc mul*[R,X,Y:array](r: var R, x: X, y: Y) =
  for i in 0..<r.len:
    r[i] = x[i] * y[i]
proc mul*(r: var Color, x: SomeNumber, y: ptr Color2) =
  mul(r[], x, y[][])
proc mul*(r: var Color, x: ptr Color2, y: ptr Color3) =
  mul(r[], x[][], y[][])

template add*[R:SomeNumber,X:SomeNumber,Y:AsFloat](r: R, x: X, y: Y) =
  r = x + eval(y)
  #r = x
template mul*[R:SomeNumber,X:AsFloat,Y:AsFloat](r: R, x: X, y: Y) =
  r = eval(x) * eval(y)
template imadd*[R:SomeNumber,X:AsFloat,Y:AsFloat](r: R, x: X, y: Y) =
  r += eval(x) * eval(y)
template imsub*[R:SomeNumber,X:AsFloat,Y:AsFloat](r: R, x: X, y: Y) =
  r -= eval(x) * eval(y)
proc imadd*[R,X,Y:array](r: var R, x: X, y: Y) =
  for i in 0..<r.len:
    r[i] += x[i] * y[i]
proc imsub*[R,X,Y:array](r: var R, x: X, y: Y) =
  for i in 0..<r.len:
    r[i] -= x[i] * y[i]

proc add*[R:array,X:array,Y:array](r: var R, x: X, y: Y) =
  for i in 0..<r.len:
    r[i] = x[i] + y[i]
proc add*(r: var Color, x: Color, y: ptr Color) =
  add(r[], x[], y[][])


type
  SiteV[V:static int] = distinct int
template `[]`*(x: SiteV): int = int(x)

when backendIsGpu:
  iterator gpuSites(n:int, V:static int): SiteV[1] =
    for s in gpuRange(n*V):
      yield SiteV[1](s)
else:
  iterator gpuSites(n:int, V:static int): SiteV[V] =
    for s in gpuRange(n):
      yield SiteV[V](s)

type
  #FieldProxy[T] = object
  #  f: T
  GpuFieldObj*[V:static int, T] = object
    n: int
    p: ptr UncheckedArray[T]
  GpuFieldExpr*[OP:static string, A] = object
    args: A
  GpuFieldExpr2*[OP:static string, A] = GpuFieldExpr[OP,A]
  GpuField*[V:static int, T] = GpuFieldObj[V,T]
  GpuField2*[V:static int, T] = GpuFieldObj[V,T]
  SomeGpuField = GpuField | GpuFieldExpr
  SomeGpuField2 = GpuField2 | GpuFieldExpr2
template gpuSites(x: GpuField): auto = gpuSites(x.n, x.V)

proc bytes*[T:GpuField](x: T): int = x.n * sizeof(T.T)

template index[T;V,L:static int](X: typedesc[GpuField[V,T]], I: typedesc[SiteV[L]]): typedesc =
  when V == L:
    ptr T
  else:
    index(T, Simd[int])
proc `[]`[T;V,L:static int](x: GpuField[V,T], i: SiteV[L]): auto =
  when V == L:
    addr x.p[][i[]]
  else:
    let s = i[] div V
    let v = i[] mod V
    x.p[][s][asSimd(v)]

proc `[]=`[T;V:static int](r: GpuField[V,T], i: SiteV[1], x: auto) =
  let s = i[] div V
  let v = i[] mod V
  r.p[][s][asSimd(v)] = x

proc `[]=`[T;V:static int](r: GpuField[V,T], i: SiteV[V], x: auto) =
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

#template toGpuField(t: typedesc[SimdD16Obj]): typedesc =
#  array[16,float]
template toGpuField[T](t: typedesc[Simd[T]]): typedesc =
  Simd[array[T.numNumbers,T.numberType]]
template toGpuField[T](t: typedesc[ComplexType[T]]): typedesc =
  ComplexType[toGpuField(T)]
template toGpuField[N:static int; T](t: typedesc[VectorArray[N,T]]): typedesc =
  VectorArray[N,toGpuField(T)]
template toGpuField[N,M:static int; T](t: typedesc[MatrixArray[N,M,T]]): typedesc =
  MatrixArray[N,M,toGpuField(T)]
template toGpuField[T](t: typedesc[Color[T]]): typedesc =
  Color[toGpuField(T)]
template toGpuField[V:static int, T](t: typedesc[Field[V,T]]): typedesc =
  GpuField[V,toGpuField(T)]

type
  CgField*[T] = CpuGpu[T, toGpuField(T)]
  CgFld[C:Field,G:GpuField] = CpuGpu[C,G]
  CgFld2[C:Field,G:GpuField] = CpuGpu[C,G]
template `[]`*(x: CgFld, i: int): auto = x.cpu[i]
template numberType*[T:CgFld](x: T): typedesc = numberType(T.C)
template `:=`*(r: CgFld, x: SomeNumber) =
  r.cpu := x
template gpuSites(x: CgFld): auto = gpuSites(x.gpu.n, x.gpu.V)

proc destroy*[C:Field,G:GpuField](x: var CpuGpu[C,G]) =
  when backendIsGpu:
    if x.gpu.p != nil:
      gpuFree(x.gpu.p)
proc newCgField(c: Field): auto =
  type T = typeof c
  var r: CgField[T]
  r.cpu = c
  r.gpu.n = c.l.nSitesOuter
  when backendIsGpu:
    r.gpu.p = cast[type r.gpu.p](gpuMalloc(r.gpu.bytes))
  else:
    r.gpu.p = cast[type r.gpu.p](addr r.cpu[0])
  r
proc CgColorVector(lo: Layout): auto = newCgField(lo.ColorVector())
proc CgColorMatrix(lo: Layout): auto = newCgField(lo.ColorMatrix())
proc CgColorVectorS(lo: Layout): auto = newCgField(lo.ColorVectorS())
proc CgColorMatrixS(lo: Layout): auto = newCgField(lo.ColorMatrixS())
proc CgColorVectorD(lo: Layout): auto = newCgField(lo.ColorVectorD())
proc CgColorMatrixD(lo: Layout): auto = newCgField(lo.ColorMatrixD())

template toGpu*(g: var GpuField, x: CgFld, cpy: bool) =
  when backendIsGpu:
    if cpy:
      gpuMemCpyToGPU(g.p, addr x.cpu[0], g.bytes)

template getGpu*(x: CgFld, g: GpuField): auto = g

template fromGpu*(x: CgFld, g: GpuField, cpy: bool) =
  when backendIsGpu:
    if cpy:
      gpuMemCpyToCPU(addr x.cpu[0], g.p, g.bytes)

template `*`(x: SomeNumber, y: CgFld): auto = x * y.cpu
template `*`(x: CgFld, y: CgFld2): auto = x.cpu * y.cpu
template `+`(x: SomeField, y: CgFld): auto = x + y.cpu
template `:=`(x: var CgFld, y: SomeField) =
  x.cpu := y
template `:=`(x: var CgFld, y: CgFld) =
  x.cpu := y.cpu

macro exp2string(x:untyped):auto =
  #echo x.treerepr
  #echo x.repr
  #var s = repr fixBracket symToIdent x
  var s = repr symToIdent x
  #echo x.repr
  #echo s
  let n = skipWhitespace(s)
  result = newlit s[n..^1].replace("\n"," ")
  #echo result

template check(r,eqn: untyped) =
  r.cpuReadWrite
  block:
    var s = r.cpu.newOneOf
    threads:
      r := 0
      threadBarrier()
      eqn
      threadBarrier()
      s := r.cpu
      threadBarrier()
      r := 0
    onGpu:
      eqn
    threads:
      s -= r.cpu
      echo s.norm2
  r.cpuUnused

var bs = newSeq[string](0)
var br = newSeq[array[3,string]](0)
var repin = 0
template bench(r,fps,bps,mm,eqn: untyped) =
  check(r, eqn)
  block:
    let N = r.cpu.l.nSitesOuter
    let V = r.cpu.l.V
    let vol = N*V
    var b = newBench()
    case repin
    of 0:
      benchSingle(b):
        let nrep = b.nrep
        onGpu(vol):
          for rep in 0..<nrep:
            eqn
    of 1:
      benchSingle(b):
        let nrep = b.nrep
        for rep in 0..<nrep:
          onGpu(vol):
            eqn
    else:
      r.cpuReadWrite
      benchSingle(b):
        let nrep = b.nrep
        for rep in 0..<nrep:
          onGpu(vol):
            eqn
      r.cpuUnused
    let memMB = vol * mm.float / (1024.0*1024.0)
    let gb = vol * bps.float * b.perNs
    let gf = vol * fps.float * b.perNs
    #echo &"{memMB:8.3f} {gb:8.3f} {gf:8.3f}"
    #echo memMB|(8,-3), " ", gb|(9,-3), " ", gf|(8,-3)
    bs.add exp2string(eqn)
    br.add [memMB|(8,-3), gb|(9,-3), gf|(8,-3)]

proc test(lat:auto, double:static bool=false) =
  let minl = intParam("minl",0)
  if lat[0] < minl: return
  let maxl = intParam("maxl",0)
  if maxl>0:
    if lat[0] > maxl: return
  var lo = newLayout(lat)

  when double:
    type Real = float
    template newCV(): auto = lo.CgColorVectorD()
    template newCM(): auto = lo.CgColorMatrixD()
    #template newHF: untyped = lo.HalfFermion()
    #template newDF: untyped = lo.DiracFermion()
  else:
    type Real = float32
    template newCV(): auto = lo.CgColorVectorS()
    template newCM(): auto = lo.CgColorMatrixS()

  var v1 = newCV()
  var v2 = newCV()
  #var v3 = newCV()
  var m1 = newCM()
  #var m2 = newCM()
  #var m3 = newCM()
  #var m4 = newCM()
  #var m5 = newCM()
  #var h1 = newHF()
  #var d1 = newDF()
  #var d2 = newDF()
  const nc = v1[0].len
  let sf = sizeof(numberType(v1[0]))
  let nc2 = 2*nc
  let vb = nc2*sf
  let mb = nc*vb
  let mvf = (2*nc2-1)*nc2
  echo "Float type: ", $(v1.numberType)
  macro xfer(x: varargs[untyped]): auto =
    result = newNimNode(nnkStmtList)
    for n in x:
      result.add quote do:
        `n`.copyToGpu
        `n`.cpuUnused
  proc reset =
    threads:
      m1 := 2
      v1 := 1
      v2 := 0
      #v3 := 0
    xfer(m1, v1, v2)
  #echo "done setup"
  #var nbench = 0
  template ioGpu:untyped = declared(inOnGpu)

  reset()
  template job1(v2,v1:untyped) =
    when ioGpu:
      #v2 := Real(0.5)*v2 + v1
      for s in gpuRange(v2.n):
        v2.p[s] = Real(0.5)*v2.p[s] + v1.p[s]
      #for s in gpuRange(v2.n):
      #  for ic in 0..<v2.p[0].getNc:
      #    v2.p[s][ic] = Real(0.5)*v2.p[s][ic] + v1.p[s][ic]
    else:
      v2.cpu := Real(0.5)*v2.cpu + v1.cpu
  bench(v2, 2*nc2, 3*vb, 2*vb):
    job1(v2, v1);  ## loop over outer sites

  reset()
  bench(v2, 2*nc2, 3*vb, 2*vb):
    v2 := Real(0.5)*v2 + v1

  reset()
  bench(v2, mvf, mb+2*vb, mb+2*vb):
    v2 := m1 * v1

qexInit()
echo "rank ", myRank, "/", nRanks
threads:
  echo "thread ", threadNum, "/", numThreads

var bss = newSeq[typeof bs](0)
var brs = newSeq[typeof br](0)
template runtests(double:static bool) =
  bs.setLen(0)
  br.setLen(0)
  test([4,4,4,4], double)
  test([8,8,8,8], double)
  test([12,12,12,12], double)
  test([16,16,16,16], double)
  test([24,24,24,24], double)
  test([32,32,32,32], double)
  test([40,40,40,40], double)
  #test([48,48,48,48], double)
  bss.add bs
  brs.add br
template runt(double:static bool) =
  repin = 0
  runtests(double)
  repin = 1
  runtests(double)
  repin = 2
  runtests(double)
runt(false)
let bss32 = bss
let brs32 = brs
runt(true)

proc echoResult(p: string, s: seq, r: auto) =
  let sd = s[0].deduplicate
  for t in 0..<sd.len:
    echo "      MB      GB/s     ", p, "  ", sd[t]
    for i in 0..<r[0].len:
      if s[0][i] == sd[t]:
        echo r[0][i][0], " ", r[0][i][1], " ", r[1][i][1], " ", r[2][i][1]
echoResult("float32", bss32, brs32)
echoResult("float64", bss, brs)

qexFinalize()
