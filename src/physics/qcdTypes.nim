import base
export base
#import alignedMem
#export alignedMem
import layout
#import basicOps
#export basicOps
import simd
export simd
import macros
#import metaUtils
#import threading
#import strutils
#import globals
#import comms
import field
export field
#import complexConcept
#export complexConcept
#import complexType
#export complexType
#import matrixConcept
#export matrixConcept
import maths
export maths
import maths/types
export types
import base/wrapperTypes
import color
export color
import spin
export spin
getOptimPragmas()

const nc {.intDefine.} = getDefaultNc()
static: echo "INFO: qcdTypes default Nc ", nc
const ns = 4
const nh = 2
setType(Svec0, "SimdS" & $VLEN)
setType(Dvec0, "SimdD" & $VLEN)
type
  SDvec = Svec0 | Dvec0
  #SComplex* = AsComplex[tuple[re,im:float32]]
  #SComplexV* = AsComplex[tuple[re,im:Svec0]]
  SComplex* = ComplexType[float32]
  SComplexV* = ComplexType[Svec0]
  #SComplex* = Complex[float32,float32]
  #SComplexV* = Complex[Svec0,Svec0]
  SColorVector* = Color[VectorArray[nc,SComplex]]
  SColorVectorV* = Color[VectorArray[nc,SComplexV]]
  SColorMatrix* = Color[MatrixArray[nc,nc,SComplex]]
  SColorMatrixV* = Color[MatrixArray[nc,nc,SComplexV]]

  SLatticeReal* = Field[1,float32]
  #SLatticeRealG*[V:static[int]] = Field[V,Svec0]
  SLatticeRealV* = Field[VLEN,Svec0]
  SLatticeComplex* = Field[1,SComplex]
  SLatticeComplexV* = Field[VLEN,SComplexV]
  SLatticeColorVector* = Field[1,SColorVector]
  SLatticeColorVectorV* = Field[VLEN,SColorVectorV]
  SLatticeColorMatrix* = Field[1,SColorMatrix]
  SLatticeColorMatrixV* = Field[VLEN,SColorMatrixV]

  #DComplex* = AsComplex[tuple[re,im:float64]]
  #DComplexV* = AsComplex[tuple[re,im:Dvec0]]
  DComplex* = ComplexType[float64]
  DComplexV* = ComplexType[Dvec0]
  #DComplex* = Complex[float64,float64]
  #DComplexV* = Complex[Dvec0,Dvec0]
  DColorVector* = Color[VectorArray[nc,DComplex]]
  DColorVectorV* = Color[VectorArray[nc,DComplexV]]
  DColorMatrix* = Color[MatrixArray[nc,nc,DComplex]]
  DColorMatrixV* = Color[MatrixArray[nc,nc,DComplexV]]

  ColorMatrixN*[n:static[int],T] = Color[MatrixArray[n,n,T]]

  DLatticeReal* = Field[1,float64]
  DLatticeRealV* = Field[VLEN,Dvec0]
  DLatticeComplex* = Field[1,DComplex]
  DLatticeComplexV* = Field[VLEN,DComplexV]
  DLatticeColorVector* = Field[1,DColorVector]
  DLatticeColorVectorV* = Field[VLEN,DColorVectorV]
  DLatticeColorMatrix* = Field[1,DColorMatrix]
  DLatticeColorMatrixV* = Field[VLEN,DColorMatrixV]

  SDiracFermion* = Spin[VectorArray[ns,SColorVector]]
  SDiracFermionV* = Spin[VectorArray[ns,SColorVectorV]]
  SHalfFermion* = Spin[VectorArray[nh,SColorVector]]
  SHalfFermionV* = Spin[VectorArray[nh,SColorVectorV]]
  SLatticeDiracFermion* = Field[1,SDiracFermion]
  SLatticeDiracFermionV* = Field[VLEN,SDiracFermionV]
  SLatticeHalfFermion* = Field[1,SHalfFermion]
  SLatticeHalfFermionV* = Field[VLEN,SHalfFermionV]

  DDiracFermion* = Spin[VectorArray[ns,DColorVector]]
  DDiracFermionV* = Spin[VectorArray[ns,DColorVectorV]]
  DHalfFermion* = Spin[VectorArray[nh,DColorVector]]
  DHalfFermionV* = Spin[VectorArray[nh,DColorVectorV]]
  DLatticeDiracFermion* = Field[1,DDiracFermion]
  DLatticeDiracFermionV* = Field[VLEN,DDiracFermionV]
  DLatticeHalfFermion* = Field[1,DHalfFermion]
  DLatticeHalfFermionV* = Field[VLEN,DHalfFermionV]

#template simdLength*(x:typedesc[SColorMatrixV]):untyped =
#  mixin simdLength
#  simdLength(Svec0)
#template simdLength*(x:typedesc[SColorVectorV]):untyped = simdLength(Svec0)
#template simdLength*(x:SColorVectorV):untyped = simdLength(Svec0)
#template simdLength*(x:SColorMatrixV):untyped = simdLength(Svec0)
#template simdLength*(x:typedesc[DColorMatrixV]):untyped = simdLength(Dvec0)
#template simdLength*(x:typedesc[DColorVectorV]):untyped =
#  mixin simdLength
#  simdLength(Dvec0)
#template simdLength*(x:DColorVectorV):untyped = simdLength(Dvec0)
#template simdLength*(x:DColorMatrixV):untyped = simdLength(Dvec0)
#template simdLength*(x:AsComplex):untyped = simdLength(x.re)
#template simdLength*(x:Complex):untyped = simdLength(x.re)
#template simdLength*(x:AsMatrix):untyped = simdLength(x[0,0])
#template simdLength*(x:AsVector):untyped = simdLength(x[0])

#template nVectors(x:SColorVectorV):untyped = 2*nc
#template nVectors(x:SColorMatrixV):untyped = 2*nc*nc
#template nVectors(x:DColorVectorV):untyped = 2*nc
#template nVectors(x:DColorMatrixV):untyped = 2*nc*nc
template nVectors(x:Svec0):untyped = 1
template nVectors(x:Dvec0):untyped = 1
template nVectors(x:AsComplex):untyped = 2*nVectors(x.re)
#template nVectors(x:DComplexV):untyped = 2*nVectors(x.re)
template nVectors*(x:AsVector):untyped = x.len*nVectors(x[0])
template nVectors*(x:AsMatrix):untyped = x.nrows*x.ncols*nVectors(x[0,0])

template simdType*(x:tuple):untyped = simdType(x[0])
template simdType*(x:array):untyped = simdType(x[x.low])
template simdType*(x:AsComplex):untyped = simdType(x[])
template simdType*(x:ComplexObj):untyped = simdType(x.re)
#template simdType*(x:DComplexV):untyped = simdType(x.re)
template simdType*(x:AsVector):untyped = simdType(x[0])
#template simdType*(x:AsMatrix):untyped = simdType(x[])
template simdType*(x:AsMatrix):untyped = simdType(x[0,0])


template `*`*(x: Color, y: Spin): untyped =
  staticTraceReturn timesColorSpin:
    asSpin(x * y[])


#template trace*(x:SComplexV):untyped = x
#proc simdSum*(x:SComplexV):SComplex = complexConcept.map(result, simdSum, x)
#proc simdSum*(x:DComplexV):DComplex = complexConcept.map(result, simdSum, x)
#template simdSum*(x:ToDouble):untyped = toDouble(simdSum(x[]))
#template simdSum*(x: AsComplex): untyped =
#  let tSimdSum = simdSum(x[])
#  asComplex(tSimdSum)
#template simdSum*(x:Complex):untyped = simdSum(x[])
#template simdSum*(xx:tuple):untyped =
#  lets(x,xx):
#    map(x, simdSum)
#template simdSum*(x: ComplexObj): untyped =
#  #mapComplexObj(x, simdReduce)
#  map(x, simdReduce)
#template rankSum*(x:AsComplex) =
#  #mixin qmpSum
#  rankSum(x[])
#template `/`*(x:SComplex|DComplex; y:SomeNumber):untyped =
#  mixin divd
#  var r{.noInit.}:type(x)
#  echoType: x
#  divd(r, x, y)
#  r
#proc `/`*(x:SComplex|DComplex; y:SomeNumber):auto =
#  mixin divd
#  var r{.noInit.}:type(x)
#  #echoType: x
#  divd(r, x, y)
#  r

#template assign*(r: DComplexV, x: SomeNumber) = assignCU(r, x)
#template assign*(r: DColorVectorV, x: SomeNumber) = assignVS(r, x)
#template assign*(r: DColorMatrixV, x: SomeNumber) = assignMS(r, x)
#template `*`*(x: float64, y: DColorVectorV): untyped = mul(x, y)
#template `*`*(x: ToDouble[Dvec0], y: Dvec0): untyped = `*`(x[],y)

proc assign*(r:var SomeNumber; m:Masked[SDvec]) =
  var i = 0
  var b = m.mask
  while b != 0:
    if (b and 1) != 0:
      r = m.pobj[][i]
      break
    b = b shr 1
    i.inc
template `:=`*(r: SomeNumber; m: Masked[SDvec]) = assign(r, m)
proc assign*(m:Masked[SDvec], x:SomeNumber) =
  var i = 0
  var b = m.mask
  while b != 0:
    if (b and 1) != 0:
      m.pobj[][i] = x
    b = b shr 1
    i.inc
proc `:=`*(m:Masked[SDvec], x:SomeNumber) = assign(m, x)
proc assign*(m:Masked[SDvec], x:SDvec) =
  var i = 0
  var b = m.mask
  while b != 0:
    if (b and 1) != 0:
      m.pobj[][i] = x[i]
    b = b shr 1
    i.inc
proc `:=`*(m:Masked[SDvec], x:SDvec) = assign(m,x)
proc assign*(m:Masked[SDvec], x:Masked[SDvec]) =
  ## Only works for the same number of unmasked bits,
  ## and assign those from RHS to LHS in sequence.
  var
    i,j = 0
    b = m.mask
    c = x.mask
  while b != 0:
    if (b and 1) != 0:
      while c != 0:
        let p = (c and 1) != 0
        if p: m.pobj[][i] = x.pobj[][j]
        c = c shr 1
        j.inc
        if p: break
    b = b shr 1
    i.inc
proc `:=`*(m:Masked[SDvec], x:Masked[SDvec]) = assign(m,x)
proc mul*(m:Masked[SDvec]; x:SDvec; y:SomeNumber) =
  var i = 0
  var b = m.mask
  #echo b
  while b != 0:
    if (b and 1) != 0:
      #echo i
      m.pobj[][i] = x[i] * (type(m.pobj[][i]))(y)
    b = b shr 1
    i.inc
#[
proc mul*(m:Masked[SDvec]; x:SomeNumber): auto =
  mul(m, m.pobj[], x)
  m
proc `*`*(m:Masked[SDvec]; x:SomeNumber): auto =
  mul(m, m.pobj[], x)
  m
]#
template `*`*(m: Masked[SDvec]; x: SomeNumber): untyped =
  var t_mulMasked = m
  var t_mulMasked2 = t_mulMasked.pobj[] * x
  maskedObj(t_mulMasked2, t_mulMasked.mask)

proc imul*(m:Masked[SDvec]; x:SomeNumber) =
  var t = m[]
  imul(t, x)
  assign(m, t)
proc `*=`*(m:Masked[SDvec]; x:SomeNumber) = m.imul x
#proc assign*(m:Masked[SComplexV], y:int) =
#  let p = m.pobj[].re.addr
#  assign(Masked[type(p[])](pobj:p,mask:m.mask), y)
#proc mul*(r:Masked[SComplexV], x:SComplexV, y:int) =
#  let prr = r.pobj[].re.addr
#  var mrr = Masked[type(prr[])](pobj:prr,mask:r.mask)
#  mul(mrr, x.re, y)
#  let pri = r.pobj[].im.addr
#  var mri = Masked[type(pri[])](pobj:pri,mask:r.mask)
#  mul(mri, x.im, y)
proc norm2*(m:Masked[SDvec]):auto =
  var r:type(m.pobj[][0])
  var i = 0
  var b = m.mask
  while b != 0:
    if (b and 1) != 0:
      let t = m.pobj[][i]
      r += t*t
    b = b shr 1
    inc i
  r
proc norm2*(r:var SomeNumber; m:Masked[SDvec]) =
  let t = norm2(m)
  r = (type(r))(t)
proc inorm2*(r:var SomeNumber; m:Masked[SDvec]) =
  let t = norm2(m)
  r += (type(r))(t)

template permX*[T](r:T; prm:int; x0:T) =
  block:
    let xp = getPtr x0; template x:untyped = xp[]
    const n = x.numNumbers div x.simdLength
    let rr = cast[ptr array[n,simdType(r)]](addr r)
    let xx = cast[ptr array[n,simdType(x)]](addr x)
    template loop(f:untyped) =
      when compiles(f(rr[0], xx[0])):
        forStatic i, 0, n-1: f(rr[i], xx[i])
    case prm
    of 0: loop(assign)
    of 1: loop(perm1)
    of 2: loop(perm2)
    of 4: loop(perm4)
    of 8: loop(perm8)
    else: discard
proc perm*[T](r0: var T; prm: int; x0: T) {.alwaysInline.} =
  #mixin assign, perm1, perm2, perm4, perm8
  const n = x0.numNumbers div x0.simdLength
  var r = cast[ptr array[n,simdType(r0)]](addr r0)
  let x = cast[ptr array[n,simdType(x0)]](unsafeaddr x0)
  template loop(f:untyped) =
    #mixin assign, perm1, perm2, perm4, perm8
    when compiles(f(r[][0], x[][0])):
      forStatic i, 0, n-1: f(r[][i], x[][i])
  case prm
  #of 0: loop(assign)  # doesn't work
  of 1: loop(perm1)
  of 2: loop(perm2)
  of 4: loop(perm4)
  of 8: loop(perm8)
  else: discard
template perm2*[T](r: var T; prm: int; x: T) =
  const n = x.nVectors
  let rr = cast[ptr array[n,simdType(r)]](r.addr)
  var xt = x
  let xx = cast[ptr array[n,simdType(xt)]](addr xt)
  forStatic i, 0, n-1:
    rr[i] = perm(xx[i], prm)
#proc pack*(r:ptr auto; l:ptr auto; pck:int; x:PackTypes) {.inline.} =
proc pack*(r:ptr auto; l:ptr auto; pck:int; x:auto) {.inline.} =
  mixin simdLength
  if pck==0:
    #const n = x.nVectors
    const n = x.numNumbers div x.simdLength
    let rr = cast[ptr array[n,array[simdLength(x),type(r[])]]](r)
    let xx = cast[ptr array[n,simdType(x)]](unsafeAddr(x))
    for i in 0..<n:
      assign(rr[i], xx[i])
  else:
    #const n = x.nVectors
    const n = x.numNumbers div x.simdLength
    const vl2 = x.simdLength div 2
    let rr = cast[ptr array[n,array[vl2,type(r[])]]](r)
    let ll = cast[ptr array[n,array[vl2,type(l[])]]](l)
    let xx = cast[ptr array[n,simdType(x)]](unsafeAddr(x))
    template loop(f:untyped):untyped =
      forStatic i, 0, n.pred: f(rr[i], xx[i], ll[i])
    case pck
    of  1: loop(packp1)
    of -1: loop(packm1)
    of  2: loop(packp2)
    of -2: loop(packm2)
    of  4: loop(packp4)
    of -4: loop(packm4)
    of  8: loop(packp8)
    of -8: loop(packm8)
    else: discard
#proc pack*(r:ptr char; pck:int; x:PackTypes) =
proc pack*(r:ptr char; pck:int; x: auto) =
  if pck==0:
    const n = x.nVectors
    let rr = cast[ptr array[n,simdType(r)]](r)
    let xx = cast[ptr array[n,simdType(x)]](unsafeAddr(x))
    rr[] = xx[]
  else:
    const n = x.nVectors
    const vl2 = x.simdLength div 2
    let rr = cast[ptr array[n,array[vl2,numberType(x)]]](r)
    let xx = cast[ptr array[n,simdType(x)]](unsafeAddr(x))
    template loop(f:untyped):untyped =
      forStatic i, 0, n.pred: f(rr[i], xx[i])
    case pck
    of  1: loop(packp1)
    of -1: loop(packm1)
    of  2: loop(packp2)
    of -2: loop(packm2)
    of  4: loop(packp4)
    of -4: loop(packm4)
    of  8: loop(packp8)
    of -8: loop(packm8)
    else: discard
proc blend*(r:var auto; x:ptr char; b:ptr char; blnd:int) {.inline.} =
  #const n = r.nVectors
  const n = r.numNumbers div r.simdLength
  #const n2 = n div 2
  const stride = r.simdLength div 2
  var rr = cast[ptr array[n,simdType(r)]](r.addr)
  let xx = cast[ptr array[n,array[stride,numberType(r)]]](x)
  let bb = cast[ptr array[n,array[stride,numberType(r)]]](b)
  template loop(f:untyped):untyped =
    forStatic i, 0, n.pred: f(rr[i], xx[i], bb[i])
  case blnd
  of  1: loop(blendp1)
  of -1: loop(blendm1)
  of  2: loop(blendp2)
  of -2: loop(blendm2)
  of  4: loop(blendp4)
  of -4: loop(blendm4)
  of  8: loop(blendp8)
  of -8: loop(blendm8)
  else: discard

template newDComplexV*[V:static[int]](l: Layout[V]): auto =
  ComplexType[`SimdD V`]()

proc newField*[V:static[int]](l: Layout[V], T:typedesc): Field[V,T] =
  result.new(l)

proc ColorMatrix*(l: Layout, n: static[int]): auto =
  type C = type(l.newDComplexV)
  type CM = ColorMatrixN[n,C]
  result = l.newField(CM)

macro makeConstructors(x: untyped): untyped =
  template mp(f,r,rslt: untyped) =
    proc f*(l: Layout): r =
      # ## create `r`
      new(rslt, l)
  let f = $x
  let r = "Lattice" & f
  result = newStmtList()
  result.add getAst mp(ident(f&"S"), ident("S"&r&"V"), ident"result")
  result.add getAst mp(ident(f&"D"), ident("D"&r&"V"), ident"result")
  result.add getAst mp(ident(f), ident(defPrec&r&"V"), ident"result")
  # non-Simd versions
  result.add getAst mp(ident(f&"S1"), ident("S"&r), ident"result")
  result.add getAst mp(ident(f&"D1"), ident("D"&r), ident"result")
  result.add getAst mp(ident(f&"1"), ident(defPrec&r), ident"result")
  #echo result.repr

makeConstructors(Real)
makeConstructors(Complex)
#makeConstructors(complex)
makeConstructors(ColorVector)
makeConstructors(ColorMatrix)
makeConstructors(HalfFermion)
makeConstructors(DiracFermion)

#proc ColorVectorS*(l: Layout): SLatticeColorVectorV = result.new(l)
#proc ColorVectorD*(l: Layout): DLatticeColorVectorV = result.new(l)
#proc ColorVector*(l: Layout): auto = ColorVectorD(l)
#proc ColorMatrixS*(l: Layout): SLatticeColorMatrixV = result.new(l)
#proc ColorMatrixD*(l: Layout): DLatticeColorMatrixV = result.new(l)
#proc ColorMatrix*(l: Layout): auto = ColorMatrixD(l)
#proc DiracFermionS*(l: Layout): SLatticeDiracFermionV = result.new(l)
#proc DiracFermionD*(l: Layout): DLatticeDiracFermionV = result.new(l)
#proc DiracFermion*(l: Layout): auto = DiracFermionD(l)

when isMainModule:
  import times
  import qex
  qexInit()
  echo "rank ", myRank, "/", nRanks
  #destructors.newSeq(0)
  #var lat = [4,4,2,2]
  #var lat = [4,4,4,4]
  #var lat = [8,8,4,4]
  var lat = [8,8,8,8]
  #var lat = [16,16,8,8]
  #var lat = [16,16,16,16]
  #var lat = [32,32,16,16]
  var lo = newLayout(lat)
  #layout.makeShift(0,1)
  #layout.makeShift(3,-2,"even")
  var v1 = lo.ColorVector()
  var v2 = lo.ColorVector()
  var m1 = lo.ColorMatrix()

  echo threadNum, "/", numThreads
  v1 := 1
  v2 := 0
  m1 := 0
  echo "done init"

  threads:
    m1["odd"] := 1
    v2["odd"] := v1
    threadBarrier()
    echo m1.norm2/lo.nSites.float
    echo v2.norm2/lo.nSites.float
  #var ns:array[8,string]
  #for i in 0..<numThreads: ns[i] = $i
  let nrep = int(1e9/lo.nSites.float)
  #let nrep = 1
  let t0 = epochTime()
  threads:
    for i in 1..nrep:
      #mul(v2[0], m1[0], v1[0])
      #let m10 = m1[0]
      #let v10 = v1[0]
      #let t = m10 * v10
      #assign(v2[0], m1[0] * v1[0])
      #v2[0] := m1[0] * v1[0]
      v2 := m1 * v1
      #mul(m1, v1, v2)
      #for e in r.all:
      #echo ns[threadNum],"/",numThreads
      #for e in allX(v2,threadNum,numThreads):
      #for e in all(v2):
      #  mul(m1.s[e], v1.s[e], v2.s[e])
  let t1 = epochTime()
  echo "time: ", (t1-t0)
  echo v1.s[0][0]
  echo v1.s[lo.nSitesOuter-1][0]
  echo v2.s[lo.nSitesOuter-1][0]
  echo "mflops: ", (66e-6*lo.nSites.float*nrep.float)/(t1-t0)

  echo("mem: (used+free)/total: (", getOccupiedMem(), "+", getFreeMem(), ")/",
       getTotalMem())

  echo GC_getStatistics()
  GC_fullCollect()
  echo GC_getStatistics()
  #echo "destructors: ", destructors.len
  #for f in destructors: f()
  #echo GC_getStatistics()
  #GC_fullCollect()
  #echo GC_getStatistics()
  v1 = nil
  v2 = nil
  m1 = nil
  echo GC_getStatistics()
  GC_fullCollect()
  echo GC_getStatistics()
  echo("mem: (used+free)/total: (", getOccupiedMem(), "+", getFreeMem(), ")/",
       getTotalMem())

  qexFinalize()
