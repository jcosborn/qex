import alignedMem
export alignedMem
import layout
import basicOps
import simd
export simd
import macros
import metaUtils
import threading
import strutils
import globals
import comms
import field
import complexConcept
export complexConcept
import matrixConcept
export matrixConcept

#var destructors:seq[proc()]

const nc = 3
setType(Svec0, "SimdS" & $VLEN)
setType(Dvec0, "SimdD" & $VLEN)
type
  SDvec = Svec0 | Dvec0
  SComplex* = AsComplex[tuple[re,im:float32]]
  SComplexV* = AsComplex[tuple[re,im:Svec0]]
  SColorVector* = VectorArray[nc,SComplex]
  SColorVectorV* = VectorArray[nc,SComplexV]
  SColorMatrix* = MatrixArray[nc,nc,SComplex]
  SColorMatrixV* = MatrixArray[nc,nc,SComplexV]
  SLatticeReal* = Field[1,float32]
  SLatticeRealV* = Field[VLEN,Svec0]
  SLatticeComplex* = Field[1,SComplex]
  SLatticeComplexV* = Field[VLEN,SComplexV]
  SLatticeColorVector* = Field[1,SColorVector]
  SLatticeColorVectorV* = Field[VLEN,SColorVectorV]
  SLatticeColorMatrix* = Field[1,SColorMatrix]
  SLatticeColorMatrixV* = Field[VLEN,SColorMatrixV]
  DComplex* = AsComplex[tuple[re,im:float64]]
  DComplexV* = AsComplex[tuple[re,im:Dvec0]]
  DColorVector* = VectorArray[nc,DComplex]
  DColorVectorV* = VectorArray[nc,DComplexV]
  DColorMatrix* = MatrixArray[nc,nc,DComplex]
  DColorMatrixV* = MatrixArray[nc,nc,DComplexV]
  DLatticeReal* = Field[1,float64]
  DLatticeRealV* = Field[VLEN,Dvec0]
  DLatticeComplex* = Field[1,DComplex]
  DLatticeComplexV* = Field[VLEN,DComplexV]
  DLatticeColorVector* = Field[1,DColorVector]
  DLatticeColorVectorV* = Field[VLEN,DColorVectorV]
  DLatticeColorMatrix* = Field[1,DColorMatrix]
  DLatticeColorMatrixV* = Field[VLEN,DColorMatrixV]

template numberType*[T1,T2](x:tuple[re:T1,im:T2]):untyped = numberType(T1)
template numberType*[I,T](x:array[I,T]):untyped = numberType(T)
template numberType*[T](x:AsComplex[T]):untyped = numberType(T)
template numberType*[T](x:AsVector[T]):untyped = numberType(T)
template numberType*[T](x:AsMatrix[T]):untyped = numberType(T)
#template numberType*[N,T](x:typedesc[array[N,T]]):untyped = numberType(T)
template numberType*[T](x:typedesc[AsVector[T]]):untyped = numberType(T)
template numberType*[T](x:typedesc[AsMatrix[T]]):untyped = numberType(T)

template simdLength*(x:typedesc[SColorMatrixV]):expr = simdLength(Svec0)
template simdLength*(x:typedesc[SColorVectorV]):expr = simdLength(Svec0)
template simdLength*(x:SColorVectorV):expr = simdLength(Svec0)
template simdLength*(x:SColorMatrixV):expr = simdLength(Svec0)
template simdLength*(x:typedesc[DColorMatrixV]):expr = simdLength(Dvec0)
template simdLength*(x:typedesc[DColorVectorV]):expr = simdLength(Dvec0)
template simdLength*(x:DColorVectorV):expr = simdLength(Dvec0)
template simdLength*(x:DColorMatrixV):expr = simdLength(Dvec0)
template simdLength*(x:AsComplex):expr = simdLength(x.re)
template simdLength*(x:AsMatrix):expr = simdLength(x[0,0])

#template nVectors(x:SColorVectorV):expr = 2*nc
#template nVectors(x:SColorMatrixV):expr = 2*nc*nc
#template nVectors(x:DColorVectorV):expr = 2*nc
#template nVectors(x:DColorMatrixV):expr = 2*nc*nc
template nVectors(x:Svec0):expr = 1
template nVectors(x:Dvec0):expr = 1
template nVectors(x:AsComplex):expr = 2*nVectors(x.re)
template nVectors(x:AsVector):expr = x.len*nVectors(x[0])
template nVectors(x:AsMatrix):expr = x.nrows*x.ncols*nVectors(x[0,0])

template simdType*(x:tuple):expr = simdType(x[0])
template simdType*(x:array):expr = simdType(x[x.low])
template simdType*(x:AsComplex):expr = simdType(x[])
template simdType*(x:AsVector):expr = simdType(x[])
template simdType*(x:AsMatrix):expr = simdType(x[])


#import complexConcept
#export complexConcept
#declareComplex(SComplex)
#declareComplex(SComplexV)


#proc `$`*(x:Masked[SComplexV]):string =
#  result = "(" & $x.re & "," & $x.im & ")"

template trace*(x:SComplexV):expr = x
proc simdSum*(x:SComplexV):SComplex = complexConcept.map(result, simdSum, x)
proc simdSum*(x:DComplexV):DComplex = complexConcept.map(result, simdSum, x)
#template rankSum*(x:AsComplex) =
#  #mixin qmpSum
#  rankSum(x[])
#template `/`*(x:SComplex|DComplex; y:SomeNumber):expr =
#  mixin divd
#  var r{.noInit.}:type(x)
#  echoType: x
#  divd(r, x, y)
#  r
proc `/`*(x:SComplex|DComplex; y:SomeNumber):auto =
  mixin divd
  var r{.noInit.}:type(x)
  #echoType: x
  divd(r, x, y)
  r

proc assign*(r:var SomeNumber; m:Masked[SDvec]) =
  var i = 0
  var b = m.mask
  while b != 0:
    if (b and 1) != 0:
      r = m.pobj[][i]
      break
    b = b shr 1
    i.inc
proc assign*(m:Masked[SDvec], x:SomeNumber) =
  var i = 0
  var b = m.mask
  while b != 0:
    if (b and 1) != 0:
      m.pobj[][i] = x
    b = b shr 1
    i.inc
proc assign*(m:Masked[SDvec], x:SDvec) =
  var i = 0
  var b = m.mask
  while b != 0:
    if (b and 1) != 0:
      m.pobj[][i] = x[i]
    b = b shr 1
    i.inc
proc mul*(m:Masked[SDvec]; x:SDvec; y:int) =
  var i = 0
  var b = m.mask
  #echo b
  while b != 0:
    if (b and 1) != 0:
      #echo i
      m.pobj[][i] = x[i] * (type(m.pobj[][i]))(y)
    b = b shr 1
    i.inc
proc imul*(m:Masked[SDvec]; x:int) =
  var t = m[]
  imul(t, x)
  assign(m, t)
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

#template isScalar*(x:float32):expr = true
#template isScalar*(x:Scalar):expr = true
#template isScalar*(x:SComplexV):expr = true
#template isVector*(x:SColorVectorV):expr = true
#template mvLevel*(x:SColorVectorV):expr = 1
#template isMatrix*(x:SColorMatrixV):expr = true
#template mvLevel*(x:SColorMatrixV):expr = 1
#template nrows*(x:SColorMatrixV):expr = nc
#template ncols*(x:SColorMatrixV):expr = nc
#template `[]`*(x:SColorMatrixV; i,j:int):untyped = x[i][j]
#template `[]=`*(x:SColorMatrixV; i,j:int, y:untyped):untyped = x[i][j] = y
#template `[]`*(x:SColorMatrix; i,j:int):untyped = x[i][j]
#template `[]=`*(x:SColorMatrix; i,j:int, y:untyped):untyped = x[i][j] = y

#import matrixConcept
#export matrixConcept

#template assign*(r:var SColorMatrixV, x:SomeNumber):untyped = assign(r, x.toScalar)

proc prefetch*(x:ptr AsComplex) {.inline.} =
  prefetch(addr(x[].re))
  prefetch(addr(x[].im))
proc prefetch*(x:ptr AsVector) {.inline.} =
  for i in 0..<x[].len:
    prefetch(addr(x[][i]))

#type PackTypes = SColorVectorV | SColorMatrixV | DColorVectorV | DColorMatrixV
type PackTypes = any
proc perm*[T](r:var T; prm:int; x:T) {.inline.} =
  const n = x.nVectors
  let rr = cast[ptr array[n,simdType(r)]](r.addr)
  let xx = cast[ptr array[n,simdType(x)]](unsafeAddr(x))
  template loop(f:untyped):untyped =
    when compiles(f(rr[0], xx[0])):
      forStatic i, 0, n-1: f(rr[i], xx[i])
  case prm
  of 0: loop(assign)
  of 1: loop(perm1)
  of 2: loop(perm2)
  of 4: loop(perm4)
  of 8: loop(perm8)
  else: discard
proc pack*(r:ptr any; l:ptr any; pck:int; x:PackTypes) {.inline.} =
  if pck==0:
    const n = x.nVectors
    let rr = cast[ptr array[n,array[simdLength(x),type(r[])]]](r)
    let xx = cast[ptr array[n,simdType(x)]](unsafeAddr(x))
    for i in 0..<n:
      assign(rr[i], xx[i])
  else:
    const n = x.nVectors
    const vl2 = x.simdLength div 2
    let rr = cast[ptr array[n,array[vl2,type(r[])]]](r)
    let ll = cast[ptr array[n,array[vl2,type(l[])]]](l)
    let xx = cast[ptr array[n,simdType(x)]](unsafeAddr(x))
    template loop(f:untyped):untyped =
      forStatic i, 0, <n: f(rr[i], xx[i], ll[i])
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
proc pack*(r:ptr char; pck:int; x:PackTypes) =
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
      forStatic i, 0, <n: f(rr[i], xx[i])
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
proc blend*(r:var any; x:ptr char; b:ptr char; blnd:int) {.inline.} =
  const n = r.nVectors
  const n2 = n div 2
  const stride = r.simdLength div 2
  var rr = cast[ptr array[n,simdType(r)]](r.addr)
  let xx = cast[ptr array[n,array[stride,numberType(r)]]](x)
  let bb = cast[ptr array[n,array[stride,numberType(r)]]](b)
  template loop(f:untyped):untyped =
    forStatic i, 0, <n: f(rr[i], xx[i], bb[i])
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


proc ColorVector*(l:Layout):SLatticeColorVectorV = result.new(l)
proc ColorMatrix*(l:Layout):SLatticeColorMatrixV = result.new(l)
#proc ColorVector*(l:Layout):DLatticeColorVectorV = result.new(l)
#proc ColorMatrix*(l:Layout):DLatticeColorMatrixV = result.new(l)

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
