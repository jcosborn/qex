import qex
import mgblocks, mgargs

type
  MgColorVectorV*[N:static[int]] = Color[VectorArray[N,SComplexV]]
  LatticeMgColorVectorV*[N:static[int]] = Field[VLEN,MgColorVectorV[N]]

  WmgVectorsV*[N:static[int]] = array[N,SDiracFermionV]
  LatticeWmgVectorsV*[N:static[int]] = Field[VLEN,WmgVectorsV[N]]

template assign*(r: var WmgVectorsV, x: WmgVectorsV) =
  for i in 0..<r.len:
    r[i] := x[i]

proc wmgZero*(r: Field) =
  let n1 = r[0].len - 1
  for e in r:
    for i in 0..n1:
      r[e][i] := 0

proc wmgZero*(r: Field, n: int) =
  for e in r:
    r[e][n] := 0

proc wmgExtract*(r: Field, x: Field2, n: int) =
  for e in r:
    r[e] := x[e][n]

proc wmgInsert*(r: Field, x: Field2, n: int) =
  for e in r:
    r[e][n] := x[e]

#proc simdSum*[TR,TI](x: ComplexObj[TR,TI]): auto =
#  newComplexObj(simdSum(x.re),simdSum(x.im))

#proc simdSum*[T](x: ComplexProxy[T]): auto =
#  newComplexProxy(simdSum(x))

proc simdSum*[N,T](x: array[N,T]): auto =
  var r: array[N, type(simdSum(x[0]))]
  for i in x.low .. x.high:
    r[i] = simdSum x[i]
  r

#proc mgdot(x: Field, y: Field2): array[x[0].N,float] =
proc wmgdot*(x: Field, y: Field2): auto =
  const n = x[0].N
  var s: array[n,type(dot(x[0][0],y[0]))]
  for e in x:
    for i in 0..<n:
      s[i] += dot(x[e][i], y[e])
  var t = simdSum s
  #FIXME threadRankSum t
  t

# assumes vectors in t are normalized
proc wmgProject*(x: var Field, t: MgTransfer) =
  var s = wmgdot(t.v, x)
  let n = s.len
  #let fv = t.mgb.fine.physVol
  #let cv = t.mgb.coarse.physVol
  #let f = cv.float/fv.float
  #let f = 1.0/cv.float
  #for i in 0..<n:
  #  s[i] *= f
  #echo s
  for e in x:
    for i in 0..<n:
      x[e] -= s[i] * t.v[e][i]

proc wmgNormalize*(x: var Field) =
  let t = x.norm2
  let s = 1/sqrt(t)
  x := s*x

proc wmgInvsqrt*(x: Field, i: int) =
  for e in x:
    let t = x[e][i].re
    x[e] := 0
    for j in 0..<t.numNumbers:
      if t[j] > 0.0:
        let s = 1/sqrt(t[j])
        x[e][i].re[j] = s

proc wmgBlockProject*(x: var Field, t: MgTransfer, f: Field2, c: Field3) =
  t.restrict(c, x)
  echo "c2: ", c.norm2
  t.prolong(f, c)
  x -= f

  t.restrict(c, x)
  echo "c2: ", c.norm2

proc wmgBlockNormalizeInsert*(t: var MgTransfer, x: var Field, i: int,
                              f: Field2, c: Field3) =
  t.restrict(c, x)
  echo "c2: ", c.norm2
  t.prolong(f, c)
  x -= f

  t.restrict(c, x)
  echo "c2: ", c.norm2

  var v = t.v
  #v.wmgzero()
  v.wmginsert(x, i)
  t.restrict(c, x)
  c.wmginvsqrt(i)
  t.prolong(f, c)
  t.v.wmginsert(f, i)
  echo "f2: ", f.norm2
  x := f

  let x2 = x.norm2
  t.restrict(c, x)
  t.prolong(x, c)
  echo "x2: ", x2, "   ", x.norm2
  x := f

proc wmgBlockProjectInsert*(t: var MgTransfer, x: var Field, i: int,
                            f: Field2, c: Field3) =
  t.restrict(c, x)
  echo "c2: ", c.norm2
  t.prolong(f, c)
  x -= f

  t.restrict(c, x)
  echo "c2: ", c.norm2

  var v = t.v
  #v.wmgzero()
  v.wmginsert(x, i)
  t.restrict(c, x)
  c.wmginvsqrt(i)
  t.prolong(f, c)
  t.v.wmginsert(f, i)
  echo "f2: ", f.norm2
  x := f

  let x2 = x.norm2
  t.restrict(c, x)
  t.prolong(x, c)
  echo "x2: ", x2, "   ", x.norm2
  x := f

when isMainModule:
  qexInit()
  #let latF = [16,16,16,16]
  let latF = [8,8,8,8]
  let loF = newLayout(latF)

  #let latC = [8,8,8,8]
  let latC = [4,4,4,4]
  echo loF.rankGeom
  echo loF.innerGeom
  let loC = newLayout(latC, loF.V, loF.rankGeom, loF.innerGeom)

  let b = newMgBlock(loF, loC)

  const nmgv1 {.intDefine.} = 2
  var rv,pv: LatticeWmgVectorsV[nmgv1]
  rv.new(loF)
  pv.new(loF)
  var fv = loF.DiracFermionS()
  var fv2 = loF.DiracFermionS()
  var rs = newRNGField(RngMilc6, loF, 987654321)
  let r = newMgTransfer(b, rv)
  let p = newMgTransfer(b, pv)
  var cv: LatticeMgColorVectorV[nmgv1]
  cv.new(loC)

  fv := 0
  for i in 0..<nmgv1:
    rv.mginsert(fv, i)
  pv := rv

  #fv.normalize
  #echo "fv all:  ", fv.norm2
  #echo "fv even: ", fv.even.norm2
  #echo "fv odd:  ", fv.odd.norm2

  for i in 0..<nmgv1:
    gaussian(fv, rs)
    echo "gaussian: ", fv.even.norm2
    cv := 0
    r.restrict(cv, fv)
    fv2 := 0
    p.prolong(fv2, cv)
    echo "fv2: ", fv2.even.norm2
    fv -= fv2
    echo "fv:  ", fv.even.norm2
    rv.mginsert(fv, i)
    pv := rv

    cv := 0
    r.restrict(cv, fv)
    cv.invsqrt(i)
    fv2 := 0
    p.prolong(fv2, cv)
    echo "fv2: ", fv2.even.norm2
    rv.mginsert(fv2, i)
    pv := rv

  for i in 0..<nmgv1:
    mgextract(fv, rv, i)
    echo "fv even: ", fv.even.norm2
    cv := 0
    r.restrict(cv, fv)
    #echo "cv all:  ", cv.norm2
    #echo "cv even: ", cv.even.norm2
    #echo "cv odd:  ", cv.odd.norm2
    fv2 := 0
    p.prolong(fv2, cv)
    echo "fv2: ", fv2.norm2

  qexFinalize()

