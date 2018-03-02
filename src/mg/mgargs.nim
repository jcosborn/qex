import qex
import mgblocks

type
  MgRestrictor*[VF,VC: static[int]; R] = object
    mgb: MgBlock[VF,VC]
    rv: Field[VF,R]
  MgProlongator*[VF,VC: static[int]; P] = object
    mgb: MgBlock[VF,VC]
    pv: Field[VF,P]

proc newMgRestrictor*[VF,VC: static[int]; R](mgb: MgBlock[VF,VC],
                                  rv: Field[VF,R]): MgRestrictor[VF,VC,R] =
  result.mgb = mgb
  result.rv = rv

proc newMgProlongator*[VF,VC: static[int]; P](mgb: MgBlock[VF,VC],
                                  pv: Field[VF,P]): MgProlongator[VF,VC,P] =
  result.mgb = mgb
  result.pv = pv

template paritySitesV*(f: Field, par: int): untyped =
  var s: Subset
  s.paritySubset(f.l, par)
  f[s]

proc apply*[VF,VC:static[int];R](r: MgRestrictor[VF,VC,R],
                                 fc: var Field, ff: Field2, par=0) =
  let nv = r.rv[0].len
  #var s = 0.0
  #var s: type(dot(r.rv[0][0],ff[0]))
  for fsv in ff.paritySitesV(par):
    for i in 0..<nv:
      let t = dot(r.rv[fsv][i],ff[fsv])
      #s += t
      for j in 0..<VF:
        let cs = r.mgb.csites[fsv*VF+j]
        #if i==0: echo "fsv: ", fsv, "  j: ", j, "  cs: ", cs
        let csv = cs div VC
        let cj = cs mod VC
        fc[csv][i].re[cj] = fc[csv][i].re[cj] + t.re[j]
        fc[csv][i].im[cj] = fc[csv][i].im[cj] + t.im[j]
        #let tt = fc[csv][i].re[cj] + t.re[j]
        #fc[csv][i].re[cj] = tt
        #if i==0: s += t*t
  #echo "s: ", s

proc apply*[VF,VC:static[int];P](p: MgProlongator[VF,VC,P],
                                 ff: var Field, fc: Field2, par=0) =
  let nv = p.pv[0].len
  for fsv in ff.paritySitesV(par):
    for i in 0..<nv:
      var t: type(dot(p.pv[fsv][i],ff[fsv]))
      for j in 0..<VF:
        let cs = p.mgb.csites[fsv*VF+j]
        let csv = cs div VC
        let cj = cs mod VC
        #let tt = fc[csv][i].re[cj]
        #t.re[j] = tt
        t.re[j] = fc[csv][i].re[cj]
        t.im[j] = fc[csv][i].im[cj]
      ff[fsv] += t * p.pv[fsv][i]

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

  const nmgv1 = 2
  type
    SMgR1V* = array[nmgv1,SDiracFermionV]
    SMgP1V* = array[nmgv1,SDiracFermionV]
    SLatticeMgR1V* = Field[loF.V,SMgR1V]
    SLatticeMgP1V* = Field[loF.V,SMgP1V]

    SMgColorVector1V* = Color[VectorArray[nmgv1,SComplexV]]
    SLatticeMgVector1V* = Field[loC.V,SMgColorVector1V]

  template assign(r: var SMgR1V, x: SMgR1V) =
    for i in 0..<r.len:
      r[i] := x[i]
  proc mgextract(r: var Field, x: Field2, n: int) =
    for e in r:
      r[e] := x[e][n]
  proc mginsert(r: var Field, x: Field2, n: int) =
    for e in r:
      r[e][n] := x[e]
  proc mgdot(x: Field, y: Field2): array[nmgv1,float] =
    var s: array[nmgv1,type(dot(x[0][0],y[0]))]
    for e in x:
      for i in 0..<nmgv1:
        s[i] += dot(x[e][i], y[e])
  proc normalize(x: var Field) =
    let t = x.norm2
    let s = 1/sqrt(t)
    x := s*x
  proc invsqrt(x: var Field, i: int) =
    for e in x:
      let t = x[e][i].re
      x[e] := 0
      for j in 0..<t.numNumbers:
        if t[j] > 0.0:
          let s = 1/sqrt(t[j])
          x[e][i].re[j] = s

  var rv: SLatticeMgR1V
  var pv: SLatticeMgP1V
  rv.new(loF)
  pv.new(loF)
  var fv = loF.DiracFermionS()
  var fv2 = loF.DiracFermionS()
  var rs = newRNGField(RngMilc6, loF, 987654321)
  let r = newMgRestrictor(b, rv)
  let p = newMgProlongator(b, pv)
  var cv: SLatticeMgVector1V
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
    r.apply(cv, fv)
    fv2 := 0
    p.apply(fv2, cv)
    echo "fv2: ", fv2.even.norm2
    fv -= fv2
    echo "fv:  ", fv.even.norm2
    rv.mginsert(fv, i)
    pv := rv

    cv := 0
    r.apply(cv, fv)
    cv.invsqrt(i)
    fv2 := 0
    p.apply(fv2, cv)
    echo "fv2: ", fv2.even.norm2
    rv.mginsert(fv2, i)
    pv := rv

  for i in 0..<nmgv1:
    mgextract(fv, rv, i)
    echo "fv even: ", fv.even.norm2
    cv := 0
    r.apply(cv, fv)
    #echo "cv all:  ", cv.norm2
    #echo "cv even: ", cv.even.norm2
    #echo "cv odd:  ", cv.odd.norm2
    fv2 := 0
    p.apply(fv2, cv)
    echo "fv2: ", fv2.norm2

  qexFinalize()

