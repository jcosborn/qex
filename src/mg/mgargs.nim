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
  for fsv in ff.paritySitesV(par):
    for i in 0..<nv:
      let t = dot(r.rv[fsv][i],ff[fsv])
      for j in 0..<VF:
        let cs = r.mgb.csites[fsv*VF+j]
        let csv = cs div VC
        let cj = cs mod VC
        fc[csv][i].re[cj] = fc[csv][i].re[cj] + t.re[j]
        fc[csv][i].im[cj] = fc[csv][i].im[cj] + t.im[j]

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
        t.re[j] = fc[csv][i].re[cj]
        t.im[j] = fc[csv][i].im[cj]
      ff[fsv] += t * p.pv[fsv][i]

when isMainModule:
  qexInit()
  let latF = [16,16,16,16]
  let loF = newLayout(latF)

  let latC = [8,8,8,8]
  let loC = newLayout(latC, loF.V, loF.rankGeom, loF.innerGeom)

  let b = newMgBlock(loF, loC)

  const nmgv1 = 20
  type
    SMgR1V* = array[nmgv1,SDiracFermionV]
    SMgP1V* = array[nmgv1,SDiracFermionV]
    SLatticeMgR1V* = Field[loF.V,SMgR1V]
    SLatticeMgP1V* = Field[loF.V,SMgP1V]

    SMgColorVector1V* = Color[VectorArray[nmgv1,SComplexV]]
    SLatticeMgVector1V* = Field[loC.V,SMgColorVector1V]

  var rv: SLatticeMgR1V
  var pv: SLatticeMgP1V
  rv.new(loF)
  pv.new(loF)

  var rs = newRNGField(RngMilc6, lo, 987654321)
  gaussian(rv, rs)

  let r = newMgRestrictor(b, rv)
  let p = newMgProlongator(b, pv)
  var fv = loF.DiracFermionS()
  var cv: SLatticeMgVector1V
  cv.new(loC)


  r.apply(cv, fv)
  p.apply(fv, cv)

  qexFinalize()



#[

]#
