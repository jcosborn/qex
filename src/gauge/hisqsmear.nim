import qex
import physics/[hisqLinks]
import gauge
import gauge/[fat7l,fat7lderiv]

export hisqLinks

proc newHISQ*(lepage: float = 0.0; naik: float = 1.0): HisqCoefs =
  result = HisqCoefs(naik: -naik/24.0)
  result.fat7first.setHisqFat7(lepage,0.0)
  result.fat7second.setHisqFat7(2.0-lepage,naik)

proc smearGetForce*[T](
    self: HisqCoefs; 
    u: T; 
    su,sul: T;
    displayPerformance: bool = false
  ): proc(dsdu: var T; dsdsu,dsdsul: T) =
  mixin projectU,projectUderiv
  let
    lo = u[0].l
    fat7l1 = self.fat7first
    fat7l2 = self.fat7second
    naik = self.naik
  var
    v = newOneOf(u)
    w = newOneOf(u)
    info: PerfInfo
  
  # Smear
  v.makeImpLinks(u,fat7l1,info) # First fat7
  threads: # Unitary projection
    for mu in 0..<w.len: 
      for s in w[mu]: w[mu][s].projectU(v[mu][s])
  makeImpLinks(su,w,fat7l2,sul,w,naik,info) # Second fat7

  # Chain rule - retains a reference to u,su,sul
  proc smearedForce(dsdu: var T; dsdsu,dsdsul: T) =
    var 
      dsdx_dxdw = newOneOf(dsdu)
      dsdx_dxdw_dwdv = newOneOf(dsdu)
    dsdx_dxdw.fat7lderiv(dsdsu,su,fat7l2,dsdsul,sul,naik,info) # Second fat7
    threads: # Unitary projection
      for mu in 0..<dsdx_dxdw_dwdv.len:
        for s in dsdx_dxdw_dwdv[mu]:
          dsdx_dxdw_dwdv[mu][s].projectUderiv(w[mu][s],v[mu][s],dsdx_dxdw[mu][s]) 
    dsdu.fat7lderiv(dsdx_dxdw_dwdv,u,fat7l1,info) # First fat7
  
  if displayPerformance: echo $(info)
  return smearedForce

if isMainModule:
  qexInit()
  let 
    defaultLat = @[8,8,8,8]
    hisq = newHISQ()
  defaultSetup()
  var
    sg = lo.newGauge()
    sgl = lo.newGauge()
    f = lo.newGauge()
    ff = lo.newGauge()
    ffl = lo.newGauge()
  g.random
  echo f.plaq
  echo sg.plaq
  echo sgl.plaq
  hisq.smear(g,sg,sgl)
  echo "--"
  echo sg.plaq
  echo sgl.plaq
  ff.gaugeForce(sg)
  ffl.gaugeForce(sgl)
  var force = hisq.smearGetForce(g,sg,sgl)
  f.force(ff,ffl)
  echo "--"
  echo f.plaq
  echo sg.plaq
  echo sgl.plaq
  qexFinalize()