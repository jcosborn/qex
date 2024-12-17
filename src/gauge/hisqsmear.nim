import qex
import gauge
import gauge/[fat7l,fat7lderiv]

export hisqLinks

proc fat7lDeriv(
    mid: auto, 
    gauge: auto, 
    deriv: auto, 
    coef: Fat7lCoefs,
    llgauge: auto, 
    llderiv: auto, 
    naik: float,
    perf: var PerfInfo
  ) =
  var (f,fll) = (newOneOf(deriv),newOneOf(llderiv))
  fat7lderiv(f,gauge,deriv,coef,fll,llgauge,llderiv,naik,perf)
  threads:
    for mu in 0..<mid.len:
      for s in mid[mu]: mid[mu][s] := f[mu][s] + fll[mu][s]

proc fat7lDeriv(mid: auto; gauge,deriv: auto; coef: Fat7lCoefs; perf: var PerfInfo) =
  fat7lderiv(mid,gauge,deriv,coef,mid,gauge,deriv,0.0,perf)

proc projectU[T](v: auto; u: T) =
  threads:
    for mu in 0..<u.len: 
      for s in u[mu]: v[mu][s].projectU(u[mu][s])

proc projectUderiv[T](dvdu: auto; v,u: T; chain: T) =
  threads:
    for mu in 0..<chain.len:
      for s in chain[mu]: 
        dvdu[mu][s].projectUderiv(v[mu][s],u[mu][s],chain[mu][s])

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
  w.projectU(v) # Unitary projection
  makeImpLinks(su,w,fat7l2,sul,w,naik,info) # Second fat7

  # Chain rule - retains a reference to u,su,sul
  proc smearedForce(dsdu: var T; dsdsu,dsdsul: T) =
    var (t1,t2) = (newOneOf(dsdu),newOneOf(dsdu))
    t1.fat7lderiv(w,dsdsu,fat7l2,w,dsdsul,naik,info) # Second fat7
    t2.projectUderiv(w,v,t1) # Unitary projection
    dsdu.fat7lderiv(u,t2,fat7l1,info) # First fat7
    if displayPerformance: echo $(info)
  
  # Display performance (if requested) and return
  if displayPerformance: echo $(info)
  return smearedForce