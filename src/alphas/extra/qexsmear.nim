import base
import layout
import field
import gauge
import gauge/[fat7l, fat7lderiv]
import alphaslinks, alphasproject

export alphaslinks
export alphasproject

proc asqtadDeriv[T](
    deriv: auto, 
    gauge,mid: T,  
    coef: Fat7lCoefs,
    llgauge,llmid: T, 
    naik: float,
    perf: var PerfInfo
  ) =
  var (f,fll) = (newOneOf(mid), newOneOf(mid))
  fat7lderiv(f,gauge,mid,coef,fll,llgauge,llmid,naik,perf)
  threads:
    for mu in 0..<mid.len:
      for s in deriv[mu]: deriv[mu][s] := f[mu][s] + fll[mu][s]

proc fat7Deriv[T](
    deriv: auto,
    gauge,mid: T,
    coef: Fat7lCoefs,
    perf: var PerfInfo
  ) = deriv.fat7lDeriv(gauge,mid,coef,perf)

proc project[T](self: UnitaryProjection, v: auto; u: T; regulate: bool) =
  threads:
    for mu in 0..<u.len:
      for s in u[mu]: self.projectU(v[mu][s], u[mu][s], regulate)

proc projectDeriv[T](
  self: UnitaryProjection, 
  dvdu: auto; 
  v, u: T; 
  chain: T,
  regulate: bool
) =
  threads:
    for mu in 0..<chain.len:
      for s in chain[mu]:
        self.projectUderiv(dvdu[mu][s], v[mu][s], u[mu][s], chain[mu][s], regulate)

template rephase(g: auto) =
  threadBarrier()
  g.setBC
  threadBarrier()
  g.stagPhase
  threadBarrier()

proc smearGetForce*[T](
    self: HisqCoefs; 
    u: T; 
    su,sul: T;
    displayPerformance: bool = false,
    regulate: bool = false
  ): proc(dsdu: var T; dsdsu,dsdsul: T) =
  mixin project, projectDeriv
  let
    lo = u[0].l
    fat7l1 = self.fat7first
    fat7l2 = self.fat7second
    naik = self.naik
  var
    g = newOneOf(u)
    v = newOneOf(u)
    w = newOneOf(u)
    info: PerfInfo
  
  # retain copy of u link, but rephased
  threads: 
    g := u
    threadBarrier()
    g.rephase

  # Smear
  v.makeImpLinks(g,fat7l1,info) # First fat7
  self.projection.project(w, v, regulate = regulate) # Unitary projection
  makeImpLinks(su,w,fat7l2,sul,w,naik,info) # Second fat7

  # Chain rule - retains a reference to u,su,sul
  proc smearedForce(dsdu: var T; dsdsu,dsdsul: T) =
    var t = newOneOf(dsdu)
    t.asqtadDeriv(w,dsdsu,fat7l2,w,dsdsul,naik,info) # Second fat7
    for mu in 0..<t.len:
      echo mu, " ", simdSum(trace(t[mu]))
    self.projection.projectDeriv(t, w, v, t, regulate = regulate)
    dsdu.fat7Deriv(g,t,fat7l1,info) # First fat7
    if displayPerformance: 
      echo &"forceSmear: {info.secs:.5f}s {1e-9*info.flops/info.secs:.3f}Gf/s"
  
  # Display performance (if requested) and return
  if displayPerformance: 
    echo &"linkSmear: {info.secs:.5f}s {1e-9*info.flops/info.secs:.3f}Gf/s"
  info.clear
  return smearedForce

when isMainModule:
  import gauge/gaugeAction, strformat
  qexInit()
  let
    defaultLat = @[8,8,8,8]
  var
    (lo, g, r) = setupLattice(defaultLat)
    dg = lo.newGauge()
    g2 = lo.newGauge()
    fd = lo.newGauge()
    fl = lo.newGauge()
    ll = lo.newGauge()
    fl2 = lo.newGauge()
    ll2 = lo.newGauge()
    fc = lo.newGauge()
    lc = lo.newGauge()
    ch = lo.newGauge()
    chl = lo.newGauge()
    gc = GaugeActionCoeffs(plaq:1.0)
    eps = floatParam("eps", 1e-6)
    hisq = newHISQ()

  g.random r
  g.stagPhase
  dg.gaussian r
  ch.gaussian r
  chl.gaussian r
  for mu in 0..<g.len:
    for e in g[mu]:
      dg[mu][e] *= eps
      g2[mu][e] := g[mu][e] + dg[mu][e]

  proc check(da: float, tol: float) =
    var ds = 0.0
    for mu in 0..3:
      ds += redot(dg[mu], fd[mu])
    let r = (da-ds)/da
    echo &"  da {da}  ds {ds}  rel {r}"
    if abs(r) > tol*eps:
      echo &"> ERROR rel error |{r}| > {tol*eps}"

  proc checkG =
    echo "Checking redot deriv"
    var a, a2 = 0.0
    for mu in 0..<fd.len:
      a += redot(ch[mu], g[mu])
      a2 += redot(ch[mu], g2[mu])
      fd[mu] := ch[mu]
    check(a2-a, 1)
  checkG()

  proc checkL(name: string, tol: float) =
    echo "Checking ", name
    resetTimers()
    let f = hisq.smearGetForce(g, fl, ll)
    let f2 = hisq.smearGetForce(g2, fl2, ll2)
    var a, a2 = 0.0
    for mu in 0..<fd.len:
      a += redot(ch[mu], fl[mu]) + redot(chl[mu], ll[mu])
      a2 += redot(ch[mu], fl2[mu]) + redot(chl[mu], ll2[mu])
      fd[mu] := 0
      fc[mu] := ch[mu]
      lc[mu] := chl[mu]
    f(fd, fc, lc)
    f2(fd, fc, lc)  # combine forces for better accuracy
    check(2.0*(a2-a), tol)  # 2.0 due to combined forces

  hisq.fat7first.oneLink = 1.0
  hisq.fat7first.threeStaple = 0.0
  hisq.fat7first.fiveStaple = 0.0
  hisq.fat7first.sevenStaple = 0.0
  hisq.fat7first.lepage = 0.0
  hisq.fat7second.oneLink = 1.0
  hisq.fat7second.threeStaple = 0.0
  hisq.fat7second.fiveStaple = 0.0
  hisq.fat7second.sevenStaple = 0.0
  hisq.fat7second.lepage = 0.0
  hisq.naik = 0.0
  checkL("oneLink", 1)

  hisq.fat7first.setHisqFat7(0.0, 0.0)
  hisq.fat7second.setHisqFat7(0.0, 0.0)
  hisq.naik = 0.0
  checkL("all", 1)

  hisq = newHisq()
  checkL("all", 1)

  echoProf()
  qexFinalize()