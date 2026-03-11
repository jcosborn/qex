import base, layout, gauge, fat7l, smearutil
export PerfInfo

const keepProj {.boolDefine.} = true
when keepProj:
  static: echo "hypsmear: keeping projected fields"
else:
  static: echo "hypsmear: NOT keeping projected fields"

type HypCoefs* = object
  alpha1*: float
  alpha2*: float
  alpha3*: float

type HypScaleFactors = object
  alp1, alp2, alp3: float
  ma1, ma2, ma3: float

proc `$`*(c: HypCoefs): string =
  result = "Hyp{\n"
  result &= "  alpha1: " & $c.alpha1 & "\n"
  result &= "  alpha2: " & $c.alpha2 & "\n"
  result &= "  alpha3: " & $c.alpha3 & "\n"
  result &= "}"

template proj(x: auto) =
  for e in x:
    x[e].projectU x[e]

template proj(r: auto, x: auto) =
  for e in r:
    r[e].projectU x[e]

template projDeriv(r: auto, x: auto, c: auto) =
  for i in r:
    r[i].projectUDeriv(x[i], c[i])

template projDeriv(r: auto, u, x: auto, c: auto) =
  for i in r:
    r[i].projectUDeriv(u[i], x[i], c[i])

template projVJP(r: auto, x: auto, c: auto) =
  ## Alias for projDeriv (kernel VJP).
  projDeriv(r, x, c)

template projVJP(r: auto, u, x: auto, c: auto) =
  ## Alias for projDeriv (kernel VJP) with explicit projected link.
  projDeriv(r, u, x, c)

template projHVP(r: auto, x: auto, c: auto, dx: auto) =
  ## Gauge-field wrapper for the directional second derivative of the unitary
  ## projection backprop.  Applies projectUHVP site-by-site.
  for i in r:
    r[i].projectUHVP(x[i], c[i], dx[i])

template projHVPu(r: auto, u, x: auto, c: auto, dx: auto) =
  ## Variant taking both the projected link u and the pre-projection x.
  for i in r:
    r[i].projectUHVPu(u[i], x[i], c[i], dx[i])

template projHVP(r: auto, x: auto, c: auto, dx: auto, dc: auto) =
  ## Variant that also propagates a chain tangent dc alongside dx.
  for i in r:
    r[i].projectUHVP(x[i], c[i], dx[i], dc[i])

template projHVPu(r: auto, u, x: auto, c: auto, dx: auto, dc: auto) =
  ## Variant taking both projected link u and upstream tangent dc.
  for i in r:
    r[i].projectUHVPu(u[i], x[i], c[i], dx[i], dc[i])

template projJVP(r: auto, u, x: auto, dx: auto) =
  ## Forward tangent of unitary projection at field level.
  for i in r:
    r[i].projectUJVP(u[i], x[i], dx[i])

template projJVP(r: auto, x: auto, dx: auto) =
  ## Forward tangent without providing u.
  for i in r:
    r[i].projectUJVP(x[i], dx[i])

template projVJPChain(chainbar: auto, u, x: auto, rbar: auto) =
  ## Adjoint of projVJP with respect to the chain input.
  ## Given rbar (the adjoint seed), computes chainbar.
  for i in chainbar:
    chainbar[i].projectUVJPChain(u[i], x[i], rbar[i])

template projVJPChain(chainbar: auto, x: auto, rbar: auto) =
  ## Variant without providing u.
  for i in chainbar:
    chainbar[i].projectUVJPChain(x[i], rbar[i])

template projHVPVJP_dx(dxbar: auto, u, x: auto, c: auto, rbar: auto) =
  ## Adjoint of projHVP with respect to dx.
  ## Given rbar (the adjoint seed), computes dxbar.
  for i in dxbar:
    dxbar[i].projectUHVPVJP_dx(u[i], x[i], c[i], rbar[i])

template projHVPVJP_dx(dxbar: auto, x: auto, c: auto, rbar: auto) =
  ## Variant without providing u (computes u = projectU(x) internally).
  for i in dxbar:
    dxbar[i].projectUHVPVJP_dx(x[i], c[i], rbar[i])

template projHVPVJP_dc(dcbar: auto, u, x: auto, c: auto, rbar: auto) =
  ## Adjoint of projHVP with respect to dc (chain tangent).
  ## Given rbar (the adjoint seed), computes dcbar.
  for i in dcbar:
    dcbar[i].projectUHVPVJP_dc(u[i], x[i], c[i], rbar[i])

template projHVPVJP_dc(dcbar: auto, x: auto, c: auto, rbar: auto) =
  ## Variant without providing u (computes u = projectU(x) internally).
  for i in dcbar:
    dcbar[i].projectUHVPVJP_dc(x[i], c[i], rbar[i])

template forPairs(mu, nu: untyped; body: untyped) =
  for mu in 0..<4:
    for nu in 0..<4:
      if nu != mu:
        body

template forNu(mu, nu: untyped; body: untyped) =
  for nu in 0..<4:
    if nu != mu:
      body

template forThird(mu, nu, a, b: untyped; body: untyped) =
  for a in 0..<4:
    if a != mu and a != nu:
      let b = 1 + 2 + 3 - mu - nu - a
      discard b
      body

# L[mu][nu] = P( (1-a1)*g[mu] + 0.5*a1 SS(g[nu],g[mu]) )
# L2[mu][nu] = P( (1-a2)*g[mu] + 0.25*a2 sum{a,b!=mu,nu} SS(L[a][b],L[mu][b]) )
# fl[mu] = P( (1-a3)*g[mu] + a3/6 sum{nu!=mu} SS(L2[nu][mu],L2[mu][nu]) )
#proc smear*(coef: HypCoefs, gf: auto, fl: auto, ht: HypTemps,
#            info: var PerfInfo) =
proc smearGetForce*[G](coef: HypCoefs, gf: G, fl: G,
            info: var PerfInfo):auto =
  ## Note that the resulting proc, smearedForce, holds a reference to the input gauge gf.
  ## The correctness of the algorithm depends on gf remaining the same.
  ## On the contrary, any changes to the smeared gauge fl would have no effects to the force calculation.
  let lo = gf[0].l
  var hs = lo.newHypSmear(coef)
  hs.smear(gf, fl)
  discard info
  proc smearedForce(f, chain: G) =
    hs.smearVJP(f, chain)
  smearedForce

proc smearPriv[G](coef: HypCoefs, gf: G, fl: G, info: var PerfInfo) {.codegenDecl:
    "__attribute__((noinline)) $# $#$#".} =
  # Avoid inlining or other compiler optimizations
  # in order to guarantee the change of the stack pointer,
  # such that Nim's GC is able to collect the memory.
  {.emit: "asm (\"\");".}
  let lo = gf[0].l
  var hs = lo.newHypSmear(coef)
  hs.smear(gf, fl)
  discard info
proc smear*[G](coef: HypCoefs, gf: G, fl: G, info: var PerfInfo) =
  ## Try our best to release memory here.
  ## Sometimes it still requires a GC after this function returns.
  coef.smearPriv(gf, fl, info)
  qexGC()

#proc smear*(c: HypCoefs, gf: auto, fl: auto, info: var PerfInfo) =
#  var t = newHypTemps(gf)
#  smear(c, gf, fl, t, info)

proc smear*(c: HypCoefs, g: auto, fl: auto) =
  var info: PerfInfo
  c.smear(g, fl, info)

#proc deriv*(coef: HypCoefs, gf: auto, fl: auto, info: var PerfInfo) =
#  ## Compatibility wrapper around the shared forward HYP builder.
#  coef.smear(gf, fl, info)

# Preallocated HYP smearing environment, following stoutsmear style
type
  HypForwardState[L,F,T] = object
    ## Persistent state produced by `smear` and consumed by derivative passes.
    gf*: array[4,F]
    ## Dedicated shifters for pair terms built directly from the original gauge
    ## links. These buffers stay stable after `smear`, so later passes do not
    ## need to restore any hidden cached field pointers.
    gaugePairShift*: array[4,array[4,Shifter[F,T]]]
    ## Separate shifters for projected Level-1/Level-2 pair terms. Keeping
    ## these distinct from `gaugePairShift` avoids mutating a single shifter set
    ## between gauge and projected-link phases.
    stagePairShift*: array[4,array[4,Shifter[F,T]]]
    sm1*: array[4,Shifter[F,T]]
    l1x*, l2x*: FieldArray[2,L.V,T]
    when keepProj:
      l1*, l2*: FieldArray[2,L.V,T]
    flx*: FieldArray[1,L.V,T]

  HypScratch[L,F,T] = object
    ## Mutable work buffers used by reverse/tangent passes.
    ## Post-projection chain values for the output, Level 1, and Level 2
    ## reverse passes. These are scratch, not cached forward state.
    outputChain*: FieldArray[1,L.V,T]
    level1Chain*, level2Chain*: FieldArray[2,L.V,T]
    ## Pre-projection chain values saved specifically for second-order
    ## projection adjoints.
    preProjLevel1Chain*, preProjLevel2Chain*: FieldArray[2,L.V,T]
    ## Reusable forward shifters indexed by direction. `chainShift` is used for
    ## active chain fields, while `seedShift` stays available as a second
    ## scratch pool when the higher-order adjoints need an extra shifted seed.
    chainShift*: array[4,Shifter[F,T]]
    seedShift*: array[4,Shifter[F,T]]
    tm1*, tm2*: F

  HypSmear*[L,F,T] = object
    ## Reusable buffers for HYP smearing and its derivative.
    ##
    ## The forward state is intentionally separated from the mutable scratch so
    ## it is clear which fields are part of the cached smear result and which
    ## are transient workspaces for reverse/tangent evaluations.
    coef*: HypCoefs
    state*: HypForwardState[L,F,T]
    work*: HypScratch[L,F,T]

  HypTangentState[GaugeTangent, OutputTangent, PairField, PairShift] = object
    dgEff*: GaugeTangent
    dl1x*, dl1*: PairField
    dl2x*, dl2*: PairField
    dflx*: OutputTangent
    dsl2*: PairShift

  HypHvpWorkspace[GaugeField, PairField] = object
    dOutputChain*: GaugeField
    dLevel1Chain*, dLevel2Chain*: PairField
    dPreProjLevel1Chain*, dPreProjLevel2Chain*: PairField

  HypChainAdjointWorkspace[GaugeField, PairField] = object
    outputChainBar*: GaugeField
    level1ChainBar*, level2ChainBar*: PairField

  HypHvpAdjointWorkspace[GaugeField, PairField] = object
    dflxbar*, dfcbar*, dgEffbar*: GaugeField
    dfl1bar*, dfl2bar*: PairField
    dl1bar*, dl2bar*: PairField
    dl1xbar*, dl2xbar*: PairField
    dfl1barL2*, dfl2barL3*: PairField

  HypLevel3HvpScratch[PairField, ShiftField] = object
    shiftedOutputChain*: ShiftField
    acc1*, acc2*: PairField

  HypLevel3AdjointScratch[GaugeField] = object
    dcbar*: GaugeField

proc scaleFactors[L,F,T](hs: HypSmear[L,F,T]): HypScaleFactors {.inline.} =
  result = HypScaleFactors(
    alp1: hs.coef.alpha1 / 2.0,
    alp2: hs.coef.alpha2 / 4.0,
    alp3: hs.coef.alpha3 / 6.0,
    ma1: 1 - hs.coef.alpha1,
    ma2: 1 - hs.coef.alpha2,
    ma3: 1 - hs.coef.alpha3,
  )

proc newHypSmear*(l: Layout, coef: HypCoefs): auto =
  ## Allocate a reusable environment for HYP smearing and derivatives
  type F = type(l.newGauge[0])
  type T = type(l.newGauge[0][0])
  type L = type(l)
  var hs: HypSmear[L,F,T]
  hs.coef = coef
  hs.state.l1x = newFieldArray2(l,F,[4,4],mu!=nu)
  hs.state.l2x = newFieldArray2(l,F,[4,4],mu!=nu)
  when keepProj:
    hs.state.l1 = newFieldArray2(l,F,[4,4],mu!=nu)
    hs.state.l2 = newFieldArray2(l,F,[4,4],mu!=nu)
  hs.state.flx = newFieldArray(l,F,4)
  hs.work.outputChain = newFieldArray(l,F,4)
  hs.work.level1Chain = newFieldArray2(l,F,[4,4],mu!=nu)
  hs.work.level2Chain = newFieldArray2(l,F,[4,4],mu!=nu)
  hs.work.preProjLevel1Chain = newFieldArray2(l,F,[4,4],mu!=nu)
  hs.work.preProjLevel2Chain = newFieldArray2(l,F,[4,4],mu!=nu)
  hs.work.tm1.new(l)
  hs.work.tm2.new(l)
  for mu in 0..<4:
    hs.state.sm1[mu] = newShifter(hs.work.outputChain[mu], mu, -1)
    hs.work.chainShift[mu] = newShifter(hs.work.outputChain[mu], mu, 1)
    hs.work.seedShift[mu] = newShifter(hs.work.outputChain[mu], mu, 1)
    for nu in 0..<4:
      if nu != mu:
        hs.state.gaugePairShift[mu][nu] = newShifter(hs.work.outputChain[mu], nu, 1)
        hs.state.stagePairShift[mu][nu] = newShifter(hs.work.outputChain[mu], nu, 1)
  hs

proc initProjectedLevelsImpl(hs: var HypSmear; l1, l2: var auto;
                             populateFromCache: static bool) =
  when keepProj:
    l1 = hs.state.l1
    l2 = hs.state.l2
  else:
    let lo = hs.state.gf[0].l
    type Fld = typeof(hs.state.gf[0])
    l1 = newFieldArray2(lo, Fld, [4,4], mu!=nu)
    l2 = newFieldArray2(lo, Fld, [4,4], mu!=nu)
    when populateFromCache:
      threads:
        forPairs(mu, nu):
          l1[mu,nu].proj hs.state.l1x[mu,nu]
          l2[mu,nu].proj hs.state.l2x[mu,nu]

proc initProjectedLevels(hs: var HypSmear; l1, l2: var auto) =
  initProjectedLevelsImpl(hs, l1, l2, true)

proc buildGaugeTangent(gf, dg: auto): auto =
  let lo = gf[0].l
  var dgEff = lo.newGauge
  threads:
    for mu in 0..<gf.len:
      dgEff[mu] := dg[mu] * gf[mu]
  dgEff

proc initShiftedTangents(f: auto): auto =
  var s: array[4,array[4,Shifter[type(f[0]),type(f[0][0])]]]
  forPairs(mu, nu):
    s[mu][nu] = newShifter(f[mu], nu, 1)
  threads:
    forPairs(mu, nu):
      discard s[mu][nu] ^*! f[mu]
  s

proc initForwardProjectedLevels(hs: var HypSmear; l1, l2: var auto) =
  ## Allocate the projected Level-1/Level-2 buffers used while rebuilding the
  ## forward HYP state. Unlike `initProjectedLevels`, this does not populate the
  ## buffers from cached state.
  initProjectedLevelsImpl(hs, l1, l2, false)

template projectAndScalePairChain(chain, projected, xfield, preProjChain,
                                  alp, derivSink, ma: untyped;
                                  accumulateDirect: static bool) =
  threads:
    forPairs(mu, nu):
      preProjChain[mu,nu] := chain[mu,nu]
      chain[mu,nu].projVJP(projected[mu,nu], xfield[mu,nu], preProjChain[mu,nu])

  threads:
    forPairs(mu, nu):
      when accumulateDirect:
        derivSink[mu] += ma * chain[mu,nu]
      chain[mu,nu] *= alp

template clearReverseChains(outputChain, level1Chain, level2Chain: auto) =
  threads:
    for mu in 0..<4:
      outputChain[mu] := 0
      forNu(mu, nu):
        level1Chain[mu,nu] := 0
        level2Chain[mu,nu] := 0

template clearTangentWorkspace(dflx, dl1x, dl2x, dl1, dl2: auto) =
  threads:
    for mu in 0..<4:
      dflx[mu] := 0
      forNu(mu, nu):
        dl1x[mu,nu] := 0
        dl2x[mu,nu] := 0
        dl1[mu,nu] := 0
        dl2[mu,nu] := 0

template clearHvpWorkspace(dderiv, dOutputChain, dLevel1Chain, dLevel2Chain: auto) =
  threads:
    for mu in 0..<4:
      dderiv[mu] := 0
      dOutputChain[mu] := 0
      forNu(mu, nu):
        dLevel1Chain[mu,nu] := 0
        dLevel2Chain[mu,nu] := 0

template clearChainAdjointWorkspace(outputChainBar, chainbar,
                                    level1ChainBar, level2ChainBar: auto) =
  threads:
    for mu in 0..<4:
      outputChainBar[mu] := 0
      chainbar[mu] := 0
      forNu(mu, nu):
        level1ChainBar[mu,nu] := 0
        level2ChainBar[mu,nu] := 0

template clearHvpAdjointWorkspace(gbar, dflxbar, dfcbar, dgEffbar,
                                  dfl1bar, dfl2bar, dl1bar, dl2bar,
                                  dl1xbar, dl2xbar: auto) =
  threads:
    for mu in 0..<4:
      gbar[mu] := 0
      dflxbar[mu] := 0
      dfcbar[mu] := 0
      dgEffbar[mu] := 0
      forNu(mu, nu):
        dfl1bar[mu,nu] := 0
        dfl2bar[mu,nu] := 0
        dl1bar[mu,nu] := 0
        dl2bar[mu,nu] := 0
        dl1xbar[mu,nu] := 0
        dl2xbar[mu,nu] := 0

proc newHvpWorkspace(gf: auto): auto =
  let lo = gf[0].l
  type Fld = typeof(gf[0])
  var dOutputChain = lo.newGauge
  var dLevel1Chain = newFieldArray2(lo, Fld, [4,4], mu!=nu)
  var dLevel2Chain = newFieldArray2(lo, Fld, [4,4], mu!=nu)
  var dPreProjLevel1Chain = newFieldArray2(lo, Fld, [4,4], mu!=nu)
  var dPreProjLevel2Chain = newFieldArray2(lo, Fld, [4,4], mu!=nu)
  HypHvpWorkspace[type(dOutputChain), type(dLevel1Chain)](
    dOutputChain: dOutputChain,
    dLevel1Chain: dLevel1Chain,
    dLevel2Chain: dLevel2Chain,
    dPreProjLevel1Chain: dPreProjLevel1Chain,
    dPreProjLevel2Chain: dPreProjLevel2Chain,
  )

proc newChainAdjointWorkspace(chainbar, gf: auto): auto =
  let lo = gf[0].l
  type Fld = typeof(gf[0])
  var outputChainBar = lo.newGauge
  var level1ChainBar = newFieldArray2(lo, Fld, [4,4], mu!=nu)
  var level2ChainBar = newFieldArray2(lo, Fld, [4,4], mu!=nu)
  clearChainAdjointWorkspace(outputChainBar, chainbar, level1ChainBar, level2ChainBar)
  HypChainAdjointWorkspace[type(outputChainBar), type(level1ChainBar)](
    outputChainBar: outputChainBar,
    level1ChainBar: level1ChainBar,
    level2ChainBar: level2ChainBar,
  )

proc newHvpAdjointWorkspace(gbar, gf: auto): auto =
  let lo = gf[0].l
  type Fld = typeof(gf[0])
  var dflxbar = lo.newGauge
  var dfcbar = lo.newGauge
  var dgEffbar = lo.newGauge
  var dfl1bar = newFieldArray2(lo, Fld, [4,4], mu!=nu)
  var dfl2bar = newFieldArray2(lo, Fld, [4,4], mu!=nu)
  var dl1bar = newFieldArray2(lo, Fld, [4,4], mu!=nu)
  var dl2bar = newFieldArray2(lo, Fld, [4,4], mu!=nu)
  var dl1xbar = newFieldArray2(lo, Fld, [4,4], mu!=nu)
  var dl2xbar = newFieldArray2(lo, Fld, [4,4], mu!=nu)
  clearHvpAdjointWorkspace(gbar, dflxbar, dfcbar, dgEffbar,
                           dfl1bar, dfl2bar, dl1bar, dl2bar, dl1xbar, dl2xbar)
  var dfl1barL2 = newFieldArray2(lo, Fld, [4,4], mu!=nu)
  var dfl2barL3 = newFieldArray2(lo, Fld, [4,4], mu!=nu)
  HypHvpAdjointWorkspace[type(dflxbar), type(dfl1bar)](
    dflxbar: dflxbar,
    dfcbar: dfcbar,
    dgEffbar: dgEffbar,
    dfl1bar: dfl1bar,
    dfl2bar: dfl2bar,
    dl1bar: dl1bar,
    dl2bar: dl2bar,
    dl1xbar: dl1xbar,
    dl2xbar: dl2xbar,
    dfl1barL2: dfl1barL2,
    dfl2barL3: dfl2barL3,
  )

proc newLevel3HvpScratch(gf, dOutputChain: auto): auto =
  let lo = gf[0].l
  type Fld = typeof(gf[0])
  var shiftedOutputChain: array[4, array[4, Shifter[Fld,type(gf[0][0])]]]
  forPairs(mu, nu):
    shiftedOutputChain[mu][nu] = newShifter(dOutputChain[mu], nu, 1)
  var acc1 = newFieldArray2(lo, Fld, [4,4], mu!=nu)
  var acc2 = newFieldArray2(lo, Fld, [4,4], mu!=nu)
  HypLevel3HvpScratch[type(acc1), type(shiftedOutputChain)](
    shiftedOutputChain: shiftedOutputChain,
    acc1: acc1,
    acc2: acc2,
  )

proc newLevel3AdjointScratch(gf: auto): auto =
  let lo = gf[0].l
  var dcbar = lo.newGauge
  threads:
    for mu in 0..<4:
      dcbar[mu] := 0
  HypLevel3AdjointScratch[type(dcbar)](dcbar: dcbar)

template primeLevel3PairShifts(hs, l2: auto) =
  threads:
    forPairs(mu, nu):
      discard hs.state.stagePairShift[mu][nu] ^*! l2[mu,nu]

template forLevel3Pairs(hs, l2, pl1, pl2: untyped; body: untyped) =
  ## Prime the Level-3 pair shifters once, then iterate over the paired links
  ## as `(pl1, pl2) = (l2[nu,mu], l2[mu,nu])`.
  primeLevel3PairShifts(hs, l2)

  threads:
    forPairs(mu, nu):
      let pl1 = l2[nu,mu]
      let pl2 = l2[mu,nu]
      body

template buildLevel3ReverseChains(hs, l2, outputChain, level2Chain, tm1, tm2: auto) =
  forLevel3Pairs(hs, l2, pl1, pl2):
    threadBarrier()
    discard hs.work.chainShift[nu] ^*! outputChain[mu]
    threadBarrier()
    symStapleVJP(level2Chain[nu,mu], level2Chain[mu,nu],
                 pl1, pl2,
                 hs.state.stagePairShift[nu][mu], hs.state.stagePairShift[mu][nu],
                 outputChain[mu], hs.work.chainShift[nu],
                 tm1, tm2, hs.state.sm1[nu], hs.state.sm1[mu])
    threadBarrier()

template buildForwardPairStage(Stage: static int, hs, pairSrc, outX, outProj,
                               ma, alp, tm1: untyped) =
  when Stage == 1:
    threads:
      forPairs(mu, nu):
        discard hs.state.gaugePairShift[mu][nu] ^*! hs.state.gf[mu]

  threads:
    forPairs(mu, nu):
      outX[mu,nu] := ma * hs.state.gf[mu]
      when Stage == 1:
        symStaple(outX[mu,nu], alp, hs.state.gf[nu], hs.state.gf[mu],
                  hs.state.gaugePairShift[nu][mu], hs.state.gaugePairShift[mu][nu],
                  tm1, hs.state.sm1[nu])
        threadBarrier()
      else:
        forThird(mu, nu, a, b):
          let pl1 = pairSrc[a,b]
          let pl2 = pairSrc[mu,b]
          threadBarrier()
          discard hs.state.stagePairShift[nu][mu] ^*! pl1
          discard hs.state.stagePairShift[mu][a] ^*! pl2
          threadBarrier()
          symStaple(outX[mu,nu], alp, pl1, pl2,
                    hs.state.stagePairShift[nu][mu], hs.state.stagePairShift[mu][a],
                    tm1, hs.state.sm1[a])
          threadBarrier()
      outProj[mu,nu].proj outX[mu,nu]

template buildTangentPairStage(Stage: static int, hs, pairSrc, tangentSrc,
                               shiftedTangents, preProjOut, tangentOut,
                               preProjBase, projBase, dgEff, ma, alp, tm1: untyped) =
  when Stage == 1:
    threads:
      forPairs(mu, nu):
        discard hs.state.gaugePairShift[mu][nu] ^*! hs.state.gf[mu]

  threads:
    forPairs(mu, nu):
      preProjOut[mu,nu] := ma * dgEff[mu]
      when Stage == 1:
        symStapleJVP(preProjOut[mu,nu], alp,
                     hs.state.gf[nu], hs.state.gf[mu],
                     hs.state.gaugePairShift[nu][mu], hs.state.gaugePairShift[mu][nu],
                     tangentSrc[nu], tangentSrc[mu],
                     shiftedTangents[nu][mu], shiftedTangents[mu][nu],
                     tm1, hs.state.sm1[nu])
        threadBarrier()
      else:
        forThird(mu, nu, a, b):
          let pl1 = pairSrc[a,b]
          let pl2 = pairSrc[mu,b]
          let dpl1 = tangentSrc[a,b]
          let dpl2 = tangentSrc[mu,b]
          threadBarrier()
          discard hs.state.stagePairShift[nu][mu] ^*! pl1
          discard hs.state.stagePairShift[mu][a] ^*! pl2
          discard shiftedTangents[mu][nu][a][0] ^*! dpl1
          discard shiftedTangents[mu][nu][a][1] ^*! dpl2
          threadBarrier()
          symStapleJVP(preProjOut[mu,nu], alp,
                       pl1, pl2, hs.state.stagePairShift[nu][mu], hs.state.stagePairShift[mu][a],
                       dpl1, dpl2,
                       shiftedTangents[mu][nu][a][0], shiftedTangents[mu][nu][a][1],
                       tm1, hs.state.sm1[a])
          threadBarrier()
      tangentOut[mu,nu].projJVP(projBase[mu,nu], preProjBase[mu,nu], preProjOut[mu,nu])

template applyReversePairStage(Stage: static int, hs, pairSrc, chainIn, dst,
                               tm1, tm2: untyped) =
  threads:
    forPairs(mu, nu):
      when Stage == 1:
        discard hs.work.chainShift[nu] ^* chainIn[mu,nu]
        symStapleVJP(dst[nu], dst[mu],
                     hs.state.gf[nu], hs.state.gf[mu],
                     hs.state.gaugePairShift[nu][mu], hs.state.gaugePairShift[mu][nu],
                     chainIn[mu,nu], hs.work.chainShift[nu],
                     tm1, tm2, hs.state.sm1[nu], hs.state.sm1[mu])
        threadBarrier()
      else:
        forThird(mu, nu, a, b):
          let pl1 = pairSrc[a,b]
          let pl2 = pairSrc[mu,b]
          threadBarrier()
          discard hs.state.stagePairShift[nu][mu] ^*! pl1
          discard hs.state.stagePairShift[mu][a] ^*! pl2
          discard hs.work.chainShift[a] ^*! chainIn[mu,nu]
          threadBarrier()
          symStapleVJP(dst[a,b], dst[mu,b],
                       pl1, pl2,
                       hs.state.stagePairShift[nu][mu], hs.state.stagePairShift[mu][a],
                       chainIn[mu,nu], hs.work.chainShift[a],
                       tm1, tm2, hs.state.sm1[a], hs.state.sm1[mu])
          threadBarrier()

template applyReversePairStageChainAdjoint(Stage: static int, hs, pairSrc,
                                           dst, seedBar, tm1, tm2: untyped) =
  threads:
    forPairs(mu, nu):
      when Stage == 1:
        symStapleVJPChain(dst[mu,nu],
                          hs.state.gf[nu], hs.state.gf[mu],
                          hs.state.gaugePairShift[nu][mu], hs.state.gaugePairShift[mu][nu],
                          seedBar[nu], seedBar[mu],
                          hs.work.chainShift[nu], hs.work.seedShift[mu],
                          tm1, tm2, hs.state.sm1[nu])
        threadBarrier()
      else:
        forThird(mu, nu, a, b):
          let pl1 = pairSrc[a,b]
          let pl2 = pairSrc[mu,b]
          threadBarrier()
          discard hs.state.stagePairShift[nu][mu] ^*! pl1
          discard hs.state.stagePairShift[mu][a] ^*! pl2
          threadBarrier()
          symStapleVJPChain(dst[mu,nu],
                            pl1, pl2,
                            hs.state.stagePairShift[nu][mu], hs.state.stagePairShift[mu][a],
                            seedBar[a,b], seedBar[mu,b],
                            hs.work.chainShift[a], hs.work.seedShift[mu],
                            tm1, tm2, hs.state.sm1[a])
          threadBarrier()

template accumulateAndScaleGauge(dst, src, ma, alp: auto) =
  threads:
    for mu in 0..<4:
      dst[mu] += ma * src[mu]
      src[mu] *= alp

template accumulateAndScalePairs(dst, src, ma, alp: auto) =
  threads:
    forPairs(mu, nu):
      dst[mu] += ma * src[mu,nu]
      src[mu,nu] *= alp

template addScaledGauge(dst, src, scale: auto) =
  threads:
    for mu in 0..<4:
      dst[mu] += scale * src[mu]

template addScaledPairs(dst, src, scale: auto) =
  threads:
    forPairs(mu, nu):
      dst[mu,nu] += scale * src[mu,nu]

template applyPairProjectionHvpPhase(chainBar, dChainBar, projected, xfield,
                                     preProjChain, tangent, derivHvp,
                                     ma, alp: untyped) =
  applyProjectionHVP(chainBar, dChainBar, projected, xfield, preProjChain, tangent)
  accumulateAndScalePairs(derivHvp, chainBar, ma, alp)

template applyPairStapleHvpPhase(Stage: static int, hs, pairSrc, tangentSrc,
                                 chainSrc, dChainIn, dst, tm1, tm2: untyped) =
  block:
    let lo = hs.state.gf[0].l
    type Fld = typeof(hs.state.gf[0])
    when Stage == 1:
      var shiftedPairTangents: array[4, array[4, array[2, Shifter[Fld,type(hs.state.gf[0][0])]]]]
      var shiftedChainTangents: array[4, array[4, Shifter[Fld,type(hs.state.gf[0][0])]]]
      forPairs(mu, nu):
        shiftedPairTangents[mu][nu][0] = newShifter(tangentSrc[nu], mu, 1)
        shiftedPairTangents[mu][nu][1] = newShifter(tangentSrc[mu], nu, 1)
        shiftedChainTangents[mu][nu] = newShifter(dChainIn[mu,nu], nu, 1)
      var acc1 = newFieldArray2(lo, Fld, [4,4], mu!=nu)
      var acc2 = newFieldArray2(lo, Fld, [4,4], mu!=nu)

      threads:
        forPairs(mu, nu):
          discard hs.work.chainShift[nu] ^* chainSrc[mu,nu]
          discard shiftedPairTangents[mu][nu][0] ^*! tangentSrc[nu]
          discard shiftedPairTangents[mu][nu][1] ^*! tangentSrc[mu]
          discard shiftedChainTangents[mu][nu] ^*! dChainIn[mu,nu]
          symStapleHVP(dst[nu], dst[mu], acc1[mu,nu], acc2[mu,nu],
                       hs.state.gf[nu], hs.state.gf[mu],
                       hs.state.gaugePairShift[nu][mu], hs.state.gaugePairShift[mu][nu],
                       chainSrc[mu,nu], hs.work.chainShift[nu],
                       tangentSrc[nu], tangentSrc[mu], dChainIn[mu,nu],
                       shiftedPairTangents[mu][nu][0], shiftedPairTangents[mu][nu][1],
                       shiftedChainTangents[mu][nu],
                       tm1, tm2, hs.state.sm1[nu], hs.state.sm1[mu])
          threadBarrier()
    else:
      var shiftedPairTangents: array[4, array[4, array[4, array[2, Shifter[Fld,type(hs.state.gf[0][0])]]]]]
      var shiftedChainTangents: array[4, array[4, array[4, Shifter[Fld,type(hs.state.gf[0][0])]]]]
      forPairs(mu, nu):
        forThird(mu, nu, a, b):
          shiftedPairTangents[mu][nu][a][0] = newShifter(tangentSrc[a,b], mu, 1)
          shiftedPairTangents[mu][nu][a][1] = newShifter(tangentSrc[mu,b], a, 1)
          shiftedChainTangents[mu][nu][a] = newShifter(dChainIn[mu,nu], a, 1)
      var acc1: array[4, array[4, array[4, Fld]]]
      var acc2: array[4, array[4, array[4, Fld]]]
      forPairs(mu, nu):
        forThird(mu, nu, a, b):
          acc1[mu][nu][a] = lo.newField(type(hs.state.gf[0][0]))
          acc2[mu][nu][a] = lo.newField(type(hs.state.gf[0][0]))

      threads:
        forPairs(mu, nu):
          forThird(mu, nu, a, b):
            let pl1 = pairSrc[a,b]
            let pl2 = pairSrc[mu,b]
            threadBarrier()
            discard hs.state.stagePairShift[nu][mu] ^*! pl1
            discard hs.state.stagePairShift[mu][a] ^*! pl2
            discard hs.work.chainShift[a] ^*! chainSrc[mu,nu]
            discard shiftedPairTangents[mu][nu][a][0] ^*! tangentSrc[a,b]
            discard shiftedPairTangents[mu][nu][a][1] ^*! tangentSrc[mu,b]
            discard shiftedChainTangents[mu][nu][a] ^*! dChainIn[mu,nu]
            threadBarrier()
            symStapleHVP(dst[a,b], dst[mu,b],
                         acc1[mu][nu][a], acc2[mu][nu][a],
                         pl1, pl2,
                         hs.state.stagePairShift[nu][mu], hs.state.stagePairShift[mu][a],
                         chainSrc[mu,nu], hs.work.chainShift[a],
                         tangentSrc[a,b], tangentSrc[mu,b], dChainIn[mu,nu],
                         shiftedPairTangents[mu][nu][a][0], shiftedPairTangents[mu][nu][a][1],
                         shiftedChainTangents[mu][nu][a],
                         tm1, tm2, hs.state.sm1[a], hs.state.sm1[mu])
            threadBarrier()

template backpropPairHvpOutputPhase(Stage: static int, hs, pairSrc, chainSrc,
                                    seedBar, pairBarOut, chainBar,
                                    nextChainBar, projected, xfield,
                                    preProjChain, alp, tm1, tm2: untyped) =
  block:
    let lo = hs.state.gf[0].l
    type Fld = typeof(hs.state.gf[0])
    var dcbar = newFieldArray2(lo, Fld, [4,4], mu!=nu)
    threads:
      forPairs(mu, nu):
        dcbar[mu,nu] := 0

    threads:
      forPairs(mu, nu):
        when Stage == 1:
          discard hs.work.chainShift[nu] ^* chainSrc[mu,nu]
          symStapleHVPVJP(pairBarOut[nu], pairBarOut[mu], dcbar[mu,nu],
                          hs.state.gf[nu], hs.state.gf[mu],
                          hs.state.gaugePairShift[nu][mu], hs.state.gaugePairShift[mu][nu],
                          chainSrc[mu,nu], hs.work.chainShift[nu],
                          seedBar[nu], seedBar[mu],
                          hs.work.seedShift[nu], hs.work.seedShift[mu],
                          tm1, tm2, hs.state.sm1[nu], hs.state.sm1[mu])
          threadBarrier()
        else:
          forThird(mu, nu, a, b):
            let pl1 = pairSrc[a,b]
            let pl2 = pairSrc[mu,b]
            threadBarrier()
            discard hs.state.stagePairShift[nu][mu] ^*! pl1
            discard hs.state.stagePairShift[mu][a] ^*! pl2
            discard hs.work.chainShift[a] ^*! chainSrc[mu,nu]
            threadBarrier()
            symStapleHVPVJP(pairBarOut[a,b], pairBarOut[mu,b], dcbar[mu,nu],
                            pl1, pl2,
                            hs.state.stagePairShift[nu][mu], hs.state.stagePairShift[mu][a],
                            chainSrc[mu,nu], hs.work.chainShift[a],
                            seedBar[a,b], seedBar[mu,b],
                            hs.work.seedShift[a], hs.work.seedShift[mu],
                            tm1, tm2, hs.state.sm1[a], hs.state.sm1[mu])
            threadBarrier()

    addScaledPairs(chainBar, dcbar, alp)

    threads:
      forPairs(mu, nu):
        nextChainBar[mu,nu].projHVPVJP_dc(
          projected[mu,nu], xfield[mu,nu], preProjChain[mu,nu],
          chainBar[mu,nu]
        )

template backpropPairTangentPhase(Stage: static int, hs, pairSrc, tangentXBar,
                                  pairBarOut, dgEffbar, ma, alp, tm1, tm2: untyped) =
  threads:
    forPairs(mu, nu):
      dgEffbar[mu] += ma * tangentXBar[mu,nu]

  threads:
    forPairs(mu, nu):
      when Stage == 1:
        symStapleJVPVJP(dgEffbar[nu], dgEffbar[mu], alp,
                        hs.state.gf[nu], hs.state.gf[mu],
                        hs.state.gaugePairShift[nu][mu], hs.state.gaugePairShift[mu][nu],
                        tangentXBar[mu,nu], hs.work.seedShift[nu],
                        tm1, tm2, hs.state.sm1[nu], hs.state.sm1[mu])
        threadBarrier()
      else:
        forThird(mu, nu, a, b):
          let pl1 = pairSrc[a,b]
          let pl2 = pairSrc[mu,b]
          threadBarrier()
          discard hs.state.stagePairShift[nu][mu] ^*! pl1
          discard hs.state.stagePairShift[mu][a] ^*! pl2
          threadBarrier()
          symStapleJVPVJP(pairBarOut[a,b], pairBarOut[mu,b], alp,
                          pl1, pl2,
                          hs.state.stagePairShift[nu][mu], hs.state.stagePairShift[mu][a],
                          tangentXBar[mu,nu], hs.work.seedShift[a],
                          tm1, tm2, hs.state.sm1[a], hs.state.sm1[mu])
          threadBarrier()

proc buildTangentState[L,F,T](hs: var HypSmear[L,F,T], dg: auto): auto =
  ## Build the shared forward tangent state used by `smearJVP` and `smearHVP`.
  let sc = hs.scaleFactors
  var l1: type(hs.state.l1x)
  var l2: type(hs.state.l2x)
  initProjectedLevels(hs, l1, l2)
  let hsp = hs.addr

  let lo = hs.state.gf[0].l
  type Fld = typeof(hs.state.gf[0])
  let dgEff = buildGaugeTangent(hs.state.gf, dg)
  let ds1 = initShiftedTangents(dgEff)

  var dl1x = newFieldArray2(lo, Fld, [4,4], mu!=nu)
  var dl2x = newFieldArray2(lo, Fld, [4,4], mu!=nu)
  var dflx = newFieldArray(lo, Fld, 4)
  var dl1 = newFieldArray2(lo, Fld, [4,4], mu!=nu)
  var dl2 = newFieldArray2(lo, Fld, [4,4], mu!=nu)

  clearTangentWorkspace(dflx, dl1x, dl2x, dl1, dl2)

  buildTangentPairStage(1, hsp, l1, dgEff, ds1, dl1x, dl1, hsp.state.l1x, l1, dgEff, sc.ma1, sc.alp1, hsp.work.tm1)

  var dsl1_l2: array[4, array[4, array[4, array[2, Shifter[Fld,type(hs.state.gf[0][0])]]]]]
  forPairs(mu, nu):
    forThird(mu, nu, a, b):
      dsl1_l2[mu][nu][a][0] = newShifter(dl1[a,b], mu, 1)
      dsl1_l2[mu][nu][a][1] = newShifter(dl1[mu,b], a, 1)

  buildTangentPairStage(2, hsp, l1, dl1, dsl1_l2, dl2x, dl2, hsp.state.l2x, l2, dgEff, sc.ma2, sc.alp2, hsp.work.tm1)

  var dsl2: array[4,array[4,Shifter[Fld,type(hs.state.gf[0][0])]]]
  forPairs(mu, nu):
    dsl2[mu][nu] = newShifter(dl2[mu,nu], nu, 1)
  threads:
    forPairs(mu, nu):
      discard dsl2[mu][nu] ^*! dl2[mu,nu]
      discard hsp.state.stagePairShift[mu][nu] ^*! l2[mu,nu]

  threads:
    for mu in 0..<4:
      dflx[mu] := sc.ma3 * dgEff[mu]
      forNu(mu, nu):
        let pl1 = l2[nu,mu]
        let pl2 = l2[mu,nu]
        let dpl1 = dl2[nu,mu]
        let dpl2 = dl2[mu,nu]
        symStapleJVP(dflx[mu], sc.alp3,
                        pl1, pl2, hsp.state.stagePairShift[nu][mu], hsp.state.stagePairShift[mu][nu],
                        dpl1, dpl2, dsl2[nu][mu], dsl2[mu][nu],
                        hsp.work.tm1, hsp.state.sm1[nu])
        threadBarrier()

  HypTangentState[typeof(dgEff), typeof(dflx), typeof(dl1x), typeof(dsl2)](
    dgEff: dgEff,
    dl1x: dl1x,
    dl1: dl1,
    dl2x: dl2x,
    dl2: dl2,
    dflx: dflx,
    dsl2: dsl2,
  )

template accumulateProjectionAdjoint(dxbar, x, chain0, tangentBar, hvpBar: auto) =
  threads:
    forPairs(mu, nu):
      var tmp {.noinit.}: type(dxbar[mu,nu][0])
      for i in dxbar[mu,nu]:
        tmp.projectUVJP(x[mu,nu][i], tangentBar[mu,nu][i])
        dxbar[mu,nu][i] := tmp

  threads:
    forPairs(mu, nu):
      var tmp {.noinit.}: type(dxbar[mu,nu][0])
      for i in dxbar[mu,nu]:
        tmp.projectUHVPVJP_dx(x[mu,nu][i], chain0[mu,nu][i], hvpBar[mu,nu][i])
        dxbar[mu,nu][i] += tmp

template applyPairProjectionChainAdjoint(chainBar, projected, xfield: auto) =
  let lo = projected[0,1].l
  type Fld = typeof(projected[0,1])
  var chainBarTmp = newFieldArray2(lo, Fld, [4,4], mu!=nu)
  threads:
    forPairs(mu, nu):
      chainBarTmp[mu,nu] := chainBar[mu,nu]
      chainBar[mu,nu] := 0
      chainBar[mu,nu].projVJPChain(projected[mu,nu], xfield[mu,nu], chainBarTmp[mu,nu])

template applyProjectionHVP(dst, dchain0, u, x, chain0, dx: auto) =
  ## Apply the projection HVP for a full level while preserving the incoming
  ## chain tangent in `dchain0`.
  threads:
    forPairs(mu, nu):
      dchain0[mu,nu] := dst[mu,nu]
      dst[mu,nu] := 0
      dst[mu,nu].projHVPu(u[mu,nu], x[mu,nu], chain0[mu,nu],
                          dx[mu,nu], dchain0[mu,nu])

template applyOutputProjectionHVP(dst, x, chain0, dx, dchain: auto) =
  threads:
    for mu in 0..<4:
      if dchain.len > 0:
        dst[mu].projHVP(x[mu], chain0[mu], dx[mu], dchain[mu])
      else:
        dst[mu].projHVP(x[mu], chain0[mu], dx[mu])

template seedSmearHVPAdjoints(dfl1bar, dfl2bar, dfcbar, derivbar,
                              ma1, ma2, ma3: auto) =
  threads:
    for mu in 0..<4:
      dfcbar[mu] += ma3 * derivbar[mu]
      forNu(mu, nu):
        dfl1bar[mu,nu] += ma1 * derivbar[mu]
        dfl2bar[mu,nu] += ma2 * derivbar[mu]

template buildReverseState(hsp, chain, deriv: untyped; accumulateDirect: static bool) =
  ## Build the reverse HYP chain state, optionally accumulating the direct
  ## ma3/ma2/ma1 contributions into `deriv`.
  let sc = hsp[].scaleFactors
  var l1: type(hsp.state.l1x)
  var l2: type(hsp.state.l2x)
  initProjectedLevels(hsp[], l1, l2)

  clearReverseChains(hsp.work.outputChain, hsp.work.level1Chain, hsp.work.level2Chain)

  threads:
    for mu in 0..<4:
      hsp.work.outputChain[mu].projVJP(hsp.state.flx[mu], chain[mu])
      when accumulateDirect:
        deriv[mu] := sc.ma3 * hsp.work.outputChain[mu]
      hsp.work.outputChain[mu] *= sc.alp3

  buildLevel3ReverseChains(hsp, l2, hsp.work.outputChain, hsp.work.level2Chain,
                           hsp.work.tm1, hsp.work.tm2)
  projectAndScalePairChain(hsp.work.level2Chain, l2, hsp.state.l2x,
                           hsp.work.preProjLevel2Chain, sc.alp2, deriv, sc.ma2,
                           accumulateDirect)
  applyReversePairStage(2, hsp, l1, hsp.work.level2Chain, hsp.work.level1Chain,
                        hsp.work.tm1, hsp.work.tm2)
  projectAndScalePairChain(hsp.work.level1Chain, l1, hsp.state.l1x,
                           hsp.work.preProjLevel1Chain, sc.alp1, deriv, sc.ma1,
                           accumulateDirect)

proc smear*(hs: var HypSmear, gf: auto, fl: auto) =
  static: doAssert(type(gf[0]) is typeof(hs.state.gf[0]))
  ## Apply HYP smearing; intermediates stored in hs for derivative
  let sc = hs.scaleFactors
  for mu in 0..<4:
    hs.state.gf[mu] = gf[mu]
  var
    cachedGauge = hs.state.gf
    l1x = hs.state.l1x
    l2x = hs.state.l2x
    flx = hs.state.flx
    tm1 = hs.work.tm1
  var l1: typeof(l1x)
  var l2: typeof(l2x)
  initForwardProjectedLevels(hs, l1, l2)
  let hsp = hs.addr

  buildForwardPairStage(1, hsp, l1, l1x, l1, sc.ma1, sc.alp1, tm1)
  buildForwardPairStage(2, hsp, l1, l2x, l2, sc.ma2, sc.alp2, tm1)
  threads:
    forPairs(mu, nu):
      discard hsp.state.stagePairShift[mu][nu] ^*! l2[mu,nu]
  threads:
    for mu in 0..<4:
      flx[mu] := sc.ma3 * cachedGauge[mu]
      forNu(mu, nu):
        let pl1 = l2[nu,mu]
        let pl2 = l2[mu,nu]
        symStaple(flx[mu], sc.alp3, pl1, pl2,
                  hsp.state.stagePairShift[nu][mu], hsp.state.stagePairShift[mu][nu],
                  tm1, hsp.state.sm1[nu])
        threadBarrier()
      fl[mu].proj flx[mu]

proc smearVJP*(hs: var HypSmear, deriv: auto, chain: auto) =
  static:
    doAssert(type(deriv[0]) is typeof(hs.state.gf[0]))
    doAssert(type(chain[0]) is typeof(hs.state.gf[0]))
  ## Backpropagate derivative through hs.smear
  let hsp = hs.addr
  buildReverseState(hsp, chain, deriv, true)
  applyReversePairStage(1, hsp, hsp.work.level1Chain, hsp.work.level1Chain,
                        deriv, hsp.work.tm1, hsp.work.tm2)

template applyOutputProjectionJvp(dst, projected, preProj, tangent: auto) =
  let lo = preProj[0].l
  var tangentTmp = lo.newGauge
  threads:
    for mu in 0..<4:
      tangentTmp[mu] := tangent[mu]
      dst[mu].projJVP(projected[mu], preProj[mu], tangentTmp[mu])

proc smearJVP*[G](hs: var HypSmear, dsg: G, dg: G, sg: G) =
  ## Forward tangent of HYP smearing.
  ##
  ## Computes: dsg = d/dε|_{ε=0} [ HYP(exp(ε dg) g) ]
  ##
  ## Under multiplicative perturbation exp(ε·dg)·gf, this computes the
  ## tangent of the projected smeared output field.
  ##
  ## Parameters:
  ##   dsg: output - the tangent of the smeared field (projected)
  ##   dg:  input  - the Lie algebra perturbation direction
  ##   sg:  input  - the projected smeared field from hs.smear(g, sg)
  ##
  ## NOTE: hs.smear(g, sg) must be called first to populate internal state.
  static:
    doAssert(type(dsg[0]) is typeof(hs.state.gf[0]))
    doAssert(type(dg[0])  is typeof(hs.state.gf[0]))
  let hsp = hs.addr
  let tangents = hs.buildTangentState(dg)

  applyOutputProjectionJvp(dsg, sg, hsp.state.flx, tangents.dflx)

template applyLevel3StapleHvpPhase(hs, dLevel2Chain, l2, dl2, dsl2,
                                   dOutputChain: untyped) =
  block:
    var scratch = newLevel3HvpScratch(hs.state.gf, dOutputChain)

    forLevel3Pairs(hs, l2, pl1, pl2):
        discard hs.work.chainShift[nu] ^*! hs.work.outputChain[mu]
        discard scratch.shiftedOutputChain[mu][nu] ^*! dOutputChain[mu]
        threadBarrier()
        symStapleHVP(dLevel2Chain[nu,mu], dLevel2Chain[mu,nu],
                     scratch.acc1[mu,nu], scratch.acc2[mu,nu],
                     pl1, pl2, hs.state.stagePairShift[nu][mu], hs.state.stagePairShift[mu][nu],
                     hs.work.outputChain[mu], hs.work.chainShift[nu],
                     dl2[nu,mu], dl2[mu,nu], dOutputChain[mu],
                     dsl2[nu][mu], dsl2[mu][nu], scratch.shiftedOutputChain[mu][nu],
                     hs.work.tm1, hs.work.tm2, hs.state.sm1[nu], hs.state.sm1[mu])
        threadBarrier()

proc smearHVP*[G](hs: var HypSmear, dderiv: G, chain: G, dg: G,
                  dchain: G = @[]) =
  ## Second derivative of the uncontracted HYP reverse pass.
  ##
  ## Computes: dderiv = d/dε|_{ε=0} [ smearVJP(exp(ε dg) g, chain + ε dchain) ]
  ##
  ## When dchain is provided, this includes the variation of the chain field,
  ## enabling full chain rule computation for T(U) = <p, smearVJP(U, gaugeActionDeriv(HYP(U)))>.
  ##
  ## This is the pure forward-over-reverse Hessian–vector product for the HYP
  ## smear. Any outer contraction, such as contractProjectTAH(g, ·), is handled
  ## by the caller instead of being mixed into this API.
  static:
    doAssert(type(dderiv[0]) is typeof(hs.state.gf[0]))
    doAssert(type(chain[0])  is typeof(hs.state.gf[0]))
    doAssert(type(dg[0])     is typeof(hs.state.gf[0]))

  let sc = hs.scaleFactors
  var l1: type(hs.state.l1x)
  var l2: type(hs.state.l2x)
  initProjectedLevels(hs, l1, l2)
  let hsp = hs.addr

  ## Phase 1: analytical forward tangents
  let tangents = hs.buildTangentState(dg)
  let dgEff = tangents.dgEff
  let dl1x = tangents.dl1x
  let dl1 = tangents.dl1
  let dl2x = tangents.dl2x
  let dl2 = tangents.dl2
  let dflx = tangents.dflx
  let dsl2 = tangents.dsl2

  var work = newHvpWorkspace(hs.state.gf)
  clearHvpWorkspace(dderiv, work.dOutputChain, work.dLevel1Chain, work.dLevel2Chain)

  ## Phase 2: shared reverse-chain state
  buildReverseState(hsp, chain, nil, false)

  ## Phase 3: second-order adjoints through the three HYP levels
  applyOutputProjectionHVP(work.dOutputChain, hsp.state.flx, chain, dflx, dchain)
  accumulateAndScaleGauge(dderiv, work.dOutputChain, sc.ma3, sc.alp3)
  applyLevel3StapleHvpPhase(hsp, work.dLevel2Chain, l2, dl2, dsl2, work.dOutputChain)

  # Level 2: projection on l2x and staples built from l1
  applyPairProjectionHvpPhase(work.dLevel2Chain, work.dPreProjLevel2Chain, l2, hsp.state.l2x,
                              hsp.work.preProjLevel2Chain, dl2x, dderiv, sc.ma2, sc.alp2)
  applyPairStapleHvpPhase(2, hsp, l1, dl1, hsp.work.level2Chain, work.dLevel2Chain,
                          work.dLevel1Chain, hsp.work.tm1, hsp.work.tm2)

  # Level 1: projection on l1x and staples built directly from gf
  applyPairProjectionHvpPhase(work.dLevel1Chain, work.dPreProjLevel1Chain, l1, hsp.state.l1x,
                              hsp.work.preProjLevel1Chain, dl1x, dderiv, sc.ma1, sc.alp1)
  applyPairStapleHvpPhase(1, hsp, dgEff, dgEff, hsp.work.level1Chain, work.dLevel1Chain,
                          dderiv, hsp.work.tm1, hsp.work.tm2)

proc smearVJPChain*[G](hs: var HypSmear, chainbar: G, derivbar: G) =
  ## Adjoint of `smearVJP` with respect to `chain`.
  ## Computes `chainbar` such that:
  ##   ⟨derivbar, smearVJP(g, chain)⟩ = ⟨chainbar, chain⟩
  ## using the cached forward HYP state and the reverse-chain phases in reverse order.
  static:
    doAssert(type(chainbar[0]) is typeof(hs.state.gf[0]))
    doAssert(type(derivbar[0]) is typeof(hs.state.gf[0]))

  let sc = hs.scaleFactors
  var l1: type(hs.state.l1x)
  var l2: type(hs.state.l2x)
  initProjectedLevels(hs, l1, l2)
  let hsp = hs.addr

  # Phase 1: seed the adjoint workspaces
  var work = newChainAdjointWorkspace(chainbar, hsp.state.gf)

  # Phase 2: backpropagate the Level 1 reverse pair stage
  applyReversePairStageChainAdjoint(1, hsp, l1, work.level1ChainBar, derivbar, hsp.work.tm1, hsp.work.tm2)

  # Phase 3: apply the Level 1 direct term and projection-chain adjoint
  threads:
    forPairs(mu, nu):
      work.level1ChainBar[mu,nu] *= sc.alp1
      work.level1ChainBar[mu,nu] += sc.ma1 * derivbar[mu]
  applyPairProjectionChainAdjoint(work.level1ChainBar, l1, hsp.state.l1x)

  # Phase 4: backpropagate the Level 2 reverse pair stage
  applyReversePairStageChainAdjoint(2, hsp, l1, work.level2ChainBar, work.level1ChainBar, hsp.work.tm1, hsp.work.tm2)

  # Phase 5: apply the Level 2 direct term and projection-chain adjoint
  threads:
    forPairs(mu, nu):
      work.level2ChainBar[mu,nu] *= sc.alp2
      work.level2ChainBar[mu,nu] += sc.ma2 * derivbar[mu]
  applyPairProjectionChainAdjoint(work.level2ChainBar, l2, hsp.state.l2x)

  # Phase 6: backpropagate the Level 3 reverse pair stage
  forLevel3Pairs(hsp, l2, pl1, pl2):
    symStapleVJPChain(work.outputChainBar[mu],
                      pl1, pl2, hsp.state.stagePairShift[nu][mu], hsp.state.stagePairShift[mu][nu],
                      work.level2ChainBar[nu,mu], work.level2ChainBar[mu,nu],
                      hsp.work.chainShift[nu], hsp.work.seedShift[mu],
                      hsp.work.tm1, hsp.work.tm2, hsp.state.sm1[nu])
    threadBarrier()

  # Phase 7: apply the Level 3 direct term and finish at the output projection
  threads:
    for mu in 0..<4:
      work.outputChainBar[mu] *= sc.alp3
      work.outputChainBar[mu] += sc.ma3 * derivbar[mu]
  threads:
    for mu in 0..<4:
      chainbar[mu].projVJPChain(hsp.state.flx[mu], work.outputChainBar[mu])

template backpropLevel3HvpOutputPhase(hs, dl2bar, dOutputChainBar,
                                      dLevel2ChainBarFromL3, l2, alp: untyped) =
  block:
    var scratch = newLevel3AdjointScratch(hs.state.gf)

    forLevel3Pairs(hs, l2, pl1, pl2):
      discard hs.work.chainShift[nu] ^*! hs.work.outputChain[mu]
      threadBarrier()
      symStapleHVPVJP(dl2bar[nu,mu], dl2bar[mu,nu], scratch.dcbar[mu],
                      pl1, pl2, hs.state.stagePairShift[nu][mu], hs.state.stagePairShift[mu][nu],
                      hs.work.outputChain[mu], hs.work.chainShift[nu],
                      dLevel2ChainBarFromL3[nu,mu], dLevel2ChainBarFromL3[mu,nu],
                      hs.work.seedShift[nu], hs.work.seedShift[mu],
                      hs.work.tm1, hs.work.tm2, hs.state.sm1[nu], hs.state.sm1[mu])
      threadBarrier()

    addScaledGauge(dOutputChainBar, scratch.dcbar, alp)

template backpropLevel3TangentPhase(hs, dgEffbar, dl2bar, dOutputTangentBar,
                                    l2, ma, alp: untyped) =
  block:
    threads:
      for mu in 0..<4:
        dgEffbar[mu] += ma * dOutputTangentBar[mu]

    forLevel3Pairs(hs, l2, pl1, pl2):
      symStapleJVPVJP(dl2bar[nu,mu], dl2bar[mu,nu], alp,
                      pl1, pl2, hs.state.stagePairShift[nu][mu], hs.state.stagePairShift[mu][nu],
                      dOutputTangentBar[mu], hs.work.seedShift[nu],
                      hs.work.tm1, hs.work.tm2, hs.state.sm1[nu], hs.state.sm1[mu])
      threadBarrier()

proc smearHVPVJP*[G](hs: var HypSmear, gbar: G, derivbar: G, chain: G) =
  ## Adjoint of the `dg` branch of `smearHVP`.
  ## Computes `gbar` such that:
  ##   <derivbar, smearHVP(dg, dchain=0)> = <gbar, dg*g>
  ## with the upstream chain held fixed. This reverses the forward-over-reverse
  ## HVP path, including the dependence of the cached HYP state on the gauge field.
  static:
    doAssert(type(gbar[0]) is typeof(hs.state.gf[0]))
    doAssert(type(derivbar[0]) is typeof(hs.state.gf[0]))
    doAssert(type(chain[0]) is typeof(hs.state.gf[0]))

  let sc = hs.scaleFactors
  var l1: type(hs.state.l1x)
  var l2: type(hs.state.l2x)
  initProjectedLevels(hs, l1, l2)
  let hsp = hs.addr
  buildReverseState(hsp, chain, nil, false)

  var work = newHvpAdjointWorkspace(gbar, hsp.state.gf)

  # Phase 1: seed the adjoints from the output contribution of `smearHVP`
  seedSmearHVPAdjoints(work.dfl1bar, work.dfl2bar, work.dfcbar, derivbar, sc.ma1, sc.ma2, sc.ma3)

  # Phase 2: reverse the Level 1 and Level 2 HVP output paths
  backpropPairHvpOutputPhase(1, hsp, l1, hsp.work.level1Chain, derivbar, work.dgEffbar,
                             work.dfl1bar, work.dfl1barL2,
                             l1, hsp.state.l1x, hsp.work.preProjLevel1Chain, sc.alp1,
                             hsp.work.tm1, hsp.work.tm2)

  backpropPairHvpOutputPhase(2, hsp, l1, hsp.work.level2Chain, work.dfl1barL2,
                             work.dl1bar, work.dfl2bar, work.dfl2barL3,
                             l2, hsp.state.l2x, hsp.work.preProjLevel2Chain, sc.alp2,
                             hsp.work.tm1, hsp.work.tm2)

  # Phase 3: reverse the Level 3 HVP output path back into the tangent state
  backpropLevel3HvpOutputPhase(hsp, work.dl2bar, work.dfcbar, work.dfl2barL3, l2, sc.alp3)
  threads:
    for mu in 0..<4:
      work.dflxbar[mu].projHVPVJP_dx(hsp.state.flx[mu], chain[mu], work.dfcbar[mu])
  backpropLevel3TangentPhase(hsp, work.dgEffbar, work.dl2bar, work.dflxbar, l2, sc.ma3, sc.alp3)

  # Phase 4: reverse the Level 2 projection and tangent contributions
  accumulateProjectionAdjoint(work.dl2xbar, hsp.state.l2x, hsp.work.preProjLevel2Chain, work.dl2bar, work.dfl2bar)
  backpropPairTangentPhase(2, hsp, l1, work.dl2xbar, work.dl1bar, work.dgEffbar, sc.ma2, sc.alp2,
                           hsp.work.tm1, hsp.work.tm2)

  # Phase 5: reverse the Level 1 projection and tangent contributions
  accumulateProjectionAdjoint(work.dl1xbar, hsp.state.l1x, hsp.work.preProjLevel1Chain, work.dl1bar, work.dfl1bar)
  backpropPairTangentPhase(1, hsp, hsp.state.gf, work.dl1xbar, work.dl1xbar, work.dgEffbar, sc.ma1, sc.alp1,
                           hsp.work.tm1, hsp.work.tm2)

  # Phase 6: return the adjoint with respect to the multiplicative tangent `dg*g`
  threads:
    for mu in 0..<4:
      gbar[mu] := work.dgEffbar[mu]
