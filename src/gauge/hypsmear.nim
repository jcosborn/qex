import base, layout, gauge, fat7l, smearutil
export PerfInfo

const keepProj {.boolDefine.} = true
when keepProj:
  static: echo "hypsmear: keeping projected fields"
else:
  static: echo "hypsmear: NOT keeping projected fields"

const
  psL1 = true
  psL2 = false

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
  tic()
  type lcm = type(gf[0])
  let lo = gf[0].l
  proc newlcm: lcm = result.new(gf[0].l)
  var
    l1x = newFieldArray2(lo,lcm,[4,4],mu!=nu)
    l2x = newOneOf(l1x)
    flx = newFieldArray(lo,lcm,4)
    tm1: lcm
    sm1: array[4,Shifter[lcm,type(gf[0][0])]]
    s1: array[4,array[4,Shifter[lcm,type(gf[0][0])]]]
    #nflop = 61632.0
    #dtime = 0.0
  when keepProj:
    var
      l1 = newOneOf(l1x)
      l2 = newOneOf(l1x)
  else:
    var
      lp1 = newlcm()
      lp2 = newlcm()

  tm1 = newlcm()
  for mu in 0..<4:
    sm1[mu] = newShifter(gf[mu], mu, -1)
    for nu in 0..<4:
      if nu!=mu:
        s1[mu][nu] = newShifter(gf[mu], nu, 1)
  threads:
    for mu in 0..<4:
      for nu in 0..<4:
        if nu!=mu:
          discard s1[mu][nu] ^*! gf[mu]

  let
    alp1 = coef.alpha1 / 2.0
    alp2 = coef.alpha2 / 4.0
    alp3 = coef.alpha3 / 6.0
    ma1 = 1 - coef.alpha1
    ma2 = 1 - coef.alpha2
    ma3 = 1 - coef.alpha3

  toc("prep")
  threads:
    for mu in 0..<4:
      for nu in 0..<4:
        if nu!=mu:
          l1x[mu,nu] := ma1 * gf[mu]
          symStaple(l1x[mu,nu], alp1, gf[nu], gf[mu],
                    s1[nu][mu], s1[mu][nu], tm1, sm1[nu])
          when keepProj:
            l1[mu,nu].proj l1x[mu,nu]
    toc("1")

    for mu in 0..<4:
      for nu in 0..<4:
        if nu!=mu:
          l2x[mu,nu] := ma2 * gf[mu]
          for a in 0..<4:
            if a!=mu and a!=nu:
              let b = 1+2+3-mu-nu-a
              when keepProj:
                template lp1:untyped = l1[a,b]
                template lp2:untyped = l1[mu,b]
              else:
                lp1.proj l1x[a,b]
                lp2.proj l1x[mu,b]
              threadBarrier()
              discard s1[nu][mu] ^*! lp1
              discard s1[mu][a] ^*! lp2
              threadBarrier()
              symStaple(l2x[mu,nu], alp2, lp1, lp2,
                        s1[nu][mu], s1[mu][a], tm1, sm1[a])
          when keepProj:
            l2[mu,nu].proj l2x[mu,nu]
    toc("2")

    when keepProj:
      for mu in 0..<4:
        for nu in 0..<4:
          if nu!=mu:
            threadBarrier()
            discard s1[mu][nu] ^*! l2[mu,nu]

    for mu in 0..<4:
      flx[mu] := ma3 * gf[mu]
      for nu in 0..<4:
        if nu!=mu:
          when keepProj:
            template lp1:untyped = l2[nu,mu]
            template lp2:untyped = l2[mu,nu]
          else:
            lp1.proj l2x[nu,mu]
            lp2.proj l2x[mu,nu]
            threadBarrier()
            discard s1[nu][mu] ^*! lp1
            discard s1[mu][nu] ^*! lp2
            threadBarrier()
          symStaple(flx[mu], alp3, lp1, lp2,
                    s1[nu][mu], s1[mu][nu], tm1, sm1[nu])
      fl[mu].proj flx[mu]
  toc("threads end")

  proc smearedForce(f,chain:G) =
    tic("smearedF")
    # fₓₚₜ ← chainₘₖₕ d/dUₓₚₜ^*[Vₘₖₕ(U)^*] + chainₘₕₖ^* d/dUₓₚₜ^*[Vₘₕₖ(U)]
    var
      fl1 = newFieldArray2(lo,lcm,[4,4],mu!=nu)
      fl2 = newOneOf(fl1)
      fc = newFieldArray(lo,lcm,4)
      fs: array[4,Shifter[lcm,type(gf[0][0])]]
      tm2: lcm
    tm2 = newlcm()
    for mu in 0..<4:
      fs[mu] = newShifter(fc[mu], mu, 1)
    toc("prep")

    threads:
      for mu in 0..<4:
        for nu in 0..<4:
          if nu!=mu:
            fl1[mu,nu] := 0
            fl2[mu,nu] := 0

      # proj flx → fl, fc ← chain
      for mu in 0..<4:
        fc[mu].projDeriv(flx[mu], chain[mu])
      # link (gf, l2) → flx, (f, fl2) ← fc
      for mu in 0..<4:
        f[mu] := ma3 * fc[mu]
        fc[mu] *= alp3
      when keepProj:
        for mu in 0..<4:
          for nu in 0..<4:
            if nu!=mu:
              threadBarrier()
              discard s1[mu][nu] ^*! l2[mu,nu]
      for mu in 0..<4:
        for nu in 0..<4:
          if nu!=mu:
            when keepProj:
              template lp1:untyped = l2[nu,mu]
              template lp2:untyped = l2[mu,nu]
            else:
              lp1.proj l2x[nu,mu]
              lp2.proj l2x[mu,nu]
              threadBarrier()
              discard s1[nu][mu] ^*! lp1
              discard s1[mu][nu] ^*! lp2
            discard fs[nu] ^*! fc[mu]
            threadBarrier()
            symStapleDeriv(fl2[nu,mu], fl2[mu,nu],
                           lp1, lp2, s1[nu][mu], s1[mu][nu],
                           fc[mu], fs[nu], tm1, tm2, sm1[nu], sm1[mu])
      toc("1")

      # proj l2x → l2, fl2 ← fl2
      for mu in 0..<4:
        for nu in 0..<4:
          if nu!=mu:
            when keepProj:
              fl2[mu,nu].projDeriv(l2[mu,nu], l2x[mu,nu], fl2[mu,nu])
            else:
              fl2[mu,nu].projDeriv(l2x[mu,nu], fl2[mu,nu])
      # link (gf, l1) → l2x, (f, fl1) ← fl2
      for mu in 0..<4:
        for nu in 0..<4:
          if nu!=mu:
            f[mu] += ma2 * fl2[mu,nu]
            fl2[mu,nu] *= alp2
      for mu in 0..<4:
        for nu in 0..<4:
          if nu!=mu:
            for a in 0..<4:
              if a!=mu and a!=nu:
                let b = 1+2+3-mu-nu-a
                when keepProj:
                  template lp1:untyped = l1[a,b]
                  template lp2:untyped = l1[mu,b]
                else:
                  lp1.proj l1x[a,b]
                  lp2.proj l1x[mu,b]
                threadBarrier()
                discard s1[nu][mu] ^*! lp1
                discard s1[mu][a] ^*! lp2
                discard fs[a] ^*! fl2[mu,nu]
                threadBarrier()
                symStapleDeriv(fl1[a,b], fl1[mu,b],
                               lp1, lp2, s1[nu][mu], s1[mu][a],
                               fl2[mu,nu], fs[a], tm1, tm2, sm1[a], sm1[mu])
      toc("2")

      # proj l1x → l1, fl1 ← fl1
      for mu in 0..<4:
        for nu in 0..<4:
          if nu!=mu:
            when keepProj:
              fl1[mu,nu].projDeriv(l1[mu,nu], l1x[mu,nu], fl1[mu,nu])
            else:
              fl1[mu,nu].projDeriv(l1x[mu,nu], fl1[mu,nu])
            discard s1[mu][nu] ^*! gf[mu]
      # link gf → l1, f ← fl1
      for mu in 0..<4:
        for nu in 0..<4:
          if nu!=mu:
            f[mu] += ma1 * fl1[mu,nu]
            fl1[mu,nu] *= alp1
      for mu in 0..<4:
        for nu in 0..<4:
          if nu!=mu:
            discard fs[nu] ^* fl1[mu,nu]
            symStapleDeriv(f[nu], f[mu],
                           gf[nu], gf[mu], s1[nu][mu], s1[mu][nu],
                           fl1[mu,nu], fs[nu], tm1, tm2, sm1[nu], sm1[mu])
    toc("end")

  toc("end")
  smearedForce

proc smearPriv[G](coef: HypCoefs, gf: G, fl: G, info: var PerfInfo) {.codegenDecl:
    "__attribute__((noinline)) $# $#$#".} =
  # Avoid inlining or other compiler optimizations
  # in order to guarantee the change of the stack pointer,
  # such that Nim's GC is able to collect the memory.
  {.emit: "asm (\"\");".}
  var f = coef.smearGetForce(gf, fl, info)
  f = nil
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
  HypForwardState[V:static[int],F,T] = object
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
    tm1*: F
    l1x*, l2x*: FieldArray[2,V,T]
    when keepProj:
      l1*, l2*: FieldArray[2,V,T]
    flx*: FieldArray[1,V,T]

  ## When reused `FieldArray` scratch from these bundles is touched inside
  ## `threads`, access it through `.addr`/`p[]` aliases. Do not capture a
  ## local `FieldArray` value directly in a `threads` block; under `refc`
  ## that pattern caused runtime crashes in derivative-path code.
  HypReverseStateScratch[V:static[int],F,T] = ref object
    outputChain*: FieldArray[1,V,T]
    level1Chain*, level2Chain*: FieldArray[2,V,T]
    preProjLevel1Chain*, preProjLevel2Chain*: FieldArray[2,V,T]

  HypKernelScratch[F,T] = ref object
    chainShift*: array[4,Shifter[F,T]]
    seedShift*: array[4,Shifter[F,T]]
    ## Scratch scalar/matrix temp paired with the reverse/tangent kernel layer.
    tm2*: F

  HypChainAdjointScratch[V:static[int],F,T] = ref object
    pairProjectionTmp*: FieldArray[2,V,T]
    outputChainBar*: seq[F]
    level1ChainBar*, level2ChainBar*: FieldArray[2,V,T]

  HypTangentScratch[V:static[int],F,T] = ref object
    tangentDgEff*: seq[F]
    tangentDflx*: FieldArray[1,V,T]
    tangentDl1x*, tangentDl2x*, tangentDl1*, tangentDl2*: FieldArray[2,V,T]
    tangentDs1*, tangentDsl2*: array[4,array[4,Shifter[F,T]]]
    tangentDsl1L2*: array[4, array[4, array[4, array[2, Shifter[F,T]]]]]

  HypProjectedTmpScratch[V:static[int],F,T] = ref object
    projectedLevel1Tmp*, projectedLevel2Tmp*: FieldArray[2,V,T]

  HypHvpScratch[V:static[int],F,T] = ref object
    dOutputChain*: seq[F]
    dLevel1Chain*, dLevel2Chain*: FieldArray[2,V,T]
    dPreProjLevel1Chain*, dPreProjLevel2Chain*: FieldArray[2,V,T]
    hvpL1Acc1*, hvpL1Acc2*: FieldArray[2,V,T]
    hvpL2Acc1*, hvpL2Acc2*: array[4, array[4, array[4, F]]]
    hvpL1ShiftedPairTangents*: array[4, array[4, array[2, Shifter[F,T]]]]
    hvpL1ShiftedChainTangents*: array[4, array[4, Shifter[F,T]]]
    hvpL2ShiftedPairTangents*: array[4, array[4, array[4, array[2, Shifter[F,T]]]]]
    hvpL2ShiftedChainTangents*: array[4, array[4, array[4, Shifter[F,T]]]]
    level3Hvp*: HypLevel3HvpScratch[FieldArray[2,V,T], array[4, array[4, Shifter[F,T]]]]

  HypHvpAdjointScratch[V:static[int],F,T] = ref object
    dflxbar*, dfcbar*, dgEffbar*: seq[F]
    dfl1bar*, dfl2bar*: FieldArray[2,V,T]
    dl1bar*, dl2bar*: FieldArray[2,V,T]
    dl1xbar*, dl2xbar*: FieldArray[2,V,T]
    dfl1barL2*, dfl2barL3*: FieldArray[2,V,T]
    hvpVjpDcbar*: FieldArray[2,V,T]
    level3Adjoint*: HypLevel3AdjointScratch[seq[F]]

  HypSmear*[V:static[int],F,T] = object
    ## Reusable buffers for HYP smearing and its derivative.
    ##
    ## The forward state is intentionally separated from the lazy scratch
    ## bundles so plain `smear` keeps a small footprint. Each derivative family
    ## allocates its own bundle on first use and then reuses it thereafter.
    coef*: HypCoefs
    state*: HypForwardState[V,F,T]
    reverse*: HypReverseStateScratch[V,F,T]
    kernel*: HypKernelScratch[F,T]
    chainAdjoint*: HypChainAdjointScratch[V,F,T]
    tangent*: HypTangentScratch[V,F,T]
    projectedTmp*: HypProjectedTmpScratch[V,F,T]
    hvp*: HypHvpScratch[V,F,T]
    hvpAdjoint*: HypHvpAdjointScratch[V,F,T]

  ProjectedLevels[PairField] = object
    l1*, l2*: PairField

  HypLevel3HvpScratch[PairField, ShiftField] = object
    shiftedOutputChain*: ShiftField
    acc1*, acc2*: PairField

  HypLevel3AdjointScratch[GaugeField] = object
    dcbar*: GaugeField

proc scaleFactors[V:static[int],F,T](hs: HypSmear[V,F,T]): HypScaleFactors {.inline.} =
  result = HypScaleFactors(
    alp1: hs.coef.alpha1 / 2.0,
    alp2: hs.coef.alpha2 / 4.0,
    alp3: hs.coef.alpha3 / 6.0,
    ma1: 1 - hs.coef.alpha1,
    ma2: 1 - hs.coef.alpha2,
    ma3: 1 - hs.coef.alpha3,
  )

template pairStageScales(hs, stage: untyped): untyped =
  block:
    let sc = hs.scaleFactors
    when stage == psL1:
      (ma: sc.ma1, alp: sc.alp1)
    else:
      (ma: sc.ma2, alp: sc.alp2)

proc outputStageScales(hs: HypSmear): auto {.inline.} =
  let sc = hs.scaleFactors
  result = (ma: sc.ma3, alp: sc.alp3)

template allocShiftedTangents(s, f: untyped) =
  forPairs(mu, nu):
    s[mu][nu] = newShifter(f[mu], nu, 1)

template allocThirdShiftedTangents(s, f: untyped) =
  forPairs(mu, nu):
    forThird(mu, nu, a, b):
      s[mu][nu][a][0] = newShifter(f[a,b], mu, 1)
      s[mu][nu][a][1] = newShifter(f[mu,b], a, 1)

template primeShiftedTangentsPtr(sp, fp: untyped) =
  threads:
    forPairs(mu, nu):
      discard sp[][mu][nu] ^*! fp[][mu]

proc ensureReverseStateScratch[V:static[int],F,T](hs: var HypSmear[V,F,T]) =
  if hs.reverse != nil:
    return
  new hs.reverse
  let lo = hs.state.flx[0].l
  type GaugeField = typeof(hs.state.flx[0])
  hs.reverse.outputChain = newFieldArray(lo, GaugeField, 4)
  hs.reverse.level1Chain = newFieldArray2(lo, GaugeField, [4,4], mu!=nu)
  hs.reverse.level2Chain = newFieldArray2(lo, GaugeField, [4,4], mu!=nu)
  hs.reverse.preProjLevel1Chain = newFieldArray2(lo, GaugeField, [4,4], mu!=nu)
  hs.reverse.preProjLevel2Chain = newFieldArray2(lo, GaugeField, [4,4], mu!=nu)

proc ensureKernelScratch[V:static[int],F,T](hs: var HypSmear[V,F,T]) =
  if hs.kernel != nil:
    return
  new hs.kernel
  hs.kernel.tm2.new(hs.state.flx[0].l)
  for mu in 0..<4:
    hs.kernel.chainShift[mu] = newShifter(hs.state.flx[mu], mu, 1)
    hs.kernel.seedShift[mu] = newShifter(hs.state.flx[mu], mu, 1)

proc ensureChainAdjointScratch[V:static[int],F,T](hs: var HypSmear[V,F,T]) =
  if hs.chainAdjoint != nil:
    return
  new hs.chainAdjoint
  let lo = hs.state.flx[0].l
  type GaugeField = typeof(hs.state.flx[0])
  hs.chainAdjoint.pairProjectionTmp = newFieldArray2(lo, GaugeField, [4,4], mu!=nu)
  hs.chainAdjoint.outputChainBar = lo.newGauge
  hs.chainAdjoint.level1ChainBar = newFieldArray2(lo, GaugeField, [4,4], mu!=nu)
  hs.chainAdjoint.level2ChainBar = newFieldArray2(lo, GaugeField, [4,4], mu!=nu)

proc ensureTangentScratch[V:static[int],F,T](hs: var HypSmear[V,F,T]) =
  if hs.tangent != nil:
    return
  new hs.tangent
  let lo = hs.state.flx[0].l
  type GaugeField = typeof(hs.state.flx[0])
  hs.tangent.tangentDgEff = lo.newGauge
  hs.tangent.tangentDflx = newFieldArray(lo, GaugeField, 4)
  hs.tangent.tangentDl1x = newFieldArray2(lo, GaugeField, [4,4], mu!=nu)
  hs.tangent.tangentDl2x = newFieldArray2(lo, GaugeField, [4,4], mu!=nu)
  hs.tangent.tangentDl1 = newFieldArray2(lo, GaugeField, [4,4], mu!=nu)
  hs.tangent.tangentDl2 = newFieldArray2(lo, GaugeField, [4,4], mu!=nu)
  allocShiftedTangents(hs.tangent.tangentDs1, hs.tangent.tangentDgEff)
  allocThirdShiftedTangents(hs.tangent.tangentDsl1L2, hs.tangent.tangentDl1)
  for mu in 0..<4:
    for nu in 0..<4:
      if nu != mu:
        hs.tangent.tangentDsl2[mu][nu] = newShifter(hs.tangent.tangentDl2[mu,nu], nu, 1)

proc ensureProjectedTmpScratch[V:static[int],F,T](hs: var HypSmear[V,F,T]) =
  if hs.projectedTmp != nil:
    return
  new hs.projectedTmp
  let lo = hs.state.flx[0].l
  type GaugeField = typeof(hs.state.flx[0])
  hs.projectedTmp.projectedLevel1Tmp = newFieldArray2(lo, GaugeField, [4,4], mu!=nu)
  hs.projectedTmp.projectedLevel2Tmp = newFieldArray2(lo, GaugeField, [4,4], mu!=nu)

proc ensureHvpScratch[V:static[int],F,T](hs: var HypSmear[V,F,T]) =
  if hs.hvp != nil:
    return
  new hs.hvp
  let lo = hs.state.flx[0].l
  type GaugeField = typeof(hs.state.flx[0])
  hs.hvp.dOutputChain = lo.newGauge
  hs.hvp.dLevel1Chain = newFieldArray2(lo, GaugeField, [4,4], mu!=nu)
  hs.hvp.dLevel2Chain = newFieldArray2(lo, GaugeField, [4,4], mu!=nu)
  hs.hvp.dPreProjLevel1Chain = newFieldArray2(lo, GaugeField, [4,4], mu!=nu)
  hs.hvp.dPreProjLevel2Chain = newFieldArray2(lo, GaugeField, [4,4], mu!=nu)
  hs.hvp.level3Hvp.acc1 = newFieldArray2(lo, GaugeField, [4,4], mu!=nu)
  hs.hvp.level3Hvp.acc2 = newFieldArray2(lo, GaugeField, [4,4], mu!=nu)
  hs.hvp.hvpL1Acc1 = newFieldArray2(lo, GaugeField, [4,4], mu!=nu)
  hs.hvp.hvpL1Acc2 = newFieldArray2(lo, GaugeField, [4,4], mu!=nu)
  forPairs(mu, nu):
    hs.hvp.level3Hvp.shiftedOutputChain[mu][nu] = newShifter(hs.hvp.dOutputChain[mu], nu, 1)
    hs.hvp.hvpL1ShiftedPairTangents[mu][nu][0] = newShifter(hs.state.flx[nu], mu, 1)
    hs.hvp.hvpL1ShiftedPairTangents[mu][nu][1] = newShifter(hs.state.flx[mu], nu, 1)
    hs.hvp.hvpL1ShiftedChainTangents[mu][nu] = newShifter(hs.hvp.dLevel1Chain[mu,nu], nu, 1)
    forThird(mu, nu, a, b):
      hs.hvp.hvpL2Acc1[mu][nu][a] = lo.newField(T)
      hs.hvp.hvpL2Acc2[mu][nu][a] = lo.newField(T)
      hs.hvp.hvpL2ShiftedPairTangents[mu][nu][a][0] = newShifter(hs.state.l1x[a,b], mu, 1)
      hs.hvp.hvpL2ShiftedPairTangents[mu][nu][a][1] = newShifter(hs.state.l1x[mu,b], a, 1)
      hs.hvp.hvpL2ShiftedChainTangents[mu][nu][a] = newShifter(hs.hvp.dLevel2Chain[mu,nu], a, 1)

proc ensureHvpAdjointScratch[V:static[int],F,T](hs: var HypSmear[V,F,T]) =
  if hs.hvpAdjoint != nil:
    return
  new hs.hvpAdjoint
  let lo = hs.state.flx[0].l
  type GaugeField = typeof(hs.state.flx[0])
  hs.hvpAdjoint.dflxbar = lo.newGauge
  hs.hvpAdjoint.dfcbar = lo.newGauge
  hs.hvpAdjoint.dgEffbar = lo.newGauge
  hs.hvpAdjoint.dfl1bar = newFieldArray2(lo, GaugeField, [4,4], mu!=nu)
  hs.hvpAdjoint.dfl2bar = newFieldArray2(lo, GaugeField, [4,4], mu!=nu)
  hs.hvpAdjoint.dl1bar = newFieldArray2(lo, GaugeField, [4,4], mu!=nu)
  hs.hvpAdjoint.dl2bar = newFieldArray2(lo, GaugeField, [4,4], mu!=nu)
  hs.hvpAdjoint.dl1xbar = newFieldArray2(lo, GaugeField, [4,4], mu!=nu)
  hs.hvpAdjoint.dl2xbar = newFieldArray2(lo, GaugeField, [4,4], mu!=nu)
  hs.hvpAdjoint.dfl1barL2 = newFieldArray2(lo, GaugeField, [4,4], mu!=nu)
  hs.hvpAdjoint.dfl2barL3 = newFieldArray2(lo, GaugeField, [4,4], mu!=nu)
  hs.hvpAdjoint.hvpVjpDcbar = newFieldArray2(lo, GaugeField, [4,4], mu!=nu)
  hs.hvpAdjoint.level3Adjoint.dcbar = lo.newGauge

proc newHypSmear*(l: Layout, coef: HypCoefs): auto =
  ## Allocate a reusable environment for HYP smearing and derivatives
  type F = type(l.newGauge[0])
  type T = type(l.newGauge[0][0])
  type L = type(l)
  var hs: HypSmear[L.V,F,T]
  hs.coef = coef
  hs.state.l1x = newFieldArray2(l,F,[4,4],mu!=nu)
  hs.state.l2x = newFieldArray2(l,F,[4,4],mu!=nu)
  when keepProj:
    hs.state.l1 = newFieldArray2(l,F,[4,4],mu!=nu)
    hs.state.l2 = newFieldArray2(l,F,[4,4],mu!=nu)
  hs.state.flx = newFieldArray(l,F,4)
  hs.state.tm1.new(l)
  for mu in 0..<4:
    hs.state.sm1[mu] = newShifter(hs.state.flx[mu], mu, -1)
    for nu in 0..<4:
      if nu != mu:
        hs.state.gaugePairShift[mu][nu] = newShifter(hs.state.flx[mu], nu, 1)
        hs.state.stagePairShift[mu][nu] = newShifter(hs.state.flx[mu], nu, 1)
  hs

proc initProjectedLevelsImpl(hs: var HypSmear; l1, l2: var auto;
                             populateFromCache: static bool) =
  when keepProj:
    l1 = hs.state.l1
    l2 = hs.state.l2
  else:
    hs.ensureProjectedTmpScratch()
    l1 = hs.projectedTmp.projectedLevel1Tmp
    l2 = hs.projectedTmp.projectedLevel2Tmp
    when populateFromCache:
      let hsp = hs.addr
      let l1p = l1.addr
      let l2p = l2.addr
      threads:
        forPairs(mu, nu):
          l1p[][mu,nu].proj hsp.state.l1x[mu,nu]
          l2p[][mu,nu].proj hsp.state.l2x[mu,nu]

proc initProjectedLevels(hs: var HypSmear; l1, l2: var auto) =
  initProjectedLevelsImpl(hs, l1, l2, true)

proc projectedLevels(hs: var HypSmear): auto =
  var levels: ProjectedLevels[type(hs.state.l1x)]
  initProjectedLevels(hs, levels.l1, levels.l2)
  levels

proc buildGaugeTangent(dst, gf, dg: auto) =
  let dstp = dst.addr
  threads:
    for mu in 0..<gf.len:
      dstp[][mu] := dg[mu] * gf[mu]

proc initForwardProjectedLevels(hs: var HypSmear; l1, l2: var auto) =
  ## Allocate the projected Level-1/Level-2 buffers used while rebuilding the
  ## forward HYP state. Unlike `initProjectedLevels`, this does not populate the
  ## buffers from cached state.
  initProjectedLevelsImpl(hs, l1, l2, false)

proc forwardProjectedLevels(hs: var HypSmear): auto =
  var levels: ProjectedLevels[type(hs.state.l1x)]
  initForwardProjectedLevels(hs, levels.l1, levels.l2)
  levels

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
  let dflxp = dflx.addr
  let dl1xp = dl1x.addr
  let dl2xp = dl2x.addr
  let dl1p = dl1.addr
  let dl2p = dl2.addr
  threads:
    for mu in 0..<4:
      dflxp[][mu] := 0
      forNu(mu, nu):
        dl1xp[][mu,nu] := 0
        dl2xp[][mu,nu] := 0
        dl1p[][mu,nu] := 0
        dl2p[][mu,nu] := 0

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

template clearHvpAdjointIntermediates(dfl1barL2, dfl2barL3: auto) =
  threads:
    forPairs(mu, nu):
      dfl1barL2[mu,nu] := 0
      dfl2barL3[mu,nu] := 0

template clearLevel3AdjointScratch(dcbar: auto) =
  threads:
    for mu in 0..<4:
      dcbar[mu] := 0

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

proc hypOutputStageStapleVJP*(hs: var HypSmear, l2, outputChain: auto;
                              level2Chain: auto) =
  hs.ensureKernelScratch()
  let hsp = hs.addr
  forLevel3Pairs(hsp[], l2, pl1, pl2):
    threadBarrier()
    discard hsp[].kernel.chainShift[nu] ^*! outputChain[mu]
    threadBarrier()
    symStapleVJP(level2Chain[nu,mu], level2Chain[mu,nu],
                 pl1, pl2,
                 hsp[].state.stagePairShift[nu][mu], hsp[].state.stagePairShift[mu][nu],
                 outputChain[mu], hsp[].kernel.chainShift[nu],
                 hsp[].state.tm1, hsp[].kernel.tm2, hsp[].state.sm1[nu], hsp[].state.sm1[mu])
    threadBarrier()

proc hypOutputStageChainVJP*(hs: var HypSmear, l2, seedBar: auto;
                             outputChainBar: auto) =
  hs.ensureKernelScratch()
  let hsp = hs.addr
  forLevel3Pairs(hsp[], l2, pl1, pl2):
    symStapleVJPChain(outputChainBar[mu],
                      pl1, pl2, hsp[].state.stagePairShift[nu][mu], hsp[].state.stagePairShift[mu][nu],
                      seedBar[nu,mu], seedBar[mu,nu],
                      hsp[].kernel.chainShift[nu], hsp[].kernel.seedShift[mu],
                      hsp[].state.tm1, hsp[].kernel.tm2, hsp[].state.sm1[nu])
    threadBarrier()

template primePairStageForwardShifts(hs, stage: untyped) =
  when stage == psL1:
    threads:
      forPairs(mu, nu):
        discard hs.state.gaugePairShift[mu][nu] ^*! hs.state.gf[mu]

template forPairStageForwardTerms(hs, stage, pairSrc,
                                  mu, nu, pl1, pl2, sh1, sh2, smDir: untyped;
                                  body: untyped) =
  when stage == psL1:
    let pl1 = hs.state.gf[nu]
    let pl2 = hs.state.gf[mu]
    let sh1 = hs.state.gaugePairShift[nu][mu]
    let sh2 = hs.state.gaugePairShift[mu][nu]
    let smDir = nu
    body
  else:
    forThird(mu, nu, a, b):
      let pl1 = pairSrc[a,b]
      let pl2 = pairSrc[mu,b]
      let sh1 = hs.state.stagePairShift[nu][mu]
      let sh2 = hs.state.stagePairShift[mu][a]
      let smDir = a
      threadBarrier()
      discard sh1 ^*! pl1
      discard sh2 ^*! pl2
      threadBarrier()
      body
      threadBarrier()

template hypPairStageForward(stage, hs, pairSrc, outX, outProj: untyped) =
  block:
    let hsp = hs.addr
    let sc = pairStageScales(hsp[], stage)
    primePairStageForwardShifts(hsp[], stage)

    threads:
      forPairs(mu, nu):
        outX[mu,nu] := sc.ma * hsp[].state.gf[mu]
        forPairStageForwardTerms(hsp[], stage, pairSrc, mu, nu, pl1, pl2, sh1, sh2, smDir):
          symStaple(outX[mu,nu], sc.alp, pl1, pl2,
                    sh1, sh2, hsp[].state.tm1, hsp[].state.sm1[smDir])
        when stage == psL1:
          threadBarrier()
        outProj[mu,nu].proj outX[mu,nu]

## Pair-stage operator family split into explicit Level-1/Level-2 entry points.
proc hypPairStageL1*(hs: var HypSmear; outX, outProj: auto) =
  hypPairStageForward(psL1, hs, hs.state.l1x, outX, outProj)

proc hypPairStageL2*(hs: var HypSmear, pairSrc: auto;
                     outX, outProj: auto) =
  hypPairStageForward(psL2, hs, pairSrc, outX, outProj)

template hypPairStageJVPPtr(stage, hsp, pairSrc, tangentSrc, shiftedTangents,
                         preProjBase, projBase, dgEff, preProjOut, tangentOut: untyped) =
  block:
    let dgEffp = dgEff.addr
    let preProjOutp = preProjOut.addr
    let tangentOutp = tangentOut.addr
    let sc = pairStageScales(hsp[], stage)
    primePairStageForwardShifts(hsp[], stage)
    threads:
      forPairs(mu, nu):
        preProjOutp[][mu,nu] := sc.ma * dgEffp[][mu]
        when stage == psL1:
          symStapleJVP(preProjOutp[][mu,nu], sc.alp,
                       hsp[].state.gf[nu], hsp[].state.gf[mu],
                       hsp[].state.gaugePairShift[nu][mu], hsp[].state.gaugePairShift[mu][nu],
                       tangentSrc[nu], tangentSrc[mu],
                       shiftedTangents[nu][mu], shiftedTangents[mu][nu],
                       hsp[].state.tm1, hsp[].state.sm1[nu])
          threadBarrier()
        else:
          forThird(mu, nu, a, b):
            let pl1 = pairSrc[a,b]
            let pl2 = pairSrc[mu,b]
            let dpl1 = tangentSrc[a,b]
            let dpl2 = tangentSrc[mu,b]
            threadBarrier()
            discard hsp[].state.stagePairShift[nu][mu] ^*! pl1
            discard hsp[].state.stagePairShift[mu][a] ^*! pl2
            discard shiftedTangents[mu][nu][a][0] ^*! dpl1
            discard shiftedTangents[mu][nu][a][1] ^*! dpl2
            threadBarrier()
            symStapleJVP(preProjOutp[][mu,nu], sc.alp,
                         pl1, pl2, hsp[].state.stagePairShift[nu][mu], hsp[].state.stagePairShift[mu][a],
                         dpl1, dpl2,
                         shiftedTangents[mu][nu][a][0], shiftedTangents[mu][nu][a][1],
                         hsp[].state.tm1, hsp[].state.sm1[a])
            threadBarrier()
        tangentOutp[][mu,nu].projJVP(projBase[mu,nu], preProjBase[mu,nu], preProjOutp[][mu,nu])

template hypPairStageJVPL1Ptr(hsp, tangentSrc, shiftedTangents,
                              preProjBase, projBase, dgEff,
                              preProjOut, tangentOut: untyped) =
  hypPairStageJVPPtr(psL1, hsp, hsp[].state.l1x, tangentSrc, shiftedTangents,
                     preProjBase, projBase, dgEff, preProjOut, tangentOut)

template hypPairStageJVPL2Ptr(hsp, pairSrc, tangentSrc, shiftedTangents,
                              preProjBase, projBase, dgEff,
                              preProjOut, tangentOut: untyped) =
  hypPairStageJVPPtr(psL2, hsp, pairSrc, tangentSrc, shiftedTangents,
                     preProjBase, projBase, dgEff, preProjOut, tangentOut)

template hypPairStageJVPL1*(hs, tangentSrc, shiftedTangents,
                            preProjBase, projBase, dgEff,
                            preProjOut, tangentOut: untyped) =
  hypPairStageJVPL1Ptr(hs.addr, tangentSrc, shiftedTangents,
                       preProjBase, projBase, dgEff, preProjOut, tangentOut)

template hypPairStageJVPL2*(hs, pairSrc, tangentSrc, shiftedTangents,
                            preProjBase, projBase, dgEff,
                            preProjOut, tangentOut: untyped) =
  hypPairStageJVPL2Ptr(hs.addr, pairSrc, tangentSrc, shiftedTangents,
                       preProjBase, projBase, dgEff, preProjOut, tangentOut)

template hypPairStageStapleVJP(stage, hs, pairSrc, chainIn, dst: untyped) =
  block:
    let hsp = hs.addr
    threads:
      forPairs(mu, nu):
        when stage == psL1:
          discard hsp[].kernel.chainShift[nu] ^* chainIn[mu,nu]
          symStapleVJP(dst[nu], dst[mu],
                       hsp[].state.gf[nu], hsp[].state.gf[mu],
                       hsp[].state.gaugePairShift[nu][mu], hsp[].state.gaugePairShift[mu][nu],
                       chainIn[mu,nu], hsp[].kernel.chainShift[nu],
                       hsp[].state.tm1, hsp[].kernel.tm2, hsp[].state.sm1[nu], hsp[].state.sm1[mu])
          threadBarrier()
        else:
          forThird(mu, nu, a, b):
            let pl1 = pairSrc[a,b]
            let pl2 = pairSrc[mu,b]
            threadBarrier()
            discard hsp[].state.stagePairShift[nu][mu] ^*! pl1
            discard hsp[].state.stagePairShift[mu][a] ^*! pl2
            discard hsp[].kernel.chainShift[a] ^*! chainIn[mu,nu]
            threadBarrier()
            symStapleVJP(dst[a,b], dst[mu,b],
                         pl1, pl2,
                         hsp[].state.stagePairShift[nu][mu], hsp[].state.stagePairShift[mu][a],
                         chainIn[mu,nu], hsp[].kernel.chainShift[a],
                         hsp[].state.tm1, hsp[].kernel.tm2, hsp[].state.sm1[a], hsp[].state.sm1[mu])
            threadBarrier()

proc hypPairStageStapleVJPL1*(hs: var HypSmear, chainIn: auto; dst: auto) =
  hs.ensureKernelScratch()
  hypPairStageStapleVJP(psL1, hs, hs.state.l1x, chainIn, dst)

proc hypPairStageStapleVJPL2*(hs: var HypSmear, pairSrc, chainIn: auto;
                              dst: auto) =
  hs.ensureKernelScratch()
  hypPairStageStapleVJP(psL2, hs, pairSrc, chainIn, dst)

template hypPairStageChainVJP(stage, hs, pairSrc, seedBar, dst: untyped) =
  block:
    let hsp = hs.addr
    threads:
      forPairs(mu, nu):
        when stage == psL1:
          symStapleVJPChain(dst[mu,nu],
                            hsp[].state.gf[nu], hsp[].state.gf[mu],
                            hsp[].state.gaugePairShift[nu][mu], hsp[].state.gaugePairShift[mu][nu],
                            seedBar[nu], seedBar[mu],
                            hsp[].kernel.chainShift[nu], hsp[].kernel.seedShift[mu],
                            hsp[].state.tm1, hsp[].kernel.tm2, hsp[].state.sm1[nu])
          threadBarrier()
        else:
          forThird(mu, nu, a, b):
            let pl1 = pairSrc[a,b]
            let pl2 = pairSrc[mu,b]
            threadBarrier()
            discard hsp[].state.stagePairShift[nu][mu] ^*! pl1
            discard hsp[].state.stagePairShift[mu][a] ^*! pl2
            threadBarrier()
            symStapleVJPChain(dst[mu,nu],
                              pl1, pl2,
                              hsp[].state.stagePairShift[nu][mu], hsp[].state.stagePairShift[mu][a],
                              seedBar[a,b], seedBar[mu,b],
                              hsp[].kernel.chainShift[a], hsp[].kernel.seedShift[mu],
                              hsp[].state.tm1, hsp[].kernel.tm2, hsp[].state.sm1[a])
            threadBarrier()

proc hypPairStageChainVJPL1*(hs: var HypSmear, seedBar: auto; dst: auto) =
  hs.ensureKernelScratch()
  hypPairStageChainVJP(psL1, hs, hs.state.l1x, seedBar, dst)

proc hypPairStageChainVJPL2*(hs: var HypSmear, pairSrc, seedBar: auto;
                             dst: auto) =
  hs.ensureKernelScratch()
  hypPairStageChainVJP(psL2, hs, pairSrc, seedBar, dst)

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

template hypPairStageHVP(stage, hs, pairSrc, tangentSrc, chainSrc, dChainIn, dst: untyped) =
  block:
    let hsp = hs.addr
    let tangentSrcp = tangentSrc.addr
    when stage == psL1:
      let shiftedPairTangentsp = hsp[].hvp.hvpL1ShiftedPairTangents.addr
      let shiftedChainTangentsp = hsp[].hvp.hvpL1ShiftedChainTangents.addr
      let acc1p = hsp[].hvp.hvpL1Acc1.addr
      let acc2p = hsp[].hvp.hvpL1Acc2.addr

      threads:
        forPairs(mu, nu):
          discard hsp[].kernel.chainShift[nu] ^* chainSrc[mu,nu]
          discard shiftedPairTangentsp[][mu][nu][0] ^*! tangentSrcp[][nu]
          discard shiftedPairTangentsp[][mu][nu][1] ^*! tangentSrcp[][mu]
          discard shiftedChainTangentsp[][mu][nu] ^*! dChainIn[mu,nu]
          symStapleHVP(dst[nu], dst[mu], acc1p[][mu,nu], acc2p[][mu,nu],
                       hsp[].state.gf[nu], hsp[].state.gf[mu],
                       hsp[].state.gaugePairShift[nu][mu], hsp[].state.gaugePairShift[mu][nu],
                       chainSrc[mu,nu], hsp[].kernel.chainShift[nu],
                       tangentSrcp[][nu], tangentSrcp[][mu], dChainIn[mu,nu],
                       shiftedPairTangentsp[][mu][nu][0], shiftedPairTangentsp[][mu][nu][1],
                       shiftedChainTangentsp[][mu][nu],
                       hsp[].state.tm1, hsp[].kernel.tm2, hsp[].state.sm1[nu], hsp[].state.sm1[mu])
          threadBarrier()
    else:
      let shiftedPairTangentsp = hsp[].hvp.hvpL2ShiftedPairTangents.addr
      let shiftedChainTangentsp = hsp[].hvp.hvpL2ShiftedChainTangents.addr
      let acc1p = hsp[].hvp.hvpL2Acc1.addr
      let acc2p = hsp[].hvp.hvpL2Acc2.addr

      threads:
        forPairs(mu, nu):
          forThird(mu, nu, a, b):
            let pl1 = pairSrc[a,b]
            let pl2 = pairSrc[mu,b]
            threadBarrier()
            discard hsp[].state.stagePairShift[nu][mu] ^*! pl1
            discard hsp[].state.stagePairShift[mu][a] ^*! pl2
            discard hsp[].kernel.chainShift[a] ^*! chainSrc[mu,nu]
            discard shiftedPairTangentsp[][mu][nu][a][0] ^*! tangentSrcp[][a,b]
            discard shiftedPairTangentsp[][mu][nu][a][1] ^*! tangentSrcp[][mu,b]
            discard shiftedChainTangentsp[][mu][nu][a] ^*! dChainIn[mu,nu]
            threadBarrier()
            symStapleHVP(dst[a,b], dst[mu,b],
                         acc1p[][mu][nu][a], acc2p[][mu][nu][a],
                         pl1, pl2,
                         hsp[].state.stagePairShift[nu][mu], hsp[].state.stagePairShift[mu][a],
                         chainSrc[mu,nu], hsp[].kernel.chainShift[a],
                         tangentSrcp[][a,b], tangentSrcp[][mu,b], dChainIn[mu,nu],
                         shiftedPairTangentsp[][mu][nu][a][0], shiftedPairTangentsp[][mu][nu][a][1],
                         shiftedChainTangentsp[][mu][nu][a],
                         hsp[].state.tm1, hsp[].kernel.tm2, hsp[].state.sm1[a], hsp[].state.sm1[mu])
            threadBarrier()

proc hypPairStageHVPL1*(hs: var HypSmear, tangentSrc, chainSrc, dChainIn: auto;
                        dst: auto) =
  hs.ensureHvpScratch()
  hypPairStageHVP(psL1, hs, hs.state.l1x, tangentSrc, chainSrc, dChainIn, dst)

proc hypPairStageHVPL2*(hs: var HypSmear, pairSrc, tangentSrc, chainSrc, dChainIn: auto;
                        dst: auto) =
  hs.ensureHvpScratch()
  hypPairStageHVP(psL2, hs, pairSrc, tangentSrc, chainSrc, dChainIn, dst)

template hypPairStageHVPVJP(stage, hs, pairSrc, chainSrc, seedBar,
                            projected, xfield, preProjChain,
                            pairBarOut, chainBar, nextChainBar: untyped) =
  block:
    let hsp = hs.addr
    let sc = pairStageScales(hsp[], stage)
    let dcbarp = hsp[].hvpAdjoint.hvpVjpDcbar.addr
    threads:
      forPairs(mu, nu):
        dcbarp[][mu,nu] := 0

    threads:
      forPairs(mu, nu):
        when stage == psL1:
          discard hsp[].kernel.chainShift[nu] ^* chainSrc[mu,nu]
          symStapleHVPVJP(pairBarOut[nu], pairBarOut[mu], dcbarp[][mu,nu],
                          hsp[].state.gf[nu], hsp[].state.gf[mu],
                          hsp[].state.gaugePairShift[nu][mu], hsp[].state.gaugePairShift[mu][nu],
                          chainSrc[mu,nu], hsp[].kernel.chainShift[nu],
                          seedBar[nu], seedBar[mu],
                          hsp[].kernel.seedShift[nu], hsp[].kernel.seedShift[mu],
                          hsp[].state.tm1, hsp[].kernel.tm2, hsp[].state.sm1[nu], hsp[].state.sm1[mu])
          threadBarrier()
        else:
          forThird(mu, nu, a, b):
            let pl1 = pairSrc[a,b]
            let pl2 = pairSrc[mu,b]
            threadBarrier()
            discard hsp[].state.stagePairShift[nu][mu] ^*! pl1
            discard hsp[].state.stagePairShift[mu][a] ^*! pl2
            discard hsp[].kernel.chainShift[a] ^*! chainSrc[mu,nu]
            threadBarrier()
            symStapleHVPVJP(pairBarOut[a,b], pairBarOut[mu,b], dcbarp[][mu,nu],
                            pl1, pl2,
                            hsp[].state.stagePairShift[nu][mu], hsp[].state.stagePairShift[mu][a],
                            chainSrc[mu,nu], hsp[].kernel.chainShift[a],
                            seedBar[a,b], seedBar[mu,b],
                            hsp[].kernel.seedShift[a], hsp[].kernel.seedShift[mu],
                            hsp[].state.tm1, hsp[].kernel.tm2, hsp[].state.sm1[a], hsp[].state.sm1[mu])
            threadBarrier()

    addScaledPairs(chainBar, dcbarp[], sc.alp)

    threads:
      forPairs(mu, nu):
        nextChainBar[mu,nu].projHVPVJP_dc(
          projected[mu,nu], xfield[mu,nu], preProjChain[mu,nu],
          chainBar[mu,nu]
        )

proc hypPairStageHVPVJPL1*(hs: var HypSmear, chainSrc, seedBar,
                           projected, xfield, preProjChain: auto;
                           gaugeBarOut, chainBar, nextChainBar: auto) =
  hs.ensureKernelScratch()
  hs.ensureHvpAdjointScratch()
  hypPairStageHVPVJP(psL1, hs, hs.state.l1x, chainSrc, seedBar,
                     projected, xfield, preProjChain,
                     gaugeBarOut, chainBar, nextChainBar)

proc hypPairStageHVPVJPL2*(hs: var HypSmear, pairSrc, chainSrc,
                           seedBar, projected, xfield, preProjChain: auto;
                           pairBarOut, chainBar, nextChainBar: auto) =
  hs.ensureKernelScratch()
  hs.ensureHvpAdjointScratch()
  hypPairStageHVPVJP(psL2, hs, pairSrc, chainSrc, seedBar,
                     projected, xfield, preProjChain,
                     pairBarOut, chainBar, nextChainBar)

template hypPairStageJVPVJP(stage, hs, pairSrc, tangentXBar, pairBarOut, dgEffbar: untyped) =
  block:
    let hsp = hs.addr
    let sc = pairStageScales(hsp[], stage)
    threads:
      forPairs(mu, nu):
        dgEffbar[mu] += sc.ma * tangentXBar[mu,nu]

    threads:
      forPairs(mu, nu):
        when stage == psL1:
          symStapleJVPVJP(dgEffbar[nu], dgEffbar[mu], sc.alp,
                          hsp[].state.gf[nu], hsp[].state.gf[mu],
                          hsp[].state.gaugePairShift[nu][mu], hsp[].state.gaugePairShift[mu][nu],
                          tangentXBar[mu,nu], hsp[].kernel.seedShift[nu],
                          hsp[].state.tm1, hsp[].kernel.tm2, hsp[].state.sm1[nu], hsp[].state.sm1[mu])
          threadBarrier()
        else:
          forThird(mu, nu, a, b):
            let pl1 = pairSrc[a,b]
            let pl2 = pairSrc[mu,b]
            threadBarrier()
            discard hsp[].state.stagePairShift[nu][mu] ^*! pl1
            discard hsp[].state.stagePairShift[mu][a] ^*! pl2
            threadBarrier()
            symStapleJVPVJP(pairBarOut[a,b], pairBarOut[mu,b], sc.alp,
                            pl1, pl2,
                            hsp[].state.stagePairShift[nu][mu], hsp[].state.stagePairShift[mu][a],
                            tangentXBar[mu,nu], hsp[].kernel.seedShift[a],
                            hsp[].state.tm1, hsp[].kernel.tm2, hsp[].state.sm1[a], hsp[].state.sm1[mu])
            threadBarrier()

proc hypPairStageJVPVJPL1*(hs: var HypSmear, tangentXBar: auto;
                           dgEffbar: auto) =
  hs.ensureKernelScratch()
  hypPairStageJVPVJP(psL1, hs, hs.state.l1x, tangentXBar, hs.state.l1x, dgEffbar)

proc hypPairStageJVPVJPL2*(hs: var HypSmear, pairSrc, tangentXBar: auto;
                           pairBarOut, dgEffbar: auto) =
  hs.ensureKernelScratch()
  hypPairStageJVPVJP(psL2, hs, pairSrc, tangentXBar, pairBarOut, dgEffbar)

## Final output-stage operator family built from L2 pair fields.
proc hypOutputStage*(hs: var HypSmear, l2: auto; flx, fl: auto) =
  let hsp = hs.addr
  let sc = outputStageScales(hsp[])
  threads:
    forPairs(mu, nu):
      discard hsp[].state.stagePairShift[mu][nu] ^*! l2[mu,nu]
  threads:
    for mu in 0..<4:
      flx[mu] := sc.ma * hsp[].state.gf[mu]
      forNu(mu, nu):
        let pl1 = l2[nu,mu]
        let pl2 = l2[mu,nu]
        symStaple(flx[mu], sc.alp, pl1, pl2,
                  hsp[].state.stagePairShift[nu][mu], hsp[].state.stagePairShift[mu][nu],
                  hsp[].state.tm1, hsp[].state.sm1[nu])
        threadBarrier()
      fl[mu].proj flx[mu]

template hypOutputStageJVPPtr(hsp, l2, dl2, dsl2, dgEff, dflx: untyped) =
  let dsl2p = dsl2.addr
  let dgEffp = dgEff.addr
  let dflxp = dflx.addr
  let sc = outputStageScales(hsp[])
  threads:
    forPairs(mu, nu):
      discard hsp[].state.stagePairShift[mu][nu] ^*! l2[mu,nu]

  threads:
    for mu in 0..<4:
      dflxp[][mu] := sc.ma * dgEffp[][mu]
      forNu(mu, nu):
        let pl1 = l2[nu,mu]
        let pl2 = l2[mu,nu]
        let dpl1 = dl2[nu,mu]
        let dpl2 = dl2[mu,nu]
        symStapleJVP(dflxp[][mu], sc.alp,
                     pl1, pl2, hsp[].state.stagePairShift[nu][mu], hsp[].state.stagePairShift[mu][nu],
                     dpl1, dpl2, dsl2p[][nu][mu], dsl2p[][mu][nu],
                     hsp[].state.tm1, hsp[].state.sm1[nu])
        threadBarrier()

template hypOutputStageJVP*(hs, l2, dl2, dsl2, dgEff, dflx: untyped) =
  hypOutputStageJVPPtr(hs.addr, l2, dl2, dsl2, dgEff, dflx)

proc buildTangentState[L,F,T](hs: var HypSmear[L,F,T], levels: auto, dg: auto) =
  ## Build the shared forward tangent state used by `smearJVP` and `smearHVP`.
  hs.ensureTangentScratch()
  let hsp = hs.addr

  let dgEffp = hs.tangent.tangentDgEff.addr
  buildGaugeTangent(dgEffp[], hsp[].state.gf, dg)
  let ds1p = hs.tangent.tangentDs1.addr
  primeShiftedTangentsPtr(ds1p, dgEffp)

  let dl1xp = hs.tangent.tangentDl1x.addr
  let dl2xp = hs.tangent.tangentDl2x.addr
  let dl1p = hs.tangent.tangentDl1.addr
  let dl2p = hs.tangent.tangentDl2.addr
  let dsl1L2p = hs.tangent.tangentDsl1L2.addr
  let dsl2p = hs.tangent.tangentDsl2.addr
  let dflxp = hs.tangent.tangentDflx.addr
  clearTangentWorkspace(dflxp[], dl1xp[], dl2xp[], dl1p[], dl2p[])

  hypPairStageJVPL1Ptr(hsp, dgEffp[], ds1p[], hsp[].state.l1x, levels.l1,
                       dgEffp[], dl1xp[], dl1p[])

  hypPairStageJVPL2Ptr(hsp, levels.l1, dl1p[], dsl1L2p[],
                       hsp[].state.l2x, levels.l2, dgEffp[], dl2xp[], dl2p[])

  threads:
    forPairs(mu, nu):
      discard dsl2p[][mu][nu] ^*! dl2p[][mu,nu]

  hypOutputStageJVPPtr(hsp, levels.l2, dl2p[], dsl2p[], dgEffp[], dflxp[])

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

template applyPairProjectionChainAdjoint(hs, chainBar, projected, xfield: auto) =
  let hsp = hs.addr
  let chainBarp = chainBar.addr
  let chainBarTmpp = hsp[].chainAdjoint.pairProjectionTmp.addr
  threads:
    forPairs(mu, nu):
      chainBarTmpp[][mu,nu] := chainBarp[][mu,nu]
      chainBarp[][mu,nu] := 0
      chainBarp[][mu,nu].projVJPChain(projected[mu,nu], xfield[mu,nu], chainBarTmpp[][mu,nu])

template applyProjectionHVP(dst, dchain0, u, x, chain0, dx: auto) =
  ## Apply the projection HVP for a full level while preserving the incoming
  ## chain tangent in `dchain0`.
  let dxp = dx.addr
  threads:
    forPairs(mu, nu):
      dchain0[mu,nu] := dst[mu,nu]
      dst[mu,nu] := 0
      dst[mu,nu].projHVPu(u[mu,nu], x[mu,nu], chain0[mu,nu],
                          dxp[][mu,nu], dchain0[mu,nu])

template applyOutputProjectionHVPPtr(dst, x, chain0, dxp, dchain: auto) =
  threads:
    for mu in 0..<4:
      if dchain.len > 0:
        dst[mu].projHVP(x[mu], chain0[mu], dxp[][mu], dchain[mu])
      else:
        dst[mu].projHVP(x[mu], chain0[mu], dxp[][mu])

template applyOutputProjectionHVP(dst, x, chain0, dx, dchain: auto) =
  applyOutputProjectionHVPPtr(dst, x, chain0, dx.addr, dchain)

template seedSmearHVPAdjoints(dfl1bar, dfl2bar, dfcbar, derivbar,
                              ma1, ma2, ma3: auto) =
  threads:
    for mu in 0..<4:
      dfcbar[mu] += ma3 * derivbar[mu]
      forNu(mu, nu):
        dfl1bar[mu,nu] += ma1 * derivbar[mu]
        dfl2bar[mu,nu] += ma2 * derivbar[mu]

template buildReverseState(hsp, levels, chain, deriv: untyped; accumulateDirect: static bool) =
  ## Build the reverse HYP chain state, optionally accumulating the direct
  ## ma3/ma2/ma1 contributions into `deriv`.
  let sc = hsp[].scaleFactors

  clearReverseChains(hsp[].reverse.outputChain, hsp[].reverse.level1Chain, hsp[].reverse.level2Chain)

  threads:
    for mu in 0..<4:
      hsp[].reverse.outputChain[mu].projVJP(hsp[].state.flx[mu], chain[mu])
      when accumulateDirect:
        deriv[mu] := sc.ma3 * hsp[].reverse.outputChain[mu]
      hsp[].reverse.outputChain[mu] *= sc.alp3

  hypOutputStageStapleVJP(hsp[], levels.l2, hsp[].reverse.outputChain, hsp[].reverse.level2Chain)
  projectAndScalePairChain(hsp[].reverse.level2Chain, levels.l2, hsp[].state.l2x,
                           hsp[].reverse.preProjLevel2Chain, sc.alp2, deriv, sc.ma2,
                           accumulateDirect)
  hypPairStageStapleVJPL2(hsp[], levels.l1, hsp[].reverse.level2Chain, hsp[].reverse.level1Chain)
  projectAndScalePairChain(hsp[].reverse.level1Chain, levels.l1, hsp[].state.l1x,
                           hsp[].reverse.preProjLevel1Chain, sc.alp1, deriv, sc.ma1,
                           accumulateDirect)

proc smear*(hs: var HypSmear, gf: auto, fl: auto) =
  static: doAssert(type(gf[0]) is typeof(hs.state.gf[0]))
  ## Apply HYP smearing; intermediates stored in hs for derivative
  for mu in 0..<4:
    hs.state.gf[mu] = gf[mu]
  var
    l1x = hs.state.l1x
    l2x = hs.state.l2x
    flx = hs.state.flx
  var levels = hs.forwardProjectedLevels()

  hypPairStageL1(hs, l1x, levels.l1)
  hypPairStageL2(hs, levels.l1, l2x, levels.l2)
  hypOutputStage(hs, levels.l2, flx, fl)

proc smearVJP*(hs: var HypSmear, deriv: auto, chain: auto) =
  static:
    doAssert(type(deriv[0]) is typeof(hs.state.gf[0]))
    doAssert(type(chain[0]) is typeof(hs.state.gf[0]))
  ## Backpropagate derivative through hs.smear
  let hsp = hs.addr
  hs.ensureKernelScratch()
  hs.ensureReverseStateScratch()
  let levels = hs.projectedLevels()
  buildReverseState(hsp, levels, chain, deriv, true)
  hypPairStageStapleVJPL1(hsp[], hsp[].reverse.level1Chain, deriv)

template hypOutputProjectionJVPPtr(dst, projected, preProj, tangentp: untyped) =
  threads:
    for mu in 0..<4:
      dst[mu].projJVP(projected[mu], preProj[mu], tangentp[][mu])

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
  let levels = hs.projectedLevels()
  hs.buildTangentState(levels, dg)

  hypOutputProjectionJVPPtr(dsg, sg, hsp.state.flx, hsp[].tangent.tangentDflx.addr)

template hypOutputStageHVPPtr(hsp, l2, dl2, dsl2p, dOutputChain, dLevel2Chain: untyped) =
  block:
    let dl2p = dl2.addr
    let scratchp = hsp[].hvp.level3Hvp.addr

    forLevel3Pairs(hsp[], l2, pl1, pl2):
      discard hsp[].kernel.chainShift[nu] ^*! hsp[].reverse.outputChain[mu]
      discard scratchp[].shiftedOutputChain[mu][nu] ^*! dOutputChain[mu]
      threadBarrier()
      symStapleHVP(dLevel2Chain[nu,mu], dLevel2Chain[mu,nu],
                   scratchp[].acc1[mu,nu], scratchp[].acc2[mu,nu],
                   pl1, pl2, hsp[].state.stagePairShift[nu][mu], hsp[].state.stagePairShift[mu][nu],
                   hsp[].reverse.outputChain[mu], hsp[].kernel.chainShift[nu],
                   dl2p[][nu,mu], dl2p[][mu,nu], dOutputChain[mu],
                   dsl2p[][nu][mu], dsl2p[][mu][nu], scratchp[].shiftedOutputChain[mu][nu],
                   hsp[].state.tm1, hsp[].kernel.tm2, hsp[].state.sm1[nu], hsp[].state.sm1[mu])
      threadBarrier()

template hypOutputStageHVP*(hs, l2, dl2, dsl2, dOutputChain, dLevel2Chain: untyped) =
  hypOutputStageHVPPtr(hs.addr, l2, dl2, dsl2, dOutputChain, dLevel2Chain)

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
  hs.ensureKernelScratch()
  hs.ensureReverseStateScratch()
  hs.ensureHvpScratch()
  let levels = hs.projectedLevels()
  let hsp = hs.addr

  ## Phase 1: analytical forward tangents
  hs.buildTangentState(levels, dg)

  let hvp = hsp[].hvp
  clearHvpWorkspace(dderiv, hvp.dOutputChain, hvp.dLevel1Chain, hvp.dLevel2Chain)

  ## Phase 2: shared reverse-chain state
  buildReverseState(hsp, levels, chain, nil, false)

  ## Phase 3: second-order adjoints through the three HYP levels
  applyOutputProjectionHVPPtr(hvp.dOutputChain, hsp.state.flx, chain,
                              hsp[].tangent.tangentDflx.addr, dchain)
  accumulateAndScaleGauge(dderiv, hvp.dOutputChain, sc.ma3, sc.alp3)
  hypOutputStageHVPPtr(hsp, levels.l2, hsp[].tangent.tangentDl2, hsp[].tangent.tangentDsl2.addr,
                       hvp.dOutputChain, hvp.dLevel2Chain)

  # Level 2: projection on l2x and staples built from l1
  applyPairProjectionHvpPhase(hvp.dLevel2Chain, hvp.dPreProjLevel2Chain, levels.l2, hsp.state.l2x,
                              hsp[].reverse.preProjLevel2Chain, hsp[].tangent.tangentDl2x, dderiv, sc.ma2, sc.alp2)
  hypPairStageHVPL2(hsp[], levels.l1, hsp[].tangent.tangentDl1,
                    hsp[].reverse.level2Chain, hvp.dLevel2Chain, hvp.dLevel1Chain)

  # Level 1: projection on l1x and staples built directly from gf
  applyPairProjectionHvpPhase(hvp.dLevel1Chain, hvp.dPreProjLevel1Chain, levels.l1, hsp.state.l1x,
                              hsp[].reverse.preProjLevel1Chain, hsp[].tangent.tangentDl1x, dderiv, sc.ma1, sc.alp1)
  hypPairStageHVPL1(hsp[], hsp[].tangent.tangentDgEff, hsp[].reverse.level1Chain, hvp.dLevel1Chain, dderiv)

proc smearVJPChain*[G](hs: var HypSmear, chainbar: G, derivbar: G) =
  ## Adjoint of `smearVJP` with respect to `chain`.
  ## Computes `chainbar` such that:
  ##   ⟨derivbar, smearVJP(g, chain)⟩ = ⟨chainbar, chain⟩
  ## using the cached forward HYP state and the reverse-chain phases in reverse order.
  static:
    doAssert(type(chainbar[0]) is typeof(hs.state.gf[0]))
    doAssert(type(derivbar[0]) is typeof(hs.state.gf[0]))

  let sc = hs.scaleFactors
  hs.ensureChainAdjointScratch()
  let levels = hs.projectedLevels()
  let hsp = hs.addr

  # Phase 1: seed the adjoint workspaces
  clearChainAdjointWorkspace(hsp[].chainAdjoint.outputChainBar, chainbar,
                             hsp[].chainAdjoint.level1ChainBar, hsp[].chainAdjoint.level2ChainBar)

  # Phase 2: backpropagate the Level 1 reverse pair stage
  hypPairStageChainVJPL1(hsp[], derivbar, hsp[].chainAdjoint.level1ChainBar)

  # Phase 3: apply the Level 1 direct term and projection-chain adjoint
  threads:
    forPairs(mu, nu):
      hsp[].chainAdjoint.level1ChainBar[mu,nu] *= sc.alp1
      hsp[].chainAdjoint.level1ChainBar[mu,nu] += sc.ma1 * derivbar[mu]
  applyPairProjectionChainAdjoint(hsp[], hsp[].chainAdjoint.level1ChainBar, levels.l1, hsp.state.l1x)

  # Phase 4: backpropagate the Level 2 reverse pair stage
  hypPairStageChainVJPL2(hsp[], levels.l1,
                         hsp[].chainAdjoint.level1ChainBar,
                         hsp[].chainAdjoint.level2ChainBar)

  # Phase 5: apply the Level 2 direct term and projection-chain adjoint
  threads:
    forPairs(mu, nu):
      hsp[].chainAdjoint.level2ChainBar[mu,nu] *= sc.alp2
      hsp[].chainAdjoint.level2ChainBar[mu,nu] += sc.ma2 * derivbar[mu]
  applyPairProjectionChainAdjoint(hsp[], hsp[].chainAdjoint.level2ChainBar, levels.l2, hsp.state.l2x)

  # Phase 6: backpropagate the Level 3 reverse pair stage
  hypOutputStageChainVJP(hsp[], levels.l2,
                         hsp[].chainAdjoint.level2ChainBar,
                         hsp[].chainAdjoint.outputChainBar)

  # Phase 7: apply the Level 3 direct term and finish at the output projection
  threads:
    for mu in 0..<4:
      hsp[].chainAdjoint.outputChainBar[mu] *= sc.alp3
      hsp[].chainAdjoint.outputChainBar[mu] += sc.ma3 * derivbar[mu]
  threads:
    for mu in 0..<4:
      chainbar[mu].projVJPChain(hsp.state.flx[mu], hsp[].chainAdjoint.outputChainBar[mu])

proc hypOutputStageHVPVJP*(hs: var HypSmear, dLevel2ChainBarFromL3, l2: auto;
                           dl2bar, dOutputChainBar: auto) =
  hs.ensureKernelScratch()
  hs.ensureHvpAdjointScratch()
  block:
    let hsp = hs.addr
    let sc = outputStageScales(hsp[])
    let scratchp = hsp[].hvpAdjoint.level3Adjoint.addr
    clearLevel3AdjointScratch(scratchp[].dcbar)

    forLevel3Pairs(hsp[], l2, pl1, pl2):
      discard hsp[].kernel.chainShift[nu] ^*! hsp[].reverse.outputChain[mu]
      threadBarrier()
      symStapleHVPVJP(dl2bar[nu,mu], dl2bar[mu,nu], scratchp[].dcbar[mu],
                      pl1, pl2, hsp[].state.stagePairShift[nu][mu], hsp[].state.stagePairShift[mu][nu],
                      hsp[].reverse.outputChain[mu], hsp[].kernel.chainShift[nu],
                      dLevel2ChainBarFromL3[nu,mu], dLevel2ChainBarFromL3[mu,nu],
                      hsp[].kernel.seedShift[nu], hsp[].kernel.seedShift[mu],
                      hsp[].state.tm1, hsp[].kernel.tm2, hsp[].state.sm1[nu], hsp[].state.sm1[mu])
      threadBarrier()

    addScaledGauge(dOutputChainBar, scratchp[].dcbar, sc.alp)

proc hypOutputStageJVPVJP*(hs: var HypSmear, dOutputTangentBar, l2: auto;
                           dgEffbar, dl2bar: auto) =
  hs.ensureKernelScratch()
  block:
    let hsp = hs.addr
    let sc = outputStageScales(hsp[])
    threads:
      for mu in 0..<4:
        dgEffbar[mu] += sc.ma * dOutputTangentBar[mu]

    forLevel3Pairs(hsp[], l2, pl1, pl2):
      symStapleJVPVJP(dl2bar[nu,mu], dl2bar[mu,nu], sc.alp,
                      pl1, pl2, hsp[].state.stagePairShift[nu][mu], hsp[].state.stagePairShift[mu][nu],
                      dOutputTangentBar[mu], hsp[].kernel.seedShift[nu],
                      hsp[].state.tm1, hsp[].kernel.tm2, hsp[].state.sm1[nu], hsp[].state.sm1[mu])
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
  hs.ensureKernelScratch()
  hs.ensureReverseStateScratch()
  hs.ensureHvpAdjointScratch()
  let levels = hs.projectedLevels()
  let hsp = hs.addr
  buildReverseState(hsp, levels, chain, nil, false)

  let adj = hsp[].hvpAdjoint
  clearHvpAdjointWorkspace(gbar, adj.dflxbar, adj.dfcbar, adj.dgEffbar,
                           adj.dfl1bar, adj.dfl2bar, adj.dl1bar, adj.dl2bar,
                           adj.dl1xbar, adj.dl2xbar)
  clearHvpAdjointIntermediates(adj.dfl1barL2, adj.dfl2barL3)

  # Phase 1: seed the adjoints from the output contribution of `smearHVP`
  seedSmearHVPAdjoints(adj.dfl1bar, adj.dfl2bar, adj.dfcbar, derivbar, sc.ma1, sc.ma2, sc.ma3)

  # Phase 2: reverse the Level 1 and Level 2 HVP output paths
  hypPairStageHVPVJPL1(hsp[], hsp[].reverse.level1Chain, derivbar,
                       levels.l1, hsp.state.l1x, hsp[].reverse.preProjLevel1Chain,
                       adj.dgEffbar, adj.dfl1bar, adj.dfl1barL2)

  hypPairStageHVPVJPL2(hsp[], levels.l1, hsp[].reverse.level2Chain, adj.dfl1barL2,
                       levels.l2, hsp.state.l2x, hsp[].reverse.preProjLevel2Chain,
                       adj.dl1bar, adj.dfl2bar, adj.dfl2barL3)

  # Phase 3: reverse the Level 3 HVP output path back into the tangent state
  hypOutputStageHVPVJP(hsp[], adj.dfl2barL3, levels.l2, adj.dl2bar, adj.dfcbar)
  threads:
    for mu in 0..<4:
      adj.dflxbar[mu].projHVPVJP_dx(hsp.state.flx[mu], chain[mu], adj.dfcbar[mu])
  hypOutputStageJVPVJP(hsp[], adj.dflxbar, levels.l2, adj.dgEffbar, adj.dl2bar)

  # Phase 4: reverse the Level 2 projection and tangent contributions
  accumulateProjectionAdjoint(adj.dl2xbar, hsp.state.l2x, hsp[].reverse.preProjLevel2Chain, adj.dl2bar, adj.dfl2bar)
  hypPairStageJVPVJPL2(hsp[], levels.l1, adj.dl2xbar, adj.dl1bar, adj.dgEffbar)

  # Phase 5: reverse the Level 1 projection and tangent contributions
  accumulateProjectionAdjoint(adj.dl1xbar, hsp.state.l1x, hsp[].reverse.preProjLevel1Chain, adj.dl1bar, adj.dfl1bar)
  hypPairStageJVPVJPL1(hsp[], adj.dl1xbar, adj.dgEffbar)

  # Phase 6: return the adjoint with respect to the multiplicative tangent `dg*g`
  threads:
    for mu in 0..<4:
      gbar[mu] := adj.dgEffbar[mu]
