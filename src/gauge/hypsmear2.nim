import base, layout, gauge, hypsmear, comms/halo, bitops
getOptimPragmas()

#nflop = 61632.0
# l1[mu,nu] : g[mu](0,nu,-nu) g[nu](0,mu)(0,-nu)
# l2[mu,nu] : g[mu](0) l1[mu,b](a,-a) l1[a,b](0,mu)(0,-a)
# fl[mu] : g[mu](0) l2[mu,nu](nu,-nu) l2[nu,mu](0,mu)(0,-nu)
# l2[mu,nu] : g[mu](0)+(0,b,-b)(a,-a)
#             g[a](0,mu)(0,-a)(0,b,-b)
#             g[b](0,mu)(0,-b)(0,a,-a)
# g[0]: (0,+-1,+-1,+-1) (-1,1,+-1,+-1) (-1,+-1,1,+-1) (-1,+-1,+-1,1)

type
  HypTemps*[L,F,T] = object
    gf: array[4,F]
    hl: HaloLayout[L]
    hm: HaloMap[L]
    l1x: FieldArray[L.V,T]
    l2x: FieldArray[L.V,T]
    l1: FieldArray[L.V,T]
    l2: FieldArray[L.V,T]
    flx: FieldArray[L.V,T]
    hgf: array[4,Halo[L,F,T]]
    h1x: array[4,array[4,Halo[L,F,T]]]
    h2x: array[4,array[4,Halo[L,F,T]]]
    h1: array[4,array[4,Halo[L,F,T]]]
    h2: array[4,array[4,Halo[L,F,T]]]
    info: PerfInfo

template proj(x,y: auto) = x.projectU y
proc proj[T](x: T): T {.noInit,alwaysInline.} = result.proj x
template projDeriv(r: auto, x: auto, c: auto) = r.projectUDeriv(x, c)
template projDeriv(r: auto, u, x: auto, c: auto) = r.projectUDeriv(u, x, c)

proc init*[L,F,T,G](ht: var HypTemps[L,F,T], gf: G, comm: Comm) =
  tic "HypTemps init"
  static: doAssert(type(gf[0]) is F)
  static: doAssert(type(gf[0].l) is L)
  static: doAssert(type(gf[0][0]) is T)
  let lo = gf[0].l
  doAssert(lo.nDim == 4)
  ht.hl = lo.makeHaloLayout([1,1,1,1],[1,1,1,1])
  toc "makeHaloLayout"
  ht.l1x = newFieldArray2(lo,F,[4,4],mu!=nu)
  ht.l2x = newFieldArray2(lo,F,[4,4],mu!=nu)
  ht.l1 = newFieldArray2(lo,F,[4,4],mu!=nu)
  ht.l2 = newFieldArray2(lo,F,[4,4],mu!=nu)
  ht.flx = newFieldArray(lo,F,4)
  for mu in 0..<4:
    ht.gf[mu] = gf[mu]
    ht.hgf[mu] = makeHalo(ht.hl, gf[mu])
    for nu in 0..<4:
      if nu == mu: continue
      ht.h1x[mu][nu] = makeHalo(ht.hl, ht.l1x[mu,nu])
      ht.h2x[mu][nu] = makeHalo(ht.hl, ht.l2x[mu,nu])
      ht.h1[mu][nu] = makeHalo(ht.hl, ht.l1[mu,nu])
      ht.h2[mu][nu] = makeHalo(ht.hl, ht.l2[mu,nu])
  toc "makeHalo"
  const hmtype = 0
  when hmtype == 0:
    var offsets = newSeq[array[4,int32]](0)
    for k in 0..15:  # use all corners for simplicity
      var t = [1'i32,1,1,1]
      for i in 0..<4:
        if k.testBit i: t[i] = -1
      offsets.add t
    ht.hm = ht.hl.makeHaloMap(comm, offsets)
  elif hmtype == 1:
    var offsets = newSeq[array[4,int32]](0)
    for k in 0..15:
      var t = [1'i32,1,1,1]
      for i in 0..<4:
        if k.testBit i: t[i] = -1
      if k != 0 and k != 15: offsets.add t
      for mu in 0..<4:
        if t[mu] == 1:
          t[mu] = 0
          offsets.add t
          t[mu] = 1
    ht.hm = hl.makeHaloMap(comm, offsets)
  toc "makeHaloMap"
proc newHypTemps*[G](gf: G): auto =
  type
    L = type gf[0].l
    F = type gf[0]
    T = type gf[0][0]
  var ht: HypTemps[L,F,T]
  let comm = getDefaultComm()
  ht.init(gf, comm)
  ht

proc symstaple(s: var auto, a: float, nu0,mu0,nu1,nu2,mu1,nu3: auto) =
  let t0 = mu0 * nu1.adj
  let t1 = nu0 * t0
  s += a*t1
  let t2 = mu1 * nu3
  let t3 = nu2.adj * t2
  s += a*t3

proc symstaple(s: var auto, a: float, x,y: auto, i,fnu,fmu,bnu,fmubnu: SomeInteger,
               p: static bool = false) =
  when p:
    template prj(m: auto): auto = proj(m)
  else:
    template prj(m: auto): auto = m
  let xi = prj x[i]
  let yfnu = prj y[fnu]
  let xfmu = prj x[fmu]
  let xbnu = prj x[bnu]
  let ybnu = prj y[bnu]
  let xfmubnu = prj x[fmubnu]
  var t = xi * yfnu * xfmu.adj
  t += xbnu.adj * ybnu * xfmubnu
  s += a*t
template symstaplep(s: var auto, a: float, x,y: auto, i,fnu,fmu,bnu,fmubnu: SomeInteger) =
  symstaple(s, a, x,y, i,fnu,fmu,bnu,fmubnu, true)

proc symderiv(s: var auto, x,y,cx,cy: auto, i,fnu,fmu,bnu,fmubnu: SomeInteger,
              p: static bool = false): int =
  mixin projectUflops
  when p:
    template prj(m: auto): auto = proj(m)
  else:
    template prj(m: auto): auto = m
  if fnu>=0 and fmu>=0:
    let xi = prj x[i]
    let yfnu = prj y[fnu]
    let xfmu = prj x[fmu]
    s += cx[i] * yfnu * xfmu.adj
    s += xi * cy[fnu] * xfmu.adj
    s += xi * yfnu * cx[fmu].adj
    result += 3*projectUflops(3) + 3*(18+2*3*66)
  if bnu>=0 and fmubnu>=0:
    let xbnu = prj x[bnu]
    let ybnu = prj y[bnu]
    let xfmubnu = prj x[fmubnu]
    s += cx[bnu].adj * ybnu * xfmubnu
    s += xbnu.adj * cy[bnu] * xfmubnu
    s += xbnu.adj * ybnu * cx[fmubnu]
    result += 3*projectUflops(3) + 3*(18+2*3*66)
template symderivp(s: var auto, x,y,cx,cy: auto, i,fnu,fmu,bnu,fmubnu: SomeInteger): int =
  symderiv(s, x,y,cx,cy, i,fnu,fmu,bnu,fmubnu, true)

proc symderiv3(s: var auto, x,y,cx,cy: auto, i,fnu,fmu,bnu,fmubnu,nOut: SomeInteger,
               p: static bool = false): int =
  mixin projectUflops
  when p:
    template prj(m: auto): auto = proj(m)
  else:
    template prj(m: auto): auto = m
  if fnu>=0 and fmu>=0:
    let xi = prj x[i]
    let yfnu = prj y[fnu]
    let xfmu = prj x[fmu]
    result += 3*projectUflops(3)
    if i<nOut:
      s += cx[i] * yfnu * xfmu.adj
      result += 18+2*3*66
    if fnu<nOut:
      s += xi * cy[fnu] * xfmu.adj
      result += 18+2*3*66
    if fmu<nOUt:
      s += xi * yfnu * cx[fmu].adj
      result += 18+2*3*66
  if bnu>=0 and fmubnu>=0:
    let xbnu = prj x[bnu]
    let ybnu = prj y[bnu]
    result += 2*projectUflops(3)
    if bnu<nOut:
      let xfmubnu = prj x[fmubnu]
      s += cx[bnu].adj * ybnu * xfmubnu
      s += xbnu.adj * cy[bnu] * xfmubnu
      result += projectUflops(3) + 2*(18+2*3*66)
    if fmubnu<nOut:
      s += xbnu.adj * ybnu * cx[fmubnu]
      result += 18+2*3*66
template symderiv3p(s: var auto, x,y,cx,cy: auto, i,fnu,fmu,bnu,fmubnu,nOut: SomeInteger): int =
  symderiv3(s, x,y,cx,cy, i,fnu,fmu,bnu,fmubnu,nOut, true)

proc symderiv2(s: var auto, x,y,cx,cy: auto, i,fnu,fmu,bnu,fmubnu: SomeInteger): int =
  template prj(m: auto): auto = proj(m)
  if fnu>=0 and fmu>=0:
    let xi = prj x[i]
    let yfnu = prj y[fnu]
    let xfmu = prj x[fmu]
    s += cx[i] * yfnu * xfmu.adj
    s += xi * cy[fnu] * xfmu.adj
    s += xi * yfnu * cx[fmu].adj
    inc result, 3
  if bnu>=0 and fmubnu>=0:
    let xbnu = prj x[bnu]
    let ybnu = prj y[bnu]
    let xfmubnu = prj x[fmubnu]
    s += cx[bnu].adj * ybnu * xfmubnu
    s += xbnu.adj * cy[bnu] * xfmubnu
    s += xbnu.adj * ybnu * cx[fmubnu]
    inc result, 3

#[
# x: side links, y: middle links
# snu: c xmu ynu'
# snu: y xmu cnu'
# snu: c' x-mu ynu-mu
# snu: y' x-mu cnu-mu
# smu: x cnu xmu'
# smu: x-nu' c-nu xmu-nu
# x, xmu, x-mu, x-nu, xmu-nu; y, ynu, ynu-mu; c, cnu, c-nu
proc symderiv(f: var auto, a: float, nu0,mu0,nu1,nu2,mu1,nu3,c,cnu: auto) =
  let t0 = mu0 * nu1.adj
  let t1 = nu0 * t0
  f += a*t1
  let t2 = mu1 * nu3
  let t3 = nu2.adj * t2
  f += a*t3
]#

proc smear*[L,F,T,G](ht: HypTemps[L,F,T], coef: HypCoefs, fl: G) =
  tic "HypTemps smear"
  static: doAssert(type(fl[0]) is F)
  static: doAssert(type(fl[0].l) is L)
  static: doAssert(type(fl[0][0]) is T)
  let lo = fl[0].l
  doAssert(lo.nDim == 4)
  let comm = getDefaultComm()
  for mu in 0..<4:
    ht.hgf[mu].halo := 0
    ht.hgf[mu].update ht.hm, comm
  toc "update"
  const V = L.V
  var
    gf = ht.gf
    hl = ht.hl
    l1x = ht.l1x
    l2x = ht.l2x
    flx = ht.flx
    hgf = ht.hgf
    h1x = ht.h1x
    h2x = ht.h2x
    h1 = ht.h1
    h2 = ht.h2
  let
    alp1 = coef.alpha1 / 2.0
    alp2 = coef.alpha2 / 4.0
    alp3 = coef.alpha3 / 6.0
    ma1 = 1 - coef.alpha1
    ma2 = 1 - coef.alpha2
    ma3 = 1 - coef.alpha3
  let nhalo = hl.nExt
  threads:
    let nc = gf[0][0].nrows
    let staplesFlops = float((4*(6*nc+2*(nc-1))+4*2)*nc*nc*V)
    let siteFlops1 = staplesFlops + float(2*nc*nc*V)
    let siteFlops2 = 2*staplesFlops+float((2*nc*nc+12*projectUflops(nc))*V)
    let siteFlops3 = 3*staplesFlops+float((2*nc*nc+19*projectUflops(nc))*V)
    var flops = 0.0
    # h1x[mu][nu] mu,nu,-nu,mu-nu
    tfor i, 0..<nhalo:
      for mu in 0..<4:
        #for nu in 0..<4:
        #  if nu == mu: continue
        #  h1x[mu][nu][i] = 0
        let fmu = hl.neighborFwd[mu][i]
        if fmu<0: continue
        for nu in 0..<4:
          if nu == mu: continue
          let fnu = hl.neighborFwd[nu][i]
          if fnu<0: continue
          let bnu = hl.neighborBck[nu][i]
          if bnu<0: continue
          let fmubnu = hl.neighborBck[nu][fmu]
          if fmubnu<0: continue
          h1x[mu][nu][i] = ma1 * hgf[mu][i]
          symstaple(h1x[mu][nu][i], alp1, hgf[nu],hgf[mu], i,fnu,fmu,bnu,fmubnu)
          h1[mu][nu][i].proj h1x[mu][nu][i]
          flops += siteFlops1
    #threadBarrier()
    flops.threadSum
    toc("1",flops=flops)
    flops = 0
    # h2x[mu][nu]: h1x[a][b]: 0,mu,-a,mu-a  h1x[mu][b]: a,-a
    tfor i, 0..<nhalo:
      for mu in 0..<4:
        #for nu in 0..<4:
        #  if nu == mu: continue
        #  h2x[mu][nu][i] = 0
        let fmu = hl.neighborFwd[mu][i]
        if fmu<0: continue
        for nu in 0..<4:
          if nu == mu: continue
          h2x[mu][nu][i] = ma2 * hgf[mu][i]
          for a in 0..<4:
            if a == mu or a == nu: continue
            let b = 1+2+3-mu-nu-a
            let fa = hl.neighborFwd[a][i]
            if fa<0: continue
            let ba = hl.neighborBck[a][i]
            if ba<0: continue
            let fmuba = hl.neighborBck[a][fmu]
            if fmuba<0: continue
            #symstaplep(h2x[mu][nu][i], alp2, h1x[a][b],h1x[mu][b], i,fa,fmu,ba,fmuba)
            symstaple(h2x[mu][nu][i], alp2, h1[a][b],h1[mu][b], i,fa,fmu,ba,fmuba)
          h2[mu][nu][i].proj h2x[mu][nu][i]
          flops += siteFlops2
    #threadBarrier()
    flops.threadSum
    toc("2",flops=flops)
    flops = 0
    # flx[mu]: h2x[nu][mu]: 0,mu,-nu,mu-nu  h2x[mu][nu]: nu,-nu
    for i in lo:
      for mu in 0..<4:
        flx[mu][i] = ma3 * gf[mu][i]
        let fmu = hl.neighborFwd[mu][i]
        for nu in 0..<4:
          if nu == mu: continue
          let fnu = hl.neighborFwd[nu][i]
          let bnu = hl.neighborBck[nu][i]
          let fmubnu = hl.neighborBck[nu][fmu]
          #symstaplep(flx[mu][i], alp3, h2x[nu][mu],h2x[mu][nu], i,fnu,fmu,bnu,fmubnu)
          symstaple(flx[mu][i], alp3, h2[nu][mu],h2[mu][nu], i,fnu,fmu,bnu,fmubnu)
        fl[mu][i].proj flx[mu][i]
        flops += siteFlops3
    #threadBarrier()
    flops.threadSum
    toc("3",flops=flops)
  toc("threads end")

proc smear2*[G](coef: HypCoefs, gf: G, fl: G, info: var PerfInfo) =
  tic "smear2"
  type
    F = type gf[0]
    L = type gf[0].l
    T = type gf[0][0]
  let comm = getDefaultComm()
  var ht: HypTemps[L,F,T]
  ht.init gf, comm
  toc "init"
  ht.smear coef, fl
  toc "smear"

proc force*[L,F,T,G,C](ht: HypTemps[L,F,T], coef: HypCoefs, f: G, chain: C) =
  tic "HypTemps force"
  static: doAssert(type(f[0]) is F)
  static: doAssert(type(f[0].l) is L)
  static: doAssert(type(f[0][0]) is T)
  let lo = f[0].l
  doAssert(lo.nDim == 4)
  let comm = getDefaultComm()
  const V = L.V
  var
    gf = ht.gf
    hl = ht.hl
    l1x = ht.l1x
    l2x = ht.l2x
    flx = ht.flx
    hgf = ht.hgf
    h1x = ht.h1x
    h2x = ht.h2x
    h1 = ht.h1
    h2 = ht.h2
    fc = newFieldArray(lo,F,4)
    fl1 = newFieldArray2(lo,F,[4,4],mu!=nu)
    fl2 = newFieldArray2(lo,F,[4,4],mu!=nu)
    hf: array[4,Halo[L,F,T]]
    hfc: array[4,Halo[L,F,T]]
    hl1: array[4,array[4,Halo[L,F,T]]]
    hl2: array[4,array[4,Halo[L,F,T]]]
  let
    alp1 = coef.alpha1 / 2.0
    alp2 = coef.alpha2 / 4.0
    alp3 = coef.alpha3 / 6.0
    ma1 = 1 - coef.alpha1
    ma2 = 1 - coef.alpha2
    ma3 = 1 - coef.alpha3
  for mu in 0..<4:
    hf[mu] = makeHalo(ht.hl, f[mu])
    hfc[mu] = makeHalo(ht.hl, fc[mu])
    for nu in 0..<4:
      if nu == mu: continue
      hl1[mu][nu] = makeHalo(ht.hl, fl1[mu,nu])
      hl2[mu][nu] = makeHalo(ht.hl, fl2[mu,nu])
  let nhalo = hl.nExt
  let nOut = hl.nOut
  threads:
    let nc = gf[0][0].nrows
    var flops = 0.0
    for mu in 0..<4:
      for i in lo:
        fc[mu][i].projDeriv(flx[mu][i], chain[mu][i])
        f[mu][i] = ma3 * fc[mu][i]
        fc[mu][i] *= alp3
      tfor i, nOut..<nhalo:
        hf[mu][i] = 0
        hfc[mu][i] = 0
    flops.threadSum
    toc "fc"
    flops = 0
    tfor i, 0..<nhalo:
      for mu in 0..<4:
        let fmu = hl.neighborFwd[mu][i]
        for nu in 0..<4:
          if nu==mu: continue
          hl2[mu][nu][i] = 0
          let fnu = hl.neighborFwd[nu][i]
          let bnu = hl.neighborBck[nu][i]
          let fmubnu = if fmu<0: -1 else: hl.neighborBck[nu][fmu]
          #let flp = symderiv3p(hl2[mu][nu][i], h2x[nu][mu],h2x[mu][nu],hfc[nu],hfc[mu], i,fnu,fmu,bnu,fmubnu,nOut)
          let flp = symderiv3(hl2[mu][nu][i], h2[nu][mu],h2[mu][nu],hfc[nu],hfc[mu], i,fnu,fmu,bnu,fmubnu,nOut)
          if flp > 0:
            hl2[mu][nu][i].projDeriv(h2x[mu][nu][i], hl2[mu][nu][i])
            hf[mu][i] += ma2 * hl2[mu][nu][i]
            hl2[mu][nu][i] *= alp2
            flops += V*(flp+54+2*projectUflops(nc)) # guess for projDeriv
    flops.threadSum
    toc "3", flops=flops
    flops = 0
    tfor i, 0..<nhalo:
      for mu in 0..<4:
        let fmu = hl.neighborFwd[mu][i]
        var flp = 0
        for nu in 0..<4:
          if nu == mu: continue
          hl1[mu][nu][i] = 0
          #var t {.noInit.}: evalType(hl1[mu][nu][i])
          #t := 0
          for a in 0..<4:
            if a == mu or a == nu: continue
            let b = 1+2+3-mu-nu-a
            let fa = hl.neighborFwd[a][i]
            let ba = hl.neighborBck[a][i]
            let fmuba = if fmu<0: -1 else: hl.neighborBck[a][fmu]
            #flp += symderivp(hl1[mu][nu][i], h1x[a][nu],h1x[mu][nu],hl2[a][b],hl2[mu][b], i,fa,fmu,ba,fmuba)
            flp += symderiv(hl1[mu][nu][i], h1[a][nu],h1[mu][nu],hl2[a][b],hl2[mu][b], i,fa,fmu,ba,fmuba)
            #flp += symderiv(t, h1[a][nu],h1[mu][nu],hl2[a][b],hl2[mu][b], i,fa,fmu,ba,fmuba)
          if flp > 0:
            hl1[mu][nu][i].projDeriv(h1x[mu][nu][i], hl1[mu][nu][i])
            hf[mu][i] += ma1 * hl1[mu][nu][i]
            hl1[mu][nu][i] *= alp1
            #t.projDeriv(h1x[mu][nu][i], t)
            #hf[mu][i] += ma1 * t
            #hl1[mu][nu][i] = alp1 * t
            flops += V*(flp+54+2*projectUflops(nc)) # guess for projDeriv
    flops.threadSum
    toc "2", flops=flops
    flops = 0
    tfor i, 0..<nhalo:
      for mu in 0..<4:
        let fmu = hl.neighborFwd[mu][i]
        for nu in 0..<4:
          if nu == mu: continue
          let fnu = hl.neighborFwd[nu][i]
          let bnu = hl.neighborBck[nu][i]
          let fmubnu = if fmu<0: -1 else: hl.neighborBck[nu][fmu]
          let flp = symderiv(hf[mu][i], hgf[nu],hgf[mu],hl1[nu][mu],hl1[mu][nu], i,fnu,fmu,bnu,fmubnu)
          flops += V*flp
    flops.threadSum
    toc "1", flops=flops
  for mu in 0..<4:
    hf[mu].updateRev ht.hm, comm
  toc("end")

when isMainModule:
  import qex
  qexInit()
  tic("main")
  var defaultLat = @[4,4,4,4]
  defaultSetup()
  var seed = 987654321'u
  #var rng = newRngField(lo, RngMilc6, seed)
  var rng = newRngField(lo, MRG32k3a, seed)
  g.gaussian rng
  toc "gaussian"
  echo 6.0 * g.plaq
  toc "plaq"

  var info: PerfInfo
  var coef: HypCoefs
  coef.alpha1 = 0.4
  coef.alpha2 = 0.5
  coef.alpha3 = 0.5
  echo coef

  var fl = lo.newGauge()
  var fl2 = lo.newGauge()

  proc testSmear =
    resetTimers()
    tic "start"
    coef.smear(g, fl, info)
    toc "smear"
    coef.smear2(g, fl2, info)
    toc "smear2"
    echoTimers()
    echo fl.plaq
    echo fl2.plaq
  #testSmear()

  proc testForce =
    var f = lo.newGauge()
    var f2 = lo.newGauge()
    var c = lo.newGauge()
    c.gaussian rng
    block:
      let fn = coef.smearGetForce(g, fl, info)
      fn(f, c)
    resetTimers()
    block:
      tic()
      let fn = coef.smearGetForce(g, fl, info)
      toc "smear"
      fn(f, c)
      toc "force"
    block:
      tic()
      let ht = newHypTemps(g)
      toc "newHypTemps"
      freezeTimers()
      ht.smear(coef, fl2)
      ht.force(coef, f2, c)
      thawTimers()
      block:
        tic()
        ht.smear(coef, fl2)
        toc "HTsmear"
        ht.force(coef, f2, c)
        toc "HTforce"
    echoTimers()
    echo fl.plaq
    echo fl2.plaq
    echo f.plaq
    echo f2.plaq
    var d = fl[0].newOneOf
    var d2fl = newSeq[float](fl.len)
    var d2f = newSeq[float](f.len)
    for mu in 0..<fl.len:
      d := fl2[mu] - fl[mu]
      d2fl[mu] = sqrt(d.norm2/fl[mu].norm2)
      d := f2[mu] - f[mu]
      d2f[mu] = sqrt(d.norm2/f[mu].norm2)
    echo "error smear: ", d2fl
    echo "error force: ", d2f

  testForce()

  qexFinalize()
