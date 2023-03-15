import base, layout, gauge, hypsmear, comms/halo, bitops
getOptimPragmas()

#nflop = 61632.0
#dtime = 0.0
# l1[mu,nu] : g[mu](0,nu,-nu) g[nu](0,mu)(0,-nu)
# l2[mu,nu] : g[mu](0) l1[mu,b](a,-a) l1[a,b](0,mu)(0,-a)
# fl[mu] : g[mu](0) l2[mu,nu](nu,-nu) l2[nu,mu](0,mu)(0,-nu)
# l2[mu,nu] : g[mu](0)+(0,b,-b)(a,-a)
#             g[a](0,mu)(0,-a)(0,b,-b)
#             g[b](0,mu)(0,-b)(0,a,-a)
# g[0]: (0,+-1,+-1,+-1) (-1,1,+-1,+-1) (-1,+-1,1,+-1) (-1,+-1,+-1,1)

template proj(x,y: auto) = x.projectU y
proc proj[T](x: T): T {.noInit,alwaysInline.} = result.proj x

proc smear2*[G](coef: HypCoefs, gf: G, fl: G, info: var PerfInfo) =
  tic "smear"
  let lo = gf[0].l
  doAssert(lo.nDim == 4)
  let hl = lo.makeHaloLayout([1,1,1,1],[1,1,1,1])
  toc "makeHaloLayout"
  let comm = getDefaultComm()
  type HM = HaloMap[type lo]
  var hm: array[4,HM]
  const hmtype = 0
  when hmtype == 0:
    var offsets = newSeq[array[4,int32]](0)
    for k in 0..15:  # use all corners for simplicity
      var t = [1'i32,1,1,1]
      for i in 0..3:
        if k.testBit i: t[i] = -1
      offsets.add t
    hm[0] = hl.makeHaloMap(comm, offsets)
    for i in 1..3: hm[i] = hm[0]
  elif hmtype == 1:
    var offsets = newSeq[array[4,int32]](0)
    for k in 0..15:  # use all corners for simplicity
      var t = [1'i32,1,1,1]
      for i in 0..3:
        if k.testBit i: t[i] = -1
      if k != 0 and k != 15: offsets.add t
      for mu in 0..3:
        if t[mu] == 1:
          t[mu] = 0
          offsets.add t
          t[mu] = 1
    hm[0] = hl.makeHaloMap(comm, offsets)
    for i in 1..3: hm[i] = hm[0]
  else:
    discard
  toc "makeHaloMap"
  type LCM = type gf[0]
  type H = type makeHalo(hm[0], gf[0])
  var
    l1x = newFieldArray2(lo,LCM,[4,4],mu!=nu)
    l2x = newFieldArray2(lo,LCM,[4,4],mu!=nu)
    flx = newFieldArray(lo,LCM,4)
    hgf: array[4,H]
    h1x: array[4,array[4,H]]
    h2x: array[4,array[4,H]]
  for mu in 0..3:
    hgf[mu] = makeHalo(hm[mu], gf[mu])
    for nu in 0..3:
      if nu == mu: continue
      h1x[mu][nu] = makeHalo(hm[0], l1x[mu,nu])
      h2x[mu][nu] = makeHalo(hm[0], l2x[mu,nu])
  toc "makeHalo"
  for mu in 0..3:
    hgf[mu].update comm
  toc "update"
  let
    alp1 = coef.alpha1 / 2.0
    alp2 = coef.alpha2 / 4.0
    alp3 = coef.alpha3 / 6.0
    ma1 = 1 - coef.alpha1
    ma2 = 1 - coef.alpha2
    ma3 = 1 - coef.alpha3
  let nhalo = hl.nExt
  #for mu in 0..3:
  #  for nu in 0..3:
  #    if nu == mu: continue
  #    l1x[mu,nu] := 0
  #toc "zero"
  threads:
    let nc = gf[0][0].nrows
    let staplesFlops = float((4*(6*nc+2*(nc-1))+4*2)*nc*nc*lo.V)
    let siteFlops1 = staplesFlops + float(2*nc*nc*lo.V)
    let siteFlops2 = 2*staplesFlops+float((2*nc*nc+12*projectUflops(nc))*lo.V)
    let siteFlops3 = 3*staplesFlops+float((2*nc*nc+19*projectUflops(nc))*lo.V)
    var flops = 0.0
    tfor i, 0..<nhalo:
      for mu in 0..3:
        let fmu = hl.neighborFwd[mu][i]
        if fmu<0: continue
        for nu in 0..3:
          if nu == mu: continue
          let fnu = hl.neighborFwd[nu][i]
          if fnu<0: continue
          let bnu = hl.neighborBck[nu][i]
          if bnu<0: continue
          let fmubnu = hl.neighborBck[nu][fmu]
          if fmubnu<0: continue
          h1x[mu][nu][i] = ma1 * hgf[mu][i]
          let t0 = hgf[mu][fnu] * hgf[nu][fmu].adj
          let t1 = hgf[nu][i] * t0
          h1x[mu][nu][i] += alp1*t1
          let t2 = hgf[mu][bnu] * hgf[nu][fmubnu]
          let t3 = hgf[nu][bnu].adj * t2
          h1x[mu][nu][i] += alp1*t3
          flops += siteFlops1
    #threadBarrier()
    flops.threadSum
    toc("1",flops=flops)
    flops = 0
    tfor i, 0..<nhalo:
      for mu in 0..3:
        let fmu = hl.neighborFwd[mu][i]
        if fmu<0: continue
        for nu in 0..3:
          if nu == mu: continue
          h2x[mu][nu][i] = ma2 * hgf[mu][i]
          for a in 0..3:
            if a == mu or a == nu: continue
            let b = 1+2+3-mu-nu-a
            let fa = hl.neighborFwd[a][i]
            if fa<0: continue
            let ba = hl.neighborBck[a][i]
            if ba<0: continue
            let fmuba = hl.neighborBck[a][fmu]
            if fmuba<0: continue
            let la0 = proj h1x[a][b][i]
            let lmu0 = proj h1x[mu][b][fa]
            let la1 = proj h1x[a][b][fmu]
            let t0 = lmu0 * la1.adj
            let t1 = la0 * t0
            h2x[mu][nu][i] += alp2*t1
            let la2 = proj h1x[a][b][ba]
            let lmu1 = proj h1x[mu][b][ba]
            let la3 = proj h1x[a][b][fmuba]
            let t2 = lmu1 * la3
            let t3 = la2.adj * t2
            h2x[mu][nu][i] += alp2*t3
          flops += siteFlops2
    #threadBarrier()
    flops.threadSum
    toc("2",flops=flops)
    flops = 0
    for i in lo:
      for mu in 0..3:
        flx[mu][i] = ma3 * gf[mu][i]
        let fmu = hl.neighborFwd[mu][i]
        for nu in 0..3:
          if nu == mu: continue
          let fnu = hl.neighborFwd[nu][i]
          let bnu = hl.neighborBck[nu][i]
          let fmubnu = hl.neighborBck[nu][fmu]
          let lnu0 = proj l2x[nu,mu][i]
          let lmu0 = proj h2x[mu][nu][fnu]
          let lnu1 = proj h2x[nu][mu][fmu]
          let t0 = lmu0 * lnu1.adj
          let t1 = lnu0 * t0
          flx[mu][i] += alp3*t1
          let lnu2 = proj h2x[nu][mu][bnu]
          let lmu1 = proj h2x[mu][nu][bnu]
          let lnu3 = proj h2x[nu][mu][fmubnu]
          let t2 = lmu1 * lnu3
          let t3 = lnu2.adj * t2
          flx[mu][i] += alp3*t3
        fl[mu][i].proj flx[mu][i]
        flops += siteFlops3
    flops.threadSum
    toc("3",flops=flops)
  toc("threads end")

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

  resetTimers()
  toc "start"
  coef.smear(g, fl, info)
  toc "smear"
  coef.smear2(g, fl2, info)
  toc "smear2"
  echoTimers()

  echo fl.plaq
  echo fl2.plaq

  qexFinalize()
