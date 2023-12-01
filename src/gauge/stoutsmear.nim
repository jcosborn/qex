import base
import layout
import gauge

type StoutSmear*[G] = object
  alpha*: float
  gf*,f*,expaf*,ds*,cg*:G
  # add other temporaries

proc newStoutSmear*(l:Layout, alpha:float):auto =
  # we only save a reference in gf, so no allocation to gf here.
  type G = type(l.newGauge)
  StoutSmear[G](alpha:alpha, f:l.newGauge, expaf:l.newGauge, ds:l.newGauge, cg:l.newGauge)

proc smear*[G](ss:var StoutSmear, gf:G, fl:G) =
  const nc = gf[0][0].nrows.float
  let
    alpha = -ss.alpha*nc  # negative from gaugeForce, and nc compensate force normalization
    f = ss.f
    expaf = ss.expaf
    ds = ss.ds
  ss.gf = gf

  gaugeActionDeriv(GaugeActionCoeffs(plaq:1.0), gf, ds)
  threads:
    for mu in 0..<f.len:
      for e in f[mu]:
        let s = gf[mu][e]*ds[mu][e].adj
        var t{.noinit.}: evalType(f[mu][e])
        t.projectTAH s
        f[mu][e] := t
        t := exp(alpha*t)
        expaf[mu][e] := t
        fl[mu][e] := t*gf[mu][e]

proc inverse*[G](ss:var StoutSmear, gf:G, fl:G, rdf2req=1e-24, maxIter=1000, verbose=false):auto =
  ## perform backward smear, fl:input, gf:output, gf and fl must be distinct
  ## return the number of iterations and relative force norm2
  ## note: the resulting ss does not have the necessary fields for smearDeriv
  ## This uses the fixed point iteration from eq. 5.6 and 5.7 in Trivializing Maps, the Wilson Flow and the HMC Algorithm, Luscher, 2010
  ## However, since we update the gauge links all at once, the bound |ε|<⅛ (eq. 5.5) is only a necessary condition.
  ## The iteration may diverge even if |ε|<⅛, depending on the gauge links.
  const nc = gf[0][0].nrows.float
  let
    alpha = ss.alpha*nc  # backward cancels negative from gaugeForce, and nc compensate force normalization
    f = ss.f
    ds = ss.ds
  ss.gf = gf

  threads:
    for mu in 0..<gf.len:
      gf[mu] := fl[mu]
      f[mu] := 0

  var iter = 0
  var rdf2:float
  var df2o = -1.0
  while iter<maxIter:
    inc iter
    gaugeActionDeriv(GaugeActionCoeffs(plaq:1.0), gf, ds)
    threads:
      var df2 = 0.0
      var f2 = 0.0
      for mu in 0..<f.len:
        for e in ds[mu]:
          let s = gf[mu][e]*ds[mu][e].adj
          var t{.noinit.}: evalType(ds[mu][e])
          t.projectTAH s
          df2 += norm2(t-f[mu][e]).simdSum
          f2 += t.norm2.simdSum
          f[mu][e] := t
          t := exp(alpha*t)
          gf[mu][e] := t*fl[mu][e]
      threadBarrier()
      threadRankSum df2
      threadRankSum f2
      threadMaster:
        rdf2 = df2/f2
        if verbose:
          echo iter," ",df2," ",f2," ",rdf2
        if df2o>=0 and df2o<df2:
          # this signals the iterative solver is diverging
          qexWarn "df^2 increased: iter ",iter," current ",df2," previous ",df2o
        df2o = df2
    if verbose:
      gf.echoPlaq
    if rdf2<rdf2req:
      break
  return (iter:iter, rdf2:rdf2)

# backpropagation
# z = f(g(h(x)))
# dz/dx_i = df/dg_j dg_j/dh_k dh_k/dx_i
# Dg_j = df/dg_j, Dh_k = Dg_j dg_j/dh_k, Dx_i = Dh_k dh_k/dx_i
# convention: z = 1/2 Tr(g s† + s g†), dz/dg_ij = s†_ji = D[g]†

proc gaugeForceDeriv*[G](gf:G, deriv:G, chain:G, f:G, cg:G) =
  ## gf: gauge field in
  ## deriv: derivative out
  ## chain: back propagated derivative in
  ## f: action derivative in
  ## cg: scratch space for intermediate field
  # f_x = T_a ReTr(T_a g_x sum(staples))
  # deriv_x ← d/d(g_x)[f_x] chain_x† + sum_y d/d(g_x)[f_y] chain_y†
  let nd = gf[0].l.nDim

  # d/d(g_x)[f_x] c_x† = d/d(g_x)[T_a ReTr(T_a g_x B)] c_x†
  #                    = 1/2 d/d(g_x)_ij[(g_x B T_a + c.c.)] T_akl c_x†_lk
  #                    = 1/2 (B T_a)†_ji Tr(T_a c_x†)
  # redefine the chain, add the complex conjugate back
  # D[g_x]† = (B projTAH(c_x†))†, B: x+ν→x
  # D[g_x] = projTAH(c_x) B†
  threads:
    for mu in 0..<nd:
      for e in gf[mu]:
        var t{.noinit.}: evalType(chain[mu][e])
        t.projectTAH chain[mu][e]
        deriv[mu][e] := t*f[mu][e]
        cg[mu][e] := t.adj*gf[mu][e]

  # d/d(g_x)[f_y] c_y† = d/d(g_x)[T_a ReTr(T_a g_y A g_x B)] c_y†
  #                    = 1/2 d/d(g_x)_ij[(g_x B T_a g_y A + c.c.)] T_akl c_y†_lk
  #                    = 1/2 (B T_a g_y A)†_ji Tr(T_a c_y†)
  # redefine the chain, add the complex conjugate back
  # D[g_x]† = (B projTAH(c_y†) g_y A)†, A: y+μ→x, B: x+ν→y
  # D[g_x] = A† (projTAH(c_y†) g_y)†  B†, staples: x → y+μ →[g_y† cTAH_μ] y → x+ν, with cTAH_μ insertion
  const nc = gf[0][0].nrows
  let cp = 1.0 / float(nc)    # normalization follows gaugeForce
  let
    t = newTransporters(gf, gf[0], 1)
    td = newTransporters(gf, gf[0], -1)
    ct = newTransporters(cg, cg[0], 1)
    ctd = newTransporters(cg, cg[0], -1)
  threads:
    for mu in 0..<nd:
      for nu in 0..<nd:
        if nu==mu: continue
        discard t[nu] ^* gf[mu]
        shiftExpr(t[mu].sb, deriv[mu][ir] += cp * t[nu].field[ir]*adj(it), cg[nu][ix])
        discard t[nu] ^* cg[mu]
        shiftExpr(t[mu].sb, deriv[mu][ir] += cp * t[nu].field[ir]*adj(it), gf[nu][ix])
        discard ct[nu] ^* gf[mu]
        shiftExpr(t[mu].sb, deriv[mu][ir] += cp * ct[nu].field[ir]*adj(it), gf[nu][ix])
        deriv[mu] += cp * td[nu] ^* t[mu] ^* cg[nu]
        deriv[mu] += cp * td[nu] ^* ct[mu] ^* gf[nu]
        deriv[mu] += cp * ctd[nu] ^* t[mu] ^* gf[nu]

proc smearDeriv*[G](ss:var StoutSmear, deriv:G, chain:G) =
  ## gf: gauge field in
  ## deriv: derivative out
  ## chain: back propagated derivative in
  ## chain and deriv much be distinct fields
  const nc = chain[0][0].nrows.float
  let
    alpha = -ss.alpha*nc  # negative from gaugeForce, and nc compensate force normalization
    gf = ss.gf
    f = ss.f
    expaf = ss.expaf

  # convention:
  # D[exp(a f)]† = d/d(exp(a f))[exp(a f) g] c†
  # D[exp(a f)]† = g c†
  # D[a f]† = d/d(a f)[exp(a f)] D[exp(a f)]†
  # D[f]† = d/d(f)[a f] D[a f]†
  threads:
    for mu in 0..<deriv.len:
      for e in deriv[mu]:
        deriv[mu][e] := alpha*expDeriv(alpha*f[mu][e], chain[mu][e]*gf[mu][e].adj)
  gaugeForceDeriv(ss.gf, deriv, deriv, ss.ds, ss.cg)

  # D[g]† = d/d(g)[exp(a f) g] c† = c† exp(a f)
  threads:
    for mu in 0..<deriv.len:
      for e in deriv[mu]:
        deriv[mu][e] += expaf[mu][e].adj*chain[mu][e]

when isMainModule:
  import qex
  import os

  proc reunit(g:auto) =
    tic()
    threads:
      let d = g.checkSU
      threadBarrier()
      echo "unitary deviation avg: ",d.avg," max: ",d.max
      g.projectSU
      threadBarrier()
      let dd = g.checkSU
      echo "new unitary deviation avg: ",dd.avg," max: ",dd.max
    toc("reunit")

  proc mplaq(g:auto, label="") =
    tic()
    let
      pl = g.plaq
      nl = pl.len div 2
      ps = pl[0..<nl].sum * 2.0
      pt = pl[nl..^1].sum * 2.0
    echo "MEASplaq ",label," ss: ",ps,"  st: ",pt,"  tot: ",0.5*(ps+pt)
    toc("plaq")

  qexinit()
  tic()
  letParam:
    gaugefile = ""
    savefile =
      if gaugefile.len > 0:
        gaugefile & ".smear.lime"
      else:
        "config.smear.lime"
    lat =
      if fileExists(gaugefile):
        getFileLattice gaugefile
      else:
        if gaugefile.len > 0:
          qexWarn "Nonexistent gauge file: ", gaugefile
        @[8,8,8,8]
    steps = @[0.1]    # a list of smearing steps
    backward:bool = 0    # inverse flow, from the last to the first in steps
    maxiter = 1000
    r2req = 1e-24
    verbose:bool = 0
    reunitarize:bool = 1
    showTimers:bool = 0
    timerWasteRatio = 0.05
    timerEchoDropped:bool = 0
    timerExpandRatio = 0.05
    verboseGCStats:bool = 0
    verboseTimer:bool = 0

  echoParams()
  echo "rank ", myRank, "/", nRanks
  threads: echo "thread ", threadNum, "/", numThreads

  DropWasteTimerRatio = timerWasteRatio
  VerboseGCStats = verboseGCStats
  VerboseTimer = verboseTimer

  let lo = lat.newLayout
  var gf = lo.newGauge
  if fileExists(gaugefile):
    tic("load")
    if 0 != gf.loadGauge gaugefile:
      qexError "failed to load gauge file: ", gaugefile
    qexLog "loaded gauge from file: ", gaugefile," secs: ",getElapsedTime()
    toc("read")
    if reunitarize:
      gf.reunit
    toc("reunit")
  else:
    gf.random
  gf.mplaq

  if backward:
    var fg = lo.newGauge
    for i in countdown(len(steps)-1, 0):
      var ss = lo.newStoutSmear(steps[i])
      let (iter, r2) = ss.inverse(fg, gf, maxIter=maxiter, rdf2req=r2req, verbose=verbose)
      if iter>=maxIter:
        qexWarn "maximum iteration count reached"
      qexLog "inverse ",i," t ",steps[i]," iter ",iter," r2 ",r2
      for mu in 0..<fg.len:
        let f = fg[mu]
        fg[mu] = gf[mu]    # use = to copy ref only
        gf[mu] = f
      gf.mplaq $i
  else:
    for i in 0..<len(steps):
      var ss = lo.newStoutSmear(steps[i])
      ss.smear(gf, gf)
      gf.mplaq $i

  block:
    tic("save")
    if 0 != gf.saveGauge(savefile):
      qexError "Failed to save gauge to file: ",savefile
    qexLog "saved gauge to file: ",savefile," secs: ",getElapsedTime()
    toc("done")

  toc("done")
  qexFinalize()
