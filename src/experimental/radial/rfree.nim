## WP-K free-limit campaign (Tier 1, arXiv:2510.03085): the deterministic headline
## numbers and figure data of the fermion and gauge sectors.
##
## Everything here is free-field, so both sectors block-diagonalize over temporal
## Matsubara modes and every quantity is dense linear algebra per mode:
##
## Fermion (antiperiodic, k = (2n+1) pi/Lt).  The spatial 2N_V x 2N_V Wilson block
## D_spatial comes from denseDw on an nt = 1 lattice with parts = dwSpace (the same
## construction twilson pins against the full dense operator to 5.8e-14), and
##   D_W(k) = D_spatial + kappa'_y [(1 - e^{-ik}) P_up + (1 - e^{+ik}) P_dn].
## Per mode: X = D_W(k) - M, H = X^dag X, zheev gives H^{-1/2}, D_ov(k) = 1 + X H^{-1/2},
## one zgesv column solve, and
##   G^{(1,1)}(t) = 1/(abar a_t) (1/Lt) sum_k e^{ikt} [D_ov(k)^{-1}]_{(src,up),(src,up)}
## with src a 5-fold vertex (V.2).  The raw Wilson propagator (Fig. 10) is the same
## sum with [D_W(k)^{-1}] instead.
##
## Gauge (periodic, k = 2 pi n/Lt).  Per mode the quadratic form (IV.24) is the
## (N_E + N_V) x (N_E + N_V) Hermitian matrix M(k) = sum_p beta_p conj(w_p) w_p^T over
## the plaquette rows w: faces are the incidence rows (spatial only) and the temporal
## plaquette of edge a->b is [(1 - e^{ik}) at the edge, +1 at site b, -1 at site a].
## ker M(k) is the momentum-k gauge orbit (plus the flat temporal mode at k = 0); it is
## regularized away with sigma G(k) G(k)^dag (+ tau u u^T/N_V at k = 0), the same
## "add G G^dag" trick as ops/gaugeact.RegOp, exact for the gauge-invariant face
## sources.  Then, averaged over all source triangles f (weights 1/A_f^2, (V.12)-(V.13)),
##   G_g(t) = 1/(g^2 N_F) (1/Lt) sum_k e^{ikt} sum_f [c_f^dag M(k)^+ c_f]/A_f^2 ,
## using F(-k) = F(k) so only k in [0, pi] is computed.
##
## Both pipelines are validated in-app before the campaign runs: the fermion against
## the exact dense D_ov (ops/overlap.denseOv) and dense D_W on a full L=1, nt=6
## lattice (the Wilson check pins the sign of k in e^{ikt}); the gauge against the
## real-space regSolve pseudo-inverse (ops/gaugeact, WP-G's validated route) and
## against WP-G's pinned exact-area L=1, Lt=120, T=16 correlator table.
##
## All couplings are the tree-wide exact-kite-area convention (doc/06 "THE COUPLING
## CONVENTION"): core/lattice.kap and newBeta's default gcExactArea.
##
## The mode-block helpers at the top are exported for rspec.nim (the spectra
## figures); the campaign itself runs under `when isMainModule`.

import base
import std/[math, os, strformat, strutils, tables, times]
import eigens/lapack
import eigens/linalgFuncs
import core/analytic
import meas/[dataio, fit, observables]

# --- dense helpers (shared with rspec) ---------------------------------------

proc zsolve*(a: var seq[Complex64], b: var seq[Complex64], n, nrhs: int) =
  ## Solve a x = b with LAPACK zgesv.  `a` (column-major n x n) is overwritten by
  ## its LU factors, `b` (n x nrhs) by the solution.
  var
    nn = fint(n)
    nr = fint(nrhs)
    ipiv = newSeq[fint](n)
    info = fint(0)
  zgesv(addr nn, addr nr, cast[ptr dcomplex](addr a[0]), addr nn,
        addr ipiv[0], cast[ptr dcomplex](addr b[0]), addr nn, addr info)
  doAssert info == 0, "zgesv: info = " & $info

proc zmm(ta, tb: cstring, n: int, a, b: seq[Complex64], c: var seq[Complex64]) =
  ## c = op(a) op(b) for square column-major n x n matrices (BLAS zgemm).
  var
    nn = fint(n)
    one = dcomplex(re: 1.0, im: 0.0)
    zero = dcomplex(re: 0.0, im: 0.0)
  zgemm(ta, tb, addr nn, addr nn, addr nn, addr one,
        cast[ptr dcomplex](unsafeAddr a[0]), addr nn,
        cast[ptr dcomplex](unsafeAddr b[0]), addr nn,
        addr zero, cast[ptr dcomplex](addr c[0]), addr nn)

proc eigvals*(a: var seq[Complex64], nd: int): seq[Complex64] =
  ## General complex eigenvalues via zgeev; `a` is destroyed.
  result = newSeq[Complex64](nd)
  zgeigs(cast[ptr float64](addr a[0]), cast[ptr float64](addr result[0]), nd)

# --- fermion Matsubara blocks -------------------------------------------------

proc spatialDense*(sph: Sphere, at: float, flatKap = false):
    tuple[l: Lat, a: seq[Complex64]] =
  ## The U = 0 spatial Wilson operator, 2 N_V x 2 N_V column-major, row index
  ## 2*site + spin.  nt = 1 because the spatial part does not see the time
  ## direction; the returned Lat still carries kappa'(at) for the mode blocks.
  ## `flatKap` overrides kappa with the flat identity (l*_1 + l*_2)/abar of (IV.2)
  ## (the slide-8-legend convention, doc/06); the default is the exact kite area
  ## that core/lattice now builds in.
  let l = newLat(sph, 1, at)
  if flatKap:
    for e in 0..<sph.ne: l.kap[e] = sph.edges[e].dual/sph.abar
  (l, denseDw(l, newGauge(l), 0.0, dwSpace))

proc modeBlockInto*(l: Lat, dsp: seq[Complex64], k: float,
                    x: var seq[Complex64]) =
  ## x = D_W(k) = D_spatial + kappa'_y (1 - cos k + i sigma_3 sin k), i.e.
  ## kappa'_y (1 - e^{-ik}) on spin up and kappa'_y (1 - e^{+ik}) on spin down
  ## (the construction twilson validates against the full dense operator).
  let
    nv = l.sph.nv
    nd = 2*nv
  if x.len != dsp.len: x = newSeq[Complex64](dsp.len)
  copyMem(addr x[0], unsafeAddr dsp[0], dsp.len*sizeof(Complex64))
  for y in 0..<nv:
    let
      kt = l.kapT[y]
      dr = kt*(1.0 - cos(k))
      di = kt*sin(k)
    x[(2*y) + nd*(2*y)] += complex64(dr, di)
    x[(2*y + 1) + nd*(2*y + 1)] += complex64(dr, -di)

proc modeBlock*(l: Lat, dsp: seq[Complex64], k: float): seq[Complex64] =
  modeBlockInto(l, dsp, k, result)

proc ovFromDwInto*(dw: var seq[Complex64], nd: int, bigM: float,
                   h, w, g: var seq[Complex64], ev: var seq[float]) =
  ## dw (= D_W(k) on entry) is replaced by D_ov(k) = 1 + X (X^dag X)^{-1/2},
  ## X = D_W(k) - bigM, via zheev on H = X^dag X.  h, w, g, ev are scratch.
  if h.len != nd*nd: h = newSeq[Complex64](nd*nd)
  if w.len != nd*nd: w = newSeq[Complex64](nd*nd)
  if g.len != nd*nd: g = newSeq[Complex64](nd*nd)
  if ev.len != nd: ev = newSeq[float](nd)
  for i in 0..<nd: dw[i + nd*i] -= complex64(bigM, 0.0)      # dw = X
  zmm("C", "N", nd, dw, dw, h)                               # h = X^dag X
  zeigs(cast[ptr float64](addr h[0]), addr ev[0], nd)        # h <- eigenvectors V
  for j in 0..<nd:
    let s = 1.0/sqrt(ev[j])
    for i in 0..<nd: w[i + nd*j] = s*h[i + nd*j]             # w = V E^{-1/2}
  zmm("N", "C", nd, w, h, g)                                 # g = V E^{-1/2} V^dag
  zmm("N", "N", nd, dw, g, w)                                # w = X g
  copyMem(addr dw[0], addr w[0], nd*nd*sizeof(Complex64))
  for i in 0..<nd: dw[i + nd*i] += complex64(1.0, 0.0)       # dw = 1 + X g

proc ovFromDw*(dw: seq[Complex64], nd: int, bigM = 1.0): seq[Complex64] =
  ## D_ov(k) from D_W(k) (fresh copy, convenience for rspec).
  result = dw
  var h, w, g: seq[Complex64]
  var ev: seq[float]
  ovFromDwInto(result, nd, bigM, h, w, g, ev)

type FreeProp* = object
  gov*: seq[float]     ## overlap G^(1,1)(t), t = i at, one period
  gw*: seq[float]      ## raw Wilson G^(1,1)(t) (empty unless requested)
  imOv*, imW*: float   ## max |Im G| / max |Re G|, a reality check on the mode sum

proc ovPropModes*(sph: Sphere, ltime: int, at: float, src: int,
                  wilson = false, bigM = 1.0): FreeProp =
  ## Free overlap (and optionally raw Wilson) temporal propagator by the
  ## antiperiodic Matsubara decomposition (module header).  Normalization (V.2):
  ## G_lat = [D^{-1}]/(abar a_t).
  let
    (l, dsp) = spatialDense(sph, at)
    nd = 2*sph.nv
    s2 = 2*src
  var
    accOv = newSeq[Complex64](ltime)
    accW = newSeq[Complex64](ltime)
    x, h, w, g: seq[Complex64]
    ev: seq[float]
    aw = newSeq[Complex64](nd*nd)
    rhs = newSeq[Complex64](nd)
  for n in 0..<ltime:
    let k = PI*float(2*n + 1)/float(ltime)
    modeBlockInto(l, dsp, k, x)
    if wilson:
      copyMem(addr aw[0], addr x[0], nd*nd*sizeof(Complex64))
      for i in 0..<nd: rhs[i] = complex64(0.0, 0.0)
      rhs[s2] = complex64(1.0, 0.0)
      zsolve(aw, rhs, nd, 1)
      let z = rhs[s2]
      for t in 0..<ltime:
        accW[t] += complex64(cos(k*float(t)), sin(k*float(t)))*z
    ovFromDwInto(x, nd, bigM, h, w, g, ev)                   # x <- D_ov(k)
    for i in 0..<nd: rhs[i] = complex64(0.0, 0.0)
    rhs[s2] = complex64(1.0, 0.0)
    zsolve(x, rhs, nd, 1)
    let z = rhs[s2]
    for t in 0..<ltime:
      accOv[t] += complex64(cos(k*float(t)), sin(k*float(t)))*z
  let sc = 1.0/(sph.abar*at*float(ltime))
  result.gov = newSeq[float](ltime)
  var remax, immax = 0.0
  for t in 0..<ltime:
    result.gov[t] = sc*accOv[t].re
    remax = max(remax, abs(accOv[t].re))
    immax = max(immax, abs(accOv[t].im))
  result.imOv = immax/remax
  if wilson:
    result.gw = newSeq[float](ltime)
    remax = 0.0
    immax = 0.0
    for t in 0..<ltime:
      result.gw[t] = sc*accW[t].re
      remax = max(remax, abs(accW[t].re))
      immax = max(immax, abs(accW[t].im))
    result.imW = immax/remax

# --- gauge Matsubara blocks ---------------------------------------------------

proc gaugeModeInto*(sph: Sphere, bt: Beta, k, sig, tau: float,
                    a: var seq[Complex64]) =
  ## a = M(k) + sig G(k) G(k)^dag (+ tau u u^T/N_V at k = 0): the per-mode gauge
  ## quadratic form, kernel-regularized (module header).  Dimension N_E + N_V,
  ## spatial edges first.  M(k)[p,q] = sum_rows beta conj(w_p) w_q, so
  ## M(-k) = conj(M(k)) and the mode sum needs only k in [0, pi].
  let
    ne = sph.ne
    n = ne + sph.nv
    ek = complex64(cos(k), sin(k))
    onec = complex64(1.0, 0.0)
  if a.len != n*n: a = newSeq[Complex64](n*n)
  else:
    for i in 0..<a.len: a[i] = complex64(0.0, 0.0)
  for f in 0..<sph.nf:
    let
      fc = sph.faces[f]
      bf = bt.face[f]
    for i in 0..2:
      for j in 0..2:
        a[fc.e[i] + n*fc.e[j]] += complex64(bf*float(fc.s[i]*fc.s[j]), 0.0)
  for e in 0..<ne:
    let
      be = bt.edge[e]
      idx = [e, ne + sph.edges[e].b, ne + sph.edges[e].a]
      w = [onec - ek, onec, -onec]
    for p in 0..2:
      for q in 0..2:
        a[idx[p] + n*idx[q]] += be*conjugate(w[p])*w[q]
  # gauge-orbit regularization: column g_y is the gauge shift of lambda = delta_y,
  # +-1 on the incident spatial edges and (e^{ik} - 1) on the temporal link at y
  let gt = ek - onec
  for y in 0..<sph.nv:
    let deg = sph.nbe[y].len
    var
      idx = newSeq[int](deg + 1)
      v = newSeq[Complex64](deg + 1)
    for j in 0..<deg:
      let e = sph.nbe[y][j]
      idx[j] = e
      v[j] = complex64(if sph.edges[e].b == y: 1.0 else: -1.0, 0.0)
    idx[deg] = ne + y
    v[deg] = gt
    for p in 0..deg:
      for q in 0..deg:
        a[idx[p] + n*idx[q]] += sig*v[p]*conjugate(v[q])
  if abs(1.0 - ek.re) + abs(ek.im) < 1e-14:                  # k = 0: flat mode
    let tv = tau/float(sph.nv)
    for y1 in 0..<sph.nv:
      for y2 in 0..<sph.nv:
        a[(ne + y1) + n*(ne + y2)] += complex64(tv, 0.0)

type GaugeCorr* = object
  g*: seq[float]       ## G_g(t), t = i at, one period
  imax*: float         ## max |Im c^dag A^{-1} c| over modes and sources
  secs*: float

proc gaugeCorrModes*(lv, ltime: int, at, g2: float): GaugeCorr =
  ## G_g(t) = (1/g^2) <J^t(t) J^t(0)> averaged over all source triangles, from the
  ## per-mode pseudo-inverse; couplings newBeta (gcExactArea, the paper's).
  let
    t0 = epochTime()
    sph = newSphere(lv)
    lat = newLat(sph, ltime, at)
    bt = newBeta(lat, g2)
    ne = sph.ne
    n = ne + sph.nv
    nf = sph.nf
  var s = 0.0
  for x in bt.face: s += x
  for x in bt.edge: s += x
  s /= float(bt.face.len + bt.edge.len)          # newRegOp's default scale
  var
    a: seq[Complex64]
    b = newSeq[Complex64](n*nf)
  result.g = newSeq[float](ltime)
  doAssert ltime mod 2 == 0
  for m in 0..(ltime div 2):
    let k = 2.0*PI*float(m)/float(ltime)
    gaugeModeInto(sph, bt, k, s, s, a)
    for i in 0..<b.len: b[i] = complex64(0.0, 0.0)
    for f in 0..<nf:
      let fc = sph.faces[f]
      for i in 0..2: b[fc.e[i] + n*f] += complex64(float(fc.s[i]), 0.0)
    zsolve(a, b, n, nf)
    var fk = 0.0
    for f in 0..<nf:
      let fc = sph.faces[f]
      var z = complex64(0.0, 0.0)
      for i in 0..2: z += float(fc.s[i])*b[fc.e[i] + n*f]
      fk += z.re/(bt.afac[f]*bt.afac[f])
      result.imax = max(result.imax, abs(z.im))
    let wt = if m == 0 or 2*m == ltime: 1.0 else: 2.0
    for t in 0..<ltime:
      result.g[t] += wt*cos(k*float(t))*fk
  let sc = 1.0/(g2*float(nf)*float(ltime))
  for t in 0..<ltime: result.g[t] *= sc
  result.secs = epochTime() - t0

# ==============================================================================
#  the campaign
# ==============================================================================

when isMainModule:
  qexInit()
  freezeTimers()

  letParam:
    fLev = 8           ## fermion sector: L = 1, 2, 4, ... up to this level
    gLev = 8           ## gauge sector likewise
    doCheck = 1        ## run the dense / real-space oracle self-checks first
    doFermion = 1
    doGauge = 1
    nmaxScan = 0       ## echo res(nmax) curves for the fermion fit variants
    g2inv = 20.0       ## gauge 1/g^2 (paper: 20; G_g is exactly g-independent)
    nmaxHi = 60        ## upper end of the n_max search range
    fitLo = 4.0        ## (V.6) plateau window [fitLo, fitHi), the paper's 4 <= t < 8
    fitHi = 8.0

  installStandardParams()
  echoParams()
  processHelpParam()

  const outDir = currentSourcePath().parentDir.parentDir.parentDir.parentDir /
                 "output" / "radial" / "free"
  createDir outDir

  # Per-(L, Lt, T) correlator cache: every mode-sum result is written as soon as
  # it is computed and reloaded on a rerun, so an interrupted campaign resumes at
  # the point where it stopped (the L = 8 fermion points cost ~15 min each).
  # Delete output/radial/free/cache/ to force a full recompute.
  const cacheDir = outDir / "cache"
  createDir cacheDir

  proc cachedProp(sph: Sphere, lv, lt: int, tt, at: float, src: int,
                  wilson: bool): tuple[p: FreeProp, secs: float,
                                       cached: bool] =
    let path = cacheDir / &"f_L{lv}_Lt{lt}_T{int(tt)}.tsv"
    if fileExists(path):
      let r = readTsv(path)
      if r.cols.len >= 1 and r.cols[0].len == lt and
         (not wilson or r.cols.len >= 2):
        result.p.gov = r.cols[0]
        if r.cols.len >= 2: result.p.gw = r.cols[1]
        result.p.imOv = parseFloat(r.meta.getOrDefault("imOv", "0"))
        result.p.imW = parseFloat(r.meta.getOrDefault("imW", "0"))
        result.cached = true
        return
    let t0 = epochTime()
    result.p = ovPropModes(sph, lt, at, src, wilson)
    result.secs = epochTime() - t0
    if wilson:
      writeTsv(path, {"imOv": &"{result.p.imOv:.17g}",
                      "imW": &"{result.p.imW:.17g}", "src": $src,
                      "secs": &"{result.secs:.1f}"},
               ["gov", "gw"], [result.p.gov, result.p.gw])
    else:
      writeTsv(path, {"imOv": &"{result.p.imOv:.17g}", "src": $src,
                      "secs": &"{result.secs:.1f}"},
               ["gov"], [result.p.gov])

  proc cachedGauge(lv, lt: int, tt, at, g2: float):
      tuple[c: GaugeCorr, cached: bool] =
    let path = cacheDir / &"g_L{lv}_Lt{lt}_T{int(tt)}.tsv"
    if fileExists(path):
      let r = readTsv(path)
      if r.cols.len >= 1 and r.cols[0].len == lt:
        result.c.g = r.cols[0]
        result.c.imax = parseFloat(r.meta.getOrDefault("imax", "0"))
        result.cached = true
        return
    result.c = gaugeCorrModes(lv, lt, at, g2)
    writeTsv(path, {"imax": &"{result.c.imax:.17g}", "g2": &"{g2:.17g}",
                    "secs": &"{result.c.secs:.1f}"},
             ["gg"], [result.c.g])

  proc levels(top: int): seq[int] =
    var l = 1
    while l <= top:
      result.add l
      l *= 2

  proc levIdx(lv: int): int =
    ## Index of L in (1, 2, 4, 8) for the published-value tables.
    case lv
    of 1: 0
    of 2: 1
    of 4: 2
    else: 3

  # published values (doc/02 section 8)
  const
    pubF1 = 0.953918                ## fermion Delta_0, L=1, Lt=168, T=16
    pubFc = "0.999998(34)"          ## fermion Delta_0^cont
    pubFn = [6, 10, 19, 32]         ## fermion n_max, L=1,2,4,8
    pubFr = [0.028, 0.012, 0.0039, 0.038]
    pubG1 = 1.33242                 ## gauge Delta_0, L=1, Lt=120, T=16
    pubGc = "1.41409(18)"           ## gauge Delta_0^cont
    pubGn = [3, 8, 18, 35]          ## gauge n_max, L=1,2,4,8
    pubGr = [0.0031, 0.0023, 0.0031, 0.0037]
    # WP-G's pinned exact-area L=1, Lt=120, T=16 correlator (doc/06 T1.5a/T1.5b)
    refIdx = [0, 4, 8, 15, 23, 30, 38, 45]
    refG = [1.93333729, 4.66030750e-1, 1.30684453e-1, 1.99691311e-2,
            3.40566260e-3, 8.59960736e-4, 1.94526577e-4, 5.49935238e-5]

  var master: seq[array[3, string]]    # quantity | published | ours (+dev)

  proc addRow(q, pub, ours: string) = master.add [q, pub, ours]

  # ---------------------------------------------------------------------------
  # self-checks: the mode pipelines against exact dense / real-space oracles
  # ---------------------------------------------------------------------------

  if doCheck != 0:
    echo ""
    echo "================ self-checks (dense and real-space oracles) ================"
    block:                                   # fermion: L=1, nt=6, at=0.4, dense
      let
        sph = newSphere(1)
        nt = 6
        at = 0.4
        l = newLat(sph, nt, at)
        u0 = newGauge(l)
        src = fiveFoldSite(sph)
        nd = 2*l.nsite
        sc = 1.0/(sph.abar*at)
        p = ovPropModes(sph, nt, at, src, wilson = true)
      # dense overlap via ops/overlap.denseOv (independent code path)
      let o = newOv(l, 1.0, newRat(1.0, 50.0, 11), 1e-24, 1e-20, 10000)
      var dov = denseOv(o, u0)
      var rhs = newSeq[Complex64](nd)
      rhs[2*sIdx(l, src, 0)] = complex64(1.0, 0.0)
      zsolve(dov, rhs, nd, 1)
      var eOv = 0.0
      var gs = 0.0
      for t in 0..<nt: gs = max(gs, abs(p.gov[t]))
      for t in 0..<nt:
        eOv = max(eOv, abs(p.gov[t] - sc*rhs[2*sIdx(l, src, t)].re)/gs)
      # dense raw Wilson: pins the sign of k in e^{ikt} (D_W is T-asymmetric)
      var dw = denseDw(l, u0, 0.0)
      for i in 0..<nd: rhs[i] = complex64(0.0, 0.0)
      rhs[2*sIdx(l, src, 0)] = complex64(1.0, 0.0)
      zsolve(dw, rhs, nd, 1)
      var eW = 0.0
      gs = 0.0
      for t in 0..<nt: gs = max(gs, abs(p.gw[t]))
      for t in 0..<nt:
        eW = max(eW, abs(p.gw[t] - sc*rhs[2*sIdx(l, src, t)].re)/gs)
      echo &"  fermion modes vs dense D_ov (L=1, nt=6):  max rel dev = {eOv:.3e}" &
           &"   (Im/Re of the mode sum: {p.imOv:.2e})"
      echo &"  fermion modes vs dense D_W  (L=1, nt=6):  max rel dev = {eW:.3e}" &
           &"   (Im/Re: {p.imW:.2e})"
      doAssert eOv < 1e-9, "fermion overlap mode pipeline disagrees with denseOv"
      doAssert eW < 1e-9, "fermion Wilson mode pipeline disagrees with denseDw"
    block:                                   # gauge: L=1, nt=16, real-space regSolve
      let
        lv = 1
        nt = 16
        at = 0.25
        g2 = 1.0/g2inv
        sph = newSphere(lv)
        lat = newLat(sph, nt, at)
        bt = newBeta(lat, g2)
        c = gaugeCorrModes(lv, nt, at, g2)
      var reg = newRegOp(lat, bt)
      var
        b = newGauge(lat)
        x = newGauge(lat)
        gr = newSeq[float](nt)
      for f in 0..<sph.nf:
        triSource(lat, b, f, 0)
        let info = regSolve(lat, x, b, reg, 1e-26, 400000)
        doAssert info.converged
        let ia = 1.0/(bt.afac[f]*bt.afac[f])
        for t in 0..<nt: gr[t] += ia*plaqSpatial(lat, x, f, t)
      let sc = 1.0/(g2*float(sph.nf))
      var e = 0.0
      for t in 0..<nt:
        gr[t] *= sc
        e = max(e, abs(c.g[t]/gr[t] - 1.0))
      echo &"  gauge modes vs real-space regSolve (L=1, nt=16): max rel dev = {e:.3e}" &
           &"   (contraction Im max: {c.imax:.2e})"
      doAssert e < 1e-8, "gauge mode pipeline disagrees with the regSolve oracle"

  # ---------------------------------------------------------------------------
  # A. fermion sector
  # ---------------------------------------------------------------------------

  type FRow = object
    lv, lt: int
    at, d0, ed0, chi: float

  proc delta0Fit(g: seq[float], at, T: float): PlateauFit =
    ## effMass + (V.6) plateau fit over [fitLo, fitHi), rgauge's estimator.
    let half = int(round(0.5*T/at))
    var gh = newSeq[float](half + 1)
    for i in 0..half: gh[i] = g[i]
    let m = effMass(gh, at, T)
    result = plateauFit(m, int(round(fitLo/at)), min(int(round(fitHi/at)), m.len), at)
    if result.status != fitOk:
      raise newException(ValueError, "delta0Fit: " & $result.status)

  proc normDev(g: seq[float], at, T: float, lo, hi: float,
               model: proc(t: float): float): float =
    ## max relative deviation of g from the analytic model over t in [lo, hi].
    for i in 0..<g.len:
      let t = at*float(i)
      if t >= lo and t <= hi:
        result = max(result, abs(g[i]/model(t) - 1.0))

  var
    fGrid: seq[FRow]
    fNorm: seq[(int, float)]         # (L, Fig-7 normalization dev)
    fNmax: seq[(int, int, float, int, int, float, float)]
      # (L, nmax_rel, res/dof_rel, nmax_abs, nmax_log, res/dof_abs, C_rel)
    f1d0: PlateauFit
    tFermion = 0.0

  if doFermion != 0:
    let tf0 = epochTime()
    echo ""
    echo "================ A. fermion sector =========================================="
    const
      lt12 = 168
      t12 = 12.0
      t16 = 16.0
      ntGrid = [96, 120, 144, 168]
    let at12 = t12/float(lt12)
    for lv in levels(fLev):
      let
        sph = newSphere(lv)
        src = fiveFoldSite(sph)
        tl0 = epochTime()
      # --- T = 12, Lt = 168: Fig 7, Fig 10 (L=1), n_max fits ---
      let (p12, secs12, cached12) = cachedProp(sph, lv, lt12, t12, at12, src,
                                               wilson = (lv == 1))
      block:
        var cols = @[newSeq[float](0), newSeq[float](0), newSeq[float](0),
                     newSeq[float](0), newSeq[float](0)]
        for i in 0..<lt12:
          let t = at12*float(i)
          cols[0].add t
          cols[1].add p12.gov[i]
          cols[2].add fermionGPeriodic(t, t12)
          cols[3].add (if i > 0: fermionG(t) else: NaN)
          cols[4].add (if i > 0: p12.gov[i]/fermionGPeriodic(t, t12) else: NaN)
        writeTsv(outDir / &"fig7_L{lv}.tsv",
                 {"source": "arXiv:2510.03085 Fig. 7", "lattice": &"L{lv}",
                  "Lt": $lt12, "T": $t12, "at": &"{at12:.17g}",
                  "abar": &"{sph.abar:.17g}", "src": $src,
                  "normalization": "G_lat = D_ov^{-1}/(abar at), (V.2)"},
                 ["t", "G_lat", "G_periodic", "G_inf", "ratio"], cols)
      let nd = normDev(p12.gov, at12, t12, 2.0, 6.0,
                       proc(t: float): float = fermionGPeriodic(t, t12))
      fNorm.add (lv, nd)
      echo &"--- L = {lv}: Fig 7 (T=12, Lt=168)  max |G_lat/G_analytic - 1| over" &
           &" t in [2,6] = {nd:.3e}   (mode-sum Im/Re {p12.imOv:.1e}, " &
           &"{secs12:.1f} s{(if cached12: \", cached\" else: \"\")})"
      if lv == 1:
        var cols = @[newSeq[float](0), newSeq[float](0), newSeq[float](0)]
        for i in 0..<lt12:
          cols[0].add at12*float(i)
          cols[1].add p12.gov[i]
          cols[2].add p12.gw[i]
        writeTsv(outDir/"fig10_L1.tsv",
                 {"source": "arXiv:2510.03085 Fig. 10", "lattice": "L1",
                  "Lt": $lt12, "T": $t12, "at": &"{at12:.17g}",
                  "note": "overlap preserves the t -> T-t fold, raw Wilson breaks it"},
                 ["t", "G_overlap", "G_wilson"], cols)
        var vOv = 0.0
        var vW = 0.0
        var gsO = 0.0
        var gsW = 0.0
        for i in 1..<lt12:
          gsO = max(gsO, abs(p12.gov[i]))
          gsW = max(gsW, abs(p12.gw[i]))
        for i in 1..<lt12:
          vOv = max(vOv, abs(p12.gov[i] - p12.gov[lt12 - i])/gsO)
          vW = max(vW, abs(p12.gw[i] - p12.gw[lt12 - i])/gsW)
        echo &"    Fig 10: fold violation max|G(t)-G(T-t)|/max|G|:  overlap" &
             &" {vOv:.3e}   Wilson {vW:.3e}   (ratio {vW/vOv:.1e})"
        addRow("Fig10 fold violation (Wilson / overlap)", "T-violation visible",
               &"{vW:.2e} / {vOv:.2e}")
      # --- n_max fits (T1.4f/g), T = 12 ---
      # The t = 0 point of the (1,1) OVERLAP propagator is pure Ginsparg-Wilson
      # contact: (IV.17) forces Re[D_ov^{-1}]_{xx} = 1/2 exactly per diagonal
      # entry, and the T-odd physical part vanishes at t = 0, so
      # G(0) abar a_t = 1/2 to roundoff (measured and echoed below).  It carries
      # no spectral information, so the fit uses the Lt - 1 points t > 0
      # (fitting it raw drags the minimum down; "subtracting the contact" leaves
      # an exact zero, which a relative weight cannot hold).  PRIMARY convention:
      # (V.9) with C free and the relative weight w = 1/g^2 (nmaxFit's default,
      # WP-J).  Absolute and log weightings are echoed alongside, and so is the
      # relative residual evaluated AT the published n_max -- see the WP-K status
      # entry: the published (n_max, residual/DOF) pairs lie on our relative
      # residual curve, but they are not its minimum.
      block:
        var
          gn0 = p12.gov[1..^1]
          tsn0 = newSeq[float](lt12 - 1)
        for i in 1..<lt12: tsn0[i-1] = at12*float(i)
        let modeln0 = proc(nmax: int): seq[float] =
          result = newSeq[float](tsn0.len)
          for i in 0..<tsn0.len: result[i] = fermionGPeriodic(tsn0[i], t12, nmax)
        proc logFit(g: seq[float], md: proc(nmax: int): seq[float]):
            tuple[nmax: int, resDof: float] =
          var rlog = Inf
          for x in g:
            if x <= 0.0: return (0, NaN)
          for nm in 1..nmaxHi:
            let m = md(nm)
            var s = 0.0
            for i in 0..<g.len: s += ln(g[i]) - ln(m[i])
            let lc = s/float(g.len)
            var r = 0.0
            for i in 0..<g.len:
              let d = ln(g[i]) - ln(m[i]) - lc
              r += d*d
            if r < rlog:
              rlog = r
              result.nmax = nm
          result.resDof = rlog/float(g.len - 2)
        proc relAt(g: seq[float], md: proc(nmax: int): seq[float],
                   nm: int): float =
          ## relative res/dof of (V.9) at a fixed truncation
          let m = md(nm)
          var aa, bb = 0.0
          for i in 0..<g.len:
            let w = 1.0/(g[i]*g[i])
            aa += w*g[i]*m[i]
            bb += w*m[i]*m[i]
          let cc = aa/bb
          var r = 0.0
          for i in 0..<g.len:
            let d = (g[i] - cc*m[i])/g[i]
            r += d*d
          r/float(g.len - 2)
        var ones = newSeq[float](gn0.len)
        for i in 0..<gn0.len: ones[i] = 1.0
        let
          ip = levIdx(lv)
          fr = nmaxFit(gn0, modeln0, 1..nmaxHi)               # relative (primary)
          fa = nmaxFit(gn0, modeln0, 1..nmaxHi, ones)         # absolute
          lg = logFit(gn0, modeln0)                           # log
          rpub = relAt(gn0, modeln0, pubFn[ip])
        fNmax.add (lv, fr.nmax, fr.resDof, fa.nmax, lg.nmax, rpub, fr.c)
        echo &"    T1.4f/g n_max (T=12, {gn0.len} points t>0, dof {fr.dof}):"
        echo &"      relative {fr.nmax}  res/dof {fr.resDof:.3e} (C {fr.c:.4f})" &
             &"   absolute {fa.nmax} ({fa.resDof:.2e})  log {lg.nmax}" &
             &" ({lg.resDof:.2e})"
        echo &"      published {pubFn[ip]} (res/dof {pubFr[ip]});  our relative" &
             &" res/dof AT n_max = {pubFn[ip]}: {rpub:.4f}"
        echo &"      GW contact purity: G(0) abar at - 1/2 = " &
             &"{p12.gov[0]*sph.abar*at12 - 0.5:.3e}"
        if nmaxScan != 0:
          # residual curves under every (data window, model, weighting) variant,
          # to identify the paper's convention rather than guess it
          var ts = newSeq[float](lt12)
          for i in 0..<lt12: ts[i] = at12*float(i)
          let model = proc(nmax: int): seq[float] =
            result = newSeq[float](ts.len)
            for i in 0..<ts.len: result[i] = fermionGPeriodic(ts[i], t12, nmax)
          var half = newSeq[float](0)          # t in [at, T/2]
          var tsh = newSeq[float](0)
          for i in 1..int(round(0.5*t12/at12)):
            half.add p12.gov[i]
            tsh.add ts[i]
          let modelh = proc(nmax: int): seq[float] =
            result = newSeq[float](tsh.len)
            for i in 0..<tsh.len: result[i] = fermionGPeriodic(tsh[i], t12, nmax)
          let modelhi = proc(nmax: int): seq[float] =
            result = newSeq[float](tsh.len)
            for i in 0..<tsh.len: result[i] = fermionG(tsh[i], nmax)
          proc curve(tag: string, g: seq[float],
                     md: proc(nmax: int): seq[float], e: seq[float]) =
            var line = &"      scan {tag}:"
            for nm in 1..min(nmaxHi, 40):
              let m = md(nm)
              var aa, bb = 0.0
              for i in 0..<g.len:
                let w = if e.len == 0: 1.0/(g[i]*g[i]) else: 1.0/(e[i]*e[i])
                aa += w*g[i]*m[i]
                bb += w*m[i]*m[i]
              let cc = aa/bb
              var r = 0.0
              for i in 0..<g.len:
                let w = if e.len == 0: 1.0/(g[i]*g[i]) else: 1.0/(e[i]*e[i])
                let d = g[i] - cc*m[i]
                r += w*d*d
              line.add &" {r/float(g.len - 2):.2e}"
            echo line
          var eh: seq[float]
          curve("full rel      ", gn0, modeln0, eh)
          var eabs = newSeq[float](gn0.len)
          for i in 0..<gn0.len: eabs[i] = 1.0
          curve("full abs      ", gn0, modeln0, eabs)
          curve("half rel      ", half, modelh, eh)
          var ehabs = newSeq[float](half.len)
          for i in 0..<half.len: ehabs[i] = 1.0
          curve("half abs      ", half, modelh, ehabs)
          curve("half rel inf  ", half, modelhi, eh)
          var esq = newSeq[float](half.len)     # w = 1/g -> e = sqrt(g)
          for i in 0..<half.len: esq[i] = sqrt(half[i])
          curve("half w=1/g    ", half, modelh, esq)
          proc logCurve(tag: string, g: seq[float],
                        md: proc(nmax: int): seq[float]) =
            var line = &"      scan {tag}:"
            for nm in 1..min(nmaxHi, 40):
              let m = md(nm)
              var s = 0.0
              for i in 0..<g.len: s += ln(g[i]) - ln(m[i])
              let lc = s/float(g.len)
              var r = 0.0
              for i in 0..<g.len:
                let d = ln(g[i]) - ln(m[i]) - lc
                r += d*d
              line.add &" {r/float(g.len - 2):.2e}"
            echo line
          logCurve("full log      ", gn0, modeln0)
          logCurve("half log      ", half, modelh)
          logCurve("raw  log      ", p12.gov, model)
          # C fixed to 1 (no normalization fit), relative residual
          proc c1Curve(tag: string, g: seq[float],
                       md: proc(nmax: int): seq[float]) =
            var line = &"      scan {tag}:"
            for nm in 1..min(nmaxHi, 40):
              let m = md(nm)
              var r = 0.0
              for i in 0..<g.len:
                let d = (g[i] - m[i])/g[i]
                r += d*d
              line.add &" {r/float(g.len - 2):.2e}"
            echo line
          c1Curve("full rel C=1  ", gn0, modeln0)
          # Delta_eff-based selection: no normalization freedom at all
          block:
            let nh = int(round(0.5*t12/at12))
            var gh = newSeq[float](nh + 1)
            for i in 0..nh: gh[i] = p12.gov[i]
            let md = effMass(gh, at12, t12)
            var l1 = &"      scan dEff abs      :"
            var l2 = &"      scan dEff rel      :"
            for nm in 1..min(nmaxHi, 40):
              var ga = newSeq[float](nh + 1)
              for i in 0..nh: ga[i] = fermionGPeriodic(at12*float(max(i, 1)), t12, nm)
              ga[0] = ga[1]
              let mm = effMass(ga, at12, t12)
              var ra, rr = 0.0
              for i in 1..<md.len:
                let d = md[i] - mm[i]
                ra += d*d
                rr += d*d/(md[i]*md[i])
              l1.add &" {ra/float(md.len - 2):.2e}"
              l2.add &" {rr/float(md.len - 2):.2e}"
            echo l1
            echo l2
      # --- T = 16 grid: Fig 8, T1.4c, T1.4e ---
      for lt in ntGrid:
        let
          at = t16/float(lt)
          (p, _, _) = cachedProp(sph, lv, lt, t16, at, src, false)
          fit = delta0Fit(p.gov, at, t16)
        fGrid.add FRow(lv: lv, lt: lt, at: at, d0: fit.d0, ed0: fit.ed0,
                       chi: fit.chi2dof)
        echo &"    T=16 Lt={lt:3d}: Delta_0 = {fit.d0:.6f} +- {fit.ed0:.6f}" &
             &"   chi2/dof = {fit.chi2dof:.2e}  (status {fit.status})"
        if lv == 1 and lt == 168: f1d0 = fit
        if lt == 168:                                  # Fig 8 at Lt = 168
          let half = int(round(0.5*t16/at))
          var gh = newSeq[float](half + 1)
          for i in 0..half: gh[i] = p.gov[i]
          let m = effMass(gh, at, t16)
          var ga = newSeq[float](half + 1)
          for i in 0..half: ga[i] = fermionGPeriodic(at*float(max(i, 1)), t16)
          ga[0] = ga[1]                                # contact: index 0 unused
          let ma = effMass(ga, at, t16)
          var cols = @[newSeq[float](0), newSeq[float](0), newSeq[float](0)]
          for i in 1..<m.len:
            cols[0].add at*float(i)
            cols[1].add m[i]
            cols[2].add ma[i]
          writeTsv(outDir / &"fig8_L{lv}.tsv",
                   {"source": "arXiv:2510.03085 Fig. 8", "lattice": &"L{lv}",
                    "Lt": $lt, "T": $t16, "at": &"{at:.17g}"},
                   ["t", "Delta_eff", "Delta_eff_analytic"], cols)
      echo &"    L = {lv} done in {epochTime() - tl0:.1f} s"
    # T1.4c
    echo ""
    echo &"--- T1.4c  fermion Delta_0 (L=1, Lt=168, T=16, window [{fitLo},{fitHi})):"
    echo &"    ours = {f1d0.d0:.6f} +- {f1d0.ed0:.6f}   published = {pubF1}" &
         &"   dev = {f1d0.d0 - pubF1:+.6f} ({100.0*(f1d0.d0/pubF1 - 1.0):+.4f}%)"
    addRow("fermion Delta_0(L=1, Lt=168, T=16)", $pubF1,
           &"{f1d0.d0:.6f}  (dev {f1d0.d0 - pubF1:+.1e})")
    # T1.4e: (V.7) over the paper's grid
    if fLev >= 4:
      proc planeF(minL: int, lts: openArray[int]): PlaneFit =
        var as2, at2, y, e: seq[float]
        for r in fGrid:
          if r.lv >= minL and r.lt in lts:
            let sph = newSphere(r.lv)
            as2.add sph.abar*sph.abar
            at2.add r.at*r.at
            y.add r.d0
            e.add 1.0
        contFit2(as2, at2, y, e)
      let
        pc = planeF(2, [120, 144, 168])
        ps = planeF(2, [96, 120, 144, 168])
        scC = sqrt(pc.chi2/float(pc.dof))
        scS = sqrt(ps.chi2/float(ps.dof))
      echo ""
      echo &"--- T1.4e  (V.7) fermion Delta_0^cont, central grid L>=2 x " &
           &"Lt in (120,144,168):"
      echo &"    Delta_0^cont = {pc.a:.6f} +- {pc.ea*scC:.6f}   c_s = {pc.cs:+.5f}" &
           &"   c_t = {pc.ct:+.5f}   rms resid = {scC:.2e}   (dof {pc.dof})"
      echo &"    systematic (+ Lt=96): {ps.a:.6f} +- {ps.ea*scS:.6f}" &
           &"   shift = {ps.a - pc.a:+.6f}"
      echo &"    published = {pubFc}   exact = 1"
      addRow("fermion Delta_0^cont (V.7)", pubFc,
             &"{pc.a:.6f} +- {pc.ea*scC:.6f}  [syst {ps.a - pc.a:+.6f}]")
      # Fig 9 data: the grid with both projections onto the central fit
      var cl, clt, cas2, cat2, cd0, ce, cps, cpt: seq[float]
      for r in fGrid:
        let sph = newSphere(r.lv)
        cl.add float(r.lv)
        clt.add float(r.lt)
        cas2.add sph.abar*sph.abar
        cat2.add r.at*r.at
        cd0.add r.d0
        ce.add r.ed0
        cps.add r.d0 - pc.ct*r.at*r.at             # vs abar^2 at a_t^2 removed
        cpt.add r.d0 - pc.cs*sph.abar*sph.abar     # vs a_t^2 at abar^2 removed
      writeTsv(outDir/"fig9_fermion_scaling.tsv",
               {"source": "arXiv:2510.03085 Fig. 9", "T": $t16,
                "fit": "central grid L in (2,4,8) x Lt in (120,144,168)",
                "Delta0cont": &"{pc.a:.17g}", "err": &"{pc.ea*scC:.17g}",
                "cs": &"{pc.cs:.17g}", "ct": &"{pc.ct:.17g}",
                "published": pubFc},
               ["L", "Lt", "abar2", "at2", "Delta0", "eDelta0",
                "proj_s", "proj_t"],
               [cl, clt, cas2, cat2, cd0, ce, cps, cpt])
    tFermion = epochTime() - tf0
    echo &"    fermion sector total: {tFermion:.1f} s"

  # ---------------------------------------------------------------------------
  # B. gauge sector
  # ---------------------------------------------------------------------------

  var
    gGrid: seq[FRow]
    gNorm: seq[(int, float)]
    gNmax: seq[(int, int, float, int, int, float, float)]
    g1d0: PlateauFit
    tGauge = 0.0

  if doGauge != 0:
    let tg0 = epochTime()
    echo ""
    echo "================ B. gauge sector ============================================"
    const
      lt12 = 120
      t12 = 12.0
      t16 = 16.0
      ntGrid = [48, 64, 96, 120]
    let
      g2 = 1.0/g2inv
      at12 = t12/float(lt12)
    for lv in levels(gLev):
      let
        sph = newSphere(lv)
        tl0 = epochTime()
      # --- T = 12, Lt = 120: Fig 11 and the n_max fits ---
      let (c12, _) = cachedGauge(lv, lt12, t12, at12, g2)
      block:
        let half = lt12 div 2
        var gh = newSeq[float](half + 1)
        for i in 0..half: gh[i] = c12.g[i]
        let m = effMass(gh, at12, t12)
        var cols = @[newSeq[float](0), newSeq[float](0), newSeq[float](0),
                     newSeq[float](0)]
        for i in 0..<lt12:
          let t = at12*float(i)
          cols[0].add t
          cols[1].add c12.g[i]
          cols[2].add gaugeGPeriodic(t, t12)
          cols[3].add (if i < m.len: m[i] else: NaN)
        writeTsv(outDir / &"fig11_L{lv}.tsv",
                 {"source": "arXiv:2510.03085 Fig. 11", "lattice": &"L{lv}",
                  "Lt": $lt12, "T": $t12, "at": &"{at12:.17g}",
                  "g2inv": &"{g2inv:.17g}"},
                 ["t", "G_g", "G_analytic", "Delta_eff"], cols)
      let nd = normDev(c12.g, at12, t12, 2.0, 6.0,
                       proc(t: float): float = gaugeGPeriodic(t, t12))
      gNorm.add (lv, nd)
      echo &"--- L = {lv}: Fig 11 (T=12, Lt=120)  max |G_g/G_analytic - 1| over" &
           &" t in [2,6] = {nd:.3e}   ({c12.secs:.1f} s, Im max {c12.imax:.1e})"
      block:                                     # T1.5e/f
        # The gauge correlator has no contact subtlety (the lattice t = 0 value
        # and the truncated model are both finite and physical), so all Lt
        # points enter, DOF = Lt - 2, matching the paper.  Weightings as in the
        # fermion fit; relative is primary.
        var ts = newSeq[float](lt12)
        for i in 0..<lt12: ts[i] = at12*float(i)
        let model = proc(nmax: int): seq[float] =
          result = newSeq[float](ts.len)
          for i in 0..<ts.len: result[i] = gaugeGPeriodic(ts[i], t12, nmax)
        let fr = nmaxFit(c12.g, model, 1..nmaxHi)
        var ones = newSeq[float](lt12)
        for i in 0..<lt12: ones[i] = 1.0
        let fa = nmaxFit(c12.g, model, 1..nmaxHi, ones)
        var nlog = 0
        var rlog = Inf
        for nm in 1..nmaxHi:
          let m = model(nm)
          var s = 0.0
          for i in 0..<lt12: s += ln(c12.g[i]) - ln(m[i])
          let lc = s/float(lt12)
          var r = 0.0
          for i in 0..<lt12:
            let d = ln(c12.g[i]) - ln(m[i]) - lc
            r += d*d
          if r < rlog:
            rlog = r
            nlog = nm
        let ip = levIdx(lv)
        var rpub: float
        block:
          let m = model(pubGn[ip])
          var aa, bb = 0.0
          for i in 0..<lt12:
            let w = 1.0/(c12.g[i]*c12.g[i])
            aa += w*c12.g[i]*m[i]
            bb += w*m[i]*m[i]
          let cc = aa/bb
          var r = 0.0
          for i in 0..<lt12:
            let d = (c12.g[i] - cc*m[i])/c12.g[i]
            r += d*d
          rpub = r/float(lt12 - 2)
        gNmax.add (lv, fr.nmax, fr.resDof, fa.nmax, nlog, rpub, fr.c)
        echo &"    T1.5e/f n_max (T=12, {lt12} points, dof {fr.dof}):  " &
             &"relative {fr.nmax}  res/dof {fr.resDof:.3e} (C {fr.c:.4f})   " &
             &"absolute {fa.nmax}  log {nlog}"
        echo &"      published {pubGn[ip]} (res/dof {pubGr[ip]});  our relative" &
             &" res/dof AT n_max = {pubGn[ip]}: {rpub:.4f}"
      # --- T = 16 grid: T1.5b, T1.5d, Fig 12 ---
      for lt in ntGrid:
        let
          at = t16/float(lt)
          (c, _) = cachedGauge(lv, lt, t16, at, g2)
          fit = delta0Fit(c.g, at, t16)
        gGrid.add FRow(lv: lv, lt: lt, at: at, d0: fit.d0, ed0: fit.ed0,
                       chi: fit.chi2dof)
        echo &"    T=16 Lt={lt:3d}: Delta_0 = {fit.d0:.6f} +- {fit.ed0:.6f}" &
             &"   chi2/dof = {fit.chi2dof:.2e}   ({c.secs:.1f} s)"
        if lv == 1 and lt == 120:
          g1d0 = fit
          var e = 0.0
          for i in 0..<refIdx.len:
            e = max(e, abs(c.g[refIdx[i]]/refG[i] - 1.0))
          echo &"    WP-G pinned exact-area table (8 t-values): max rel dev = {e:.3e}"
          doAssert e < 5e-7, "gauge pipeline disagrees with WP-G's pinned correlator"
      echo &"    L = {lv} done in {epochTime() - tl0:.1f} s"
    # T1.5b
    echo ""
    echo &"--- T1.5b  gauge Delta_0 (L=1, Lt=120, T=16, window [{fitLo},{fitHi})):"
    echo &"    ours = {g1d0.d0:.6f} +- {g1d0.ed0:.6f}   published = {pubG1}" &
         &"   dev = {g1d0.d0 - pubG1:+.6f} ({100.0*(g1d0.d0/pubG1 - 1.0):+.4f}%)"
    addRow("gauge Delta_0(L=1, Lt=120, T=16)", $pubG1,
           &"{g1d0.d0:.6f}  (dev {g1d0.d0 - pubG1:+.1e})")
    # T1.5d
    if gLev >= 4:
      proc planeG(minL, maxL: int): PlaneFit =
        var as2, at2, y, e: seq[float]
        for r in gGrid:
          if r.lv >= minL and r.lv <= maxL:
            let sph = newSphere(r.lv)
            as2.add sph.abar*sph.abar
            at2.add r.at*r.at
            y.add r.d0
            e.add 1.0
        contFit2(as2, at2, y, e)
      let
        pc = planeG(2, 8)
        scC = sqrt(pc.chi2/float(pc.dof))
      echo ""
      echo &"--- T1.5d  (V.7) gauge Delta_0^cont, paper grid L in (2..{gLev}) x " &
           &"Lt in (48,64,96,120):"
      echo &"    Delta_0^cont = {pc.a:.6f} +- {pc.ea*scC:.6f}   c_s = {pc.cs:+.5f}" &
           &"   c_t = {pc.ct:+.5f}   rms resid = {scC:.2e}   (dof {pc.dof})"
      echo &"    published = {pubGc}   exact sqrt(2) = {sqrt(2.0):.6f}" &
           &"   dev = {100.0*(pc.a/sqrt(2.0) - 1.0):+.4f}%"
      if gLev >= 8:                               # WP-G's grid, for comparison
        let
          pw = planeG(1, 4)
          scW = sqrt(pw.chi2/float(pw.dof))
        echo &"    (WP-G's L in (1,2,4) grid from this pipeline: " &
             &"{pw.a:.6f} +- {pw.ea*scW:.6f}; WP-G measured 1.414535(286))"
      addRow("gauge Delta_0^cont (V.7)", pubGc & ", exact 1.414214",
             &"{pc.a:.6f} +- {pc.ea*scC:.6f}")
      var cl, clt, cas2, cat2, cd0, ce, cps, cpt: seq[float]
      for r in gGrid:
        let sph = newSphere(r.lv)
        cl.add float(r.lv)
        clt.add float(r.lt)
        cas2.add sph.abar*sph.abar
        cat2.add r.at*r.at
        cd0.add r.d0
        ce.add r.ed0
        cps.add r.d0 - pc.ct*r.at*r.at
        cpt.add r.d0 - pc.cs*sph.abar*sph.abar
      writeTsv(outDir/"fig12_gauge_scaling.tsv",
               {"source": "arXiv:2510.03085 Fig. 12", "T": $t16,
                "fit": &"paper grid L in (2..{gLev}) x Lt in (48,64,96,120)",
                "Delta0cont": &"{pc.a:.17g}", "err": &"{pc.ea*scC:.17g}",
                "cs": &"{pc.cs:.17g}", "ct": &"{pc.ct:.17g}",
                "published": pubGc},
               ["L", "Lt", "abar2", "at2", "Delta0", "eDelta0",
                "proj_s", "proj_t"],
               [cl, clt, cas2, cat2, cd0, ce, cps, cpt])
    tGauge = epochTime() - tg0
    echo &"    gauge sector total: {tGauge:.1f} s"

  # ---------------------------------------------------------------------------
  # master table
  # ---------------------------------------------------------------------------

  for (lv, nd) in fNorm:
    addRow(&"Fig 7 normalization, L={lv} (max rel dev, t in [2,6])",
           "agreement incl. normalization", &"{nd:.3e}")
  for r in fNmax:
    let ip = levIdx(r[0])
    addRow(&"fermion n_max(L={r[0]})", $pubFn[ip],
           &"{r[1]} (rel; abs {r[3]}, log {r[4]}); res/dof {r[2]:.2e};" &
           &" res@pub {r[5]:.4f} [pub {pubFr[ip]}]")
  for (lv, nd) in gNorm:
    addRow(&"Fig 11 normalization, L={lv} (max rel dev, t in [2,6])",
           "agreement incl. normalization", &"{nd:.3e}")
  for r in gNmax:
    let ip = levIdx(r[0])
    addRow(&"gauge n_max(L={r[0]})", $pubGn[ip],
           &"{r[1]} (rel; abs {r[3]}, log {r[4]}); res/dof {r[2]:.2e};" &
           &" res@pub {r[5]:.4f} [pub {pubGr[ip]}]")
  addRow("wall clock fermion / gauge",
         "-", &"{tFermion:.0f} s / {tGauge:.0f} s")

  echo ""
  echo "================ MASTER TABLE (quantity | published | ours) ================="
  for r in master:
    echo &"  {r[0]:<52} | {r[1]:<28} | {r[2]}"
  echo ""
  echo &"wrote {outDir}"

  processSaveParams()
  writeParamFile()
  qexFinalize()
