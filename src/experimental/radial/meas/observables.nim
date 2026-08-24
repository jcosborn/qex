## Measurements on S^2 x R  (WP-I): fermion propagator, condensate, conserved
## link currents, scalar correlators, gluonic operators and the Wilson-loop
## GEVP basis.
##
## Normative reference: doc/07-observables.md (all sections),
## doc/04-interfaces.md section 14.
##
## THE single current kernel (doc/07 section 1.1): the conserved current is
## J_l = dS_F/dtheta_l = Psibar K_l Psi with K_l = dD_ov/dtheta_l, and the only
## implementation of K_l in the tree is `ovGradient` (WP-F).  Every current
## measurement here reaches K_l through `linkCurrent` (two ovGradient calls),
## and the connected correlator is estimated in the factorized form
##   tr[K_A S K_B S] = E_{eta,xi} <eta, K_A S xi> <xi, K_B S eta>
## (eta, xi independent unit-covariance Gaussian noises), which needs ONLY
## sequential solves and pullbacks -- no separately coded tangent, no adjoint
## tangent.  The dense oracle `denseOvDeriv` (tests only) is the same rational
## formula evaluated in exact linear algebra and is pinned against ovGradient.
##
## Propagator conventions: S(m) = (D_ov + m)^{-1} (additive mass),
##   S b       = (D^dag D)^{-1} D^dag b        -- `propSolve` (adjoint FIRST)
##   S^dag b   = D (D^dag D)^{-1} b            -- `propSolveDag`
## Getting the order wrong gives D^dag (D^dag D)^{-1} = D^{-1} D^{-dag} D^dag
## only for normal D -- the test pins `propSolve` against the dense inverse.
##
## Dense helpers are column-major, dimension 2*nsite, row index 2*sIdx + spin,
## exactly like `denseDw`/`denseOv`, and are for tests/small lattices only.

import std/[math, complex]
import ../core/lattice
import ../core/spinor
import ../ops/overlap
import ../ops/gaugeact
import harmonics
import eigens/linalgFuncs

export overlap, gaugeact, harmonics

# --- small utilities ---------------------------------------------------------

func fiveFoldSite*(sph: Sphere): int =
  ## The first 5-fold (original icosahedron) site: present at every refinement
  ## level, so continuum extrapolations compare the same physical point.
  for y in 0..<sph.nv:
    if sph.nbr[y].len == 5: return y
  doAssert false, "no 5-fold site found"

proc propSolve*(o: Ov, x: var Spin, b: Spin, u: Gauge, mass = 0.0): CgInfo =
  ## x = S b = (D^dag D)^{-1} D^dag b, D = D_ov + mass.
  var rhs = newSpin(o.l.nsite)
  applyOvAdj(o, rhs, b, u, mass)
  solveNormal(o, x, rhs, u, mass)

proc propSolveDag*(o: Ov, x: var Spin, b: Spin, u: Gauge, mass = 0.0): CgInfo =
  ## x = S^dag b = D (D^dag D)^{-1} b.
  var y = newSpin(o.l.nsite)
  result = solveNormal(o, y, b, u, mass)
  applyOv(o, x, y, u, mass)

# --- fermion propagator ------------------------------------------------------

proc propagatorT*(o: Ov, u: Gauge, mass: float, src: int): seq[Spinor] =
  ## G(t) at the coincident spatial site: the solution of D(mass) x = b for a
  ## unit source in spinor component 0 at (src, t = 0), read back at (src, t).
  ## `src` should be a 5-fold site (`fiveFoldSite`).
  let b = pointSource(o.l.nsite, sIdx(o.l, src, 0), 0)
  var x = newSpin(o.l.nsite)
  let ci = propSolve(o, x, b, u, mass)
  doAssert ci.converged, "propagatorT: solve failed"
  result = newSeq[Spinor](o.l.nt)
  for t in 0..<o.l.nt: result[t] = x[sIdx(o.l, src, t)]

# --- chiral condensate -------------------------------------------------------

proc condensatePS*(o: Ov, u: Gauge, mass: float, nnoise: int,
                   r: var Threefry4x64): tuple[v, e: float] =
  ## GW contact-subtracted condensate (doc/07 section 3):
  ##   Sigma = Re tr[(1 - D_ov/2)(D_ov + mass)^{-1}] / (2 nsite)
  ## by Gaussian volume noise, <eta eta^dag> = 1 per complex component.
  ## Returns the noise mean and the standard error of the mean.
  let n = o.l.nsite
  var
    eta = newSpin(n)
    y = newSpin(n)
    z = newSpin(n)
    s = 0.0
    s2 = 0.0
  for k in 0..<nnoise:
    eta.gaussian r
    let ci = propSolve(o, y, eta, u, mass)   # y = (D + m)^{-1} eta
    doAssert ci.converged, "condensatePS: solve failed"
    applyOv(o, z, y, u, 0.0)                 # z = D_ov y (massless D_ov)
    let v = (redot(eta, y) - 0.5*redot(eta, z))/float(2*n)
    s += v
    s2 += v*v
  result.v = s/float(nnoise)
  let varm = max(0.0, s2/float(nnoise) - result.v*result.v)
  result.e = sqrt(varm/float(max(1, nnoise - 1)))

proc condensateDense*(o: Ov, u: Gauge, mass: float): float =
  ## Exact: (1/N) sum_k Re[(1 - lambda_k/2)/(lambda_k + mass)] over the dense
  ## D_ov spectrum, N = 2 nsite.  The trace of the same rational function of
  ## D_ov, so it equals `condensatePS`'s expectation up to the rational/solver
  ## error of the iterative path.
  let
    nd = 2*o.l.nsite
    a = denseOv(o, u)
  var m = a
  var ev = newSeq[Complex64](nd)
  zgeigs(cast[ptr float64](addr m[0]), cast[ptr float64](addr ev[0]), nd)
  for k in 0..<nd:
    let z = (1.0 - 0.5*ev[k])/(ev[k] + mass)
    result += z.re
  result /= float(nd)

# --- dense oracles (tests only) ----------------------------------------------

proc zmm(a, b: seq[Complex64], n: int): seq[Complex64] =
  ## c = a b, column-major.
  result = newSeq[Complex64](n*n)
  for j in 0..<n:
    for k in 0..<n:
      let bkj = b[k + n*j]
      if bkj.re != 0.0 or bkj.im != 0.0:
        for i in 0..<n:
          result[i + n*j] += a[i + n*k]*bkj

proc zmmAdjL(a, b: seq[Complex64], n: int): seq[Complex64] =
  ## c = a^dag b, column-major.
  result = newSeq[Complex64](n*n)
  for j in 0..<n:
    for i in 0..<n:
      var s = complex64(0.0, 0.0)
      for k in 0..<n: s += conjugate(a[k + n*i])*b[k + n*j]
      result[i + n*j] = s

proc denseS*(o: Ov, u: Gauge, mass = 0.0): seq[Complex64] =
  ## Exact dense S = (D_ov + mass)^{-1} via A^{-1} = (A^dag A)^{-1} A^dag with
  ## the Hermitian A^dag A eigendecomposed by zheev.  Tests only.
  let nd = 2*o.l.nsite
  var a = denseOv(o, u)
  for i in 0..<nd: a[i + nd*i] += complex64(mass, 0.0)
  var h = zmmAdjL(a, a, nd)
  var ev = newSeq[float](nd)
  zeigs(cast[ptr float64](addr h[0]), addr ev[0], nd)   # h <- eigenvectors V
  let w = zmm(a, h, nd)                                 # w = A V
  # A^{-1} = V diag(1/ev) (A V)^dag
  result = newSeq[Complex64](nd*nd)
  for j in 0..<nd:
    for k in 0..<nd:
      let f = conjugate(w[j + nd*k])/ev[k]
      if f.re != 0.0 or f.im != 0.0:
        for i in 0..<nd:
          result[i + nd*j] += h[i + nd*k]*f

proc denseDwDeriv(l: Lat, u, du: Gauge): seq[Complex64] =
  ## Dense delta D_W[du], column by column through applyDwDeriv.
  let
    n = l.nsite
    nd = 2*n
  result = newSeq[Complex64](nd*nd)
  var
    b = newSpin(n)
    col = newSpin(n)
  for j in 0..<nd:
    b.zero
    b[j shr 1][j and 1] = complex64(1.0, 0.0)
    applyDwDeriv(l, col, b, u, du)
    for i in 0..<nd: result[i + nd*j] = col[i shr 1][i and 1]

proc denseOvDeriv*(o: Ov, u: Gauge, du: Gauge): seq[Complex64] =
  ## Exact dense tangent delta D_ov[du] of the RATIONAL overlap operator --
  ## the same formula ovGradient pulls back (doc/04 section 10), evaluated in
  ## exact linear algebra:
  ##   delta D_ov = dX R(H) + X dR,  dR = -sum_j r_j G_j dH G_j,
  ##   dH = dX^dag X + X^dag dX,     G_j = (H + q_j)^{-1},
  ## with H = X^dag X eigendecomposed.  Tests only; pinned against ovGradient
  ## by contraction in tmeas.
  let
    l = o.l
    nd = 2*l.nsite
    x = denseDw(l, u, o.m)
    dx = denseDwDeriv(l, u, du)
  var v = zmmAdjL(x, x, nd)              # H
  var ev = newSeq[float](nd)
  zeigs(cast[ptr float64](addr v[0]), addr ev[0], nd)   # v <- eigenvectors
  let
    xv = zmm(x, v, nd)
    dxv = zmm(dx, v, nd)
    p = zmmAdjL(dxv, xv, nd)             # (dX V)^dag (X V)
  # wt = V^dag dH V = p + p^dag ; drt_{kl} = -wt_{kl} sum_j r_j/((ev_k+q_j)(ev_l+q_j))
  var drt = newSeq[Complex64](nd*nd)
  for j in 0..<nd:
    for i in 0..<nd:
      var s = 0.0
      for q in 0..<o.rat.npole:
        s += o.rat.res[q]/((ev[i] + o.rat.pole[q])*(ev[j] + o.rat.pole[q]))
      drt[i + nd*j] = -s*(p[i + nd*j] + conjugate(p[j + nd*i]))
  # rh = V diag(R(ev)) V^dag
  var rh = newSeq[Complex64](nd*nd)
  for j in 0..<nd:
    for k in 0..<nd:
      let f = ratValue(o.rat, ev[k])*conjugate(v[j + nd*k])
      for i in 0..<nd:
        rh[i + nd*j] += v[i + nd*k]*f
  # dr = V drt V^dag
  let q1 = zmm(v, drt, nd)
  var dr = newSeq[Complex64](nd*nd)
  for j in 0..<nd:
    for k in 0..<nd:
      let f = conjugate(v[j + nd*k])
      for i in 0..<nd:
        dr[i + nd*j] += q1[i + nd*k]*f
  result = zmm(dx, rh, nd)
  let t2 = zmm(x, dr, nd)
  for i in 0..<nd*nd: result[i] += t2[i]

# --- the current kernel, measurement face ------------------------------------

proc tsliceForm*(l: Lat, w: openArray[float], t: int): Gauge =
  ## The one-form du with du^t_{(y,t)} = w_y on the temporal links of slice t
  ## and zero elsewhere: the direction whose tangent delta D_ov[du] is the
  ## w-weighted temporal current insertion at slice t.
  result = newGauge(l)
  for y in 0..<l.sph.nv:
    result.t[tIdx(l, y, t)] = w[y]

proc linkCurrent*(o: Ov, u: Gauge, left, right: Spin): tuple[re, im: Gauge] =
  ## <left, K_l right> per link, K_l = dD_ov/dtheta_l, through the ONE kernel:
  ## ovGradient gives 2 Re<left, K_l right>; a second call with i*left gives
  ## 2 Im.  result.re.s[l] + i*result.im.s[l] = <left, K_l right> (same for .t).
  let n = o.l.nsite
  result.re = newGauge(o.l)
  result.im = newGauge(o.l)
  ovGradient(o, result.re, left, right, u, 0.5)
  var il = newSpin(n)
  for i in 0..<n:
    for c in 0..1:
      il[i][c] = complex64(-left[i][c].im, left[i][c].re)
  ovGradient(o, result.im, il, right, u, 0.5)

proc tsliceAmp*(l: Lat, g: tuple[re, im: Gauge], w: openArray[float]):
    seq[Complex64] =
  ## A(t) = sum_y w_y <left, K_{(y,t)} right> from a `linkCurrent` result:
  ## the temporal-link current amplitude per time slice with site weights w
  ## (e.g. w_y = Y_lm(pos_y); doc/07 1.3 -- the dual area lives inside J_link,
  ## so the projection weight is Y_lm alone and the l = 0 weight w = 1 is the
  ## exactly conserved total charge).
  result = newSeq[Complex64](l.nt)
  for t in 0..<l.nt:
    var s = complex64(0.0, 0.0)
    for y in 0..<l.sph.nv:
      let i = tIdx(l, y, t)
      s += w[y]*complex64(g.re.t[i], g.im.t[i])
    result[t] = s

proc wardChargeScan*(o: Ov, u: Gauge, mass: float, v0, ta, tb: int):
    tuple[c: seq[Complex64], jump: Complex64] =
  ## The exact measurement-level Ward/charge test (doc/07 section 1.1).
  ## With x = S e_a (source at (v0, ta, spin 0)) and r = S^dag e_b (sink at
  ## (v0, tb, spin 0)), C(t) = sum_y <r, K_{(y,t)} x> = (S Q(t) S)_{ba} is the
  ## total-charge insertion in a fermion line.  Q(t) - Q(t') = i(Lam D - D Lam)
  ## for the slice-band indicator Lam (a pure gauge variation), so
  ##   C(t) - C(t+1) = i S_{ba} (delta_{t+1,ta} - delta_{t+1,tb}):
  ## C is EXACTLY constant on each arc between ta and tb and jumps by
  ## -+ i S_{ba} at them -- the fermion line carries unit charge.  `jump`
  ## returns i S_{ba} = i x[(v0,tb)][0], the predicted jump.  Holds at any
  ## mass (it is gauge variance, not GW).
  let
    l = o.l
    ba = pointSource(l.nsite, sIdx(l, v0, ta), 0)
    bb = pointSource(l.nsite, sIdx(l, v0, tb), 0)
  var
    x = newSpin(l.nsite)
    rr = newSpin(l.nsite)
  var ci = propSolve(o, x, ba, u, mass)
  doAssert ci.converged
  ci = propSolveDag(o, rr, bb, u, mass)
  doAssert ci.converged
  let g = linkCurrent(o, u, rr, x)
  var w = newSeq[float](l.sph.nv)
  for y in 0..<w.len: w[y] = 1.0
  result.c = tsliceAmp(l, g, w)
  let s = x[sIdx(l, v0, tb)][0]
  result.jump = complex64(-s.im, s.re)          # i S_{ba}

# --- connected and disconnected current correlators ---------------------------

type CurrentSample* = object
  ## One factorized sample of the connected trace and one-noise traces:
  ##   a[k] = <eta, K_{(slice k)} S xi>,  b[k] = <xi, K_{(slice k)} S eta>,
  ##   d[k] = <eta, K_{(slice k)} S eta>
  ## flattened over k = iop*nt + t for the requested site-weight vectors.
  ## E[a[k1] b[k2]] = tr[K_{k1} S K_{k2} S]   (the connected contraction),
  ## E[d[k]]        = tr[K_{k} S]             (the disconnected trace).
  a*, b*, d*: seq[Complex64]

proc currentSample*(o: Ov, u: Gauge, mass: float, w: openArray[seq[float]],
                    r: var Threefry4x64): CurrentSample =
  ## One noise pair (eta, xi): 2 sequential solves + 3 linkCurrent calls.
  ## Estimators are unbiased for independent unit-covariance complex noise.
  let
    l = o.l
    n = l.nsite
    nt = l.nt
  var
    eta = newSpin(n)
    xi = newSpin(n)
    xeta = newSpin(n)
    xxi = newSpin(n)
  eta.gaussian r
  xi.gaussian r
  var ci = propSolve(o, xeta, eta, u, mass)
  doAssert ci.converged
  ci = propSolve(o, xxi, xi, u, mass)
  doAssert ci.converged
  let
    ga = linkCurrent(o, u, eta, xxi)    # <eta, K S xi>
    gb = linkCurrent(o, u, xi, xeta)    # <xi, K S eta>
    gd = linkCurrent(o, u, eta, xeta)   # <eta, K S eta>
  result.a = newSeq[Complex64](w.len*nt)
  result.b = newSeq[Complex64](w.len*nt)
  result.d = newSeq[Complex64](w.len*nt)
  for iop in 0..<w.len:
    let
      aa = tsliceAmp(l, ga, w[iop])
      bb = tsliceAmp(l, gb, w[iop])
      dd = tsliceAmp(l, gd, w[iop])
    for t in 0..<nt:
      result.a[iop*nt + t] = aa[t]
      result.b[iop*nt + t] = bb[t]
      result.d[iop*nt + t] = dd[t]

proc currentCorrConn*(samples: openArray[CurrentSample], k1, k2: int):
    tuple[v, e: Complex64] =
  ## Mean and standard error of tr[K_{k1} S K_{k2} S] over the samples,
  ## from the factorized products a[k1]*b[k2].  The physical connected
  ## correlator of doc/07 1.2 per 4-component flavor is
  ##   C_conn = -2 Re tr[K_A S K_B S].
  var
    s = complex64(0.0, 0.0)
    s2r = 0.0
    s2i = 0.0
  for sm in samples:
    let p = sm.a[k1]*sm.b[k2]
    s += p
    s2r += p.re*p.re
    s2i += p.im*p.im
  let n = float(samples.len)
  result.v = s/n
  let
    vr = max(0.0, s2r/n - result.v.re*result.v.re)
    vi = max(0.0, s2i/n - result.v.im*result.v.im)
    m = float(max(1, samples.len - 1))
  result.e = complex64(sqrt(vr/m), sqrt(vi/m))

proc currentTraceDisc*(samples: openArray[CurrentSample], k1, k2: int):
    tuple[v, e: float] =
  ## Disconnected piece 2Re tr[K_{k1} S] * 2Re tr[K_{k2} S] from unbiased
  ## cross products of DIFFERENT noise samples (same-sample products are
  ## biased by the noise variance).  Keep nsamples modest; the noise on this
  ## is what the deck's slide 12 complains about.
  var
    s = 0.0
    s2 = 0.0
    n = 0
  for i in 0..<samples.len:
    for j in 0..<samples.len:
      if i == j: continue
      let p = 4.0*samples[i].d[k1].re*samples[j].d[k2].re
      s += p
      s2 += p*p
      inc n
  result.v = s/float(n)
  let varm = max(0.0, s2/float(n) - result.v*result.v)
  result.e = sqrt(varm/float(max(1, n - 1)))   # correlated pairs: optimistic

# --- scalar correlators (doc/07 section 3) ------------------------------------

proc scalarCorrDense*(o: Ov, u: Gauge, mass = 0.0): tuple[ps, fs: seq[float]] =
  ## Connected timeslice correlators of sigma_PS and sigma_FS from the dense
  ## propagator, averaged over the source slice (tests only).  With slice
  ## projectors P_t and S = (D_ov + m)^{-1}, Sb = S^dag:
  ##   C_PS(dt) = -(1/nt) sum_t ( tr[P_{t+dt} S P_t S] + tr[P_{t+dt} Sb P_t Sb] )
  ##   C_FS(dt) = -(1/nt) sum_t ( tr[P_{t+dt} S P_t S]
  ##                              + tr[P_{t+dt} (Sb-1) P_t (Sb-1)] )
  ## from sigma_FS = eta^dag xi - xi^dag (1 - D_ov^dag) eta and
  ## (1 - D_ov^dag) Sb = Sb - 1 at mass 0 (the GW contact subtraction).
  ## At mass 0 and dt != 0 the two are IDENTICAL configuration by
  ## configuration; at dt = 0 they differ by the contact term.
  let
    l = o.l
    nv = l.sph.nv
    nt = l.nt
    nd = 2*l.nsite
    s = denseS(o, u, mass)
  result.ps = newSeq[float](nt)
  result.fs = newSeq[float](nt)
  # site-diagonal spin blocks of S restricted to slice pairs
  template idx(v, t, c: int): int = 2*(v + nv*t) + c
  for dt in 0..<nt:
    var
      aps = 0.0
      afs = 0.0
    for t1 in 0..<nt:
      let t2 = (t1 + dt) mod nt
      var
        z1 = complex64(0.0, 0.0)   # tr[P2 S P1 S]
        z2 = complex64(0.0, 0.0)   # tr[P2 (Sdag-1) P1 (Sdag-1)]
      for x in 0..<nv:
        for y in 0..<nv:
          for c in 0..1:
            for cp in 0..1:
              let
                i2 = idx(x, t2, c)
                i1 = idx(y, t1, cp)
                sxy = s[i2 + nd*i1]                       # S_{x2, y1}
                syx = s[i1 + nd*i2]                       # S_{y1, x2}
              z1 += sxy*syx
              # (Sdag - 1)_{x2,y1} = conj(S_{y1,x2}) - delta
              var
                a = conjugate(syx)
                b = conjugate(sxy)
              if i2 == i1:
                a -= complex64(1.0, 0.0)
                b -= complex64(1.0, 0.0)
              z2 += a*b
      aps += -2.0*z1.re
      afs += -z1.re - z2.re
    result.ps[dt] = aps/float(nt)
    result.fs[dt] = afs/float(nt)

proc scalarCorrPoint*(o: Ov, u: Gauge, mass: float, v0, t0: int):
    tuple[ps, fs: seq[float]] =
  ## Source-resolved scalar correlators from a point at (v0, t0): the same
  ## contractions with P_{t1} replaced by the projector on the two spin
  ## components at (v0, t0).  4 solves: S e_c and S^dag e_c for c = 0, 1.
  ##   ps[dt] = -sum_{x in slice t0+dt} sum_{c,c'} 2 Re[ S_{xc,y0c'} S_{y0c',xc} ]
  ##   fs[dt] = same with the second factor pair from (S^dag - 1).
  ## Identical for dt != 0 at mass 0 (assert in tests); at dt = 0 they differ
  ## by the GW contact term.
  let
    l = o.l
    n = l.nsite
    nv = l.sph.nv
    nt = l.nt
    y0 = sIdx(l, v0, t0)
  var
    col: array[2, Spin]
    rw: array[2, Spin]
  for c in 0..1:
    col[c] = newSpin(n)
    rw[c] = newSpin(n)
    let b = pointSource(n, y0, c)
    var ci = propSolve(o, col[c], b, u, mass)
    doAssert ci.converged
    ci = propSolveDag(o, rw[c], b, u, mass)
    doAssert ci.converged
  result.ps = newSeq[float](nt)
  result.fs = newSeq[float](nt)
  for dt in 0..<nt:
    let t2 = (t0 + dt) mod nt
    var
      z1 = complex64(0.0, 0.0)
      z2 = complex64(0.0, 0.0)
    for x in 0..<nv:
      let i = sIdx(l, x, t2)
      for c in 0..1:
        for cp in 0..1:
          let
            sxy = col[cp][i][c]                # S_{xc, y0cp}
            syx = conjugate(rw[cp][i][c])      # S_{y0cp, xc}
          z1 += sxy*syx
          var
            a = conjugate(syx)                 # (Sdag)_{xc,y0cp}
            b = conjugate(sxy)                 # (Sdag)_{y0cp,xc}
          if i == y0 and c == cp:
            a -= complex64(1.0, 0.0)
            b -= complex64(1.0, 0.0)
          z2 += a*b
    result.ps[dt] = -2.0*z1.re
    result.fs[dt] = -z1.re - z2.re

# --- gluonic operators (doc/07 section 4) --------------------------------------

proc jtopProject*(l: Lat, u: Gauge, lh, mh: int): seq[float] =
  ## O^top_lm(t) = sum_f Y_lm(cc_f) Theta_f(t)  (the area cancels, doc/07 4.1).
  result = newSeq[float](l.nt)
  for t in 0..<l.nt:
    var s = 0.0
    for f in 0..<l.sph.nf:
      s += ylm(lh, mh, l.sph.faces[f].cc)*plaqSpatial(l, u, f, t)
    result[t] = s

proc f2Project*(l: Lat, u: Gauge, lh, mh: int): seq[float] =
  ## O^{F^2}_lm(t) = sum_f Y_lm(cc_f) Theta_f(t)^2 / A_f
  ##              + sum_e Y_lm(mid_e) Theta_e(t)^2 * 2 A_e / (l_e a_t)^2.
  ## RAW: the l = 0 channel carries <O>^2 -- subtract the ensemble mean
  ## (vacuum subtraction) at analysis time, per doc/07 4.1.
  result = newSeq[float](l.nt)
  for t in 0..<l.nt:
    var s = 0.0
    for f in 0..<l.sph.nf:
      let th = plaqSpatial(l, u, f, t)
      s += ylm(lh, mh, l.sph.faces[f].cc)*th*th/l.sph.faces[f].area
    for e in 0..<l.sph.ne:
      let
        ed = l.sph.edges[e]
        mid = unit(l.sph.pos[ed.a] + l.sph.pos[ed.b])
        th = plaqTemporal(l, u, e, t)
        den = ed.len*l.at
      s += ylm(lh, mh, mid)*th*th*2.0*ed.area/(den*den)
    result[t] = s

# --- Wilson-loop shapes (doc/07 4.2; see doc/06 WP-I for the change record) ----

type LoopShape* = enum
  ## The GEVP basis.  All shapes are exactly gauge-invariant flux sums and
  ## exactly I_h-covariant pseudoscalars, so the Y_lm projection is clean.
  ## Shapes 5-7 of doc/07 4.2 are realized as the temporal plaquettes/rectangles
  ## COMBINED round a face (or rhombus) with the face orientation signs -- the
  ## raw per-edge temporal plaquette carries the arbitrary canonical-orientation
  ## sign of the edge, which is not I_h-covariant and would poison the
  ## projection -- and in the time-reflection-EVEN second-difference form, so
  ## the correlator matrix is symmetric (a one-sided difference makes
  ## C_ij(dt) != C_ji(dt) and an indefinite symmetrized C(t0); measured, see
  ## doc/06 WP-I).  The combinations equal exact second time differences of
  ## the spatial fluxes -- a legitimate enlargement of the variational space.
  ##
  ## CAUTION (measured): at L = 1 the fixed-l projection collapses every
  ## spatial shape onto ONE operator (T1 has multiplicity 1 in the 20-face
  ## permutation rep), so the independent l = 1 basis at L = 1 is
  ## {lsTri, lsTPlaq, lsTRect2}; the full 7 needs L >= 2.
  lsTri        ## 1: elementary spatial triangle, Theta_f
  lsRhomb      ## 2: two triangles sharing an edge, Theta_f0 + Theta_f1
  lsStar       ## 3: the faces round a site, sum_{f in star(y)} Theta_f
  lsQuad       ## 4 (L>=2): face + its 3 edge neighbours (the L/2 triangle)
  lsTPlaq      ## 5: temporal plaquettes round a face, even difference:
               ##    Theta_f(t-1) - 2 Theta_f(t) + Theta_f(t+1)
  lsTRect2     ## 6: extent-2 temporal rectangles round a face:
               ##    Theta_f(t-2) - 2 Theta_f(t) + Theta_f(t+2)
  lsTRhomb     ## 7: temporal plaquettes round a rhombus (two-link spatial
               ##    sides), even difference of the rhombus flux

func loopCount*(l: Lat, sh: LoopShape): int =
  case sh
  of lsTri, lsQuad, lsTPlaq, lsTRect2: l.sph.nf
  of lsRhomb, lsTRhomb: l.sph.ne
  of lsStar: l.sph.nv

func loopCenter*(sph: Sphere, sh: LoopShape, i: int): Vec3 =
  case sh
  of lsTri, lsQuad, lsTPlaq, lsTRect2: sph.faces[i].cc
  of lsRhomb, lsTRhomb: unit(sph.pos[sph.edges[i].a] + sph.pos[sph.edges[i].b])
  of lsStar: sph.pos[i]

func otherFace(sph: Sphere, e, f: int): int =
  let ed = sph.edges[e]
  if ed.f[0] == f: ed.f[1] else: ed.f[0]

func loopFlux*(l: Lat, u: Gauge, sh: LoopShape, i, t: int): float =
  ## Theta_C for loop `i` of shape `sh` on slice t (temporal shapes span t..).
  case sh
  of lsTri:
    result = plaqSpatial(l, u, i, t)
  of lsRhomb:
    let ed = l.sph.edges[i]
    result = plaqSpatial(l, u, ed.f[0], t) + plaqSpatial(l, u, ed.f[1], t)
  of lsStar:
    for f in l.sph.nbf[i]: result += plaqSpatial(l, u, f, t)
  of lsQuad:
    result = plaqSpatial(l, u, i, t)
    for k in 0..2:
      result += plaqSpatial(l, u, otherFace(l.sph, l.sph.faces[i].e[k], i), t)
  of lsTPlaq:
    let fc = l.sph.faces[i]
    for k in 0..2:
      result += float(fc.s[k])*(plaqTemporal(l, u, fc.e[k], t-1) -
                                plaqTemporal(l, u, fc.e[k], t))
  of lsTRect2:
    let fc = l.sph.faces[i]
    for k in 0..2:
      result += float(fc.s[k])*(plaqTemporal(l, u, fc.e[k], t-2) +
                                plaqTemporal(l, u, fc.e[k], t-1) -
                                plaqTemporal(l, u, fc.e[k], t) -
                                plaqTemporal(l, u, fc.e[k], t+1))
  of lsTRhomb:
    let ed = l.sph.edges[i]
    for f in ed.f:
      let fc = l.sph.faces[f]
      for k in 0..2:
        result += float(fc.s[k])*(plaqTemporal(l, u, fc.e[k], t-1) -
                                  plaqTemporal(l, u, fc.e[k], t))

proc loopProject*(l: Lat, u: Gauge, sh: LoopShape, lh, mh: int): seq[float] =
  ## O_{sh,lm}(t) = sum_i Y_lm(center_i) Theta_C(i, t).
  result = newSeq[float](l.nt)
  let nc = loopCount(l, sh)
  for t in 0..<l.nt:
    var s = 0.0
    for i in 0..<nc:
      s += ylm(lh, mh, loopCenter(l.sph, sh, i))*loopFlux(l, u, sh, i, t)
    result[t] = s

proc loopOps*(l: Lat, u: Gauge, shapes: openArray[LoopShape],
              lh, mh: int): seq[seq[float]] =
  ## result[k][t] = the Y_lm-projected shape-k operator per time slice.
  result = newSeq[seq[float]](shapes.len)
  for k in 0..<shapes.len:
    result[k] = loopProject(l, u, shapes[k], lh, mh)

# incidence sources for the exact (deterministic) Gaussian correlators

proc addTriSource(l: Lat, b: var Gauge, f, t: int, w: float) =
  let fc = l.sph.faces[f]
  for k in 0..2:
    b.s[eIdx(l, fc.e[k], t)] += w*float(fc.s[k])

proc addTPlaqSource(l: Lat, b: var Gauge, e, t: int, w: float) =
  let ed = l.sph.edges[e]
  b.s[eIdx(l, e, t)] += w
  b.t[tIdx(l, ed.b, t)] += w
  b.s[eIdx(l, e, t+1)] -= w
  b.t[tIdx(l, ed.a, t)] -= w

proc loopSource*(l: Lat, b: var Gauge, sh: LoopShape, lh, mh, t: int) =
  ## b = the link-incidence vector of O_{sh,lm}(t), so that
  ## <O_{sh,lm}(t) O_{sh',l'm'}(t')> = b^T Mtilde^{-1} b'.  Always a
  ## combination of plaquette incidence rows, hence exactly transverse
  ## (in range(M)) -- regSolve applies.
  b.zero
  let nc = loopCount(l, sh)
  for i in 0..<nc:
    let w = ylm(lh, mh, loopCenter(l.sph, sh, i))
    case sh
    of lsTri:
      addTriSource(l, b, i, t, w)
    of lsRhomb:
      addTriSource(l, b, l.sph.edges[i].f[0], t, w)
      addTriSource(l, b, l.sph.edges[i].f[1], t, w)
    of lsStar:
      for f in l.sph.nbf[i]: addTriSource(l, b, f, t, w)
    of lsQuad:
      addTriSource(l, b, i, t, w)
      for k in 0..2:
        addTriSource(l, b, otherFace(l.sph, l.sph.faces[i].e[k], i), t, w)
    of lsTPlaq:
      let fc = l.sph.faces[i]
      for k in 0..2:
        addTPlaqSource(l, b, fc.e[k], t-1, w*float(fc.s[k]))
        addTPlaqSource(l, b, fc.e[k], t, -w*float(fc.s[k]))
    of lsTRect2:
      let fc = l.sph.faces[i]
      for k in 0..2:
        addTPlaqSource(l, b, fc.e[k], t-2, w*float(fc.s[k]))
        addTPlaqSource(l, b, fc.e[k], t-1, w*float(fc.s[k]))
        addTPlaqSource(l, b, fc.e[k], t, -w*float(fc.s[k]))
        addTPlaqSource(l, b, fc.e[k], t+1, -w*float(fc.s[k]))
    of lsTRhomb:
      let ed = l.sph.edges[i]
      for f in ed.f:
        let fc = l.sph.faces[f]
        for k in 0..2:
          addTPlaqSource(l, b, fc.e[k], t-1, w*float(fc.s[k]))
          addTPlaqSource(l, b, fc.e[k], t, -w*float(fc.s[k]))

proc jtopCorrExact*(l: Lat, bt: Beta, lh: int,
                    r2req = 1e-26, maxits = 200000): seq[seq[float]] =
  ## Exact free-theory correlator matrix of the J_top projections:
  ##   result[dt][m*(2l+1) + m'] = <O_lm(dt) O_lm'(0)> = b_lm(dt)^T M^+ b_lm'(0)
  ## via regSolve (the heatbath covariance is exactly M^+).  Deterministic --
  ## this is the machine-precision icosahedral-degeneracy object of doc/07
  ## section 2 and the oracle for the Monte-Carlo pipeline.
  let nm = 2*lh + 1
  result = newSeq[seq[float]](l.nt)
  for dt in 0..<l.nt: result[dt] = newSeq[float](nm*nm)
  var op = newRegOp(l, bt)
  var
    b = newGauge(l)
    x = newGauge(l)
  for mp in 0..<nm:
    loopSource(l, b, lsTri, lh, mp - lh, 0)
    let ci = regSolve(l, x, b, op, r2req, maxits)
    doAssert ci.converged, "jtopCorrExact: regSolve failed"
    for dt in 0..<l.nt:
      for m in 0..<nm:
        var s = 0.0
        for f in 0..<l.sph.nf:
          s += ylm(lh, m - lh, l.sph.faces[f].cc)*plaqSpatial(l, x, f, dt)
        result[dt][m*nm + mp] = s

proc loopCorrExact*(l: Lat, bt: Beta, shapes: openArray[LoopShape], lh, mh: int,
                    r2req = 1e-26, maxits = 200000): seq[seq[seq[float]]] =
  ## Exact free-theory correlator matrices of the Y_lm-projected loop basis:
  ##   result[dt][i][j] = <O_i(dt) O_j(0)>,
  ## one regSolve per shape.  The GEVP on this is the deterministic reference
  ## for the Monte-Carlo loop GEVP.
  result = newSeq[seq[seq[float]]](l.nt)
  for dt in 0..<l.nt:
    result[dt] = newSeq[seq[float]](shapes.len)
    for i in 0..<shapes.len: result[dt][i] = newSeq[float](shapes.len)
  var op = newRegOp(l, bt)
  var
    b = newGauge(l)
    x = newGauge(l)
  for j in 0..<shapes.len:
    loopSource(l, b, shapes[j], lh, mh, 0)
    let ci = regSolve(l, x, b, op, r2req, maxits)
    doAssert ci.converged, "loopCorrExact: regSolve failed"
    for dt in 0..<l.nt:
      for i in 0..<shapes.len:
        var s = 0.0
        let nc = loopCount(l, shapes[i])
        for c in 0..<nc:
          s += ylm(lh, mh, loopCenter(l.sph, shapes[i], c))*
               loopFlux(l, x, shapes[i], c, dt)
        result[dt][i][j] = s
