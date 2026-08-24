## Two-component Wilson-Dirac operator on S^2 x R  (WP-E).
##
## Normative reference: doc/02-formulation.md section 3 -- `applyDw` is Eq. (IV.1)
## written out; doc/04-interfaces.md section 7.
##
##   [D_W psi]_{y1,t} = sum_{y2 nn y1} kappa_{y1y2}
##       [ -1/2 (1 - e^a_{y1y2}(y1) sigma_a) U_{y1,t;y2,t} Omega_{y1y2} psi_{y2,t}
##         + 1/2 psi_{y1,t} ]
##     + kappa'_{y1} [ -1/2 (1 - sigma_3) U_{y1,t;y1,t+1} psi_{y1,t+1}
##                     -1/2 (1 + sigma_3) U_{y1,t;y1,t-1} psi_{y1,t-1} + psi_{y1,t} ]
##
## Conventions, all fixed here and used by every other work package:
##
## * Link angles are non-compact reals.  `Gauge.s[eIdx(e,t)]` is theta on the edge
##   in its canonical orientation a -> b, so the hop a -> b (row b, column a) carries
##   U = exp(+i theta) and the reverse hop carries its conjugate.  `Gauge.t[tIdx(y,t)]`
##   is theta on the temporal link (y,t) -> (y,t+1), so U_{y,t+1;y,t} = exp(+i theta).
##   Both are periodic in t.  Under `gaugeTransform` with `alpha`,
##   theta_e -> theta_e + alpha_b - alpha_a and theta^t_{y,t} -> theta^t + alpha_{y,t+1}
##   - alpha_{y,t}, matching the plaquettes of doc/02 section 5.
## * Fermions are antiperiodic in t: the temporal hop across the t = nt-1 -> 0 seam
##   carries an explicit -1.  Both directions of that link get it, which is what makes
##   the naive/Wilson split below exactly (anti)hermitian.
## * The spin connection is `Omega_ab = edges[e].sgn * expIsig3(edges[e].omega)` for the
##   canonical direction; the reverse hop uses its adjoint.  The Z2 `sgn` is the pole cut
##   (antiperiodicity in phi), so it is carried by the operator, not by the boundary.
##
## Performance: every apply is allocation-free.  `dst` must already have length
## `l.nsite` and must not alias `src`.

import std/[math, complex]
import ../core/lattice
import ../core/spinor

export lattice, spinor

type
  Gauge* = object
    ## Non-compact U(1) link angles: `s` per spatial link per time slice (index
    ## `eIdx(e,t)`), `t` per site per temporal link (index `tIdx(v,t)`).
    s*: seq[float]
    t*: seq[float]

  DwPart* = enum
    dwSpatial      ## the kappa sum over spatial neighbours of (IV.1)
    dwTemporal     ## the kappa' terms of (IV.1)
    dwNaive        ## the e^a sigma_a and sigma_3 pieces: the antihermitian C of (IV.4)
    dwWilson       ## the 1 and the diagonal pieces: the hermitian B of (IV.4)
  DwParts* = set[DwPart]

const
  dwAll*: DwParts = {dwSpatial, dwTemporal, dwNaive, dwWilson}
  dwC*: DwParts = {dwSpatial, dwTemporal, dwNaive}    ## antihermitian part of (IV.4)
  dwB*: DwParts = {dwSpatial, dwTemporal, dwWilson}   ## hermitian part of (IV.4)
  dwSpace*: DwParts = {dwSpatial, dwNaive, dwWilson}  ## the whole spatial operator
  dwTime*: DwParts = {dwTemporal, dwNaive, dwWilson}  ## the whole temporal operator

proc newGauge*(l: Lat): Gauge =
  ## All link angles zero, i.e. U = 1.
  Gauge(s: newSeq[float](l.sph.ne*l.nt), t: newSeq[float](l.sph.nv*l.nt))

proc nlinkS*(l: Lat): int = l.sph.ne*l.nt
proc nlinkT*(l: Lat): int = l.sph.nv*l.nt

proc zero*(u: var Gauge) =
  for i in 0..<u.s.len: u.s[i] = 0.0
  for i in 0..<u.t.len: u.t[i] = 0.0

# ---------------------------------------------------------------------------
# the kernel

type DwMode = enum dmFwd, dmAdj, dmDeriv

proc dwKernel(l: Lat, dst: var Spin, src: Spin, u, du: Gauge, m: float,
              parts: DwParts, mode: DwMode) =
  ## One pass over the lattice.  `dmAdj` assembles the conjugate transpose block by
  ## block from the *reverse* hop's own data (not from an identity relating the two
  ## tangents), so <x, D y> = <D^dag x, y> holds to roundoff independently of how
  ## accurately the geometry satisfies the tetrad hypothesis.  `dmDeriv` differentiates
  ## with respect to the link angles in the direction `du`, which kills every term
  ## that does not carry a U.
  let
    sph = l.sph
    nv = sph.nv
    ne = sph.ne
    nt = l.nt
    ew = if dwWilson in parts: 1.0 else: 0.0
    en = if dwNaive in parts: 1.0 else: 0.0
    doSpat = dwSpatial in parts
    doTemp = dwTemporal in parts
    adj = mode == dmAdj
    der = mode == dmDeriv
  doAssert dst.len == l.nsite, "dst must be sized for the lattice"
  doAssert src.len == l.nsite, "src must be sized for the lattice"
  for t in 0..<nt:
    let
      o = nv*t
      oe = ne*t
      tp = if t + 1 == nt: 0 else: t + 1
      tm = if t == 0: nt - 1 else: t - 1
      op = nv*tp
      om = nv*tm
      sf = if t + 1 == nt: -1.0 else: 1.0    ## antiperiodic seam, forward hop
      sb = if t == 0: -1.0 else: 1.0         ## antiperiodic seam, backward hop
    for y1 in 0..<nv:
      let i1 = o + y1
      var r: Spinor
      if doSpat:
        var dg = 0.0
        for k in 0..<sph.nbr[y1].len:
          let
            y2 = sph.nbr[y1][k]
            e = sph.nbe[y1][k]
            ed = sph.edges[e]
            kap = l.kap[e]
            d = if y1 == ed.a: 1.0 else: -1.0
            th = u.s[oe + e]
            hw = 0.5*ed.omega
            g = -0.5*kap*ed.sgn
            # U_{y1y2} Omega_{y1y2} = sgn diag(e^{i d (w/2 - th)}, e^{-i d (w/2 + th)})
            p0 = d*(hw - th)
            p1 = -d*(hw + th)
          dg += kap
          var
            f0 = complex64(g*cos(p0), g*sin(p0))
            f1 = complex64(g*cos(p1), g*sin(p1))
          if der:
            let q = -d*du.s[oe + e]           ## dU/dtheta = -i d dtheta U
            f0 = complex64(-q*f0.im, q*f0.re)
            f1 = complex64(-q*f1.im, q*f1.re)
          # e^a of the hop: at y1 for D, at y2 for D^dag
          var e0, e1: float
          if (d > 0.0) xor adj:
            e0 = ed.ea[0]
            e1 = ed.ea[1]
          else:
            e0 = -ed.eb[0]
            e1 = -ed.eb[1]
          let
            w = complex64(en*e0, en*e1)
            v = src[o + y2]
          if adj:
            # (f0, f1) . [ (ew - en e^a sigma_a) v ]
            let
              b0 = ew*v[0] - conjugate(w)*v[1]
              b1 = ew*v[1] - w*v[0]
            r[0] += f0*b0
            r[1] += f1*b1
          else:
            # (ew - en e^a sigma_a) [ (f0 v0, f1 v1) ]
            let
              a0 = f0*v[0]
              a1 = f1*v[1]
            r[0] += ew*a0 - conjugate(w)*a1
            r[1] += ew*a1 - w*a0
        if not der:
          let dh = 0.5*ew*dg
          r[0] += dh*src[i1][0]
          r[1] += dh*src[i1][1]
      if doTemp:
        let
          kt = l.kapT[y1]
          thp = u.t[o + y1]                  ## theta on (y,t) -> (y,t+1)
          thm = u.t[om + y1]                 ## theta on (y,t-1) -> (y,t)
        var
          zp = complex64(cos(thp), -sin(thp))   ## U_{y,t;y,t+1}
          zm = complex64(cos(thm), sin(thm))    ## U_{y,t;y,t-1}
        if der:
          let
            qp = -du.t[o + y1]
            qm = du.t[om + y1]
          zp = complex64(-qp*zp.im, qp*zp.re)
          zm = complex64(-qm*zm.im, qm*zm.re)
        # -1/2 (1 -+ sigma_3) split into Wilson (1) and naive (sigma_3);
        # the adjoint swaps the two projectors.
        var cf0, cf1, cb0, cb1: float
        if adj:
          cf0 = -0.5*(ew + en); cf1 = -0.5*(ew - en)
          cb0 = -0.5*(ew - en); cb1 = -0.5*(ew + en)
        else:
          cf0 = -0.5*(ew - en); cf1 = -0.5*(ew + en)
          cb0 = -0.5*(ew + en); cb1 = -0.5*(ew - en)
        let
          vp = src[op + y1]
          vm = src[om + y1]
          af = kt*sf
          ab = kt*sb
        r[0] += (af*cf0)*(zp*vp[0]) + (ab*cb0)*(zm*vm[0])
        r[1] += (af*cf1)*(zp*vp[1]) + (ab*cb1)*(zm*vm[1])
        if not der:
          let dh = kt*ew
          r[0] += dh*src[i1][0]
          r[1] += dh*src[i1][1]
      if m != 0.0 and not der:
        r[0] -= m*src[i1][0]
        r[1] -= m*src[i1][1]
      dst[i1] = r

proc applyDw*(l: Lat, dst: var Spin, src: Spin, u: Gauge, m = 0.0, parts = dwAll) =
  ## dst = (D_W - m) src, exactly Eq. (IV.1).  `parts` restricts the sum; the mass is
  ## always subtracted, so pass m = 0 with a restricted `parts`.
  dwKernel(l, dst, src, u, u, m, parts, dmFwd)

proc applyDwAdj*(l: Lat, dst: var Spin, src: Spin, u: Gauge, m = 0.0, parts = dwAll) =
  ## dst = (D_W - m)^dag src.
  dwKernel(l, dst, src, u, u, m, parts, dmAdj)

proc applyDwDeriv*(l: Lat, dst: var Spin, src: Spin, u: Gauge, du: Gauge,
                   parts = dwAll) =
  ## dst = (delta D_W)[du] src, the tangent of `applyDw` in the link angles.
  dwKernel(l, dst, src, u, du, 0.0, parts, dmDeriv)

# ---------------------------------------------------------------------------
# adjoint-mode derivative

proc dwPullback*(l: Lat, f: var Gauge, left, right: Spin, u: Gauge,
                 scale = 1.0, add = false, parts = dwAll) =
  ## f_link += scale * d[ 2 Re <left, D_W right> ] / d theta_link.
  ##
  ## Every term of (IV.1) that carries a link is of the form z(theta) = c exp(i s theta)
  ## with s = +-1, so d(2 Re z)/dtheta = -2 s Im z; the whole routine is that identity
  ## applied hop by hop.  Contracted with any `du`, this reproduces
  ## 2 Re <left, applyDwDeriv[du] right>.
  let
    sph = l.sph
    nv = sph.nv
    ne = sph.ne
    nt = l.nt
    ew = if dwWilson in parts: 1.0 else: 0.0
    en = if dwNaive in parts: 1.0 else: 0.0
    doSpat = dwSpatial in parts
    doTemp = dwTemporal in parts
  doAssert f.s.len == ne*nt
  doAssert f.t.len == nv*nt
  if not add:
    for i in 0..<f.s.len: f.s[i] = 0.0
    for i in 0..<f.t.len: f.t[i] = 0.0
  for t in 0..<nt:
    let
      o = nv*t
      oe = ne*t
      tp = if t + 1 == nt: 0 else: t + 1
      op = nv*tp
      sf = if t + 1 == nt: -1.0 else: 1.0
    if doSpat:
      for e in 0..<ne:
        let
          ed = sph.edges[e]
          kap = l.kap[e]
          th = u.s[oe + e]
          hw = 0.5*ed.omega
          g = -0.5*kap*ed.sgn
          xa = left[o + ed.a]
          xb = left[o + ed.b]
          va = right[o + ed.a]
          vb = right[o + ed.b]
        # hop b <- a: d = -1, U = e^{+i th}; e^a is the tangent at b toward a = -eb
        var
          p0 = -(hw - th)
          p1 = hw + th
          f0 = complex64(g*cos(p0), g*sin(p0))
          f1 = complex64(g*cos(p1), g*sin(p1))
          w = complex64(-en*ed.eb[0], -en*ed.eb[1])
          a0 = f0*va[0]
          a1 = f1*va[1]
          h0 = ew*a0 - conjugate(w)*a1
          h1 = ew*a1 - w*a0
          z = conjugate(xb[0])*h0 + conjugate(xb[1])*h1
        var g1 = -2.0*z.im                     ## s = +1
        # hop a <- b: d = +1, U = e^{-i th}; e^a is the tangent at a toward b = ea
        p0 = hw - th
        p1 = -(hw + th)
        f0 = complex64(g*cos(p0), g*sin(p0))
        f1 = complex64(g*cos(p1), g*sin(p1))
        w = complex64(en*ed.ea[0], en*ed.ea[1])
        a0 = f0*vb[0]
        a1 = f1*vb[1]
        h0 = ew*a0 - conjugate(w)*a1
        h1 = ew*a1 - w*a0
        z = conjugate(xa[0])*h0 + conjugate(xa[1])*h1
        g1 += 2.0*z.im                         ## s = -1
        f.s[oe + e] += scale*g1
    if doTemp:
      for y in 0..<nv:
        let
          kt = l.kapT[y]
          th = u.t[o + y]
          zp = complex64(cos(th), -sin(th))    ## U_{y,t;y,t+1}, s = -1
          zm = complex64(cos(th), sin(th))     ## U_{y,t+1;y,t}, s = +1
          xt = left[o + y]
          xp = left[op + y]
          vt = right[o + y]
          vp = right[op + y]
          cf0 = -0.5*(ew - en)
          cf1 = -0.5*(ew + en)
          cb0 = -0.5*(ew + en)
          cb1 = -0.5*(ew - en)
          af = kt*sf
          z1 = conjugate(xt[0])*((af*cf0)*(zp*vp[0])) +
               conjugate(xt[1])*((af*cf1)*(zp*vp[1]))
          z2 = conjugate(xp[0])*((af*cb0)*(zm*vt[0])) +
               conjugate(xp[1])*((af*cb1)*(zm*vt[1]))
        f.t[o + y] += scale*(2.0*z1.im - 2.0*z2.im)

# ---------------------------------------------------------------------------
# gauge transformations

proc gaugeTransform*(l: Lat, u: var Gauge, alpha: openArray[float]) =
  ## theta_e -> theta_e + alpha_b - alpha_a,  theta^t_{y,t} -> + alpha_{y,t+1} - alpha_{y,t}.
  ## `alpha` is indexed like a spinor field, `sIdx(v,t)`.
  let
    sph = l.sph
    nv = sph.nv
    ne = sph.ne
    nt = l.nt
  doAssert alpha.len == l.nsite
  for t in 0..<nt:
    let
      o = nv*t
      oe = ne*t
      op = nv*(if t + 1 == nt: 0 else: t + 1)
    for e in 0..<ne:
      u.s[oe + e] += alpha[o + sph.edges[e].b] - alpha[o + sph.edges[e].a]
    for y in 0..<nv:
      u.t[o + y] += alpha[op + y] - alpha[o + y]

proc spinGaugeTransform*(l: Lat, x: var Spin, alpha: openArray[float]) =
  ## psi_x -> exp(i alpha_x) psi_x.
  doAssert alpha.len == l.nsite
  doAssert x.len == l.nsite
  for i in 0..<l.nsite:
    let
      c = cos(alpha[i])
      s = sin(alpha[i])
    for a in 0..1:
      let
        re = x[i][a].re
        im = x[i][a].im
      x[i][a].re = c*re - s*im
      x[i][a].im = c*im + s*re

# ---------------------------------------------------------------------------
# volume-normalized kernel

proc hatScale*(l: Lat): seq[float] =
  ## abar / sqrt(A_y), one per sphere site: the diagonal similarity that turns D_W into
  ## the overlap kernel Dhat_W = abar^2 A^{-1/2} D_W A^{-1/2}.
  result = newSeq[float](l.sph.nv)
  for y in 0..<l.sph.nv:
    result[y] = l.sph.abar/sqrt(l.sph.area[y])

template hatBody(l: Lat, dst: var Spin, src: Spin, work: var Spin,
                 m: float, body: untyped) =
  let sph = l.sph
  doAssert work.len == l.nsite
  doAssert dst.len == l.nsite
  for y in 0..<sph.nv:
    let s = sph.abar/sqrt(sph.area[y])
    for t in 0..<l.nt:
      let i = y + sph.nv*t
      work[i][0] = s*src[i][0]
      work[i][1] = s*src[i][1]
  body
  for y in 0..<sph.nv:
    let s = sph.abar/sqrt(sph.area[y])
    for t in 0..<l.nt:
      let i = y + sph.nv*t
      dst[i][0] = s*dst[i][0] - m*src[i][0]
      dst[i][1] = s*dst[i][1] - m*src[i][1]

proc applyDwHat*(l: Lat, dst: var Spin, src: Spin, u: Gauge, work: var Spin,
                 m = 0.0, parts = dwAll) =
  ## dst = (Dhat_W - m) src with Dhat_W = abar^2 A^{-1/2} D_W A^{-1/2}, A = diag(A_y).
  ## `work` is a caller-owned scratch field of length `l.nsite`; nothing is allocated.
  ## The conjugation is by a real diagonal, so the inner product stays plain, the
  ## adjoint is the same routine with `applyDwAdj` inside, and a mass is site
  ## independent in physical units.  eig(Dhat_W) = abar * eig(D_continuum).
  hatBody(l, dst, src, work, m):
    applyDw(l, dst, work, u, 0.0, parts)

proc applyDwHatAdj*(l: Lat, dst: var Spin, src: Spin, u: Gauge, work: var Spin,
                    m = 0.0, parts = dwAll) =
  ## dst = (Dhat_W - m)^dag src.
  hatBody(l, dst, src, work, m):
    applyDwAdj(l, dst, work, u, 0.0, parts)

# ---------------------------------------------------------------------------
# dense assembly (tests and small-lattice spectra only)

proc denseDw*(l: Lat, u: Gauge, m = 0.0, parts = dwAll): seq[Complex64] =
  ## Column-major dense matrix of dimension 2*nsite, row index 2*sIdx + spin.
  ## Assembled block by block out of `Mat2` algebra -- deliberately a second,
  ## independent coding of (IV.1), so comparing it with `applyDw` has content.
  let
    sph = l.sph
    nv = sph.nv
    ne = sph.ne
    nt = l.nt
    nd = 2*l.nsite
    ew = if dwWilson in parts: 1.0 else: 0.0
    en = if dwNaive in parts: 1.0 else: 0.0
  result = newSeq[Complex64](nd*nd)
  template put(ri, ci: int, mm: Mat2) =
    let
      rr = 2*ri
      cc = 2*ci
    for a in 0..1:
      for b in 0..1:
        result[(rr + a) + nd*(cc + b)] += mm[a][b]
  for t in 0..<nt:
    let
      o = nv*t
      oe = ne*t
      tp = if t + 1 == nt: 0 else: t + 1
      tm = if t == 0: nt - 1 else: t - 1
      op = nv*tp
      om = nv*tm
      sf = if t + 1 == nt: -1.0 else: 1.0
      sb = if t == 0: -1.0 else: 1.0
    for y1 in 0..<nv:
      let i1 = o + y1
      if dwSpatial in parts:
        for k in 0..<sph.nbr[y1].len:
          let
            y2 = sph.nbr[y1][k]
            e = sph.nbe[y1][k]
            ed = sph.edges[e]
            kap = l.kap[e]
            d = if y1 == ed.a: 1.0 else: -1.0
            th = u.s[oe + e]
            eh = if d > 0.0: ed.ea else: [-ed.eb[0], -ed.eb[1]]
            omg = cr(ed.sgn)*expIsig3(d*ed.omega)        ## Omega_{y1y2}
            up = complex64(cos(d*th), -sin(d*th))        ## U_{y1y2}
            pr = (cr ew)*id2 - (cr en)*esig(eh)
          put(i1, o + y2, (cr(-0.5*kap)*up)*(pr*omg))
          put(i1, i1, cr(0.5*ew*kap)*id2)
      if dwTemporal in parts:
        let
          kt = l.kapT[y1]
          thp = u.t[o + y1]
          thm = u.t[om + y1]
          fp = (cr(-0.5*ew))*id2 + (cr(0.5*en))*sig3     ## -1/2 (1 - sigma_3)
          bp = (cr(-0.5*ew))*id2 - (cr(0.5*en))*sig3     ## -1/2 (1 + sigma_3)
        put(i1, op + y1, (cr(kt*sf)*complex64(cos(thp), -sin(thp)))*fp)
        put(i1, om + y1, (cr(kt*sb)*complex64(cos(thm), sin(thm)))*bp)
        put(i1, i1, cr(kt*ew)*id2)
      if m != 0.0:
        put(i1, i1, cr(-m)*id2)
