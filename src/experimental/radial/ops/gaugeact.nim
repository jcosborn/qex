## Non-compact U(1) gauge action on S^2 x R: the action (IV.24) with the free-limit
## couplings (IV.26), its analytic force, the gauge zero modes, an exact heatbath and
## the pseudo-inverse propagator (V.16)-(V.17).
##
## Normative reference: doc/02-formulation.md section 5 and section 7.2,
## doc/04-interfaces.md section 11.
##
##   S = sum_{tri,t} (beta_tri/2) Theta_tri(t)^2 + sum_{e,t} (beta_e/2) Theta_e(t)^2
##
## with Theta_tri the oriented sum of the three spatial link angles round a face and
## Theta_e the temporal plaquette of the spatial edge e = (a -> b),
##   Theta_e(t) = theta_e(t) + theta^t_b(t) - theta_e(t+1) - theta^t_a(t).
## Both gauge fields are PERIODIC in t (doc/04 section 4).
##
## The action is exactly Gaussian, S = theta^T M theta / 2 with M = C^T W C, C the
## plaquette-link incidence matrix and W = diag(beta).  Consequences used throughout:
## `gaugeForce` IS M theta; the heatbath is a Gaussian draw plus one solve; the
## gradient flow is linear.
##
## ker M is the gauge orbit d(alpha) PLUS one extra flat direction, the uniform
## temporal (Polyakov) mode theta^t = const, which is not a gauge mode because a
## gauge function periodic in t can only produce temporal shifts summing to zero
## round the time circle.  `projectGauge` removes the first, `projectKernel` both.

import std/math
import ../core/lattice
import wilson
import solve
import rng/threefry4x64

export lattice, threefry4x64
export CgInfo
export Gauge, newGauge, zero, gaugeTransform

## `Gauge` and `newGauge` come from `ops/wilson.nim` (doc/04 section 7) -- there is
## exactly one gauge-field type in the tree.  Note `wilson.gaugeTransform` is the
## same operation as `gradient` below, in place and with `+=`.

type
  GeomConv* = enum
    ## Which O(abar^2) transcription of the geometric weights of (IV.26) to use.
    ## The paper (section III) says the flat one is "equally possible" and differs
    ## at O(a^2).  gcExactArea is the paper's production convention -- both published
    ## Delta_0 values pin it to six digits (doc/06 "THE COUPLING CONVENTION") -- and
    ## is the default; `Edge.area` now stores this exact kite area.  gcGeodesic
    ## (flat area with geodesic lengths) reproduces the slide-8 spectrum legends.
    gcGeodesic     ## geodesic l, spherical excess A_tri, A_e = l (l*_1 + l*_2)/2
    gcExactArea    ## as gcGeodesic but A_e = the exact spherical diamond area
    gcFlat         ## chord l, planar A_tri, in-plane l*

  Beta* = object
    ## Precomputed plaquette couplings (IV.26).  Hoisting the divisions out of the
    ## force kernel matters: `gaugeForce` is the CG inner loop.
    face*: seq[float]        ## beta_tri, one per face
    edge*: seq[float]        ## beta_l, one per spatial edge
    afac*: seq[float]        ## the A_tri that went into beta_tri (J^t normalization)
    g2*: float
    conv*: GeomConv

func nlink*(l: Lat): int = (l.sph.ne + l.sph.nv)*l.nt
  ## Total number of link degrees of freedom.
func slink*(l: Lat): int = l.sph.ne*l.nt
  ## Flat index range of the spatial links; temporal links follow.

# --- vector space -----------------------------------------------------------

proc `:=`*(x: var Gauge, y: Gauge) =
  for i in 0..<x.s.len: x.s[i] = y.s[i]
  for i in 0..<x.t.len: x.t[i] = y.t[i]

proc axpy*(x: var Gauge, a: float, y: Gauge) =
  for i in 0..<x.s.len: x.s[i] += a*y.s[i]
  for i in 0..<x.t.len: x.t[i] += a*y.t[i]

proc axpby*(x: var Gauge, a: float, y: Gauge, b: float) =
  ## x = a*y + b*x, the CG search-direction update.
  for i in 0..<x.s.len: x.s[i] = a*y.s[i] + b*x.s[i]
  for i in 0..<x.t.len: x.t[i] = a*y.t[i] + b*x.t[i]

proc scale*(x: var Gauge, a: float) =
  for i in 0..<x.s.len: x.s[i] *= a
  for i in 0..<x.t.len: x.t[i] *= a

func dot*(x, y: Gauge): float =
  for i in 0..<x.s.len: result += x.s[i]*y.s[i]
  for i in 0..<x.t.len: result += x.t[i]*y.t[i]

func norm2*(x: Gauge): float = dot(x, x)

func toSeq*(u: Gauge): seq[float] =
  ## Flat link vector: spatial links first (index eIdx), then temporal (slink + tIdx).
  result = newSeq[float](u.s.len + u.t.len)
  for i in 0..<u.s.len: result[i] = u.s[i]
  for i in 0..<u.t.len: result[u.s.len + i] = u.t[i]

proc fromSeq*(l: Lat, v: openArray[float]): Gauge =
  result = newGauge(l)
  for i in 0..<result.s.len: result.s[i] = v[i]
  for i in 0..<result.t.len: result.t[i] = v[result.s.len + i]

proc unitSource*(l: Lat, b: var Gauge, link: int) =
  ## b = e_link in the flat link indexing of `toSeq`.
  b.zero
  if link < b.s.len: b.s[link] = 1.0
  else: b.t[link - b.s.len] = 1.0

# --- couplings --------------------------------------------------------------

proc flatWeights(sph: Sphere): tuple[fa, el, ed: seq[float]] =
  ## Chordal/planar counterparts of Face.area, Edge.len and Edge.dual.
  ## arXiv:2510.03085 section III says the flat (simplicial) convention is
  ## "equally possible" and differs at O(a^2); this makes that variant testable.
  ## Planar circumcenter O = A + (|AC|^2 (n x AB) - |AB|^2 (n x AC))/(2|n|^2),
  ## n = AB x AC; l* is the in-plane signed distance from O to the edge, positive
  ## on the side of the third vertex.
  result.fa = newSeq[float](sph.nf)
  result.el = newSeq[float](sph.ne)
  result.ed = newSeq[float](sph.ne)
  for f in 0..<sph.nf:
    let
      fc = sph.faces[f]
      a = sph.pos[fc.v[0]]
      b = sph.pos[fc.v[1]]
      c = sph.pos[fc.v[2]]
      ab = b - a
      ac = c - a
      n = cross(ab, ac)
      n2 = dot(n, n)
      o = a + (1.0/(2.0*n2))*(dot(ac, ac)*cross(n, ab) - dot(ab, ab)*cross(n, ac))
      nn = unit n
    result.fa[f] = 0.5*sqrt(n2)
    for i in 0..2:
      let
        p = sph.pos[fc.v[i]]
        q = sph.pos[fc.v[(i+1) mod 3]]
        z = sph.pos[fc.v[(i+2) mod 3]]
        d = q - p
        le = norm d
      var u = (1.0/le)*cross(nn, d)
      if dot(u, z - p) < 0.0: u = -1.0*u
      result.el[fc.e[i]] = le
      result.ed[fc.e[i]] = result.ed[fc.e[i]] + dot(o - p, u)

proc newBeta*(l: Lat, g2: float, conv = gcExactArea): Beta =
  ## (IV.26): beta_tri = a_t/(g2 A_tri), beta_l = 2 A_l/(g2 l^2 a_t).
  result.g2 = g2
  result.conv = conv
  result.face = newSeq[float](l.sph.nf)
  result.edge = newSeq[float](l.sph.ne)
  result.afac = newSeq[float](l.sph.nf)
  if conv == gcFlat:
    let w = flatWeights(l.sph)
    for f in 0..<l.sph.nf:
      result.afac[f] = w.fa[f]
      result.face[f] = l.at/(g2*w.fa[f])
    for e in 0..<l.sph.ne:
      let ae = 0.5*w.el[e]*w.ed[e]
      result.edge[e] = 2.0*ae/(g2*w.el[e]*w.el[e]*l.at)
  else:
    for f in 0..<l.sph.nf:
      result.afac[f] = l.sph.faces[f].area
      result.face[f] = l.betaFace(f, g2)
    for e in 0..<l.sph.ne:
      if conv == gcExactArea:
        result.edge[e] = l.betaEdge(e, g2)          # Edge.area is the exact kite area
      else:
        let ed = l.sph.edges[e]                     # flat A_e = l dual/2, slide-8 legend
        result.edge[e] = ed.dual/(g2*ed.len*l.at)

# --- plaquettes -------------------------------------------------------------

func plaqSpatial*(l: Lat, u: Gauge, f, t: int): float =
  ## Theta_tri(f, t), the oriented sum round the spatial triangle.
  let fc = l.sph.faces[f]
  for i in 0..2:
    result += float(fc.s[i])*u.s[eIdx(l, fc.e[i], t)]

func plaqTemporal*(l: Lat, u: Gauge, e, t: int): float =
  ## Theta_e(t) = theta_e(t) + theta^t_b(t) - theta_e(t+1) - theta^t_a(t).
  let ed = l.sph.edges[e]
  u.s[eIdx(l, e, t)] + u.t[tIdx(l, ed.b, t)] -
    u.s[eIdx(l, e, t+1)] - u.t[tIdx(l, ed.a, t)]

func jtop*(l: Lat, u: Gauge, f, t: int): float =
  ## (V.12): J^t_lat = Theta_tri / A_tri.
  plaqSpatial(l, u, f, t)/l.sph.faces[f].area

func jtop*(l: Lat, u: Gauge, b: Beta, f, t: int): float =
  ## Same, with the area convention that produced `b`.
  plaqSpatial(l, u, f, t)/b.afac[f]

# --- action and force -------------------------------------------------------

proc gaugeActionParts*(l: Lat, u: Gauge, b: Beta): tuple[sp, tp: float] =
  ## (IV.24) split into the spatial (magnetic) and temporal (electric) sums.
  let
    sph = l.sph
    ne = sph.ne
    nv = sph.nv
    nt = l.nt
  var ss = 0.0
  for f in 0..<sph.nf:
    let
      fc = sph.faces[f]
      e0 = fc.e[0]
      e1 = fc.e[1]
      e2 = fc.e[2]
      s0 = float(fc.s[0])
      s1 = float(fc.s[1])
      s2 = float(fc.s[2])
    var a = 0.0
    for t in 0..<nt:
      let
        o = ne*t
        th = s0*u.s[e0+o] + s1*u.s[e1+o] + s2*u.s[e2+o]
      a += th*th
    ss += 0.5*b.face[f]*a
  var st = 0.0
  for e in 0..<ne:
    let
      ia = sph.edges[e].a
      ib = sph.edges[e].b
    var a = 0.0
    for t in 0..<nt:
      let
        o = ne*t
        o1 = if t+1 == nt: 0 else: ne*(t+1)
        ov = nv*t
        th = u.s[e+o] + u.t[ib+ov] - u.s[e+o1] - u.t[ia+ov]
      a += th*th
    st += 0.5*b.edge[e]*a
  (ss, st)

proc gaugeAction*(l: Lat, u: Gauge, b: Beta): float =
  let p = gaugeActionParts(l, u, b)
  p.sp + p.tp

proc gaugeAction*(l: Lat, u: Gauge, g2: float): float =
  gaugeAction(l, u, newBeta(l, g2))

proc gaugeForce*(l: Lat, f: var Gauge, u: Gauge, b: Beta) =
  ## f = dS/dtheta.  The action is quadratic, so this is exactly M u.
  let
    sph = l.sph
    ne = sph.ne
    nv = sph.nv
    nt = l.nt
  f.zero
  for fa in 0..<sph.nf:
    let
      fc = sph.faces[fa]
      bf = b.face[fa]
      e0 = fc.e[0]
      e1 = fc.e[1]
      e2 = fc.e[2]
      s0 = float(fc.s[0])
      s1 = float(fc.s[1])
      s2 = float(fc.s[2])
    for t in 0..<nt:
      let
        o = ne*t
        w = bf*(s0*u.s[e0+o] + s1*u.s[e1+o] + s2*u.s[e2+o])
      f.s[e0+o] += s0*w
      f.s[e1+o] += s1*w
      f.s[e2+o] += s2*w
  for e in 0..<ne:
    let
      be = b.edge[e]
      ia = sph.edges[e].a
      ib = sph.edges[e].b
    for t in 0..<nt:
      let
        o = ne*t
        o1 = if t+1 == nt: 0 else: ne*(t+1)
        ov = nv*t
        w = be*(u.s[e+o] + u.t[ib+ov] - u.s[e+o1] - u.t[ia+ov])
      f.s[e+o] += w
      f.t[ib+ov] += w
      f.s[e+o1] -= w
      f.t[ia+ov] -= w

proc gaugeForce*(l: Lat, f: var Gauge, u: Gauge, g2: float) =
  gaugeForce(l, f, u, newBeta(l, g2))

proc mDiagonal*(l: Lat, d: var Gauge, b: Beta) =
  ## The diagonal of M, link by link.  Used to size the gradient-flow step:
  ## lambda_max(M) is a small multiple of max_i M_ii (Gershgorin with the
  ## plaquette row sums), and an explicit RK step needs h lambda_max = O(1).
  let
    sph = l.sph
    ne = sph.ne
    nv = sph.nv
    nt = l.nt
  d.zero
  for f in 0..<sph.nf:
    for t in 0..<nt:
      let o = ne*t
      for i in 0..2: d.s[sph.faces[f].e[i]+o] += b.face[f]
  for e in 0..<ne:
    let
      be = b.edge[e]
      ia = sph.edges[e].a
      ib = sph.edges[e].b
    for t in 0..<nt:
      let
        o = ne*t
        o1 = if t+1 == nt: 0 else: ne*(t+1)
        ov = nv*t
      d.s[e+o] += be
      d.s[e+o1] += be
      d.t[ib+ov] += be
      d.t[ia+ov] += be

proc mDiagMax*(l: Lat, b: Beta): float =
  var d = newGauge(l)
  mDiagonal(l, d, b)
  for x in d.s: result = max(result, x)
  for x in d.t: result = max(result, x)

template applyM*(l: Lat, dst: var Gauge, src: Gauge, b: Beta) =
  ## The Gaussian kernel of the action; identical to `gaugeForce` by construction.
  gaugeForce(l, dst, src, b)

# --- zero modes -------------------------------------------------------------

proc gradient*(l: Lat, p: var Gauge, alpha: openArray[float]) =
  ## p = d alpha, the pure-gauge direction:
  ##   p_e(t)   = alpha_b(t) - alpha_a(t)
  ##   p^t_v(t) = alpha_v(t+1) - alpha_v(t)
  ## `alpha` is indexed by sIdx(v, t).
  let
    sph = l.sph
    ne = sph.ne
    nv = sph.nv
    nt = l.nt
  if p.s.len != ne*nt: p = newGauge(l)
  for t in 0..<nt:
    let
      ov = nv*t
      ov1 = if t+1 == nt: 0 else: nv*(t+1)
      o = ne*t
    for e in 0..<ne:
      p.s[e+o] = alpha[sph.edges[e].b + ov] - alpha[sph.edges[e].a + ov]
    for v in 0..<nv:
      p.t[v+ov] = alpha[v+ov1] - alpha[v+ov]

proc divergence*(l: Lat, d: var seq[float], p: Gauge) =
  ## d = d^dagger p, the adjoint of `gradient` for the plain Euclidean products:
  ##   d_v(t) = sum_{e: b=v} p_e - sum_{e: a=v} p_e + p^t_v(t-1) - p^t_v(t).
  let
    sph = l.sph
    ne = sph.ne
    nv = sph.nv
    nt = l.nt
  if d.len != nv*nt: d = newSeq[float](nv*nt)
  for i in 0..<d.len: d[i] = 0.0
  for t in 0..<nt:
    let
      ov = nv*t
      ovm = if t == 0: nv*(nt-1) else: nv*(t-1)
      o = ne*t
    for e in 0..<ne:
      let x = p.s[e+o]
      d[sph.edges[e].b + ov] += x
      d[sph.edges[e].a + ov] -= x
    for v in 0..<nv:
      d[v+ov] += p.t[v+ovm] - p.t[v+ov]

proc laplace*(l: Lat, dst: var seq[float], src: openArray[float]) =
  ## d^dagger d alpha, the graph Laplacian of S^2 x S^1_t.  Positive semidefinite
  ## with kernel exactly the constants (the graph is connected).
  let
    sph = l.sph
    nv = sph.nv
    nt = l.nt
  if dst.len != nv*nt: dst = newSeq[float](nv*nt)
  for t in 0..<nt:
    let
      ov = nv*t
      ovp = if t+1 == nt: 0 else: nv*(t+1)
      ovm = if t == 0: nv*(nt-1) else: nv*(t-1)
    for v in 0..<nv:
      let x = src[v+ov]
      var s = 2.0*x - src[v+ovp] - src[v+ovm]
      for w in sph.nbr[v]: s += x - src[w+ov]
      dst[v+ov] = s

proc removeMean(x: var seq[float]) =
  var s = 0.0
  for v in x: s += v
  s /= float(x.len)
  for i in 0..<x.len: x[i] -= s

proc projectGauge*(l: Lat, p: var Gauge, r2req = 1e-24, maxits = 10000): CgInfo =
  ## Remove the gauge-orbit component of `p`: solve (d^dagger d) alpha = d^dagger p
  ## by CG with the constant mode projected out at every step, then p -= d alpha.
  ## Afterwards d^dagger p = the CG residual, so the tolerance is a direct statement
  ## about the transversality of the result.
  let n = l.sph.nv*l.nt
  var b = newSeq[float](n)
  divergence(l, b, p)
  removeMean b
  var b2 = 0.0
  for v in b: b2 += v*v
  if b2 == 0.0: return CgInfo(iters: 0, r2: 0.0, converged: true)
  var
    al = newSeq[float](n)
    r = b
    q = b
    ap = newSeq[float](n)
    r2 = b2
    its = 0
  let stop = r2req*b2
  while r2 > stop and its < maxits:
    laplace(l, ap, q)
    removeMean ap
    var pap = 0.0
    for i in 0..<n: pap += q[i]*ap[i]
    if pap <= 0.0: break
    let a = r2/pap
    var r2n = 0.0
    for i in 0..<n:
      al[i] += a*q[i]
      r[i] -= a*ap[i]
      r2n += r[i]*r[i]
    let bt = r2n/r2
    for i in 0..<n: q[i] = r[i] + bt*q[i]
    r2 = r2n
    inc its
  laplace(l, ap, al)
  removeMean ap
  var t2 = 0.0
  for i in 0..<n:
    let d = b[i] - ap[i]
    t2 += d*d
  var da = newGauge(l)
  gradient(l, da, al)
  axpy(p, -1.0, da)
  CgInfo(iters: its, r2: t2/b2, converged: t2 <= 1.001*r2req*b2)

proc projectFlat*(l: Lat, p: var Gauge) =
  ## Remove the uniform temporal (Polyakov) direction theta^t = const, theta^s = 0.
  ## It costs no action -- sum_t Theta_e(t) = W_b - W_a only sees differences of the
  ## temporal Wilson lines -- and is orthogonal to every gauge direction, because
  ## sum_t (alpha_v(t+1) - alpha_v(t)) = 0.
  var s = 0.0
  for x in p.t: s += x
  s /= float(p.t.len)
  for i in 0..<p.t.len: p.t[i] -= s

proc projectKernel*(l: Lat, p: var Gauge, r2req = 1e-24, maxits = 10000): CgInfo =
  ## Project out all of ker M = range(d) (+) span(uniform temporal).
  result = projectGauge(l, p, r2req, maxits)
  projectFlat(l, p)

# --- solvers ----------------------------------------------------------------

proc cgM*(l: Lat, x: var Gauge, b: Gauge, bt: Beta,
          r2req = 1e-20, maxits = 20000): CgInfo =
  ## CG for M x = b from x = 0.  `b` must lie in range(M): the Krylov space then
  ## never touches the kernel, which is exactly the (V.16)-(V.17) prescription.
  if x.s.len != b.s.len: x = newGauge(l)
  let b2 = norm2(b)
  if b2 == 0.0:
    x.zero
    return CgInfo(iters: 0, r2: 0.0, converged: true)
  var
    r = newGauge(l)
    q = newGauge(l)
    ap = newGauge(l)
  x.zero
  r := b
  q := b
  var
    r2 = b2
    its = 0
  let stop = r2req*b2
  while r2 > stop and its < maxits:
    gaugeForce(l, ap, q, bt)
    let pap = dot(q, ap)
    if pap <= 0.0: break
    let a = r2/pap
    axpy(x, a, q)
    axpy(r, -a, ap)
    let r2n = norm2(r)
    axpby(q, 1.0, r, r2n/r2)
    r2 = r2n
    inc its
  gaugeForce(l, ap, x, bt)
  var t2 = 0.0
  for i in 0..<b.s.len:
    let d = b.s[i] - ap.s[i]
    t2 += d*d
  for i in 0..<b.t.len:
    let d = b.t[i] - ap.t[i]
    t2 += d*d
  CgInfo(iters: its, r2: t2/b2, converged: t2 <= 1.001*r2req*b2)

proc pseudoSolve*(l: Lat, x: var Gauge, b: Gauge, bt: Beta,
                  r2req = 1e-20, maxits = 20000): tuple[proj, sol: CgInfo] =
  ## x = Mtilde^{-1} b, the double-CG pseudo-inverse of (V.16)-(V.17):
  ##   b' = Mtilde^{-1}(M b)   -- project the source onto range(M)
  ##   x  = Mtilde^{-1} b'     -- the actual solve
  ## Both CGs start from zero, so neither Krylov space contains the kernel.
  var mb = newGauge(l)
  gaugeForce(l, mb, b, bt)
  var bp = newGauge(l)
  result.proj = cgM(l, bp, mb, bt, r2req, maxits)
  result.sol = cgM(l, x, bp, bt, r2req, maxits)

type RegOp* = object
  ## M regularized on its kernel:  A = M + sig * d d^dagger + tau * P P^T/|P|^2.
  ## Both added blocks vanish on range(M) = ker(M)^perp and are positive definite on
  ## ker(M), so A is symmetric POSITIVE DEFINITE on the whole link space while
  ## A^{-1} b = Mtilde^{-1} b exactly for every b orthogonal to the kernel -- which
  ## every gauge-invariant source is.  This is the "add G G^dagger" regularization,
  ## and it is what makes tight tolerances reachable: the literal (V.16)-(V.17)
  ## double CG has no control over the kernel component that roundoff injects into
  ## the Krylov space, and past r2 ~ 1e-26 that component takes over and diverges
  ## (measured: at r2req = 1e-28 the L=1, nt=120 solve blows up).
  bt*: Beta
  sig*, tau*: float
  d: seq[float]
  g: Gauge

proc newRegOp*(l: Lat, bt: Beta, sig = 0.0, tau = 0.0): RegOp =
  ## sig/tau default to the mean plaquette coupling, which puts the regularized
  ## block on the same scale as M and leaves the condition number alone.
  var s = 0.0
  for x in bt.face: s += x
  for x in bt.edge: s += x
  s /= float(bt.face.len + bt.edge.len)
  RegOp(bt: bt, sig: (if sig > 0.0: sig else: s), tau: (if tau > 0.0: tau else: s),
        d: newSeq[float](l.sph.nv*l.nt), g: newGauge(l))

proc applyReg*(l: Lat, o: var RegOp, dst: var Gauge, src: Gauge) =
  gaugeForce(l, dst, src, o.bt)
  divergence(l, o.d, src)
  gradient(l, o.g, o.d)
  axpy(dst, o.sig, o.g)
  var p = 0.0
  for x in src.t: p += x
  p = o.tau*p/float(src.t.len)
  for i in 0..<dst.t.len: dst.t[i] += p

proc regSolve*(l: Lat, x: var Gauge, b: Gauge, o: var RegOp,
               r2req = 1e-24, maxits = 100000): CgInfo =
  ## CG on the regularized operator.  For b in range(M) this returns Mtilde^{-1} b
  ## and is stable to the roundoff floor.
  if x.s.len != b.s.len: x = newGauge(l)
  let b2 = norm2(b)
  if b2 == 0.0:
    x.zero
    return CgInfo(iters: 0, r2: 0.0, converged: true)
  var
    r = newGauge(l)
    q = newGauge(l)
    ap = newGauge(l)
  x.zero
  r := b
  q := b
  var
    r2 = b2
    its = 0
  let stop = r2req*b2
  while r2 > stop and its < maxits:
    applyReg(l, o, ap, q)
    let pap = dot(q, ap)
    if pap <= 0.0: break
    let a = r2/pap
    axpy(x, a, q)
    axpy(r, -a, ap)
    let r2n = norm2(r)
    axpby(q, 1.0, r, r2n/r2)
    r2 = r2n
    inc its
  applyReg(l, o, ap, x)
  var t2 = 0.0
  for i in 0..<b.s.len:
    let d = b.s[i] - ap.s[i]
    t2 += d*d
  for i in 0..<b.t.len:
    let d = b.t[i] - ap.t[i]
    t2 += d*d
  CgInfo(iters: its, r2: t2/b2, converged: t2 <= 1.001*r2req*b2)

proc gaugePropagator*(l: Lat, g2: float, srcLink: int,
                      r2req = 1e-20, maxits = 20000): seq[float] =
  ## Column `srcLink` of <theta_m theta_n> = (Mtilde^{-1})_{mn}, flat link indexing.
  let bt = newBeta(l, g2)
  var b = newGauge(l)
  unitSource(l, b, srcLink)
  var x = newGauge(l)
  discard pseudoSolve(l, x, b, bt, r2req, maxits)
  x.toSeq

# --- exact heatbath ---------------------------------------------------------

proc heatbathSource*(l: Lat, b: var Gauge, bt: Beta, r: var Threefry4x64) =
  ## b = C^T W^{1/2} xi with xi ~ N(0,1) per plaquette, so <b b^T> = C^T W C = M
  ## exactly, and b lies in range(M) = range(C^T) by construction.
  let
    sph = l.sph
    ne = sph.ne
    nv = sph.nv
    nt = l.nt
  b.zero
  for fa in 0..<sph.nf:
    let
      fc = sph.faces[fa]
      w = sqrt(bt.face[fa])
    for t in 0..<nt:
      let
        o = ne*t
        x = w*r.gaussian
      for i in 0..2: b.s[fc.e[i]+o] += float(fc.s[i])*x
  for e in 0..<ne:
    let
      w = sqrt(bt.edge[e])
      ia = sph.edges[e].a
      ib = sph.edges[e].b
    for t in 0..<nt:
      let
        o = ne*t
        o1 = if t+1 == nt: 0 else: ne*(t+1)
        ov = nv*t
        x = w*r.gaussian
      b.s[e+o] += x
      b.t[ib+ov] += x
      b.s[e+o1] -= x
      b.t[ia+ov] -= x

proc heatbath*(l: Lat, u: var Gauge, bt: Beta, r: var Threefry4x64,
               r2req = 1e-20, maxits = 20000): CgInfo =
  ## Exact Gaussian sampling of the transverse modes.  theta = M^+ b with
  ## <b b^T> = M gives <theta theta^T> = M^+ M M^+ = M^+, the exact free-field
  ## covariance, with the kernel of M set to zero.  No HMC, no accept/reject.
  var b = newGauge(l)
  heatbathSource(l, b, bt, r)
  cgM(l, u, b, bt, r2req, maxits)

proc heatbath*(l: Lat, u: var Gauge, g2: float, r: var Threefry4x64) =
  discard heatbath(l, u, newBeta(l, g2), r)

# --- observable sources -----------------------------------------------------

proc triSource*(l: Lat, b: var Gauge, f, t: int) =
  ## b = the incidence row of the spatial plaquette (f, t), so that
  ## <Theta_f(t) Theta_f'(t')> = b^T Mtilde^{-1} b'.
  b.zero
  let fc = l.sph.faces[f]
  for i in 0..2: b.s[eIdx(l, fc.e[i], t)] += float(fc.s[i])
