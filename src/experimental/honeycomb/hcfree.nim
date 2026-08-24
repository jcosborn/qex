## Free (U = 1) Wilson-type Dirac operator in momentum space, for the 16-cell
## honeycomb (D4*) and for the ordinary hypercubic lattice.
##
## This module is deliberately dependency-light: `std/[math, complex]` plus
## `hcgeom` for the direction tables.  No QEX `Layout`, no fields, no threads --
## everything here is 4x4 / 8x8 linear algebra and scalar momentum-space
## arithmetic.  See doc/FORMULATION.md section 5 for the derivations.
##
## Conventions (FORMULATION section 7): `a = 1`, mu = 0..3 with mu = 3 the time
## direction, gamma matrices in the DeGrand-Rossi basis used by QEX
## (`src/physics/spinOld.nim`, `gamma1..gamma4` here called `gammaMat[0..3]`).
##
## The operator
## ------------
## For a set of unit neighbour vectors `n_i` with prefactor `pref`
## (`1/6` for the 24 honeycomb directions, `1/2` for the 8 cubic ones),
##
##   D(p)   = m + M(p) + i gamma.K(p)
##   M(p)   = r*pref * sum_i ( 1 - cos(p.n_i) )
##   K_mu(p)= pref   * sum_i n_{mu,i} sin(p.n_i)
##
## Since `(gamma.K)^2 = |K|^2`, the eigenvalues are `m + M +- i|K|`, each doubly
## degenerate -- exact, no eigensolver needed.
##
## Brillouin zone
## --------------
## The honeycomb site set `Z^4 u (Z+1/2)^4` has reciprocal lattice `2*pi*D4`
## (covolume `2*(2*pi)^4`), so its BZ is *twice* the cubic one and an
## `Ns^3 x Nt` cell lattice carries `2*Ns^3*Nt` momenta.  A complete set is the
## cubic set with the range of one component doubled; `momenta` doubles `p_0`,
## `momentaAlt` doubles `p_3`.  (NB: FORMULATION section 5.3 states the BZ
## volume as `(2*pi)^4/2`; that is a typo, it is `2*(2*pi)^4`.)

import std/[math, complex]
import hcgeom

type
  HcLat* = enum        ## which lattice
    lCubic = "cubic",
    lHoneycomb = "16cell"

  Vec4* = array[nDim, float]
  Mat4* = array[4, array[4, Complex64]]
  Mat8* = array[8, array[8, Complex64]]

const
  cubicDirsF* = block:
    var a: array[nAxis, Vec4]
    for i in 0..<nAxis: a[i] = toFloat(allDirs[i])
    a
  hcDirsF* = block:
    var a: array[nDirs, Vec4]
    for i in 0..<nDirs: a[i] = toFloat(allDirs[i])
    a

func nNeighbours*(lat: HcLat): int {.inline.} =
  if lat == lCubic: nAxis else: nDirs

func pref*(lat: HcLat): float {.inline.} =
  ## the `1/6` (honeycomb) / `1/2` (cubic) of FORMULATION 5.1
  if lat == lCubic: 0.5 else: 1.0/6.0

func dirsOf*(lat: HcLat): seq[Vec4] =
  if lat == lCubic:
    for d in cubicDirsF: result.add d
  else:
    for d in hcDirsF: result.add d

func dot*(a, b: Vec4): float {.inline.} =
  a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3]

func norm2*(a: Vec4): float {.inline.} = dot(a, a)

# ---------------------------------------------------------------------------
# M(p) and K(p): direct sum over the neighbour vectors
# ---------------------------------------------------------------------------

func freeM*(p: Vec4; lat: HcLat; r = 1.0): float =
  ## Wilson term `M(p) = r*pref*sum_i (1-cos(p.n_i))`, by direct summation.
  var s = 0.0
  if lat == lCubic:
    for n in cubicDirsF: s += 1.0 - cos(dot(p, n))
  else:
    for n in hcDirsF: s += 1.0 - cos(dot(p, n))
  r*pref(lat)*s

func freeK*(p: Vec4; lat: HcLat): Vec4 =
  ## `K_mu(p) = pref*sum_i n_{mu,i} sin(p.n_i)`, by direct summation.
  template acc(tab: untyped) =
    for n in tab:
      let s = sin(dot(p, n))
      for mu in 0..<nDim: result[mu] += n[mu]*s
  if lat == lCubic: acc(cubicDirsF) else: acc(hcDirsF)
  let f = pref(lat)
  for mu in 0..<nDim: result[mu] *= f

# ---------------------------------------------------------------------------
# closed forms (FORMULATION 5.3)
# ---------------------------------------------------------------------------

func freeMClosed*(p: Vec4; lat: HcLat; r = 1.0): float =
  ## cubic:     `M = r sum_mu (1 - cos p_mu)`
  ## honeycomb: `M = 4r - (r/6)[ 2 sum_mu cos p_mu + 16 prod_mu cos(p_mu/2) ]`
  if lat == lCubic:
    var s = 0.0
    for mu in 0..<nDim: s += 1.0 - cos(p[mu])
    r*s
  else:
    var sc = 0.0
    var pc = 1.0
    for mu in 0..<nDim:
      sc += cos(p[mu])
      pc *= cos(0.5*p[mu])
    4.0*r - (r/6.0)*(2.0*sc + 16.0*pc)

func freeKClosed*(p: Vec4; lat: HcLat): Vec4 =
  ## cubic:     `K_mu = sin p_mu`
  ## honeycomb: `K_mu = (1/6)[ 2 sin p_mu + 8 sin(p_mu/2) prod_{nu/=mu} cos(p_nu/2) ]`
  if lat == lCubic:
    for mu in 0..<nDim: result[mu] = sin(p[mu])
  else:
    for mu in 0..<nDim:
      var pc = 1.0
      for nu in 0..<nDim:
        if nu != mu: pc *= cos(0.5*p[nu])
      result[mu] = (2.0*sin(p[mu]) + 8.0*sin(0.5*p[mu])*pc)/6.0

func absK*(p: Vec4; lat: HcLat): float =
  let k = freeKClosed(p, lat)
  sqrt(norm2(k))

func freeEigs*(p: Vec4; lat: HcLat; r = 1.0; m = 0.0): array[2, Complex64] =
  ## The two distinct eigenvalues `m + M(p) +- i|K(p)|` of `D(p)`; each is
  ## doubly degenerate, so the 4x4 spectrum is `{e[0],e[0],e[1],e[1]}`.
  let
    mm = m + freeMClosed(p, lat, r)
    kk = absK(p, lat)
  [complex64(mm, kk), complex64(mm, -kk)]

# ---------------------------------------------------------------------------
# gamma matrices (DeGrand-Rossi, as in QEX src/physics/spinOld.nim)
# ---------------------------------------------------------------------------

const
  z0 = complex64(0.0, 0.0)
  z1 = complex64(1.0, 0.0)
  zi = complex64(0.0, 1.0)
  n1 = complex64(-1.0, 0.0)
  ni = complex64(0.0, -1.0)

const
  gammaMat*: array[4, Mat4] = [
    # gamma1
    [[ z0, z0, z0, zi ],
     [ z0, z0, zi, z0 ],
     [ z0, ni, z0, z0 ],
     [ ni, z0, z0, z0 ]],
    # gamma2
    [[ z0, z0, z0, n1 ],
     [ z0, z0, z1, z0 ],
     [ z0, z1, z0, z0 ],
     [ n1, z0, z0, z0 ]],
    # gamma3
    [[ z0, z0, zi, z0 ],
     [ z0, z0, z0, ni ],
     [ ni, z0, z0, z0 ],
     [ z0, zi, z0, z0 ]],
    # gamma4 (time)
    [[ z0, z0, z1, z0 ],
     [ z0, z0, z0, z1 ],
     [ z1, z0, z0, z0 ],
     [ z0, z1, z0, z0 ]]]

  gamma5Mat*: Mat4 =
    [[ z1, z0, z0, z0 ],
     [ z0, z1, z0, z0 ],
     [ z0, z0, n1, z0 ],
     [ z0, z0, z0, n1 ]]

func gammaDot*(v: Vec4): Mat4 =
  ## `gamma.v = sum_mu v_mu gamma_mu`
  for mu in 0..<nDim:
    for i in 0..3:
      for j in 0..3:
        result[i][j] = result[i][j] + v[mu]*gammaMat[mu][i][j]

# --- small dense complex linear algebra ------------------------------------

func mmul*[N: static int](a, b: array[N, array[N, Complex64]]):
    array[N, array[N, Complex64]] =
  for i in 0..<N:
    for j in 0..<N:
      var s = z0
      for k in 0..<N: s = s + a[i][k]*b[k][j]
      result[i][j] = s

func adj*[N: static int](a: array[N, array[N, Complex64]]):
    array[N, array[N, Complex64]] =
  for i in 0..<N:
    for j in 0..<N: result[i][j] = conjugate(a[j][i])

func tr*[N: static int](a: array[N, array[N, Complex64]]): Complex64 =
  result = z0
  for i in 0..<N: result = result + a[i][i]

func maxAbsDiff*[N: static int](a, b: array[N, array[N, Complex64]]): float =
  for i in 0..<N:
    for j in 0..<N: result = max(result, abs(a[i][j] - b[i][j]))

proc det*[N: static int](a0: array[N, array[N, Complex64]]): Complex64 =
  ## complex LU with partial pivoting
  var a = a0
  var d = z1
  for k in 0..<N:
    var pv = k
    for i in k+1..<N:
      if abs(a[i][k]) > abs(a[pv][k]): pv = i
    if pv != k:
      swap(a[pv], a[k])
      d = -d
    if a[k][k] == z0: return z0
    d = d*a[k][k]
    for i in k+1..<N:
      let f = a[i][k]/a[k][k]
      for j in k..<N: a[i][j] = a[i][j] - f*a[k][j]
  d

# ---------------------------------------------------------------------------
# the operator as a matrix
# ---------------------------------------------------------------------------

func freeD4*(p: Vec4; lat: HcLat; r = 1.0; m = 0.0): Mat4 =
  ## `D(p) = (m + M(p)) 1 + i gamma.K(p)` at a *true* momentum `p`.
  let
    mm = m + freeMClosed(p, lat, r)
    k = freeKClosed(p, lat)
  result = gammaDot(k)
  for i in 0..3:
    for j in 0..3:
      result[i][j] = zi*result[i][j]
    result[i][i] = result[i][i] + mm

func freeD4Direct*(p: Vec4; lat: HcLat; r = 1.0; m = 0.0): Mat4 =
  ## The same operator assembled the way the position-space code does it:
  ## `D = (m + 4r) + (1/(2 or 6)) sum_i (gamma.n_i - r) e^{i p.n_i}`.
  let
    dirs = dirsOf(lat)
    f = pref(lat)
  for i in 0..3: result[i][i] = complex64(m + 4.0*r)
  for n in dirs:
    let
      gn = gammaDot(n)
      ph = complex64(cos(dot(p, n)), sin(dot(p, n)))
    for i in 0..3:
      for j in 0..3:
        var e = gn[i][j]
        if i == j: e = e - r
        result[i][j] = result[i][j] + f*(e*ph)

func freeD8*(k: Vec4; r = 1.0; m = 0.0): Mat8 =
  ## The 16-cell operator at a *cell* momentum `k`, as an 8x8 matrix in
  ## (spin x sublattice); index = 2*spin + sub, sub = 0 (A) / 1 (B).
  ##
  ##   D_AA = D_BB = (m+4r) + (1/6) sum_mu [ 2 i gamma_mu sin k_mu - 2 r cos k_mu ]
  ##   D_BA = (1/6) sum_delta (gamma.d(delta) - r) e^{ i k.delta}
  ##   D_AB = (1/6) sum_delta (gamma.d(delta) - r) e^{-i k.deltabar},  deltabar = delta xor 15
  var daa: Mat4
  for mu in 0..<nDim:
    let
      s = sin(k[mu])
      c = cos(k[mu])
    for i in 0..3:
      for j in 0..3:
        daa[i][j] = daa[i][j] + (2.0*s/6.0)*(zi*gammaMat[mu][i][j])
      daa[i][i] = daa[i][i] - complex64(2.0*r*c/6.0)
  for i in 0..3: daa[i][i] = daa[i][i] + complex64(m + 4.0*r)

  var dba, dab: Mat4
  for delta in 0..15:
    let
      d = toFloat(deltaVec(delta))
      gd = gammaDot(d)
      dbar = delta xor 15
    var kd = 0.0
    var kdb = 0.0
    for mu in 0..<nDim:
      if ((delta shr mu) and 1) == 1: kd += k[mu]
      if ((dbar shr mu) and 1) == 1: kdb += k[mu]
    let
      ph1 = complex64(cos(kd), sin(kd))
      ph2 = complex64(cos(kdb), -sin(kdb))
    for i in 0..3:
      for j in 0..3:
        var e = gd[i][j]
        if i == j: e = e - r
        dba[i][j] = dba[i][j] + (1.0/6.0)*(e*ph1)
        dab[i][j] = dab[i][j] + (1.0/6.0)*(e*ph2)

  for i in 0..3:
    for j in 0..3:
      result[2*i+0][2*j+0] = daa[i][j]
      result[2*i+1][2*j+1] = daa[i][j]
      result[2*i+0][2*j+1] = dab[i][j]
      result[2*i+1][2*j+0] = dba[i][j]

func gamma5x8*(): Mat8 =
  ## `gamma5 (x) 1` in the (spin x sublattice) ordering used by `freeD8`
  for i in 0..3:
    for j in 0..3:
      result[2*i+0][2*j+0] = gamma5Mat[i][j]
      result[2*i+1][2*j+1] = gamma5Mat[i][j]

# ---------------------------------------------------------------------------
# momentum enumeration
# ---------------------------------------------------------------------------

iterator momenta*(lat: HcLat; ns, nt: int; antiperiodicTime = true): Vec4 =
  ## All `Ns^3*Nt` (cubic) / `2*Ns^3*Nt` (16-cell) momenta of an `Ns^3 x Nt`
  ## lattice.  Spatial `p_j = 2 pi n_j / Ns`, time
  ## `p_3 = 2 pi (n_3 + 1/2)/Nt` (antiperiodic) or `2 pi n_3/Nt` (periodic).
  ## For the 16-cell the range of `n_0` is doubled, `p_0 in [0, 4 pi)`.
  let
    n0max = if lat == lCubic: ns else: 2*ns
    sh = if antiperiodicTime: 0.5 else: 0.0
  var p: Vec4
  for n0 in 0..<n0max:
    p[0] = 2.0*PI*n0.float/ns.float
    for n1 in 0..<ns:
      p[1] = 2.0*PI*n1.float/ns.float
      for n2 in 0..<ns:
        p[2] = 2.0*PI*n2.float/ns.float
        for n3 in 0..<nt:
          p[3] = 2.0*PI*(n3.float + sh)/nt.float
          yield p

iterator momentaAlt*(lat: HcLat; ns, nt: int; antiperiodicTime = true): Vec4 =
  ## The same complete set of momenta, but with the *time* range doubled
  ## instead of `p_0` (16-cell only; identical to `momenta` for the cubic
  ## lattice).  Both are valid sets of coset representatives of the momentum
  ## lattice modulo the reciprocal lattice `2 pi D4`; summing any periodic
  ## function of `p` over either must give the same answer.
  let
    n3max = if lat == lCubic: nt else: 2*nt
    sh = if antiperiodicTime: 0.5 else: 0.0
  var p: Vec4
  for n0 in 0..<ns:
    p[0] = 2.0*PI*n0.float/ns.float
    for n1 in 0..<ns:
      p[1] = 2.0*PI*n1.float/ns.float
      for n2 in 0..<ns:
        p[2] = 2.0*PI*n2.float/ns.float
        for n3 in 0..<n3max:
          p[3] = 2.0*PI*(n3.float + sh)/nt.float
          yield p

iterator cellMomenta*(ns, nt: int; antiperiodicTime = true): Vec4 =
  ## The `Ns^3*Nt` *cell* momenta at which `freeD8` is the full 8x8 operator.
  let sh = if antiperiodicTime: 0.5 else: 0.0
  var k: Vec4
  for n0 in 0..<ns:
    k[0] = 2.0*PI*n0.float/ns.float
    for n1 in 0..<ns:
      k[1] = 2.0*PI*n1.float/ns.float
      for n2 in 0..<ns:
        k[2] = 2.0*PI*n2.float/ns.float
        for n3 in 0..<nt:
          k[3] = 2.0*PI*(n3.float + sh)/nt.float
          yield k

func branchMomenta*(k: Vec4): array[2, Vec4] =
  ## The two true momenta belonging to one cell momentum: `k` and `k + 2 pi e_3`
  ## (they differ by an element of `2 pi Z^4` that is *not* in `2 pi D4`).
  result[0] = k
  result[1] = k
  result[1][3] += 2.0*PI

# ---------------------------------------------------------------------------
# transfer-matrix modes  (used by the free-energy / pressure calculation)
# ---------------------------------------------------------------------------
#
# At fixed spatial momentum `q`, write `x = cos(theta)` with
#   cubic     theta = p_3      (period 2 pi)
#   honeycomb theta = p_3/2    (period 2 pi; p_3 has period 4 pi)
# Then `M(p)^2 + |K(p)|^2` is a polynomial `P(x)` -- degree 1 for the cubic
# lattice, degree 3 for the honeycomb.  Its roots `x_i = cosh(E_i)` give the
# transfer-matrix eigenvalues `v_i = e^{-E_i}`, `|v_i| <= 1`.
#
# It is numerically much better to solve for `y = 1/x` : the reversed
# polynomial has leading coefficient `P(0) = B0 >= 100/9`, whereas the leading
# coefficient of `P` itself is `(16/9) prod_j cos(q_j/2)`, which vanishes
# whenever some `q_j = pi`.  From `y` one gets `v = y/(1 + sqrt(1-y^2))`
# directly, with the principal square root; that automatically selects
# `|v| <= 1` and has no cancellation.

func polyB*(q: array[3, float]; lat: HcLat): seq[float] =
  ## Coefficients `B_k` of `P(x) = M^2 + |K|^2 = sum_k B_k x^k`, `r = 1`, `m = 0`,
  ## with `x = cos p_3` (cubic) or `x = cos(p_3/2)` (honeycomb).
  if lat == lCubic:
    var ms = 0.0
    var ks = 0.0
    for j in 0..2:
      ms += 1.0 - cos(q[j])
      ks += sin(q[j])^2
    result = @[ms*ms + ks + 2.0*ms + 2.0, -2.0*(ms + 1.0)]
  else:
    var c, s, cq, sq: array[3, float]
    for j in 0..2:
      c[j] = cos(0.5*q[j])
      s[j] = sin(0.5*q[j])
      cq[j] = 2.0*c[j]*c[j] - 1.0
      sq[j] = 2.0*s[j]*c[j]
    let cc = c[0]*c[1]*c[2]
    var alpha = 13.0/3.0
    for j in 0..2: alpha -= cq[j]/3.0
    var a, g: array[3, float]
    for j in 0..2:
      a[j] = 2.0*sq[j]
      g[j] = 8.0*s[j]*c[(j+1) mod 3]*c[(j+2) mod 3]
    var aa, ag, gg = 0.0
    for j in 0..2:
      aa += a[j]*a[j]
      ag += a[j]*g[j]
      gg += g[j]*g[j]
    result = @[alpha*alpha + aa/36.0 + (16.0/9.0)*cc*cc,
               -(16.0/3.0)*alpha*cc + ag/18.0 + (16.0/9.0)*cc,
               (16.0/3.0)*cc*cc - (4.0/3.0)*alpha + gg/36.0 + 4.0/9.0,
               (16.0/9.0)*cc]

func polyEval*(b: openArray[float]; x: float): float =
  result = 0.0
  for i in countdown(b.high, 0): result = result*x + b[i]

proc cubicRoots*(a3, a2, a1, a0: float): array[3, Complex64] =
  ## Roots of `a3 y^3 + a2 y^2 + a1 y + a0` (requires `a3 != 0`), by Cardano /
  ## trigonometric solution followed by Newton polishing.
  let
    a = a2/a3
    b = a1/a3
    c = a0/a3
    q = (a*a - 3.0*b)/9.0
    r = (2.0*a*a*a - 9.0*a*b + 27.0*c)/54.0
    d = r*r - q*q*q
  if d < 0.0:
    let
      th = arccos(max(-1.0, min(1.0, r/sqrt(q*q*q))))
      sq = -2.0*sqrt(q)
    result[0] = complex64(sq*cos(th/3.0) - a/3.0)
    result[1] = complex64(sq*cos((th + 2.0*PI)/3.0) - a/3.0)
    result[2] = complex64(sq*cos((th - 2.0*PI)/3.0) - a/3.0)
  else:
    let sd = sqrt(d)
    var aa = -copySign(pow(abs(r) + sd, 1.0/3.0), r)
    let bb = if aa != 0.0: q/aa else: 0.0
    result[0] = complex64(aa + bb - a/3.0)
    result[1] = complex64(-0.5*(aa + bb) - a/3.0, 0.5*sqrt(3.0)*(aa - bb))
    result[2] = conjugate(result[1])
  for i in 0..2:
    var z = result[i]
    for it in 0..2:
      let
        f = ((a3*z + a2)*z + a1)*z + a0
        fp = (3.0*a3*z + 2.0*a2)*z + a1
      if fp == z0: break
      z = z - f/fp
    result[i] = z

func vFromY(y: Complex64): Complex64 {.inline.} =
  ## `x = 1/y = cosh E`, `v = e^{-E}`.  With the principal square root
  ## `Re sqrt(1-y^2) >= 0`, so `|v| <= 1` automatically.
  y/(1.0 + sqrt(1.0 - y*y))

proc modeV*(q: array[3, float]; lat: HcLat): seq[Complex64] =
  ## Transfer-matrix eigenvalues `v_i = e^{-E_i}` with `|v_i| <= 1`, one per
  ## mode: 1 for the cubic lattice, 3 for the 16-cell.  `E_i` is conjugate to
  ## `theta` (= `p_3` for the cubic lattice, `p_3/2` for the 16-cell), so the
  ## physical energies are `omega_i = E_i` and `2 E_i` respectively.
  let b = polyB(q, lat)
  if lat == lCubic:
    # P(x) = B0 + B1 x, B1 = -2(Ms+1) < 0, root x = -B0/B1 = A/B >= 1.
    # d = 1 - B/A is computed directly to avoid cancellation.
    let
      aa = b[0]
      bb = -b[1]
      d = (aa - bb)/aa            # = (Ms^2+Ks^2)/A, exact and >= 0
      s = sqrt(d*(2.0 - d))
    result = @[complex64((1.0 - d)/(1.0 + s))]
  else:
    let ys = cubicRoots(b[0], b[1], b[2], b[3])   # reversed cubic: y = 1/x
    result = @[vFromY(ys[0]), vFromY(ys[1]), vFromY(ys[2])]

func modeExponent*(lat: HcLat; nt: int): int {.inline.} =
  ## Number of Matsubara points per period of `theta`: `Nt` for the cubic
  ## lattice, `2*Nt` for the 16-cell (where `theta = p_3/2`).
  if lat == lCubic: nt else: 2*nt

proc pressureIntegrand*(q: array[3, float]; lat: HcLat; nt: int): float =
  ## `sum_i ln|1 + v_i^{N'}|` with `N' = modeExponent`.  The free-energy
  ## difference is
  ##   O(Nt) = Nt^3 * < pressureIntegrand >_q ,
  ## the average being `int d^3q/(2 pi)^3` over `[0,2 pi)^3` for *both*
  ## lattices (for the 16-cell the doubling of the BZ has been moved into the
  ## time direction, see `momentaAlt`).
  let np = modeExponent(lat, nt)
  for v in modeV(q, lat):
    var z = complex64(1.0, 0.0)
    for i in 1..np: z = z*v
    result += 0.5*ln((1.0 + z.re)^2 + z.im^2)

const contPressure* = 7.0*PI*PI/720.0
  ## `7 pi^2 / 720` -- the continuum `(p(T)-p(0))/T^4` per fermionic dof.

when isMainModule:
  echo "free Wilson-Dirac operator, r = 1, m = 0"
  for lat in [lCubic, lHoneycomb]:
    let pi4: Vec4 = [PI, PI, PI, PI]
    echo "  ", lat, ": neighbours = ", nNeighbours(lat), ", pref = ", pref(lat)
    echo "    M(pi,pi,pi,pi)     = ", freeM(pi4, lat)
    echo "    max |K| (analytic) = ",
      (if lat == lCubic: absK([0.5*PI, 0.5*PI, 0.5*PI, 0.5*PI], lat)
       else: absK([2.0*arccos(0.5*(sqrt(3.0)-1.0)), 0.0, 0.0, 0.0], lat))
