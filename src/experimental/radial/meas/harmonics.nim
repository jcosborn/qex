## Real spherical harmonics on the icosahedral S^2 lattice, the two quadratures
## (site and face), and the icosahedral rotation group  (WP-I).
##
## Normative reference: doc/07-observables.md sections 1.3 and 2,
## doc/04-interfaces.md section 14.
##
## Y_lm are the REAL harmonics (doc/07 section 1.3: everything stays real),
## orthonormal on the unit sphere, sum_m Y_lm(v)^2 = (2l+1)/(4pi).  They are
## hard-coded polynomials in the components of the unit vector for l <= 4 --
## QEX has no Legendre code and the associated recurrences buy nothing at
## this degree.  `legendreP` is provided as the independent cross-check
## (addition theorem) the tests use.
##
## Quadratures:
##   siteProject: sum_y A_y Y_lm(pos_y) x_y     -- the lattice measure of doc/07 1.3
##   faceProject: sum_f Y_lm(cc_f) x_f          -- NO area weight: for J^t_lat =
##     Theta/A_tri the area cancels (doc/07 section 4.1), so the caller passes
##     the flux Theta_f itself
##
## Icosahedral protection (doc/07 section 2): every orbit of the icosahedral
## group is a spherical 5-design (the first nontrivial I_h-invariant harmonic
## is l = 6), so at L = 1 both quadratures integrate Y_lm Y_l'm' EXACTLY for
## l + l' <= 5; at any L, Schur's lemma makes each l <= 2 Gram block exactly
## proportional to the identity and every block between inequivalent irreps
## exactly zero.  `icosaGroup` supplies the 60 rotations and the site/face
## permutations they induce, for the exact single-configuration degeneracy
## tests.

import std/math
import ../core/geom
export geom

func legendreP*(l: int, x: float): float =
  ## Legendre polynomial by the three-term recurrence; the independent oracle
  ## for `ylm` through the addition theorem
  ##   sum_m Y_lm(u) Y_lm(v) = (2l+1)/(4pi) P_l(u.v).
  if l == 0: return 1.0
  if l == 1: return x
  var
    pm = 1.0
    p = x
  for k in 2..l:
    let pn = (float(2*k - 1)*x*p - float(k - 1)*pm)/float(k)
    pm = p
    p = pn
  p

func ylm*(l, m: int, v: Vec3): float =
  ## Real orthonormal spherical harmonic, l <= 4, -l <= m <= l, `v` a unit
  ## vector.  Standard real forms: m > 0 ~ cos(m phi), m < 0 ~ sin(|m| phi).
  let
    x = v[0]
    y = v[1]
    z = v[2]
  case l
  of 0:
    result = 0.5/sqrt(PI)
  of 1:
    let c = sqrt(3.0/(4.0*PI))
    case m
    of -1: result = c*y
    of 0: result = c*z
    of 1: result = c*x
    else: doAssert false, "ylm: bad m"
  of 2:
    case m
    of -2: result = 0.5*sqrt(15.0/PI)*x*y
    of -1: result = 0.5*sqrt(15.0/PI)*y*z
    of 0: result = 0.25*sqrt(5.0/PI)*(3.0*z*z - 1.0)
    of 1: result = 0.5*sqrt(15.0/PI)*x*z
    of 2: result = 0.25*sqrt(15.0/PI)*(x*x - y*y)
    else: doAssert false, "ylm: bad m"
  of 3:
    case m
    of -3: result = 0.25*sqrt(35.0/(2.0*PI))*y*(3.0*x*x - y*y)
    of -2: result = 0.5*sqrt(105.0/PI)*x*y*z
    of -1: result = 0.25*sqrt(21.0/(2.0*PI))*y*(5.0*z*z - 1.0)
    of 0: result = 0.25*sqrt(7.0/PI)*z*(5.0*z*z - 3.0)
    of 1: result = 0.25*sqrt(21.0/(2.0*PI))*x*(5.0*z*z - 1.0)
    of 2: result = 0.25*sqrt(105.0/PI)*(x*x - y*y)*z
    of 3: result = 0.25*sqrt(35.0/(2.0*PI))*x*(x*x - 3.0*y*y)
    else: doAssert false, "ylm: bad m"
  of 4:
    case m
    of -4: result = 0.75*sqrt(35.0/PI)*x*y*(x*x - y*y)
    of -3: result = 0.75*sqrt(35.0/(2.0*PI))*y*(3.0*x*x - y*y)*z
    of -2: result = 0.75*sqrt(5.0/PI)*x*y*(7.0*z*z - 1.0)
    of -1: result = 0.75*sqrt(5.0/(2.0*PI))*y*z*(7.0*z*z - 3.0)
    of 0: result = (3.0/16.0)*sqrt(1.0/PI)*(35.0*z*z*z*z - 30.0*z*z + 3.0)
    of 1: result = 0.75*sqrt(5.0/(2.0*PI))*x*z*(7.0*z*z - 3.0)
    of 2: result = (3.0/8.0)*sqrt(5.0/PI)*(x*x - y*y)*(7.0*z*z - 1.0)
    of 3: result = 0.75*sqrt(35.0/(2.0*PI))*x*(x*x - 3.0*y*y)*z
    of 4: result = (3.0/16.0)*sqrt(35.0/PI)*(x*x*x*x - 6.0*x*x*y*y + y*y*y*y)
    else: doAssert false, "ylm: bad m"
  else:
    doAssert false, "ylm: only l <= 4 is coded"

func siteProject*(sph: Sphere, x: openArray[float], l, m: int): float =
  ## sum_y A_y Y_lm(pos_y) x_y over one time slice (`x` has length sph.nv).
  ## A_y is the dual-cell area, the lattice measure of doc/07 section 1.3.
  for y in 0..<sph.nv:
    result += sph.area[y]*ylm(l, m, sph.pos[y])*x[y]

func faceProject*(sph: Sphere, x: openArray[float], l, m: int): float =
  ## sum_f Y_lm(cc_f) x_f over one time slice (`x` has length sph.nf).
  ## No area weight: for x = Theta_f this is the (V.12) projection
  ## sum_f A_f Y_lm (Theta_f/A_f) with the area cancelled (doc/07 section 4.1).
  for f in 0..<sph.nf:
    result += ylm(l, m, sph.faces[f].cc)*x[f]

func siteGram*(sph: Sphere, l1, l2: int): seq[float] =
  ## G[(l1 block) x (l2 block)] with G_{mm'} = sum_y A_y Y_{l1 m} Y_{l2 m'},
  ## row-major (2l1+1) x (2l2+1).  Equals delta_{l1 l2} delta_{mm'} up to the
  ## quadrature error; the identity-proportionality within each l <= 2 block
  ## and the vanishing of blocks between inequivalent I_h irreps are EXACT.
  let
    n1 = 2*l1 + 1
    n2 = 2*l2 + 1
  result = newSeq[float](n1*n2)
  for y in 0..<sph.nv:
    let
      a = sph.area[y]
      p = sph.pos[y]
    for i in 0..<n1:
      let yi = ylm(l1, i - l1, p)
      for j in 0..<n2:
        result[i*n2 + j] += a*yi*ylm(l2, j - l2, p)

func faceGram*(sph: Sphere, l1, l2: int): seq[float] =
  ## Same Gram on the face quadrature, sum_f A_tri Y Y' at the dual points.
  let
    n1 = 2*l1 + 1
    n2 = 2*l2 + 1
  result = newSeq[float](n1*n2)
  for f in 0..<sph.nf:
    let
      a = sph.faces[f].area
      p = sph.faces[f].cc
    for i in 0..<n1:
      let yi = ylm(l1, i - l1, p)
      for j in 0..<n2:
        result[i*n2 + j] += a*yi*ylm(l2, j - l2, p)

# --- icosahedral rotation group ----------------------------------------------

type IcosaGroup* = object
  ## The 60 rotations of the icosahedral rotation group I and the permutations
  ## they induce on the lattice.  `rot[g]` is the matrix as three ROWS, so the
  ## image of v is [dot(rot[g][0],v), dot(rot[g][1],v), dot(rot[g][2],v)].
  ## Built numerically from the 5-fold axes; every permutation is verified to
  ## be a bijection with |R pos - pos[perm]| below `tol` (the residual is the
  ## measured equivariance of the lattice construction and is reported).
  rot*: seq[array[3, Vec3]]
  site*: seq[seq[int]]       ## site[g][y] = image site of y under rotation g
  face*: seq[seq[int]]       ## face[g][f] = image face (rotations preserve the
                             ## CCW orientation, so fluxes map with sign +1)
  siteDev*: float            ## max |R pos[y] - pos[site[g][y]]| over all g, y
  faceDev*: float            ## max |R cc[f] - cc[face[g][f]]| over all g, f

func rotate(r: array[3, Vec3], v: Vec3): Vec3 =
  [dot(r[0], v), dot(r[1], v), dot(r[2], v)]

func nearest(pts: seq[Vec3], v: Vec3): tuple[i: int, d: float] =
  result.d = 1e300
  for k in 0..<pts.len:
    let d = norm(pts[k] - v)
    if d < result.d:
      result.d = d
      result.i = k

func frameAt(a, b: Vec3): array[3, Vec3] =
  ## Orthonormal frame (columns as rows here): a, the tangent at a toward b,
  ## and their cross product.
  let t = unit(b - dot(a, b)*a)
  [a, t, cross(a, t)]

proc icosaGroup*(sph: Sphere, tol = 1e-6): IcosaGroup =
  ## All 60 rotations: the reference 5-fold site and one of its neighbours can
  ## go to any (5-fold site, neighbour) pair, and that fixes the rotation.
  var five: seq[int]
  for y in 0..<sph.nv:
    if sph.nbr[y].len == 5: five.add y
  doAssert five.len == 12, "icosaGroup: expected 12 five-fold sites"
  let
    a0 = five[0]
    b0 = sph.nbr[a0][0]
    f0 = frameAt(sph.pos[a0], sph.pos[b0])
  var ccs = newSeq[Vec3](sph.nf)
  for f in 0..<sph.nf: ccs[f] = sph.faces[f].cc
  for a in five:
    for b in sph.nbr[a]:
      let f1 = frameAt(sph.pos[a], sph.pos[b])
      # R = f1 * f0^T: rows of R are sum_k f1[k][i] f0[k][j] assembled per row.
      var r: array[3, Vec3]
      for i in 0..2:
        for j in 0..2:
          r[i][j] = f1[0][i]*f0[0][j] + f1[1][i]*f0[1][j] + f1[2][i]*f0[2][j]
      var
        sp = newSeq[int](sph.nv)
        fp = newSeq[int](sph.nf)
        used = newSeq[bool](sph.nv)
        usedF = newSeq[bool](sph.nf)
      for y in 0..<sph.nv:
        let (i, d) = nearest(sph.pos, rotate(r, sph.pos[y]))
        doAssert d < tol, "icosaGroup: site set not closed under rotation"
        doAssert not used[i], "icosaGroup: site map not a bijection"
        used[i] = true
        sp[y] = i
        result.siteDev = max(result.siteDev, d)
      for f in 0..<sph.nf:
        let (i, d) = nearest(ccs, rotate(r, ccs[f]))
        doAssert d < tol, "icosaGroup: face set not closed under rotation"
        doAssert not usedF[i], "icosaGroup: face map not a bijection"
        usedF[i] = true
        fp[f] = i
        result.faceDev = max(result.faceDev, d)
      result.rot.add r
      result.site.add sp
      result.face.add fp
  doAssert result.rot.len == 60
