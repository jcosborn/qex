## Shared value types for the radial (S^2 x R) QED3 code.
##
## Plain value semantics throughout: no QEX Layout/Field, no expression templates.
## Complex64 has the same memory layout as LAPACK's dcomplex, so a seq[Complex64]
## can be passed to zgeev/zheev with a cast.

import std/[complex, math]
export complex

type
  Vec3* = array[3, float]
  Spinor* = array[2, Complex64]           ## two-component Dirac spinor
  Spin* = seq[Spinor]                     ## spinor field, flat index v + nv*t
  Mat2* = array[2, array[2, Complex64]]   ## 2x2 complex matrix

template c*(x, y: float): Complex64 = complex64(x, y)
template cr*(x: float): Complex64 = complex64(x, 0.0)

const
  id2*: Mat2 = [[cr 1.0, cr 0.0], [cr 0.0, cr 1.0]]
  sig1*: Mat2 = [[cr 0.0, cr 1.0], [cr 1.0, cr 0.0]]
  sig2*: Mat2 = [[cr 0.0, c(0.0, -1.0)], [c(0.0, 1.0), cr 0.0]]
  sig3*: Mat2 = [[cr 1.0, cr 0.0], [cr 0.0, cr(-1.0)]]

# --- Vec3 -------------------------------------------------------------------

func dot*(a, b: Vec3): float = a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
func cross*(a, b: Vec3): Vec3 =
  [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]]
func norm*(a: Vec3): float = sqrt(dot(a, a))
func unit*(a: Vec3): Vec3 =
  let s = 1.0/norm(a)
  [s*a[0], s*a[1], s*a[2]]
func `+`*(a, b: Vec3): Vec3 = [a[0]+b[0], a[1]+b[1], a[2]+b[2]]
func `-`*(a, b: Vec3): Vec3 = [a[0]-b[0], a[1]-b[1], a[2]-b[2]]
func `*`*(s: float, a: Vec3): Vec3 = [s*a[0], s*a[1], s*a[2]]

func geodesic*(a, b: Vec3): float =
  ## Arc length on the unit sphere. atan2 form is accurate for small separations.
  arctan2(norm(cross(a, b)), dot(a, b))

func tangent*(a, b: Vec3): Vec3 =
  ## Unit tangent at `a` of the geodesic a -> b.
  let d = b - dot(a, b)*a
  unit d

# --- Mat2 / Spinor ----------------------------------------------------------

func `*`*(m: Mat2, x: Spinor): Spinor =
  [m[0][0]*x[0] + m[0][1]*x[1], m[1][0]*x[0] + m[1][1]*x[1]]

func `*`*(m, n: Mat2): Mat2 =
  for i in 0..1:
    for j in 0..1:
      result[i][j] = m[i][0]*n[0][j] + m[i][1]*n[1][j]

func `*`*(z: Complex64, m: Mat2): Mat2 =
  for i in 0..1:
    for j in 0..1: result[i][j] = z*m[i][j]

func `+`*(m, n: Mat2): Mat2 =
  for i in 0..1:
    for j in 0..1: result[i][j] = m[i][j] + n[i][j]

func `-`*(m, n: Mat2): Mat2 =
  for i in 0..1:
    for j in 0..1: result[i][j] = m[i][j] - n[i][j]

func adj*(m: Mat2): Mat2 =
  for i in 0..1:
    for j in 0..1: result[i][j] = conjugate(m[j][i])

func expIsig3*(w: float): Mat2 =
  ## Omega = exp(i sigma3 w/2) = diag(e^{iw/2}, e^{-iw/2}).  The factor 1/2 is built in.
  let h = 0.5*w
  [[c(cos h, sin h), cr 0.0], [cr 0.0, c(cos h, -sin h)]]

func esig*(e: array[2, float]): Mat2 =
  ## e^a sigma_a for a spatial tangent-frame vector (e^1, e^2).
  [[cr 0.0, c(e[0], -e[1])], [c(e[0], e[1]), cr 0.0]]

func `+`*(x, y: Spinor): Spinor = [x[0]+y[0], x[1]+y[1]]
func `-`*(x, y: Spinor): Spinor = [x[0]-y[0], x[1]-y[1]]
func `*`*(a: float, x: Spinor): Spinor = [a*x[0], a*x[1]]
func `*`*(a: Complex64, x: Spinor): Spinor = [a*x[0], a*x[1]]
func sdot*(x, y: Spinor): Complex64 = conjugate(x[0])*y[0] + conjugate(x[1])*y[1]
func snorm2*(x: Spinor): float = abs2(x[0]) + abs2(x[1])
