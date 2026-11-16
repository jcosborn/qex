## Double-double reference for the complex matrix exponential, used by
## texpaccuracy: hi + lo pairs built from the error free transforms twoSum and
## Dekker's twoProd, about 2^-104 relative rounding per operation. The transforms
## need IEEE semantics; the production build uses -ffast-math, which may drop
## the error terms, and GCC keeps -ffp-contract=fast without it, fusing the
## Dekker split into FMAs. This module alone is compiled without both, and
## `exact` checks at run time that the terms survive.
{.localPassC: "-fno-fast-math -ffp-contract=off".}
import std/math

type
  DD = object            # hi + lo, |lo| <= ulp(hi)/2
    hi, lo: float64
  CDD = object
    re, im: DD
  Mat[N: static int] = array[N, array[N, CDD]]
  FMat*[N: static int] = array[N, array[N, tuple[re, im: float64]]]

proc twoSum(a, b: float64): DD =
  let s = a + b
  let bb = s - a
  DD(hi: s, lo: (a - (s - bb)) + (b - bb))

proc quick(a, b: float64): DD =    # |a| >= |b|
  let s = a + b
  DD(hi: s, lo: b - (s - a))

proc twoProd(a, b: float64): DD =
  # Dekker split by 2^27+1: a = ah + al with 26 bit halves, a*b exact as sum.
  let p = a*b
  let ca = 134217729.0*a
  let ah = ca - (ca - a)
  let al = a - ah
  let cb = 134217729.0*b
  let bh = cb - (cb - b)
  let bl = b - bh
  DD(hi: p, lo: ((ah*bh - p) + ah*bl + al*bh) + al*bl)

proc `+`(a, b: DD): DD =
  let s = twoSum(a.hi, b.hi)
  let t = twoSum(a.lo, b.lo)
  let u = quick(s.hi, s.lo + t.hi)
  quick(u.hi, u.lo + t.lo)
proc `-`(a: DD): DD = DD(hi: -a.hi, lo: -a.lo)
proc `-`(a, b: DD): DD = a + (-b)
proc `*`(a, b: DD): DD =
  let p = twoProd(a.hi, b.hi)
  quick(p.hi, p.lo + (a.hi*b.lo + a.lo*b.hi))
proc `/`(a: DD, b: float64): DD =
  # q1 = a/b rounded; the remainder a - q1 b is exact in dd, so q = q1 + r/b.
  let q1 = a.hi / b
  let r = a - twoProd(q1, b)
  quick(q1, (r.hi + r.lo) / b)

var operands {.volatile.} = [1e-20, 1.0 + 1.0/1073741824.0]
  # read at run time, so `exact` sees the arithmetic of the reference

proc exact*(): bool =
  ## False when the compiler dropped or fused the error terms: twoSum keeps
  ## 1e-20 next to 1, and with e = 2^-30 the dd cube (1 + e)^3 is exactly
  ## (1 + 3e) + (3e^2 + e^3).
  const e = 1.0/1073741824.0
  let s = twoSum(1.0, operands[0])
  let x = DD(hi: operands[1])
  let c = x*x*x
  s.lo == operands[0] and c.hi == 1.0 + 3*e and c.lo == 3*e*e + e*e*e

proc `+`(a, b: CDD): CDD = CDD(re: a.re + b.re, im: a.im + b.im)
proc `-`(a, b: CDD): CDD = CDD(re: a.re - b.re, im: a.im - b.im)
proc `*`(a, b: CDD): CDD =
  CDD(re: a.re*b.re - a.im*b.im, im: a.re*b.im + a.im*b.re)
proc conj(a: CDD): CDD = CDD(re: a.re, im: -a.im)
proc norm2(a: CDD): DD = a.re*a.re + a.im*a.im

proc toMat[N: static int](x: FMat[N]): Mat[N] =
  for i in 0..<N:
    for j in 0..<N:
      result[i][j] = CDD(re: DD(hi: x[i][j].re), im: DD(hi: x[i][j].im))
proc eye(N: static int): Mat[N] =
  for i in 0..<N: result[i][i].re.hi = 1.0
proc `+`[N: static int](a, b: Mat[N]): Mat[N] =
  for i in 0..<N:
    for j in 0..<N: result[i][j] = a[i][j] + b[i][j]
proc `-`[N: static int](a, b: Mat[N]): Mat[N] =
  for i in 0..<N:
    for j in 0..<N: result[i][j] = a[i][j] - b[i][j]
proc `*`[N: static int](a, b: Mat[N]): Mat[N] =
  for i in 0..<N:
    for j in 0..<N:
      for k in 0..<N: result[i][j] = result[i][j] + a[i][k]*b[k][j]
proc `/`[N: static int](a: Mat[N], s: float64): Mat[N] =
  for i in 0..<N:
    for j in 0..<N: result[i][j] = CDD(re: a[i][j].re/s, im: a[i][j].im/s)
proc adj[N: static int](a: Mat[N]): Mat[N] =
  for i in 0..<N:
    for j in 0..<N: result[i][j] = conj(a[j][i])
proc norm2[N: static int](a: Mat[N]): DD =
  for i in 0..<N:
    for j in 0..<N: result = result + norm2(a[i][j])
proc det[N: static int](a: Mat[N]): CDD =
  when N == 1: a[0][0]
  elif N == 2: a[0][0]*a[1][1] - a[0][1]*a[1][0]
  else:
    (a[0][0]*a[1][1] - a[0][1]*a[1][0])*a[2][2] +
    (a[0][2]*a[1][0] - a[0][0]*a[1][2])*a[2][1] +
    (a[0][1]*a[1][2] - a[0][2]*a[1][1])*a[2][0]

proc expRef*[N: static int](x: FMat[N]): Mat[N] =
  ## exp(x) = T(Y)^(2^s), Y = x/2^s with ||Y||_F <= 1/16, T = sum_{j<=20} Y^j/j!.
  ## Truncation: (1/16)^21/21! < 1e-45. Rounding: < 5 N 2^-104 per product and
  ## term, 20 terms, then s squarings each doubling the error present, so
  ## < 2^s (20 + s) 5 N 2^-104 < 1e-26 for s <= 8 (||x||_F <= 16), below 1e-10
  ## eps of float64.
  var n = 0.0
  for i in 0..<N:
    for j in 0..<N: n += x[i][j].re*x[i][j].re + x[i][j].im*x[i][j].im
  var s = 0
  while n > 1.0/256.0:
    n *= 0.25
    inc s
  let y = toMat(x)/float(1 shl s)
  var t = eye(N)
  result = t
  for k in 1..20:
    t = (t*y)/float(k)
    result = result + t
  for i in 0..<s: result = result*result

proc measure*[N: static int](f: FMat[N], r: Mat[N]): tuple[fwd, unit, det: float64] =
  ## fwd = ||f - r||_F/||r||_F, unit = ||f^dag f - 1||_F, det = |det f - det r|,
  ## all in dd arithmetic; det r = exp(tr x) = 1 for traceless x.
  let m = toMat(f)
  result.fwd = sqrt(norm2(m - r).hi/norm2(r).hi)
  result.unit = sqrt(norm2(adj(m)*m - eye(N)).hi)
  result.det = sqrt(norm2(det(m) - det(r)).hi)
