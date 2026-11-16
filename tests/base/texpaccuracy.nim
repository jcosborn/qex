import qex
from maths/matexp import expAH
import testutils, expaccuracyRef
import std/[algorithm, math, strformat]

qexInit()

# Accuracy and precision of exp (maths/matrixFunctions: degree 4 Taylor of
# X/2^20, then 20 squarings) and expAH (maths/matexp: X/2^k with
# ||X/2^k||_F <= 1/4, degree 12 Taylor, then k squarings) on anti-Hermitian X,
# traceless for N > 1, against the double-double reference of expaccuracyRef.
# Errors are in units of eps of the element type. Every matrix is evaluated
# once with scalar elements and once in a lane of the lattice SIMD type; the
# float32 kernels receive the float32 rounding of the same matrices and the
# reference is taken from that rounded input.
#
# Range of interest: ||X||_F at the call sites (mean .. max over the links of
# an 8^4 Wilson beta 6 configuration; HMC momenta P from randomTAH have
# E||P||_F^2 = N^2 - 1, so ||P||_F ~ 2.8 and up to 7 on 4 x 32^4 links):
#   gauge fixing, polish steps   eps 0.1..0.2 times the remaining gradient   1e-4 .. 2e-3
#   gradient flow                eps 0.005..0.1, (15/36) eps N ||force||      0.008 .. 0.3
#   HMC momentum update          ||eps P||, eps = c tau/steps = 0.03..0.12    0.09 .. 0.8
#   stout smearing               ||rho N Q||, rho 0.1..0.125                  0.4 .. 0.8
#   warm start                   ||s P||, s = 0.4..0.5                        1 .. 3
# Taken as 1e-4 <= ||X||_F <= 2; the grid extends to 8.
#
# Error model. E = exp(X) is unitary, ||E||_F = sqrt(N). Both kernels form
# E = 1 + r from r = expm1(X/2^k) and k squarings r <- r (r + 2), each rounding
# relative to ||r||. A squaring doubles the error d already in r,
# d -> (2 + r) d + d r ~ 2 E d, and adds its own rounding c eps ||r||, so
#   ||E_k - E||_F ~ c eps sum_j 2^(k-j) ||r_j||_F ~ (p + k) c eps ||X||_F,
# p stages for the polynomial; 2^(k-j) ||r_j|| ~ ||X|| for small X, and once
# ||r|| saturates at 2 sqrt(N) the doubling goes on, so the error stays linear
# in ||X||. Relative to ||E||_F, plus eps/2 from rounding the diagonal of 1 + r:
#   err = ||E_k - E||_F/||E||_F ~ eps/2 + (p + k) c eps ||X||_F/sqrt(N).
# exp: k = 20, p ~ 1 (degree 4 at ||X||/2^20). expAH: k = 0 for ||X|| <= 1/4,
# else k = ceil(log2(4 ||X||)), 5 at ||X|| = 8; p ~ 3 (degree 12 at norm
# 1/8..1/4). The measured maxima follow this with c = 0.15..0.3 (roundings are
# random, not worst case; means are 4 to 5 times smaller): expAH's error is 2 to
# 3 times smaller than exp's above ||X|| = 1/4 and both sit at the eps/4 floor
# below 0.1. N = 1 is the scalar complex exp (libm sin, cos), error below eps/2;
# in SIMD lanes x86 gcc -Ofast calls glibc's libmvec instead, accurate to 4 ulp,
# which is 2 eps for cos in [1/2, 1), so the N = 1 bounds start at 5/2 eps.
#
# Checked per norm band at its upper norm, margin 1.4 or more over the maxima
# measured with clang -Ofast -ffast-math on AVX (2 x 516 samples per norm for N = 3):
#   exp:    err <= (1/2 + 6 ||X||/sqrt(N)) eps, 5/2 in place of 1/2 for N = 1
#   expAH:  err <= (1/2 + 3 ||X||/sqrt(N)) eps, likewise
#   unitarity: U = E + D, U'U - 1 = E'D + D'E + D'D, so ||U'U - 1||_F <= 2 sqrt(N) err + O(D^2)
#   determinant: det U - 1 = tr(E'D) + O(D^2), so |det U - 1| <= N err
#   scalar against SIMD lane: under -ffast-math the compiler orders the scalar
#   and the vector arithmetic differently, so the two results are not bitwise
#   equal (the table shows how often they are); both are measured against the
#   reference and their distance is at most the sum of their errors:
#   ||scalar - lane||_F/||E||_F <= 2 err.
# Inside the range of interest expAH's maximum and mean error must not exceed
# exp's in any band, up to eps/16 for ties at the floor, where both kernels
# round the same 1 + r.

type
  CM[N: static int, T] = MatrixArray[N, N, ComplexType[T]]
  Stat = object
    fwd: seq[float]          # ||f(X) - E||_F/||E||_F in eps, scalar and SIMD evaluations
    unit, det, simd: float   # maxima in eps: ||U'U - 1||_F, |det U - 1|, ||scalar - lane||_F/||E||_F
    same: int                # SIMD lanes bitwise equal to the scalar result
  Band = object
    lo, hi: float            # norms in (lo, hi]
    n: int                   # matrices
    e, a: Stat               # exp, expAH

const
  roi = (lo: 1e-4, hi: 2.0)
  edges = [1e-6, 1e-4, 1e-3, 1e-2, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0]
  tie = 1.0/16                # allowance for the expAH against exp comparison at the floor

# ---- sampling ----

var rs: RngMilc6
rs.seed(13, 987654321)

proc toF[N: static int](x: CM[N, float64]): FMat[N] =
  for i in 0..<N:
    for j in 0..<N: result[i][j] = (x[i,j].re, x[i,j].im)

proc unit[N: static int](x: CM[N, float64]): FMat[N] =
  var u = x
  u *= 1.0/sqrt(norm2(u))
  toF(u)

proc tah(N: static int): CM[N, float64] =
  var g: CM[N, float64]
  gaussian(g, rs)
  result.projectTAH g

proc cases(N: static int, nrand: int): seq[FMat[N]] =
  ## Unit Frobenius norm, anti-Hermitian (traceless for N > 1) by construction:
  ## the structured spectra first, then nrand random matrices.
  when N == 1:
    for s in [1.0, -1.0]:
      var m: FMat[1]
      m[0][0] = (0.0, s)
      result.add m
  else:
    var d: CM[N, float64]                 # diagonal i (1, .., N-1, -N(N-1)/2)
    d := 0
    for i in 0..<N-1: d[i,i].im := float(i+1)
    d[N-1,N-1].im := -float(N*(N-1) div 2)
    result.add unit(d)
    when N == 3:
      var b: CM[3, float64]               # one SU(2) block, third eigenvalue 0
      b := 0
      b[0,0].im := 0.3
      b[1,1].im := -0.3
      b[0,1].re := 0.4
      b[0,1].im := 0.5
      b[1,0].re := -0.4
      b[1,0].im := 0.5
      result.add unit(b)
      for delta in [0.0, 1e-5]:           # spectrum i (1+delta, 1-delta, -2), rotated
        var s, m: CM[3, float64]
        s := 0
        s[0,0].im := 1.0 + delta
        s[1,1].im := 1.0 - delta
        s[2,2].im := -2.0
        let v = expAH(tah(3))
        m.projectTAH(v*s*v.adj)
        result.add unit(m)
    for c in 0..<nrand: result.add unit(tah(N))

proc grid(): seq[float] =
  ## 8 norms per decade inside the range of interest, 2 per decade below it,
  ## 4 per decade above up to 8, and both sides of the expAH branch points
  ## 2^k/4 (k squarings below, k + 1 above).
  for j in -12 .. -9: result.add pow(10.0, j/2)      # 1e-6 .. 3e-5
  for j in -32 .. 2: result.add pow(10.0, j/8)       # 1e-4 .. 1.78
  for j in 2 .. 3: result.add pow(10.0, j/4)         # 3.16, 5.62
  result.add 8.0
  for k in 0..4:
    let b = pow(2.0, float(k))/4
    result.add b*(1 - 1e-3)
    result.add b*(1 + 1e-3)
  result.sort

proc scaled[N: static int](c: FMat[N], s: float): FMat[N] =
  for i in 0..<N:
    for j in 0..<N: result[i][j] = (s*c[i][j].re, s*c[i][j].im)

# ---- evaluation in the element type T, scalar or SIMD ----

proc pack[N: static int](xs: openArray[FMat[N]], T: typedesc): CM[N, T] =
  ## Lane k holds xs[k] rounded to numberType(T).
  for i in 0..<N:
    for j in 0..<N:
      when T is SomeFloat:
        result[i,j].re := T(xs[0][i][j].re)
        result[i,j].im := T(xs[0][i][j].im)
      else:
        var re, im: array[simdLength(T), numberType(T)]
        for k in 0..<re.len:
          re[k] = numberType(T)(xs[k][i][j].re)
          im[k] = numberType(T)(xs[k][i][j].im)
        result[i,j].re := re
        result[i,j].im := im

proc lane(x: auto, k: int): auto =
  ## Lane k of a matrix with scalar or SIMD elements, as float64 (a scalar is its own lane 0).
  var r: FMat[x.nrows]
  for i in 0..<x.nrows:
    for j in 0..<x.ncols: r[i][j] = (float64(x[i,j].re[k]), float64(x[i,j].im[k]))
  r

proc rounded[N: static int](x: FMat[N], S: typedesc): FMat[N] = lane(pack([x], S), 0)

proc eval[N: static int](xs: openArray[FMat[N]], T: typedesc): tuple[e, a: seq[FMat[N]]] =
  ## exp and expAH of every xs[i] rounded to numberType(T), evaluated in T, as float64.
  ## Spare lanes of the last vector repeat the last matrix.
  const L = simdLength(T)
  var i = 0
  while i < xs.len:
    var vs: array[L, FMat[N]]
    for k in 0..<L: vs[k] = xs[min(i + k, xs.len - 1)]
    let x = pack(vs, T)
    let e = exp(x)
    let a = expAH(x)
    for k in 0..<min(L, xs.len - i):
      result.e.add lane(e, k)
      result.a.add lane(a, k)
    i += L

# ---- measurement ----

proc dist[N: static int](a, b: FMat[N]): float =
  ## ||a - b||_F/sqrt(N) in float64.
  for i in 0..<N:
    for j in 0..<N:
      result += (a[i][j].re - b[i][j].re)^2 + (a[i][j].im - b[i][j].im)^2
  sqrt(result/float(N))

proc add(s: var Stat, scalar, simd: tuple[fwd, unit, det: float], d, eps: float) =
  for m in [scalar, simd]:
    s.fwd.add m.fwd/eps
    s.unit = max(s.unit, m.unit/eps)
    s.det = max(s.det, m.det/eps)
  s.simd = max(s.simd, d/eps)
  if d == 0: inc s.same

proc bandOf(x: float): int =
  for i in 1..<edges.len:
    if x <= edges[i]: return i - 1
  edges.len - 2

proc run[N: static int](cs: seq[FMat[N]], T: typedesc): seq[Band] =
  ## T is the SIMD element type of the lattice fields, numberType(T) the scalar type.
  type S = numberType(T)
  let eps = float(epsilon(S))
  for i in 1..<edges.len: result.add Band(lo: edges[i-1], hi: edges[i])
  for nrm in grid():
    let b = bandOf(nrm)
    var xs: seq[FMat[N]]
    for c in cs: xs.add scaled(c, nrm)
    let s = eval(xs, S)
    let v = eval(xs, T)
    for i, x in xs:
      let r = expRef(rounded(x, S))
      result[b].e.add(measure(s.e[i], r), measure(v.e[i], r), dist(s.e[i], v.e[i]), eps)
      result[b].a.add(measure(s.a[i], r), measure(v.a[i], r), dist(s.a[i], v.a[i]), eps)
      inc result[b].n

# ---- bounds and report ----

proc floorEps(N: int): float = (if N == 1: 2.5 else: 0.5)
proc tolExp(N: int, nrm: float): float = floorEps(N) + 6.0*nrm/sqrt(float(N))
proc tolAH(N: int, nrm: float): float = floorEps(N) + 3.0*nrm/sqrt(float(N))
proc mean(xs: seq[float]): float = sum(xs)/float(xs.len)
proc pct(xs: seq[float], p: float): float = sorted(xs)[int(p*float(xs.len - 1))]

proc report(name: string, N: int, bs: seq[Band]) =
  echo &"\n{name}, N = {N}: errors in units of eps from the scalar and the SIMD evaluation of n matrices"
  echo "      ||X||_F           n    exp max   mean    p99  expAH max   mean    p99" &
       "  |U'U-1| exp  expAH  |detU-1| exp  expAH  scalar-SIMD exp  expAH  bitwise equal"
  for b in bs:
    echo &"  {b.lo:>7g} .. {b.hi:<6g} {b.n:5}  {max(b.e.fwd):7.2f} {mean(b.e.fwd):6.2f} {pct(b.e.fwd, 0.99):6.2f}" &
         &"  {max(b.a.fwd):7.2f} {mean(b.a.fwd):6.2f} {pct(b.a.fwd, 0.99):6.2f}" &
         &"  {b.e.unit:10.2f} {b.a.unit:6.2f}  {b.e.det:10.2f} {b.a.det:6.2f}" &
         &"  {b.e.simd:12.2f} {b.a.simd:6.2f}  {100*(b.e.same + b.a.same) div (2*b.n):12}%"

proc verdict(name: string, N: int, b: Band): tuple[ok: bool, line: string] =
  ## The bounds of the header at the band's upper norm; the expAH against exp
  ## comparison only inside the range of interest.
  let te = tolExp(N, b.hi)
  let ta = tolAH(N, b.hi)
  let emax = max(b.e.fwd)
  let amax = max(b.a.fwd)
  var bad: seq[string]
  if emax > te: bad.add "exp error"
  if amax > ta: bad.add "expAH error"
  if b.e.unit > 2*sqrt(float(N))*te or b.a.unit > 2*sqrt(float(N))*ta: bad.add "unitarity"
  if b.e.det > float(N)*te or b.a.det > float(N)*ta: bad.add "determinant"
  if b.e.simd > 2*te or b.a.simd > 2*ta: bad.add "scalar-SIMD"
  var s = &"{name} N={N}, norms {b.lo:g}..{b.hi:g}: expAH max error {amax:.2f} eps (bound {ta:.2f}), exp {emax:.2f} eps (bound {te:.2f})"
  if b.hi > roi.lo and b.lo < roi.hi:
    if amax <= emax + tie and mean(b.a.fwd) <= mean(b.e.fwd) + tie: s &= ", expAH <= exp"
    else: bad.add "expAH > exp"
  s &= (if bad.len == 0: ": ok" else: ": FAIL " & $bad)
  (bad.len == 0, s)

proc overall(name: string, N: int, bs: seq[Band]): string =
  var e, a: seq[float]
  for b in bs:
    if b.hi > roi.lo and b.lo < roi.hi:
      e.add b.e.fwd
      a.add b.a.fwd
  &"{name} N={N}, range of interest {roi.lo:g}..{roi.hi:g}: expAH max error {max(a):.2f} eps, mean {mean(a):.3f}; exp max {max(e):.2f} eps, mean {mean(e):.3f}"

# ---- tests ----

let c3 = cases(3, 512)
let c2 = cases(2, 512)
let c1 = cases(1, 0)

suite "exp and expAH on anti-Hermitian matrices":
  test "reference arithmetic keeps the error terms":
    check exact()
  template run(N: static int, T: typedesc) =
    test "N=" & $N & " " & $numberType(T) & " and " & $T:
      let name = $numberType(T) & ", scalar and SIMD (" & $simdLength(T) & " lanes)"
      let bs = run((when N == 3: c3 elif N == 2: c2 else: c1), T)
      report(name, N, bs)
      for b in bs:
        let (ok, line) = verdict($numberType(T), N, b)
        echo line
        check ok
      echo overall($numberType(T), N, bs)
  run(3, Dvec0)
  run(3, Svec0)
  run(2, Dvec0)
  run(2, Svec0)
  run(1, Dvec0)
  run(1, Svec0)

qexFinalize()
