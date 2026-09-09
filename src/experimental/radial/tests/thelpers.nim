## Fixtures and dense oracles shared by the radial test suites: seeded random
## fields, norms, dense spinor application, and the dense incidence matrix
## M = C^T W C of the gauge action with its eigen-decomposition -- the oracle of
## tgauge and tflow, written from doc/02 section 5 independently of
## ops/gaugeact's loops.

import std/[math, complex]
import ../core/dense
import ../core/spinor
import ../ops/gaugeact

export dense

# --- random fixtures ------------------------------------------------------------

proc rng*(sed: int): Threefry4x64 =
  result.seedIndep(sed, 0)

proc randGauge*(l: Lat, r: var Threefry4x64, amp = 1.0): Gauge =
  result = newGauge(l)
  for i in 0..<result.s.len: result.s[i] = amp*r.gaussian
  for i in 0..<result.t.len: result.t[i] = amp*r.gaussian

proc randGauge*(l: Lat, sed: int, amp = 1.0): Gauge =
  var r = rng(sed)
  randGauge(l, r, amp)

proc randSpin*(n: int, r: var Threefry4x64): Spin =
  result = newSpin(n)
  result.gaussian r

proc randSpin*(n, sed: int): Spin =
  var r = rng(sed)
  randSpin(n, r)

proc randReal*(n: int, r: var Threefry4x64, amp = 1.0): seq[float] =
  result = newSeq[float](n)
  for i in 0..<n: result[i] = amp*r.gaussian

proc randReal*(n, sed: int, amp = 1.0): seq[float] =
  var r = rng(sed)
  randReal(n, r, amp)

# --- norms ----------------------------------------------------------------------

func maxAbs*(x: openArray[float]): float =
  for v in x: result = max(result, abs v)

func maxAbs*(a: openArray[Complex64]): float =
  for z in a: result = max(result, abs z)

func maxDiff*(a, b: openArray[float]): float =
  for i in 0..<a.len: result = max(result, abs(a[i] - b[i]))

proc reldiff*(x, y: Spin): float =
  ## |x - y| / |y|
  var d = 0.0
  for i in 0..<x.len:
    for c in 0..1: d += abs2(x[i][c] - y[i][c])
  sqrt(d/norm2(y))

# --- dense spinor algebra ---------------------------------------------------------

proc denseApply*(a: seq[Complex64], nd: int, dst: var Spin, src: Spin, dag = false) =
  ## dst = A src (or A^dag src) for the column-major matrix `a`, spinor flat
  ## index i = 2*site + comp.
  for i in 0..<nd:
    var s = complex64(0.0, 0.0)
    for j in 0..<nd:
      let m = if dag: conjugate(a[j + nd*i]) else: a[i + nd*j]
      s += m*src[j shr 1][j and 1]
    dst[i shr 1][i and 1] = s

proc symEig*(m: seq[float], n: int): seq[float] =
  ## Ascending eigenvalues of a real symmetric matrix given row-major n x n.
  var a = newSeq[Complex64](n*n)
  for i in 0..<n:
    for j in 0..<n:
      a[i + n*j] = complex64(0.5*(m[i*n + j] + m[j*n + i]), 0.0)
  heig(a, n)

# --- dense gauge-action oracle ------------------------------------------------------

type DenseM* = object
  n*: int                    ## number of links
  m*: seq[float]             ## column-major M = C^T W C
  nplaq*: int

proc plaqRow(l: Lat, kind, idx, t: int, row: var seq[float]) =
  ## Dense incidence row of a plaquette.  kind 0 = spatial triangle `idx`,
  ## kind 1 = temporal plaquette of spatial edge `idx`.  Written literally from
  ## doc/02 section 5, independently of gaugeact's loops.
  for i in 0..<row.len: row[i] = 0.0
  let ns = l.sph.ne*l.nt
  if kind == 0:
    let fc = l.sph.faces[idx]
    for i in 0..2:
      row[eIdx(l, fc.e[i], t)] += float(fc.s[i])
  else:
    let ed = l.sph.edges[idx]
    row[eIdx(l, idx, t)] += 1.0
    row[ns + tIdx(l, ed.b, t)] += 1.0
    row[eIdx(l, idx, t+1)] -= 1.0
    row[ns + tIdx(l, ed.a, t)] -= 1.0

proc newDenseM*(l: Lat, b: Beta): DenseM =
  result.n = nlink(l)
  result.m = newSeq[float](result.n*result.n)
  var row = newSeq[float](result.n)
  var nz: seq[int]
  for kind in 0..1:
    let nk = if kind == 0: l.sph.nf else: l.sph.ne
    for idx in 0..<nk:
      let w = if kind == 0: b.face[idx] else: b.edge[idx]
      for t in 0..<l.nt:
        plaqRow(l, kind, idx, t, row)
        nz.setLen 0
        for i in 0..<result.n:
          if row[i] != 0.0: nz.add i
        for i in nz:
          for j in nz:
            result.m[i + result.n*j] += w*row[i]*row[j]
        inc result.nplaq

proc act*(d: DenseM, x: openArray[float]): float =
  ## x^T M x / 2
  for j in 0..<d.n:
    var s = 0.0
    for i in 0..<d.n: s += d.m[i + d.n*j]*x[i]
    result += 0.5*s*x[j]

proc mul*(d: DenseM, x: openArray[float]): seq[float] =
  ## M x
  result = newSeq[float](d.n)
  for j in 0..<d.n:
    let xj = x[j]
    if xj == 0.0: continue
    for i in 0..<d.n: result[i] += d.m[i + d.n*j]*xj

type DenseEig* = tuple[v: seq[Complex64], w: seq[float]]

proc eigen*(d: DenseM): DenseEig =
  ## Real symmetric eigenproblem through the Hermitian LAPACK path (zheev).
  var h = newSeq[Complex64](d.n*d.n)
  for i in 0..<d.n*d.n: h[i] = complex64(d.m[i], 0.0)
  let ew = heig(h, d.n)
  (h, ew)

proc pinv*(d: DenseM, ev: DenseEig, tol: float): seq[float] =
  ## Moore-Penrose pseudo-inverse, kernel dropped.
  result = newSeq[float](d.n*d.n)
  for k in 0..<d.n:
    if abs(ev.w[k]) <= tol: continue
    let iw = 1.0/ev.w[k]
    for j in 0..<d.n:
      let c = iw*ev.v[j + d.n*k].re
      if c == 0.0: continue
      for i in 0..<d.n: result[i + d.n*j] += c*ev.v[i + d.n*k].re

proc expm*(d: DenseM, ev: DenseEig, s: float): seq[float] =
  ## exp(-M s), column major.
  result = newSeq[float](d.n*d.n)
  for k in 0..<d.n:
    let e = exp(-s*ev.w[k])
    for j in 0..<d.n:
      let c = e*ev.v[j + d.n*k].re
      for i in 0..<d.n: result[i + d.n*j] += c*ev.v[i + d.n*k].re

proc expmv*(d: DenseM, ev: DenseEig, s: float, x: openArray[float]): seq[float] =
  ## exp(-M s) x, formed mode by mode so no dense matrix is ever squared.
  result = newSeq[float](d.n)
  for k in 0..<d.n:
    var c = 0.0
    for i in 0..<d.n: c += ev.v[i + d.n*k].re*x[i]
    c *= exp(-s*ev.w[k])
    for i in 0..<d.n: result[i] += c*ev.v[i + d.n*k].re
