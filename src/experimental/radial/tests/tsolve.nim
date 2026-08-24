#RUNCMD $RUN1
import std/[unittest, math, complex, strformat]
import eigens/linalgFuncs
import ../core/spinor
import ../ops/solve

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))

# A dense Hermitian positive-definite test operator with a prescribed spectrum.
# `a = h diag(w) h^dag`, with `h` the eigenvectors of a random Hermitian matrix
# (unitary to roundoff).  `(v, ev)` is a *fresh* eigendecomposition of the
# assembled `a` and supplies the reference inverse, so the oracle inverts the
# matrix the solver actually sees.  LAPACK column-major: m[i + nd*j] = M_ij.

type Dense = object
  nd: int
  a: seq[Complex64]         ## the operator
  v: seq[Complex64]         ## its eigenvectors, columns
  ev: seq[float]            ## its eigenvalues

proc newDense(nd: int, w: openArray[float], sed: int): Dense =
  var r: Threefry4x64
  r.seedIndep(sed, 0)
  var h = newSeq[Complex64](nd*nd)
  for j in 0..<nd:
    h[j + nd*j] = complex64(r.gaussian, 0.0)
    for i in (j+1)..<nd:
      let z = complex64(r.gaussian, r.gaussian)
      h[i + nd*j] = z
      h[j + nd*i] = conjugate(z)
  var ew = newSeq[float](nd)
  zeigs(cast[ptr float64](addr h[0]), addr ew[0], nd)   # h <- eigenvectors

  result.nd = nd
  result.a = newSeq[Complex64](nd*nd)
  for i in 0..<nd:
    for k in 0..i:
      var s = complex64(0.0, 0.0)
      for j in 0..<nd: s += w[j]*h[i + nd*j]*conjugate(h[k + nd*j])
      result.a[i + nd*k] = s
      result.a[k + nd*i] = conjugate(s)   # exactly Hermitian by construction

  result.v = result.a
  result.ev = newSeq[float](nd)
  zeigs(cast[ptr float64](addr result.v[0]), addr result.ev[0], nd)

proc apply(d: Dense, dst: var Spin, src: Spin) =
  let nd = d.nd
  for k in 0..<nd:
    var sr = 0.0
    var si = 0.0
    for j in 0..<nd:
      let m = d.a[k + nd*j]
      let s = src[j shr 1][j and 1]
      sr += m.re*s.re - m.im*s.im
      si += m.re*s.im + m.im*s.re
    dst[k shr 1][k and 1] = complex64(sr, si)

proc op(d: Dense): proc(dst: var Spin, src: Spin) =
  result = proc(dst: var Spin, src: Spin) = d.apply(dst, src)

proc exact(d: Dense, b: Spin, s: float): Spin =
  ## (A + s)^{-1} b from the eigendecomposition of the assembled A.
  let nd = d.nd
  var c = newSeq[Complex64](nd)          # c_j = <v_j, b> / (ev_j + s)
  for j in 0..<nd:
    var zr = 0.0
    var zi = 0.0
    for i in 0..<nd:
      let m = d.v[i + nd*j]
      let x = b[i shr 1][i and 1]
      zr += m.re*x.re + m.im*x.im
      zi += m.re*x.im - m.im*x.re
    let f = 1.0/(d.ev[j] + s)
    c[j] = complex64(f*zr, f*zi)
  result = newSpin(nd div 2)
  for i in 0..<nd:
    var zr = 0.0
    var zi = 0.0
    for j in 0..<nd:
      let m = d.v[i + nd*j]
      zr += m.re*c[j].re - m.im*c[j].im
      zi += m.re*c[j].im + m.im*c[j].re
    result[i shr 1][i and 1] = complex64(zr, zi)

proc logSpace(lo, hi: float, n: int): seq[float] =
  result = newSeq[float](n)
  for i in 0..<n:
    result[i] = lo*pow(hi/lo, float(i)/float(n-1))

proc reldiff2(x, y: Spin): float =
  var d = x
  axpy(d, -1.0, y)
  norm2(d)/norm2(y)

# ---------------------------------------------------------------------------
# easy operator: dimension 120, spectrum in [0.05, 5], cond 100
const nsEasy = 60
const ndEasy = 2*nsEasy
let dEasy = newDense(ndEasy, logSpace(0.05, 5.0, ndEasy), 20240921)
let opEasy = op(dEasy)

var rb: Threefry4x64
rb.seedIndep(4242, 0)
var bEasy = newSpin(nsEasy)
bEasy.gaussian rb

# hard operator: dimension 400, spectrum in [1e-8, 1], cond 1e8.  With shifts
# spanning 1e6 the base system A + 1e-6 still has cond ~1e6.
const nsHard = 200
const ndHard = 2*nsHard
let dHard = newDense(ndHard, logSpace(1e-8, 1.0, ndHard), 777)
let opHard = op(dHard)
let sHard = logSpace(1e-6, 1.0, 8)

var rh: Threefry4x64
rh.seedIndep(99, 0)
var bHard = newSpin(nsHard)
bHard.gaussian rh

suite "solve":

  test "cgSolve vs dense inverse":
    let r2req = 1e-26
    var x = newSpin(nsEasy)
    let info = cgSolve(x, bEasy, r2req, 5000, opEasy)
    let xe = dEasy.exact(bEasy, 0.0)
    let e2 = reldiff2(x, xe)
    echo &"  cgSolve: iters={info.iters} r2={info.r2:.3e} |dx|/|x|={sqrt(e2):.3e}"
    check info.converged
    check info.r2 <= 1.001*r2req
    check sqrt(e2) < 1e-9

  test "T1.3h multishift vs independent single-shift":
    let shifts = logSpace(1e-3, 3.0, 15)
    let r2req = 1e-26
    var xs: seq[Spin]
    let mi = cgmSolve(xs, bEasy, shifts, r2req, 5000, opEasy)
    echo &"  cgmSolve: iters={mi.iters} refined={mi.refined} refits={mi.refits}"
    check mi.converged
    var worst = 0.0
    var worstRes = 0.0
    for j in 0..<shifts.len:
      var x = newSpin(nsEasy)
      let ci = cgSolve(x, bEasy, r2req, 5000,
                       proc(dst: var Spin, src: Spin) =
                         dEasy.apply(dst, src)
                         axpy(dst, shifts[j], src))
      check ci.converged
      let e2 = reldiff2(xs[j], x)
      let xe = dEasy.exact(bEasy, shifts[j])
      let a2 = reldiff2(xs[j], xe)
      echo &"    s={shifts[j]:.4e} its={ci.iters:4d} r2m={mi.r2[j]:.3e} r2s={ci.r2:.3e}" &
           &" |dx|^2/|x|^2={e2:.3e} vs-dense={sqrt(a2):.3e}"
      if e2 > worst: worst = e2
      if mi.r2[j] > worstRes: worstRes = mi.r2[j]
      check mi.r2[j] <= 1.001*r2req
    echo &"  worst multishift-vs-single |dx|^2/|x|^2 = {worst:.3e}"
    check worst < 1e-18

  test "iteration count tracks the smallest shift":
    let shifts = logSpace(1e-3, 3.0, 15)
    let r2req = 1e-26
    var xs: seq[Spin]
    let mi = cgmSolve(xs, bEasy, shifts, r2req, 5000, opEasy)
    var x = newSpin(nsEasy)
    let ci = cgSolve(x, bEasy, r2req, 5000,
                     proc(dst: var Spin, src: Spin) =
                       dEasy.apply(dst, src)
                       axpy(dst, shifts[0], src))
    echo &"  cgm iters={mi.iters}  smallest-shift cg iters={ci.iters}"
    check mi.refined == 0
    check mi.iters <= ci.iters + 2

  test "determinism":
    let shifts = logSpace(1e-3, 3.0, 15)
    # 1e-24 converges cleanly; 1e-32 is below the roundoff floor and drives the
    # refinement branch, so both paths are covered.
    for r2req in [1e-24, 1e-32]:
      var xs1, xs2: seq[Spin]
      let i1 = cgmSolve(xs1, bEasy, shifts, r2req, 5000, opEasy)
      let i2 = cgmSolve(xs2, bEasy, shifts, r2req, 5000, opEasy)
      echo &"  r2req={r2req:.0e} iters={i1.iters} refined={i1.refined} refits={i1.refits}"
      check i1.iters == i2.iters
      check i1.refined == i2.refined
      check i1.refits == i2.refits
      for j in 0..<shifts.len:
        check i1.r2[j] == i2.r2[j]
        for i in 0..<nsEasy:
          for c in 0..1:
            check xs1[j][i][c].re == xs2[j][i][c].re
            check xs1[j][i][c].im == xs2[j][i][c].im

  test "hard case: cond 1e8, shift ratio 1e6":
    ## The recurrence estimate must hold up when the base system is genuinely
    ## hard.  r2req = 1e-16 is about two decades above the base's roundoff
    ## floor (measured: it stops meeting 1.001*r2req somewhere near 1e-18).
    let r2req = 1e-16
    var xs: seq[Spin]
    let mi = cgmSolve(xs, bHard, sHard, r2req, 20000, opHard)
    echo &"  iters={mi.iters} refined={mi.refined} refits={mi.refits} conv={mi.converged}"
    var mxPre = 0.0
    for j in 0..<sHard.len:
      let e2 = reldiff2(xs[j], dHard.exact(bHard, sHard[j]))
      echo &"    s={sHard[j]:.3e} r2/req={mi.r2[j]/r2req:.4f} vs-dense={sqrt(e2):.3e}"
      if mi.r2pre[j]/r2req > mxPre: mxPre = mi.r2pre[j]/r2req
      check mi.r2[j] <= 1.001*r2req
      check sqrt(e2) < 1e-7
    echo &"  worst recurrence residual before refinement = {mxPre:.4f} x r2req"
    check mi.converged
    # Pins the shift-freezing fix: without it the auxiliary z recurrence
    # underflows to 0/0 and NaN-poisons the far shifts, which then all need
    # refinement.  With it the estimate is good to better than 0.1%.
    check mi.refined == 0

  test "refinement fires when the request is below the roundoff floor":
    ## The one regime where the recurrence estimate really does lie: r2req far
    ## under the operator's attainable residual.  The fallback fires, improves
    ## the worst shifts, and cgmSolve reports honest failure instead of
    ## accepting the estimate.
    let r2req = 1e-24
    var xs: seq[Spin]
    let mi = cgmSolve(xs, bHard, sHard, r2req, 20000, opHard)
    echo &"  iters={mi.iters} refined={mi.refined} refits={mi.refits} conv={mi.converged}"
    var mxPre = 0.0
    var mxPost = 0.0
    for j in 1..<sHard.len:
      echo &"    s={sHard[j]:.3e} pre/req={mi.r2pre[j]/r2req:.4e} post/req={mi.r2[j]/r2req:.4e}"
      if mi.r2pre[j] > mxPre: mxPre = mi.r2pre[j]
      if mi.r2[j] > mxPost: mxPost = mi.r2[j]
    check mi.refined > 0            ## the fallback path runs
    check mxPre > 1.001*r2req       ## the recurrence really did miss
    check mxPost < mxPre            ## and the re-solve improved it
    check not mi.converged          ## but the floor is the floor: report it
