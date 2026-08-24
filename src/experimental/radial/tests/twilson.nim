#RUNCMD $RUN1
## WP-E acceptance tests: doc/03-targets.md T1.2a - T1.2h.
## Also writes the Fig. 4 / Fig. 5 data of arXiv:2510.03085 to output/radial/wilson/.

import std/[algorithm, math, complex, os, strformat, unittest]
import base/alignedMem
import eigens/linalgFuncs
import ../core/analytic
import ../meas/dataio
import ../ops/wilson

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))

const
  outDir = currentSourcePath().parentDir.parentDir.parentDir.parentDir.parentDir /
           "output" / "radial" / "wilson"

# --- helpers ----------------------------------------------------------------

proc randGauge(l: Lat, sed: int, amp = 1.0): Gauge =
  result = newGauge(l)
  var r: Threefry4x64
  r.seedIndep(sed, 0)
  for i in 0..<result.s.len: result.s[i] = amp*r.gaussian
  for i in 0..<result.t.len: result.t[i] = amp*r.gaussian

proc randSpin(n, sed: int): Spin =
  result = newSpin(n)
  var r: Threefry4x64
  r.seedIndep(sed, 0)
  result.gaussian r

proc randReal(n, sed: int, amp = 1.0): seq[float] =
  result = newSeq[float](n)
  var r: Threefry4x64
  r.seedIndep(sed, 0)
  for i in 0..<n: result[i] = amp*r.gaussian

proc denseApply(a: seq[Complex64], nd: int, dst: var Spin, src: Spin, dag = false) =
  ## dst = A src (or A^dag src) for the column-major matrix `a`.
  for i in 0..<nd:
    var s = complex64(0.0, 0.0)
    for j in 0..<nd:
      let m = if dag: conjugate(a[j + nd*i]) else: a[i + nd*j]
      s += m*src[j shr 1][j and 1]
    dst[i shr 1][i and 1] = s

proc maxAbs(a: seq[Complex64]): float =
  for z in a: result = max(result, abs(z))

proc antiHermErr(a: seq[Complex64], nd: int): float =
  ## max |A + A^dag|
  for i in 0..<nd:
    for j in 0..<nd:
      result = max(result, abs(a[i + nd*j] + conjugate(a[j + nd*i])))

proc hermErr(a: seq[Complex64], nd: int): float =
  ## max |A - A^dag|
  for i in 0..<nd:
    for j in 0..<nd:
      result = max(result, abs(a[i + nd*j] - conjugate(a[j + nd*i])))

proc reldiff(x, y: Spin): float =
  ## |x - y| / |y|
  var d = 0.0
  for i in 0..<x.len:
    for c in 0..1:
      let z = x[i][c] - y[i][c]
      d += abs2(z)
  sqrt(d/norm2(y))

proc eigvals(a: var seq[Complex64], nd: int): seq[Complex64] =
  ## General complex eigenvalues; `a` is destroyed.
  result = newSeq[Complex64](nd)
  zgeigs(cast[ptr float64](addr a[0]), cast[ptr float64](addr result[0]), nd)

proc byIm(x, y: Complex64): int =
  if x.im < y.im: -1 elif x.im > y.im: 1
  elif x.re < y.re: -1 elif x.re > y.re: 1 else: 0

# --- fixtures ---------------------------------------------------------------

let
  sph1 = newSphere(1)
  sph2 = newSphere(2)
  lat1 = newLat(sph1, 6, 0.4)        ## small enough for a dense 144x144
  lat2 = newLat(sph2, 4, 0.3)
  u1 = randGauge(lat1, 20250101, 0.7)
  u2 = randGauge(lat2, 20250102, 0.7)

suite "T1.2a  naive C antihermitian, Wilson B hermitian":

  test "C + C^dag = 0 and B - B^dag = 0":
    for (nm, l, u) in [("L1 nt6", lat1, u1), ("L2 nt4", lat2, u2)]:
      let nd = 2*l.nsite
      let c = denseDw(l, u, 0.0, dwC)
      let b = denseDw(l, u, 0.0, dwB)
      let
        ec = antiHermErr(c, nd)/maxAbs(c)
        eb = hermErr(b, nd)/maxAbs(b)
      echo &"  {nm}: |C+C^dag|/|C| = {ec:.3e}   |B-B^dag|/|B| = {eb:.3e}"
      check ec < 1e-13
      check eb < 1e-13

  test "C + B = D_W exactly":
    for (nm, l, u) in [("L1 nt6", lat1, u1), ("L2 nt4", lat2, u2)]:
      let nd = 2*l.nsite
      let
        c = denseDw(l, u, 0.0, dwC)
        b = denseDw(l, u, 0.0, dwB)
        d = denseDw(l, u, 0.0, dwAll)
      var e = 0.0
      for i in 0..<nd*nd: e = max(e, abs(d[i] - b[i] - c[i]))
      echo &"  {nm}: max |D - B - C| = {e:.3e}"
      check e < 1e-14*maxAbs(d)

  test "C is pure imaginary spectrum, B is pure real":
    ## Consequence of the split quoted below (IV.4).
    let l = lat1
    let nd = 2*l.nsite
    var c = denseDw(l, u1, 0.0, dwC)
    var b = denseDw(l, u1, 0.0, dwB)
    let
      ec = eigvals(c, nd)
      eb = eigvals(b, nd)
    var mr = 0.0
    var mi = 0.0
    for z in ec: mr = max(mr, abs(z.re))
    for z in eb: mi = max(mi, abs(z.im))
    echo &"  max |Re eig C| = {mr:.3e}   max |Im eig B| = {mi:.3e}"
    check mr < 1e-12
    check mi < 1e-12

suite "T1.2b  explicit adjoint":

  test "<x, D y> = <D^dag x, y>":
    for (nm, l, u) in [("L1 nt6", lat1, u1), ("L2 nt4", lat2, u2)]:
      let
        x = randSpin(l.nsite, 991)
        y = randSpin(l.nsite, 992)
      var dy = newSpin(l.nsite)
      var dx = newSpin(l.nsite)
      for m in [0.0, 1.0]:
        applyDw(l, dy, y, u, m)
        applyDwAdj(l, dx, x, u, m)
        let
          a = dot(x, dy)
          b = dot(dx, y)
          sc = sqrt(norm2(x)*norm2(dy))
          e = abs(a - b)/sc
        echo &"  {nm} m={m}: |<x,Dy> - <D^dag x,y>|/|x||Dy| = {e:.3e}"
        check e < 1e-13

  test "applyDwAdj equals the conjugate transpose of denseDw":
    for (nm, l, u) in [("L1 nt6", lat1, u1), ("L2 nt4", lat2, u2)]:
      let nd = 2*l.nsite
      let a = denseDw(l, u, 0.0)
      let x = randSpin(l.nsite, 993)
      var p = newSpin(l.nsite)
      var q = newSpin(l.nsite)
      applyDwAdj(l, p, x, u)
      denseApply(a, nd, q, x, dag = true)
      let e = reldiff(p, q)
      echo &"  {nm}: |D^dag x - A^dag x|/|A^dag x| = {e:.3e}"
      check e < 1e-13

suite "T1.2c  U(1) gauge covariance":

  test "D(theta + d alpha) e^{i alpha} = e^{i alpha} D(theta)":
    for (nm, l, u) in [("L1 nt6", lat1, u1), ("L2 nt4", lat2, u2)]:
      # the seam and the temporal wrap must both be exercised
      var nsgn = 0
      for e in l.sph.edges:
        if e.sgn < 0.0: inc nsgn
      check nsgn > 0
      let alpha = randReal(l.nsite, 4711, 1.3)
      var ug = u
      gaugeTransform(l, ug, alpha)
      # the temporal seam link really did move
      var seam = 0.0
      for y in 0..<l.sph.nv:
        seam = max(seam, abs(ug.t[y + l.sph.nv*(l.nt-1)] - u.t[y + l.sph.nv*(l.nt-1)]))
      check seam > 1e-3
      let x = randSpin(l.nsite, 4712)
      var xg = x
      spinGaugeTransform(l, xg, alpha)
      var a = newSpin(l.nsite)
      var b = newSpin(l.nsite)
      applyDw(l, a, x, u, 0.3)            ## D(theta) x
      spinGaugeTransform(l, a, alpha)     ## e^{i alpha} D(theta) x
      applyDw(l, b, xg, ug, 0.3)          ## D(theta + d alpha) e^{i alpha} x
      let e = reldiff(a, b)
      echo &"  {nm}: |dD|/|D| = {e:.3e}  (squared: {e*e:.3e}), sgn<0 edges = {nsgn}"
      check e*e < 1e-20                   ## T1.2c read as a relative squared residual
      check e < 1e-14                     ## the honest amplitude-level statement

  test "a pure gauge on a Z2 seam edge is still covariant":
    ## Isolate the phi seam: alpha supported on one endpoint of a sgn = -1 edge.
    let l = lat2
    var esgn = -1
    for i, e in l.sph.edges:
      if e.sgn < 0.0:
        esgn = i
        break
    check esgn >= 0
    var alpha = newSeq[float](l.nsite)
    alpha[l.sph.edges[esgn].a] = 0.9
    alpha[l.sph.edges[esgn].b + l.sph.nv*(l.nt-1)] = -1.4
    var ug = u2
    gaugeTransform(l, ug, alpha)
    let x = randSpin(l.nsite, 4713)
    var xg = x
    spinGaugeTransform(l, xg, alpha)
    var a = newSpin(l.nsite)
    var b = newSpin(l.nsite)
    applyDw(l, a, x, u2)
    spinGaugeTransform(l, a, alpha)
    applyDw(l, b, xg, ug)
    let e = reldiff(a, b)
    echo &"  seam edge {esgn}: |dD|/|D| = {e:.3e}"
    check e < 1e-14

suite "T1.2d  matrix-free apply equals the dense assembly":

  test "applyDw == denseDw":
    for (nm, l, u) in [("L1 nt6", lat1, u1), ("L2 nt4", lat2, u2)]:
      let nd = 2*l.nsite
      for m in [0.0, 0.7]:
        for (pn, p) in [("all", dwAll), ("C", dwC), ("B", dwB),
                        ("space", dwSpace), ("time", dwTime)]:
          let a = denseDw(l, u, m, p)
          let x = randSpin(l.nsite, 5150)
          var q = newSpin(l.nsite)
          var r = newSpin(l.nsite)
          applyDw(l, q, x, u, m, p)
          denseApply(a, nd, r, x)
          let e = reldiff(q, r)
          if pn == "all":
            echo &"  {nm} m={m} {pn}: |dD x|/|D x| = {e:.3e}"
          check e < 1e-12

suite "T1.2e  antiperiodic temporal mode":

  test "half-integer Matsubara eigenvalue kappa'(1 - e^{-+ik})":
    ## With U = 0 and psi_{y,t} = e^{ikt} chi_pm, chi_pm the sigma_3 eigenspinors,
    ## the temporal operator of (IV.1) gives, for the antiperiodic extension
    ## f(t+nt) = -f(t), i.e. k = (2n+1) pi / nt,
    ##    D_T psi = kappa'_y [1 - cos k +- i sin k] psi = kappa'_y (1 - e^{-+ik}) psi.
    for (nm, l) in [("L1 nt6", lat1), ("L2 nt4", lat2)]:
      let
        u0 = newGauge(l)
        nv = l.sph.nv
      var worst = 0.0
      var x = newSpin(l.nsite)
      var y = newSpin(l.nsite)
      for n in 0..<l.nt:
        let k = PI*float(2*n + 1)/float(l.nt)
        for c in 0..1:
          let sg = if c == 0: 1.0 else: -1.0     ## sigma_3 eigenvalue
          x.zero
          for t in 0..<l.nt:
            let ph = complex64(cos(k*float(t)), sin(k*float(t)))
            for v in 0..<nv: x[v + nv*t][c] = ph
          applyDw(l, y, x, u0, 0.0, dwTime)
          let lam = complex64(1.0 - cos(k), sg*sin(k))
          for i in 0..<l.nsite:
            let ex = (l.kapT[i mod nv]*lam)*x[i][c]
            worst = max(worst, abs(y[i][c] - ex))
            worst = max(worst, abs(y[i][1 - c]))
      echo &"  {nm}: max |D_T psi - kappa'(1 - e^-+ik) psi| = {worst:.3e}"
      check worst < 1e-12

  test "integer Matsubara modes are NOT eigenvectors":
    ## The explicit -1 on the t = nt-1 -> 0 seam is what makes this fail; without it
    ## the test above would hold at k = 2 pi n / nt instead.
    let l = lat1
    let
      u0 = newGauge(l)
      nv = l.sph.nv
    var x = newSpin(l.nsite)
    var y = newSpin(l.nsite)
    let k = 2.0*PI/float(l.nt)
    x.zero
    for t in 0..<l.nt:
      let ph = complex64(cos(k*float(t)), sin(k*float(t)))
      for v in 0..<nv: x[v + nv*t][0] = ph
    applyDw(l, y, x, u0, 0.0, dwTime)
    let lam = complex64(1.0 - cos(k), sin(k))
    var dev = 0.0
    for i in 0..<l.nsite:
      dev = max(dev, abs(y[i][0] - (l.kapT[i mod nv]*lam)*x[i][0]))
    echo &"  integer-Matsubara deviation = {dev:.3e} (must be O(kappa'))"
    check dev > 0.1*l.kapT[0]

suite "derivatives":

  test "applyDwDeriv vs centered finite differences":
    for (nm, l, u) in [("L1 nt6", lat1, u1), ("L2 nt4", lat2, u2)]:
      let
        du = randGauge(l, 8801, 1.0)
        x = randSpin(l.nsite, 8802)
      var an = newSpin(l.nsite)
      applyDwDeriv(l, an, x, u, du)
      var up = newGauge(l)
      var um = newGauge(l)
      var a = newSpin(l.nsite)
      var b = newSpin(l.nsite)
      var best = 1e300
      var bestEps = 0.0
      for e in [1e-3, 1e-4, 1e-5, 3e-6, 1e-6]:
        for i in 0..<u.s.len:
          up.s[i] = u.s[i] + e*du.s[i]
          um.s[i] = u.s[i] - e*du.s[i]
        for i in 0..<u.t.len:
          up.t[i] = u.t[i] + e*du.t[i]
          um.t[i] = u.t[i] - e*du.t[i]
        applyDw(l, a, x, up)
        applyDw(l, b, x, um)
        for i in 0..<l.nsite:
          for c in 0..1:
            a[i][c] = (0.5/e)*(a[i][c] - b[i][c])
        let r = reldiff(a, an)
        if r < best:
          best = r
          bestEps = e
      echo &"  {nm}: best |dD_fd - dD|/|dD| = {best:.3e} at eps = {bestEps:.0e}"
      check best < 1e-10

  test "pullback vs tangent contraction":
    for (nm, l, u) in [("L1 nt6", lat1, u1), ("L2 nt4", lat2, u2)]:
      let
        du = randGauge(l, 8803, 1.0)
        left = randSpin(l.nsite, 8804)
        right = randSpin(l.nsite, 8805)
      var dd = newSpin(l.nsite)
      applyDwDeriv(l, dd, right, u, du)
      let rhs = 2.0*redot(left, dd)
      var f = newGauge(l)
      dwPullback(l, f, left, right, u)
      var lhs = 0.0
      for i in 0..<f.s.len: lhs += f.s[i]*du.s[i]
      for i in 0..<f.t.len: lhs += f.t[i]*du.t[i]
      let e = abs(lhs - rhs)/max(abs(rhs), 1e-300)
      echo &"  {nm}: <f,du> = {lhs:.12e}  2Re<x,dD y> = {rhs:.12e}  rel = {e:.3e}"
      check e < 1e-12

  test "pullback vs finite differences, link by link":
    let l = lat1
    let
      u = u1
      left = randSpin(l.nsite, 8806)
      right = randSpin(l.nsite, 8807)
    var f = newGauge(l)
    dwPullback(l, f, left, right, u)
    var up = u
    var a = newSpin(l.nsite)
    var b = newSpin(l.nsite)
    let eps = 1e-5
    var worst = 0.0
    var sc = 0.0
    for i in 0..<f.s.len: sc = max(sc, abs(f.s[i]))
    for i in 0..<f.t.len: sc = max(sc, abs(f.t[i]))
    for i in countup(0, f.s.len - 1, 7):
      up.s[i] = u.s[i] + eps
      applyDw(l, a, right, up)
      up.s[i] = u.s[i] - eps
      applyDw(l, b, right, up)
      up.s[i] = u.s[i]
      let fd = (2.0*redot(left, a) - 2.0*redot(left, b))/(2.0*eps)
      worst = max(worst, abs(fd - f.s[i]))
    for i in countup(0, f.t.len - 1, 5):
      up.t[i] = u.t[i] + eps
      applyDw(l, a, right, up)
      up.t[i] = u.t[i] - eps
      applyDw(l, b, right, up)
      up.t[i] = u.t[i]
      let fd = (2.0*redot(left, a) - 2.0*redot(left, b))/(2.0*eps)
      worst = max(worst, abs(fd - f.t[i]))
    echo &"  worst |f_fd - f| = {worst:.3e}  (scale {sc:.3e})"
    check worst < 1e-8*sc

  test "add and scale semantics of dwPullback":
    let l = lat1
    let
      left = randSpin(l.nsite, 8808)
      right = randSpin(l.nsite, 8809)
    var f1 = newGauge(l)
    var f2 = newGauge(l)
    dwPullback(l, f1, left, right, u1, 2.0)
    dwPullback(l, f2, left, right, u1, 0.5)
    dwPullback(l, f2, left, right, u1, 1.5, add = true)
    var e = 0.0
    for i in 0..<f1.s.len: e = max(e, abs(f1.s[i] - f2.s[i]))
    for i in 0..<f1.t.len: e = max(e, abs(f1.t[i] - f2.t[i]))
    echo &"  |scale 2 - (0.5 + 1.5)| = {e:.3e}"
    check e < 1e-12

suite "volume-normalized kernel":

  test "Dhat = S D S with S = abar/sqrt(A_y)":
    for (nm, l, u) in [("L1 nt6", lat1, u1), ("L2 nt4", lat2, u2)]:
      let
        sc = hatScale(l)
        nv = l.sph.nv
        x = randSpin(l.nsite, 6001)
      var w = newSpin(l.nsite)
      var a = newSpin(l.nsite)
      var b = newSpin(l.nsite)
      applyDwHat(l, a, x, u, w, 0.25)
      # reference: scale, apply, scale, subtract
      var y = newSpin(l.nsite)
      for i in 0..<l.nsite:
        let s = sc[i mod nv]
        y[i][0] = s*x[i][0]
        y[i][1] = s*x[i][1]
      applyDw(l, b, y, u)
      for i in 0..<l.nsite:
        let s = sc[i mod nv]
        b[i][0] = s*b[i][0] - 0.25*x[i][0]
        b[i][1] = s*b[i][1] - 0.25*x[i][1]
      let e = reldiff(a, b)
      echo &"  {nm}: |Dhat - S D S| / |.| = {e:.3e}"
      check e < 1e-14

  test "the hat adjoint is the adjoint":
    for (nm, l, u) in [("L1 nt6", lat1, u1), ("L2 nt4", lat2, u2)]:
      let
        x = randSpin(l.nsite, 6002)
        y = randSpin(l.nsite, 6003)
      var w = newSpin(l.nsite)
      var dy = newSpin(l.nsite)
      var dx = newSpin(l.nsite)
      applyDwHat(l, dy, y, u, w, 0.9)
      applyDwHatAdj(l, dx, x, u, w, 0.9)
      let
        a = dot(x, dy)
        b = dot(dx, y)
        e = abs(a - b)/sqrt(norm2(x)*norm2(dy))
      echo &"  {nm}: |<x,Dhat y> - <Dhat^dag x,y>| rel = {e:.3e}"
      check e < 1e-13

  test "gauge covariance of the hat operator":
    let l = lat2
    let alpha = randReal(l.nsite, 6004, 1.1)
    var ug = u2
    gaugeTransform(l, ug, alpha)
    let x = randSpin(l.nsite, 6005)
    var xg = x
    spinGaugeTransform(l, xg, alpha)
    var w = newSpin(l.nsite)
    var a = newSpin(l.nsite)
    var b = newSpin(l.nsite)
    applyDwHat(l, a, x, u2, w, 1.0)
    spinGaugeTransform(l, a, alpha)
    applyDwHat(l, b, xg, ug, w, 1.0)
    let e = reldiff(a, b)
    echo &"  |dDhat|/|Dhat| = {e:.3e}"
    check e < 1e-14

suite "performance discipline":

  test "applyDw allocates nothing":
    let l = lat2
    let u = u2
    var x = randSpin(l.nsite, 7001)
    var y = newSpin(l.nsite)
    var w = newSpin(l.nsite)
    var f = newGauge(l)
    let du = randGauge(l, 7002, 1.0)
    # warm up every code path once
    applyDw(l, y, x, u)
    applyDwAdj(l, y, x, u)
    applyDwDeriv(l, y, x, u, du)
    applyDwHat(l, y, x, u, w)
    applyDwHatAdj(l, y, x, u, w)
    dwPullback(l, f, x, y, u)
    let raw0 = getRawMemAllocated()
    let occ0 = getOccupiedMem()
    var acc = 0.0
    for k in 0..<64:
      applyDw(l, y, x, u, 0.3)
      applyDwAdj(l, x, y, u, 0.3)
      applyDwDeriv(l, y, x, u, du)
      applyDwHat(l, x, y, u, w, 0.1)
      applyDwHatAdj(l, y, x, u, w, 0.1)
      dwPullback(l, f, x, y, u, 1e-9, add = true)
      let s = 1.0/sqrt(norm2(x) + norm2(y))
      scale(x, s)
      scale(y, s)
      acc += s
    let raw1 = getRawMemAllocated()
    let occ1 = getOccupiedMem()
    echo &"  raw {raw0} -> {raw1}   occupied {occ0} -> {occ1}   acc = {acc:.6e}"
    check raw1 == raw0
    check occ1 == occ0
    check acc != 0.0

# --- free-field spectra -----------------------------------------------------
#
# With U = 0 the operator is invariant under translations in t, so it block
# diagonalizes on the antiperiodic Matsubara modes e^{i k t}, k = (2n+1) pi / nt:
#     D(k) = D_spatial + kappa'_y (1 - cos k + i sigma_3 sin k).
# That turns the 2*nv*nt dense problem into nt problems of size 2*nv, which is what
# makes the L = 4, Lt = 24 panel of Fig. 4 (a 7776 x 7776 matrix) cheap.  The
# decomposition is checked against a full dense diagonalization below.

proc spatialDense(sph: Sphere, at: float): tuple[l: Lat, a: seq[Complex64]] =
  ## The U = 0 spatial operator, 2*nv x 2*nv.  `nt = 1` because the spatial part does
  ## not see the time direction; the returned Lat still carries kappa'(at).
  let l = newLat(sph, 1, at)
  (l, denseDw(l, newGauge(l), 0.0, dwSpace))

proc matsuBlock(l: Lat, dsp: seq[Complex64], k: float): seq[Complex64] =
  result = dsp
  let
    nv = l.sph.nv
    nd = 2*nv
  for y in 0..<nv:
    let
      kt = l.kapT[y]
      dr = kt*(1.0 - cos(k))
      di = kt*sin(k)
    result[(2*y) + nd*(2*y)] += complex64(dr, di)
    result[(2*y + 1) + nd*(2*y + 1)] += complex64(dr, -di)

proc freeSpectrum(sph: Sphere, nt: int, at: float): seq[Complex64] =
  ## Full spectrum of the free D_W on sph x nt slices, Matsubara mode by Matsubara mode.
  let (l, dsp) = spatialDense(sph, at)
  result = newSeqOfCap[Complex64](2*sph.nv*nt)
  for n in 0..<nt:
    var b = matsuBlock(l, dsp, PI*float(2*n + 1)/float(nt))
    for z in eigvals(b, 2*sph.nv): result.add z

proc hatDense(l: Lat, a: seq[Complex64]): seq[Complex64] =
  ## S A S with S = diag(abar/sqrt(A_y)): the overlap normalization Dhat_W.
  result = a
  let
    nv = l.sph.nv
    nd = 2*l.nsite
    sc = hatScale(l)
  for i in 0..<nd:
    for j in 0..<nd:
      result[i + nd*j] = (sc[(i div 2) mod nv]*sc[(j div 2) mod nv])*result[i + nd*j]

proc volDense(l: Lat, a: seq[Complex64]): seq[Complex64] =
  ## diag(volw/volbar) A -- the literal reading of doc/02 section 3.2.
  result = a
  let
    nv = l.sph.nv
    nd = 2*l.nsite
  for i in 0..<nd:
    let w = l.volw[(i div 2) mod nv]/l.volbar
    for j in 0..<nd:
      result[i + nd*j] = w*result[i + nd*j]

proc byAbs(x, y: Complex64): int =
  if abs(x) < abs(y): -1 elif abs(x) > abs(y): 1 else: byIm(x, y)

proc doubler(ev: seq[Complex64], cut: float): tuple[re, im: float] =
  ## Smallest Re over the eigenvalues that sit within `cut` of the positive real axis:
  ## the lattice image of the "lambda real and nonzero" doubler point of (IV.8).
  result = (1e300, 0.0)
  for z in ev:
    if abs(z.im) <= cut and z.re > cut and z.re < result.re:
      result = (z.re, abs(z.im))

const towerMult = [4, 8, 12]     ## 2 (l+1) eigenvalues per sign, l = 0, 1, 2

suite "T1.2h  free spatial spectrum":

  test "volume-rescaled spectrum approaches +- i (l+1)":
    ## Dhat_W = abar^2 A^{-1/2} D_W A^{-1/2} has eig = abar * eig(D_continuum), so the
    ## continuum check is on eig(Dhat_spatial)/abar.  The residual is dominated by the
    ## Wilson term's O(abar) real part; the imaginary part converges as O(abar^2).
    var prev = newSeq[float](towerMult.len)
    var prevIm = newSeq[float](towerMult.len)
    for i in 0..<towerMult.len:
      prev[i] = 1e300
      prevIm[i] = 1e300
    for lev in [1, 2, 4]:
      let sph = newSphere(lev)
      let (l, a) = spatialDense(sph, 0.2)
      var h = hatDense(l, a)
      var ev = eigvals(h, 2*sph.nv)
      for i in 0..<ev.len: ev[i] = ev[i]/sph.abar
      ev.sort(byAbs)
      var k = 0
      for li in 0..<towerMult.len:
        let
          mult = towerMult[li]
          target = float(li + 1)
        var
          err = 0.0
          eim = 0.0
          npos = 0
        for j in k..<(k + mult):
          let z = ev[j]
          if z.im > 0.0: inc npos
          let w = complex64(z.re, z.im - (if z.im > 0.0: target else: -target))
          err = max(err, abs(w))
          eim = max(eim, abs(abs(z.im) - target))
        # the group must be separated from the next one (L=1 has no modes past l=2)
        let gap = if k + mult < ev.len: abs(ev[k + mult]) - abs(ev[k + mult - 1])
                  else: NaN
        echo &"  L={lev} l={li}: mult {mult} ({npos} with Im>0)  " &
             &"max|lam -+ i(l+1)| = {err:.5f}  max||Im|-(l+1)| = {eim:.5f}  gap = {gap:.4f}"
        check npos == mult div 2                  ## conjugate-symmetric multiplet
        if k + mult < ev.len: check gap > 0.02
        check err < prev[li]                      ## the error decreases with L
        check eim < prevIm[li]
        prev[li] = err
        prevIm[li] = eim
        k += mult

  test "the doc/02 section 3.2 weight diag(volw/volbar) does not converge":
    ## (IV.11) says abar a_t D_lat = deltaV D_cont, so D_cont = abar diag(1/A_y) D_lat
    ## and the weight is deltaVbar/deltaV, not deltaV/deltaVbar.  Recorded because
    ## doc/02 section 3.2 states the reciprocal; the numbers below settle it.
    for lev in [1, 2, 4]:
      let sph = newSphere(lev)
      let (l, a) = spatialDense(sph, 0.2)
      var h = hatDense(l, a)
      var v = volDense(l, a)
      var eh = eigvals(h, 2*sph.nv)
      var evv = eigvals(v, 2*sph.nv)
      for i in 0..<eh.len: eh[i] = eh[i]/sph.abar
      eh.sort(byAbs)
      evv.sort(byAbs)
      var ih = 0.0
      var iv = 0.0
      for j in 0..<4:
        ih = max(ih, abs(abs(eh[j].im) - 1.0))
        iv = max(iv, abs(abs(evv[j].im) - 1.0))
      echo &"  L={lev}: l=0 ||Im|-1|  hat/abar = {ih:.5f}   volw/volbar = {iv:.5f}"
      check ih < iv

suite "T1.2f, T1.2g  Fig. 4 and Fig. 5":

  test "the Matsubara decomposition reproduces the full dense spectrum":
    let
      sph = sph1
      nt = 6
      at = 0.4
      l = newLat(sph, nt, at)
    var d = denseDw(l, newGauge(l), 0.0)
    var full = eigvals(d, 2*l.nsite)
    var blk = freeSpectrum(sph, nt, at)
    check full.len == blk.len
    # greedy nearest-neighbour matching: the eigenvalues are 8- and 12-fold degenerate
    # to the last bit, so a lexicographic sort does not align the two lists.
    var used = newSeq[bool](full.len)
    var e = 0.0
    for z in blk:
      var best = 1e300
      var bi = -1
      for j in 0..<full.len:
        if used[j]: continue
        let d = abs(z - full[j])
        if d < best:
          best = d
          bi = j
      used[bi] = true
      e = max(e, best)
    echo &"  L1 nt6: max |lam_dense - lam_matsubara| = {e:.3e} over {full.len} modes"
    check e < 1e-10

  test "Fig. 4 and Fig. 5 data, and the doubler comparison":
    createDir(outDir)
    const
      tt = 4.0
      cases = [(1, 24), (2, 24), (4, 24), (2, 16), (2, 48)]
      kflat = 1.0/sqrt(3.0)
    var
      cl, clt, cabar, cat, ckt, cdb, cdim, cflat, cmaxre, cfmaxre: seq[float]
    for (lev, nt) in cases:
      let
        sph = newSphere(lev)
        at = tt/float(nt)
        ktflat = 0.5*sqrt(3.0)*sph.abar/at
        ev = freeSpectrum(sph, nt, at)
      # the raw spectrum, in the units of (IV.8)
      var re = newSeq[float](ev.len)
      var im = newSeq[float](ev.len)
      var idx = newSeq[float](ev.len)
      for i, z in ev:
        re[i] = z.re
        im[i] = z.im
        idx[i] = float(i div (2*sph.nv))
      writeTsv(outDir/(&"fig4_L{lev}_Lt{nt}.tsv"),
               {"source": "arXiv:2510.03085 Fig. 4", "lattice": &"L{lev}",
                "Lt": $nt, "T": $tt, "at": &"{at:.17g}",
                "abar": &"{sph.abar:.17g}", "normalization": "raw D_lat of (IV.1)"},
               ["re", "im", "matsubara"], [re, im, idx])
      # matching flat spectrum, (IV.8)
      let fl = flatSpectrum(kflat, ktflat, 16, 16, nt)
      var fre = newSeq[float](fl.len)
      var fim = newSeq[float](fl.len)
      for i, z in fl:
        fre[i] = z.re
        fim[i] = z.im
      writeTsv(outDir/(&"fig5_flat_L{lev}_Lt{nt}.tsv"),
               {"source": "arXiv:2510.03085 Fig. 5, Eq. (IV.8)",
                "kappa": &"{kflat:.17g}", "kappaT": &"{ktflat:.17g}",
                "grid": "16x16x" & $nt, "abar": &"{sph.abar:.17g}", "at": &"{at:.17g}"},
               ["re", "im"], [fre, fim])
      # the doubler, from the spatial (k_t -> 0) operator, where 4 kappa lives
      let (_, sa) = spatialDense(sph, at)
      var sv = sa
      let sev = eigvals(sv, 2*sph.nv)
      let cutd = 0.4*sph.abar
      let (dre, dim) = doubler(sev, cutd)
      var mre = 0.0
      for z in ev: mre = max(mre, z.re)
      var fmre = 0.0
      for z in fl: fmre = max(fmre, z.re)
      let dflat = min(4.0*kflat, 2.0*ktflat)
      # the temporal doubler is exact: the k_t = pi mode of the temporal operator has
      # eigenvalue 2 kappa'_y, site by site (see the T1.2e suite).
      var ktlo = 1e300
      var kthi = 0.0
      var ktbar = 0.0
      for y in 0..<sph.nv:
        let kt = sph.area[y]/(sph.abar*at)
        ktlo = min(ktlo, kt)
        kthi = max(kthi, kt)
        ktbar += kt
      ktbar /= float(sph.nv)
      echo &"  L={lev} Lt={nt} at={at:.5f} abar={sph.abar:.5f} kappa'flat={ktflat:.4f}"
      echo &"    spatial doubler = {dre:.6f} (|Im| = {dim:.4f}, cut {cutd:.3f})" &
           &"   flat 4 kappa = {4.0*kflat:.6f}   rel = {abs(dre/(4.0*kflat) - 1.0):.3e}"
      echo &"    temporal doubler 2 kappa'_y in [{2.0*ktlo:.5f}, {2.0*kthi:.5f}]" &
           &" mean {2.0*ktbar:.5f}   flat 2 kappa' = {2.0*ktflat:.5f}" &
           &"   rel(mean) = {abs(ktbar/ktflat - 1.0):.3e}"
      echo &"    max Re: curved {mre:.5f}   flat grid {fmre:.5f}" &
           &"   rel = {abs(mre/fmre - 1.0):.3e}" &
           &"   (exact 4.5k+2k' = {4.5*kflat + 2.0*ktflat:.5f}; the 16x16 grid misses the" &
           &" K point where 4.5 kappa lives)"
      # 0.025: the exact-kappa convention (doc/06 "THE COUPLING CONVENTION") shifts
      # the L=1 spatial doubler to 2.13% from the flat estimate (was 1.51% flat-kappa);
      # (IV.10)'s alpha = 0.9 safety factor absorbs far more than this.
      check abs(dre/(4.0*kflat) - 1.0) < 0.025
      check abs(ktbar/ktflat - 1.0) < 0.02
      check abs(mre/fmre - 1.0) < 0.06
      cl.add float(lev)
      clt.add float(nt)
      cabar.add sph.abar
      cat.add at
      ckt.add ktflat
      cdb.add dre
      cdim.add dim
      cflat.add dflat
      cmaxre.add mre
      cfmaxre.add fmre
    writeTsv(outDir/"doublers.tsv",
             {"source": "WP-E, doc/03-targets T1.2f/T1.2g",
              "doubler_def": "min Re over eig(D_spatial) with |Im| <= 0.4 abar",
              "flat_doubler": "min(4 kappa, 2 kappa') of (IV.8), kappa = 1/sqrt(3)"},
             ["L", "Lt", "abar", "at", "kappaTflat", "doubler", "doublerIm",
              "flatDoubler", "maxReCurved", "maxReFlat"],
             [cl, clt, cabar, cat, ckt, cdb, cdim, cflat, cmaxre, cfmaxre])
    check fileExists(outDir/"doublers.tsv")
    echo "  wrote ", outDir

suite "normalization cross-check":

  test "applyDwHat equals the dense S D S used for the spectra":
    ## Closes the loop between the matrix-free kernel WP-F will call and the dense
    ## matrix the T1.2h spectra are read off.
    for (nm, l, u) in [("L1 nt6", lat1, u1), ("L2 nt4", lat2, u2)]:
      let
        nd = 2*l.nsite
        a = hatDense(l, denseDw(l, u, 0.0))
        x = randSpin(l.nsite, 6100)
      var w = newSpin(l.nsite)
      var p = newSpin(l.nsite)
      var q = newSpin(l.nsite)
      applyDwHat(l, p, x, u, w)
      denseApply(a, nd, q, x)
      let e = reldiff(p, q)
      echo &"  {nm}: |Dhat x - (S A S) x| / |.| = {e:.3e}"
      check e < 1e-13

  test "Ward identity: delta_{d alpha} D = i (alpha D - D alpha)":
    ## The infinitesimal form of T1.2c, and an independent check of applyDwDeriv:
    ## D(theta + eps d alpha) = e^{i eps alpha} D(theta) e^{-i eps alpha} differentiates
    ## to delta D = i [alpha, D] = i (alpha D - D alpha).  NOTE the sign: doc/04 section
    ## 15 item 6 writes i (D alpha - alpha D), which is the opposite link convention
    ## (theta_e -> theta_e + alpha_a - alpha_b).  Ours is doc/02 section 5's,
    ## theta_e -> theta_e + alpha_b - alpha_a, and it is the one that makes the
    ## plaquettes of (IV.24) invariant.
    for (nm, l, u) in [("L1 nt6", lat1, u1), ("L2 nt4", lat2, u2)]:
      let
        alpha = randReal(l.nsite, 3141, 1.0)
        x = randSpin(l.nsite, 2718)
      var da = newGauge(l)
      gaugeTransform(l, da, alpha)          ## da = d alpha, since da started at zero
      var p = newSpin(l.nsite)
      applyDwDeriv(l, p, x, u, da)
      var ax = x
      for i in 0..<l.nsite:
        ax[i][0] = alpha[i]*x[i][0]
        ax[i][1] = alpha[i]*x[i][1]
      var q = newSpin(l.nsite)
      var s = newSpin(l.nsite)
      applyDw(l, q, ax, u)                  ## D alpha x
      applyDw(l, s, x, u)                   ## D x
      for i in 0..<l.nsite:
        for c in 0..1:
          let z = alpha[i]*s[i][c] - q[i][c]
          q[i][c] = complex64(-z.im, z.re)  ## i (alpha D - D alpha) x
      let e = reldiff(p, q)
      echo &"  {nm}: |delta D x - i[alpha,D] x| / |.| = {e:.3e}"
      check e < 1e-13

  test "min |Dhat_W - 1| in the free limit (the T2.1 free column)":
    ## Not a WP-E acceptance criterion -- reported because the published free-limit
    ## values 1.154 / 1.010 / 0.965 (doc/03 T2.1) are a direct, deterministic test of
    ## which normalization goes into X = D_W - M.  a_t = 0.2, M = 1.
    const at = 0.2
    for lev in [1, 2, 4]:
      let sph = newSphere(lev)
      var line = &"  L={lev} abar={sph.abar:.5f} abar/at={sph.abar/at:.3f}:"
      for nt in [60, 80]:
        let (l, dsp) = spatialDense(sph, at)
        var dhat = 1e300
        var draw = 1e300
        let sc = hatScale(l)
        for n in 0..<nt:
          var b = matsuBlock(l, dsp, PI*float(2*n + 1)/float(nt))
          var h = b
          let nd = 2*sph.nv
          for i in 0..<nd:
            for j in 0..<nd:
              h[i + nd*j] = (sc[i div 2]*sc[j div 2])*h[i + nd*j]
          for z in eigvals(b, nd): draw = min(draw, abs(z - 1.0))
          for z in eigvals(h, nd): dhat = min(dhat, abs(z - 1.0))
        line.add &"  Lt={nt}: hat {dhat:.4f} raw {draw:.4f}"
      echo line
    # cheap L = 1 scan in Lt, to see whether the published 1.154 is an L_t effect
    block:
      let sph = newSphere(1)
      let (l, dsp) = spatialDense(sph, at)
      var line = "  L=1 raw min|D_W - 1| vs Lt:"
      for nt in [12, 16, 20, 24, 32, 48, 60, 80, 120, 168]:
        var d = 1e300
        for n in 0..<nt:
          var b = matsuBlock(l, dsp, PI*float(2*n + 1)/float(nt))
          for z in eigvals(b, 2*sph.nv): d = min(d, abs(z - 1.0))
        line.add &" {nt}:{d:.4f}"
      echo line
