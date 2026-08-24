#RUNCMD $RUN1
## WP-F acceptance tests: doc/03-targets.md T1.3c-T1.3f, T1.3h (operator level)
## plus the doc/04-interfaces.md section 10 contract: adjoint, gauge covariance,
## ovGradient (tangent / pullback / Ward), kernelWindow, allocation, solve counts.
##
## Small L = 1 lattices so the dense polar factor is a cheap exact oracle.  Every
## operator test runs on the free field AND a random non-compact gauge field; each
## field gets its own frozen rational window, read off the dense spec(X^dag X).

import std/[math, complex, strformat, unittest]
import base/alignedMem
import eigens/linalgFuncs
import ../ops/overlap

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))

# --- helpers -----------------------------------------------------------------

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

proc denseApply(a: seq[Complex64], nd: int, dst: var Spin, src: Spin) =
  ## dst = A src for the column-major matrix `a`.
  for i in 0..<nd:
    var s = complex64(0.0, 0.0)
    for j in 0..<nd: s += a[i + nd*j]*src[j shr 1][j and 1]
    dst[i shr 1][i and 1] = s

proc reldiff(x, y: Spin): float =
  ## |x - y| / |y|
  var d = 0.0
  for i in 0..<x.len:
    for c in 0..1: d += abs2(x[i][c] - y[i][c])
  sqrt(d/norm2(y))

proc eigvals(a: seq[Complex64], nd: int): seq[Complex64] =
  var m = a
  result = newSeq[Complex64](nd)
  zgeigs(cast[ptr float64](addr m[0]), cast[ptr float64](addr result[0]), nd)

proc denseHBounds(x: seq[Complex64], nd: int): tuple[smin, smax: float] =
  ## sigma bounds of X from the eigenvalues of the dense X^dag X.
  var h = newSeq[Complex64](nd*nd)
  for j in 0..<nd:
    for i in j..<nd:
      var s = complex64(0.0, 0.0)
      for k in 0..<nd: s += conjugate(x[k + nd*i])*x[k + nd*j]
      h[i + nd*j] = s
      h[j + nd*i] = conjugate(s)
  var ev = newSeq[float](nd)
  zeigs(cast[ptr float64](addr h[0]), addr ev[0], nd)   # ascending
  (sqrt(ev[0]), sqrt(ev[nd-1]))

proc hShiftOp(l: Lat, u: Gauge, m, s: float): proc(dst: var Spin, src: Spin) =
  ## (X^dag X + s), built directly on applyDw/applyDwAdj -- independent of Ov.
  var w = newSpin(l.nsite)
  result = proc(dst: var Spin, src: Spin) =
    applyDw(l, w, src, u, m)
    applyDwAdj(l, dst, w, u, m)
    if s != 0.0: axpy(dst, s, src)

proc assembleOv(o: Ov, u: Gauge): seq[Complex64] =
  ## The rational D_ov, column by column through applyOv.
  let
    n = o.l.nsite
    nd = 2*n
  result = newSeq[Complex64](nd*nd)
  var b = newSpin(n)
  var col = newSpin(n)
  for j in 0..<nd:
    b.zero
    b[j shr 1][j and 1] = complex64(1.0, 0.0)
    applyOv(o, col, b, u)
    for i in 0..<nd: result[i + nd*j] = col[i shr 1][i and 1]

proc conjClosure(ev: seq[Complex64]): float =
  ## max over i of the distance from conj(ev_i) to the spectrum: 0 iff the
  ## spectrum is closed under complex conjugation (as a set).
  for z in ev:
    let zc = conjugate(z)
    var best = abs(zc - ev[0])
    for w in ev: best = min(best, abs(zc - w))
    result = max(result, best)

# --- fixtures ----------------------------------------------------------------
# L = 1, nt = 6, at = 0.4: nsite 72, dense dimension 144.  abar/at = 2.77 >= 4/3
# and M = 1 < maxM = 2.078.  Windows are the dense sigma bounds padded by 5%.

type Fix = object
  name: string
  u: Gauge
  smin, smax: float          ## dense sigma bounds of X
  x: seq[Complex64]          ## dense X
  dov: seq[Complex64]        ## dense exact D_ov (the oracle)
  o31, o11: Ov

const
  r2in = 1e-26
  r2out = 1e-22
  mxit = 20000

let
  sph = newSphere(1)
  lat = newLat(sph, 6, 0.4)
  nd = 2*lat.nsite

proc makeFix(name: string, u: Gauge): Fix =
  result.name = name
  result.u = u
  result.x = denseDw(lat, u, 1.0)
  (result.smin, result.smax) = denseHBounds(result.x, nd)
  let
    r31 = newRat(0.95*result.smin, 1.05*result.smax, 31)
    r11 = newRat(0.95*result.smin, 1.05*result.smax, 11)
  result.o31 = newOv(lat, 1.0, r31, r2in, r2out, mxit)
  result.o11 = newOv(lat, 1.0, r11, r2in, r2out, mxit)
  result.dov = denseOv(result.o31, u)
  echo &"  {name}: sigma = [{result.smin:.6f}, {result.smax:.6f}]" &
       &"  cond(X^dag X) = {(result.smax/result.smin)^2:.3f}"
  echo &"    window [{r31.smin:.6f}, {r31.smax:.6f}]:" &
       &" maxRelErr(31) = {r31.maxRelErr:.3e}  maxRelErr(11) = {r11.maxRelErr:.3e}"

let fixes = [makeFix("free", newGauge(lat)),
             makeFix("random", randGauge(lat, 20260821, 0.3))]

suite "fixtures and the dense oracle":

  test "dense polar factor is unitary and on the GW circle":
    for fx in fixes:
      # V = D_ov - 1 must satisfy V^dag V = 1 exactly (oracle self-check).
      var e = 0.0
      for j in 0..<nd:
        for i in 0..<nd:
          var s = complex64(0.0, 0.0)
          for k in 0..<nd:
            var vki = fx.dov[k + nd*i]
            var vkj = fx.dov[k + nd*j]
            if k == i: vki -= complex64(1.0, 0.0)
            if k == j: vkj -= complex64(1.0, 0.0)
            s += conjugate(vki)*vkj
          if i == j: s -= complex64(1.0, 0.0)
          e = max(e, abs(s))
      echo &"  {fx.name}: max |V^dag V - 1| = {e:.3e}"
      check e < 1e-12

suite "T1.3c  rational overlap vs the exact dense polar factor":

  test "applyOv == denseOv to 5 maxRelErr + 1e-10, orders 31 and 11":
    for fx in fixes:
      for (nm, o) in [("31", fx.o31), ("11", fx.o11)]:
        var worst = 0.0
        var y = newSpin(lat.nsite)
        var z = newSpin(lat.nsite)
        for s in 0..<6:
          let v = randSpin(lat.nsite, 3000 + s)
          applyOv(o, y, v, fx.u)
          denseApply(fx.dov, nd, z, v)
          var d = 0.0
          for i in 0..<lat.nsite:
            for c in 0..1: d += abs2(y[i][c] - z[i][c])
          worst = max(worst, sqrt(d/norm2(v)))
        let bound = 5.0*o.rat.maxRelErr + 1e-10
        echo &"  {fx.name} order {nm}: worst |D_rat v - D_exact v|/|v| = {worst:.3e}" &
             &"  (bound {bound:.3e})"
        check worst < bound
        check o.stats.ok

suite "T1.3d, T1.3e  Ginsparg-Wilson circle identities":
  ## The 4-component GW relation (IV.22), Gamma calD + calD Gamma = calD Gamma calD
  ## with calD = diag(D, D^dag), reduces for BOTH Gamma = gamma_4 (offdiag 1, 1) and
  ## Gamma = gamma_5 (offdiag -i, i) to the same pair of 2-component identities:
  ##   gamma_4:  LHS = offdiag(D + D^dag, D + D^dag),
  ##             RHS = offdiag(D D^dag, D^dag D)   (top-right, bottom-left)
  ##   gamma_5:  LHS = offdiag(-i(D + D^dag), i(D + D^dag)),
  ##             RHS = offdiag(-i D D^dag, i D^dag D)
  ## i.e. (IV.22)  <=>  D + D^dag - D^dag D = 0  AND  D + D^dag - D D^dag = 0,
  ## the circle identity (equivalent to (IV.17)) and its adjoint-order partner.
  ## For the exact polar factor both hold identically; the rational violates them
  ## by 1 - x R(x)^2 = -2e - e^2, i.e. at 2 maxRelErr plus solve residuals.

  test "D + D^dag - D^dag D = 0 and D + D^dag - D D^dag = 0":
    for fx in fixes:
      let o = fx.o31
      let v = randSpin(lat.nsite, 4100)
      var a = newSpin(lat.nsite)   # D v
      var b = newSpin(lat.nsite)   # D^dag v
      var c = newSpin(lat.nsite)   # D^dag D v
      var e = newSpin(lat.nsite)   # D D^dag v
      applyOv(o, a, v, fx.u)
      applyOvAdj(o, b, v, fx.u)
      applyOvAdj(o, c, a, fx.u)
      applyOv(o, e, b, fx.u)
      var r1 = 0.0
      var r2 = 0.0
      for i in 0..<lat.nsite:
        for s in 0..1:
          r1 += abs2(a[i][s] + b[i][s] - c[i][s])
          r2 += abs2(a[i][s] + b[i][s] - e[i][s])
      let n = sqrt(norm2(v))
      r1 = sqrt(r1)/n
      r2 = sqrt(r2)/n
      echo &"  {fx.name}: |D + Ddag - DdagD|/|v| = {r1:.3e}" &
           &"   |D + Ddag - DDdag|/|v| = {r2:.3e}"
      check r1 < 1e-13
      check r2 < 1e-13
      check o.stats.ok

suite "adjoint and gauge covariance":

  test "applyOvAdj is the true adjoint":
    for fx in fixes:
      let
        x = randSpin(lat.nsite, 4201)
        y = randSpin(lat.nsite, 4202)
      var dy = newSpin(lat.nsite)
      var dx = newSpin(lat.nsite)
      for mass in [0.0, 0.13]:
        applyOv(fx.o31, dy, y, fx.u, mass)
        applyOvAdj(fx.o31, dx, x, fx.u, mass)
        let
          a = dot(x, dy)
          b = dot(dx, y)
          e = abs(a - b)/sqrt(norm2(x)*norm2(dy))
        echo &"  {fx.name} mass={mass}: |<x,Dy> - <Ddag x,y>| rel = {e:.3e}"
        check e < 1e-12

  test "gauge covariance of D_ov":
    ## D_ov(theta + d alpha) e^{i alpha} = e^{i alpha} D_ov(theta): inherited from
    ## D_W, but the shifted solves must not break it -- this catches sign errors in
    ## the shifts.  The rational window is unchanged (conjugation by a unitary).
    let fx = fixes[1]
    let alpha = randReal(lat.nsite, 4301, 1.1)
    var ug = fx.u
    gaugeTransform(lat, ug, alpha)
    let x = randSpin(lat.nsite, 4302)
    var xg = x
    spinGaugeTransform(lat, xg, alpha)
    var a = newSpin(lat.nsite)
    var b = newSpin(lat.nsite)
    applyOv(fx.o31, a, x, fx.u, 0.2)
    spinGaugeTransform(lat, a, alpha)
    applyOv(fx.o31, b, xg, ug, 0.2)
    let e = reldiff(a, b)
    echo &"  |D_ov(theta + d alpha) e^ia x - e^ia D_ov(theta) x| rel = {e:.3e}"
    check e < 1e-12

suite "T1.3h  multishift path == independent single-shift solves":

  test "applyOv vs npole plain CG solves":
    for fx in fixes:
      let o = fx.o31
      let v = randSpin(lat.nsite, 4400)
      var y = newSpin(lat.nsite)
      applyOv(o, y, v, fx.u)
      # replicate with one cgSolve per pole at the same tolerance
      var z = newSpin(lat.nsite)
      z := v
      scale(z, o.rat.cst)
      var xj = newSpin(lat.nsite)
      for j in 0..<o.rat.npole:
        let ci = cgSolve(xj, v, r2in, mxit, hShiftOp(lat, fx.u, 1.0, o.rat.pole[j]))
        check ci.converged
        axpy(z, o.rat.res[j], xj)
      var w = newSpin(lat.nsite)
      applyDw(lat, w, z, fx.u, 1.0)      # X z
      axpy(w, 1.0, v)
      let e = reldiff(y, w)
      echo &"  {fx.name}: |applyOv - single-shift assembly| rel = {e:.3e}"
      check e < 1e-12

suite "ovGradient":

  test "tangent from the delta D_ov formula vs centered finite differences":
    ## T = <left, delta D_ov right>
    ##   = <left, dX z> - sum_j r_j ( <dX t_j, X s_j> + <X t_j, dX s_j> ),
    ## every ingredient built here from applyDw/applyDwAdj/applyDwDeriv and
    ## cgmSolve directly, not through ovGradient.
    for fx in fixes:
      let
        o = fx.o31
        n = lat.nsite
        left = randSpin(n, 4501)
        right = randSpin(n, 4502)
        du = randGauge(lat, 4503, 1.0)
      # s_j, t_j, z
      var ss, tt: seq[Spin]
      let mi1 = cgmSolve(ss, right, o.rat.pole, r2in, mxit, hShiftOp(lat, fx.u, 1.0, 0.0))
      var w = newSpin(n)
      applyDwAdj(lat, w, left, fx.u, 1.0)
      let mi2 = cgmSolve(tt, w, o.rat.pole, r2in, mxit, hShiftOp(lat, fx.u, 1.0, 0.0))
      check mi1.converged and mi2.converged
      var z = newSpin(n)
      z := right
      scale(z, o.rat.cst)
      for j in 0..<o.rat.npole: axpy(z, o.rat.res[j], ss[j])
      # T
      var dd = newSpin(n)
      var xv = newSpin(n)
      applyDwDeriv(lat, dd, z, fx.u, du)
      var tan = dot(left, dd)
      for j in 0..<o.rat.npole:
        applyDw(lat, xv, ss[j], fx.u, 1.0)        # X s_j
        applyDwDeriv(lat, dd, tt[j], fx.u, du)    # dX t_j
        tan -= o.rat.res[j]*dot(dd, xv)
        applyDw(lat, xv, tt[j], fx.u, 1.0)        # X t_j
        applyDwDeriv(lat, dd, ss[j], fx.u, du)    # dX s_j
        tan -= o.rat.res[j]*dot(xv, dd)
      # FD of <left, D_ov right>
      var up = newGauge(lat)
      var um = newGauge(lat)
      var a = newSpin(n)
      var b = newSpin(n)
      var best = 1e300
      var bestEps = 0.0
      for eps in [1e-3, 3e-4, 1e-4, 3e-5, 1e-5]:
        for i in 0..<fx.u.s.len:
          up.s[i] = fx.u.s[i] + eps*du.s[i]
          um.s[i] = fx.u.s[i] - eps*du.s[i]
        for i in 0..<fx.u.t.len:
          up.t[i] = fx.u.t[i] + eps*du.t[i]
          um.t[i] = fx.u.t[i] - eps*du.t[i]
        applyOv(o, a, right, up)
        applyOv(o, b, right, um)
        let fd = (dot(left, a) - dot(left, b))/(2.0*eps)
        let r = abs(fd - tan)/abs(tan)
        if r < best:
          best = r
          bestEps = eps
      echo &"  {fx.name}: best |dF_fd - dF|/|dF| = {best:.3e} at eps = {bestEps:.0e}"
      check best < 1e-8
      # pullback vs the same tangent contraction
      var f = newGauge(lat)
      ovGradient(o, f, left, right, fx.u)
      var lhs = 0.0
      for i in 0..<f.s.len: lhs += f.s[i]*du.s[i]
      for i in 0..<f.t.len: lhs += f.t[i]*du.t[i]
      let rhs = 2.0*tan.re
      let e = abs(lhs - rhs)/abs(rhs)
      echo &"  {fx.name}: <f,du> = {lhs:.12e}  2Re<l,dD r> = {rhs:.12e}  rel = {e:.3e}"
      check e < 1e-12
      check o.stats.ok

  test "scale and add semantics":
    let fx = fixes[1]
    let
      left = randSpin(lat.nsite, 4601)
      right = randSpin(lat.nsite, 4602)
    var f1 = newGauge(lat)
    var f2 = newGauge(lat)
    ovGradient(fx.o11, f1, left, right, fx.u, 2.0)
    ovGradient(fx.o11, f2, left, right, fx.u, 0.5)
    ovGradient(fx.o11, f2, left, right, fx.u, 1.5, add = true)
    var e = 0.0
    for i in 0..<f1.s.len: e = max(e, abs(f1.s[i] - f2.s[i]))
    for i in 0..<f1.t.len: e = max(e, abs(f1.t[i] - f2.t[i]))
    echo &"  |scale 2 - (0.5 + 1.5)| = {e:.3e}"
    check e < 1e-12

  test "Ward identity: pure-gauge du gives delta D_ov = i(alpha D_ov - D_ov alpha)":
    ## Sign convention: doc/04 section 7 (WP-E), theta_e -> theta_e + alpha_b -
    ## alpha_a, hence delta_{d alpha} D = i(alpha D - D alpha) -- the opposite sign
    ## to the section 15 item 6 sketch.  D_ov inherits it exactly (every ingredient
    ## of the rational is covariant), up to solve residuals.
    for fx in fixes:
      let
        o = fx.o31
        n = lat.nsite
        left = randSpin(n, 4701)
        right = randSpin(n, 4702)
        alpha = randReal(n, 4703, 1.0)
      var da = newGauge(lat)
      gaugeTransform(lat, da, alpha)        # da = d alpha
      var f = newGauge(lat)
      ovGradient(o, f, left, right, fx.u)
      var lhs = 0.0
      for i in 0..<f.s.len: lhs += f.s[i]*da.s[i]
      for i in 0..<f.t.len: lhs += f.t[i]*da.t[i]
      # 2 Re <left, i(alpha D_ov - D_ov alpha) right>
      var dr = newSpin(n)
      applyOv(o, dr, right, fx.u)           # D right
      var ar = newSpin(n)
      for i in 0..<n:
        for c in 0..1: ar[i][c] = alpha[i]*right[i][c]
      var dar = newSpin(n)
      applyOv(o, dar, ar, fx.u)             # D (alpha right)
      var rhs = 0.0
      for i in 0..<n:
        for c in 0..1:
          let z = alpha[i]*dr[i][c] - dar[i][c]
          rhs += 2.0*(left[i][c].re*(-z.im) + left[i][c].im*z.re)
      let e = abs(lhs - rhs)/abs(rhs)
      echo &"  {fx.name}: <f,d alpha> = {lhs:.12e}  Ward rhs = {rhs:.12e}  rel = {e:.3e}"
      check e < 1e-12

suite "T1.3f  parity and time reversal (frame-independent consequences)":

  test "free D_lat spectrum is closed under complex conjugation":
    let ev = eigvals(fixes[0].x, nd)        # X = D_W - 1; conj closure is the same
    let e = conjClosure(ev)
    echo &"  max dist(conj(lambda), spec) = {e:.3e}"
    check e < 1e-10

  test "free D_ov spectrum: conjugation-closed and on the GW circle":
    let fx = fixes[0]
    # the rational operator we actually use, assembled column by column
    let a = assembleOv(fx.o31, fx.u)
    let ev = eigvals(a, nd)
    var circ = 0.0
    for z in ev: circ = max(circ, abs(abs(z - complex64(1.0, 0.0)) - 1.0))
    let cc = conjClosure(ev)
    # the exact dense polar factor: circle to machine precision
    let eve = eigvals(fx.dov, nd)
    var circE = 0.0
    for z in eve: circE = max(circE, abs(abs(z - complex64(1.0, 0.0)) - 1.0))
    let bound = 5.0*fx.o31.rat.maxRelErr + 1e-9
    echo &"  rational: max ||lam-1|-1| = {circ:.3e} (bound {bound:.3e})" &
         &"  conj closure = {cc:.3e}"
    echo &"  exact:    max ||lam-1|-1| = {circE:.3e}"
    check circ < bound
    check cc < 1e-8
    check circE < 1e-12

  test "time reversal: free overlap propagator folds about T/2, Wilson does not":
    ## Paper Fig. 10.  With the antiperiodic BC the free two-point function at a
    ## coincident spatial site satisfies G(t) = G(nt - t) (the continuum
    ## G_T(T - t) = G_T(t) fold; component (1,1), source component 0).  The overlap
    ## inherits it through D_ov^{-1} = 1 - (D_ov^dag)^{-1} (IV.17); the Wilson term
    ## breaks it at O(a_t) for D_W.
    let
      latP = newLat(sph, 8, 0.25)
      ndP = 2*latP.nsite
      uP = newGauge(latP)
      xP = denseDw(latP, uP, 1.0)
      (pmin, pmax) = denseHBounds(xP, ndP)
      oP = newOv(latP, 1.0, newRat(0.95*pmin, 1.05*pmax, 31), r2in, r2out, mxit)
    let b = pointSource(latP.nsite, 0, 0)   # site (v = 0, t = 0), component 0
    # overlap: x = D_ov^{-1} b = (Ddag D)^{-1} Ddag b
    var rhs = newSpin(latP.nsite)
    applyOvAdj(oP, rhs, b, uP)
    var x = newSpin(latP.nsite)
    let ci = solveNormal(oP, x, rhs, uP)
    check ci.converged
    # Wilson: same source, same site, via its own normal equations
    var rhsW = newSpin(latP.nsite)
    applyDwAdj(latP, rhsW, b, uP)
    var xw = newSpin(latP.nsite)
    let cw = cgSolve(xw, rhsW, 1e-24, mxit, hShiftOp(latP, uP, 0.0, 0.0))
    check cw.converged
    var gv = newSeq[Complex64](latP.nt)
    var gw = newSeq[Complex64](latP.nt)
    for t in 0..<latP.nt:
      gv[t] = x[sIdx(latP, 0, t)][0]
      gw[t] = xw[sIdx(latP, 0, t)][0]
    var scOv = 0.0
    var scW = 0.0
    var dOv = 0.0
    var dW = 0.0
    for t in 1..<latP.nt:
      scOv = max(scOv, abs(gv[t]))
      scW = max(scW, abs(gw[t]))
      dOv = max(dOv, abs(gv[t] - gv[latP.nt - t]))
      dW = max(dW, abs(gw[t] - gw[latP.nt - t]))
    echo &"  overlap: max |G(t) - G(nt-t)| / max |G| = {dOv/scOv:.3e}"
    echo &"  Wilson:  max |G(t) - G(nt-t)| / max |G| = {dW/scW:.3e}" &
         &"   (ratio {dW/scW/(dOv/scOv):.2e})"
    check dOv/scOv < 1e-7          ## rational + solver accuracy
    check dW/scW > 1e-2            ## the physical O(a_t) violation
    check dW/scW > 1e3*(dOv/scOv)  ## the Fig. 10 discriminator

suite "kernelWindow":

  test "reproduces the dense sigma bounds within its margins":
    for fx in fixes:
      let kw = kernelWindow(fx.o31, fx.u, 400)
      echo &"  {fx.name}: smin = {kw.smin:.8f} (dense {fx.smin:.8f})" &
           &"  smax = {kw.smax:.8f} (dense {fx.smax:.8f})"
      echo &"    margins [{kw.lo:.8f}, {kw.hi:.8f}]  inside = {kw.inside}"
      # Rayleigh quotients bound the extremes one-sidedly; the residual expansion
      # must cover them from the other side.
      check kw.smin >= fx.smin*(1.0 - 1e-10)
      check kw.lo <= fx.smin*(1.0 + 1e-10)
      check kw.smax <= fx.smax*(1.0 + 1e-10)
      check kw.hi >= fx.smax*(1.0 - 1e-10)
      # and the point estimates are close after 400 iterations
      check abs(kw.smin - fx.smin) < 1e-6*fx.smin
      check abs(kw.smax - fx.smax) < 1e-3*fx.smax
      check kw.inside

  test "a deliberately narrow window reports inside = false":
    let fx = fixes[0]
    let on = newOv(lat, 1.0, newRat(1.5*fx.smin, 1.05*fx.smax, 11), r2in, r2out, mxit)
    let kw = kernelWindow(on, fx.u, 400)
    echo &"  narrow window [{on.rat.smin:.6f}, {on.rat.smax:.6f}]:" &
         &" lo = {kw.lo:.6f}  inside = {kw.inside}"
    check not kw.inside

suite "solve counts and allocation":

  test "solves per call":
    let fx = fixes[1]
    let o = fx.o31
    let v = randSpin(lat.nsite, 4901)
    var y = newSpin(lat.nsite)
    var z = newSpin(lat.nsite)
    o.clearStats
    applyOv(o, y, v, fx.u)
    echo &"  applyOv:    nmulti = {o.stats.nmulti}  miters = {o.stats.miters}" &
         &"  nx = {o.stats.nx}"
    check o.stats.nmulti == 1
    o.clearStats
    applyOvAdj(o, y, v, fx.u)
    check o.stats.nmulti == 1
    o.clearStats
    applyNormal(o, y, v, fx.u)
    check o.stats.nmulti == 2
    o.clearStats
    var f = newGauge(lat)
    ovGradient(o, f, v, y, fx.u)
    echo &"  ovGradient: nmulti = {o.stats.nmulti}  miters = {o.stats.miters}" &
         &"  nx = {o.stats.nx}"
    check o.stats.nmulti == 2
    o.clearStats
    let ci = solveNormal(o, z, v, fx.u)
    echo &"  solveNormal: iters = {ci.iters}  nmulti = {o.stats.nmulti}" &
         &"  (= 2 (iters + 1) = {2*(ci.iters + 1)})"
    check ci.converged
    check o.stats.nmulti == 2*(ci.iters + 1)
    check o.stats.ok

  test "allocation regression over 64 applies":
    ## The Ov-owned path allocates nothing: work/xs/xt are preallocated and the
    ## solver closures are built once in newOv.  cgmSolve/cgSolve DO allocate their
    ## nshift+3 scratch fields per call (a documented WP-D property of the stateless
    ## section-9 signatures), all of it garbage on return; under this build's refc
    ## GC the raw occupied counter therefore drifts with the collection cadence.
    ## The assertion with teeth is live memory across GC_fullCollect: it must be
    ## identical, i.e. nothing in the apply path GROWS the live set.  The raw drift
    ## is echoed so the churn stays visible.
    let fx = fixes[1]
    let o = fx.o11
    var x = randSpin(lat.nsite, 4902)
    var y = newSpin(lat.nsite)
    var f = newGauge(lat)
    # warm up every path
    for k in 0..<4:
      applyOv(o, y, x, fx.u, 0.1)
      applyOvAdj(o, x, y, fx.u, 0.1)
      applyNormal(o, y, x, fx.u, 0.2)
      ovGradient(o, f, x, y, fx.u, 1e-9, add = true)
      let s = 1.0/sqrt(norm2(x) + norm2(y))
      scale(x, s)
      scale(y, s)
    GC_fullCollect()
    let raw0 = getRawMemAllocated()
    let occ0 = getOccupiedMem()
    var acc = 0.0
    for k in 0..<64:
      applyOv(o, y, x, fx.u, 0.1)
      applyOvAdj(o, x, y, fx.u, 0.1)
      applyNormal(o, y, x, fx.u, 0.2)
      ovGradient(o, f, x, y, fx.u, 1e-9, add = true)
      let s = 1.0/sqrt(norm2(x) + norm2(y))
      scale(x, s)
      scale(y, s)
      acc += s
    let drift = getOccupiedMem() - occ0
    GC_fullCollect()
    let raw1 = getRawMemAllocated()
    let occ1 = getOccupiedMem()
    echo &"  raw {raw0} -> {raw1}   live occupied {occ0} -> {occ1}" &
         &"   uncollected drift {drift}   acc = {acc:.6e}"
    check raw1 == raw0
    check occ1 == occ0
    check acc != 0.0
    check o.stats.ok

suite "production window report (L = 1 free field, M = 1)":

  test "sigma bounds and rational errors at a_t = 0.2, L_t = 60":
    ## Sizes the production windows (doc/03 tier 2 runs at a_t = 0.2, T = 12).
    ## Free field: X block-diagonalizes on the antiperiodic Matsubara modes,
    ## X(k) = D_spatial + kappa'(1 - cos k + i sigma_3 sin k) - M, so the sigma
    ## bounds come from 60 dense 24 x 24 problems.  Cross-check: min |eig(X)|
    ## must reproduce WP-E's raw min |D_W - 1| = 1.1234 at L = 1, Lt = 60.
    const
      at = 0.2
      nt = 60
    let
      l1 = newLat(sph, 1, at)
      nv = sph.nv
      nb = 2*nv
      dsp = denseDw(l1, newGauge(l1), 0.0, dwSpace)
    var
      hmin = 1e300
      hmax = 0.0
      xmin = 1e300
    for n in 0..<nt:
      let k = PI*float(2*n + 1)/float(nt)
      var x = dsp
      for y in 0..<nv:
        let kt = l1.kapT[y]
        x[(2*y) + nb*(2*y)] += complex64(kt*(1.0 - cos(k)) - 1.0, kt*sin(k))
        x[(2*y + 1) + nb*(2*y + 1)] += complex64(kt*(1.0 - cos(k)) - 1.0, -kt*sin(k))
      for z in eigvals(x, nb): xmin = min(xmin, abs(z))
      let (a, b) = denseHBounds(x, nb)
      hmin = min(hmin, a)
      hmax = max(hmax, b)
    let
      r31 = newRat(0.95*hmin, 1.05*hmax, 31)
      r11 = newRat(0.95*hmin, 1.05*hmax, 11)
    echo &"  sigma = [{hmin:.6f}, {hmax:.6f}]  cond(X^dag X) = {(hmax/hmin)^2:.2f}"
    echo &"  window [{r31.smin:.6f}, {r31.smax:.6f}]:" &
         &" maxRelErr(31) = {r31.maxRelErr:.3e}  maxRelErr(11) = {r11.maxRelErr:.3e}"
    # 1.1770 under the exact-kappa convention (doc/06 "THE COUPLING CONVENTION");
    # the flat-kappa value was 1.1234 (= WP-E's raw min |D_W - 1|, slide-8 legend
    # convention), and the shift +0.0536 was predicted by the Python oracle.
    echo &"  min |eig(X)| = {xmin:.6f}  (flat-kappa value was 1.1234)"
    check abs(xmin - 1.1770) < 2e-4
    check hmin > 0.0 and hmax > hmin
