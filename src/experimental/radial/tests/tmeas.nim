#RUNCMD $RUN1
## WP-I acceptance tests: harmonics (Gram matrices, icosahedral protection),
## propagator/condensate against dense oracles, the conserved-current Ward
## test, sigma_PS vs sigma_FS, the gluonic J_top pipeline (exact + Monte
## Carlo), the 7-shape Wilson-loop basis and the GEVP.
##
## References: doc/07-observables.md, doc/04-interfaces.md section 14,
## doc/06-status.md WP-I section (measured numbers recorded there).

import std/[math, complex, strformat, unittest, algorithm]
import eigens/linalgFuncs
import ../meas/observables
import ../meas/gevp
import ../meas/fit
import ../ops/flow

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))

# --- helpers -------------------------------------------------------------------

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

proc denseHBounds(x: seq[Complex64], nd: int): tuple[smin, smax: float] =
  var h = newSeq[Complex64](nd*nd)
  for j in 0..<nd:
    for i in j..<nd:
      var s = complex64(0.0, 0.0)
      for k in 0..<nd: s += conjugate(x[k + nd*i])*x[k + nd*j]
      h[i + nd*j] = s
      h[j + nd*i] = conjugate(s)
  var ev = newSeq[float](nd)
  zeigs(cast[ptr float64](addr h[0]), addr ev[0], nd)
  (sqrt(ev[0]), sqrt(ev[nd-1]))

proc symEig(m: seq[float], n: int): seq[float] =
  ## Ascending eigenvalues of a real symmetric matrix (row-major n x n).
  var a = newSeq[Complex64](n*n)
  for i in 0..<n:
    for j in 0..<n:
      a[i + n*j] = complex64(0.5*(m[i*n + j] + m[j*n + i]), 0.0)
  result = newSeq[float](n)
  zeigs(cast[ptr float64](addr a[0]), addr result[0], n)

proc clusters(ev: seq[float]): tuple[k: int, s1, s2, gap: float] =
  ## Best split of the sorted eigenvalues into two clusters: the split index
  ## with the largest gap; s1/s2 are the within-cluster spreads.
  var e = ev
  e.sort
  var best = 0.0
  var kb = 1
  for k in 1..<e.len:
    let g = e[k] - e[k-1]
    if g > best:
      best = g
      kb = k
  result.k = kb
  result.gap = best
  result.s1 = e[kb-1] - e[0]
  result.s2 = e[e.len-1] - e[kb]

func offDiagId(g: seq[float], n: int): tuple[off, spread, mean: float] =
  ## For an n x n block: max |offdiag|, diagonal spread, diagonal mean.
  var
    dmin = 1e300
    dmax = -1e300
  for i in 0..<n:
    for j in 0..<n:
      if i == j:
        dmin = min(dmin, g[i*n + j])
        dmax = max(dmax, g[i*n + j])
      else:
        result.off = max(result.off, abs(g[i*n + j]))
      result.mean += (if i == j: g[i*n + j]/float(n) else: 0.0)
  result.spread = dmax - dmin

func maxAbs(g: seq[float]): float =
  for x in g: result = max(result, abs(x))

# ==============================================================================
suite "harmonics: addition theorem and quadrature Gram matrices":

  test "ylm vs legendreP addition theorem, l <= 4":
    var r: Threefry4x64
    r.seedIndep(101, 0)
    var worst = 0.0
    for rep in 0..<20:
      var u = [r.gaussian, r.gaussian, r.gaussian].unit
      var v = [r.gaussian, r.gaussian, r.gaussian].unit
      for l in 0..4:
        var s = 0.0
        for m in -l..l: s += ylm(l, m, u)*ylm(l, m, v)
        let p = float(2*l + 1)/(4.0*PI)*legendreP(l, dot(u, v))
        worst = max(worst, abs(s - p))
    echo &"  max |sum_m Y(u)Y(v) - (2l+1)/(4pi) P_l(u.v)| = {worst:.3e}"
    check worst < 1e-13

  test "Gram matrices at L = 1, 2, 4: exact protections and O(a^2) values":
    for lev in [1, 2, 4]:
      let sph = newSphere(lev)
      echo &"  L = {lev}:"
      for (nm, gr) in [("site", siteGram), ("face", faceGram)]:
        # diagonal blocks
        for l in 0..4:
          let
            n = 2*l + 1
            g = gr(sph, l, l)
            (off, spr, mean) = offDiagId(g, n)
          if l <= 2:
            # Schur: exactly proportional to the identity on any lattice
            check off < 1e-12
            check spr < 1e-12
            echo &"    {nm} l={l}: diag = {mean:.15f}  (|off| {off:.1e}, spread {spr:.1e})"
            if lev == 1:
              check abs(mean - 1.0) < 1e-13   # icosahedral orbits are 5-designs
          else:
            let ev = symEig(g, n)
            let c = clusters(ev)
            echo &"    {nm} l={l}: eig clusters {c.k}+{n-c.k}: " &
                 &"[{ev.min:.6f} .. {ev.max:.6f}] spreads {c.s1:.1e}/{c.s2:.1e} gap {c.gap:.4f}"
            if l == 3:
              check c.k == 3 or c.k == 4          # T2 (3) + G (4)
              check max(c.s1, c.s2) < 1e-10*max(1.0, maxAbs(g))
        # cross blocks that vanish by Schur (no common I irrep):
        # (0,l<=4 except pairs with A: none), (1,2), (1,3), (1,4), (2,3)
        var worst = 0.0
        for (l1, l2) in [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3),
                         (1, 4), (2, 3)]:
          worst = max(worst, maxAbs(gr(sph, l1, l2)))
        check worst < 1e-12
        # allowed couplings between equivalent irreps: report only
        let
          c24 = maxAbs(gr(sph, 2, 4))
          c34 = maxAbs(gr(sph, 3, 4))
        echo &"    {nm} Schur-zero cross max {worst:.1e};  (2,4) {c24:.3e}  (3,4) {c34:.3e}"

# ==============================================================================
suite "GEVP on synthetic data":

  test "recovers the exact energies through zeigsgv":
    const
      nop = 4
      nt = 8
      at = 0.25
    let e = [0.7, 1.3, 2.1, 3.4]
    var v: array[nop, array[nop, float]]
    var r: Threefry4x64
    r.seedIndep(202, 0)
    for i in 0..<nop:
      for n in 0..<nop: v[i][n] = r.gaussian + (if i == n: 2.0 else: 0.0)
    var c = newSeq[seq[seq[float]]](nt)
    for t in 0..<nt:
      c[t] = newSeq[seq[float]](nop)
      for i in 0..<nop:
        c[t][i] = newSeq[float](nop)
        for j in 0..<nop:
          for n in 0..<nop:
            c[t][i][j] += v[i][n]*v[j][n]*exp(-e[n]*at*float(t))
    let hc = gevpCheck(c, 1)
    echo &"  C(t0=1): evmin {hc.evmin:.3e}  cond {hc.cond:.3e}  asym {hc.asym:.1e}"
    check hc.evmin > 0.0
    let dims = gevpDims(c, 1, at)
    var worst = 0.0
    for t in 1..<nt-1:
      for n in 0..<nop:
        worst = max(worst, abs(dims[t][n] - e[n]))
    echo &"  max |Delta_n(t) - E_n| over t in [1,{nt-2}] = {worst:.3e}"
    check worst < 1e-10

# ==============================================================================
# fermion fixture: L = 1, nt = 6, at = 0.4 (dense dimension 144), free + random

const
  ## Solver targets for the measurement fixtures.  1e-26/1e-22 (the WP-F test
  ## values) sit at the convergence-guard floor for LONG solve chains on the
  ## random window: measured with them, the outer cgSolve's recomputed-residual
  ## guard trips marginally (ok=false, mrefits=0) and the multishift refinement
  ## fires thousands of times over the noise chains.  One decade of headroom
  ## fixes both; every acceptance tolerance below is unchanged and still met
  ## with orders of magnitude to spare.
  ## (Raised again 1e-25 -> 1e-24 after the exact-kappa convention switch
  ## (doc/06 "THE COUPLING CONVENTION") nudged the roundoff floor: same marginal
  ## ok=false, mrefits=0 signature, same one-decade fix.)
  r2in = 1e-24
  r2out = 1e-20
  mxit = 20000

let
  sph1 = newSphere(1)
  latF = newLat(sph1, 6, 0.4)
  ndF = 2*latF.nsite
  v5 = fiveFoldSite(sph1)

proc makeOv(u: Gauge, order: int): Ov =
  let
    x = denseDw(latF, u, 1.0)
    (smin, smax) = denseHBounds(x, ndF)
  newOv(latF, 1.0, newRat(0.95*smin, 1.05*smax, order), r2in, r2out, mxit)

let
  uF0 = newGauge(latF)
  uFr = randGauge(latF, 20260821, 0.3)
  ovF0 = makeOv(uF0, 31)
  ovFr = makeOv(uFr, 31)

echo &"fermion fixture: L=1 nt=6 at=0.4, 5-fold site {v5}, " &
     &"windows [{ovF0.rat.smin:.4f},{ovF0.rat.smax:.4f}] / " &
     &"[{ovFr.rat.smin:.4f},{ovFr.rat.smax:.4f}]"

suite "fermion propagator and condensate vs dense oracles":

  test "propagatorT == dense inverse column (both fields, mass 0 and 0.13)":
    for (nm, u, o) in [("free", uF0, ovF0), ("random", uFr, ovFr)]:
      for mass in [0.0, 0.13]:
        o.clearStats
        let
          s = denseS(o, u, mass)
          g = propagatorT(o, u, mass, v5)
          j0 = 2*sIdx(latF, v5, 0)
        var worst = 0.0
        var scale = 0.0
        for t in 0..<latF.nt:
          for c in 0..1:
            let d = s[(2*sIdx(latF, v5, t) + c) + ndF*j0]
            worst = max(worst, abs(g[t][c] - d))
            scale = max(scale, abs(d))
        echo &"  {nm} mass={mass}: max |G_solver - G_dense| = {worst:.3e}  (scale {scale:.3e})"
        echo &"    stats: nmulti {o.stats.nmulti} miters {o.stats.miters} " &
             &"mrefits {o.stats.mrefits} ncg {o.stats.ncg} cgiters {o.stats.cgiters} ok {o.stats.ok}"
        check worst < 1e-10
        check o.stats.ok

  test "condensatePS stochastic vs condensateDense":
    var r: Threefry4x64
    r.seedIndep(303, 0)
    for (nm, u, o, mass) in [("free", uF0, ovF0, 0.08),
                             ("random", uFr, ovFr, 0.08)]:
      let
        d = condensateDense(o, u, mass)
        (v, e) = condensatePS(o, u, mass, 24, r)
        pull = abs(v - d)/max(e, 1e-30)
      echo &"  {nm} mass={mass}: dense {d:.8f}  noise {v:.8f} +- {e:.8f}  pull {pull:.2f}"
      check abs(v - d) < 6.0*e + 1e-8

suite "the current kernel: dense tangent pinned against ovGradient":

  test "2 Re<l, denseOvDeriv r> == <ovGradient, du> (both fields and masses)":
    for (nm, u, o) in [("free", uF0, ovF0), ("random", uFr, ovFr)]:
      for mass in [0.0, 0.17]:
        let
          du = randGauge(latF, 404, 1.0)
          left = randSpin(latF.nsite, 405)
          right = randSpin(latF.nsite, 406)
          dd = denseOvDeriv(o, u, du, mass)
        var rhs = complex64(0.0, 0.0)
        for i in 0..<ndF:
          var s = complex64(0.0, 0.0)
          for j in 0..<ndF: s += dd[i + ndF*j]*right[j shr 1][j and 1]
          rhs += conjugate(left[i shr 1][i and 1])*s
        var f = newGauge(latF)
        ovGradient(o, f, left, right, u, mass = mass)
        var lhs = 0.0
        for i in 0..<f.s.len: lhs += f.s[i]*du.s[i]
        for i in 0..<f.t.len: lhs += f.t[i]*du.t[i]
        let e = abs(lhs - 2.0*rhs.re)/abs(lhs)
        echo &"  {nm} mass={mass}: <f,du> = {lhs:.12e}  " &
             &"2Re<l,dD r> = {2.0*rhs.re:.12e}  rel {e:.3e}"
        check e < 1e-11

suite "Ward / charge conservation (measurement level)":

  test "total charge in a fermion line: exact plateaus and unit-charge jump":
    for (nm, u, o, mass) in [("random m=0", uFr, ovFr, 0.0),
                             ("random m=0.1", uFr, ovFr, 0.1)]:
      const
        ta = 1
        tb = 4
      o.clearStats
      let (c, jump) = wardChargeScan(o, u, mass, v5, ta, tb)
      var scale = abs(jump)
      for z in c: scale = max(scale, abs(z))
      # plateaus: t in {ta..tb-1} and {tb..nt-1, 0..ta-1}
      var flat = 0.0
      for t in ta..<tb-1: flat = max(flat, abs(c[t+1] - c[t]))
      var t = tb
      while (t + 1) mod latF.nt != ta:
        flat = max(flat, abs(c[(t+1) mod latF.nt] - c[t mod latF.nt]))
        inc t
      let
        j1 = c[ta-1] - c[ta]        # = +i S_ba
        j2 = c[tb-1] - c[tb]        # = -i S_ba
        ej = max(abs(j1 - jump), abs(j2 + jump))
      echo &"  {nm}: plateau flatness {flat/scale:.3e}  " &
           &"jump err {ej/scale:.3e}  (|jump| {abs(jump):.3e}, scale {scale:.3e})"
      check flat/scale < 1e-10
      check ej/scale < 1e-10
      check o.stats.ok

  test "l=0 current-current correlator: dense probe of the contact structure":
    ## The current-current connected correlator C(t2,t1) = -2Re tr[K(t2)SK(t1)S]
    ## with the OVERLAP current: the one-point trace tr[K(t)S] is exactly
    ## t-independent (pure gauge variance), and the deviation of C(t2,t1) from
    ## t2-independence at t2 != t1 measures the GW-smeared contact term of the
    ## non-ultralocal kernel.  Measured and reported; the exact statement
    ## asserted is the one-point invariance and the fermion-line test above.
    let
      u = uFr
      o = ovFr
      nt = latF.nt
      s = denseS(o, u, 0.0)
    var w1 = newSeq[float](sph1.nv)
    for y in 0..<w1.len: w1[y] = 1.0
    var kt = newSeq[seq[Complex64]](nt)
    for t in 0..<nt:
      kt[t] = denseOvDeriv(o, u, tsliceForm(latF, w1, t))
    # one-point traces tr[K(t) S]
    var tr1 = newSeq[Complex64](nt)
    for t in 0..<nt:
      var z = complex64(0.0, 0.0)
      for i in 0..<ndF:
        for j in 0..<ndF: z += kt[t][i + ndF*j]*s[j + ndF*i]
      tr1[t] = z
    var dev1 = 0.0
    for t in 1..<nt: dev1 = max(dev1, abs(tr1[t] - tr1[0]))
    echo &"  one-point tr[K(t)S]: value {tr1[0].re:.6e}{tr1[0].im:+.3e}i  " &
         &"max t-dev {dev1:.3e}"
    check dev1 < 1e-10*max(1.0, abs(tr1[0]))
    # two-point: C(t2; t1=0)
    let t1 = 0
    var sk = newSeq[Complex64](ndF*ndF)
    block:
      var a = newSeq[Complex64](ndF*ndF)
      # a = K(t1) S ; sk = S a
      for j in 0..<ndF:
        for k in 0..<ndF:
          let f = s[k + ndF*j]
          if f.re != 0.0 or f.im != 0.0:
            for i in 0..<ndF: a[i + ndF*j] += kt[t1][i + ndF*k]*f
      for j in 0..<ndF:
        for k in 0..<ndF:
          let f = a[k + ndF*j]
          if f.re != 0.0 or f.im != 0.0:
            for i in 0..<ndF: sk[i + ndF*j] += s[i + ndF*k]*f
    var cc = newSeq[float](nt)
    for t2 in 0..<nt:
      var z = complex64(0.0, 0.0)
      for i in 0..<ndF:
        for j in 0..<ndF: z += kt[t2][i + ndF*j]*sk[j + ndF*i]
      cc[t2] = -2.0*z.re
    var line = "  conn C(t2;0):"
    for t2 in 0..<nt: line &= &" {cc[t2]:.6e}"
    echo line
    var dev2 = 0.0
    for t2 in 1..<nt-1: dev2 = max(dev2, abs(cc[t2+1] - cc[t2]))
    echo &"  t2-dependence away from contact (t2 = 1..{nt-1}): {dev2:.3e} " &
         &"of scale {abs(cc[1]):.3e} -- the smeared GW contact term, reported not asserted"

suite "connected and disconnected current estimators vs dense":

  test "factorized stochastic estimator reproduces the dense traces":
    ## Two slice operators: the l = 0 total charge (w = 1) and the Y_10
    ## projection (w_y = Y_10(pos_y)), the actual current operators of
    ## doc/07 1.3.  Dense reference for both from denseOvDeriv + denseS.
    let
      u = uFr
      o = ovFr
      nt = latF.nt
      s = denseS(o, u, 0.0)
    var w1 = newSeq[float](sph1.nv)
    for y in 0..<w1.len: w1[y] = 1.0
    var wy = newSeq[float](sph1.nv)
    for y in 0..<wy.len: wy[y] = ylm(1, 0, sph1.pos[y])
    let wops = [w1, wy]
    # dense reference: tr[K_a(t2) S K_b(0) S] and tr[K_a(t) S]
    var kt = newSeq[seq[Complex64]](2*nt)
    for a in 0..1:
      for t in 0..<nt:
        kt[a*nt + t] = denseOvDeriv(o, u, tsliceForm(latF, wops[a], t))
    var
      dconn = newSeq[Complex64](2*nt)   # [a*nt+t2] vs source op a at t1 = 0
      done1 = newSeq[Complex64](2*nt)
    for b in 0..1:
      var a = newSeq[Complex64](ndF*ndF)
      var sk = newSeq[Complex64](ndF*ndF)
      for j in 0..<ndF:
        for k in 0..<ndF:
          let f = s[k + ndF*j]
          if f.re != 0.0 or f.im != 0.0:
            for i in 0..<ndF: a[i + ndF*j] += kt[b*nt][i + ndF*k]*f
      for j in 0..<ndF:
        for k in 0..<ndF:
          let f = a[k + ndF*j]
          if f.re != 0.0 or f.im != 0.0:
            for i in 0..<ndF: sk[i + ndF*j] += s[i + ndF*k]*f
      for t in 0..<nt:
        var z1 = complex64(0.0, 0.0)
        var z2 = complex64(0.0, 0.0)
        for i in 0..<ndF:
          for j in 0..<ndF:
            z1 += kt[b*nt + t][i + ndF*j]*sk[j + ndF*i]
            z2 += kt[b*nt + t][i + ndF*j]*s[j + ndF*i]
        dconn[b*nt + t] = z1
        done1[b*nt + t] = z2
    # stochastic samples.  The estimator error is statistical (~20 % here), so
    # the sampling Ov runs at 1e-23/1e-16: over ~8000 multishift solves on
    # noise right-hand sides the 1.001-guard of the UNREFINABLE seed system
    # (cgmSolve j = 0; refining it would reproduce the same vector, WP-D) trips
    # at 1e-25, and the outer normal-equation recompute at 1e-20 -- measured,
    # ok=false with mrefits=0 in both cases.  Solver error ~3e-12 is 10 orders
    # below the noise.
    var r: Threefry4x64
    r.seedIndep(505, 0)
    const npair = 128
    let oc = newOv(latF, 1.0, o.rat, 1e-23, 1e-16, mxit)
    var samples = newSeq[CurrentSample](npair)
    for k in 0..<npair:
      samples[k] = currentSample(oc, u, 0.0, wops, r)
    echo &"  stats over {npair} samples: nmulti {oc.stats.nmulti} mrefits {oc.stats.mrefits} " &
         &"ncg {oc.stats.ncg} ok {oc.stats.ok}"
    check oc.stats.ok
    for a in 0..1:
      var worstC = 0.0
      var worstD = 0.0
      for t in 0..<nt:
        # sink op a at t2 = t against source op a at t1 = 0 (same-op correlator)
        let (cv, ce) = currentCorrConn(samples, a*nt + t, a*nt)
        let pr = abs(cv.re - dconn[a*nt + t].re)/max(ce.re, 1e-30)
        let pi = abs(cv.im - dconn[a*nt + t].im)/max(ce.im, 1e-30)
        worstC = max(worstC, max(pr, pi))
        if t == 0 or t == 2:
          echo &"  op{a} conn t2={t}: dense {dconn[a*nt+t].re:.6e}  " &
               &"noise {cv.re:.6e} +- {ce.re:.2e}  pull {pr:.2f}"
        # one-noise trace estimator (disconnected building block)
        var
          sm = 0.0
          s2 = 0.0
        for k in 0..<npair:
          sm += samples[k].d[a*nt + t].re
          s2 += samples[k].d[a*nt + t].re*samples[k].d[a*nt + t].re
        sm /= float(npair)
        let se = sqrt(max(0.0, s2/float(npair) - sm*sm)/float(npair - 1))
        worstD = max(worstD, abs(sm - done1[a*nt + t].re)/max(se, 1e-30))
      echo &"  op{a} worst pull: connected {worstC:.2f}  one-point {worstD:.2f}  (npair {npair})"
      check worstC < 6.0
      check worstD < 6.0
    # disconnected correlator piece: cross-noise products vs the dense product
    # 2Re tr[K(t2)S] * 2Re tr[K(0)S] (a fixed-configuration statement)
    var worstX = 0.0
    for t in [0, 2]:
      let (dv, de) = currentTraceDisc(samples, t, 0)
      let dref = 4.0*done1[t].re*done1[0].re
      worstX = max(worstX, abs(dv - dref)/max(de, 1e-30))
      if t == 0:
        echo &"  disc t2={t}: dense {dref:.6e}  noise {dv:.6e} +- {de:.2e}"
    echo &"  disc worst pull {worstX:.2f} (stderr from correlated pairs: optimistic)"
    check worstX < 10.0

suite "sigma_PS vs sigma_FS":

  test "connected timeslice correlators are IDENTICAL at every dt (mass 0)":
    ## Derivation (doc/06 WP-I): with the contractions <xi eta^dag> = S,
    ## <eta xi^dag> = S^dag, the two correlators differ only through
    ## (1 - D^dag) S^dag = S^dag - 1.  At dt != 0 every extra term carries a
    ## slice overlap P_{t2} P_{t1} = 0, and at dt = 0 the would-be contact
    ## 2 sum_{x in t} Re tr S^dag_xx - 2 nv vanishes too, because the GW
    ## relation (IV.17) reads S + S^dag = 1 at mass 0, so the site-diagonal
    ## blocks obey tr(S + S^dag)_xx = 2 exactly.  So PS == FS at ALL dt,
    ## configuration by configuration -- stronger than "identical spectra".
    for (nm, u, o) in [("free", uF0, ovF0), ("random", uFr, ovFr)]:
      o.clearStats
      let (ps, fs) = scalarCorrDense(o, u, 0.0)
      var scale = 0.0
      for x in ps: scale = max(scale, abs(x))
      var dev = 0.0
      for dt in 0..<latF.nt: dev = max(dev, abs(ps[dt] - fs[dt]))
      echo &"  {nm} dense: max|PS-FS| over ALL dt = {dev:.3e} of scale {scale:.3e}"
      check dev < 1e-11*scale
      # solver path from the 5-fold vertex, checked against its dense analogue
      let (pp, pf) = scalarCorrPoint(o, u, 0.0, v5, 0)
      var devp = 0.0
      for dt in 0..<latF.nt: devp = max(devp, abs(pp[dt] - pf[dt]))
      # dense point-restricted reference for ps
      let s = denseS(o, u, 0.0)
      var worst = 0.0
      for dt in 0..<latF.nt:
        var z = complex64(0.0, 0.0)
        for x in 0..<sph1.nv:
          let i = sIdx(latF, x, dt)
          for c in 0..1:
            for cp in 0..1:
              let
                i2 = 2*i + c
                i1 = 2*sIdx(latF, v5, 0) + cp
              z += s[i2 + ndF*i1]*s[i1 + ndF*i2]
        worst = max(worst, abs(pp[dt] - (-2.0*z.re)))
      echo &"  {nm} point: max|PS-FS| (dt!=0) = {devp:.3e}; vs dense {worst:.3e}"
      check devp < 1e-9*scale
      check worst < 1e-9*scale
      check o.stats.ok

  test "finite-mass FS uses the standard-overlap contact identity":
    let mass = 0.13
    for (nm, u, o) in [("free", uF0, ovF0), ("random", uFr, ovFr)]:
      let
        d0 = denseOv(o, u)
        s = denseS(o, u, mass)
        alpha = ovMassAlpha(mass)
        beta = ovMassBeta(mass)
      # Matrix-level identity, independent of the correlator implementation:
      # (1-D0^dag)S^dag = (beta*S^dag-1)/alpha.
      var worstId = 0.0
      for j in 0..<ndF:
        for i in 0..<ndF:
          var lhs = complex64(0.0, 0.0)
          for k in 0..<ndF:
            var a = -conjugate(d0[k + ndF*i])
            if i == k: a += complex64(1.0, 0.0)
            lhs += a*conjugate(s[j + ndF*k])
          var rhs = beta*conjugate(s[j + ndF*i])
          if i == j: rhs -= complex64(1.0, 0.0)
          rhs = rhs/alpha
          worstId = max(worstId, abs(lhs - rhs))
      check worstId < 1e-11

      # The iterative point contraction must match the dense contraction made
      # from that identity at the same nonzero mass.
      let (pp, pf) = scalarCorrPoint(o, u, mass, v5, 0)
      var worstPs = 0.0
      var worstFs = 0.0
      var scale = 0.0
      for dt in 0..<latF.nt:
        var z1 = complex64(0.0, 0.0)
        var z2 = complex64(0.0, 0.0)
        for x in 0..<sph1.nv:
          let site = sIdx(latF, x, dt)
          for c in 0..1:
            for cp in 0..1:
              let
                i2 = 2*site + c
                i1 = 2*sIdx(latF, v5, 0) + cp
                sxy = s[i2 + ndF*i1]
                syx = s[i1 + ndF*i2]
              z1 += sxy*syx
              var
                a = beta*conjugate(syx)
                b = beta*conjugate(sxy)
              if i2 == i1:
                a -= complex64(1.0, 0.0)
                b -= complex64(1.0, 0.0)
              z2 += (a/alpha)*(b/alpha)
        let
          psRef = -2.0*z1.re
          fsRef = -z1.re - z2.re
        scale = max(scale, max(abs(psRef), abs(fsRef)))
        worstPs = max(worstPs, abs(pp[dt] - psRef))
        worstFs = max(worstFs, abs(pf[dt] - fsRef))
      echo &"  {nm} mass={mass}: FS identity {worstId:.3e}, " &
           &"point PS/FS vs dense {worstPs:.3e}/{worstFs:.3e}"
      check worstPs < 1e-9*scale
      check worstFs < 1e-9*scale
      check o.stats.ok

# ==============================================================================
# gluonic fixture: L = 1, nt = 60, at = 0.2, g2 R = 1, exact-area convention

let
  latG = newLat(sph1, 60, 0.2)
  btG = newBeta(latG, 1.0, gcExactArea)

suite "gluonic exact correlator: icosahedral protection (deterministic)":

  test "l = 1 and l = 2 correlator matrices are proportional to identity":
    for lh in [1, 2]:
      let
        nm = 2*lh + 1
        c = jtopCorrExact(latG, btG, lh)
      var worst = 0.0
      for dt in 1..4:
        let (off, spr, mean) = offDiagId(c[dt], nm)
        worst = max(worst, max(off, spr)/abs(mean))
      echo &"  l={lh}: max (|off|, diag spread)/diag over dt=1..4 = {worst:.3e}"
      check worst < 1e-9

  test "l = 3 splits into 3-fold + 4-fold blocks; the splitting size":
    for (lev, lat) in [(1, latG), (2, newLat(newSphere(2), 60, 0.2))]:
      let
        bt = newBeta(lat, 1.0, gcExactArea)
        c = jtopCorrExact(lat, bt, 3)
        ndt = 12
      var
        lo = newSeq[float](ndt)
        hi = newSeq[float](ndt)
        nlo0 = 0
      for dt in 1..<ndt:
        let
          ev = symEig(c[dt], 7)
          cl = clusters(ev)
        var e = ev
        e.sort
        var nlo: int
        if cl.k == 3:
          lo[dt] = (e[0] + e[1] + e[2])/3.0
          hi[dt] = (e[3] + e[4] + e[5] + e[6])/4.0
          nlo = 3
        else:
          lo[dt] = (e[0] + e[1] + e[2] + e[3])/4.0
          hi[dt] = (e[4] + e[5] + e[6])/3.0
          nlo = 4
        if dt == 1: nlo0 = nlo
        check cl.k == 3 or cl.k == 4
        check nlo == nlo0                       # no cluster crossing in dt
        let wm = (float(nlo)*lo[dt] + float(7-nlo)*hi[dt])/7.0
        check max(cl.s1, cl.s2) < 1e-8*abs(wm)
        if dt <= 2:
          echo &"  L={lev} dt={dt}: clusters {nlo}+{7-nlo}, values {lo[dt]:.6e} / {hi[dt]:.6e}, " &
               &"amplitude split {100.0*(hi[dt]-lo[dt])/wm:.2f} %"
      # the deck-comparable number: the splitting of the effective DIMENSION
      for dt in [2, 5, 8]:
        let
          dLo = -ln(lo[dt+1]/lo[dt])/lat.at
          dHi = -ln(hi[dt+1]/hi[dt])/lat.at
          dm = (float(nlo0)*dLo + float(7-nlo0)*dHi)/7.0
        echo &"  L={lev} t={float(dt)*lat.at:.1f}: Delta({nlo0}-fold) {dLo:.6f}  " &
             &"Delta({7-nlo0}-fold) {dHi:.6f}  mean {dm:.6f}  " &
             &"split {100.0*abs(dHi-dLo)/dm:.2f} %  (free tower sqrt(12) = {sqrt(12.0):.6f})"

suite "single heatbath configuration: group-averaged degeneracy (exact)":

  test "l <= 2 identity and l = 3 clusters at machine precision":
    let grp = icosaGroup(sph1)
    echo &"  group: 60 rotations, siteDev {grp.siteDev:.2e}  faceDev {grp.faceDev:.2e}"
    check grp.faceDev < 1e-9
    var r: Threefry4x64
    r.seedIndep(606, 0)
    var u = newGauge(latG)
    let ci = heatbath(latG, u, btG, r)
    check ci.converged
    # fluxes once
    var th = newSeq[seq[float]](sph1.nf)
    for f in 0..<sph1.nf:
      th[f] = newSeq[float](latG.nt)
      for t in 0..<latG.nt: th[f][t] = plaqSpatial(latG, u, f, t)
    for lh in 1..3:
      let nm = 2*lh + 1
      # Y matrix per face
      var yf = newSeq[seq[float]](nm)
      for m in 0..<nm:
        yf[m] = newSeq[float](sph1.nf)
        for f in 0..<sph1.nf: yf[m][f] = ylm(lh, m - lh, sph1.faces[f].cc)
      var cavg = newSeq[float](nm*nm)
      for g in 0..<60:
        var o = newSeq[seq[float]](nm)
        for m in 0..<nm:
          o[m] = newSeq[float](latG.nt)
          for f in 0..<sph1.nf:
            let fp = grp.face[g][f]
            for t in 0..<latG.nt: o[m][t] += yf[m][f]*th[fp][t]
        for m in 0..<nm:
          for mp in 0..<nm:
            var s = 0.0
            for t in 0..<latG.nt: s += o[m][t]*o[mp][t]
            cavg[m*nm + mp] += s/60.0
      if lh <= 2:
        let (off, spr, mean) = offDiagId(cavg, nm)
        echo &"  l={lh}: diag {mean:.6e}  (|off| {off:.2e}, spread {spr:.2e})" &
             &"  rel {max(off, spr)/abs(mean):.3e}"
        check max(off, spr) < 3e-12*abs(mean)
      else:
        let
          ev = symEig(cavg, nm)
          cl = clusters(ev)
        var e = ev
        e.sort
        echo &"  l=3: clusters {cl.k}+{7-cl.k}, spreads {cl.s1:.2e}/{cl.s2:.2e}, gap {cl.gap:.3e}"
        check cl.k == 3 or cl.k == 4
        check max(cl.s1, cl.s2) < 1e-10*max(abs(e[0]), abs(e[6]))

suite "gluonic Monte Carlo pipeline at L=1, nt=60, at=0.2, g2R=1":

  test "J_top l=1: MC vs exact correlator, effective dimension vs sqrt(2)":
    const
      ncfg = 128
      maxDt = 20
    let
      nt = latG.nt
      at = latG.at
      tt = at*float(nt)
    # exact reference (m-averaged diagonal)
    let cex0 = jtopCorrExact(latG, btG, 1)
    var cex = newSeq[float](nt)
    for dt in 0..<nt:
      for m in 0..2: cex[dt] += cex0[dt][m*3 + m]/3.0
    # ensemble
    var r: Threefry4x64
    r.seedIndep(707, 0)
    var u = newGauge(latG)
    var cmc = newSeq[seq[float]](ncfg)     # per-config correlator
    for k in 0..<ncfg:
      let ci = heatbath(latG, u, btG, r)
      doAssert ci.converged
      var o = newSeq[seq[float]](3)
      for m in 0..2: o[m] = jtopProject(latG, u, 1, m - 1)
      cmc[k] = newSeq[float](nt)
      for dt in 0..<nt:
        var s = 0.0
        for m in 0..2:
          for t1 in 0..<nt: s += o[m][(t1 + dt) mod nt]*o[m][t1]
        cmc[k][dt] = s/float(3*nt)
    # pulls MC vs exact
    var
      cm = newSeq[float](nt)
      ce = newSeq[float](nt)
    var worstPull = 0.0
    for dt in 0..<nt:
      var s = 0.0
      for k in 0..<ncfg: s += cmc[k][dt]
      cm[dt] = s/float(ncfg)
      var s2 = 0.0
      for k in 0..<ncfg: s2 += (cmc[k][dt] - cm[dt])^2
      ce[dt] = sqrt(s2/float(ncfg*(ncfg - 1)))
      if dt <= maxDt:
        worstPull = max(worstPull, abs(cm[dt] - cex[dt])/max(ce[dt], 1e-300))
    echo &"  C(0): exact {cex[0]:.6e}  MC {cm[0]:.6e} +- {ce[0]:.2e}"
    echo &"  C(5): exact {cex[5]:.6e}  MC {cm[5]:.6e} +- {ce[5]:.2e}"
    echo &"  worst |MC - exact|/sigma over dt <= {maxDt}: {worstPull:.2f}  (ncfg {ncfg})"
    check worstPull < 5.0
    # effective dimension: exact correlator (deterministic physics number)
    let dex = effMass(cex, at, tt)
    echo &"  exact Delta_eff at t = 1, 2, 3, 4: {dex[5]:.6f} {dex[10]:.6f} " &
         &"{dex[15]:.6f} {dex[20]:.6f}"
    let pf = plateauFit(dex, 15, 30, at)
    echo &"  exact plateau fit over t in [3,6): Delta_0 = {pf.d0:.6f} +- {pf.ed0:.6f} " &
         &"(sqrt(2) = {sqrt(2.0):.6f}, published L=1 exact-area value 1.33242 at T=16)"
    # MC effective dimension with delete-1 jackknife (t = 1 and t = 2)
    for tIdx in [5, 10]:
      let dmc = effMass(cm, at, tt)[tIdx]
      var js = newSeq[float](ncfg)
      for k in 0..<ncfg:
        var cj = newSeq[float](nt)
        for dt in 0..<nt: cj[dt] = (cm[dt]*float(ncfg) - cmc[k][dt])/float(ncfg - 1)
        js[k] = effMass(cj, at, tt)[tIdx]
      var jm = 0.0
      for x in js: jm += x/float(ncfg)
      var jv = 0.0
      for x in js: jv += (x - jm)^2
      let de = sqrt(jv*float(ncfg - 1)/float(ncfg))
      let dexact = dex[tIdx]
      echo &"  Delta_eff(t={float(tIdx)*at:.0f}): MC {dmc:.4f} +- {de:.4f}   exact {dexact:.6f}   " &
           &"pull {(dmc - dexact)/de:.2f};  vs sqrt(2): {(dmc - sqrt(2.0))/de:.2f} sigma"
      check abs(dmc - dexact) < 5.0*de

suite "7-shape Wilson-loop basis and GEVP":

  test "all shapes are exactly gauge invariant":
    let u = randGauge(latG, 808, 0.7)
    var ug = u
    var alpha = newSeq[float](latG.nsite)
    var r: Threefry4x64
    r.seedIndep(809, 0)
    for i in 0..<alpha.len: alpha[i] = r.gaussian
    gaugeTransform(latG, ug, alpha)
    var worst = 0.0
    for sh in [lsTri, lsRhomb, lsStar, lsQuad, lsTPlaq, lsTRect2, lsTRhomb]:
      let
        a = loopProject(latG, u, sh, 1, 0)
        b = loopProject(latG, ug, sh, 1, 0)
      var scale = 1e-300
      for t in 0..<latG.nt: scale = max(scale, abs(a[t]))
      for t in 0..<latG.nt: worst = max(worst, abs(a[t] - b[t])/scale)
    echo &"  max relative change under a random gauge transformation: {worst:.3e}"
    check worst < 1e-12

  test "L=1 collapse: every spatial shape is ONE operator after l-projection":
    ## Measured finding (doc/06 WP-I): T1 has multiplicity 1 in the 20-face
    ## permutation rep of the icosahedron, so at L = 1 the l = 1 projections
    ## of ALL spatial shapes (tri, rhombus, star, quad) are proportional to
    ## each other configuration by configuration, and the 6-shape GEVP matrix
    ## has rank 3.  Demonstrated here on one heatbath configuration.
    var r: Threefry4x64
    r.seedIndep(818, 0)
    var u = newGauge(latG)
    let ci = heatbath(latG, u, btG, r)
    check ci.converged
    let a = loopProject(latG, u, lsTri, 1, 0)
    var n2a = 0.0
    for t in 0..<latG.nt: n2a += a[t]*a[t]
    for sh in [lsRhomb, lsStar, lsQuad]:
      let b = loopProject(latG, u, sh, 1, 0)
      var ab = 0.0
      var n2b = 0.0
      for t in 0..<latG.nt:
        ab += a[t]*b[t]
        n2b += b[t]*b[t]
      let c = ab/n2a
      var res = 0.0
      for t in 0..<latG.nt: res += (b[t] - c*a[t])^2
      echo &"  {sh} = {c:.12f} * tri  (relative residual {sqrt(res/n2b):.3e})"
      check sqrt(res/n2b) < 1e-12

  test "exact loop GEVP ground state == l=1-projected correlator dimension":
    ## L = 1: the independent basis is {tri, tplaq, trect2} (see the collapse
    ## test); at L = 2 all 7 shapes including the quadruple triangle enter.
    let cex0 = jtopCorrExact(latG, btG, 1)
    var cexL1 = newSeq[float](latG.nt)
    for dt in 0..<latG.nt:
      for m in 0..2: cexL1[dt] += cex0[dt][m*3 + m]/3.0
    let latG2 = newLat(newSphere(2), 60, 0.2)
    let btG2 = newBeta(latG2, 1.0, gcExactArea)
    let cex2 = jtopCorrExact(latG2, btG2, 1)
    var cexL2 = newSeq[float](latG2.nt)
    for dt in 0..<latG2.nt:
      for m in 0..2: cexL2[dt] += cex2[dt][m*3 + m]/3.0
    ## t0 = 5: the extent-2 temporal shapes span t +- 2 slices, so C(t0) is a
    ## genuine (positive) transfer-matrix correlator only once the sink and
    ## source supports are disjoint, t0 > 4.  (At t0 = 2 the measured C(t0)
    ## is INDEFINITE, evmin ~ -0.86 -- overlapping-support contact terms.)
    ## The GEVP is compared to the correlator at t <= 2: beyond that the
    ## periodic image (T = 12) shifts the plain-ratio GEVP dimension at the
    ## 1e-3..1e-2 level while arccosh-based effMass folds it out; reported.
    for (lev, lat, bt, cex, shapes) in [
        (1, latG, btG, cexL1, @[lsTri, lsTPlaq, lsTRect2]),
        (2, latG2, btG2, cexL2,
         @[lsTri, lsRhomb, lsStar, lsQuad, lsTPlaq, lsTRect2, lsTRhomb])]:
      let c = loopCorrExact(lat, bt, shapes, 1, 0)
      let hc = gevpCheck(c, 5)
      var cs = newSeq[seq[seq[float]]](26)
      for dt in 0..<cs.len: cs[dt] = c[dt]
      let
        dims = gevpDims(cs, 5, lat.at)
        dref = effMass(cex, lat.at, lat.at*float(lat.nt))
      echo &"  L={lev} C(t0=5): evmin {hc.evmin:.3e}  evmax {hc.evmax:.3e}  " &
           &"asym {hc.asym:.2e}  rank kept {dims[6].len} of {shapes.len} " &
           &"(the channel holds that many states)"
      var line = &"  L={lev} GEVP Delta_0(t):  "
      for t in [6, 8, 10, 15, 20]: line &= &"{dims[t][0]:.6f} "
      echo line
      echo &"  L={lev} corr Delta_eff(t): {dref[6]:.6f} {dref[8]:.6f} " &
           &"{dref[10]:.6f} {dref[15]:.6f} {dref[20]:.6f}"
      var worst = 0.0
      for t in [6, 8, 10]: worst = max(worst, abs(dims[t][0] - dref[t]))
      echo &"  L={lev} max |GEVP - correlator| over t in [1.2, 2.0]: {worst:.2e}" &
           &"  (image drift at t=3, 4: {abs(dims[15][0]-dref[15]):.1e}, " &
           &"{abs(dims[20][0]-dref[20]):.1e})"
      check hc.asym < 1e-8
      check hc.evmin > -1e-10*hc.evmax     # PSD up to roundoff; rank < nshapes is physics
      check worst < 1e-3

  test "flowed-config MC loop GEVP vs the exact flowed GEVP":
    const
      ncfg = 96
      nkeep = 16          # correlator matrices kept for dt < nkeep
    let
      shapes = [lsTri, lsTPlaq, lsTRect2]   # the independent l=1 basis at L=1
      nsh = shapes.len
      s = 0.2             # MC flow time (deck: t in [0.2, 1.6])
      h = 0.1/mDiagMax(latG, btG)
    # exact flowed correlator: c_i^T e^{-sM} M^+ e^{-sM} c_j, at two flow times
    var cflow: seq[seq[seq[float]]]         # kept for the MC comparison (s = 0.2)
    for sf in [0.2, 0.6]:
      var cf = newSeq[seq[seq[float]]](nkeep)
      for dt in 0..<nkeep:
        cf[dt] = newSeq[seq[float]](nsh)
        for i in 0..<nsh: cf[dt][i] = newSeq[float](nsh)
      var op = newRegOp(latG, btG)
      var b = newGauge(latG)
      var x = newGauge(latG)
      for j in 0..<nsh:
        loopSource(latG, b, shapes[j], 1, 0, 0)
        flowRun(latG, b, btG, [sf], h, RK4CK_2N, nil)
        let ci = regSolve(latG, x, b, op, 1e-26, 200000)
        doAssert ci.converged
        flowRun(latG, x, btG, [sf], h, RK4CK_2N, nil)
        for dt in 0..<nkeep:
          for i in 0..<nsh:
            var acc = 0.0
            let nc = loopCount(latG, shapes[i])
            for c in 0..<nc:
              acc += ylm(1, 0, loopCenter(sph1, shapes[i], c))*
                     loopFlux(latG, x, shapes[i], c, dt)
            cf[dt][i][j] = acc
      let dt0 = gevpDims(cf, 5, latG.at, 0.05)
      echo &"  exact flowed GEVP at s = {sf}: Delta_0(t=1.2, 1.6, 2.2) = " &
           &"{dt0[6][0]:.6f} {dt0[8][0]:.6f} {dt0[11][0]:.6f}"
      if sf == s: cflow = cf
    # MC on flowed heatbath configurations, per-config matrices for jackknife
    var r: Threefry4x64
    r.seedIndep(909, 0)
    var u = newGauge(latG)
    var cmcK = newSeq[seq[seq[seq[float]]]](ncfg)
    for k in 0..<ncfg:
      let ci = heatbath(latG, u, btG, r)
      doAssert ci.converged
      var uf = u
      flowRun(latG, uf, btG, [s], h, RK4CK_2N, nil)
      var o = newSeq[seq[float]](nsh)
      for i in 0..<nsh: o[i] = loopProject(latG, uf, shapes[i], 1, 0)
      cmcK[k] = newSeq[seq[seq[float]]](nkeep)
      for dt in 0..<nkeep:
        cmcK[k][dt] = newSeq[seq[float]](nsh)
        for i in 0..<nsh:
          cmcK[k][dt][i] = newSeq[float](nsh)
          for j in 0..<nsh:
            var acc = 0.0
            for t1 in 0..<latG.nt:
              acc += o[i][(t1 + dt) mod latG.nt]*o[j][t1]
            cmcK[k][dt][i][j] = acc/float(latG.nt)
    proc mcMean(skip: int): seq[seq[seq[float]]] =
      result = newSeq[seq[seq[float]]](nkeep)
      let w = 1.0/float(if skip < 0: ncfg else: ncfg - 1)
      for dt in 0..<nkeep:
        result[dt] = newSeq[seq[float]](nsh)
        for i in 0..<nsh:
          result[dt][i] = newSeq[float](nsh)
          for j in 0..<nsh:
            var acc = 0.0
            for k in 0..<ncfg:
              if k != skip: acc += cmcK[k][dt][i][j]
            result[dt][i][j] = acc*w
    ## cut = 0.05: the L=1 l=1 channel is rank 1 (see the collapse test), so
    ## the 2nd/3rd C(t0) directions are exactly null for the exact correlator
    ## and pure noise for the MC one; both must be truncated at the MC noise
    ## level for the two pipelines to solve the same eigenproblem.
    const tCmp = 8   # t = 1.6, inside the deck's flow-time window
    let
      dEx = gevpDims(cflow, 5, latG.at, 0.05)
      dMc = gevpDims(mcMean(-1), 5, latG.at, 0.05)
    var js = newSeq[float](ncfg)
    for k in 0..<ncfg: js[k] = gevpDims(mcMean(k), 5, latG.at, 0.05)[tCmp][0]
    var jm = 0.0
    for x in js: jm += x/float(ncfg)
    var jv = 0.0
    for x in js: jv += (x - jm)^2
    let de = sqrt(jv*float(ncfg - 1)/float(ncfg))
    echo &"  ranks kept: exact {dEx[6].len}  MC {dMc[6].len}  of {nsh}"
    echo &"  flow s = {s}: exact flowed GEVP Delta_0(t): " &
         &"{dEx[6][0]:.6f} {dEx[8][0]:.6f} {dEx[11][0]:.6f}  " &
         &"(unflowed exact 1.329856: flow leaves the energy alone)"
    echo &"                MC    flowed GEVP Delta_0(t): " &
         &"{dMc[6][0]:.6f} {dMc[8][0]:.6f} {dMc[11][0]:.6f}   (ncfg {ncfg})"
    let pull = (dMc[tCmp][0] - dEx[tCmp][0])/de
    echo &"  at t = 1.6: MC {dMc[tCmp][0]:.4f} +- {de:.4f} (jackknife)  " &
         &"exact {dEx[tCmp][0]:.6f}  pull {pull:.2f}"
    check abs(pull) < 5.0
