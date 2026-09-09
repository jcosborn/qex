## the tree-level clover term of hcwilson (cSW) and of
## cubic.CubicWilson.
##
## 1. cSW = 0 reduces to the bare operator (norm2diff = 0);
## 2. gauge covariance of D_c to 1e-12;
## 3. gamma5-hermiticity of D_c; the added term C = D_c - D is HERMITIAN and
##    gamma5-even (both facts are what gamma5-hermiticity of D_c requires --
##    the task brief's "anti-Hermitian" is a slip, checked numerically here);
## 4. on-site locality: support(D_c delta) = support(D delta), and
##    (D_c - D) delta lives on the source site only;
## 5. THE ACCEPTANCE TEST -- constant-flux normalisation pin, both lattices:
##    on the exact Atiyah-Singer background (F_01 = f1, F_23 = f2 via
##    T = diag(1,-1,0)) the measured on-site 12x12 matrix (D_c - D)(x) must be
##      (i/2) [ f1 s1 gamma2gamma1 + f2 s2 gamma4gamma3 ] (x) T,
##    i.e. the standard tree-level -(cSW r/4) sigma_munu F_munu at cSW = 1,
##    times the exactly known clover artifact factors
##      s_hc = 4 sin(f/4)/f      (hexagon clover)
##      s_cubic = sin(f)/f       (1x1 clover)
##    -> measured/continuum-analytic reported per plane (= s, the artifact,
##    1 - s = O(1/L^4), scaling shown), measured/lattice-analytic = 1 to
##    ~1e-12, structure residual ~1e-13, and the chirality-splitting
##    cross-check: for f1 = f2 ("self-dual") the term vanishes on the
##    gamma5 = +1 spin sector, for f1 = -f2 on gamma5 = -1, identically on
##    BOTH lattices (a relative sign error between the two planes, an overall
##    sign error, or an hc/cubic convention mismatch would each break this).
## The F = 0 case (D_c == D on the unit gauge) is covered by tspectrum 1a/1b.

import std/[math, strformat, unittest, complex]
import qex except epsilon
import physics/qcdTypes
import gauge, gauge/gaugefix
import ../hcgeom
import ../hcgauge
import ../hcwilson
import ../hctopo
import ../cubic
import ../hcfree
import helpers

qexInit()

let seed = 1928374655'u64

proc spinColor(f: auto, i: int): array[12, Complex64] =
  ## the 12 spin-color components at site i, index s*3 + c
  for s in 0..<4:
    for c in 0..<nc:
      result[s*3 + c] = complex64(toF f{i}[s][c].re, toF f{i}[s][c].im)

# 12x12 complex spin (x) color matrices

type CMat12 = array[12, array[12, Complex64]]

proc kronz(z: Complex64, s: hcfree.Mat4,
           t: array[3, array[3, Complex64]]): CMat12 =
  for s1 in 0..<4:
    for c1 in 0..<3:
      for s2 in 0..<4:
        for c2 in 0..<3:
          result[s1*3+c1][s2*3+c2] = z*s[s1][s2]*t[c1][c2]

proc tMat(): array[3, array[3, Complex64]] =
  result[0][0] = complex64(1.0)
  result[1][1] = complex64(-1.0)

proc iMat3(): array[3, array[3, Complex64]] =
  for c in 0..<3: result[c][c] = complex64(1.0)

proc frob(a, b: CMat12): Complex64 =
  for i in 0..<12:
    for j in 0..<12:
      result += conjugate(a[i][j])*b[i][j]

proc fnorm(a: CMat12): float =
  sqrt(frob(a, a).re)

proc mmul12(a, b: CMat12): CMat12 =
  for i in 0..<12:
    for j in 0..<12:
      var z: Complex64
      for k in 0..<12:
        z += a[i][k]*b[k][j]
      result[i][j] = z

proc maxDiff12(a, b: CMat12): float =
  for i in 0..<12:
    for j in 0..<12:
      result = max(result, abs(a[i][j] - b[i][j]))

# chirality projectors P+- = (1 +- gamma5)/2 (x) 1_color
proc chirProj(sign: float): CMat12 =
  let g5 = kronz(complex64(1.0), gamma5Mat, iMat3())
  for i in 0..<12:
    for j in 0..<12:
      result[i][j] = 0.5*(complex64(float(i == j)) + complex64(sign)*g5[i][j])

# ---------------------------------------------------------------------------
# on-site clover matrix extraction: 12 point sources, C = (D_c - D) delta
# ---------------------------------------------------------------------------

proc measureCloverHc(cop, cop0: auto, hl: auto, i0: int, sub: int,
                     m: float): CMat12 =
  var x = newHcFermion(hl)
  var r1 = newHcFermion(hl)
  var r0 = newHcFermion(hl)
  for j in 0..<12:
    let s = j div 3
    let cc = j mod 3
    threads:
      x := 0
    if sub == 0: x.a{i0}[s][cc].re := 1.0
    else: x.b{i0}[s][cc].re := 1.0
    cop.D(r1, x, m)             # improved
    cop0.D(r0, x, m)          # bare
    let v1 = spinColor(if sub == 0: r1.a else: r1.b, i0)
    let v0 = spinColor(if sub == 0: r0.a else: r0.b, i0)
    for i in 0..<12:
      result[i][j] = v1[i] - v0[i]

proc measureCloverCubic(cop: auto, lo: auto, i0: int, m: float): CMat12 =
  var x = lo.DiracFermion()
  var r1 = lo.DiracFermion()
  var r0 = lo.DiracFermion()
  for j in 0..<12:
    let s = j div 3
    let cc = j mod 3
    threads:
      x := 0
    x{i0}[s][cc].re := 1.0
    cop.D(r1, x, m)             # improved (opens threads itself)
    cop.s.D(r0, x, m)           # bare QEX Wilson (serial call is fine)
    let v1 = spinColor(r1, i0)
    let v0 = spinColor(r0, i0)
    for i in 0..<12:
      result[i][j] = v1[i] - v0[i]

type CloverPin = object
  s1, s2: float                 ## measured/continuum-analytic per plane
  resid: float                  ## structure residual after the 2-basis fit
  npos, nneg, ncross: float     ## chirality block norms of C
  cnorm: float

proc analyzeClover(cm: CMat12, f1, f2: float): CloverPin =
  ## project the measured C onto the two analytic basis matrices
  ##   B1 = i gamma2gamma1 (x) T,  B2 = i gamma4gamma3 (x) T
  ## (analytic target: C = (f1 s1/2) B1 + (f2 s2/2) B2)
  let b1 = kronz(complex64(0.0, 1.0), mmul(gammaMat[1], gammaMat[0]), tMat())
  let b2 = kronz(complex64(0.0, 1.0), mmul(gammaMat[3], gammaMat[2]), tMat())
  let n1 = frob(b1, b1).re
  let n2 = frob(b2, b2).re
  let c1 = frob(b1, cm)/n1      # should be real (f1 s1/2)
  let c2 = frob(b2, cm)/n2
  result.s1 = 2.0*c1.re/f1
  result.s2 = 2.0*c2.re/f2
  var res = cm
  for i in 0..<12:
    for j in 0..<12:
      res[i][j] -= c1*b1[i][j] + c2*b2[i][j]
  result.cnorm = fnorm(cm)
  result.resid = fnorm(res)/result.cnorm
  let pp = chirProj(1.0)
  let pm = chirProj(-1.0)
  result.npos = fnorm(mmul12(pp, mmul12(cm, pp)))
  result.nneg = fnorm(mmul12(pm, mmul12(cm, pm)))
  result.ncross = fnorm(mmul12(pp, mmul12(cm, pm))) +
                  fnorm(mmul12(pm, mmul12(cm, pp)))

# ---------------------------------------------------------------------------

suite "hcclover (task D2 part 2)":

  test "1. cSW = 0 reduces to the bare operator":
    block:                      # honeycomb
      let hl = newLayout([4, 4, 4, 6])
      var r = hl.newRNGField(RngMilc6, seed)
      var g = newHcGauge(hl)
      threads:
        g.random r
      var c0 = newHcWilson(g, 0.0)
      var w = newHcWilson(g)    # independently built bare operator
      var x = newHcFermion(hl)
      var y1 = newHcFermion(hl)
      var y0 = newHcFermion(hl)
      threads:
        x.gaussian r
      c0.D(y1, x, 0.1)
      w.D(y0, x, 0.1)
      let d = norm2diff(y1, y0)
      echo &"  honeycomb |D_c(cSW=0) x - D x|^2 = {d:.3e}"
      ok("honeycomb cSW=0 == bare (norm2diff = 0)", d == 0.0)
      c0.Ddag(y1, x, 0.1)
      w.Ddag(y0, x, 0.1)
      ok("honeycomb cSW=0 Ddag == bare Ddag", norm2diff(y1, y0) == 0.0)
    block:                      # cubic
      let lo = newLayout(@[4, 4, 4, 8])
      var r = lo.newRNGField(RngMilc6, seed + 1)
      var g = lo.newGauge
      threads:
        g.random r
      var c0 = newCubicWilson(g, 0.0)
      var x = lo.DiracFermion()
      var y1 = lo.DiracFermion()
      var y0 = lo.DiracFermion()
      threads:
        x.gaussian r
      c0.D(y1, x, 0.1)
      c0.s.D(y0, x, 0.1)
      let d = norm2diff(y1, y0)
      echo &"  cubic |D_c(cSW=0) x - D x|^2 = {d:.3e}"
      ok("cubic cSW=0 == bare (norm2diff = 0)", d == 0.0)

  test "2. gauge covariance of D_c":
    block:                      # honeycomb
      let hl = newLayout([4, 4, 4, 6])
      let lo = hl
      var rs = lo.newRNGField(RngMilc6, seed + 2)
      var g = newHcGauge(hl)
      var g2 = newHcGauge(hl)
      threads:
        g.random rs
        g2 := g
      var c = newHcWilson(g, 1.0)
      var x = newHcFermion(hl)
      var xg = newHcFermion(hl)
      var y1 = newHcFermion(hl)
      var y2 = newHcFermion(hl)
      var yv = newHcFermion(hl)
      threads:
        x.gaussian rs
      c.D(y1, x, 0.1)
      var vA = lo.ColorMatrix(nc)
      var vB = lo.ColorMatrix(nc)
      threads:
        vA.randomSU rs
        vB.randomSU rs
      g2.gaugeTransform(vA, vB)
      var c2 = newHcWilson(g2, 1.0)
      threads:
        for e in xg.a:
          xg.a[e] := vA[e] * x.a[e]
        for e in xg.b:
          xg.b[e] := vB[e] * x.b[e]
      c2.D(y2, xg, 0.1)
      threads:
        for e in yv.a:
          yv.a[e] := vA[e] * y1.a[e]
        for e in yv.b:
          yv.b[e] := vB[e] * y1.b[e]
      let rel = sqrt(norm2diff(y2, yv)/norm2(yv))
      echo &"  honeycomb |D_c[U^g] V x - V D_c[U] x| / |.| = {rel:.3e}"
      ok(&"honeycomb D_c gauge covariant to 1e-12 ({rel:.3e})", rel < 1e-12)
    block:                      # cubic
      let lo = newLayout(@[4, 4, 4, 8])
      var rs = lo.newRNGField(RngMilc6, seed + 3)
      var g = lo.newGauge
      var g2 = lo.newGauge
      var t = lo.ColorMatrix(nc)
      threads:
        g.random rs
        t.randomSU rs
      var c = newCubicWilson(g, 1.0)
      var x = lo.DiracFermion()
      var xg = lo.DiracFermion()
      var y1 = lo.DiracFermion()
      var y2 = lo.DiracFermion()
      var yv = lo.DiracFermion()
      threads:
        x.gaussian rs
      c.D(y1, x, 0.1)
      gaugeTransform(g2, g, t)
      var c2 = newCubicWilson(g2, 1.0)
      threads:
        for e in xg:
          xg[e] := t[e] * x[e]
      c2.D(y2, xg, 0.1)
      threads:
        for e in yv:
          yv[e] := t[e] * y1[e]
      let rel = sqrt(norm2diff(y2, yv)/norm2(yv))
      echo &"  cubic |D_c[U^g] V x - V D_c[U] x| / |.| = {rel:.3e}"
      ok(&"cubic D_c gauge covariant to 1e-12 ({rel:.3e})", rel < 1e-12)

  test "3. gamma5-hermiticity; the added term is Hermitian and gamma5-even":
    block:                      # honeycomb
      let hl = newLayout([4, 4, 4, 4])
      var rs = hl.newRNGField(RngMilc6, seed + 4)
      var g = newHcGauge(hl)
      threads:
        g.random rs
      var c = newHcWilson(g, 1.0)
      var x = newHcFermion(hl)
      var y = newHcFermion(hl)
      var dy = newHcFermion(hl)
      var dx = newHcFermion(hl)
      var t1 = newHcFermion(hl)
      var t2 = newHcFermion(hl)
      threads:
        x.gaussian rs
        y.gaussian rs
      let m = 0.15
      # adjointness <x, D_c y> = <D_c^dag x, y>
      c.D(dy, y, m)
      c.Ddag(dx, x, m)
      let z1 = cdot(x, dy)
      let z2 = cdot(dx, y)
      let relAdj = abs(z1 - z2)/abs(z1)
      ok(&"honeycomb Ddag is the adjoint of D_c ({relAdj:.3e})", relAdj < 1e-12)
      # gamma5-hermiticity <x, g5 D_c g5 y> = <D_c x, y>
      threads:
        t1.applyGamma5 y
      c.D(t2, t1, m)
      threads:
        t1.applyGamma5 t2
      c.D(dy, x, m)
      let z3 = cdot(x, t1)
      let z4 = cdot(dy, y)
      let relG5 = abs(z3 - z4)/abs(z3)
      echo &"  honeycomb <x, g5 D_c g5 y> = {z3.re:.10f}{z3.im:+.10f} i"
      echo &"            <D_c x, y>       = {z4.re:.10f}{z4.im:+.10f} i"
      ok(&"honeycomb <x, g5 D_c g5 y> = <D_c x, y> ({relG5:.3e})", relG5 < 1e-12)
      # the ADDED term C = D_c - D: Hermitian, gamma5-even
      var w = newHcWilson(g)
      var cy = newHcFermion(hl)
      var cx = newHcFermion(hl)
      c.D(dy, y, m)
      w.D(t2, y, m)
      threads:
        for e in cy.a: cy.a[e] := dy.a[e] - t2.a[e]
        for e in cy.b: cy.b[e] := dy.b[e] - t2.b[e]
      c.D(dx, x, m)
      w.D(t2, x, m)
      threads:
        for e in cx.a: cx.a[e] := dx.a[e] - t2.a[e]
        for e in cx.b: cx.b[e] := dx.b[e] - t2.b[e]
      let h1 = cdot(x, cy)
      let h2 = cdot(cx, y)
      let relH = abs(h1 - h2)/abs(h1)
      let relA = abs(h1 + h2)/abs(h1)
      echo &"  added term: <x, C y> = {h1.re:.10f}{h1.im:+.10f} i"
      echo &"              <C x, y> = {h2.re:.10f}{h2.im:+.10f} i"
      echo &"  Hermitian dev {relH:.2e}; anti-Hermitian dev {relA:.2e} (O(1) as expected)"
      ok(&"added term is HERMITIAN ({relH:.2e}); NOT anti-Hermitian " &
         &"({relA:.2f} ~ 2)", relH < 1e-12 and relA > 0.5)
      # gamma5-evenness: g5 C g5 y == C y
      threads:
        t1.applyGamma5 y
      c.D(t2, t1, m)
      w.D(dx, t1, m)
      threads:
        for e in t2.a: t2.a[e] := t2.a[e] - dx.a[e]
        for e in t2.b: t2.b[e] := t2.b[e] - dx.b[e]
        threadBarrier()
        t1.applyGamma5 t2
      let relE = sqrt(norm2diff(t1, cy)/norm2(cy))
      ok(&"added term is gamma5-even ({relE:.3e})", relE < 1e-12)
    block:                      # cubic
      let lo = newLayout(@[4, 4, 4, 8])
      var rs = lo.newRNGField(RngMilc6, seed + 5)
      var g = lo.newGauge
      threads:
        g.random rs
      var c = newCubicWilson(g, 1.0)
      var x = lo.DiracFermion()
      var y = lo.DiracFermion()
      var dy = lo.DiracFermion()
      var dx = lo.DiracFermion()
      var t1 = lo.DiracFermion()
      var t2 = lo.DiracFermion()
      threads:
        x.gaussian rs
        y.gaussian rs
      let m = 0.15
      c.D(dy, y, m)
      c.Ddag(dx, x, m)
      let z1 = cdot(x, dy)
      let z2 = cdot(dx, y)
      let relAdj = abs(z1 - z2)/abs(z1)
      ok(&"cubic Ddag is the adjoint of D_c ({relAdj:.3e})", relAdj < 1e-12)
      threads:
        for e in t1: t1[e] := gamma5 * y[e]
      c.D(t2, t1, m)
      threads:
        for e in t1: t1[e] := gamma5 * t2[e]
      c.D(dy, x, m)
      let z3 = cdot(x, t1)
      let z4 = cdot(dy, y)
      let relG5 = abs(z3 - z4)/abs(z3)
      ok(&"cubic <x, g5 D_c g5 y> = <D_c x, y> ({relG5:.3e})", relG5 < 1e-12)

  test "4. on-site locality of the added term":
    block:                      # honeycomb
      let geom = @[4, 4, 4, 4]
      let hl = newLayout(geom)
      let lo = hl
      var rs = lo.newRNGField(RngMilc6, seed + 6)
      var g = newHcGauge(hl)
      threads:
        g.random rs
      var c = newHcWilson(g, 1.0)
      var cb = newHcWilson(g, 0.0)
      var x = newHcFermion(hl)
      var r1 = newHcFermion(hl)
      var r0 = newHcFermion(hl)
      let c0 = [0, 3, 0, 3]     # wrap corner (twilson test 4)
      let i0 = lo.rankIndex(c0).index
      threads:
        x := 0
      x.a{i0}[1][0].re := 1.0
      c.D(r1, x, 0.2, 0.9)
      cb.D(r0, x, 0.2, 0.9)
      var nnz1, nnz0 = 0
      var offsite = 0.0
      var onsite = 0.0
      for i in 0..<lo.nSites:
        for sub in 0..1:
          let f1 = (if sub == 0: r1.a else: r1.b)
          let f0 = (if sub == 0: r0.a else: r0.b)
          if siteNorm2(f1, i) > 1e-24: inc nnz1
          if siteNorm2(f0, i) > 1e-24: inc nnz0
          var d = 0.0
          for s in 0..<4:
            for cc in 0..<nc:
              d += abs(toF(f1{i}[s][cc].re) - toF(f0{i}[s][cc].re)) +
                   abs(toF(f1{i}[s][cc].im) - toF(f0{i}[s][cc].im))
          if i == i0 and sub == 0: onsite = d
          else: offsite = max(offsite, d)
      echo &"  honeycomb: support(D_c delta) = {nnz1} sites, support(D delta) = {nnz0}"
      echo &"  |C delta| on-site = {onsite:.3e}, max off-site = {offsite:.3e}"
      ok("honeycomb support unchanged (25 sites)", nnz1 == 25 and nnz0 == 25)
      ok("honeycomb added term strictly on-site", offsite == 0.0 and onsite > 1e-6)
    block:                      # cubic
      let lo = newLayout(@[4, 4, 4, 8])
      var rs = lo.newRNGField(RngMilc6, seed + 7)
      var g = lo.newGauge
      threads:
        g.random rs
      var c = newCubicWilson(g, 1.0)
      var x = lo.DiracFermion()
      var r1 = lo.DiracFermion()
      var r0 = lo.DiracFermion()
      let c0 = [0, 3, 0, 7]
      let i0 = lo.rankIndex(c0).index
      threads:
        x := 0
      x{i0}[1][0].re := 1.0
      c.D(r1, x, 0.2)
      c.s.D(r0, x, 0.2)
      var nnz1, nnz0 = 0
      var offsite = 0.0
      var onsite = 0.0
      for i in 0..<lo.nSites:
        if siteNorm2(r1, i) > 1e-24: inc nnz1
        if siteNorm2(r0, i) > 1e-24: inc nnz0
        var d = 0.0
        for s in 0..<4:
          for cc in 0..<nc:
            d += abs(toF(r1{i}[s][cc].re) - toF(r0{i}[s][cc].re)) +
                 abs(toF(r1{i}[s][cc].im) - toF(r0{i}[s][cc].im))
        if i == i0: onsite = d
        else: offsite = max(offsite, d)
      echo &"  cubic: support(D_c delta) = {nnz1} sites, support(D delta) = {nnz0}"
      echo &"  |C delta| on-site = {onsite:.3e}, max off-site = {offsite:.3e}"
      ok("cubic support unchanged (9 sites)", nnz1 == 9 and nnz0 == 9)
      ok("cubic added term strictly on-site", offsite == 0.0 and onsite > 1e-6)

  test "5. ACCEPTANCE: constant-flux normalisation pin, both lattices":
    # exact targets: C = (i/2)[f1 s1 g2g1 + f2 s2 g4g3] (x) T with
    # s_hc = 4 sin(f/4)/f, s_cubic = sin(f)/f
    echo "  honeycomb (n1 = n2 = 1):"
    echo "  L    f = 2pi/L^2   s_meas (= meas/cont)  s_exact       |s_meas - s_exact|  resid"
    var devHc: array[3, float]
    let hls = [4, 8, 12]
    var worstPin = 0.0
    for il in 0..<3:
      let L = hls[il]
      let hl = newLayout([L, L, L, L])
      var g = newHcGauge(hl)
      setFluxHc(g, 1, 1)
      var c = newHcWilson(g, 1.0)
      var c0 = newHcWilson(g, 0.0)
      let f1 = 2.0*PI/float(L*L)
      let i0 = hl.rankIndex([1, 2 mod L, 0, 1]).index
      let cm = measureCloverHc(c, c0, hl, i0, 0, 0.1)
      let pin = analyzeClover(cm, f1, f1)
      let sExact = 4.0*sin(f1/4.0)/f1
      devHc[il] = 1.0 - pin.s1
      worstPin = max(worstPin, abs(pin.s1 - sExact))
      worstPin = max(worstPin, abs(pin.s2 - sExact))
      echo &"  {L:2d}   {f1:.6f}     {pin.s1:.12f}    {sExact:.12f}   {abs(pin.s1-sExact):.2e}   {pin.resid:.2e}"
      ok(&"hc L={L}: measured/lattice-analytic = 1 to 1e-9 " &
         &"(dev {abs(pin.s1-sExact):.1e}), structure resid {pin.resid:.1e}",
         abs(pin.s1 - sExact) < 1e-9 and abs(pin.s2 - sExact) < 1e-9 and
         pin.resid < 1e-10)
      if L == 8:
        # sublattice B and the transition-function seam site must agree
        let cmB = measureCloverHc(c, c0, hl, i0, 1, 0.1)
        let iSeam = hl.rankIndex([L-1, 0, 0, 0]).index
        let cmS = measureCloverHc(c, c0, hl, iSeam, 0, 0.1)
        echo &"  L=8 extra: |C(sub B) - C(sub A)| = {maxDiff12(cmB, cm):.2e}, " &
             &"|C(seam) - C(interior)| = {maxDiff12(cmS, cm):.2e}"
        ok("hc clover matrix constant across sublattices and the seam",
           maxDiff12(cmB, cm) < 1e-12 and maxDiff12(cmS, cm) < 1e-12)
        # chirality splitting preview: f1 = f2 kills the gamma5 = +1 sector
        echo &"  L=8 chirality: |P+ C P+| = {pin.npos:.2e}, " &
             &"|P- C P-| = {pin.nneg:.2e}"
    echo &"  artifact 1 - s scaling: dev(4)/dev(8) = {devHc[0]/devHc[1]:.2f} " &
         &"(16 for 1/L^4), dev(8)/dev(12) = {devHc[1]/devHc[2]:.2f} (5.06)"
    ok("hc artifact scales as 1/L^4",
       devHc[0]/devHc[1] > 14.0 and devHc[0]/devHc[1] < 18.0 and
       devHc[1]/devHc[2] > 4.4 and devHc[1]/devHc[2] < 5.8)
    echo "  cubic (n1 = n2 = 1):"
    echo "  L    f = 2pi/L^2   s_meas (= meas/cont)  s_exact       |s_meas - s_exact|  resid"
    var devCb: array[3, float]
    let cls = [8, 12, 16]
    for il in 0..<3:
      let L = cls[il]
      let lo = newLayout(@[L, L, L, L])
      var g = lo.newGauge
      setFluxCubic(g, lo, 1, 1)
      var c = newCubicWilson(g, 1.0)
      let f1 = 2.0*PI/float(L*L)
      let i0 = lo.rankIndex([1, 2, 0, 1]).index
      let cm = measureCloverCubic(c, lo, i0, 0.1)
      let pin = analyzeClover(cm, f1, f1)
      let sExact = sin(f1)/f1
      devCb[il] = 1.0 - pin.s1
      worstPin = max(worstPin, abs(pin.s1 - sExact))
      echo &"  {L:2d}   {f1:.6f}     {pin.s1:.12f}    {sExact:.12f}   {abs(pin.s1-sExact):.2e}   {pin.resid:.2e}"
      ok(&"cubic L={L}: measured/lattice-analytic = 1 to 1e-9 " &
         &"(dev {abs(pin.s1-sExact):.1e}), resid {pin.resid:.1e}",
         abs(pin.s1 - sExact) < 1e-9 and abs(pin.s2 - sExact) < 1e-9 and
         pin.resid < 1e-10)
    echo &"  artifact 1 - s scaling: dev(8)/dev(12) = {devCb[0]/devCb[1]:.2f} " &
         &"((12/8)^4 = 5.06), dev(12)/dev(16) = {devCb[1]/devCb[2]:.2f} ((16/12)^4 = 3.16)"
    ok("cubic artifact scales as 1/L^4 (this also pins cubicFmunuSign = +1)",
       devCb[0]/devCb[1] > 4.4 and devCb[0]/devCb[1] < 5.8 and
       devCb[1]/devCb[2] > 2.8 and devCb[1]/devCb[2] < 3.6)
    # chirality-splitting cross-check on both lattices, both duality signs
    proc chirCase(hc: bool, n2: int): tuple[np, nm, nx, scale: float] =
      if hc:
        let hl = newLayout([8, 8, 8, 8])
        var g = newHcGauge(hl)
        setFluxHc(g, 1, n2)
        var c = newHcWilson(g, 1.0)
        var c0 = newHcWilson(g, 0.0)
        let i0 = hl.rankIndex([1, 2, 0, 1]).index
        let cm = measureCloverHc(c, c0, hl, i0, 0, 0.1)
        let f1 = 2.0*PI/64.0
        let pin = analyzeClover(cm, f1, float(n2)*f1)
        (pin.npos, pin.nneg, pin.ncross, pin.cnorm)
      else:
        let lo = newLayout(@[8, 8, 8, 8])
        var g = lo.newGauge
        setFluxCubic(g, lo, 1, n2)
        var c = newCubicWilson(g, 1.0)
        let i0 = lo.rankIndex([1, 2, 0, 1]).index
        let cm = measureCloverCubic(c, lo, i0, 0.1)
        let f1 = 2.0*PI/64.0
        let pin = analyzeClover(cm, f1, float(n2)*f1)
        (pin.npos, pin.nneg, pin.ncross, pin.cnorm)
    echo "  chirality splitting (L=8): |P+ C P+|, |P- C P-|, cross"
    var chirOk = true
    for (name, hc) in [("honeycomb", true), ("cubic", false)]:
      for n2 in [1, -1]:
        let (np, nm, nx, sc) = chirCase(hc, n2)
        echo &"  {name:9s} n2 = {n2:2d}:  {np:.3e}  {nm:.3e}  {nx:.3e}"
        if nx > 1e-12*sc: chirOk = false
        if n2 == 1:       # f1 = f2: gamma5 = +1 sector must be null
          if np > 1e-9*sc or nm < 0.5*sc: chirOk = false
        else:             # f1 = -f2: gamma5 = -1 sector must be null
          if nm > 1e-9*sc or np < 0.5*sc: chirOk = false
    ok("self-dual flux kills the gamma5 = +1 sector, anti-self-dual the " &
       "gamma5 = -1 sector, identically on both lattices", chirOk)
    ok(&"acceptance pin: worst |measured - lattice-analytic| = {worstPin:.2e}",
       worstPin < 1e-9)

qexFinalize()
