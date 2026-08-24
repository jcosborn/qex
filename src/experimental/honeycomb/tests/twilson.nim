## Task D1 test: `hcwilson.nim`, the interacting Wilson-Dirac operator on the
## 16-cell honeycomb.
##
##   export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
##   cd build_mac && make src/experimental/honeycomb/tests/twilson.nim && \
##     OMP_NUM_THREADS=4 ./bin/twilson
##
## The make-or-break check is test 2: on a unit gauge field the position-space
## operator, projected onto plane waves, must reproduce `hcfree.freeD8`'s 8x8
## (spin x sublattice) cell-momentum blocks ENTRYWISE.  Plane-wave convention
## (documented in hcwilson.nim): psiA(y) = e_s exp(i k.y) and
## psiB(y) = e_s exp(i k.y), both phased with the integer CELL coordinate y --
## no half-site offset phase on the B sublattice.
##
## QEX dot convention (established empirically in test 1): dot(x, y) conjugates
## the FIRST argument, dot(x, y) = sum conj(x)_i y_i.

import math, strformat, unittest
import std/[complex, times, os]
import qex except epsilon
import ../hcgeom
import ../hclayout
import ../hcgauge
import ../hcwilson
import ../hcfree

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))

template ok(msg: string, cond: bool) =
  let c = cond
  if c: echo "PASS: ", msg
  else: echo "FAIL: ", msg
  check c

qexInit()

const nc = getDefaultNc()
let seed = 987654321'u64

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

template toF(x: untyped): float =
  ## QEX single-site element accessors return lane proxies, not raw floats
  block:
    var v: float
    v := x
    v

proc cdot(x, y: auto): Complex64 =
  ## <x|y> with the first argument conjugated (QEX's convention, test 1)
  let z = dot(x, y)
  complex64(z.re, z.im)

proc setPlaneWave(f: auto, k: array[4, float], s, sub: int) =
  ## f := plane wave e_s exp(i k.y) on sublattice `sub` (cell-coordinate
  ## phase), zero elsewhere; single site writes, call outside threads
  threads:
    f.a := 0
    f.b := 0
  let lo = f.a.l
  for i in 0..<lo.nSites:
    var c: array[4, cint]
    lo.coord(c, (lo.myRank, i))
    var ph = 0.0
    for mu in 0..<4: ph += k[mu]*float(c[mu])
    if sub == 0:
      f.a{i}[s][0].re := cos(ph)
      f.a{i}[s][0].im := sin(ph)
    else:
      f.b{i}[s][0].re := cos(ph)
      f.b{i}[s][0].im := sin(ph)

proc measureBlock(w: auto, basis: auto, work: var HcFermion,
                  k: array[4, float], m, rw: float, dag = false): Mat8 =
  ## the 8x8 (spin x sublattice) matrix <u_i, D u_j>/V on plane waves at cell
  ## momentum k; index = 2*spin + sub, matching hcfree.freeD8
  let vol = float(basis[0].a.l.physVol)
  for j in 0..<8:
    setPlaneWave(basis[j], k, j div 2, j and 1)
  for j in 0..<8:
    if dag: w.Ddag(work, basis[j], m, rw)
    else: w.D(work, basis[j], m, rw)
    for i in 0..<8:
      let z = cdot(basis[i], work)
      result[i][j] = complex64(z.re/vol, z.im/vol)

proc maxEntryDiff(a, b: Mat8): float =
  for i in 0..<8:
    for j in 0..<8:
      result = max(result, abs(a[i][j] - b[i][j]))

proc wrapCoord(c: Cell, geom: openArray[int]): array[4, int] =
  for mu in 0..<nDim:
    result[mu] = ((c[mu] mod geom[mu]) + geom[mu]) mod geom[mu]

proc spinor(f: auto, i, c: int): array[4, Complex64] =
  ## the 4 spin components of color c at site index i
  for s in 0..<4:
    result[s] = complex64(toF f{i}[s][c].re, toF f{i}[s][c].im)

proc siteNorm2(f: auto, i: int): float =
  for s in 0..<4:
    for c in 0..<nc:
      result += toF(f{i}[s][c].re)^2 + toF(f{i}[s][c].im)^2

# ---------------------------------------------------------------------------

suite "hcwilson (task D1)":

  test "1. conventions: QEX dot conjugation; spinOld vs hcfree gammas":
    let hl = newHcLayout([4, 4, 4, 4])
    let lo = hl.lo
    var x = newHcFermion(hl)
    var y = newHcFermion(hl)
    setPlaneWave(x, [0.0, 0.0, 0.0, 0.0], 0, 0)
    threads:
      for e in y.a:
        y.a[e] := newComplex(0.0, 1.0) * x.a[e]
      y.b := 0
    let vol = float(lo.physVol)
    let d = cdot(x, y)
    echo &"  dot(x, i x)/V = {d.re/vol:.3f} + {d.im/vol:.3f} i"
    ok("dot(x, ix) = +i |x|^2: QEX dot conjugates its FIRST argument",
       abs(d.re) < 1e-12*vol and abs(d.im - vol) < 1e-12*vol)
    # spinOld's gamma1..4 match hcfree.gammaMat entrywise
    var wg = 0.0
    template cmpGamma(gs: typed, mu: int) =
      for i in 0..<4:
        for j in 0..<4:
          wg = max(wg, abs(gs[i, j].re - gammaMat[mu][i][j].re))
          wg = max(wg, abs(gs[i, j].im - gammaMat[mu][i][j].im))
    cmpGamma(gamma1, 0)
    cmpGamma(gamma2, 1)
    cmpGamma(gamma3, 2)
    cmpGamma(gamma4, 3)
    ok(&"spinOld gamma1..4 == hcfree.gammaMat ({wg:.1e})", wg == 0.0)
    # hcwilson.hopMat matches (gammaDot - r)/6 built from hcfree, all 24 dirs
    var wh = 0.0
    let rw = 0.8
    for dir in 0..<nDirs:
      let hm = hopMat(dir, rw)
      let gd = gammaDot(hcgeom.toFloat(dirVec(dir)))
      for i in 0..<4:
        for j in 0..<4:
          var z = gd[i][j]
          if i == j: z = z - complex64(rw)
          z = z/complex64(6.0)
          wh = max(wh, abs(complex64(toF hm[i, j].re, toF hm[i, j].im) - z))
    ok(&"hopMat(dir, r) == (gamma.n - r)/6 for all 24 dirs ({wh:.1e})",
       wh < 1e-15)

  test "2. free field: D == hcfree.freeD8 entrywise (4^4, periodic)":
    let geom = @[4, 4, 4, 4]
    let hl = newHcLayout(geom)
    var g = newHcGauge(hl)           # unit links
    var w = newHcWilson(g)
    var basis: array[8, typeof(newHcFermion(hl))]
    for j in 0..<8: basis[j] = newHcFermion(hl)
    var work = newHcFermion(hl)
    let ns = [
      [0, 0, 0, 0], [2, 2, 2, 2], [1, 0, 0, 0], [0, 0, 0, 1], [1, 2, 3, 0],
      [3, 3, 1, 2], [1, 1, 1, 1], [2, 0, 1, 3], [3, 1, 0, 2], [0, 3, 2, 1]]
    let mrw = [(0.0, 1.0), (0.1, 1.0), (0.3, 0.8), (-0.2, 1.0), (0.0, 0.7)]
    var worst = 0.0
    for t in 0..<ns.len:
      var k: array[4, float]
      for mu in 0..<4: k[mu] = 2.0*PI*float(ns[t][mu])/float(geom[mu])
      let (m, rw) = mrw[t mod mrw.len]
      let got = measureBlock(w, basis, work, k, m, rw)
      let want = freeD8(k, rw, m)
      let d = maxEntryDiff(got, want)
      worst = max(worst, d)
      echo &"  k = 2pi/4*({ns[t][0]},{ns[t][1]},{ns[t][2]},{ns[t][3]}), " &
           &"m = {m}, r = {rw}: max entry diff = {d:.3e}"
    ok(&"all 8x8 blocks match freeD8 entrywise, worst = {worst:.3e}",
       worst < 1e-12)
    # and the same for Ddag against the adjoint block, one momentum
    block:
      var k: array[4, float]
      for mu in 0..<4: k[mu] = 2.0*PI*float(ns[4][mu])/float(geom[mu])
      let got = measureBlock(w, basis, work, k, 0.1, 0.9, dag = true)
      let want = adj(freeD8(k, 0.9, 0.1))
      let d = maxEntryDiff(got, want)
      ok(&"Ddag block == freeD8^dag entrywise ({d:.3e})", d < 1e-12)

  test "3. free field, asymmetric geometry [4,4,2,6]":
    let geom = @[4, 4, 2, 6]
    let hl = newHcLayout(geom)
    var g = newHcGauge(hl)
    var w = newHcWilson(g)
    var basis: array[8, typeof(newHcFermion(hl))]
    for j in 0..<8: basis[j] = newHcFermion(hl)
    var work = newHcFermion(hl)
    let ns = [[0, 0, 0, 0], [1, 3, 1, 5], [2, 1, 0, 3], [3, 2, 1, 1],
              [1, 0, 1, 2]]
    var worst = 0.0
    for t in 0..<ns.len:
      var k: array[4, float]
      for mu in 0..<4: k[mu] = 2.0*PI*float(ns[t][mu])/float(geom[mu])
      let got = measureBlock(w, basis, work, k, 0.05, 1.0)
      let want = freeD8(k, 1.0, 0.05)
      worst = max(worst, maxEntryDiff(got, want))
    echo &"  worst entry diff over {ns.len} momenta = {worst:.3e}"
    ok(&"asymmetric geometry blocks match ({worst:.3e})", worst < 1e-12)

  test "3b. free field on a V = 1 (non-SIMD) layout":
    let geom = @[4, 4, 4, 4]
    let hl = newHcLayout(geom, 1)
    var g = newHcGauge(hl)
    var w = newHcWilson(g)
    var basis: array[8, typeof(newHcFermion(hl))]
    for j in 0..<8: basis[j] = newHcFermion(hl)
    var work = newHcFermion(hl)
    let ns = [[1, 2, 3, 0], [3, 3, 1, 2]]
    var worst = 0.0
    for t in 0..<ns.len:
      var k: array[4, float]
      for mu in 0..<4: k[mu] = 2.0*PI*float(ns[t][mu])/float(geom[mu])
      let got = measureBlock(w, basis, work, k, 0.1, 1.0)
      worst = max(worst, maxEntryDiff(got, freeD8(k, 1.0, 0.1)))
    echo &"  V = 1: worst entry diff = {worst:.3e}"
    ok(&"V = 1 layout matches freeD8 ({worst:.3e})", worst < 1e-12)

  test "4. locality: point source -> exactly 25 sites, correct values":
    let geom = @[4, 4, 4, 4]
    let hl = newHcLayout(geom)
    let lo = hl.lo
    var g = newHcGauge(hl)
    var w = newHcWilson(g)
    var x = newHcFermion(hl)
    var r = newHcFermion(hl)
    let m = 0.2
    let rw = 0.9
    for srcSub in 0..1:
      let c0 = [0, 3, 0, 3]              # near the boundary: exercise wraps
      let s0 = 1                         # source spin
      threads:
        x := 0
      let i0 = lo.rankIndex(c0).index
      if srcSub == 0: x.a{i0}[s0][0].re := 1.0
      else: x.b{i0}[s0][0].re := 1.0
      w.D(r, x, m, rw)
      # expected: (m+4r) at the source, (1/6)(gamma.n_i - r) e_s0 at src - n_i
      var expApos: seq[(array[4, int], array[4, Complex64])]
      var expBpos: seq[(array[4, int], array[4, Complex64])]
      block:
        var v: array[4, Complex64]
        v[s0] = complex64(m + 4.0*rw)
        if srcSub == 0: expApos.add (wrapCoord(c0, geom), v)
        else: expBpos.add (wrapCoord(c0, geom), v)
      for dir in 0..<nDirs:
        let (_, dest) = step(Site(cell: c0, sub: srcSub), opposite(dir))
        let gd = gammaDot(hcgeom.toFloat(dirVec(dir)))
        var v: array[4, Complex64]
        for s in 0..<4:
          var z = gd[s][s0]
          if s == s0: z = z - complex64(rw)
          v[s] = z/complex64(6.0)
        if dest.sub == 0: expApos.add (wrapCoord(dest.cell, geom), v)
        else: expBpos.add (wrapCoord(dest.cell, geom), v)
      ok(&"sub {srcSub}: expected support is 25 sites " &
         &"({expApos.len} A + {expBpos.len} B)",
         expApos.len + expBpos.len == 25 and
         (if srcSub == 0: expApos.len == 9 else: expBpos.len == 9))
      # scan every site of the result
      var nnz = 0
      var worstVal = 0.0
      var badSite = false
      for i in 0..<lo.nSites:
        var c: array[4, cint]
        lo.coord(c, (lo.myRank, i))
        let cc = [c[0].int, c[1].int, c[2].int, c[3].int]
        for sub in 0..1:
          let f = (if sub == 0: r.a else: r.b)
          let n2 = siteNorm2(f, i)
          # find the expectation for this site, if any
          var want: array[4, Complex64]
          var listed = false
          for (p, v) in (if sub == 0: expApos else: expBpos):
            if p == cc:
              want = v
              listed = true
              break
          if n2 > 1e-24:
            inc nnz
            if not listed: badSite = true
          # entrywise check (zero sites checked against zero, color 0)
          let got = spinor(f, i, 0)
          for s in 0..<4:
            worstVal = max(worstVal, abs(got[s] - want[s]))
          for cl in 1..<nc:                # colors != 0 must vanish
            let gc = spinor(f, i, cl)
            for s in 0..<4:
              worstVal = max(worstVal, abs(gc[s]))
      echo &"  sub {srcSub}: {nnz} nonzero sites, worst value dev = {worstVal:.3e}"
      ok(&"sub {srcSub}: support = 25 sites exactly", nnz == 25 and not badSite)
      ok(&"sub {srcSub}: all site values match (gamma.n - r)/6 columns " &
         &"({worstVal:.3e})", worstVal < 1e-14)

  test "5. gauge covariance: D[U^g](V psi) = V (D[U] psi)":
    let hl = newHcLayout([4, 4, 4, 6])
    let lo = hl.lo
    var rs = lo.newRNGField(RngMilc6, seed)
    var g = newHcGauge(hl)
    var g2 = newHcGauge(hl)
    threads:
      g.random rs
      g2 := g
    var w = newHcWilson(g)
    var x = newHcFermion(hl)
    var xg = newHcFermion(hl)
    var y1 = newHcFermion(hl)
    var y2 = newHcFermion(hl)
    var yv = newHcFermion(hl)
    threads:
      x.gaussian rs
    let m = 0.1
    let rw = 1.0
    w.D(y1, x, m, rw)
    var vA = lo.ColorMatrix(nc)
    var vB = lo.ColorMatrix(nc)
    threads:
      vA.randomSU rs
      vB.randomSU rs
    g2.gaugeTransform(vA, vB)
    var w2 = newHcWilson(g2)
    threads:
      for e in xg.a:
        xg.a[e] := vA[e] * x.a[e]
      for e in xg.b:
        xg.b[e] := vB[e] * x.b[e]
    w2.D(y2, xg, m, rw)
    threads:
      for e in yv.a:
        yv.a[e] := vA[e] * y1.a[e]
      for e in yv.b:
        yv.b[e] := vB[e] * y1.b[e]
    let rel = sqrt(norm2diff(y2, yv)/norm2(yv))
    echo &"  |D[U^g] V psi - V D[U] psi| / |V D psi| = {rel:.3e}"
    ok(&"gauge covariant to 1e-12 ({rel:.3e})", rel < 1e-12)
    # the transform is non-trivial: transformed application differs from raw
    let relRaw = sqrt(norm2diff(y2, y1)/norm2(y1))
    ok(&"transformed result differs from untransformed ({relRaw:.2f})",
       relRaw > 0.1)

  test "6. gamma5-hermiticity and the adjoint":
    let hl = newHcLayout([4, 4, 4, 4])
    let lo = hl.lo
    var rs = lo.newRNGField(RngMilc6, seed + 7)
    var g = newHcGauge(hl)
    threads:
      g.random rs
    var w = newHcWilson(g)
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
    let rw = 0.85
    # <x, D y> == <Ddag x, y>  (dot conjugates the first argument)
    w.D(dy, y, m, rw)
    w.Ddag(dx, x, m, rw)
    let z1 = cdot(x, dy)
    let z2 = cdot(dx, y)
    let relAdj = abs(z1 - z2)/abs(z1)
    echo &"  <x, D y>      = {z1.re:.12f} + {z1.im:.12f} i"
    echo &"  <Ddag x, y>   = {z2.re:.12f} + {z2.im:.12f} i"
    ok(&"Ddag is the adjoint of D ({relAdj:.3e})", relAdj < 1e-12)
    # gamma5 D gamma5 x == Ddag x, as fields
    threads:
      t1.applyGamma5 x
    w.D(t2, t1, m, rw)
    threads:
      t1.applyGamma5 t2
    let relG5 = sqrt(norm2diff(t1, dx)/norm2(dx))
    echo &"  |g5 D g5 x - Ddag x| / |Ddag x| = {relG5:.3e}"
    ok(&"gamma5 D gamma5 = Ddag ({relG5:.3e})", relG5 < 1e-12)
    # and directly: <x, g5 D g5 y> == <D x, y>
    threads:
      t1.applyGamma5 y
    w.D(t2, t1, m, rw)
    threads:
      t1.applyGamma5 t2
    w.D(dy, x, m, rw)
    let z3 = cdot(x, t1)
    let z4 = cdot(dy, y)
    let relS = abs(z3 - z4)/abs(z3)
    ok(&"<x, g5 D g5 y> = <D x, y> ({relS:.3e})", relS < 1e-12)

  test "7. antiperiodic BC: block at k equals freeD8 at k3 + pi/Nt":
    let geom = @[4, 4, 4, 4]
    let hl = newHcLayout(geom)
    var g = newHcGauge(hl)           # unit links ...
    threads:
      g.setBC                        # ... made antiperiodic
    var w = newHcWilson(g)           # (refresh built in: setBC came first)
    var basis: array[8, typeof(newHcFermion(hl))]
    for j in 0..<8: basis[j] = newHcFermion(hl)
    var work = newHcFermion(hl)
    let ns = [[0, 0, 0, 0], [1, 2, 3, 0], [0, 0, 0, 3], [2, 1, 0, 2]]
    var worst = 0.0
    for t in 0..<ns.len:
      var k: array[4, float]
      for mu in 0..<3: k[mu] = 2.0*PI*float(ns[t][mu])/float(geom[mu])
      k[3] = 2.0*PI*(float(ns[t][3]) + 0.5)/float(geom[3])   # k3 -> k3 + pi/Nt
      let got = measureBlock(w, basis, work, k, 0.1, 1.0)
      let want = freeD8(k, 1.0, 0.1)
      let d = maxEntryDiff(got, want)
      worst = max(worst, d)
      echo &"  n = ({ns[t][0]},{ns[t][1]},{ns[t][2]},{ns[t][3]}+1/2): " &
           &"max entry diff = {d:.3e}"
    ok(&"antiperiodic blocks match freeD8 at shifted k3 ({worst:.3e})",
       worst < 1e-12)
    # sanity: at an UNshifted (periodic) momentum the antiperiodic operator
    # must NOT reproduce the periodic block
    block:
      let k = [0.0, 0.0, 0.0, 2.0*PI*1.0/float(geom[3])]
      let got = measureBlock(w, basis, work, k, 0.1, 1.0)
      let want = freeD8(k, 1.0, 0.1)
      let d = maxEntryDiff(got, want)
      ok(&"setBC really changed the operator (diff at periodic k = {d:.2e})",
         d > 1e-3)

  test "8. mass linearity: D(m1) - D(m0) = (m1 - m0) * 1":
    let hl = newHcLayout([4, 4, 4, 4])
    let lo = hl.lo
    var rs = lo.newRNGField(RngMilc6, seed + 11)
    var g = newHcGauge(hl)
    threads:
      g.random rs
    var w = newHcWilson(g)
    var x = newHcFermion(hl)
    var p = newHcFermion(hl)
    var q = newHcFermion(hl)
    var t = newHcFermion(hl)
    threads:
      x.gaussian rs
    let m0 = 0.1
    let m1 = 0.7
    w.D(p, x, m0)
    w.D(q, x, m1)
    let dm = m1 - m0
    threads:
      for e in t.a:
        t.a[e] := q.a[e] - p.a[e] - dm*x.a[e]
      for e in t.b:
        t.b[e] := q.b[e] - p.b[e] - dm*x.b[e]
    let rel = sqrt(norm2(t)/norm2(x))
    echo &"  |D(m1)x - D(m0)x - (m1-m0)x| / |x| = {rel:.3e}"
    ok(&"mass linearity ({rel:.3e})", rel < 1e-13)

  test "9. timing: D on 8^4 cells":
    let hl = newHcLayout([8, 8, 8, 8])
    let lo = hl.lo
    var rs = lo.newRNGField(RngMilc6, seed + 13)
    var g = newHcGauge(hl)
    threads:
      g.random rs
    var w = newHcWilson(g)
    var x = newHcFermion(hl)
    var r = newHcFermion(hl)
    threads:
      x.gaussian rs
    w.D(r, x, 0.1)                   # warm up (also builds hop matrices)
    w.D(r, x, 0.1)
    # the machine is heavily shared: take the fastest of several batches
    let nrep = 5
    var dt = 1e30
    for batch in 0..<5:
      let t0 = epochTime()
      for it in 0..<nrep:
        w.D(r, x, 0.1)
      dt = min(dt, (epochTime() - t0)/float(nrep))
    let vol = float(hl.nCells)
    echo &"  D on {hl.nCells} cells: {1e3*dt:.2f} ms/application, best of 5 " &
         &"batches ({1e9*dt/vol:.1f} ns/cell, OMP_NUM_THREADS = " &
         getEnv("OMP_NUM_THREADS") & ")"
    var dtg = 1e30
    for batch in 0..<3:
      let t0 = epochTime()
      w.gaugeRefresh
      dtg = min(dtg, epochTime() - t0)
    echo &"  gaugeRefresh (once per configuration): {1e3*dtg:.2f} ms"
    ok("timing measured", dt > 0.0)

qexFinalize()
