## Task L test: `hclayout.nim` (cell layout + the 16-way binary-tree shift) and
## `hcgauge.nim` (24 link fields + triangle loops).
##
##   export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
##   cd build_mac && make src/experimental/honeycomb/tests/tgauge.nim && ./bin/tgauge
##
## Everything is cross-checked against a brute-force reference built from
## `hcgeom.triPath` and single-site `lo.coord` / `lo.rankIndex` indexing, which
## shares no code with the field-level implementation.

import math, tables, strformat, unittest
import qex except epsilon
import ../hcgeom
import ../hclayout
import ../hcgauge

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))

template ok(msg: string, cond: bool) =
  let c = cond
  if c: echo "PASS: ", msg
  else: echo "FAIL: ", msg
  check c

qexInit()

const nc = getDefaultNc()

# ---------------------------------------------------------------------------
# a completely independent 3x3 complex matrix reference
# ---------------------------------------------------------------------------

type Mat = array[nc, array[nc, array[2, float]]]   ## [row][col][re/im]

proc mid(): Mat =
  for i in 0..<nc: result[i][i][0] = 1.0

proc mmul(a, b: Mat): Mat =
  for i in 0..<nc:
    for j in 0..<nc:
      var re = 0.0
      var im = 0.0
      for k in 0..<nc:
        re += a[i][k][0]*b[k][j][0] - a[i][k][1]*b[k][j][1]
        im += a[i][k][0]*b[k][j][1] + a[i][k][1]*b[k][j][0]
      result[i][j][0] = re
      result[i][j][1] = im

proc mdag(a: Mat): Mat =
  for i in 0..<nc:
    for j in 0..<nc:
      result[i][j][0] = a[j][i][0]
      result[i][j][1] = -a[j][i][1]

proc reTr(a: Mat): float =
  for i in 0..<nc: result += a[i][i][0]

template toF(x: untyped): float =
  ## QEX single-site element accessors return lane proxies, not raw floats
  block:
    var v: float
    v := x
    v

proc getMat(f: auto, idx: int): Mat =
  for a in 0..<nc:
    for b in 0..<nc:
      result[a][b][0] = toF f{idx}[a, b].re
      result[a][b][1] = toF f{idx}[a, b].im

proc wrapCoord(c: Cell, geom: openArray[int]): array[4, cint] =
  for mu in 0..<nDim:
    result[mu] = cint(((c[mu] mod geom[mu]) + geom[mu]) mod geom[mu])

proc siteIndex(lo: Layout, c: Cell, geom: openArray[int]): int =
  var cc = wrapCoord(c, geom)
  let ri = lo.rankIndex(cc)
  doAssert ri.rank == lo.myRank, "tgauge is a single-rank test"
  ri.index

iterator lexCells(geom: openArray[int]): Cell =
  var c: Cell
  for x3 in 0..<geom[3]:
    for x2 in 0..<geom[2]:
      for x1 in 0..<geom[1]:
        for x0 in 0..<geom[0]:
          c = [x0, x1, x2, x3]
          yield c

# ---------------------------------------------------------------------------
# brute-force triangleSum straight from hcgeom.triPath
# ---------------------------------------------------------------------------

proc triangleSumRef(g: auto): tuple[re, im: float] =
  let
    hl = g.hl
    lo = hl.lo
    geom = hl.geom
  var sre = 0.0
  var sim = 0.0
  for y in lexCells(geom):
    for sub in 0..1:
      for t in apexTris:
        var p = mid()
        for l in triPath(Site(cell: y, sub: sub), t):
          var m = getMat(g.link(l.kind, l.idx), siteIndex(lo, l.cell, geom))
          if l.dag: m = mdag(m)
          p = mmul(p, m)
        sre += reTr(p)
        var ti = 0.0
        for i in 0..<nc: ti += p[i][i][1]
        sim += ti
  let n = float(nTriPerSite*hl.nSites*nc)
  (sre/n, sim/n)

# ---------------------------------------------------------------------------
# brute-force check of HcShift16
# ---------------------------------------------------------------------------

proc shift16Error(src: auto, s: auto, geom: openArray[int], sign: int): float =
  let lo = src.l
  var worst = 0.0
  for delta in 0..<nDiag:
    for i in 0..<lo.nSites:
      var c: array[4, cint]
      lo.coord(c, (lo.myRank, i))
      var want: Cell
      for mu in 0..<nDim:
        want[mu] = c[mu].int + sign*((delta shr mu) and 1)
      let j = siteIndex(lo, want, geom)
      let
        a = getMat(s.f[delta], i)
        b = getMat(src, j)
      for p in 0..<nc:
        for q in 0..<nc:
          worst = max(worst, abs(a[p][q][0] - b[p][q][0]))
          worst = max(worst, abs(a[p][q][1] - b[p][q][1]))
  worst

# ---------------------------------------------------------------------------

let
  geom4 = @[4, 4, 4, 6]
  seed = 987654321'u64

suite "hclayout / hcgauge":

  test "1. HcShift16 matches a brute-force gather (SIMD layout)":
    let hl = newHcLayout(geom4)
    let lo = hl.lo
    echo &"  cells = {hl.nCells}, sites = {hl.nSites}, links = {hl.nLinks}, V = {lo.V}"
    ok(&"nCells = {hl.nCells}", hl.nCells == 4*4*4*6)
    ok(&"nSites = {hl.nSites} = 2 nCells", hl.nSites == 2*hl.nCells)
    ok(&"nLinks = {hl.nLinks} = 24 nCells", hl.nLinks == 24*hl.nCells)
    var r = lo.newRNGField(RngMilc6, seed)
    var src = lo.ColorMatrix(nc)
    threads:
      src.randomSU r
    var sf = newHcShift16(src, 1)
    var sb = newHcShift16(src, -1)
    threads:
      sf.run
      sb.run
    let ef = shift16Error(src, sf, geom4, 1)
    let eb = shift16Error(src, sb, geom4, -1)
    echo &"  max |f[delta] - src(y+delta)| = {ef:.3e}"
    echo &"  max |f[delta] - src(y-delta)| = {eb:.3e}"
    ok(&"forward 16-shift exact (err {ef:.1e})", ef == 0.0)
    ok(&"backward 16-shift exact (err {eb:.1e})", eb == 0.0)
    # f[0] aliases the source, so changing the source and re-running suffices
    ok("f[0] aliases src", sf.f[0] == src)
    threads:
      src.randomSU r
      sf.run
    let ef2 = shift16Error(src, sf, geom4, 1)
    ok(&"re-running after the source changed (err {ef2:.1e})", ef2 == 0.0)
    # ... and rebinding to a different field works too
    var src2 = lo.ColorMatrix(nc)
    threads:
      src2.randomSU r
    sf.run src2                      # setSrc + run, outside threads:
    let ef3 = shift16Error(src2, sf, geom4, 1)
    ok(&"run(s, src) rebinds the source (err {ef3:.1e})", ef3 == 0.0)

  test "2. HcShift16 matches a brute-force gather (V = 1 layout)":
    let hl = newHcLayout(geom4, 1)
    let lo = hl.lo
    echo &"  cells = {hl.nCells}, V = {lo.V}"
    var r = lo.newRNGField(RngMilc6, seed)
    var src = lo.ColorMatrix(nc)
    threads:
      src.randomSU r
    var sf = newHcShift16(src, 1)
    var sb = newHcShift16(src, -1)
    threads:
      sf.run
      sb.run
    let ef = shift16Error(src, sf, geom4, 1)
    let eb = shift16Error(src, sb, geom4, -1)
    echo &"  max |f[delta] - src(y+delta)| = {ef:.3e}"
    echo &"  max |f[delta] - src(y-delta)| = {eb:.3e}"
    ok(&"forward 16-shift exact (err {ef:.1e})", ef == 0.0)
    ok(&"backward 16-shift exact (err {eb:.1e})", eb == 0.0)

  test "3. unit gauge gives triangleSum = 1":
    let hl = newHcLayout(geom4)
    var g = newHcGauge(hl)
    let s = g.triangleSum
    let si = g.triangleSumIm
    echo &"  triangleSum(1) = {s:.17g}, Im = {si:.3e}"
    ok(&"triangleSum = 1 to 1e-14 (|1-s| = {abs(1.0-s):.2e})", abs(1.0 - s) < 1e-14)
    ok(&"Im part vanishes ({abs(si):.2e})", abs(si) < 1e-14)
    let (rr, ri) = triangleSumRef(g)
    ok(&"brute-force reference agrees ({abs(rr-s):.2e})", abs(rr - s) < 1e-14)
    ok(&"reference Im vanishes ({abs(ri):.2e})", abs(ri) < 1e-14)

  test "4. field triangleSum == brute-force hcgeom.triPath reference":
    let hl = newHcLayout(geom4)
    let lo = hl.lo
    var r = lo.newRNGField(RngMilc6, seed)
    var g = newHcGauge(hl)
    threads:
      g.random r
    let s = g.triangleSum
    let si = g.triangleSumIm
    let (rr, ri) = triangleSumRef(g)
    echo &"  triangleSum       = {s:.15g}  (Im {si:.3e})"
    echo &"  triangleSum (ref) = {rr:.15g}  (Im {ri:.3e})"
    echo &"  difference        = {abs(s-rr):.3e}"
    ok(&"field and reference agree to 1e-12 (diff {abs(s-rr):.2e})",
       abs(s - rr) < 1e-12)
    ok(&"imaginary parts agree ({abs(si-ri):.2e})", abs(si - ri) < 1e-12)
    let cs = g.checkSU
    echo &"  checkSU: avg {cs.avg:.3e}, max {cs.max:.3e}"
    ok(&"random config is SU(3) (max {cs.max:.2e})", cs.max < 1e-9)

  test "5. gauge invariance of triangleSum":
    let hl = newHcLayout(geom4)
    let lo = hl.lo
    var r = lo.newRNGField(RngMilc6, seed + 1)
    var g = newHcGauge(hl)
    var g0 = newHcGauge(hl)
    threads:
      g.random r
      g0 := g
    let before = g.triangleSum
    var vA = lo.ColorMatrix(nc)
    var vB = lo.ColorMatrix(nc)
    threads:
      vA.randomSU r
      vB.randomSU r
    g.gaugeTransform(vA, vB)
    let after = g.triangleSum
    let d = abs(after - before)
    echo &"  triangleSum before = {before:.17g}"
    echo &"  triangleSum after  = {after:.17g}"
    echo &"  |difference|       = {d:.3e}"
    ok(&"gauge invariant to 1e-12 (diff {d:.2e})", d < 1e-12)
    # ... and the transformation really did change every link
    var t = g.uA[0].newOneOf
    var dn = newSeq[float](nDirs)
    let ga = g.allLinks
    let g0a = g0.allLinks
    threads:
      for k in 0..<nDirs:
        t := ga[k] - g0a[k]
        let v = t.norm2
        if threadNum == 0: dn[k] = v
    var tot = 0.0
    var minChange = 1e300
    for v in dn:
      tot += v
      minChange = min(minChange, v)
    echo &"  sum_links |U' - U|^2 = {tot:.6g}, smallest per-field = {minChange:.6g}"
    ok(&"gaugeTransform changed all 24 link fields (min {minChange:.3g})",
       minChange > 1.0)
    let cs = g.checkSU
    echo &"  checkSU after transform: avg {cs.avg:.3e}, max {cs.max:.3e}"
    ok(&"transformed links are still SU(3) (max {cs.max:.2e})", cs.max < 1e-9)
    # and the reference agrees with the transformed configuration too
    let (rr, _) = triangleSumRef(g)
    ok(&"reference agrees after the transform ({abs(rr-after):.2e})",
       abs(rr - after) < 1e-12)
    # a second, non-trivial check: transforming a *unit* configuration must
    # still give exactly 1
    var u = newHcGauge(hl)
    u.gaugeTransform(vA, vB)
    let su = u.triangleSum
    echo &"  triangleSum(gauge-transformed unit) = {su:.17g}"
    ok(&"pure gauge stays 1 ({abs(1.0-su):.2e})", abs(1.0 - su) < 1e-12)

  test "6. triangleSum of a random configuration is small":
    let hl = newHcLayout(geom4)
    let lo = hl.lo
    var r = lo.newRNGField(RngMilc6, seed + 2)
    var g = newHcGauge(hl)
    var vals: seq[float]
    var maxRe = 0.0
    var maxIm = 0.0
    for k in 0..<8:
      threads:
        g.random r
      let v = g.triangleSum
      vals.add v
      maxRe = max(maxRe, abs(v))
      maxIm = max(maxIm, abs(g.triangleSumIm))
    var mean = 0.0
    for v in vals: mean += v
    mean /= float(vals.len)
    echo &"  triangleSum over {vals.len} random configs:"
    for v in vals: echo &"    {v: .9f}"
    echo &"  mean = {mean:.3e}, max |Re| = {maxRe:.3e}, max |Im| = {maxIm:.3e}"
    # 64*nCells triangles, each ~ Re Tr(random SU(3))/3 with variance 1/(2 Nc^2)
    let sigma = sqrt(1.0/(2.0*nc*nc*float(64*hl.nCells)))
    echo &"  expected size ~ {sigma:.3e} (1 sigma)"
    ok(&"|mean| = {abs(mean):.3e} << 1", abs(mean) < 10.0*sigma)
    ok(&"every value within 5 sigma of 0 (max {maxRe:.2e})", maxRe < 5.0*sigma)
    # NOTE: the 32 apex triangles are enumerated with a *fixed* orientation, so
    # unlike the cubic plaquette sum there is no exact cancellation of the
    # imaginary part; it is only statistically zero, the same size as the real
    # part.  It *is* exactly zero for unit and pure-gauge configurations
    # (checked in tests 3 and 5).
    ok(&"|Im| is the same statistical size as |Re| ({maxIm:.2e})",
       maxIm < 5.0*sigma)

  test "7. every link of the cell torus is in exactly 8 triangles":
    let geom = @[4, 4, 4, 4]
    var tally = initTable[(int, int, int, int, int, int), int]()
    var nTri = 0
    for y in lexCells(geom):
      for sub in 0..1:
        for t in apexTris:
          inc nTri
          for l in triPath(Site(cell: y, sub: sub), t):
            let c = wrapCoord(l.cell, geom)
            let k = (l.kind.ord, l.idx, c[0].int, c[1].int, c[2].int, c[3].int)
            tally[k] = tally.getOrDefault(k) + 1
    var nCells = 1
    for x in geom: nCells *= x
    echo &"  {nTri} triangles on {nCells} cells = {nTri div nCells} per cell"
    echo &"  {tally.len} distinct links (expected {24*nCells})"
    ok(&"64 triangles per cell", nTri == 64*nCells)
    ok(&"all {24*nCells} links appear", tally.len == 24*nCells)
    var allEight = true
    for _, v in tally:
      if v != 8: allEight = false
    ok("every link is contained in exactly 8 triangles", allEight)

  test "8. HcGauge helpers: unit / warm / reunit / link accessor":
    let hl = newHcLayout(@[4, 4, 4, 4])
    let lo = hl.lo
    var r = lo.newRNGField(RngMilc6, seed + 3)
    var g = newHcGauge(hl)
    threads:
      g.warm(0.3, r)
    let sw = g.triangleSum
    echo &"  triangleSum after warm(0.3) = {sw:.6f}"
    ok(&"warm start is near but not at 1 ({sw:.4f})", sw > 0.5 and sw < 1.0)
    threads:
      g.unit
    ok("unit() restores triangleSum = 1", abs(1.0 - g.triangleSum) < 1e-14)
    threads:
      g.random r
      g.reunit
    let cs = g.checkSU
    ok(&"reunit keeps SU(3) (max {cs.max:.2e})", cs.max < 1e-13)
    var nMatch = 0
    let al = g.allLinks
    for mu in 0..<nDim:
      if g.link(lkA, mu) == al[mu]: inc nMatch
      if g.link(lkB, mu) == al[nDim+mu]: inc nMatch
    for d in 0..<nDiag:
      if g.link(lkD, d) == al[2*nDim+d]: inc nMatch
    ok(&"link() accessor consistent with allLinks ({nMatch}/24)", nMatch == 24)

  test "9. triangleSum vs reference on a V = 1 layout":
    let hl = newHcLayout(geom4, 1)
    let lo = hl.lo
    var r = lo.newRNGField(RngMilc6, seed + 4)
    var g = newHcGauge(hl)
    threads:
      g.random r
    let s = g.triangleSum
    let (rr, _) = triangleSumRef(g)
    echo &"  V = 1: triangleSum = {s:.15g}, reference = {rr:.15g}"
    ok(&"agree to 1e-12 (diff {abs(s-rr):.2e})", abs(s - rr) < 1e-12)
    ok(&"unit gauge still 1", abs(1.0 - (block:
      var u = newHcGauge(hl)
      u.triangleSum)) < 1e-14)

  test "10. triangleSum vs reference on a short (L=2) direction":
    # L = 2 makes the +delta shifts wrap around repeatedly; a good stress test
    # for the binary-tree shift and for the boundary paths of the shifters.
    let geom = @[2, 4, 4, 6]
    let hl = newHcLayout(geom)
    let lo = hl.lo
    var r = lo.newRNGField(RngMilc6, seed + 5)
    var g = newHcGauge(hl)
    threads:
      g.random r
    let s = g.triangleSum
    let (rr, _) = triangleSumRef(g)
    echo &"  geom {geom}: triangleSum = {s:.15g}, reference = {rr:.15g}"
    ok(&"agree to 1e-12 (diff {abs(s-rr):.2e})", abs(s - rr) < 1e-12)
    var src = lo.ColorMatrix(nc)
    threads:
      src.randomSU r
    var sf = newHcShift16(src, 1)
    var sb = newHcShift16(src, -1)
    threads:
      sf.run
      sb.run
    let ef = shift16Error(src, sf, geom, 1)
    let eb = shift16Error(src, sb, geom, -1)
    ok(&"forward 16-shift exact with L=2 (err {ef:.1e})", ef == 0.0)
    ok(&"backward 16-shift exact with L=2 (err {eb:.1e})", eb == 0.0)
    var vA = lo.ColorMatrix(nc)
    var vB = lo.ColorMatrix(nc)
    threads:
      vA.randomSU r
      vB.randomSU r
    g.gaugeTransform(vA, vB)
    let after = g.triangleSum
    ok(&"gauge invariant with L=2 (diff {abs(after-s):.2e})", abs(after - s) < 1e-12)

qexFinalize()
