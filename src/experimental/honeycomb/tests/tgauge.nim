## `hcgauge.nim` (cell layout, the 16-way binary-tree shift, the
## 24 link fields and the triangle loops).
##
## Everything is cross-checked against a brute-force reference built from
## `hcgeom.triPath` and single-site `lo.coord` / `lo.rankIndex` indexing, which
## shares no code with the field-level implementation.

import std/[math, strformat, unittest]
import qex except epsilon
import ../hcgeom
import ../hcgauge
import helpers

qexInit()

proc triangleSumIm(g: auto): float =
  triangleTrace(g).im/float(nTriPerSite*2*g.lo.physVol*nc)

proc triangleSumRef(g: auto): tuple[re, im: float] =
  ## brute-force triangleSum straight from hcgeom.triPath
  let
    lo = g.lo
    geom = lo.physGeom
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
        sim += imTr(p)
  let n = float(nTriPerSite*2*lo.physVol*nc)
  (sre/n, sim/n)

proc shift16Error(src: auto, s: auto, geom: openArray[int], sign: int): float =
  ## max |f[delta](y) - src(y + sign delta)| over delta and sites
  let lo = src.l
  for delta in 0..<nDiag:
    for i in 0..<lo.nSites:
      var c: array[4, cint]
      lo.coord(c, (lo.myRank, i))
      var want: Cell
      for mu in 0..<nDim:
        want[mu] = c[mu].int + sign*((delta shr mu) and 1)
      result = max(result, mmaxdiff(getMat(s.f[delta], i),
                                    getMat(src, siteIndex(lo, want, geom))))

let
  geom4 = @[4, 4, 4, 6]
  seed = 987654321'u64

proc altLayout(lo: auto, geom: seq[int]): tuple[dtri, ef, eb: float] =
  ## |triangleSum - reference| and the forward/backward 16-shift errors on lo
  var r = lo.newRNGField(RngMilc6, seed + 4)
  var g = newHcGauge(lo)
  threads:
    g.random r
  let (rr, _) = triangleSumRef(g)
  var src = lo.ColorMatrix(nc)
  threads:
    src.randomSU r
  var sf = newHcShift16(src, 1)
  var sb = newHcShift16(src, -1)
  threads:
    sf.run
    sb.run
  (abs(g.triangleSum - rr), shift16Error(src, sf, geom, 1),
   shift16Error(src, sb, geom, -1))

suite "hcgauge":

  test "1. HcShift16 matches a brute-force gather (SIMD layout)":
    let hl = newLayout(geom4)
    let lo = hl
    echo &"  cells = {hl.physVol}, sites = {2*hl.physVol}, links = {24*hl.physVol}, V = {lo.V}"
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
    sf.setSrc src2
    threads:
      sf.run
    let ef3 = shift16Error(src2, sf, geom4, 1)
    ok(&"setSrc rebinds the source (err {ef3:.1e})", ef3 == 0.0)

  test "2. unit gauge gives triangleSum = 1":
    let hl = newLayout(geom4)
    var g = newHcGauge(hl)
    let s = g.triangleSum
    let si = g.triangleSumIm
    echo &"  triangleSum(1) = {s:.17g}, Im = {si:.3e}"
    ok(&"triangleSum = 1 to 1e-14 (|1-s| = {abs(1.0-s):.2e})", abs(1.0 - s) < 1e-14)
    ok(&"Im part vanishes ({abs(si):.2e})", abs(si) < 1e-14)
    let (rr, ri) = triangleSumRef(g)
    ok(&"brute-force reference agrees ({abs(rr-s):.2e})", abs(rr - s) < 1e-14)
    ok(&"reference Im vanishes ({abs(ri):.2e})", abs(ri) < 1e-14)

  test "3. field triangleSum == brute-force hcgeom.triPath reference":
    let hl = newLayout(geom4)
    let lo = hl
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

  test "4. gauge invariance of triangleSum":
    let hl = newLayout(geom4)
    let lo = hl
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
    var dn = newSeq[float](nDirs)
    for k in 0..<nDirs:
      dn[k] = norm2diff(g.links[k], g0.links[k])
    let minChange = min(dn)
    echo &"  sum_links |U' - U|^2 = {sum(dn):.6g}, smallest per-field = {minChange:.6g}"
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

  test "5. triangleSum of a random configuration is small":
    let hl = newLayout(geom4)
    let lo = hl
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
    let mean = sum(vals)/float(vals.len)
    echo &"  triangleSum over {vals.len} random configs:"
    for v in vals: echo &"    {v: .9f}"
    echo &"  mean = {mean:.3e}, max |Re| = {maxRe:.3e}, max |Im| = {maxIm:.3e}"
    # 64*nCells triangles, each ~ Re Tr(random SU(3))/3 with variance 1/(2 Nc^2)
    let sigma = sqrt(1.0/(2.0*nc*nc*float(64*hl.physVol)))
    echo &"  expected size ~ {sigma:.3e} (1 sigma)"
    ok(&"|mean| = {abs(mean):.3e} << 1", abs(mean) < 10.0*sigma)
    ok(&"every value within 5 sigma of 0 (max {maxRe:.2e})", maxRe < 5.0*sigma)
    # NOTE: the 32 apex triangles are enumerated with a *fixed* orientation, so
    # unlike the cubic plaquette sum there is no exact cancellation of the
    # imaginary part; it is only statistically zero, the same size as the real
    # part.  It *is* exactly zero for unit and pure-gauge configurations
    # (checked in tests 2 and 4).
    ok(&"|Im| is the same statistical size as |Re| ({maxIm:.2e})",
       maxIm < 5.0*sigma)

  test "6. HcGauge helpers: unit / warm / reunit / link accessor":
    let hl = newLayout(@[4, 4, 4, 4])
    let lo = hl
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
    for mu in 0..<nDim:
      if g.link(lkA, mu) == g.links[mu]: inc nMatch
      if g.link(lkB, mu) == g.links[nDim+mu]: inc nMatch
    for d in 0..<nDiag:
      if g.link(lkD, d) == g.links[2*nDim+d]: inc nMatch
    ok(&"link() accessor consistent with links ({nMatch}/24)", nMatch == 24)

  test "7. V = 1 layout and a short (L = 2) direction":
    # L = 2 makes the +delta shifts wrap around repeatedly; V = 1 takes the
    # non-vectorised path
    let geom2 = @[2, 4, 4, 6]
    let v1 = altLayout(newLayout(geom4, 1), geom4)
    let l2 = altLayout(newLayout(geom2), geom2)
    echo &"  V = 1: |triangleSum - ref| = {v1.dtri:.2e}, 16-shift errors {v1.ef:.1e} {v1.eb:.1e}"
    echo &"  L = 2: |triangleSum - ref| = {l2.dtri:.2e}, 16-shift errors {l2.ef:.1e} {l2.eb:.1e}"
    ok(&"V = 1: field and reference agree to 1e-12 (diff {v1.dtri:.2e})", v1.dtri < 1e-12)
    ok("V = 1: forward and backward 16-shifts exact", v1.ef == 0.0 and v1.eb == 0.0)
    ok(&"L = 2: field and reference agree to 1e-12 (diff {l2.dtri:.2e})", l2.dtri < 1e-12)
    ok("L = 2: forward and backward 16-shifts exact", l2.ef == 0.0 and l2.eb == 0.0)

qexFinalize()
