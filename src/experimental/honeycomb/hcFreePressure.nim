## Free lattice pressure  O = (p(T) - p(0))/T^4  per fermionic degree of
## freedom -- reproduces slides 12-13 / paper Fig. 4.
##
##   export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
##   cd build_mac && make run src/experimental/honeycomb/hcFreePressure.nim
##
## Optional arguments:  hcFreePressure [ntMax] [maxGrid]
##
## Method
## ------
## `p = (T/V) ln Z`, `ln Z = sum_p ln det D(p)`, `det D(p) = (M^2+|K|^2)^2`
## for `r = 1, m = 0`.  At fixed spatial momentum `q`, `M^2+|K|^2` is a
## polynomial `P(x)` in `x = cos theta` with `theta = p_3` (cubic) or
## `p_3/2` (16-cell, whose BZ is doubled -- the doubling is put in the time
## direction here, see `momentaAlt`).  Writing `P(x) = L prod_i (x - x_i)`,
## `x_i = cosh E_i`, and using
##
##   prod_{n=0}^{N'-1} ( cos theta_n - x ) = 2^{-N'} ( w^{N'/2} + w^{-N'/2} )^2
##   int_0^{2pi} (dtheta/2pi) ln( cos theta - x ) = -ln 2 + ln w,   w = e^E
##
## at the antiperiodic Matsubara points `theta_n = 2 pi (n+1/2)/N'` gives the
## *exact*, singularity-free result
##
##   O(N_t) = N_t^3 * < sum_i ln | 1 + v_i(q)^{N'} | >_q ,    v_i = e^{-E_i}
##
## with `N' = N_t` (cubic) / `2 N_t` (16-cell) and `< >_q` the average
## `int d^3q/(2pi)^3` over `[0,2pi)^3`.  The `ln p^2` singularity of the
## `T = 0` term has cancelled analytically; what is left is a smooth periodic
## integrand whose only non-analyticity is the `|q|` cone at `q = 0`.  A
## uniform grid therefore converges as `a/N_g^4 + b/N_g^6 + ...` (the odd part
## of `ln(1+e^{-x})` is exactly `-x/2`, so `|q|`, `|q|^3`, ... are the only
## singular terms), which is removed by Richardson extrapolation.
##
## Continuum value per fermionic dof: `7 pi^2/720`.

import std/[algorithm, math, complex, os, strformat, strutils, times]
import hcfree

const srcDir = currentSourcePath().parentDir
let outDir = srcDir / "doc" / "plots"

# ---------------------------------------------------------------------------
# the 3d Brillouin-zone average
# ---------------------------------------------------------------------------

proc gridO*(lat: HcLat; ng: int; ntMin, ntMax: int): seq[float] =
  ## `O(N_t) = N_t^3 <h>_q` on a uniform `ng^3` grid, using the 48-fold cubic
  ## symmetry of the integrand (reflections `q_j -> -q_j` and permutations).
  doAssert ng mod 2 == 0
  let m = ng div 2
  var acc = newSeq[float](ntMax + 1)
  var cmp = newSeq[float](ntMax + 1)      # Kahan compensation
  var wtot = 0'i64
  var q: array[3, float]
  let dq = 2.0*PI/ng.float
  for n0 in 0..m:
    q[0] = dq*n0.float
    let w0 = if n0 == 0 or n0 == m: 1 else: 2
    for n1 in n0..m:
      q[1] = dq*n1.float
      let w1 = w0*(if n1 == 0 or n1 == m: 1 else: 2)
      for n2 in n1..m:
        q[2] = dq*n2.float
        var w = w1*(if n2 == 0 or n2 == m: 1 else: 2)
        if n0 == n2: discard                       # all three equal
        elif n0 == n1 or n1 == n2: w *= 3
        else: w *= 6
        wtot += w
        let wf = w.float
        let dbl = lat == lHoneycomb
        for v in modeV(q, lat):
          # z = v^{N'} for N_t = ntMin .. ntMax, built incrementally
          let step = if dbl: v*v else: v      # v^{N'} advances by this per N_t
          var z = complex64(1.0, 0.0)
          for i in 1..ntMin: z = z*step
          for nt in ntMin..ntMax:
            if nt > ntMin: z = z*step
            let m2 = z.re*z.re + z.im*z.im
            if m2 < 1e-40: break                   # and smaller for larger nt
            let t = wf*(if m2 < 1e-20: z.re
                        else: 0.5*ln((1.0 + z.re)^2 + z.im^2))
            let y = t - cmp[nt]
            let s = acc[nt] + y
            cmp[nt] = (s - acc[nt]) - y
            acc[nt] = s
  doAssert wtot == ng.int64*ng.int64*ng.int64, $wtot
  result = newSeq[float](ntMax + 1)
  let vol = ng.float*ng.float*ng.float
  for nt in ntMin..ntMax:
    result[nt] = nt.float^3*acc[nt]/vol

proc gridOBrute(lat: HcLat; ng: int; ntMin, ntMax: int): seq[float] =
  ## unreduced reference for `gridO` (used only as a self-check)
  var acc = newSeq[float](ntMax + 1)
  var q: array[3, float]
  let dq = 2.0*PI/ng.float
  for n0 in 0..<ng:
    q[0] = dq*n0.float
    for n1 in 0..<ng:
      q[1] = dq*n1.float
      for n2 in 0..<ng:
        q[2] = dq*n2.float
        for nt in ntMin..ntMax:
          acc[nt] += pressureIntegrand(q, lat, nt)
  result = newSeq[float](ntMax + 1)
  for nt in ntMin..ntMax:
    result[nt] = nt.float^3*acc[nt]/(ng.float^3)

# ---------------------------------------------------------------------------
# extrapolation helpers
# ---------------------------------------------------------------------------

proc richardson(x, y: openArray[float]; p0, dp: float): seq[seq[float]] =
  ## Richardson table for `y_i = y_exact + a_1 x_i^-p0 + a_2 x_i^-(p0+dp) + ...`
  ## Column `k` has eliminated the first `k` error terms.
  result.add @y
  var k = 0
  while result[k].len > 1:
    let p = p0 + dp*k.float
    var col: seq[float]
    for i in 1..<result[k].len:
      let
        a = pow(x[i + k], p)
        b = pow(x[i + k - 1], p)
      col.add (a*result[k][i] - b*result[k][i-1])/(a - b)
    result.add col
    inc k

proc neville(t, y: openArray[float]): float =
  ## polynomial extrapolation of `y(t)` to `t = 0`
  var c = @y
  for k in 1..t.high:
    for i in countdown(t.high, k):
      c[i] = (c[i]*(0.0 - t[i-k]) - c[i-1]*(0.0 - t[i]))/(t[i] - t[i-k])
  c[^1]

proc lstsq(a0: seq[seq[float]]; b: seq[float]): seq[float] =
  ## least squares via normal equations + Gaussian elimination.  The power
  ## basis is horribly scaled, so normalise the columns first.
  let n = a0[0].len
  var sc = newSeq[float](n)
  for j in 0..<n:
    var s = 0.0
    for k in 0..<a0.len: s += a0[k][j]*a0[k][j]
    sc[j] = if s > 0.0: 1.0/sqrt(s) else: 1.0
  var a = a0
  for k in 0..<a.len:
    for j in 0..<n: a[k][j] = a0[k][j]*sc[j]
  var m = newSeq[seq[float]](n)
  for i in 0..<n:
    m[i] = newSeq[float](n + 1)
    for j in 0..<n:
      var s = 0.0
      for k in 0..<a.len: s += a[k][i]*a[k][j]
      m[i][j] = s
    var s = 0.0
    for k in 0..<a.len: s += a[k][i]*b[k]
    m[i][n] = s
  for c in 0..<n:
    var pv = c
    for r in c+1..<n:
      if abs(m[r][c]) > abs(m[pv][c]): pv = r
    swap(m[pv], m[c])
    for r in 0..<n:
      if r == c: continue
      let f = m[r][c]/m[c][c]
      for j in c..n: m[r][j] -= f*m[c][j]
  result = newSeq[float](n)
  for i in 0..<n: result[i] = sc[i]*m[i][n]/m[i][i]

# ---------------------------------------------------------------------------

const
  c2Cubic = 248.0/147.0
  c4Cubic = 635.0/147.0
  c4Hc = 127.0/980.0
  c6Hc = 73.0/4158.0

proc seriesCubic(nt: float): float =
  1.0 + c2Cubic*PI^2/nt^2 + c4Cubic*PI^4/nt^4
proc seriesHc(nt: float): float =
  1.0 + c4Hc*PI^4/nt^4 + c6Hc*PI^6/nt^6

when isMainModule:
  var
    ntMax = 32
    maxGrid = 1024
  let av = commandLineParams()
  if av.len > 0: ntMax = parseInt(av[0])
  if av.len > 1: maxGrid = parseInt(av[1])
  const ntMin = 3
  # Two *geometric* grid families (ratio 2).  The Richardson recursion below
  # eliminates 1/ng^4, 1/ng^6, ... exactly only when the ratio ng_i/ng_{i-1} is
  # constant, so each family is extrapolated on its own; the difference between
  # the two independent answers is the error estimate.
  var famA, famB: seq[int]
  var ng = 32
  while ng <= maxGrid:
    famA.add ng
    ng *= 2
  ng = 48
  while ng <= maxGrid:
    famB.add ng
    ng *= 2
  var grids = famA & famB
  sort(grids)
  createDir(outDir)

  echo "self-check: symmetry-reduced grid sum vs brute force (ng = 12)"
  for lat in [lCubic, lHoneycomb]:
    let
      a = gridO(lat, 12, ntMin, 6)
      b = gridOBrute(lat, 12, ntMin, 6)
    var w = 0.0
    for nt in ntMin..6: w = max(w, abs(a[nt] - b[nt]))
    echo &"  {lat:8}: max |reduced - brute| = {w:.3e}"
    doAssert w < 1e-13

  # ---- the grids ----------------------------------------------------------
  var raw: array[HcLat, seq[seq[float]]]      # raw[lat][gridIndex][nt]
  for lat in [lCubic, lHoneycomb]:
    for g in grids:
      let t0 = epochTime()
      raw[lat].add gridO(lat, g, ntMin, ntMax)
      echo &"  {lat:8} ng = {g:5}  ({epochTime()-t0:6.2f} s)"

  # ---- Richardson in 1/ng^4, 1/ng^6, ... ----------------------------------
  var final, ferr: array[HcLat, seq[float]]
  for lat in [lCubic, lHoneycomb]:
    final[lat] = newSeq[float](ntMax + 1)
    ferr[lat] = newSeq[float](ntMax + 1)
  var convLines: seq[string]
  echo ""
  echo "grid convergence of O/O_cont  (rows: raw ng, then Richardson columns)"
  proc extrapFamily(lat: HcLat; fam: seq[int]; nt: int): seq[seq[float]] =
    var xg, y: seq[float]
    for g in fam:
      xg.add g.float
      y.add raw[lat][grids.find(g)][nt]
    richardson(xg, y, 4.0, 2.0)
  for lat in [lCubic, lHoneycomb]:
    for nt in ntMin..ntMax:
      let
        ta = extrapFamily(lat, famA, nt)
        tb = extrapFamily(lat, famB, nt)
      final[lat][nt] = ta[^1][0]
      ferr[lat][nt] = abs(ta[^1][0] - tb[^1][0])       # independent families
      for k in max(1, ta.high-2)..ta.high:             # and successive orders
        ferr[lat][nt] = max(ferr[lat][nt], abs(ta[k][^1] - final[lat][nt]))
      var s = &"{lat:8} nt={nt:3}  raw:"
      for i in 0..grids.high:
        s.add &"  {raw[lat][i][nt]/contPressure:.9f}"
      convLines.add s
      var s2 = &"{lat:8} nt={nt:3}  extrapolants A:"
      for k in 1..ta.high: s2.add &"  {ta[k][^1]/contPressure:.9f}"
      s2.add &"   | B:"
      for k in 1..tb.high: s2.add &"  {tb[k][^1]/contPressure:.9f}"
      convLines.add s2
      if nt in [4, 8, 12, 16, 20, 24, 32]:
        echo "  ", s
        echo "  ", s2

  # ---- the table ----------------------------------------------------------
  echo ""
  echo "  Nt     cubic O/O_cont   series      16-cell O/O_cont  series" &
       "      (grid err)"
  var f = open(outDir / "pressure.dat", fmWrite)
  f.write("# free lattice pressure  O = (p(T)-p(0))/T^4 per fermionic dof\n")
  f.write(&"# O_cont = 7 pi^2/720 = {contPressure:.12f}\n")
  f.write("# geometric grid families A = " & $famA & " and B = " & $famB & ",\n")
  f.write("# each Richardson-extrapolated in 1/ng^4, 1/ng^6, ...; the value is A,\n")
  f.write("# the |A-B| difference is the error: 1e-13..1e-10 for Nt <= 26,\n")
  f.write("# growing to 1e-9 at Nt = 32.  See pressure_conv.dat.\n")
  f.write("# Nt   O_cubic/O_cont   O_16cell/O_cont   series_cubic   series_16cell\n")
  for nt in ntMin..ntMax:
    let
      rc = final[lCubic][nt]/contPressure
      rh = final[lHoneycomb][nt]/contPressure
    f.write(&"{nt} {rc:.10f} {rh:.10f} {seriesCubic(nt.float):.10f} " &
            &"{seriesHc(nt.float):.10f}\n")
    echo &"  {nt:3}   {rc:14.8f}  {seriesCubic(nt.float):10.6f}   " &
         &"{rh:14.8f}  {seriesHc(nt.float):10.6f}   " &
         &"{ferr[lCubic][nt]/contPressure:.1e} {ferr[lHoneycomb][nt]/contPressure:.1e}"
  f.close()

  var fc = open(outDir / "pressure_conv.dat", fmWrite)
  fc.write("# raw grid values and Richardson extrapolants of O/O_cont\n")
  fc.write("# grids: " & $grids & "\n")
  for l in convLines: fc.write(l & "\n")
  fc.close()

  # analytic series on a fine Nt grid, for the plot
  var fs = open(outDir / "pressure_series.dat", fmWrite)
  fs.write("# Nt  series_cubic  series_16cell (truncated series of slide 12)\n")
  var x = 3.0
  while x <= 20.0001:
    fs.write(&"{x:.3f} {seriesCubic(x):.8f} {seriesHc(x):.8f}\n")
    x += 0.05
  fs.close()

  # ---- the asymptotic coefficients ---------------------------------------
  var dev: array[HcLat, seq[float]]           # O/O_cont - 1
  for lat in [lCubic, lHoneycomb]:
    dev[lat] = newSeq[float](ntMax + 1)
    for nt in ntMin..ntMax: dev[lat][nt] = final[lat][nt]/contPressure - 1.0

  echo ""
  echo "A. 5-node Neville extrapolation in t = 1/Nt^2 to t = 0, over windows"
  echo "   ending at Nt = nhi.  Large nhi reduces the truncation error but the"
  echo "   signal shrinks like Nt^-2k, so the numerical noise eventually wins."
  var est: array[HcLat, array[1..3, seq[float]]]   # collected estimates of c_2k
  proc cascade(lat: HcLat; kmin, nhiBest: int; want: seq[float]) =
    var d = dev[lat]
    for j in 0..<want.len:
      let k = kmin + j
      var trail: seq[string]
      var best = 0.0
      for nhi in [12, 16, 20, 24, 28, 32]:
        if nhi > ntMax: continue
        var t, y: seq[float]
        for nt in nhi-4..nhi:
          t.add 1.0/(nt.float*nt.float)
          y.add d[nt]*pow(nt.float*nt.float/(PI*PI), k.float)
        let v = neville(t, y)
        trail.add &"{v:11.7f}"
        # the 16-cell signal falls like Nt^-4, so its large-nhi Neville values
        # are noise-limited; only the cubic ones are collected for the summary
        if nhi >= 24 and lat == lCubic: est[lat][k].add v
        if nhi == nhiBest: best = v
      echo &"  {lat:8} c_{2*k} = {best:12.8f}  (target {want[j]:.8f}, " &
           &"rel {abs(best/want[j]-1.0)*100:7.3f}%)"
      echo &"                nhi = 12,16,20,24,28,32: " & trail.join(" ")
      for nt in ntMin..ntMax: d[nt] -= best*pow(PI*PI/(nt.float*nt.float), k.float)
  cascade(lCubic, 1, 32, @[c2Cubic, c4Cubic])
  cascade(lHoneycomb, 2, 20, @[c4Hc, c6Hc])

  echo ""
  echo "B. least-squares fits of O/O_cont - 1 in the basis (pi^2/Nt^2)^k"
  proc fit(lat: HcLat; kmin, nterm, nlo: int; verbose = true): seq[float] =
    var a: seq[seq[float]]
    var b: seq[float]
    for nt in nlo..ntMax:
      var row: seq[float]
      for k in kmin..<kmin+nterm: row.add pow(PI*PI/(nt.float*nt.float), k.float)
      a.add row
      b.add dev[lat][nt]
    result = lstsq(a, b)
    if verbose:
      var s = &"  {lat:8} Nt={nlo:3}..{ntMax}, k={kmin}..{kmin+nterm-1}:"
      for k in 0..<nterm: s.add &"  c_{2*(kmin+k)} = {result[k]:11.7f}"
      echo s
      if nlo == 12 and nterm >= 4:              # the well-conditioned window
        for k in 0..min(1, nterm-1):
          if kmin + k <= 3: est[lat][kmin+k].add result[k]
  echo "  cubic (targets  c_2 = 1.6870748, c_4 = 4.3197279):"
  echo "   -- literal 2-term fit at large Nt (contaminated by the c_6 term):"
  for nlo in [16, 20, 24]: discard fit(lCubic, 1, 2, nlo)
  echo "   -- multi-term fits:"
  for nterm in [4, 5, 6]:
    for nlo in [6, 8, 12]: discard fit(lCubic, 1, nterm, nlo)
  echo "  16-cell (targets  c_4 = 0.1295918, c_6 = 0.0175565):"
  echo "   -- literal 2-term fit at large Nt:"
  for nlo in [10, 14, 18]: discard fit(lHoneycomb, 2, 2, nlo)
  echo "   -- multi-term fits:"
  for nterm in [3, 4, 5]:
    for nlo in [6, 8, 12]: discard fit(lHoneycomb, 2, nterm, nlo)
  echo "  16-cell with a 1/Nt^2 term added -- c_2 must come out consistent with 0:"
  for nterm in [4, 5, 6]:
    for nlo in [6, 8, 12]: discard fit(lHoneycomb, 1, nterm, nlo)

  echo ""
  echo "C. residual test with the *published* coefficients (the sharp check):"
  echo "   R(Nt) = [O/O_cont - 1 - published terms] * Nt^p / pi^p must tend to a"
  echo "   finite constant (the next coefficient); if a published coefficient were"
  echo "   wrong, R would diverge as a power of Nt."
  echo "   Nt      cubic R (p=6)      16-cell R (p=8)"
  for nt in ntMin..ntMax:
    let x = PI*PI/(nt.float*nt.float)
    let rc = (dev[lCubic][nt] - c2Cubic*x - c4Cubic*x*x)/(x*x*x)
    let rh = (dev[lHoneycomb][nt] - c4Hc*x*x - c6Hc*x*x*x)/(x*x*x*x)
    if nt >= 6: echo &"  {nt:3}   {rc:16.7f}   {rh:16.7f}"

  echo ""
  echo "SUMMARY -- asymptotic coefficients, spread over the stable estimators"
  echo "          of A (Neville, nhi >= 24) and B (least squares, Nt = 12..32)"
  proc summarize(lat: HcLat; k: int; want: float) =
    let v = est[lat][k]
    if v.len == 0: return
    var lo = v[0]
    var hi = v[0]
    var mean = 0.0
    for x in v:
      lo = min(lo, x); hi = max(hi, x); mean += x/v.len.float
    echo &"  {lat:8} c_{2*k} = {mean:.7f} +- {0.5*(hi-lo):.7f}   " &
         &"target {want:.7f}   deviation {abs(mean/want - 1.0)*100:.3f}%   " &
         &"[{v.len} estimators, {lo:.7f} .. {hi:.7f}]"
  summarize(lCubic, 1, c2Cubic)
  summarize(lCubic, 2, c4Cubic)
  summarize(lHoneycomb, 2, c4Hc)
  summarize(lHoneycomb, 3, c6Hc)
  block:
    var lo = 1e300
    var hi = -1e300
    for nterm in [4, 5, 6]:
      let c = fit(lHoneycomb, 1, nterm, 12, verbose = false)
      lo = min(lo, c[0]); hi = max(hi, c[0])
    echo &"  16cell   c_2  = {0.5*(lo+hi):.3e} +- {0.5*(hi-lo):.3e}   " &
         "target EXACTLY ZERO (slide 12: leading 16-cell correction is O(a^4))"

  # ---- plot ---------------------------------------------------------------
  let gp = outDir / "pressure.gnuplot"
  var g = open(gp, fmWrite)
  g.write("""
# Free lattice pressure per fermionic dof, p/p_cont vs N_t
# reproduces slide 13 / Fig. 4 of Katz-Nogradi, "Lattice QCD on the 16-cell honeycomb"
set terminal pngcairo size 900,700 font "Helvetica,18"
set output "pressure.png"
set xlabel "N_t"
set ylabel "p / p_{cont}"
set xrange [3:20]
set yrange [1:4]
set xtics 2
set ytics 0.5
set key top right
set grid
plot "pressure_series.dat" u 1:2 w l lw 2 lc rgb "#d00000" notitle, \
     "pressure_series.dat" u 1:3 w l lw 2 lc rgb "#0040c0" notitle, \
     "pressure.dat" u 1:2 w p pt 5 ps 1.4 lc rgb "#d00000" t "cubic", \
     "pressure.dat" u 1:3 w p pt 7 ps 1.4 lc rgb "#0040c0" t "16-cell"
""")
  g.close()
  let rc = execShellCmd(&"cd {outDir.quoteShell} && gnuplot pressure.gnuplot")
  echo ""
  echo "wrote ", outDir / "pressure.dat", ", pressure_conv.dat, ",
       "pressure_series.dat, pressure.gnuplot"
  if rc == 0: echo "wrote ", outDir / "pressure.png"
  else: echo "gnuplot failed (rc = ", rc, ")"
