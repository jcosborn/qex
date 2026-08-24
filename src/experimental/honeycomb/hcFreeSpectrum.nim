## Free Wilson-Dirac spectrum on 16^4 -- reproduces slide 14 / paper Fig. 3.
##
##   export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
##   cd build_mac && make run src/experimental/honeycomb/hcFreeSpectrum.nim
##
## Optional arguments:  hcFreeSpectrum [ns] [nt] [per|anti]
## Default: 16^4, periodic in time (which is what makes the cubic spectrum
## reach exactly Re lambda = 8, as on the slide).
##
## Writes  doc/plots/freespec_cubic.dat, doc/plots/freespec_16cell.dat,
##         doc/plots/freespec.gnuplot  and renders  doc/plots/freespec.png.

import std/[math, complex, os, sets, strformat, strutils]
import hcfree

const srcDir = currentSourcePath().parentDir
let outDir = srcDir / "doc" / "plots"

proc dump(lat: HcLat; ns, nt: int; anti: bool; fn: string):
    tuple[n, nuniq: int, maxRe, maxIm: float] =
  ## Every eigenvalue of the free r=1, m=0 operator, deduplicated at 1e-6
  ## (the point group maps most momenta onto each other, so the raw list has
  ## a ~1000-fold redundancy; the scatter plot is identical either way).
  var seen = initHashSet[(int64, int64)]()
  var lines: seq[string]
  for p in momenta(lat, ns, nt, anti):
    let e = freeEigs(p, lat, 1.0, 0.0)
    for k in 0..1:
      result.n += 1
      result.maxRe = max(result.maxRe, e[k].re)
      result.maxIm = max(result.maxIm, abs(e[k].im))
      let key = ((e[k].re*1e6).round.int64, (e[k].im*1e6).round.int64)
      if not seen.containsOrIncl(key):
        lines.add &"{e[k].re:.10f} {e[k].im:.10f}"
  result.nuniq = lines.len
  var f = open(fn, fmWrite)
  f.write("# free Wilson-Dirac eigenvalues, r = 1, m = 0\n")
  f.write(&"# lattice = {lat}, {ns}^3 x {nt}, time BC = " &
          (if anti: "antiperiodic" else: "periodic") & "\n")
  f.write(&"# modes = {result.n}, distinct to 1e-6 = {result.nuniq}\n")
  f.write("# Re(lambda)  Im(lambda)\n")
  for l in lines: f.write(l & "\n")
  f.close()

when isMainModule:
  var
    ns = 16
    nt = 16
    anti = false
  let a = commandLineParams()
  if a.len > 0: ns = parseInt(a[0])
  if a.len > 1: nt = parseInt(a[1])
  if a.len > 2: anti = a[2].startsWith("a")
  createDir(outDir)

  echo &"free Wilson-Dirac spectrum, r = 1, m = 0, {ns}^3 x {nt}, " &
       "time BC = " & (if anti: "antiperiodic" else: "periodic")
  var files: array[2, string]
  for i, lat in [lCubic, lHoneycomb]:
    let fn = outDir / (if lat == lCubic: "freespec_cubic.dat"
                       else: "freespec_16cell.dat")
    files[i] = fn
    let r = dump(lat, ns, nt, anti, fn)
    echo &"  {lat:8}: modes = {r.n:8}  distinct = {r.nuniq:7}  " &
         &"max Re = {r.maxRe:.10f}  max |Im| = {r.maxIm:.10f}"
    echo &"            -> {fn}"

  let gp = outDir / "freespec.gnuplot"
  var f = open(gp, fmWrite)
  f.write("""
# Free Wilson-Dirac spectrum on 16^4 (r = 1, m = 0)
# reproduces slide 14 of Katz-Nogradi, "Lattice QCD on the 16-cell honeycomb"
set terminal pngcairo size 900,700 font "Helvetica,18"
set output "freespec.png"
set xlabel "Re({/Symbol l})"
set ylabel "Im({/Symbol l})"
set xrange [0:9]
set yrange [-3:3]
set xtics 1
set ytics 1
set key top left
set grid
plot "freespec_cubic.dat"  u 1:2 w p pt 7 ps 0.35 lc rgb "#d00000" t "cubic", \
     "freespec_16cell.dat" u 1:2 w p pt 7 ps 0.35 lc rgb "#0040c0" t "16-cell"
""")
  f.close()
  echo &"            -> {gp}"
  let rc = execShellCmd(&"cd {outDir.quoteShell} && gnuplot freespec.gnuplot")
  if rc == 0: echo "            -> " & (outDir / "freespec.png")
  else: echo "  gnuplot failed (rc = ", rc, "); run it by hand in ", outDir
