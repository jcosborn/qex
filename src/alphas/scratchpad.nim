import qex
import std/[sequtils]
import gauge/[gaugefix, gaugeUtils]

qexInit()
let 
  geom = [8, 8, 8, 16]
  lo = geom.newLayout()
let nd = geom.len
var 
  f = lo.newGauge()
  u = lo.newGauge()
  m = lo.ColorMatrix()

# test gauge invariance of polyakov loop
let 
  gfstop = 1e-8
  gforf = 1.75
u.random()
threads: m := 1
m.getGaugeFixTransform(u, @[0, 1, 2], gfstop, gforf, verb = 0)
f.gaugeTransform(u, m)
for l in 0..<nd: 
  var psum, fsum, absdiff: float
  psum = u.wline(repeat(l+1, geom[l])).re()
  fsum = f.wline(repeat(l+1, geom[l])).re()
  absdiff = abs(psum - fsum)/abs(psum)
  echo "[" & $l & "]" & "poly: ", psum
  echo "[" & $l & "]" & "poly fixed: ", fsum
  echo "[" & $l & "]" & "abs diff: ", absdiff

qexFinalize()