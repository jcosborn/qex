## Spinor-field linear algebra.
##
## `Spin` is a flat `seq[Spinor]` (doc/04 section 5); every routine below is
## allocation-free except the two constructors, and serial by design.
## Fields passed to a binary routine must have the same length; a mismatch
## shows up as an index error.
##
## `rng/threefry4x64` is imported directly: the `rng` umbrella module pulls in
## `field` and with it the hypercubic machinery we are avoiding.

import std/math
import types
import rng/threefry4x64
export types, threefry4x64

proc newSpin*(n: int): Spin =
  ## Zero-initialized field of `n` sites.
  newSeq[Spinor](n)

proc zero*(x: var Spin) =
  for i in 0..<x.len:
    for c in 0..1:
      x[i][c].re = 0.0
      x[i][c].im = 0.0

proc `:=`*(x: var Spin, y: Spin) =
  for i in 0..<y.len: x[i] = y[i]

proc axpy*(x: var Spin, a: float, y: Spin) =
  ## x += a*y
  for i in 0..<y.len:
    for c in 0..1:
      x[i][c].re += a*y[i][c].re
      x[i][c].im += a*y[i][c].im

proc axpy*(x: var Spin, a: Complex64, y: Spin) =
  ## x += a*y
  let ar = a.re
  let ai = a.im
  for i in 0..<y.len:
    for c in 0..1:
      let yr = y[i][c].re
      let yi = y[i][c].im
      x[i][c].re += ar*yr - ai*yi
      x[i][c].im += ar*yi + ai*yr

proc axpby*(x: var Spin, a: float, y: Spin, b: float) =
  ## x = a*y + b*x -- the CG search-direction update.
  for i in 0..<y.len:
    for c in 0..1:
      x[i][c].re = a*y[i][c].re + b*x[i][c].re
      x[i][c].im = a*y[i][c].im + b*x[i][c].im

proc scale*(x: var Spin, a: float) =
  for i in 0..<x.len:
    for c in 0..1:
      x[i][c].re *= a
      x[i][c].im *= a

proc dot*(x, y: Spin): Complex64 =
  ## sum_i conj(x_i) y_i
  var sr = 0.0
  var si = 0.0
  for i in 0..<x.len:
    for c in 0..1:
      let xr = x[i][c].re
      let xi = x[i][c].im
      let yr = y[i][c].re
      let yi = y[i][c].im
      sr += xr*yr + xi*yi
      si += xr*yi - xi*yr
  complex64(sr, si)

proc redot*(x, y: Spin): float =
  ## Re sum_i conj(x_i) y_i
  var s = 0.0
  for i in 0..<x.len:
    for c in 0..1:
      s += x[i][c].re*y[i][c].re + x[i][c].im*y[i][c].im
  s

proc norm2*(x: Spin): float =
  var s = 0.0
  for i in 0..<x.len:
    for c in 0..1:
      s += x[i][c].re*x[i][c].re + x[i][c].im*x[i][c].im
  s

proc gaussian*(x: var Spin, r: var Threefry4x64) =
  ## <|x_i|^2> = 1 per complex component, i.e. Re and Im are each N(0, 1/2).
  let s = 1.0/sqrt(2.0)
  for i in 0..<x.len:
    for c in 0..1:
      x[i][c].re = s*r.gaussian
      x[i][c].im = s*r.gaussian

proc pointSource*(n, site, comp: int): Spin =
  ## Unit source in spinor component `comp` at flat index `site`.
  result = newSpin(n)
  result[site][comp].re = 1.0
