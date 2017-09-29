import math
import base, field, layout, maths

type
  RNG* = concept var r
    r.uniform
    r.gaussian
  RNGField* = concept var r
    r[0] is RNG

when defined(FUELCompat):
  # For maximal compatibility, see below.
  proc gaussian_call2(x: AsVarComplex, a,b:float) =
    x.re = a
    x.im = b
proc gaussian*(x: AsVarComplex, r: var RNG) =
  mixin gaussian
  when defined(FUELCompat):
    # This is how QLA does it for complex types (e.g. QLA_D3_V_veq_gaussian_S).
    # Technically which one in this call gets evaluated is undefined in C.
    # Let's hope if you use the same C compiler,
    # the evaluation order turns out to be the same.
    x.gaussian_call2(gaussian(r), gaussian(r))
  else:
    x.re = gaussian(r)
    x.im = gaussian(r)
proc gaussian*(x: AsVarVector, r: var RNG) =
  forO i, 0, x.len-1:
    gaussian(x[i], r)
proc gaussian*(x: AsVarMatrix, r: var RNG) =
  forO i, 0, x.nrows-1:
    forO j, 0, x.ncols-1:
      gaussian(x[i,j], r)
proc gaussian*(v: Field, r: RNGField) =
  for i in v.l.sites:
    gaussian(v{i}, r[i])

proc uniform*(x: AsVarComplex, r: var RNG) =
  mixin uniform
  x.re = uniform(r)
  x.im = uniform(r)
proc uniform*(x: AsVarVector, r: var RNG) =
  forO i, 0, x.len-1:
    uniform(x[i], r)
proc uniform*(x: AsVarMatrix, r: var RNG) =
  forO i, 0, x.nrows-1:
    forO j, 0, x.ncols-1:
      uniform(x[i,j], r)
proc uniform*(v: Field, r: var RNGField) =
  for i in v.l.sites:
    uniform(v{i}, r[i])

proc z4*(x: AsVarComplex, r: var RNG) =
  when defined(FUELCompat):
    x.gaussian r
    var n,o {.noinit.}: float
    n := x.re
    o := x.im
    if n >= 0:
      if o >= 0:
        x.re = 1.0
        x.im = 0.0
      else:
        x.re = 0.0
        x.im = 1.0
    else:
      if o >= 0:
        x.re = -1.0
        x.im = 0.0
      else:
        x.re = 0.0
        x.im = -1.0
  else:
    let n = r.uniform
    if n < 0.5:
      if n < 0.25:
        x.re = 1.0
        x.im = 0.0
      else:
        x.re = 0.0
        x.im = 1.0
    else:
      if n < 0.75:
        x.re = -1.0
        x.im = 0.0
      else:
        x.re = 0.0
        x.im = -1.0
proc z4*(x: AsVarVector, r: var RNG) =
  forO i, 0, x.len-1: x[i].z4 r
proc z4*(x: AsVarMatrix, r: var RNG) =
  forO i, 0, x.nrows-1:
    forO i, 0, x.ncols-1:
      x[i,j].z4 r
proc z4*(x: Field, r: var RNGField) =
  for i in x.l.sites:
    x{i}.z4 r[i]

proc z2*(x: AsVarComplex, r: var RNG) =
  when defined(FUELCompat):
    x.gaussian r
    var n {.noinit.}:float
    n := x.re
    if n >= 0: x := 1.0
    else: x := -1.0
  else:
    let n = r.uniform
    x.im = 0.0
    if n < 0.5: x.re = 1.0
    else: x.re = -1.0
proc z2*(x: AsVarVector, r: var RNG) =
  forO i, 0, x.len-1: x[i].z2 r
proc z2*(x: AsVarMatrix, r: var RNG) =
  forO i, 0, x.nrows-1:
    forO i, 0, x.ncols-1:
      x[i,j].z2 r
proc z2*(x: Field, r: var RNGField) =
  for i in x.l.sites:
    x{i}.z2 r[i]

proc u1*(x: AsVarComplex, r: var RNG) =
  when defined(FUELCompat):
    x.gaussian r
    let n = x.norm2
    if n == 0:
      x.re = 1.0
      x.im = 0.0
    else:
      let s = 1.0 / sqrt n
      x.re *= s
      x.im *= s
  else:
    let n = 2.0 * PI * r.uniform
    x.re = cos n
    x.im = sin n
proc u1*(x: AsVarVector, r: var RNG) =
  forO i, 0, x.len-1: x[i].u1 r
proc u1*(x: AsVarMatrix, r: var RNG) =
  forO i, 0, x.nrows-1:
    forO i, 0, x.ncols-1:
      x[i,j].u1 r
proc u1*(x: Field, r: var RNGField) =
  for i in x.l.sites:
    x{i}.u1 r[i]

proc newRNGField*(R:typedesc[RNG], lo:Layout,
                  seed:uint64 = uint64(17^7)):auto =
  var r:Field[1,R]
  r.new(lo.physGeom.newLayout 1)
  threads:
    for j in lo.sites:
      var l = lo.coords[lo.nDim-1][j].int
      for i in countdown(lo.nDim-2, 0):
        l = l * lo.physGeom[i].int + lo.coords[i][j].int
      r[j].seed(seed, l)
  return r
