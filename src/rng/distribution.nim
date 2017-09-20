import math
import base, field, layout, maths

type
  RNG* = concept var r
    r.uniform
    r.gaussian
  RNGField* = concept var r
    r[0] is RNG

proc gaussian*(x: AsVarComplex, r: var RNG) =
  mixin gaussian
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
  # This is different from the implemetation used in FUEL
  let n = r.uniform
  if n < 0.25:
    x.re = 1.0
    x.im = 0.0
  elif n < 0.5:
    x.re = 0.0
    x.im = 1.0
  elif n < 0.75:
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
  # This is different from the implemetation used in FUEL
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
  x.gaussian r
  let n = x.norm2
  if n == 0:
    x.re = 1.0
    x.im = 0.0
  else:
    let s = 1.0 / sqrt n
    x.re *= s
    x.im *= s
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
