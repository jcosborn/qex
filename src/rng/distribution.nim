import math
import base, field, layout, maths

type
  RNG* = concept var r
    r.uniform
    r.gaussian
  RNGField* = concept var r
    r[0] is RNG

proc gaussian*(x: var AsComplex, r: var RNG) =
  mixin gaussian
  x.re = gaussian(r)
  x.im = gaussian(r)
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
proc gaussian*[T: RNGField](v: Field, r: var T) =
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

proc newRNGField*[R: RNG](rng: typedesc[R], lo: Layout,
                          seed: uint64 = uint64(17^7)): Field[1,R] =
  var r: Field[1,R]
  r.new(lo.physGeom.newLayout 1)
  threads:
    for j in lo.sites:
      var l = lo.coords[lo.nDim-1][j].int
      for i in countdown(lo.nDim-2, 0):
        l = l * lo.physGeom[i].int + lo.coords[i][j].int
      r[j].seed(seed, l)
  r
