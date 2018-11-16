import math
import base, field, layout, maths, maths/types

type
  RNG* = concept var r
    r.uniform
    r.gaussian
  RNGField* = concept r
    r[0] is RNG

when defined(FUELCompat):
  # For maximal compatibility, see below.
  proc gaussian_call2(x: var AsComplex, a,b:float) =
    x.re = a
    x.im = b
proc gaussian*(x: var AsComplex, r: var RNG) =
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
proc gaussian*[T:array](x: MaskedObj[T], r: var RNG) =
  for i in 0..<x.len:
    gaussian(x[i], r)
proc gaussian*(x: var array, r: var RNG) =
  for i in 0..<x.len:
    gaussian(x[i], r)
proc gaussian*(x: var AsVector, r: var RNG) =
  forO i, 0, x.len-1:
    gaussian(x[i], r)
proc gaussian*(x: var AsMatrix, r: var RNG) =
  forO i, 0, getConst(x.nrows-1):
    forO j, 0, getConst(x.ncols-1):
      gaussian(x[i,j], r)
template gaussian*(r: AsVar, x: untyped) =
  mixin gaussian
  var t = r[]
  gaussian(t, x)
proc gaussian*(v: Field, r: RNGField) =
  for i in v.l.sites:
    gaussian(v{i}, r{i})

proc uniform*(x: var AsComplex, r: var RNG) =
  mixin uniform
  x.re = uniform(r)
  x.im = uniform(r)
proc uniform*(x: var AsVector, r: var RNG) =
  forO i, 0, x.len-1:
    uniform(x[i], r)
proc uniform*(x: var AsMatrix, r: var RNG) =
  forO i, 0, x.nrows-1:
    forO j, 0, x.ncols-1:
      uniform(x[i,j], r)
template uniform*(r: AsVar, x: untyped) =
  mixin uniform
  var t = r[]
  uniform(t, x)
proc uniform*(v: Field, r: RNGField) =
  for i in v.l.sites:
    uniform(v{i}, r{i}[])

proc z4*(x: var AsComplex, r: var RNG) =
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
proc z4*(x: var AsVector, r: var RNG) =
  forO i, 0, x.len-1: x[i].z4 r
proc z4*(x: var AsMatrix, r: var RNG) =
  forO i, 0, x.nrows-1:
    forO j, 0, x.ncols-1:
      x[i,j].z4 r
template z4*(r: AsVar, x: untyped) =
  mixin z4
  var t = r[]
  z4(t, x)
proc z4*(x: Field, r: RNGField) =
  for i in x.l.sites:
    x{i}.z4 r{i}[]

proc z2*(x: var AsComplex, r: var RNG) =
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
proc z2*(x: var AsVector, r: var RNG) =
  forO i, 0, x.len-1: x[i].z2 r
proc z2*(x: var AsMatrix, r: var RNG) =
  forO i, 0, x.nrows-1:
    forO j, 0, x.ncols-1:
      x[i,j].z2 r
template z2*(r: AsVar, x: untyped) =
  mixin z2
  var t = r[]
  z2(t, x)
proc z2*(x: Field, r: RNGField) =
  for i in x.l.sites:
    x{i}.z2 r{i}[]

proc u1*(x: var AsComplex, r: var RNG) =
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
    let n = 2.0 * PI * r.uniform.float
    x.re = cos n
    x.im = sin n
proc u1*(x: var AsVector, r: var RNG) =
  forO i, 0, x.len-1: x[i].u1 r
proc u1*(x: var AsMatrix, r: var RNG) =
  forO i, 0, x.nrows-1:
    forO j, 0, x.ncols-1:
      x[i,j].u1 r
template u1*(r: AsVar, x: untyped) =
  mixin u1
  var t = r[]
  u1(t, x)
proc u1*(x: Field, r: RNGField) =
  for i in x.l.sites:
    x{i}.u1 r{i}[]

proc newRNGField*[R: RNG](lo: Layout, rng: typedesc[R],
                          s: uint64 = uint64(17^7)): Field[1,R] =
  var r: Field[1,rng]
  r.new(lo.physGeom.newLayout 1)
  let t = r[0]  # Workaround Nim bug (Nim needs to see the type instantiated.)
  threads:
    for j in lo.sites:
      var l = lo.coords[lo.nDim-1][j].int
      for i in countdown(lo.nDim-2, 0):
        l = l * lo.physGeom[i].int + lo.coords[i][j].int
      seed(r[j], s, l)
  r
proc newRNGField*[R: RNG](rng: typedesc[R], lo: Layout,
                          s: uint64 = uint64(17^7)): Field[1,R] =
  lo.newRNGField(rng, s)
