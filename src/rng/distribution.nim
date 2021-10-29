import math
import base, field, layout, maths, maths/types, simd
import comms/qmp

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
proc gaussian*(x: var SomeNumber, r: var RNG) =
  mixin gaussian
  x = gaussian(r)
proc gaussian*(x: var Simd, r: var RNG) =  # FIXME to set all lanes
  mixin gaussian
  x[] := gaussian(r)
proc gaussian*(x: var AsComplex, r: var RNG) =  # FIXME to set all lanes
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
proc gaussian*[T](a: openArray[T], r: RNGField) =
  for i in 0..<a.len:
    gaussian(a[i], r)

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
    x{i}.z4 r{i}

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
    x{i}.z2 r{i}

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
    x{i}.u1 r{i}

proc vonMisesWithExp[D](rng:var RNG, lambda:D):auto =
  ## sample x ~ exp(lambda*cos(x))
  ## using exponential distribution for rejection sampling
  ## exp(lambda*(b-a*x))
  var
    a {.noinit.}:D
    b {.noinit.}:D
    x {.noinit.}:D
  if lambda>1.904538388056459:
    # optimal in the limit of large lambda
    a = 1.0/sqrt(lambda)
    b = a*arcsin(a)+sqrt(1-a*a);
  else:
    # keep the line above cos(x), b-a*%pi > cos(%pi) = -1
    a = 0.7246113537767085
    b = 1.276433705732662
  a *= lambda
  while true:
    let
      r = (2.0*rng.uniform-1.0)*expm1(-PI*a)
      u = rng.uniform
      acc =
        if r<0:
          let y = log1p(r)
          x = y/a
          exp(lambda*(cos(x)-b)-y)
        else:
          let y = -log1p(-r)
          x = y/a
          exp(lambda*(cos(x)-b)+y)
    if u < acc: return x

#[
proc vonMisesQOPQDP[D](rng:var RNG, g:D):auto =
  ## from QOPQDP
  const WENSLEY_CONST = 1.05110196582237  # a*asin(a)+sqrt(1-a*a),a:1/%pi
  if g<0.01:  # simple accept/reject
    let norm = exp( g )
    while true:
      let xr = rng.uniform
      var theta = 2*xr - 1.0
      theta *= PI
      let f = exp( g * cos( theta ) )
      let r = rng.uniform
      if f > r*norm:
        return theta
  else:  # Wensley linear filter
    let norm = exp( g*WENSLEY_CONST )
    while true:
      let xr = rng.uniform
      var theta =
        if xr<0.5:
          ln( 1 + 2*( exp( g ) - 1 )*xr ) / g - 1;
        else:
          1 - ln( 1 + 2*( exp( g ) - 1 )*( 1 - xr ) ) / g;
      theta *= PI;
      let f = exp( g*( cos( theta ) + abs( theta )/PI ) ) / norm;
      let r = rng.uniform
      if f > r:
        return theta
]#
#[
proc vonMisesWithWrappedCauchy[D](rng:var RNG, k:D):auto =
  ## Best, D., & Fisher, N. (1979).
  ## Efficient Simulation of the von Mises Distribution.
  ## Journal of the Royal Statistical Society.
  ## Series C (Applied Statistics), 28(2), 152-157.
  ## doi:10.2307/2346732
  let
    t = 1.0+sqrt(1.0+4.0*k*k)
    p = (t-sqrt(2.0*t))/(2.0*k)
    r = (1.0+p*p)/(2.0*p)
  var f {.noinit.}:D
  while true:
    let
      u1 = rng.uniform
      z = cos(PI*u1)
    f = (1.0+r*z)/(r+z)
    let
      c = k*(r-f)
      u2 = rng.uniform
    if c*(2.0-c)-u2>0 or ln(c/u2)+1.0-c>0:
      break
  let
    u3 = rng.uniform
    theta = arcCos(f)
  if u3<0.5: return -theta
  else: return theta
]#

proc vonMises*[D](rng:var RNG, lambda:D):auto =
  vonMisesWithExp(rng,lambda)

proc newRNGField*[R: RNG](lo: Layout, rng: typedesc[R],
                          s: uint64 = uint64(17^7)): Field[1,R] =
  ## The seed `s` is broadcasted from rank 0.
  var ss = s
  QMP_broadcast(ss.addr, sizeof(ss).csize_t)
  var r: Field[1,rng]
  when lo.V == 1:
    r.new(lo)
  else:
    echo "#newRNGField lo:"
    r.new(lo.physGeom.newLayout(1, lo.rankGeom))
  let t = r[0]  # Workaround Nim bug (Nim needs to see the type instantiated.)
  threads:
    for j in lo.sites:
      var l = lo.coords[lo.nDim-1][j].int
      for i in countdown(lo.nDim-2, 0):
        l = l * lo.physGeom[i].int + lo.coords[i][j].int
      seedIndep(r[j], ss, l)
  r
proc newRNGField*[R: RNG](rng: typedesc[R], lo: Layout,
                          s: uint64 = uint64(17^7)): Field[1,R] =
  lo.newRNGField(rng, s)
