## Brief: RNG field operations for MILC version 6 generator
## 
## Author: James C. Osborn
## Author: Xiao-Yong Jin
## Modified: Curtis Taylor Peterson <curtistaylorpetersonwork@gmail.com>

import milcrngv6
export milcrngv6

import base 
import field 
import layout 
import simd
import maths 
import maths/[types]
import comms/[qmp]
import base/[wrapperTypes]
import physics/[color]

import math

type
  MILCRNG* = concept var r
    r.uniform
    r.agaussian
  MILCRNGField* = concept r
    r[0] is MILCRNG

proc newRNGField*[R: MILCRNG](lo: Layout, rng: typedesc[R], s: uint64 = uint64(17^7)): Field[1,R] =
  ## The seed `s` is broadcasted from rank 0.
  var ss = s
  QMP_broadcast(ss.addr, sizeof(ss).csize_t)
  var r: Field[1,rng]
  when lo.V == 1: r.new(lo)
  else:
    echo "#newRNGField lo:"
    r.new(lo.physGeom.newLayout(1, lo.rankGeom))
  threads:
    # "defined(RandCoordOrder) or not defined(RandRawOrder)" = true by default
    for j in r.l.sites:
      var l = r.l.coords[r.l.nDim-1][j].int
      for i in countdown(r.l.nDim-2, 0):
        l = l * r.l.physGeom[i].int + r.l.coords[i][j].int
      seedIndep(r[j], ss, l)
  return r

proc newRNGField*[R: MILCRNG](rng: typedesc[R], lo: Layout, s: uint64 = uint64(17^7)): Field[1,R] =
  return lo.newRNGField(rng, s)

# "defined(RandCoordOrder) or not defined(RandRawOrder)" = true by default
template mapRngField*(fn: untyped, x: untyped, r: untyped) = 
  let nd = x.l.nDim
  var c = newSeq[int32](nd)
  for i in x.l.sites:
    x.l.coord(c, i)
    let j = r.l.rankIndex(c).index
    fn(x{i}, r{j})

proc uniform*(x: var AsNumber, r: var MILCRNG) =
  mixin uniform
  x := uniform(r)

proc uniform*(x: var AsComplex, r: var MILCRNG) =
  mixin uniform
  uniform(x.re, r)
  uniform(x.im, r)

proc uniform*(x: var AsVector, r: var MILCRNG) =
  forO i, 0, x.len-1: uniform(x[i], r)

proc uniform*(x: var AsMatrix, r: var MILCRNG) =
  forO i, 0, x.nrows-1:
    forO j, 0, x.ncols-1: uniform(x[i,j], r)

template uniform*(r: AsVar, x: untyped) =
  mixin uniform
  var t = r[]
  uniform(t, x)

proc uniform*(v: Field, r: MILCRNGField) = mapRngField(uniform, v, r)

template uniform*(x: var Color, r: var untyped) =
  uniform(x[], r)

proc agaussian*(x: var SomeNumber, r: var MILCRNG) =
  mixin agaussian
  x = agaussian(r)

proc agaussian*(x: var AsNumber, r: var MILCRNG) =
  mixin agaussian
  x := agaussian(r)

proc gaussian_call2(x: var AsComplex, a,b:float) = (x.re, x.im) = (a, b)
proc agaussian*(x: var AsComplex, r: var MILCRNG) =
  mixin agaussian
  # This is how QLA does it for complex types (e.g. QLA_D3_V_veq_gaussian_S).
  # Technically which one in this call gets evaluated is undefined in C.
  # Let's hope if you use the same C compiler,
  # the evaluation order turns out to be the same.
  when numNumbers(x.re) > 1:
    static: echo "agaussian for type ", typeof(x), " not implemented"
    {.error.}
  x.gaussian_call2(agaussian(r), agaussian(r))

proc agaussian*[T:array](x: MaskedObj[T], r: var MILCRNG) =
  for i in 0..<x.len: agaussian(x[i], r)

proc agaussian*(x: var array, r: var MILCRNG) =
  for i in 0..<x.len: agaussian(x[i], r)

proc agaussian*(x: var AsVector, r: var MILCRNG) =
  forO i, 0, x.len-1: agaussian(x[i], r)

proc agaussian*(x: var AsMatrix, r: var MILCRNG) =
  forO i, 0, getConst(x.nrows-1):
    forO j, 0, getConst(x.ncols-1): agaussian(x[i,j], r)

template agaussian*(r: AsVar, x: untyped) =
  mixin agaussian
  var t = r[]
  agaussian(t, x)

proc agaussian*(v: Field, r: MILCRNGField) = mapRngField(agaussian, v, r)

proc agaussian*[T](a: openArray[T], r: MILCRNGField) =
  for i in 0..<a.len: agaussian(a[i], r)

template agaussian*(x: Color, r: untyped) =
  agaussian(x[], r)

proc randTah3(m: var auto, s: var auto) =
  let s2 = 0.70710678118654752440;  # sqrt(1/2)
  let s3 = 0.57735026918962576450;  # sqrt(1/3)
  let r3 = s2 * agaussian(s)
  let r8 = s2 * s3 * agaussian(s)
  m[0,0].set 0, r8+r3
  m[1,1].set 0, r8-r3
  m[2,2].set 0, -2*r8
  let r01 = s2 * agaussian(s)
  let r02 = s2 * agaussian(s)
  let r12 = s2 * agaussian(s)
  let i01 = s2 * agaussian(s)
  let i02 = s2 * agaussian(s)
  let i12 = s2 * agaussian(s)
  m[0,1].set  r01, i01
  m[1,0].set -r01, i01
  m[0,2].set  r02, i02
  m[2,0].set -r02, i02
  m[1,2].set  r12, i12
  m[2,1].set -r12, i12

proc randomTAH*(x: Field, r: MILCRNGField) =
  when x[0].nrows == 3: mapRngField(randTah3, x, r)
  else:
    agaussian(x,r)
    x.projectTAH

proc randomU*(x: Field, r: var MILCRNGField) =
  agaussian(x,r)
  x.projectU

proc randomSU*(x: Field, r: var MILCRNGField) =
  agaussian(x,r)
  x.projectSU

proc warmSU*(x: Field, s: float, r: var MILCRNGField) =
  randomTAH(x,r)
  for i in x:
    let t = s * x[i]
    x[i] = exp(t)

template random*(x: var Color) =
  agaussian(x[], r)

proc random*[F:Field](g: openArray[F], r: var MILCRNGField) =
  for mu in g.low..g.high:
    when g[mu][0].nrows==1: randomU(g[mu], r)
    else: randomSU(g[mu], r)

proc warm*[F:Field](g: openArray[F], s: float, r: var MILCRNGField) =
  for mu in g.low..g.high:
    when g[mu][0].nrows==1:
      agaussian(g[mu],r)
      g[mu] := (1-s) + s*g[mu]
      g[mu].projectU
    else: warmSU(g[mu],s,r)

proc random*(g: array or seq) =
  var r = newRNGField(MilcRngv6, g[0].l)
  threads: random(g,r)