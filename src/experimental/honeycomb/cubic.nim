## Cubic-lattice counterparts used by the comparison drivers and tests:
## repeated QEX stout steps and the tree-level clover Wilson operator.
##
##   D_c = D_W - (c_SW/2) sum_{a>b} gamma_a gamma_b Fhat_ab,   Fhat_ab = fmunu(g,1)[a][b] = +i a^2 F_ab T
## the same term and pair storage as hcwilson (r_w = 1 in physics/wilsonD).

import base, layout, field, maths
import physics/qcdTypes
import physics/wilsonD
import gauge
import gauge/stoutsmear
import hcgeom, hctopo, hcwilson

export stoutsmear, wilsonD

proc smearN*[G](s: var StoutSmear[G], g: G, gout: G, n: int) =
  ## n stout steps; n = 0 copies g to gout
  if n <= 0:
    threads:
      for mu in 0..<gout.len:
        gout[mu] := g[mu]
  else:
    s.smear(g, gout)
    for i in 1..<n:
      s.smear(gout, gout)

type
  CubicWilson*[W, MF] = ref object
    s*: W               ## physics/wilsonD operator (holds s.g)
    f*: array[6, MF]    ## Fhat_ab, pairIndex order
    cSW*: float

proc gaugeRefresh*(c: CubicWilson) =
  ## recompute the clover field from c.s.g; allocates (QEX fmunu)
  let fm = fmunu(c.s.g, 1)
  threads:
    for a in 1..<4:
      for b in 0..<a:
        c.f[pairIndex(a, b)] := fm[a][b]

proc newCubicWilson*[G](g: seq[G], cSW: float): auto =
  var s = newWilson(g)
  var c = CubicWilson[typeof(s), G](s: s, cSW: cSW)
  for p in 0..<6:
    c.f[p] = g[0].newOneOf
  c.gaugeRefresh
  c

proc newCubicWilson*[G, T](g: seq[G], v: T, cSW: float): auto =
  ## `v` is a prototype fermion field (any SIMD length)
  var s = newWilson(g, v)
  var c = CubicWilson[typeof(s), G](s: s, cSW: cSW)
  for p in 0..<6:
    c.f[p] = g[0].newOneOf
  c.gaugeRefresh
  c

proc applyClover(c: CubicWilson, r: auto, x: auto, cf: float) =
  var gm {.noinit.}: array[6, SpinMat]
  for p in 0..<6:
    gm[p] := cf*cloverGam[p]
  threads:
    for p in 0..<6:
      for e in r:
        r[e] += gm[p] * (c.f[p][e] * x[e])

proc D*(c: CubicWilson, r: var auto, x: auto, m: SomeNumber) =
  ## r must not alias x; opens its own `threads:`
  let rr = r
  threads:
    c.s.D(rr, x, m)
  if c.cSW != 0.0:
    applyClover(c, rr, x, -0.5*c.cSW)

proc Ddag*(c: CubicWilson, r: var auto, x: auto, m: SomeNumber) =
  let rr = r
  threads:
    c.s.Ddag(rr, x, m)
  if c.cSW != 0.0:
    applyClover(c, rr, x, -0.5*c.cSW)
