import math, unittest
import qex except epsilon
import algorithms/numdiff
import helpers
import ../[core, scalar, gauge]
import ../gauge/[types, field_ops, transport, stencil]

let grt = initGraphRuntime()
include gauge/gaugehelpers

qexInit()
echo "stencil ranks: ", nRanks
letParam:
  lat = latticeFromLocalLattice(@[4,4,8,8], nRanks)
let lo = lat.newLayout
doAssert lo.nDim >= 2, "stencil tests need at least two lattice dimensions"
var axes: seq[int]
for d, n in lo.rankGeom:
  if n > 1:
    axes.add d
let split = axes.len
for d, n in lo.rankGeom:
  if n == 1:
    axes.add d
let
  mu = axes[0]
  nu = axes[1]
echo "stencil axes: ", mu, ", ", nu, "; split dimensions: ", split
var
  rng = lo.newRNGField(Philox4x64, 1290911u64)
  g = lo.newGauge
  p = lo.newGauge
threads:
  g.random rng
  p.randomTAH rng
  for f in p:
    f *= 0.1
let
  gg = grt.toGvalue(g)
  gp = grt.toGvalue(p)

suite "graph remote stencil and reverse exchange":
  test "gather/scatter cross faces, corners, and two-step offsets":
    let f = linkField(gg, mu)
    let b = linkField(gp, nu)
    var faces, corners = false
    for off in [[1,0], [0,1], [1,-1], [-2,1]]:
      var sh = newSeq[int](lo.nDim)
      sh[mu] = off[0]
      sh[nu] = off[1]
      var remote = 0
      for d, n in sh:
        if n != 0 and lo.rankGeom[d] > 1:
          inc remote
      faces = faces or remote > 0
      corners = corners or remote > 1
      var refv = f
      for d, n in sh:
        for _ in 0..<abs(n):
          refv = shift(refv, d, if n > 0: -1 else: 1)
      norm2(gather(f, sh) - refv) :< 1e-18
      norm2(scatter(gather(f, sh), sh) - f) :< 1e-18
      (redot(gather(f, sh), b) - redot(f, scatter(b, sh))) :< 1e-9
    if nRanks > 1:
      check faces
    if split > 1:
      check corners

  test "path products and further pullbacks match unit-hop references":
    let
      a = mu+1
      b = nu+1
      paths = @[@[a,b,-a,-b], @[a,a,b,-a,-a,-b]]
    if nRanks > 1:
      check lo.rankGeom[mu] > 1
    if split > 1:
      check lo.rankGeom[nu] > 1
    proc loss(x: Ggauge, fused: bool): Gscalar =
      let ps = if fused: lineProducts(x, paths) else:
        @[wilsonLine(x, paths[0]), wilsonLine(x, paths[1])]
      # The upstream depends on x, exercising live seeds and reverse exchanges.
      redot(ps[0], linkField(x, nu)) + redot(ps[1], linkField(gp, mu))
    let
      sf = loss(gg, true)
      sr = loss(gg, false)
      df = grad(sf, gg)
      dr = grad(sr, gg)
      hf = grad(redot(df, gp), gg)
      hr = grad(redot(dr, gp), gg)
    for _ in 0..1:
      (sf - sr) :< 1e-8
      norm2(df - dr) :< 1e-18
      norm2(hf - hr) :< 1e-18
      gg.update p

qexFinalize()
