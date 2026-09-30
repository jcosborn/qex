import testutils
import qex, physics/qcdTypes

# RNGFieldV: site i of its layout holds and draws with the generator that
# newRNGField gives the site of the same coordinates.  Fields of the layout of
# the RNGFieldV and of a one-lane layout draw the numbers of newRNGField, the
# generators exactly, the gaussians to rounding: the two layouts compute them
# in different code, which -ffast-math may contract or vectorize differently.

qexInit()
let lat = latticeFromLocalLattice([4,4,4,8], nRanks)
let lo = lat.newLayout
let l1 = lo.physGeom.newLayout(1, lo.rankGeom)
const seed = 987654321

proc near(x, y: Field; eps: float): bool =
  ## |x-y|^2 <= (8 eps)^2 max(1, |x|^2, |y|^2) for numbers drawn with
  ## precision eps, far below the difference of a site or stream mismatch
  var d, a, b = 0.0
  threads:
    let t = norm2(x - y)
    let u = x.norm2
    let v = y.norm2
    threadMaster:
      d = t
      a = u
      b = v
  d <= (8*eps)*(8*eps)*max(1.0, max(a, b))

proc sameGenerators(rv: RNGFieldV; r: Field): int =
  ## sites of rv whose generator or coordinates differ from those of r
  for i in rv.l.sites:
    let q = r.l.rankIndex(rv.l.coords, i)
    var ok = q.rank == myRank and rv[i] == r[q.index]
    for d in 0..<rv.l.nDim:
      ok = ok and r.l.coords[d][q.index] == rv.l.coords[d][i]
    if not ok: inc result

template checkSites(R: typedesc) =
  var r = lo.newRNGField(R, seed)
  var rv = lo.newRNGFieldV(R, seed)
  var r1 = l1.newRNGField(R, seed)
  var rv1 = lo.newRNGFieldV(R, seed)
  var v, vv = lo.ColorVector()
  var w, ww = l1.ColorVector1()
  var g: R
  let eps = float(epsilon(typeof(gaussian(g))))  # RngMilc6 draws float32
  check sameGenerators(rv, r) == 0
  threads:
    v.gaussian r
    vv.gaussian rv
    w.gaussian r1
    ww.gaussian rv1
  check near(v, vv, eps)
  check near(w, ww, eps)
  check sameGenerators(rv, r) == 0
  check sameGenerators(rv1, r1) == 0

suite "RNGFieldV sites":
  test "RngMilc6": checkSites(RngMilc6)
  test "MRG32k3a": checkSites(MRG32k3a)
  test "Philox4x64": checkSites(Philox4x64)
  test "Threefry4x64": checkSites(Threefry4x64)

qexFinalize()
