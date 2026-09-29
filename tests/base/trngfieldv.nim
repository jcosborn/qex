import testutils
import qex, physics/qcdTypes

# RNGFieldV: site i of its layout holds and draws with the generator that
# newRNGField gives the site of the same coordinates.  Fields of the layout of
# the RNGFieldV and of a one-lane layout draw the numbers of newRNGField.

qexInit()
let lat = latticeFromLocalLattice([4,4,4,8], nRanks)
let lo = lat.newLayout
let l1 = lo.physGeom.newLayout(1, lo.rankGeom)
const seed = 987654321

proc diff(x, y: Field): float =
  var d = 0.0
  threads:
    let t = norm2(x - y)
    threadMaster: d = t
  d

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
  check sameGenerators(rv, r) == 0
  threads:
    v.gaussian r
    vv.gaussian rv
    w.gaussian r1
    ww.gaussian rv1
  check diff(v, vv) == 0.0
  check diff(w, ww) == 0.0
  check sameGenerators(rv, r) == 0
  check sameGenerators(rv1, r1) == 0

suite "RNGFieldV sites":
  test "RngMilc6": checkSites(RngMilc6)
  test "MRG32k3a": checkSites(MRG32k3a)
  test "Philox4x64": checkSites(Philox4x64)
  test "Threefry4x64": checkSites(Threefry4x64)

qexFinalize()
