import testutils
import qex, simd, simd/simdArray, comms/halo, base/hyper, sequtils

# Shifts and halo updates (corners and axes) on layouts with 2 sites per lane
# along dimensions split across ranks, outer extent 2, against the global
# coordinates of every lane; on 1, 2 and 4 ranks.  Shifts of a split
# dimension longer than its outer extent are unsupported (shiftX stops).

type SD8 = Simd[SimdArrayObj[8, float]]

proc setf(x: Field, offset: seq[int]) =
  let lo = x.l
  const vl = lo.V
  let lat = lo.physGeom
  let nd = lat.len
  for e in x:
    var t: array[vl,int]
    for i in 0..<nd:
      t = 1 + 10*t + ((offset[i] + lat[i] + lo.vcoords(i,e)) mod lat[i])
    x[e] := t

proc shiftErr(x, y, z: auto; mu, d: int): float =
  let f = newShifter(x, mu, d)
  var offs = newSeq[int](x.l.nDim)
  offs[mu] = d
  z.setf(offs)
  var r = 0.0
  threads:
    y := f ^* x
    let t = norm2(z - y)
    threadMaster: r = t
  r

proc site(lo: auto, e, k: int, o: seq[int32]): int =
  let xs = lo.vcoords(e)
  var x = newSeq[int](lo.nDim)
  for d in 0..<lo.nDim:
    x[d] = (int(xs[d][k]) + int(o[d]) + lo.physGeom[d]) mod lo.physGeom[d]
  lexIndex(x, lo.physGeom)

proc nbr(hl: auto, i: int, o: seq[int32]): int =
  result = i
  for d in 0..<o.len:
    for s in 1..abs(o[d]):
      if result < 0: return
      result = if o[d] > 0: hl.neighborFwd[d][result] else: hl.neighborBck[d][result]

proc haloErr(g: auto, hl: auto, hm: auto, offsets: seq[seq[int32]]): int =
  let lo = g.l
  const V = lo.V
  let z = newSeq[int32](lo.nDim)
  for e in 0..<lo.nSitesOuter:
    var t: array[V,float]
    for k in 0..<V: t[k] = float lo.site(e, k, z)
    g[e] := t
  let h = makeHalo(hl, g)
  h.update hm, getDefaultComm()
  for i in 0..<hl.nOut:
    for o in offsets:
      let j = hl.nbr(i, o)
      if j < 0:
        inc result
        continue
      for k in 0..<V:
        let want = float lo.site(i, k, o)
        let got = h[j][k]
        if got != want:
          if result < 3: echo "  halo site ", i, " offset ", o, " lane ", k, ": ", got, " != ", want
          inc result

qexInit()

proc run(lat, rg, ig: seq[int]) =
  let lo = newLayout(lat, 8, rg, ig)
  echo "lat ", lat, " rankGeom ", rg, " innerGeom ", ig, " outerGeom ", lo.outerGeom
  var x, y, z: Field[8, SD8]
  x.new(lo)
  y.new(lo)
  z.new(lo)
  x.setf(newSeq[int](lat.len))
  var bad = 0
  for mu in 0..<lat.len:
    let dmax = if rg[mu] > 1: min(3, lo.outerGeom[mu]) else: 3
    for d in -dmax..dmax:
      if d == 0: continue
      let r = shiftErr(x, y, z, mu, d)
      if r != 0.0:
        echo "  shift mu ", mu, " d ", d, ": ", r
        inc bad
  let nd = lat.len
  var corners, axes: seq[seq[int32]]
  for c in 0..<(1 shl nd):
    var t = newSeq[int32](nd)
    for d in 0..<nd: t[d] = if ((c shr d) and 1) == 1: -1 else: 1
    corners.add t
  for d in 0..<nd:
    for s in [-1'i32, 1'i32]:
      var t = newSeq[int32](nd)
      t[d] = s
      axes.add t
  let hl = haloLayout(lo, newSeqWith(nd, 1), newSeqWith(nd, 1))
  let hc = haloErr(x, hl, haloMap(hl, getDefaultComm(), corners), corners)
  let ha = haloErr(x, hl, haloMap(hl, getDefaultComm(), axes), axes)
  echo "  shift errors ", bad, ", halo corner errors ", hc, ", axis errors ", ha
  check(bad == 0 and hc == 0 and ha == 0)

suite "split dimensions with 2 sites per lane":
  case nRanks
  of 1:
    test "t outer 2, not split": run(@[8,8,8,4], @[1,1,1,1], @[2,2,1,2])
  of 2:
    test "t outer 4, split": run(@[8,8,8,16], @[1,1,1,2], @[2,2,1,2])
    test "t outer 2, split": run(@[8,8,8,8], @[1,1,1,2], @[2,2,1,2])
  of 4:
    test "z, t outer 4, split": run(@[8,8,8,16], @[1,1,2,2], @[2,2,1,2])
    test "t outer 2, split": run(@[8,8,8,8], @[1,1,2,2], @[2,2,1,2])
    test "z outer 2, split": run(@[8,8,8,8], @[1,1,2,2], @[2,2,2,1])
    test "t outer 2, 4 ranks along t": run(@[8,8,8,16], @[1,1,1,4], @[2,2,1,2])
  else: discard

qexFinalize()
