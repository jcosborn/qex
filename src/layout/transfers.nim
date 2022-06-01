import comms/comms
import comms/gather
import layout
import field

proc mapcoord(x2: var seq, lo2: Layout, x: seq, lo: Layout2) =
  let n = x2.len
  for i in 0..<n:
    x2[i] = x[i] mod lo2.physGeom[i].int32

# currently only works with non-vectorized layouts
proc replicateFrom*(f: Field, f2: Field2) =
  let lo = f.l
  let lo2 = f2.l
  var rl = newSeq[RecvList](0)
  var x = newSeq[int32](lo.nDim)
  var x2 = newSeq[int32](lo2.nDim)
  for i in 0..<lo.nSites:
    lo.coord(x, i)
    mapcoord(x2, lo2, x, lo)
    let ri = lo2.rankIndex(x2)
    rl.add RecvList(didx: i.int32, srank: ri.rank.int32, sidx: ri.index.int32)

  let c = getDefaultComm()
  let gm = c.makeGatherMap(rl)

  template `&&`(x: Field): untyped = cast[pointer](unsafeAddr(x[0]))
  c.gather(gm, sizeof(f[0]), &&f, &&f2)

proc remapLocalFrom*(f: Field, f2: Field2) =
  let lo = f.l
  let lo2 = f2.l
  var x = newSeq[int32](lo.nDim)
  for i in 0..<lo.nSites:
    lo.coord(x, i)
    let ri = lo2.rankIndex(x)
    assert(ri.rank==lo.myrank)
    let j = ri.index
    when lo.V == 1:
      when lo2.V == 1:
        f[i] := f2[j]
      else:
        f[i] := f2{j}
    else:
      when lo2.V == 1:
        f{i} := f2[j]
      else:
        f{i} := f2{j}

when isMainModule:
  import qex

  qexInit()
  echo "rank ", myRank, "/", nRanks
  var lat = intSeqParam("lat", @[8,8,8,16])
  var t0 = intParam("t0", lat[^1] div 4)
  var t1 = intParam("t1", (3*lat[^1]) div 4)

  var lat2 = lat
  lat2[^1] = t1-t0+1

  var lo = newLayout(lat, 1)
  var lo2 = newLayout(lat2, 1)

  var f = lo.ColorVector1()
  var f2 = lo2.ColorVector1()

  f := 1
  f2 := 2
  let fn1 = f.norm2
  let f2n1 = f2.norm2
  echo "f: ", fn1
  echo "f2: ", f2n1
  f2.replicateFrom(f)
  let fn2 = f.norm2
  let f2n2 = f2.norm2
  echo "f: ", fn2
  echo "f2: ", f2n2
  let f2n = (fn1*lat2[^1].float)/lat[^1].float
  echo "f2: ", f2n

  f2 := 2
  f.replicateFrom(f2)
  let fn3 = f.norm2
  let f2n3 = f2.norm2
  echo "f: ", fn3
  echo "f2: ", f2n3
  let fn = fn1 - f2n + f2n3
  echo "f: ", fn

  qexFinalize()
