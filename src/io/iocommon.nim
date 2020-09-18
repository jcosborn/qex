import base, layout

proc listNumIoRanks*[T](rankGeom: seq[T]): seq[int32] =
  result.newSeq(1)
  result[0] = 1
  let nd = rankGeom.len
  var r = 1'i32
  for i in countdown(nd-1,0):
    let g = int32 rankGeom[i]
    let s = getDivisors(g)
    for j in 1..<s.len:
      result.add r*s[j]
    r *= g

proc getClosestNumRanks*(rl: seq[int32], nr: int32): int32 =
  for n in rl:
    result = n
    if n >= nr:
      break

proc getIoGeom*[T](rankGeom: seq[T], nr: int32): seq[int32] =
  let nd = rankGeom.len
  result.newSeq(nd)
  var r = nr
  for i in countdown(nd-1,0):
    let n = min(int32 rankGeom[i], r)
    result[i] = n
    r = r div n

proc getIoRanks*(l: Layout, iogeom: seq[int32]): seq[int32] =
  let nr = l.nRanks
  let nd = l.nDim
  result.newSeq(nr)
  var coords = newSeq[int32](nd)
  for r in 0..<nr:
    l.coord(coords, r, 0)
    for i in 0..<nd:
      let t = (coords[i]*iogeom[i]) div l.physGeom[i]
      coords[i] = (t*l.physGeom[i]).int32 div iogeom[i]
    let ri = l.rankindex(coords)
    result[r] = ri.rank.int32

when isMainModule:
  #import qex
  #qexInit()
  let g = @[4,5,6,7]
  echo g
  let lnr = listNumIoRanks(g)
  echo lnr
  let nr = getClosestNumRanks(lnr, 20)
  echo nr
  let iog = getIoGeom(g, nr)
  echo iog
  #let l = newLayout(g,1)
  #let ior = getIoRanks(l, iog)
  #echo ior
  #qexFinalize()
