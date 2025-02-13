import tables
export tables
import base
#import qgather
import layoutTypes
export layoutTypes
import shiftX
export shiftX

import qlayout

var noSplitInnerList = newSeq[int](0)
proc noSplitInnerDim*(d: int) =
  noSplitInnerList.add d

# partition x into n blocks with geometry nx, only works with powers of 2
# dist=0: prefer to split already split direction
# dist=1: split unsplit directions first
proc partitionGeom[T](lx,nx:var openArray[T]; x:openArray[T]; n,dist:int) =
  for i,xi in x:
    nx[i] = 1
    lx[i] = x[i]
  var ww = n
  while ww > 1:
    var k = x.len-1
    while (lx[k] and 1) > 0 or (k in noSplitInnerList):
      k.dec
      if k < 0:
        echo "not enough 2's in partitioned geom:"
        for i in 0..<x.len: echo " ", $x[i]
        echo " /", n
        quit(-1)
    for i in countdown(k-1,0):
      if (lx[i] and 1) == 0 and (k notin noSplitInnerList):
        if dist == 0:
          if lx[i]>lx[k] or (lx[i]==lx[k] and nx[i]>nx[k]): k = i
        else:
          if nx[i]<nx[k] or (nx[i]==nx[k] and lx[i]>lx[k]): k = i
    nx[k] *= 2
    lx[k] = lx[k] div 2
    ww = ww div 2

# partition x into n blocks with geometry nx
# dist=0: prefer to split already split direction
# dist=1: split unsplit directions first
proc partitionGeomF[T](lx,nx:var openArray[T]; x:openArray[T]; n,dist:int) =
  for i,xi in x:
    nx[i] = 1
    lx[i] = xi
  let fs = factor(n)
  for fi in countdown(fs.len-1,0):
    let f = fs[fi]
    var k = lx.len-1
    while (lx[k] mod f) != 0:
      dec k
      if k < 0:
        echo "not enough factors of ", f, " in partitioned geom:"
        for i in 0..<x.len: echo " ", $x[i]
        echo " /", n
        quit(-1)
    for i in countdown(k-1,0):
      if (lx[i] mod f) == 0:
        if dist == 0:
          if lx[i]>lx[k] or (lx[i]==lx[k] and nx[i]>nx[k]): k = i
        else:
          if nx[i]<nx[k] or (nx[i]==nx[k] and lx[i]>lx[k]): k = i
    nx[k] *= f
    lx[k] = lx[k] div f

proc newLayoutX*(comm: Comm; lat: openArray[int]; V: static[int];
                 rg0,ig0: openArray[int]): Layout[V] =
  let nRanks = comm.size
  let myRank = comm.rank
  result.new()
  let nd = lat.len
  var rg,ig,og,lg: seq[int]
  og.newSeq(nd)
  lg.newSeq(nd)

  if rg0.len == 0:
    rg = intSeqParam("rankgeom")
  else:
    echo rg0
    rg.newSeq(nd)
    for i in 0..<nd:
      rg[i] = rg0[i]
  if rg.len == 0:
    rg.setLen(nd)
    partitionGeomF(lg, rg, lat, nRanks, 0)
  else:
    for i in 0..<nd:
      lg[i] = lat[i] div rg[i]

  echo "#physGeom:" & $(@lat)
  echo "#rankGeom:" & $rg
  echo "#localGeom:" & $lg

  if ig0.len == 0:
    ig = intSeqParam("innergeom")
  else:
    ig.newSeq(nd)
    for i in 0..<nd:
      ig[i] = ig0[i]
  if ig.len == 0:
    ig.setLen(nd)
    partitionGeom(og, ig, lg, V, 1)
  else:
    for i in 0..<nd:
      og[i] = lg[i] div ig[i]
  echo "#innerGeom:" & $ig
  echo "#outerGeom:" & $og

  for i in 0..<nd:
    if ig[i]>1 and (og[i] mod 2)==1:
      for j in 1..<nd:
        let k = (i+j) mod nd
        if ig[k]==1 and (og[k] mod 4)==0:
          ig[k] *= 2
          og[k] = og[k] div 2
          ig[i] = ig[i] div 2
          og[i] *= 2
          break
    if (ig[i]>1 and (og[i] mod 2)==1) or (ig[i]>2):
      echo "error: can't lay out inner geom"
      quit -1

  echo "#innerGeom:" & $ig
  echo "#outerGeom:" & $og
  result.lq.myrank = cint(myRank)
  result.lq.nranks = cint(nRanks)
  result.lq.nDim = cint(nd)
  proc ca(a:openArray[int]):ptr cArray[cint] =
    result = cast[ptr cArray[cint]](alloc(a.len*sizeOf(cint)))
    for i,x in a: result[i] = x.cint
  result.lq.physGeom = ca(lat)
  result.lq.rankGeom = ca(rg)
  result.lq.innerGeom = ca(ig)
  #result.lq.outerGeom = ca(og)
  #result.lq.localGeom = ca(lg)
  layoutSetupQ(result.lq)
  result.nDim = nd
  result.physGeom = @lat
  result.localGeom = lg
  result.rankGeom = rg
  result.innerGeom = ig
  result.outerGeom = og
  result.temp = newSeq[int](nd)
  result.physVol = result.lq.physVol.int
  result.nEven = result.lq.nEven.int
  result.nOdd = result.lq.nOdd.int
  result.nSites = result.lq.nSites.int
  result.nEvenOuter = result.lq.nEvenOuter.int
  result.nOddOuter = result.lq.nOddOuter.int
  result.nSitesOuter = result.lq.nSitesOuter.int
  result.nSitesInner = result.lq.nSitesInner.int
  result.comm = comm
  result.nranks = nRanks
  result.myrank = myRank
  result.shifts.init
  result.coords.newSeq(nd)
  for i in 0..<nd: result.coords[i].newSeq(result.nSites)
  var coords = newSeq[cint](nd)
  #let coa = cast[ptr cArray[cint]](addr(coords[0]))
  for i in 0..<result.nSites:
    var li = LayoutIndexQ((rank:myRank.cint,index:i.cint))
    #layoutCoordQ(result.lq.addr, coa, li.addr)
    layoutCoordQ(result.lq.addr, coords, li.addr)
    #echo coords[0]
    for d in 0..<nd: result.coords[d][i] = coords[d].int16
  result.vcoordTemp.newSeq(nd)
template newLayout*(l:openArray[int]; n:static[int], rg,ig: seq[int]):untyped =
  let comm = getDefaultComm()
  newLayoutX(comm, l, n, rg, ig)
template newLayout*(l:openArray[int]; n:static[int], rg: seq[int]):untyped =
  let comm = getDefaultComm()
  newLayoutX(comm, l, n, rg, [])
template newLayout*(l:openArray[int]; n:static[int]):untyped =
  let comm = getDefaultComm()
  newLayoutX(comm, l, n, [], [])
template newLayout*(l:openArray[int]; rg: seq[int]):untyped =
  let comm = getDefaultComm()
  newLayoutX(comm, l, VLEN, rg, [])
template newLayout*(l:openArray[int]):untyped =
  let comm = getDefaultComm()
  newLayoutX(comm, l, VLEN, [], [])

template newLayout*(c:Comm; l:openArray[int]):untyped =
  newLayoutX(c, l, VLEN, [], [])
template newLayout*(c:Comm; l:openArray[int]; n:static[int]):untyped =
  newLayoutX(c, l, n, [], [])

template `[]`*(l: Layout, i: int): untyped = l.physGeom[i]
template `[]`*(l: Layout, i: BackwardsIndex): untyped = l.physGeom[i]

proc localIndex*(l:Layout, coords: var openArray[SomeInteger]): int =
  var li:LayoutIndexQ
  layoutIndexQ(l.lq, li, coords)
  result = li.index.int

proc rankIndex*(l:Layout, coords: ptr cArray[cint]):tuple[rank,index:int] =
  var li:LayoutIndexQ
  layoutIndexQ(l.lq, li, coords)
  result = (rank:li.rank.int,index:li.index.int)
proc rankIndex*(l: Layout, coords: var openArray[cint]): tuple[rank,index:int] =
  var li: LayoutIndexQ
  #layoutIndexQ(l.lq.addr, li.addr,
  #             cast[ptr cArray[cint]](unsafeAddr(coords[0])))
  layoutIndexQ(l.lq, li, coords)
  result = (rank:li.rank.int, index:li.index.int)
proc rankIndex*(l:Layout, coords:openArray[int]):tuple[rank,index:int] =
  when compiles((const n=coords.len;n)):
    const n=coords.len
    var c:array[n,cint]
  else:
    let n=coords.len
    var c = newSeq[cint](n)
  for i in 0..<n: c[i] = coords[i].cint
  result = l.rankIndex(c)
proc rankIndex*(l:Layout, coords: ptr cint):tuple[rank,index:int] =
  #rankIndex(l, cast[ptr cArray[cint]](coords))
  var ca = cast[ptr cArray[cint]](coords)
  rankIndex(l, toOpenArray(ca, 0, l.nDim-1))
#proc coord*(l:Layout, coord:openArray[cint], ri:tuple[rank,index:int]) =
#  var li = LayoutIndexQ(rank:ri.rank.cint, index:ri.index.cint)
#  layoutCoordQ(l.lq.addr, coord, li)
proc coord*(l: Layout, coord: var openArray[cint]; rank,index: int) =
  var li: LayoutIndexQ
  li.rank = rank.cint
  li.index = index.cint
  #layoutCoordQ(l.lq.addr, cast[ptr cArray[cint]](coord[0].addr), li.addr)
  var ca = cast[ptr cArray[cint]](coord[0].addr)
  layoutCoordQ(l.lq.addr, toOpenArray(ca,coord.low,coord.high), li.addr)
proc coord*(l: Layout, crd: var auto, ri: tuple[rank,index:int]) =
  coord(l, crd, ri.rank, ri.index)
template coord*(l: Layout, crd: var typed, index: int) =
  coord(l, crd, l.myrank, index)
proc coord*(l: Layout, coord: ptr cint, ri: tuple[rank,index:cint]) =
  var ca = cast[ptr cArray[cint]](coord)
  if ri.rank == l.myrank:
    for i in 0..<l.nDim:
      ca[i] = l.coords[i][ri.index]
  else:
    var li: LayoutIndexQ
    li.rank = ri.rank.cint
    li.index = ri.index.cint
    #layoutCoordQ(l.lq.addr, cast[ptr cArray[cint]](coord), li.addr)
    layoutCoordQ(l.lq.addr, toOpenArray(ca,0,l.nDim-1), li.addr)

# from physical index to local rank and index
proc rankIndex*(lo: Layout, lex: int): tuple[rank,index:int] =
  let n = lo.nDim
  var k = lex
  for i in 0..<n:
    lo.temp[i] = (k mod lo.physGeom[i]).cint
    k = k div lo.physGeom[i]
  rankIndex(lo, lo.temp)

proc physIndex*(lo: Layout, rank,index: int): int =
  lo.coord(lo.temp, rank, index)
  let n = lo.nDim
  result = 0
  for i in countdown(n-1,0):
    result = (result*lo.physGeom[i]) + lo.temp[i]

proc rankFromRankCoords*(l: Layout, coords: ptr cArray): int =
  for i in 0..<l.nDim:
    l.temp[i] = l.localGeom[i] * coords[i]
  let ri = l.rankIndex(l.temp)
  ri.rank

proc vcoords*[V:static[int]](l:Layout[V]; i:int):seq[array[V,int16]] =
  #for d in 0..<l.nDim:
  #  for j in 0..<V:
  #    l.vcoordTemp[d][j] = l.coords[d][i*V+j]
  #l.vcoordTemp
  result.newSeq(l.nDim)
  for d in 0..<l.nDim:
    for j in 0..<V:
      result[d][j] = l.coords[d][i*V+j]

proc vcoords*[V:static[int]](l:Layout[V]; d,i:int):array[V,int16] =
  for j in 0..<V:
    result[j] = l.coords[d][i*V+j]

proc layoutSubset*(s:var Subset, l:Layout, sub:string) =
  s.low = 0
  s.high = l.nSites
  s.lowOuter = 0
  s.highOuter = l.nSitesOuter
  if sub[0]=='e':
    s.high = l.nEven
    s.highOuter = l.nEvenOuter
  elif sub[0]=='o':
    s.low = l.nEven
    s.lowOuter = l.nEvenOuter
template getSubset*(l:Layout; sub:string):Subset =
  var s{.noInit.}:Subset
  layoutSubset(s, l, sub)
  s
proc paritySubset*(s: var Subset; l: Layout; par: int) =
  if par==0:
    s.layoutSubset(l, "e")
  else:
    s.layoutSubset(l, "o")
proc timesliceSubsetsRange*(l: Layout, i0,i1: int): seq[Subset] =
  # find timeslice subsets within a range of outer sites [i0,i1] inclusive
  doAssert(l.innerGeom[^1] == 1)
  let nt = l.physGeom[^1]
  result.newSeq(nt)
  var x = newSeq[int32](l.nDim)
  for i in 0..<nt:
    result[i].lowOuter = l.nSites
    result[i].highOuter = 0
  for i in i0..i1:
    let k = i * l.V
    let t = l.coords[^1][k]
    result[t].lowOuter = min(result[t].lowOuter, i)
    result[t].highOuter = max(result[t].highOuter, i)
  for i in 0..<nt:
    var a = result[i].lowOuter
    var b = result[i].highOuter
    if a <= b:
      inc b
    else:
      a = 0
      b = 0
    result[i].lowOuter = a
    result[i].highOuter = b
    result[i].low = a * l.V
    result[i].high = b * l.V
proc timesliceSubsets*(l: Layout): array[2,seq[Subset]] =
  result[0] = l.timesliceSubsetsRange(0, l.nEvenOuter-1)
  result[1] = l.timesliceSubsetsRange(l.nEvenOuter, l.nSitesOuter-1)

template `len`*(s:Subset):untyped = s.high-s.low
template `lenOuter`*(s:Subset):untyped = s.highOuter-s.lowOuter

template singleSites*(l: Layout): untyped = 0..<l.nSites
template singleSites*(s: Subset): untyped = s.low ..< s.high

template threadRankSum*(l: Layout, a: varargs[untyped]) =
  l.comm.threadRankSum(a)
