const hdr = currentSourcePath()[0..^12] & "qlayout.h"
{. pragma: ql, header:hdr .}
{. passC:"-I." .}
{. compile:"qlayout.c" .}
{. compile:"qshifts.c" .}
{. compile:"qgather.c" .}

import tables
export tables
import base
import comms/qmp

type
  llist* {.importc,ql.} = object
    value*: pointer
    next*: ptr llist

# k = sum_i sum_j ((x[i]/d[i][j])%m[i][j])*f[i][j]
# x[i] = sum_j ((k/f[i][j])%m[i][j])*d[i][j]
# parity?
type
  LayoutQ* {.importC:"Layout",ql.} = object
    nDim*: cint
    physGeom*: ptr cArray[cint]
    rankGeom*: ptr cArray[cint]
    innerGeom*: ptr cArray[cint]      #wrap
    outerGeom*: ptr cArray[cint]      #wls
    localGeom*: ptr cArray[cint]
    physVol*: cint
    nEven*: cint
    nOdd*: cint
    nSites*: cint
    nEvenOuter*: cint
    nOddOuter*: cint
    nSitesOuter*: cint
    nSitesInner*: cint
    innerCb*: cint
    innerCbDir*: cint
    shifts*: ptr llist
    nranks*: cint
    myrank*: cint
  LayoutIndexQ* = tuple[rank,index:cint]
  SubsetQ* {.importc:"Subset",ql.} = object
    begin*: cint
    `end`*: cint
    beginOuter*: cint
    endOuter*: cint

proc layoutSetupQ*(l:ptr LayoutQ) {.importC:"layoutSetup",ql.}
proc layoutIndexQ*(l:ptr LayoutQ; li:ptr LayoutIndexQ, coords:ptr cint)
  {.importC:"layoutIndex",ql.}
proc layoutCoordQ*(l:ptr LayoutQ; coords:ptr cint, li:ptr LayoutIndexQ)
  {.importC:"layoutCoord",ql.}

type GatherIndicesQ* {.importc:"GatherIndices",ql.} = object

type ShiftIndicesQ* {.importc:"ShiftIndices",ql.} = object
  gi*: GatherIndicesQ
  disp*: ptr cArray[cint]
  sidx*: ptr cArray[cint]
  pidx*: ptr cArray[cint]
  nRecvRanks*: cint
  recvRanks*: ptr cArray[cint]
  recvRankSizes*: ptr cArray[cint]
  recvRankSizes1*: ptr cArray[cint]
  recvRankOffsets*: ptr cArray[cint]
  recvRankOffsets1*: ptr cArray[cint]
  nRecvSites*: cint
  nRecvSites1*: cint
  nRecvDests*: cint
  recvDests*: ptr cArray[cint]
  recvLocalSrcs*: ptr cArray[cint]
  recvRemoteSrcs*: ptr cArray[cint]
  nSendRanks*: cint
  sendRanks*: ptr cArray[cint]
  sendRankSizes*: ptr cArray[cint]
  sendRankSizes1*: ptr cArray[cint]
  sendRankOffsets*: ptr cArray[cint]
  sendRankOffsets1*: ptr cArray[cint]
  nSendSites*: cint
  nSendSites1*: cint
  sendSites*: ptr cArray[cint]
  vv*: cint
  #offr*: cint
  #lenr*: cint
  #nthreads*: cint
  perm*: cint
  pack*: cint
  blend*: cint
  #sqmpmem*: QMP_msgmem_t
  #smsg*: QMP_msghandle_t
  #rqmpmem*: QMP_msgmem_t
  #rmsg*: QMP_msghandle_t
  #pairmsg*: QMP_msghandle_t

type ShiftIndices* = ref object
  sq*: ShiftIndicesQ
  nRecvRanks*: int
  nRecvDests*: int
  nSendRanks*: int
  nSendSites*: int
  sendSites*: seq[int32]
  perm*: int
  pack*: int
  blend*: int
  nSitesInner*: int

type
  ShiftKey = tuple[dir,len:int;sub:string]
  #Layout[D:static[int]]* = ref object
  Layout*[V:static[int]] = ref object
    nDim*: int
    physGeom*:seq[int]
    rankGeom*:seq[int]
    innerGeom*:seq[int]
    outerGeom*:seq[int]
    localGeom*:seq[int]
    physVol*: int
    nEven*: int
    nOdd*: int
    nSites*: int
    nEvenOuter*: int
    nOddOuter*: int
    nSitesOuter*: int
    nSitesInner*: int
    nranks*: int
    myrank*: int
    lq*: LayoutQ
    shifts*: Table[ShiftKey,ShiftIndices]
    coords*: seq[seq[int16]]
    vcoordTemp*: seq[array[V,int16]]
  Subset* = object
    low*: int
    high*: int
    lowOuter*: int
    highOuter*: int

proc makeShiftKey*(d,l:int;s:string):ShiftKey = ShiftKey((dir:d,len:l,sub:s))
#proc layoutCoord*(l:Layout; coords:var openArray[int]; rank,index:int) =
#  var c = cast[ptr array[0,cint]](coords[0].addr)
#  var liq = LayoutIndexQ((rank.cint,index.cint))
#  layoutCoordQ(l.lq.addr, cast[ptr cint](c), liq.addr)
#proc layoutIndexQ*(l:ptr LayoutQ; li:ptr LayoutIndexQ, coords:ptr cint)
#  {.importC:"layoutIndex",ql.}

# partition x into n blocks with geometry nx
# dist=0: prefer to split already split direction
# dist=1: split unsplit directions first
proc partitionGeom(lx,nx:var openArray[int]; x:openArray[int]; n,dist:int) =
  for i,xi in x:
    nx[i] = 1
    lx[i] = x[i]
  var ww = n
  while ww > 1:
    var k = x.len-1
    while (lx[k] and 1) > 0:
      k.dec
      if k < 0:
        echo "not enough 2's in partitioned geom:"
        for i in 0..<x.len: echo " ", x[i]
        echo " /", n
        quit(-1)
    for i in countdown(k-1,0):
      if (lx[i] and 1) == 0:
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
proc partitionGeomF(lx,nx:var openArray[int]; x:openArray[int]; n,dist:int) =
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
        for i in 0..<x.len: echo " ", x[i]
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

proc newLayoutX*(lat:openArray[int],V:static[int]):Layout[V] =
  result.new()
  let nd = lat.len
  var rg,ig,og,lg:seq[int]
  #rg.newSeq(nd)
  #ig.newSeq(nd)
  og.newSeq(nd)
  lg.newSeq(nd)

  rg = intSeqParam("rankgeom")
  if rg.len == 0:
    rg.setLen(nd)
    partitionGeomF(lg, rg, lat, nRanks, 1)
  else:
    for i in 0..<nd:
      lg[i] = lat[i] div rg[i]

  echo "#physGeom:" & $(@lat)
  echo "#rankGeom:" & $rg
  echo "#localGeom:" & $lg

  ig = intSeqParam("innergeom")
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
    if ig[i]>1 and (og[i] mod 2)==1:
      echo "error: can't layout inner geom"
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
  layoutSetupQ(result.lq.addr)
  result.nDim = nd
  result.physGeom = @lat
  result.localGeom = lg
  result.physVol = result.lq.physVol.int
  result.nEven = result.lq.nEven.int
  result.nOdd = result.lq.nOdd.int
  result.nSites = result.lq.nSites.int
  result.nEvenOuter = result.lq.nEvenOuter.int
  result.nOddOuter = result.lq.nOddOuter.int
  result.nSitesOuter = result.lq.nSitesOuter.int
  result.nSitesInner = result.lq.nSitesInner.int
  result.nranks = nRanks
  result.myrank = myRank
  result.shifts = initTable[ShiftKey,ShiftIndices]()
  result.coords.newSeq(nd)
  for i in 0..<nd: result.coords[i].newSeq(result.nSites)
  var coords = newSeq[cint](nd)
  for i in 0..<result.nSites:
    var li = LayoutIndexQ((rank:myRank.cint,index:i.cint))
    layoutCoordQ(result.lq.addr, coords[0].addr, li.addr)
    #echo coords[0]
    for d in 0..<nd: result.coords[d][i] = coords[d].int16
  result.vcoordTemp.newSeq(nd)
template newLayout*(l:openArray[int]; n:static[int]):untyped = newLayoutX(l,n)
template newLayout*(l:openArray[int]):untyped = newLayoutX(l,VLEN)

proc rankIndex*(l:Layout, coords:openArray[cint]):tuple[rank,index:int] =
  var li:LayoutIndexQ
  layoutIndexQ(l.lq.addr, li.addr, unsafeAddr(coords[0]))
  result = (rank:li.rank.int,index:li.index.int)
proc rankIndex*(l:Layout, coords:openArray[int]):tuple[rank,index:int] =
  when compiles((const n=coords.len;n)):
    const n=coords.len
    var c:array[n,cint]
  else:
    let n=coords.len
    var c = newSeq[cint](n)
  for i in 0..<n: c[i] = coords[i].cint
  result = l.rankIndex(c)
#proc coord*(l:Layout, coord:openArray[cint], ri:tuple[rank,index:int]) =
#  var li = LayoutIndexQ(rank:ri.rank.cint, index:ri.index.cint)
#  layoutCoordQ(l.lq.addr, coord, li)
proc coord*(l:Layout, coord:var openArray[cint], ri:tuple[rank,index:int]) =
  var li: LayoutIndexQ
  li.rank = ri.rank.cint
  li.index = ri.index.cint
  layoutCoordQ(l.lq.addr, coord[0].addr, li.addr)

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
template `len`*(s:Subset):untyped = s.high-s.low
template `lenOuter`*(s:Subset):untyped = s.highOuter-s.lowOuter

proc makeShiftSubQ(si:ptr ShiftIndicesQ; l:ptr LayoutQ; d:ptr cint; s:cstring)
  {.importC:"makeShiftSub",ql.}
proc makeShiftMultiSubQ(si:ptr ptr ShiftIndicesQ; l:ptr LayoutQ;
                        d:ptr ptr cint; s:ptr cstring; ndisp:cint)
  {.importC:"makeShiftMultiSub",ql.}

# x>0 -> 2*x-1; x<= -> -2*x
#proc makeShift(l:Layout; disp:openArray[int]; sub:string="all") =
proc makeShift*(l:var Layout; dir,len:int; sub:string="all") =
  var si = ShiftIndices.new()
  #echo l.nDim
  var disp = newSeq[cint](l.nDim)
  disp[dir] = -len.cint
  makeShiftSubQ(si.sq.addr, l.lq.unsafeAddr, disp[0].addr, sub)
  let key = makeShiftKey(dir, len, sub)
  l.shifts[key] = si
  si.nRecvRanks = si.sq.nRecvRanks
  si.nRecvDests = si.sq.nRecvDests
  si.nSendRanks = si.sq.nSendRanks
  si.nSendSites = si.sq.nSendSites
  si.sendSites.newSeq(si.nSendSites)
  for i in 0..<si.nSendSites: si.sendSites[i] = si.sq.sendSites[i]
  si.perm = si.sq.perm
  si.pack = si.sq.pack
  si.blend = si.sq.blend
  si.nSitesInner = l.nSitesInner
proc getShift*(l:var Layout; dir,len:int; sub:string="all"):ShiftIndices =
  let key = makeShiftKey(dir, len, sub)
  if not hasKey(l.shifts, key):
    makeShift(l, dir, len, sub)
  result = l.shifts[key]

type
  ShiftBufQ* {.importc:"ShiftBuf",ql.} = object
    sqmpmem*: QMP_msgmem_t
    smsg*: QMP_msghandle_t
    rqmpmem*: QMP_msgmem_t
    rmsg*: QMP_msghandle_t
    pairmsg*: QMP_msghandle_t
    sbuf*: ptr cArray[char]
    rbuf*: ptr cArray[char]
    sbufSize*: cint
    rbufSize*: cint
    first*: cint
    offr*: ptr cArray[cint]
    lenr*: ptr cArray[cint]
    nthreads*: ptr cArray[cint]
  ShiftBufObj* = object
    sq*:ShiftBufQ
    lbufSize*: int
    lbuf*: ptr cArray[char]
  ShiftBuf* = ref ShiftBufObj

proc prepareShiftBufQ*(sb:ptr ShiftBufQ, si:ptr ShiftIndicesQ, esize:cint)
  {.importc:"prepareShiftBuf",ql.}
proc freeShiftBufQ*(sb:ptr ShiftBufQ)
  {.importc:"freeShiftBuf", ql.}
proc startRecvBufQ*(sb:ptr ShiftBufQ)
  {.importc:"startRecvBuf", ql.}
proc waitRecvBufQ*(sb:ptr ShiftBufQ)
  {.importc:"waitRecvBuf", ql.}
proc doneRecvBufQ*(sb:ptr ShiftBufQ)
  {.importc:"doneRecvBuf", ql.}
proc startSendBufQ*(sb:ptr ShiftBufQ)
  {.importc:"startSendBuf", ql.}
proc waitSendBufQ*(sb:ptr ShiftBufQ)
  {.importc:"waitSendBuf", ql.}

proc freeShiftBuf*(sb:ShiftBuf) =
  if sb.lbuf != nil:
    let a = unsafeaddr(sb.sq)
    #echo "freeShiftBuf: ", cast[int](a)
    freeShiftBufQ(sb.sq.addr)
    dealloc(sb.lbuf)
    sb.lbuf = nil
proc prepareShiftBuf*(sb:var ShiftBuf, si:ShiftIndices, esize:int) =
  #sb.new()
  sb.new(freeShiftBuf)
  let a = unsafeaddr(sb.sq)
  #echo "prepareShiftBuf: ", cast[int](a)
  prepareShiftBufQ(sb.sq.addr, si.sq.addr, esize.cint)
  #sb.lbufSize = (si.sq.nSendSites*si.nSitesInner-si.sq.nSendSites1)*esize
  sb.lbufSize = si.sq.nSendSites*si.nSitesInner*esize
  sb.lbuf = cast[type(sb.lbuf)](alloc(sb.lbufSize))
proc startRecvBuf*(sb:ShiftBuf) = startRecvBufQ(unsafeAddr(sb.sq))
proc waitRecvBuf*(sb:ShiftBuf) = waitRecvBufQ(unsafeAddr(sb.sq))
proc doneRecvBuf*(sb:ShiftBuf) = doneRecvBufQ(unsafeAddr(sb.sq))
proc startSendBuf*(sb:ShiftBuf) = startSendBufQ(unsafeAddr(sb.sq))
proc waitSendBuf*(sb:ShiftBuf) = waitSendBufQ(unsafeAddr(sb.sq))

when isMainModule:
  import qex
  qexInit()
  echo "rank ", myRank, "/", nRanks
  var lat = [4,4,4,4]
  var lo = newLayout(lat)
  lo.makeShift(0,1)
  lo.makeShift(3,-2,"even")
  for i in 0..<lo.nSites:
    let x = lo.vcoords(i)
  qexFinalize()
