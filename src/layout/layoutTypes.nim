import base
import tables
export tables
import qgather

type
  llist* = object
    value*: pointer
    next*: ptr llist

# k = sum_i sum_j ((x[i]/d[i][j])%m[i][j])*f[i][j]
# x[i] = sum_j ((k/f[i][j])%m[i][j])*d[i][j]
# parity?
type
  LayoutQ* = object
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
  SubsetQ* = object
    begin*: cint
    `end`*: cint
    beginOuter*: cint
    endOuter*: cint

type ShiftIndicesQ* = object
  gi*: ptr GatherIndices
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
  Layout2*[V:static[int]] = Layout[V]
  Layout3*[V:static[int]] = Layout[V]
  Subset* = object
    low*: int
    high*: int
    lowOuter*: int
    highOuter*: int

proc makeShiftKey*(d,l:int;s:string):ShiftKey =
  ShiftKey((dir:d,len:l,sub:s))

#  result.shifts = initTable[ShiftKey,ShiftIndices]()
proc init*(x: var Table) =
  x = initTable[x.A,x.B]()
