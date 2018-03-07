const hdr = currentSourcePath()[0..^11] & "qlayout.h"
{. pragma: ql, header:hdr .}
import layoutTypes
import comms/qmp
import base

type
  ShiftBufQ* = object
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
    activeRecv*: bool
    activeSend*: bool
  ShiftBuf* = ref ShiftBufObj

import qshifts

#[
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
]#

proc freeShiftBuf*(sb:ShiftBuf) =
  if sb.lbuf != nil:
    #let a = unsafeaddr(sb.sq)
    #echo "freeShiftBuf: ", cast[int](a)
    freeShiftBufQ(sb.sq.addr)
    dealloc(sb.lbuf)
    sb.lbuf = nil
proc prepareShiftBuf*(sb:var ShiftBuf, si:ShiftIndices, esize:int) =
  sb.new()
  #sb.new(freeShiftBuf)
  #let a = unsafeaddr(sb.sq)
  #echo "prepareShiftBuf: ", cast[int](a)
  prepareShiftBufQ(sb.sq.addr, si.sq.addr, esize.cint)
  #sb.lbufSize = (si.sq.nSendSites*si.nSitesInner-si.sq.nSendSites1)*esize
  sb.lbufSize = si.sq.nSendSites*si.nSitesInner*esize
  sb.lbuf = cast[type(sb.lbuf)](alloc(sb.lbufSize))
  sb.activeRecv = false
  sb.activeSend = false
template startRecvBuf*(sb: ShiftBuf) =
  if not sb.activeRecv:
    sb.activeRecv = true
    startRecvBufQ(unsafeAddr(sb.sq))
template waitRecvBuf*(sb: ShiftBuf) =
  if sb.activeRecv:
    sb.activeRecv = false
    waitRecvBufQ(unsafeAddr(sb.sq))
proc doneRecvBuf*(sb:ShiftBuf) = doneRecvBufQ(unsafeAddr(sb.sq))
template startSendBuf*(sb: ShiftBuf) =
  if not sb.activeSend:
    sb.activeSend = true
    startSendBufQ(unsafeAddr(sb.sq))
template waitSendBuf*(sb: ShiftBuf) =
  if sb.activeSend:
    sb.activeSend = false
    waitSendBufQ(unsafeAddr(sb.sq))


#proc makeShiftSubQ(si:ptr ShiftIndicesQ; l:ptr LayoutQ; d:ptr cArray[cint];
#                   s:cstring) {.importC:"makeShiftSub",ql.}
#proc makeShiftMultiSubQ(si:ptr ptr ShiftIndicesQ; l:ptr LayoutQ;
#                        d:ptr ptr cint; s:ptr cstring; ndisp:cint)
#  {.importC:"makeShiftMultiSub",ql.}

# x>0 -> 2*x-1; x<= -> -2*x
#proc makeShift(l:Layout; disp:openArray[int]; sub:string="all") =
proc makeShift*(l:var Layout; dir,len:int; sub:string="all") =
  var si = ShiftIndices.new()
  #echo l.nDim
  var disp = newSeq[cint](l.nDim)
  disp[dir] = -len.cint
  makeShiftSubQ(si.sq.addr, l.lq.unsafeAddr, cast[ptr cArray[cint]](disp[0].addr), sub)
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
