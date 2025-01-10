import qgather, layoutTypes, shiftX, base
import comms/qmp
import strutils, strformat
import bitops

const
  PAIR = true
  MAXTHREADS = 512

proc aalloc(n: SomeInteger): pointer =
  let a = 64'u
  let b = a + sizeof(uint).uint
  let x = cast[uint](allocShared(n.uint+b))
  let y = a * ((x+b) div a)
  let z = cast[ptr uint](y - sizeof(uint).uint)
  z[] = x
  cast[pointer](y)

proc afree(y: pointer) =
  let z = cast[ptr pointer](cast[uint](y) - sizeof(uint).uint)
  deallocShared(z[])

proc prepareShiftBufsQ*(sb: openArray[ptr ShiftBufQ];
                        si: openArray[ptr ShiftIndicesQ];
                        n: cint; esize: cint) =
  var
    sbs = 0
    rbs = 0
  for i in 0..<n:
    sbs += esize * si[i].nSendSites1
    rbs += esize * si[i].nRecvSites1
  for i in 0..<n:
    sb[i].sbufSize = int32 sbs
    sb[i].rbufSize = int32 rbs
    #sb[i].smsg.clear
    #sb[i].rmsg.clear
    #echo &"sb[{i}].rmsg.isEmpty: {sb[i].rmsg.isEmpty}"
    sb[i].first = 0
    sb[i].offr = cast[type(sb[i].offr)](allocShared(MAXTHREADS * sizeof(cint)))
    sb[i].lenr = cast[type(sb[i].lenr)](allocShared(MAXTHREADS * sizeof(cint)))
    sb[i].nthreads = cast[type(sb[i].nthreads)](allocShared(MAXTHREADS*sizeof(cint)))
    for j in 0..<MAXTHREADS:
      sb[i].nthreads[j] = 0
  sb[0].first = 1
  if sbs > 0:
    var sbuf = aalloc(sbs)
    #printf("sbuf: %p\n", sbuf);
    for i in 0..<n:
      sb[i].sbuf = cast[type(sb[i].sbuf)](sbuf)
      sb[i].nsend = si[i].nSendRanks
      if si[i].nSendRanks > 0:
        sb[i].sqmpmem = cast[ptr cArray[QMP_msgmem_t]](allocShared(si[i].nSendRanks*sizeof(QMP_msgmem_t)))
        sb[i].smsg = cast[ptr cArray[QMP_msghandle_t]](allocShared(si[i].nSendRanks*sizeof(QMP_msghandle_t)))
        for j in 0..<si[i].nSendRanks:
          let sbo = cast[pointer](cast[int](sbuf)+esize*si[i].sendRankOffsets1[j])
          sb[i].sqmpmem[][j] = QMP_declare_msgmem(sbo, csize_t esize*si[i].sendRankSizes1[j])
          sb[i].smsg[][j] = QMP_declare_send_to(sb[i].sqmpmem[][j], si[i].sendRanks[j], 0)
          #echo &"->{si[i].sendRanks[0]}: {sbs}"
  if rbs > 0:
    var rbuf = aalloc(rbs)
    for i in 0..<n:
      sb[i].rbuf = cast[type(sb[i].rbuf)](rbuf)
      sb[i].nrecv = si[i].nRecvRanks
      if si[i].nRecvRanks > 0:
        sb[i].rqmpmem = cast[ptr cArray[QMP_msgmem_t]](allocShared(si[i].nRecvRanks*sizeof(QMP_msgmem_t)))
        sb[i].rmsg = cast[ptr cArray[QMP_msghandle_t]](allocShared(si[i].nRecvRanks*sizeof(QMP_msghandle_t)))
        for j in 0..<si[i].nRecvRanks:
          let rbo = cast[pointer](cast[int](rbuf)+esize*si[i].recvRankOffsets1[j])
          sb[i].rqmpmem[][j] = QMP_declare_msgmem(rbo, csize_t esize*si[i].recvRankSizes1[j])
          sb[i].rmsg[][j] = QMP_declare_receive_from(sb[i].rqmpmem[][j], si[i].recvRanks[j], 0)
        #echo &"<-{si[i].recvRanks[0]}: {rbs}"
        #echo &"sb[{i}].rmsg.isEmpty: {sb[i].rmsg.isEmpty}"
  when PAIR:
    var nmsg = 0
    for i in 0..<n:
      nmsg += si[i].nSendRanks
      nmsg += si[i].nRecvRanks
    var p = newSeq[QMP_msghandle_t](nmsg)
    var nn = cint 0
    for i in 0..<n:
      for j in 0..<si[i].nRecvRanks:
        p[nn] = sb[i].rmsg[][j]
        inc nn
      for j in 0..<si[i].nSendRanks:
        p[nn] = sb[i].smsg[][j]
        inc nn
    var pairmsg: QMP_msghandle_t
    if nn > 0:
      pairmsg = QMP_declare_send_recv_pairs(p[0].addr, nn)
    for i in 0..<n:
      sb[i].pairmsg = pairmsg
      #printf("pair[%i]: %p\t%p\t%p\n",i,sb[i]->rmsg,sb[i]->smsg,sb[i]->pairmsg);
    #fflush(stdout);

proc prepareShiftBufQ*(sb: ptr ShiftBufQ; si: ptr ShiftIndicesQ; esize: cint) =
  prepareShiftBufsQ([sb], [si], 1, esize)

proc startSendBufQ*(sb: ptr ShiftBufQ) =
  #printf("send: %g\n",*(float *)(sb->sbuf));
  #echo "startSendBufQ"
  when PAIR:
    if not isEmpty sb.pairmsg:
      #echo "QMP_start"
      discard QMP_start(sb.pairmsg)
  else:
    for i in 0..<sb.nsend:
      discard QMP_start(sb.smsg[][i])

proc startRecvBufQ*(sb: ptr ShiftBufQ) =
  when not PAIR:
    for i in 0..<sb.nrecv:
      discard QMP_start(sb.rmsg[][i])

proc waitSendBufQ*(sb: ptr ShiftBufQ) =
  when not PAIR:
    for i in 0..<sb.nsend:
      discard QMP_wait(sb.smsg[][i])

proc waitRecvBufQ*(sb: ptr ShiftBufQ) =
  when PAIR:
    if not isEmpty sb.pairmsg:
      discard QMP_wait(sb.pairmsg)
  else:
    for i in 0..<sb.nrecv:
      QMP_wait(sb.rmsg[][i])
  #printf("recv: %g\n",*(float *)(sb->rbuf));

proc doneRecvBufQ*(sb: ptr ShiftBufQ) =
  when PAIR:
    if not isEmpty sb.pairmsg:
      discard QMP_clear_to_send(sb.pairmsg, QMP_CTS_READY)

proc freeShiftBufsQ*(sb: openArray[ptr ShiftBufQ]) =
  let n = sb.len
  for i in 0..<n:
    deallocShared(sb[i].offr)
    deallocShared(sb[i].lenr)
    deallocShared(sb[i].nthreads)
  #FIXME: free all messages
  when PAIR:
    for i in 0..<n:
      if sb[i].first!=0 and not isEmpty sb[i].pairmsg:
        QMP_free_msghandle(sb[i].pairmsg)
      sb[i].pairmsg.clear
      #if not isEmpty sb[i].smsg:
      #  sb[i].smsg.clear
      #  QMP_free_msgmem(sb[i].sqmpmem)
      #  sb[i].sqmpmem = nil
      #if not isEmpty sb[i].rmsg:
      #  sb[i].rmsg.clear
      #  QMP_free_msgmem(sb[i].rqmpmem)
      #  sb[i].rqmpmem = nil
  else:
    for i in 0..<n:
      discard
      #if sb[i].smsg:
      #  QMP_free_msghandle(sb[i].smsg)
      #  sb[i].smsg = nil
      #  QMP_free_msgmem(sb[i].sqmpmem)
      #  sb[i].sqmpmem = nil
      #if sb[i].rmsg:
      #  QMP_free_msghandle(sb[i].rmsg)
      #  sb[i].rmsg = nil
      #  QMP_free_msgmem(sb[i].rqmpmem)
      #  sb[i].rqmpmem = nil
  for i in 0..<n:
    if sb[i].first!=0 and sb[i].sbufSize > 0: afree(sb[i].sbuf)
    if sb[i].first!=0 and sb[i].rbufSize > 0: afree(sb[i].rbuf)

proc freeShiftBufQ*(sb: ptr ShiftBufQ) =
  freeShiftBufsQ([sb])

import qlayout, base

type
  mapargs* = object
    l*: ptr LayoutQ
    disp*: ptr cArray[cint]
    parity*: cint

proc map*(sr: var cint; si: var cint; dr: cint; di: var cint; args: pointer) =
  var ma: ptr mapargs = cast[ptr mapargs](args)
  var l: ptr LayoutQ = ma.l
  var nd: cint = l.nDim
  if di >= 0:  #  (dr,di) -> (sr,si)
    var x = newSeq[cint](nd)
    var
      dli: LayoutIndexQ
      sli: LayoutIndexQ
    dli.rank = dr
    dli.index = di
    layoutCoordQ(l, x, addr(dli))
    var y = newSeq[cint](nd)
    var p: cint = 0
    var k: cint = 0
    while k < nd:
      inc(p, x[k])
      y[k] = (x[k] - ma.disp[k] + l.physGeom[k]) mod l.physGeom[k]
      inc(k)
    if ma.parity >= 0 and (p and 1) != ma.parity:
      sr = -1
      si = -1
    else:
      layoutIndexQ(l[], sli, y)
      sr = sli.rank
      si = sli.index
  else:  #  (dr,-di-1) -> (sr,si) increasing di0 until get sr
    # search for site after or including di0 from rank sr to dr
    var di0: cint = - (di + 1)
    while di0 < l.nSites:
      var sr0: cint
      map(sr0, si, dr, di0, args)
      if sr0 == sr:
        di = di0
        return
      inc(di0)
    si = -1

# nRecvRanks (remote ranks)
# start recvs
# nSendRanks
# start sends
# local + perm
# recv buf
# nSendRanks
# - sendRanks
# - nSendPacks
# - - sendPacks
# - - nSendSites
# - - - sendSites

template SUB2PAR*(s: typed): untyped =
  (if (s)[0] == 'e': 0 else: (if (s)[0] == 'o': 1 else: -1))

proc makeGDFromShiftSubs*(gd: ptr GatherDescription; l: ptr LayoutQ;
                          disps: openArray[ptr cArray[cint]];
                          subs: openArray[cstring]; ndisps: cint) =
  var myRank = l.myrank
  var myndi = l.nSites
  var nndi = ndisps * myndi
  var args: mapargs
  args.l = l
  var sidx: ptr cArray[cint] = cast[ptr cArray[cint]](allocShared(nndi * sizeof(cint)))
  var srank: ptr cArray[cint] = cast[ptr cArray[cint]](allocShared(nndi * sizeof(cint)))
  # find shift sources
  var nRecvDests: cint = 0
  var n: cint = 0
  while n < ndisps:
    var n0 = n * myndi
    args.disp = disps[n]
    args.parity = SUB2PAR(subs[n])
    ##pragma omp parallel for reduction(+:nRecvDests)
    var di = 0
    while di < myndi:
      var
        sr: cint
        si: cint
        di0 = int32 di
      map(sr, si, myRank, di0, addr args)
      srank[n0 + di] = sr
      sidx[n0 + di] = si
      if sr != myRank: inc(nRecvDests)
      inc(di)
    inc(n)
  gd.myRank = myRank
  gd.nIndices = int32 nndi
  gd.srcRanks = srank
  gd.srcIndices = sidx
  gd.nRecvDests = nRecvDests
  # use inverse map
  var nd: cint = l.nDim
  var dispi = newSeq[cint](nd)
  args.disp = cast[ptr cArray[cint]](addr dispi[0])
  var sendSrcIndices = newSeq[cint]()
  var sendDestRanks = newSeq[cint]()
  var sendDestIndices = newSeq[cint]()
  var tlen: array[MAXTHREADS, cint]
  # find who to send to
  n=0
  while n < ndisps:
    var n0 = n * myndi
    var sp = 0
    for i in 0..<nd:
      sp += abs(disps[n][i])
      dispi[i] = - disps[n][i]
    args.parity = SUB2PAR(subs[n])
    if (sp and 1) == 1 and args.parity >= 0:
      args.parity = 1 - args.parity
    #int tid = THREADNUM;
    #int nid = NUMTHREADS;
    var tid = 0
    var nid = 1
    var sendSrcIndicesT = newSeq[cint]()
    var sendDestRanksT = newSeq[cint]()
    var sendDestIndicesT = newSeq[cint]()
    ##pragma omp for
    var di = 0'i32
    while di < myndi:
      var
        dr: cint = myRank
        sr: cint
        si: cint
      map(sr, si, dr, di, addr args)
      if sr >= 0 and si >= 0 and sr != myRank:
        if tid == 0:
          sendSrcIndices.add di
          sendDestRanks.add sr
          sendDestIndices.add int32(n0+si)
        else:
          sendSrcIndicesT.add di
          sendDestRanksT.add sr
          sendDestIndicesT.add int32(n0+si)
      inc(di)
    tlen[tid] = sendSrcIndicesT.len.cint
    #TBARRIER;
    var i0: cint = cint sendSrcIndices.len
    for i in 0..<tid:
      inc(i0, tlen[i])
    if tid == nid-1:
      var ln = i0 + sendSrcIndicesT.len
      sendSrcIndices.setLen ln
      sendDestRanks.setLen ln
      sendDestIndices.setLen ln
    #TBARRIER;
    for i in 0..<sendSrcIndicesT.len:
      sendSrcIndices[i0 + i] = sendSrcIndicesT[i]
      sendDestRanks[i0 + i] = sendDestRanksT[i]
      sendDestIndices[i0 + i] = sendDestIndicesT[i]
    # end parallel
    inc(n)
  gd.nSendIndices = cint sendSrcIndices.len
  template ARRAY_CLONE(x,y: typed) =
    x = cast[type(x)](allocShared(y.len*sizeof(type(x[0]))))
    for i in 0..<y.len: x[i] = y[i]
  ARRAY_CLONE(gd.sendSrcIndices, sendSrcIndices)
  ARRAY_CLONE(gd.sendDestRanks, sendDestRanks)
  ARRAY_CLONE(gd.sendDestIndices, sendDestIndices)

proc makeGDFromShifts*(gd: ptr GatherDescription; l: ptr LayoutQ;
                       disps: openArray[ptr cArray[cint]]; ndisps: cint) =
  var subs = newSeq[cstring](ndisps)
  var i: cint = 0
  let s = "all"
  while i < ndisps:
    subs[i] = cstring s
    inc(i)
  makeGDFromShiftSubs(gd, l, disps, subs, ndisps)

proc makeShiftMultiSubQ*(si: openArray[ptr ShiftIndicesQ];
                         l: ptr LayoutQ; disp: openArray[ptr cArray[cint]];
                         subs: openArray[cstring]; ndisp: cint) =
  var myRank = l.myrank
  var nd = l.nDim
  var vvol = l.nSitesOuter
  var gi = cast[ptr GatherIndices](allocShared(sizeof(GatherIndices)))
  for n in 0..<ndisp:
    si[n].gi = gi
    si[n].disp = cast[type(si[n].disp)](allocShared(nd*sizeof((cint))))
    for i in 0..<nd:
      si[n].disp[i] = disp[n][i]
    si[n].pidx = cast[ptr cArray[cint]](allocShared(vvol*sizeof((cint))))
    si[n].sidx = cast[ptr cArray[cint]](allocShared(vvol*sizeof((cint))))
    #si[n].sendSites = cast[ptr cArray[cint]](allocShared(vvol*sizeof((cint))))
    for i in 0..<vvol:
      si[n].pidx[i] = -1
      si[n].sidx[i] = -1
  #mapmargs args;
  #args.l = l;
  #args.disp = disp;
  #args.ndisp = ndisp;
  #makeGather(gi, mapm, &args,l->nranks,l->nranks,l->nSites*ndisp,l->myrank);
  var gd = cast[ptr GatherDescription](allocShared(sizeof((GatherDescription))))
  makeGDFromShiftSubs(gd, l, disp, subs, ndisp)
  #makeGDFromShiftSubs(gd, l,
  #                    cast[ptr carray[ptr carray[cint]]](disp[0].unsafeaddr),
  #                    cast[ptr carray[cstring]](subs[0].unsafeaddr), ndisp)
  makeGatherFromGD(gi, gd)
  var si0 = ndisp-1  # si index for send structures
  var
    perm: cint = 0
    pack: cint = 0
    packs = [0,0,0,0]
    packbits = [0,0,0,0]
    lsrcloc = newSeq[int32](l.nSites)
  si[si0].sendSites.newSeq(0)
  si[si0].sbufcount.newSeq(0)
  si[si0].lbufcount.newSeq(0)
  if gi.nSendIndices > 0:
    for i in 0..<lsrcloc.len: lsrcloc[i] = -1
    pack = gi.sendIndices[0] mod l.nSitesInner
    if pack == 0:
      var i: cint = 1
      while (i < gi.nSendIndices) and (i < l.nSitesInner) and
          (gi.sendIndices[i] == gi.sendIndices[0] + i):
        inc(i)
      pack = - (i mod l.nSitesInner)
    var ssi0 = gi.sendIndices[0] div l.nSitesInner
    var pck = 0
    var pckbits = 0
    var scount = 0
    var lcount = 0
    for i in 0..<gi.nSendIndices:
      let ss = gi.sendIndices[i]
      let ssi = ss div l.nSitesInner
      let ssv = ss mod l.nSitesInner
      let pckn = 1 shl ssv
      pck += pckn
      inc pckbits
      let ssi1 = if i+1<gi.nSendIndices: gi.sendIndices[i+1] div l.nSitesInner else: -1
      if ssi1!=ssi0:
        ssi0 = ssi1
        var imask = 0
        while true:
          if imask >= packs.len:
            echo "pck ", pck, " not in packs ", packs
            qexError()
          if pck == packs[imask]: break
          if packs[imask] == 0:
            packs[imask] = pck
            packbits[imask] = pckbits
            break
          inc imask
        si[si0].sendSites.add SendSite(maskidx:uint32 imask, site:uint32 ssi)
        si[si0].sbufcount.add int32 scount
        si[si0].lbufcount.add int32 lcount
        scount += pckbits
        #lcount += l.nSitesInner - pckbits
        for iv in 0..<l.nSitesInner:
          if not pck.testBit(iv):
            let locidx = ssi*l.nSitesInner + iv
            lsrcloc[locidx] = int32 lcount
            inc lcount
        pck = 0
        pckbits = 0
    #echo "packs: ", packs, "  ", packbits
  si[si0].packmasks = packs
  si[si0].packbits = packbits
  for i in 0..<ndisp:
    si[i].nSendRanks = 0
    si[i].nSendSites = 0
    si[i].nSendSites1 = 0
  si[si0].nSendRanks = gi.nSendRanks
  si[si0].nSendSites = si[si0].sendSites.len
  si[si0].nSendSites1 = gi.nSendIndices
  if gi.nSendRanks > 0:
    si[si0].sendRanks = gi.sendRanks
    #si[si0].sendRankSizes = cast[ptr cArray[int]](allocShared(gi.nSendRanks*sizeof(int)))
    si[si0].sendRankSizes1 = gi.sendRankSizes
    #si[si0].sendRankOffsets = cast[ptr cArray[cint]](allocShared(gi.nSendRanks*sizeof(cint)))
    si[si0].sendRankOffsets1 = gi.sendRankOffsets
    #si[si0].sendRankSizes[0] = si[si0].sendSites.len
    #si[si0].sendRankOffsets[0] = 0
    #if gi.sendRankSizes[0] != gi.nSendIndices:
    #  echo &"gi.sendRankSizes[0] {gi.sendRankSizes[0]} != gi.nSendIndices {gi.nSendIndices}"
    #  qexError()
    #if gi.sendRankOffsets[0] != 0:
    #  echo &"gi.sendRankOffsets[0] {gi.sendRankOffsets[0]} != 0"
    #  qexError()

  var
    nrsites: cint = 0
    nrdests = newSeq[cint](ndisp)
  for i in 0..<ndisp:
    nrdests[i] = 0
  for i in 0 ..< vvol * ndisp:
    var dd = i div l.nSitesOuter
    var ix = i mod l.nSitesOuter
    var k0 = i * l.nSitesInner
    var recv = 0
    var rbi = 0
    for ii in 0 ..< l.nSitesInner:
      var k = k0 + ii
      var s = gi.srcIndices[k]
      if s == -1:
        recv = -1
        break
      if s < 0:
        inc(recv)
        if rbi == 0: rbi = s
    if recv < 0:
      si[dd].pidx[ix] = -1
      si[dd].sidx[ix] = -1
    elif recv == 0:
      si[dd].pidx[ix] = gi.srcIndices[k0] div l.nSitesInner
      si[dd].sidx[ix] = gi.srcIndices[k0] div l.nSitesInner
      var p = gi.srcIndices[k0] mod l.nSitesInner
      if p != 0:
        perm = p
        si[dd].pidx[ix] = - (si[dd].pidx[ix]) - 2
    else:
      rbi = - (rbi + 2)
      rbi = (2 * rbi) div l.nSitesInner
      if pack == 0: rbi = rbi div 2
      si[dd].sidx[ix] = int32(- rbi - 2)
      #nrsites++;
      inc(nrdests[dd])
  nrsites = gi.recvSize div l.nSitesInner
  if pack != 0: nrsites = nrsites * 2
  for i in 0..<ndisp:
    si[i].nRecvRanks = 0
    si[i].nRecvSites1 = 0
  si[0].nRecvRanks = gi.nRecvRanks
  #si[0].nRecvSites = nrsites
  si[0].nRecvSites1 = gi.recvSize
  if gi.nRecvRanks > 0:
    si[0].recvRanks = gi.recvRanks
    #si[0].recvRankSizes = cast[ptr cArray[cint]](allocShared(si[0].nRecvRanks*sizeof(cint)))
    si[0].recvRankSizes1 = gi.recvRankSizes
    #si[0].recvRankOffsets = cast[ptr cArray[cint]](allocShared(si[0].nRecvRanks * sizeof(cint)))
    si[0].recvRankOffsets1 = gi.recvRankOffsets
    #si[0].recvRankSizes[0] = nrsites
    #si[0].recvRankOffsets[0] = 0
    #if gi.recvRankSizes[0] != gi.recvSize:
    #  echo &"gi.recvRankSizes[0] {gi.recvRankSizes[0]} != gi.recvSize {gi.recvSize}"
    #  qexError()
    #if gi.recvRankOffsets[0] != 0:
    #  echo &"gi.recvRankOffsets[0] {gi.recvRankOffsets[0]} != 0"
    #  qexError()
  for n in 0..<ndisp:
    si[n].nRecvDests = nrdests[n]
    if nrdests[n] > 0:
      si[n].recvDests = cast[ptr cArray[cint]](allocShared(nrdests[n]*sizeof(cint)))
      si[n].recvLocalSrcs = cast[ptr cArray[cint]](allocShared(nrdests[n]*sizeof(cint)))
      si[n].recvRemoteSrcs = cast[ptr cArray[cint]](allocShared(nrdests[n]*sizeof(cint)))
      var j = 0
      for i in 0..<vvol:
        if si[n].sidx[i] < -1:
          var k = - (si[n].sidx[i] + 2)
          si[n].recvDests[j] = int32 i
          si[n].recvRemoteSrcs[j] = k
          si[n].recvLocalSrcs[j] = 0
          for i0 in 0..<l.nSitesInner:
            var ii = n * l.nSites + i * l.nSitesInner + i0
            var gs = gi.srcIndices[ii]
            if gs >= 0:
              si[n].recvLocalSrcs[j] = gs div l.nSitesInner
              break
          inc(j)
          if j > nrdests[n]:
            echo "j($#)>nrdests[$#]($#)"%[$j,$n,$nrdests[n]]
    si[n].vv = int32 vvol
    si[n].perm = perm
    si[n].pack = pack
    si[n].blend = pack
  #[
  var rboffs = newSeq[int]()
  for i in 0..<gi.nRecvDests:
  var rbo = -1
  var rbc = 0
  for i in 0..<gi.nRecvDests:
    let rdi = gi.recvDestIndices[i]
    let rdiv = rdi div l.nSitesInner
    let rdiv1 = if i+1<gi.nRecvDests: recvDestIndices[i+1] else: -1
    if rbo < 0: rbo = gi.recvBufIndices[i]
    if rdiv1 != rdiv:
      let n = rdiv div l.nSitesOuter
      let ix = rdiv mod l.nSitesOuter
      rboffs[n].add rbo
      si[n].lbufOffset.add si[si0].lbufcount[rbc]
      inc rbc
      si[n].rbufOffset.add rbo
      rbo = -1
  ]#
  var ridx = newSeqOfCap[int32](l.nSitesInner)
  var lidx = newSeqOfCap[int32](l.nSitesInner)
  for n in 0..<ndisp:
    si[n].recvIndex.newSeq(vvol)
    var nrecv = 0
    si[n].recvmasks = [-1,-1]
    for io in 0..<vvol:
      ridx.setLen(0)
      lidx.setLen(0)
      var rmask = 0
      var rbit = 1
      var rbits = 0
      var k0 = (n*vvol+io) * l.nSitesInner
      for ii in 0 ..< l.nSitesInner:
        var k = k0 + ii
        var s = gi.srcIndices[k]
        if s == -1:
          rbits = -1
          break
        if s < 0:
          let rb = -2 - s
          ridx.add rb
          rmask += rbit
          inc rbits
        else:
          lidx.add lsrcloc[s]
        rbit *= 2
      if rbits < 0:
        discard
      elif rbits == 0:
        discard
      else:
        for k in 0..<ridx.len:
          if ridx[k] != ridx[0] + k:
            echo ridx
            qexError("recv buf indices not contiguous")
        for k in 0..<lidx.len:
          if lidx[k] != lidx[0] + k:
            echo lidx
            qexError("local buf indices not contiguous")
        #si[n].recvDests[nrecv] = int32 s
        if si[n].recvDests[nrecv] != io:
          echo &"si[n].recvDests[nrecv] {si[n].recvDests[nrecv]} != io {io}"
          qexError()
        si[n].recvRemoteSrcs[nrecv] = ridx[0]
        #if si[n].recvRemoteSrcs[nrecv] != ridx[0]:
        #  echo &"si[n].recvRemoteSrcs[nrecv] {si[n].recvRemoteSrcs[nrecv]} != ridx[0] {ridx[0]}"
        #  qexError()
        if lidx.len>0:
          si[n].recvLocalSrcs[nrecv] = lidx[0]
        #  if si[n].recvLocalSrcs[nrecv] != lidx[0]:
        #    echo &"si[n].recvLocalSrcs[nrecv] {si[n].recvLocalSrcs[nrecv]} != lidx[0] {lidx[0]}"
        #    qexError()
        #if (rmask and si[n].packmasks[0]) != 0:
        #  echo &"rmask {rmask} and si[n].packmasks[0] {si[n].packmasks[0]} != 0"
        #  qexError()
        var imask = 0
        while true:
          if imask >= si[n].recvmasks.len:
            echo &"si[n].recvmasks[{imask}] {si[n].recvmasks[imask]} != rmask {rmask}"
            qexError()
          if rmask == si[n].recvmasks[imask]: break
          if si[n].recvmasks[imask] < 0:
            si[n].recvmasks[imask] = rmask
            si[n].recvbits[imask] = rbits
            break
          inc imask
        #if (-2-si[n].sidx[io]) != nrecv: # reuse sidx for now, replace later
        #  echo &"-2-si[n].sidx[io] {-2-si[n].sidx[io]} != nrecv {nrecv}"
        #  qexError()
        #si[n].sidx[io] = int32 -2 - nrecv
        si[n].recvIndex[io] = RecvIdx(maskidx: uint32 imask, idx: uint32 nrecv)
        inc nrecv

proc makeShiftMultiQ*(si: openArray[ptr ShiftIndicesQ]; l: ptr LayoutQ;
                      disp: openArray[ptr cArray[cint]]; ndisp: cint) =
  var subs = newSeq[cstring](ndisp)
  var i: cint = 0
  while i < ndisp:
    subs[i] = "all"
    inc(i)
  makeShiftMultiSubQ(si, l, disp, subs, ndisp)

proc makeShiftQ*(si: ptr ShiftIndicesQ; l: ptr LayoutQ;
                 disp: ptr cArray[cint]) =
  makeShiftMultiQ([si], l, [disp], 1)

proc makeShiftSubQ*(si: ptr ShiftIndicesQ; l: ptr LayoutQ;
                    disp: ptr cArray[cint]; sub: cstring) =
  makeShiftMultiSubQ([si], l, [disp], [sub], 1)
