#const hdr = currentSourcePath()[0..^12] & "qlayout.h"
#{. pragma: ql, header:hdr .}
#{. passC:"-I." .}
#{. compile:"qlayout.c" .}
#{. compile:"qshifts.c" .}
#{. compile:"qgather.c" .}

import base
import strutils
import algorithm

type
  GatherDescription* = object
    myRank*: cint
    nIndices*: cint
    srcRanks*: ptr cArray[cint]
    srcIndices*: ptr cArray[cint]
    nRecvDests*: cint
    nSendIndices*: cint
    sendSrcIndices*: ptr cArray[cint]
    sendDestRanks*: ptr cArray[cint]
    sendDestIndices*: ptr cArray[cint]

  GatherIndices* = object
    gd*: ptr GatherDescription
    myRank*: cint
    nIndices*: cint
    srcIndices*: ptr cArray[cint]
    nRecvRanks*: cint
    recvRanks*: ptr cArray[cint]
    recvRankSizes*: ptr cArray[cint]
    recvRankOffsets*: ptr cArray[cint]
    recvSize*: cint
    nRecvDests*: cint
    recvDestIndices*: ptr cArray[cint]
    recvBufIndices*: ptr cArray[cint]
    nSendRanks*: cint
    sendRanks*: ptr cArray[cint]
    sendRankSizes*: ptr cArray[cint]
    sendRankOffsets*: ptr cArray[cint]
    sendSize*: cint            ##  same as nSendIndices
    nSendIndices*: cint
    sendIndices*: ptr cArray[cint]

  GatherMap* = proc(srcRank: ptr cint; srcIdx: ptr cint; dstRank: cint;
                    dstIdx: ptr cint; args: pointer)

#proc makeGatherFromGD*(gi: ptr GatherIndices; gd: ptr GatherDescription) {.
#  importc,header:hdr.}

# map(&sr,&si,dr>=0,&di>=0) -> sr,si
# map(&sr>=0,&si,dr>=0,&di0<0) ->
#   si>=0 then (sr,si)->(dr,di) and di>=-(di0+1) is smallest such di
#   si<0 if sr doesn't send to dr
# each dest site has zero or one source sites
# each source site can have any number of dest sites
# map(&sr,&si,&dr>=0,&di>=0,dn<0) -> sr,si
# map(&sr>=0,&si>=0,&dr,&di,dn>=0) -> dr,di for dest number dn
# GatherDescription
# pass in myRank, nIndices, srcRanks, srcIndices, 
#  nSendIndices, sendSrcIndices, sendDestRanks, sendDestIndices (not sorted)
# merge (add?)  + (need destIndexOffset per rank/gather)
# compose       *
# filter (on dest sites)
# -> getSendList (1 for multi)
# -> getRecvList (1 for multi)
# -> getDestList (N for multi)
# -> getSrcList (N for multi) (needs recv sites)
# getRecvInfo
# getSendInfo
#  multi

proc mergeGatherDescriptions*(gd: ptr GatherDescription;
                              gds: ptr cArray[GatherDescription]; n: cint) =
  var
    ni: cint = 0
    nsi: cint = 0
  var i: cint = 0
  while i < n:
    inc(ni, gds[i].nIndices)
    inc(nsi, gds[i].nSendIndices)
    if gds[i].myRank != gds[0].myRank:
      var myRank: cint = 0
      echo "ranks don\'t match: gds[$#].myRank($#)!=gds[0].myRank($#)"%
                  [$i, $gds[i].myRank, $gds[0].myRank]
      quit(-1)
    inc(i)
  var sr = cast[ptr carray[cint]](alloc(ni * sizeof(cint)))
  var si = cast[ptr carray[cint]](alloc(ni * sizeof(cint)))
  ni = 0
  i=0
  while i < n:
    var j: cint = 0
    while j < gds[i].nIndices:
      sr[ni] = gds[i].srcRanks[j]
      si[ni] = gds[i].srcIndices[j]
      inc(ni)
      inc(j)
    inc(i)
  var ssi = cast[ptr carray[cint]](alloc(nsi * sizeof(cint)))
  var sdr = cast[ptr carray[cint]](alloc(nsi * sizeof(cint)))
  var sdi = cast[ptr carray[cint]](alloc(nsi * sizeof(cint)))
  nsi = 0
  i=0
  while i < n:
    var j: cint = 0
    while j < gds[i].nSendIndices:
      ssi[nsi] = gds[i].sendSrcIndices[j]
      sdr[nsi] = gds[i].sendDestRanks[j]
      sdi[nsi] = gds[i].sendDestIndices[j]
      inc(nsi)
      inc(j)
    inc(i)
  gd.myRank = gds[0].myRank
  gd.nIndices = ni
  gd.srcRanks = sr
  gd.srcIndices = si
  gd.nSendIndices = nsi
  gd.sendSrcIndices = ssi
  gd.sendDestRanks = sdr
  gd.sendDestIndices = sdi

proc cyclicComp*(a: cint; b: cint; zero: cint): cint =
  var c: cint = abs(a) + abs(b) + 1
  result = a - b
  if a < zero: result += c
  if b < zero: result -= c

type Ssd = tuple[a:cint,b:cint,c:cint]
var g_gd*: ptr GatherDescription
var g_sd*: seq[Ssd]

proc sortSd*(a: Ssd; b: Ssd): int =
  #var pa = cast[ptr carray[cint]](a)
  #var pb = cast[ptr carray[cint]](b)
  var mr = g_gd.myRank
  var ra = a[0]
  var rb = b[0]
  var rr = cyclicComp(ra, rb, mr)
  if rr == 0:
    var ia = a[2]
    var ib = b[2]
    rr = ia - ib
  return rr

proc sortSd2*(a: cint; b: cint): int =
  #var pa: cint = 3 * (cast[ptr cint](a)[])
  #var pb: cint = 3 * (cast[ptr cint](b)[])
  var mr: cint = g_gd.myRank
  var ra: cint = g_sd[a].a
  var rb: cint = g_sd[b].a
  var rr: cint = cyclicComp(ra, rb, mr)
  if rr == 0:
    var ia: cint = g_sd[a].b
    var ib: cint = g_sd[b].b
    rr = ia - ib
    if rr == 0:
      var ja: cint = g_sd[a].c
      var jb: cint = g_sd[b].c
      rr = ja - jb
  return rr

# uses:  myRank, nIndices, srcRanks, srcIndices, 

proc makeRecvInfo*(gi: ptr GatherIndices; gd: ptr GatherDescription) =
  var mr: cint = gd.myRank
  var ni: cint = gd.nIndices
  gi.gd = gd
  gi.myRank = mr
  gi.nIndices = ni
  gi.srcIndices = nil
  gi.nRecvRanks = 0
  gi.recvRanks = nil
  gi.recvRankSizes = nil
  gi.recvRankOffsets = nil
  gi.recvSize = 0
  gi.nRecvDests = 0
  gi.recvDestIndices = nil
  gi.recvBufIndices = nil
  if ni == 0: return
  when 1:
    # sort: srcRank, destIndex
    # perm(p): srcRank, srcIndex, destIndex  : a[p[i]] < a[p[j]] (i<j)
    # add firstDestIndex[p[i]] = p[i-?]
    # rdi[i] = destIndex[i]
    # if(firstDestIndex
    # rbi[i] = j
    # sd: srcRanks, srcIndices, destIndices
    #var sd = cast[ptr carray[cint]](alloc(3 * gd.nRecvDests * sizeof(cint)))
    var si = cast[ptr carray[cint]](alloc(ni * sizeof(cint)))
    var nsd: cint = 0
    var sd = newSeq[Ssd](gd.nRecvDests)
    ##pragma omp parallel for
    var i: cint = 0
    while i < ni:
      var r: cint = gd.srcRanks[i]
      if r < 0:
        si[i] = -1
      elif r == mr:
        #local
        si[i] = gd.srcIndices[i]
      else:
        # remote
        ##pragma omp critical
        #sd[nsd] = r
        #sd[nsd + 1] = gd.srcIndices[i]
        #sd[nsd + 2] = i
        #inc(nsd, 3)
        sd[nsd] = (r, gd.srcIndices[i], i)
        nsd += 1
      inc(i)
    #nsd = nsd div 3
    g_gd = gd
    #qsort(sd, nsd, 3 * sizeof((cint)), sortSd)
    sd.setLen(nsd)
    sort(sd, sortSd)
    #var p = cast[ptr cint = myalloc(nsd * sizeof((int)))
    var p = newSeq[cint](nsd)
    i=0
    while i < nsd:
      p[i] = i
      inc(i)
    g_sd = sd
    #qsort(p, nsd, sizeof((int)), sortSd2)
    sort(p, sortSd2)
    #ARRAY_CREATE(int, recvRanks)
    #ARRAY_CREATE(int, recvRankCounts)
    var recvRanks = newSeq[cint]()
    var recvRankCounts = newSeq[cint]()
    # check for duplicate source indices
    block:
      var rr: cint = -1
      var rrc: cint = 0
      var k0: cint = -1
      var r0: cint = -1
      var i0: cint = -1
      i=0
      while i < nsd:
        var k: cint = p[i]
        #var k3: cint = 3 * k
        #var sri: cint = sd[k3]
        #var sii: cint = sd[k3 + 1]
        var sri = sd[k].a
        var sii = sd[k].b
        if sri != r0 or sii != i0:
          k0 = k
          r0 = sri
          i0 = sii
          if r0 != rr:
            if i > 0:
              #ARRAY_APPEND(int, recvRankCounts, rrc)
              recvRankCounts.add rrc
            #ARRAY_APPEND(int, recvRanks, r0)
            recvRanks.add r0
            rr = r0
            rrc = 0
          inc(rrc)
        #sd[k3 + 1] = k0
        sd[k].b = k0
        inc(i)
      #ARRAY_APPEND(int, recvRankCounts, rrc)
      recvRankCounts.add rrc
    var rsize: cint = 0
    var rr = cast[ptr carray[cint]](alloc(recvRanks.len * sizeof(cint)))
    var rrs = cast[ptr carray[cint]](alloc(recvRanks.len * sizeof(cint)))
    var rro = cast[ptr carray[cint]](alloc(recvRanks.len * sizeof(cint)))
    i=0
    while i < recvRanks.len:
      rr[i] = recvRanks[i]
      rrs[i] = recvRankCounts[i]
      rro[i] = rsize
      inc(rsize, recvRankCounts[i])
      inc(i)
    # finish si, rbuf
    var rdi = cast[ptr carray[cint]](alloc(nsd * sizeof(cint)))
    var rbi = cast[ptr carray[cint]](alloc(nsd * sizeof(cint)))
    ##pragma omp parallel for
    var j: cint = 0
    i=0
    while i < nsd:
      #var i3: cint = 3 * i
      #var di: cint = sd[i3 + 2]
      #var fdi: cint = sd[i3 + 1]
      var di = sd[i][2]
      var fdi = sd[i][1]
      var bi: cint = j
      if fdi == i: inc(j)
      else:
        #bi = sd[3 * fdi + 1]
        bi = sd[fdi][1]
      #sd[i3 + 1] = bi
      sd[i][1] = bi
      rdi[i] = di
      rbi[i] = bi
      si[di] = - bi - 2
      inc(i)
    gi.srcIndices = si
    gi.nRecvRanks = cint recvRanks.len
    gi.recvRanks = rr
    gi.recvRankSizes = rrs
    gi.recvRankOffsets = rro
    gi.recvSize = rsize
    gi.nRecvDests = nsd
    gi.recvDestIndices = rdi
    gi.recvBufIndices = rbi
    #free(sd)
    #free(p)
    #free(recvRanks)
    #free(recvRankCounts)
  #[
  else: 
    ARRAY_CREATE(int, recvRanks)
    ARRAY_CREATE(int, recvRankCounts)
    ARRAY_CREATE(int, rbufSrcRankIndices)
    ARRAY_CREATE(int, rbufSrcIndices)
    ARRAY_CREATE(int, rbufRankIndices)
    ARRAY_CREATE(int, recvDestIndices)
    ARRAY_CREATE(int, recvBufIndices)
    var si: ptr cint = myalloc(ni * sizeof((int)))
    var i: cint = 0
    while i < ni: 
      var r: cint = gd.srcRanks[i]
      if r < 0: 
        si[i] = - 1
      elif r == mr: 
        #local
        si[i] = gd.srcIndices[i]
      else: 
        # remote
        # check if new rank
        var ri: cint = 0
        while ri < nrecvRanks and r != recvRanks[ri]: inc(ri)
        if ri == nrecvRanks: 
          ARRAY_APPEND(int, recvRanks, r)
          ARRAY_APPEND(int, recvRankCounts, 0)
        var sii: cint = gd.srcIndices[i]
        var j: cint = 0
        while j < nrbufSrcRankIndices and
            (rbufSrcRankIndices[j] != ri or rbufSrcIndices[j] != sii): 
          inc(j)
        # if not found, add to rbuf
        if j == nrbufSrcRankIndices: 
          ARRAY_APPEND(int, rbufSrcRankIndices, ri)
          ARRAY_APPEND(int, rbufSrcIndices, sii)
          ARRAY_APPEND(int, rbufRankIndices, recvRankCounts[ri])
          inc(recvRankCounts[ri])
        ARRAY_APPEND(int, recvDestIndices, i)
        #ARRAY_APPEND(int, recvBufIndices, rbufRankIndices[j]);
        ARRAY_APPEND(int, recvBufIndices, j)
      inc(i)
    # sort ranks with indexing array
    var p: ptr cint = myalloc(nrecvRanks * sizeof((int)))
    var i: cint = 0
    while i < nrecvRanks: 
      p[i] = i
      inc(i)
    var i: cint = 0
    while i < nrecvRanks: 
      var ri: cint = i
      var rv: cint = recvRanks[p[i]]
      var j: cint = i + 1
      while j < nrecvRanks: 
        if cyclicComp(recvRanks[p[j]], rv, mr) <= 0: 
          ri = j
          rv = recvRanks[p[j]]
        inc(j)
      var l: cint = p[i]
      p[i] = p[ri]
      p[ri] = l
      inc(i)
    var pinv: ptr cint = myalloc(nrecvRanks * sizeof((int)))
    var i: cint = 0
    while i < nrecvRanks: 
      pinv[p[i]] = i
      inc(i)
    var rsize: cint = 0
    var rr: ptr cint = myalloc(nrecvRanks * sizeof((int)))
    var rrs: ptr cint = myalloc(nrecvRanks * sizeof((int)))
    var rro: ptr cint = myalloc(nrecvRanks * sizeof((int)))
    var i: cint = 0
    while i < nrecvRanks: 
      var ri: cint = p[i]
      rr[i] = recvRanks[ri]
      rrs[i] = recvRankCounts[ri]
      rro[i] = rsize
      inc(rsize, recvRankCounts[ri])
      inc(i)
    # finish si, rbuf
    var rdi: ptr cint = myalloc(nrecvDestIndices * sizeof((int)))
    var rbi: ptr cint = myalloc(nrecvDestIndices * sizeof((int)))
    var i: cint = 0
    while i < nrecvDestIndices: 
      var di: cint = recvDestIndices[i]
      var j: cint = recvBufIndices[i]
      var ri: cint = rbufSrcRankIndices[j]
      var bi: cint = rbufRankIndices[j]
      rdi[i] = di
      rbi[i] = rro[pinv[ri]] + bi
      si[di] = - rbi[i] - 2
      inc(i)
    #gi->gd = gd;
    #gi->myRank = mr;
    #gi->nIndices = ni;
    gi.srcIndices = si
    gi.nRecvRanks = nrecvRanks
    gi.recvRanks = rr
    gi.recvRankSizes = rrs
    gi.recvRankOffsets = rro
    gi.recvSize = rsize
    gi.nRecvDests = nrecvDestIndices
    gi.recvDestIndices = rdi
    gi.recvBufIndices = rbi
    free(recvRanks)
    free(recvRankCounts)
    free(rbufSrcRankIndices)
    free(rbufSrcIndices)
    free(rbufRankIndices)
    free(recvDestIndices)
    free(recvBufIndices)
    free(p)
    free(pinv)
  ]#

proc sortsendSrc*(a: cint; b: cint): int =
  #var pa: cint = cast[ptr cint](a)[]
  #var pb: cint = cast[ptr cint](b)[]
  var mr = g_gd.myRank
  var ra = g_gd.sendDestRanks[a]
  var rb = g_gd.sendDestRanks[b]
  var rr = cyclicComp(ra, rb, mr)
  if rr == 0:
    var ia = g_gd.sendSrcIndices[a]
    var ib = g_gd.sendSrcIndices[b]
    rr = ia - ib
    if rr == 0:
      var ja = g_gd.sendDestIndices[a]
      var jb = g_gd.sendDestIndices[b]
      rr = ja - jb
  return rr

proc sortsend*(a: cint; b: cint): int =
  #var pa: cint = cast[ptr cint](a)[]
  #var pb: cint = cast[ptr cint](b)[]
  if a<0:
    if b<0: return 0
    else: return 1
  if b<0: return -1
  var mr = g_gd.myRank
  var ra = g_gd.sendDestRanks[a]
  var rb = g_gd.sendDestRanks[b]
  var rr = cyclicComp(ra, rb, mr)
  if rr == 0:
    var ia = g_gd.sendDestIndices[a]
    var ib = g_gd.sendDestIndices[b]
    rr = ia - ib
  return rr

# uses: myRank, nSendIndices, sendSrcIndices, sendDestRanks, sendDestIndices

proc makeSendInfo*(gi: ptr GatherIndices; gd: ptr GatherDescription) =
  gi.nSendRanks = 0
  gi.sendRanks = nil
  gi.sendRankSizes = nil
  gi.sendRankOffsets = nil
  gi.sendSize = 0
  gi.nSendIndices = 0
  gi.sendIndices = nil
  var n: cint = gd.nSendIndices
  if n == 0: return
  when true:
    #var p: ptr cint = myalloc(n * sizeof((int)))
    var p = newSeq[cint](n)
    var i: cint = 0
    while i < n:
      p[i] = i
      inc(i)
    g_gd = gd
    #qsort(p, n, sizeof((int)), sortsendSrc)
    sort(p, sortSendSrc)
    #ARRAY_CREATE(int, sendRanks)
    #ARRAY_CREATE(int, sendRankCounts)
    var sendRanks = newSeq[cint]()
    var sendRankCounts = newSeq[cint]()
    var ndup: cint = 0
    block:
      var sr: cint = -1
      var src: cint = 0
      var r0: cint = -1
      var i0: cint = -1
      i=0
      while i < n:
        var pi: cint = p[i]
        var sdr: cint = gd.sendDestRanks[pi]
        var ssi: cint = gd.sendSrcIndices[pi]
        if sdr != r0 or ssi != i0:
          r0 = sdr
          i0 = ssi
          if r0 != sr:
            if i > 0:
              #ARRAY_APPEND(int, sendRankCounts, src)
              sendRankCounts.add src
            #ARRAY_APPEND(int, sendRanks, r0)
            sendRanks.add r0
            sr = r0
            src = 0
          inc(src)
        else:
          p[i] = p[ndup]
          p[ndup] = -1
          inc(ndup)
        inc(i)
      #ARRAY_APPEND(int, sendRankCounts, src)
      sendRankCounts.add src
    var ssize: cint = 0
    var sr = cast[ptr carray[cint]](alloc(sendRanks.len * sizeof(cint)))
    var srs = cast[ptr carray[cint]](alloc(sendRanks.len * sizeof(cint)))
    var sro = cast[ptr carray[cint]](alloc(sendRanks.len * sizeof(cint)))
    i=0
    while i < sendRanks.len:
      sr[i] = sendRanks[i]
      srs[i] = sendRankCounts[i]
      sro[i] = ssize
      inc(ssize, sendRankCounts[i])
      inc(i)
    var nsend: cint = n - ndup
    #qsort(p + ndup, nsend, sizeof((int)), sortsend)
    sort(p, sortsend)
    var si = cast[ptr carray[cint]](alloc(nsend * sizeof(cint)))
    i=0
    while i < nsend:
      #var k: cint = p[ndup + i]
      var k = p[i]
      si[i] = gd.sendSrcIndices[k]
      inc(i)
    gi.nSendRanks = cint sendRanks.len
    gi.sendRanks = sr
    gi.sendRankSizes = srs
    gi.sendRankOffsets = sro
    gi.sendSize = nsend
    gi.nSendIndices = nsend
    gi.sendIndices = si
    #free(p)
    #free(sendRanks)
    #free(sendRankCounts)
  #[
  else: 
    var p: ptr cint = myalloc(n * sizeof((int)))
    var i: cint = 0
    while i < n: 
      p[i] = i
      inc(i)
    var mr: cint = gd.myRank
    var i: cint = 0
    while i < n: 
      var k: cint = i
      var sdr: cint = gd.sendDestRanks[p[i]]
      var sdi: cint = gd.sendDestIndices[p[i]]
      var j: cint = i + 1
      while j < n: 
        var rj: cint = gd.sendDestRanks[p[j]]
        var ij: cint = gd.sendDestIndices[p[j]]
        if ((rj == sdr) and (ij < sdi)) or (cyclicComp(rj, sdr, mr) < 0): 
          k = j
          sdr = rj
          sdi = ij
        inc(j)
      var l: cint = p[i]
      p[i] = p[k]
      p[k] = l
      inc(i)
    ARRAY_CREATE(int, sendRanks)
    ARRAY_CREATE(int, sendRankSizes)
    ARRAY_CREATE(int, sendRankOffsets)
    ARRAY_CREATE(int, sendIndices)
    ARRAY_APPEND(int, sendRanks, gd.sendDestRanks[p[0]])
    ARRAY_APPEND(int, sendRankOffsets, 0)
    var i: cint = 0
    while i < n: 
      var r: cint = gd.sendDestRanks[p[i]]
      if r != sendRanks[nsendRanks - 1]: 
        ARRAY_APPEND(int, sendRankSizes, 
                     nsendIndices - sendRankOffsets[nsendRankOffsets - 1])
        ARRAY_APPEND(int, sendRanks, r)
        ARRAY_APPEND(int, sendRankOffsets, nsendIndices)
      var ssi: cint = gd.sendSrcIndices[p[i]]
      var j: cint = sendRankOffsets[nsendRankOffsets - 1]
      while j < nsendIndices and ssi != sendIndices[j]: inc(j)
      if j == nsendIndices: 
        ARRAY_APPEND(int, sendIndices, ssi)
      inc(i)
    ARRAY_APPEND(int, sendRankSizes, 
                 nsendIndices - sendRankOffsets[nsendRankOffsets - 1])
    gi.nSendRanks = nsendRanks
    gi.sendRanks = myalloc(nsendRanks * sizeof((int)))
    ARRAY_CLONE(int, gi.sendRanks, sendRanks)
    ARRAY_CLONE(int, gi.sendRankSizes, sendRankSizes)
    ARRAY_CLONE(int, gi.sendRankOffsets, sendRankOffsets)
    gi.sendSize = nsendIndices
    gi.nSendIndices = nsendIndices
    ARRAY_CLONE(int, gi.sendIndices, sendIndices)
    free(p)
    free(sendRanks)
    free(sendRankSizes)
    free(sendRankOffsets)
    free(sendIndices)
  ]#

proc makeGD*(gd: ptr GatherDescription; map: GatherMap; args: pointer;
             nSrcRanks: cint; nDstRanks: cint; myndi: cint; myRank: cint) =
  var sidx = cast[ptr carray[cint]](alloc(myndi * sizeof(cint)))
  var srank = cast[ptr carray[cint]](alloc(myndi * sizeof(cint)))
  # find shift sources
  var di: cint = 0
  while di < myndi:
    var
      sr: cint
      si: cint
      di0: cint = di
    map(addr(sr), addr(si), myRank, addr(di0), args)
    srank[di] = sr
    sidx[di] = si
    inc(di)
  gd.myRank = myRank
  gd.nIndices = myndi
  gd.srcRanks = srank
  gd.srcIndices = sidx
  #ARRAY_CREATE(int, sendSrcIndices)
  #ARRAY_CREATE(int, sendDestRanks)
  #ARRAY_CREATE(int, sendDestIndices)
  var sendSrcIndices = newSeq[cint]()
  var sendDestRanks = newSeq[cint]()
  var sendDestIndices = newSeq[cint]()
  # find who to send to
  var dr: cint = 0
  while dr < nDstRanks:
    if dr == myRank: continue
    var
      sr: cint = myRank
      si: cint
      di: cint = -1
    map(addr(sr), addr(si), dr, addr(di), args)
    while si >= 0:
      #ARRAY_APPEND(int, sendSrcIndices, si)
      #ARRAY_APPEND(int, sendDestRanks, dr)
      #ARRAY_APPEND(int, sendDestIndices, di)
      sendSrcIndices.add si
      sendDestRanks.add dr
      sendDestIndices.add di
      di = - di - 2
      map(addr(sr), addr(si), dr, addr(di), args)
    inc(dr)
  gd.nSendIndices = cint sendSrcIndices.len
  template ARRAY_CLONE(x,y: typed) =
    x = cast[type(x)](alloc(y.len*sizeof(type(x[0]))))
    for i in 0..<y.len: x[i] = y[i]
  ARRAY_CLONE(gd.sendSrcIndices, sendSrcIndices)
  ARRAY_CLONE(gd.sendDestRanks, sendDestRanks)
  ARRAY_CLONE(gd.sendDestIndices, sendDestIndices)
  #free(sendSrcIndices)
  #free(sendDestRanks)
  #free(sendDestIndices)

proc makeGathersFromGDs*(gi: ptr carray[ptr GatherIndices];
                         gd: ptr carray[ptr GatherDescription]; n: cint) =
  var i: cint = 0
  while i < n:
    makeRecvInfo(gi[i], gd[i])
    makeSendInfo(gi[i], gd[i])
    inc(i)

proc makeGatherFromGD*(gi: ptr GatherIndices; gd: ptr GatherDescription) =
  makeRecvInfo(gi, gd)
  makeSendInfo(gi, gd)

#[
proc makeGather*(gi: ptr GatherIndices; map: ptr GatherMap; args: pointer;
                 nSrcRanks: cint; nDstRanks: cint; myndi: cint;
                 myRank: cint) =
  var gd: ptr GatherDescription = alloc(sizeof(GatherDescription))
  makeGD(gd, map, args, nSrcRanks, nDstRanks, myndi, myRank)
  makeRecvInfo(gi, gd)
  makeSendInfo(gi, gd)
]#
