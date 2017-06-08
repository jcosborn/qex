import base
#import threading
#import comms
import layout
#import matrixConcept
#import complexConcept
#import times
import macros
#import qcdTypes
#import qmp
#import stdUtils
import field
#import future
import strUtils
#import metaUtils

type ShiftB*[T] = object
  subset*: Subset
  si*: ShiftIndices
  sb*: ShiftBuf
  size*: int

template shiftBType*(x:SomeField):untyped = ShiftB[x[0].type]

template initShiftB*(s:ShiftB; l:Layout; t:typedesc;
                     dir,len:int; sub="all"):untyped =
  if threadNum==0:
    s.subset.layoutSubset(l, sub)
    s.si = l.getShift(dir, len, sub)
    s.size = sizeOf(t) div l.nSitesInner
    prepareShiftBuf(s.sb, s.si, s.size)
template initShiftB*(s:ShiftB; x:SomeField; dir,len:int; sub="all"):untyped =
  if threadNum==0:
    #template l:untyped = x.l
    s.subset.layoutSubset(x.l, sub)
    s.si = x.l.getShift(dir, len, sub)
    s.size = sizeOf(x[0]) div x.l.nSitesInner
    prepareShiftBuf(s.sb, s.si, s.size)

template createShiftB*(x:SomeField; dir,len:int; sub="all"):untyped =
  var s:ShiftB[x[0].type]
  s.initShiftB(x, dir,len, sub)
  s

template createGlobalShiftB*(v:untyped; x:SomeField;
                             dir,len:int; sub="all"):untyped =
  var v{.global.}:ShiftB[x[0].type]
  v.initShiftB(x, dir,len, sub)
  threadBarrier()

proc createShiftBufs*(n:int; x:any; ln=1; sub="all"):auto =
  var s = newSeq[ShiftB[x[0].type]](n)
  for i in 0..<n:
    s[i].initShiftB(x, i, ln, sub)
  result = s
proc createShiftBufs*(x:any; ln=1; sub="all"):auto =
  let n = x.l.nDim
  var s = newSeq[ShiftB[x[0].type]](n)
  for i in 0..<n:
    s[i].initShiftB(x, i, ln, sub)
  result = s

#proc init*(s:var ShiftB; ;
#           dir,len:int; sub="all") =

template startSB*(s:ShiftB; e:untyped) =
  mixin assign, `[]`, numberType
  if threadNum == 0:
    if s.si.nRecvRanks > 0:
      #echoRank "startRecvBuf"
      startRecvBuf(s.sb)
  if s.si.nSendRanks > 0:
    if s.si.pack == 0:
      let b = cast[ptr cArray[s.T]](s.sb.sq.sbuf)
      #echo "sendSites: ", s.si.nSendSites
      tFor i, 0..<s.si.nSendSites:
        let ix{.inject.} = s.si.sendSites[i]
        assign(b[i], e)
        #echoAll myrank, " ", i, " ", ix, " ", b[][i]
    else:
      type F = numberType(s.T)
      let stride = sizeOf(s.T) div (2*sizeof(F))
      let b = cast[ptr cArray[F]](s.sb.sq.sbuf)
      let l = cast[ptr cArray[F]](s.sb.lbuf)
      tFor i, 0..<s.si.nSendSites:
        let ix{.inject.} = s.si.sendSites[i]
        let j = stride * i
        pack(b[j].addr, l[j].addr, s.si.pack, e)
        #echoAll myrank, " ", i, " ", ix, " "
    t0wait()
    if threadNum == 0:
      #echoRank "send: ", cast[ptr float32](s.sb.sq.sbuf)[]
      startSendBuf(s.sb)

template isLocal*(s: ShiftB; i: int): bool =
  s.si.sq.pidx[i] != -1

template prefetchSB*(ss:ShiftB; ii:int; e1x:untyped):untyped =
  subst(s,ss,i,ii,e1,e1x):
    block:
      let imax = s.subset.highOuter
      if i<imax:
        let k1 = s.si.sq.pidx[i]
        if k1 >= 0:
          let ix{.inject.} = k1
          prefetch(addr(e1))
        elif k1 + 2 <= 0:
          let ix{.inject.} = -(k1 + 2)
          prefetch(addr(e1))

template localSB*(ss:ShiftB; ii:int; e1x,e2x:untyped):untyped {.dirty.} =
  bind subst
  subst(s,ss,i,ii,e1,e1x,e2,e2x,k1,_):
    block:
      let k1 = s.si.sq.pidx[i]
      var ix{.noInit.}: int
      if k1 >= 0:
        ix = k1
        template it: untyped = e2
        e1
      elif k1 + 2 <= 0:
        ix = -(k1 + 2)
        var it{.noInit.}: type(e2)
        perm(it, s.si.perm, e2)
        e1

proc boundaryOffsetSB*(s:ShiftB) =
  var ti0 = threadDivideLow(s.subset.lowOuter, s.subset.highOuter)
  var ti1 = threadDivideHigh(s.subset.lowOuter, s.subset.highOuter)
  var i0 = 0
  var step = (s.si.nRecvDests+1) div 2
  template search(i0,ti0:untyped):untyped =
    while true:
      while i0<s.si.nRecvDests and s.si.sq.recvDests[i0]<ti0:
        i0 += step
        step = (step+1) div 2
      if i0>s.si.nRecvDests: i0 = s.si.nRecvDests
      if i0==0 or s.si.sq.recvDests[i0-1]<ti0: break
      while i0>0 and s.si.sq.recvDests[i0-1]>=ti0:
        i0 -= step
        step = (step+1) div 2
      if i0<0: i0 = 0
  search(i0, ti0)
  var i1 = i0
  step = (i1+s.si.nRecvDests+1) div 2
  search(i1, ti1)
  s.sb.sq.offr[threadNum] = cint(i0)
  s.sb.sq.lenr[threadNum] = cint(i1)
  s.sb.sq.nthreads[threadNum] = numThreads.cint

template boundaryWaitSB*(s:ShiftB, e:untyped):untyped =
  if s.si.nRecvDests > 0:
    e
    if s.si.nRecvRanks > 0:
      if threadNum == 0:
        waitRecvBuf(s.sb)

template boundarySyncSB*():untyped =
  twait0()

template boundaryGetSB*(ss:ShiftB; irr:untyped; e:untyped):untyped =
  subst(s,ss,ir,irr,rr,_,i,_,k2,_,stride,_,itt,_):
    if s.si.nRecvDests > 0:
      let rr = cast[ptr cArray[s.T]](s.sb.sq.rbuf)
      if s.si.blend == 0:
        let i = -2 - s.si.sq.sidx[ir]
        let k2 = s.si.sq.recvRemoteSrcs[i]
        #echo "blend0: ", i, " ir: ", ir, " k2: ", k2
        subst(it,rr[k2]):
          e
      else:
        let stride = sizeof(s.T) div 2
        let i = -2 - s.si.sq.sidx[ir]
        let k2 = s.si.sq.recvRemoteSrcs[i]
        #echo "blendb: ", irr, " sidx: ", s.si.sq.sidx[irr].int
        var itt{.noInit.}: s.T  # should be load1(s.T)?
        blend(itt, s.sb.lbuf[stride*i].addr,
              s.sb.sq.rbuf[stride*k2].addr, s.si.blend)
        subst(it,itt):
          e

template getSB*(ss:ShiftB; ii:int; e1x,e2x:untyped):untyped =
  if isLocal(ss, ii):
    localSB(ss, ii, e1x, e2x)
  else:
    boundaryGetSB(ss, ii, e1x)

template boundarySB*(s:ShiftB; e:untyped):untyped =
  var needBoundary = false
  boundaryWaitSB(s): needBoundary = true
  if needBoundary:
    boundarySyncSB()
    if s.si.nRecvDests > 0:
      if s.sb.sq.nthreads[threadNum] != numThreads: boundaryOffsetSB(s)
      let ti0 = s.sb.sq.offr[threadNum]
      let ti1 = s.sb.sq.lenr[threadNum]
      if s.si.blend == 0:
        let rr = cast[ptr cArray[s.T]](s.sb.sq.rbuf)
        for i in ti0..<ti1:
          let irr = s.si.sq.recvDests[i]
          let k2 = s.si.sq.recvRemoteSrcs[i]
          #echo "blend0: ", i, " ir: ", irr, " k2: ", k2
          subst(ir,irr,it,rr[k2]):
            e
      else:
        let stride = sizeof(s.T) div 2
        for i in ti0..<ti1:
          let irr = s.si.sq.recvDests[i]
          let k2 = s.si.sq.recvRemoteSrcs[i]
          #echo "blendb: ", irr, " sidx: ", s.si.sq.sidx[irr].int
          var itt{.noInit.}: s.T  # should be load1(s.T)?
          blend(itt, s.sb.lbuf[stride*i].addr,
                s.sb.sq.rbuf[stride*k2].addr, s.si.blend)
          subst(ir,irr,it,itt):
            e
      if s.si.nRecvRanks > 0:
        if threadNum == 0:
          doneRecvBuf(s.sb)
    if s.si.nSendRanks > 0:
      if threadNum == 0:
        waitSendBuf(s.sb)

template shiftExpr*(sb: ShiftB, er, es) {.dirty.} =
  startSB(sb, es)
  for ir in sb.subset:
    localSB(sb, ir, er, es)
  boundarySB(sb, er)



type Shift*[V: static[int]; T] = object
  src*: Field[V,T]
  dest*: Field[V,T]
  dir*,len*: int
  sub*: string
  subset*: Subset
  si*: ShiftIndices
  sb*: ShiftBuf
  size*: int

template createShift*(x:Field):untyped =
  Shift[x.V,x.T](src:x)

proc init*[V:static[int],T](s:var Shift[V,T]; dest:Field[V,T];
                            dir,len:int; sub="all") =
  if threadNum == 0:
    let si = dest.l.getShift(dir, len, sub)
    s.dest = dest
    s.dir = dir
    s.len = len
    s.sub = sub
    s.subset.layoutSubset(dest.l, sub)
    s.si = si
    s.size = sizeOf(T) div dest.l.nSitesInner
    prepareShiftBuf(s.sb, s.si, s.size)
    #echo myrank, ": size: ", size
    #echo myrank, ": rbuf: ", cast[int](sb.sq.rbuf)

#template assign = discard
#template assign(x:any) = discard
proc start*[V:static[int],T](s:var Shift[V,T]; src:Field[V,T]) =
  mixin assign, numberType
  template si:untyped = s.si
  template sb:untyped = s.sb
  if threadNum == 0:
    s.src = src
    if si.nRecvRanks > 0:
      #echo myrank, ": startRecvBuf"
      startRecvBuf(sb)
  #echo "test3"
  if si.nSendRanks > 0:
    if si.pack == 0:
      let b = cast[ptr cArray[T]](sb.sq.sbuf)
      tFor i, 0..<si.nSendSites:
        let k = si.sendSites[i]
        assign(b[i], src[k])
    else:
      type F = numberType(T)
      let stride = sizeOf(T) div (2*sizeof(F))
      let b = cast[ptr cArray[F]](s.sb.sq.sbuf)
      let l = cast[ptr cArray[F]](s.sb.lbuf)
      tFor i, 0..<si.nSendSites:
        let k = si.sendSites[i]
        let j = stride * i
        #pack(sb.sq.sbuf[stride*i].addr, si.pack, src[k])
        pack(b[j].addr, l[j].addr, si.pack, src[k])
    t0wait()
    if threadNum == 0:
      #echo "startSendBuf"
      startSendBuf(sb)

proc local*(s:Shift) =
  template l:untyped = s.dest.l
  template si:untyped = s.si
  template sb:untyped = s.sb
  tFor i, s.subset.lowOuter..<s.subset.highOuter:
    let k1 = si.sq.pidx[i]
    #echo i, " ", k1
    if k1 >= 0:
      #echo "test6"
      assign(s.dest.s[i], s.src.s[k1])
    elif k1 + 2 <= 0:
      let k2 = -(k1 + 2)
      #echo i, " ", k2, " ", si.perm
      perm(s.dest.s[i], si.perm, s.src.s[k2])

proc boundary*(s:var Shift) =
  mixin blend
  template l:untyped = s.dest.l
  template si:untyped = s.si
  template sb:untyped = s.sb
  if si.nRecvDests > 0:
    if si.nRecvRanks > 0:
      if threadNum == 0:
        #echo "waitRecvBuf"
        waitRecvBuf(sb)
      twait0()
    if sb.sq.nthreads[threadNum] != numThreads:
      var ti0 = threadDivideLow(s.subset.lowOuter, s.subset.highOuter)
      var ti1 = threadDivideHigh(s.subset.lowOuter, s.subset.highOuter)
      var i0 = 0
      var step = (si.nRecvDests+1) div 2
      template search(i0,ti0:untyped):untyped =
        while true:
          while i0<si.nRecvDests and si.sq.recvDests[i0]<ti0:
            i0 += step
            step = (step+1) div 2
          if i0>si.nRecvDests: i0 = si.nRecvDests
          if i0==0 or si.sq.recvDests[i0-1]<ti0: break
          while i0>0 and si.sq.recvDests[i0-1]>=ti0:
            i0 -= step
            step = (step+1) div 2
          if i0<0: i0 = 0
      search(i0, ti0)
      var i1 = i0
      step = (i1+si.nRecvDests+1) div 2
      search(i1, ti1)
      sb.sq.offr[threadNum] = cint(i0)
      sb.sq.lenr[threadNum] = cint(i1)
      sb.sq.nthreads[threadNum] = numThreads.cint

    let ti0 = sb.sq.offr[threadNum]
    let ti1 = sb.sq.lenr[threadNum]
    #echo myrank, ": here1"
    if si.blend == 0:
      let rr = cast[ptr array[0,s.T]](sb.sq.rbuf)
      for i in ti0..<ti1:
        let k0 = si.sq.recvDests[i]
        let k2 = si.sq.recvRemoteSrcs[i]
        #echo "blend0: ", i, " ir: ", k0, " k2: ", k2
        #echo myrank, ": ", i, " ", k0, " ", k2
        #echo myrank, ": rr: ", cast[int](rr)
        #echo myrank, ": rbuf[0]: ", sb.sq.rbuf[0].ord
        #echo myrank, ": rr[k2][0][0]: ", rr[k2][0][0]
        s.dest.s[k0] = rr[k2]
    else:
      let stride = sizeof(s.T) div 2
      for i in ti0..<ti1:
        let k0 = si.sq.recvDests[i]
        #let k1 = si.sq.recvLocalSrcs[i]
        let k2 = si.sq.recvRemoteSrcs[i]
        #echo myrank, " ", i, " ", k0, " ", k1, " ", k2, " ", si.blend
        #blend(s.dest.s[k0], s.src.s[k1], sb.sq.rbuf[stride*k2].addr, si.blend)
        blend(s.dest.s[k0], sb.lbuf[stride*i].addr, sb.sq.rbuf[stride*k2].addr, si.blend)
    #echo myrank, ": here2"
    if si.nRecvRanks > 0:
      if threadNum == 0:
        doneRecvBuf(sb)
  if si.nSendRanks > 0:
    if threadNum == 0:
      waitSendBuf(sb)
  #QMP_barrier()
  threadBarrier()
  if threadNum == 0:
    freeShiftBuf(sb)

proc shift*(dest:var Field; dir,len:int; sub:string; src:Field) =
  const v = dest.V
  var s{.global.}:Shift[v,dest.T]
  #threadBarrier()
  s.init(dest, dir, len, sub)
  threadBarrier()  # wait for init
  s.start(src)
  threadBarrier()  # wait for src to be set
  s.local()
  #threadBarrier()
  s.boundary()
  threadBarrier()  # wait for everyone
proc shift*(dest:var Field; dir,len:int; src:Field) =
  shift(dest, dir, len, "all", src)

type
  Transporter*[U,F,T] = object
    link*: U
    field*: F
    sb*: ShiftB[T]
    len: int
  Shifter*[F,T] = Transporter[void,F,T]

#proc field*(t: Transporter): auto = t.field
#proc `link=`*(r,x: Transporter) = t.field

proc newTransporters*[U,F](u: seq[U], f: F, len: int, sub="all"): auto =
  var r: seq[Transporter[type(u[0]),F,type(f[0])]]
  let nd = u.len
  r.newSeq(nd)
  for mu in 0..<nd:
    r[mu].link = u[mu]
    r[mu].field = f.newOneOf
    #r[mu].field := 0
    r[mu].sb.initShiftB(f, mu, len, sub)
    r[mu].len = len
  r

proc newShifter*[F](f: F, dir,len: int, sub="all"): auto =
  var r: Shifter[F,type(f[0])]
  r.field = f.newOneOf
  #r.field := 0
  r.sb.initShiftB(f, dir, len, sub)
  r.len = len
  r

proc newShifters*[F](f: F, len: int, sub="all"): auto =
  var r: seq[Transporter[void,F,type(f[0])]]
  let nd = f.l.nDim
  r.newSeq(nd)
  for mu in 0..<nd:
    #r[mu].link = u[mu]
    r[mu].field = f.newOneOf
    #r[mu].field := 0
    r[mu].sb.initShiftB(f, mu, len, sub)
    r[mu].len = len
  r

proc `^*`*(x: Transporter, y: any): auto =
  mixin mul, load1, adj
  var r = x.field
  when compiles(x.link):
    if x.len >= 0:
      startSB(x.sb, y[ix])
      for ir in x.sb.subset:
        localSB(x.sb, ir, mul(r[ir], x.link[ir], it), load1(y[ix]))
      boundarySB(x.sb, mul(r[ir], x.link[ir], it))
    else:
      startSB(x.sb, x.link[ix].adj*y[ix])
      for ir in x.sb.subset:
        localSB(x.sb, ir, assign(r[ir], it), x.link[ix].adj*y[ix])
      boundarySB(x.sb, assign(r[ir], it))
  else:
    if x.len >= 0:
      startSB(x.sb, y[ix])
      for ir in x.sb.subset:
        localSB(x.sb, ir, assign(r[ir], it), load1(y[ix]))
      boundarySB(x.sb, assign(r[ir], it))
    else:
      startSB(x.sb, y[ix])
      for ir in x.sb.subset:
        localSB(x.sb, ir, assign(r[ir], it), y[ix])
      boundarySB(x.sb, assign(r[ir], it))
  threadBarrier()
  r
#template `()`*(x: Transporter, y: untyped): untyped = x ^* y

when isMainModule:
  import qex
  import physics/qcdTypes
  qexInit()
  echo "rank ", myRank, "/", nRanks
  echo threadNum, "/", numThreads, "/", numThreads
  #var lat = [4,4,4,4]
  var lat = [8,8,8,8]
  var lo = newLayout(lat)
  #lo.makeShift(0,1)
  #lo.makeShift(1,1)
  #lo.makeShift(2,1)
  #lo.makeShift(3,1)
  #layout.makeShift(3,-2,"even")
  var v1 = lo.ColorVector()
  var v2 = lo.ColorVector()
  var v3 = lo.ColorVector()
  let (sub1,sub2) = ("all","all")
  #let sub2 = "even"
  #let sub1 = "odd"
  proc lex(v,off: any) =
    for e in v.all:
      let lo = v.l
      let x = lo.vcoords(e)
      #v1.s[e] := e + 1
      var aa:array[lo.V,float32]
      for i in 0..<lo.V:
        template ff(n,r: untyped): untyped =
          r*10 + `mod`(x[n][i]+off[n]+lo.localGeom[n],lo.localGeom[n])
        #aa[i] = (((x[0][i]*10+x[1][i])*10+x[2][i])*10+x[3][i]).float32
        aa[i] = ff(3,ff(2,ff(1,ff(0,0)))).float32
      #var aa:array[lo.V,int]
      #aa = x[3]*10
      #let aa = ((x[3]*10+x[2])*10+x[1])*10+x[0]
      #assign(v1.s[e][0].re, aa)
      v[e] := 0
      v[e][0].re := aa
      if e==0:
        #echo aa
        echo v[e][0].re
        #echo v1[e][0].im
        #echo v1[e][0]
  threads:
    #echo v1[0].isVector
    v1 := -1
    threadBarrier()
    lex(v1, [0,0,0,0])
    #v1 := 1
    #v3 := v1
    threadBarrier()
    if threadNum==0:
      echo myRank, ": ", v1[0][0]
    for dir in 0..<lat.len:
      var disp = [0,0,0,0]
      v2 := -2
      shift(v2, dir,1, sub1, v1)
      echo v2[0][0].re
      disp[dir] = 1
      v3 := -3
      lex(v3, disp)
      echo "d2: ", (v2-v3).norm2
      for e in v1[sub1]:
        if (v2[e]-v3[e]).norm2.simdSum > 0:
          echo myrank, "\t", e, "\t", v2[e]
          echo myrank, "\t", e, "\t", v3[e]
      v3[sub2] := -1
      shift(v3, dir,-1, sub2, v2)
      threadMaster:
        echo v3.s[0][0]
      #for e in v1.all:
      #  echo myrank, "\t", e, "\t", v1.s[e][0]
      #  echo myrank, "\t", e, "\t", v3.s[e][0]
      let d2 = norm2(v1-v3)
      #if myRank==0:
      #  threadSingle:
      #    echo d2
      printf("d2: %g\n", d2)

  type st = type(v1[0])
  var
    sf*:array[4,ShiftB[st]]
    sb*:array[4,ShiftB[st]]
    ss = "even"
  for mu in 0..<4:
    initShiftB(sf[mu], lo, st, mu, 1, ss)
    initShiftB(sb[mu], lo, st, mu,-1, ss)
  threads:
    for mu in 0..<4:
      startSB(sf[mu], v1[ix])
      for ir in v2[ss]:
        localSB(sf[mu], ir, v2[ir] := it, v1[ix])
      boundarySB(sf[mu], v2[ir] := it)
      shift(v3, mu,1, ss, v1)
      onSubset ss:
        echo norm2(v2)
        echo norm2(v3)
        let d2 = norm2(v2-v3)
        echo d2

  #for e in v1.all:
  #  if (v1[e]-v3[e]).norm2.simdSum > 0:
  #    let x = lo.vcoords(e)
  #    echo x
  #    echo myrank, "\t", e, "\t", v1.s[e][0]
  #    echo myrank, "\t", e, "\t", v3.s[e][0]

  threads:
    shiftExpr(sf[0], v2[ir]:=it, v1[ix])

  import gauge
  var g = lo.newGauge
  let t = newTransporters(g, v1, 1)
  let td = newTransporters(g, v1, -1)
  threads:
    v2 := 0
    v2 += t[0] ^* v1
    v2 := t[0] ^* t[1] ^* t[2] ^* v1
    v2 := t[1] ^* t[0] ^* td[1] ^* v1

  qexFinalize()
