import base
import layout
import field
import io
#import qex
#import stdUtils
import times
#import profile
import os
import strUtils, sequtils
import maths, rng, physics/qcdTypes

import std/[hashes, tables]

#[
type
  GroupKind* = enum
    gkU, gkSU, gkHerm, gkAntiHerm  # traceless, real/complex
  Gauge*[T] = object
    u*: seq[T]
    n*: int
    group*: GroupKind
]#


proc newGauge*(l: Layout): auto =
  let nd = l.nDim
  result = newSeq[type(l.ColorMatrix())](nd)
  for i in 0..<nd:
    result[i] = l.ColorMatrix()
    result[i] := 1

proc newGaugeS*(l: Layout): auto =
  let nd = l.nDim
  result = newSeq[type(l.ColorMatrixS())](nd)
  for i in 0..<nd:
    result[i] = l.ColorMatrixS()
    result[i] := 1

proc newGaugeS*[T](g: seq[T]): auto =
  let nd = g.len
  result = newSeq[type(g[0].l.ColorMatrixS())](nd)
  for i in 0..<nd:
    result[i] = g[0].l.ColorMatrixS()
    result[i] := g[i]

proc newOneOf*[T](x: seq[T]): seq[T] =
  result.newSeq(x.len)
  for i in 0..<x.len:
    result[i] = newOneOf(x[i])

proc loadGauge*[T](g:openArray[T]; fn:string):int =
  var rd = g[0].l.newReader(fn)
  if rd.status!=0: return rd.status
  echo "loading gauge file: " & fn
  echo " file MD: " & rd.fileMetadata
  echo " rec MD: " & rd.recordMetadata
  echo " rec date: " & rd.recordDate
  #echo rd.datatype
  #echo rd.precision
  #echo rd.colors
  #echo rd.spins
  #echo rd.typesize
  #echo rd.datacount
  var t0 = epochTime()
  rd.read(g)
  t0 = epochTime() - t0
  if rd.status!=0: return rd.status
  echo " load seconds: " & $t0
  rd.close()
  if rd.status!=0: return rd.status
  return 0

const defFileMd = "<?xml version=\"1.0\"?>\n<note>generated by QEX</note>\n"
const defRecMd = "<?xml version=\"1.0\"?>\n<note>gauge configuration</note>\n"
proc saveGauge*[T](g:openArray[T]; fn:string; prec=""; filemd=defFileMd;
                   recordmd=defRecMd):int =
  var wr = g[0].l.newWriter(fn, filemd)
  if wr.status!=0: return wr.status
  var t0 = epochTime()
  wr.write(g, recordmd, prec)
  t0 = epochTime() - t0
  if wr.status!=0: return wr.status
  echo " save seconds: " & $t0
  wr.close
  if wr.status!=0: return wr.status
  return 0

proc setBC*(g: openArray[Field]) =
  let gt = g[3]
  tfor i, 0..<gt.l.nSites:
    #let e = i div gt.l.nSitesInner
    if gt.l.coords[3][i] == gt.l.physGeom[3]-1:
      gt{i} *= -1
      #echoAll isMatrix(gt{i})
      #echoAll i, " ", gt[e][0,0]

#[
template makeShifts(f:untyped):untyped {.dirty.} =
  template f(mu:int; b:untyped):untyped =
    tc = (tc+1) mod nf
    shiftFwd(tf[tc], mu, b)
    tf[tc]
  template `f "d"`(mu:int; b:untyped):untyped {.dirty.} =
    tc = (tc+1) mod nf
    shiftBck(tf[tc], mu, b)
    tf[tc]
template makeTransporters(f,g:untyped):untyped {.dirty.} =
  template f(mu:int; b:untyped):untyped =
    tc = (tc+1) mod nf
    transportFwd(tf[tc], g[mu], mu, b)
    tf[tc]
  template `f "d"`(mu:int; b:untyped):untyped {.dirty.} =
    tc = (tc+1) mod nf
    transportBck(tf[tc], g[mu], mu, b)
    tf[tc]
template setupTransporters(n,fld,nff:untyped):untyped =
  var sbf{.inject.} = createShiftBufs(n, fld, 1, "all")
  var sbb{.inject.} = createShiftBufs(n, fld, -1, "all")
  template SFT(r,mu,f:untyped):untyped =
    sbf[mu].startSB(f[ix])
    for ir in r.all:
      sbf[mu].localSB(ir, r[ir]:=it, f[ix])
    for ir in r.all:
      sbf[mu].boundarySB(ir, r[ir]:=it)
  proc shiftFwd(r,mu,f:auto) = SFT(r, mu, f); threadBarrier()
  #template TFT(r,g,mu,f:untyped):untyped =
  proc TFT(r,g,mu,f:auto) =
    sbf[mu].startSB(f[ix])
    for ir in r.all:
      sbf[mu].localSB(ir, r[ir]:=g[ir]*it, f[ix])
    sbf[mu].boundarySB(r[ir]:=g[ir]*it)
  proc transportFwd(r,g,mu,f:auto) = TFT(r, g, mu, f); threadBarrier()
  #template SBT(r,mu,f:untyped):untyped =
  proc SBT(r,mu,f:auto) =
    sbb[mu].startSB(f[ix])
    for ir in r.all:
      sbb[mu].localSB(ir, r[ir]:=it, f[ix])
    sbb[mu].boundarySB(r[ir]:=it)
  proc shiftBck(r,mu,f:auto) = SBT(r, mu, f); threadBarrier()
  #template TBT(r,g,mu,f:untyped):untyped =
  proc TBT(r,g,mu,f:auto) =
    sbb[mu].startSB(g[ix].adj*f[ix])
    for ir in r.all:
      sbb[mu].localSB(ir, r[ir]:=it, g[ix].adj*f[ix])
    sbb[mu].boundarySB(r[ir]:=it)
  proc transportBck(r,g,mu,f:auto) = TBT(r, g, mu, f); threadBarrier()
  var nf{.inject.} = nff
  var tf{.inject.} = newSeq[type(fld)](nf)
  for i in 0..<nf: tf[i] = newOneOf(fld)
  #var tc{.inject.} = 0

# a * pathProduct( u[nu] ^* v[mu] ^* u[nu].dag ) + b*
# s[mu] = a_mu s[mu] + f_mu_nu Unu Vmu Unu^+ + b_mu_nu Unu^+ Vmu Unu
proc staples*[T,F,B](ss,uu,vv:openArray[T];ff:openArray[F];bb:openArray[B]) =
  let nd = ss.len
  let s = cast[ptr cArray[T]](unsafeAddr(ss[0]))
  let u = cast[ptr cArray[T]](unsafeAddr(uu[0]))
  let v = cast[ptr cArray[T]](unsafeAddr(vv[0]))
  let f = cast[ptr cArray[F]](unsafeAddr(ff[0]))
  let b = cast[ptr cArray[B]](unsafeAddr(bb[0]))
  setupTransporters(nd, u[0], 10)
  makeShifts(S)
  makeTransporters(U, u)
  makeTransporters(V, v)

  threads:
    var tc = 0
    for mu in 0..<nd:
      for nu in 0..<nd:
        if nu==mu: continue
        s[mu] += f[mu][nu] * U(nu, V(mu, Sd(nu, u[nu]))) +
                 b[mu][nu] * Ud(nu, V(mu, u[nu]))
        #s[mu] += f[mu][nu] * ( U(nu, v[mu]) * Sd(mu, u[nu])) +
        #         b[mu][nu] * Ud(nu, V(mu, u[nu]))
]#

proc plaq*[T](uu: openArray[T]): auto =
  mixin mul, load1, createShiftBufs
  tic()
  template getIp(mu,nu: int): int = ((mu*(mu-1)) div 2) + nu
  let u = cast[ptr cArray[T]](unsafeAddr(uu[0]))
  let lo = u[0].l
  let nd = lo.nDim
  let nc = u[0][0].ncols
  var sf = newSeq[type(createShiftBufs(u[0],1,"all"))](nd)
  for i in 0..<nd-1:
    sf[i] = createShiftBufs(u[0], 1, "all")
  sf[nd-1].newSeq(nd)
  for i in 0..<nd-1: sf[nd-1][i] = sf[i][i]
  let np = (nd*(nd-1)) div 2
  var pl = newSeq[float64](np)
  toc("plaq setup")
  threads:
    tic()
    #var plt = newSeq[float64](np)
    var plt: array[6,float64]
    var umunu,unumu: type(load1(u[0][0]))
    for mu in 0..<nd:
      for nu in 0..<nd:
        if mu != nu:
          startSB(sf[mu][nu], u[mu][ix])
    toc("plaq start shifts")
    for ir in u[0]:
      for mu in 1..<nd:
        for nu in 0..<mu:
          if isLocal(sf[mu][nu],ir) and isLocal(sf[nu][mu],ir):
            localSB(sf[mu][nu], ir, mul(unumu,u[nu][ir],it), u[mu][ix])
            localSB(sf[nu][mu], ir, mul(umunu,u[mu][ir],it), u[nu][ix])
            let ip = getIp(mu,nu)
            let dt = redot(umunu,unumu)
            plt[ip] += simdSum(dt)
    toc("plaq local")
    var needBoundary = false
    for mu in 0..<nd:
      for nu in 0..<nd:
        if mu != nu:
          boundaryWaitSB(sf[mu][nu]): needBoundary = true
    toc("plaq wait")
    if needBoundary:
      boundarySyncSB()
      for ir in u[0]:
        for mu in 1..<nd:
          for nu in 0..<mu:
            if not isLocal(sf[mu][nu],ir) or not isLocal(sf[nu][mu],ir):
              if isLocal(sf[mu][nu], ir):
                localSB(sf[mu][nu], ir, mul(unumu,u[nu][ir],it), u[mu][ix])
              else:
                boundaryGetSB(sf[mu][nu], ir):
                  mul(unumu, u[nu][ir], it)
              if isLocal(sf[nu][mu], ir):
                localSB(sf[nu][mu], ir, mul(umunu,u[mu][ir],it), u[nu][ix])
              else:
                boundaryGetSB(sf[nu][mu], ir):
                  mul(umunu, u[mu][ir], it)
              let ip = getIp(mu,nu)
              let dt = redot(umunu,unumu)
              plt[ip] += simdSum(dt)
    toc("plaq boundary")
    threadSum(plt)
    if threadNum == 0:
      for i in 0..<pl.len:
        pl[i] = plt[i]/(lo.physVol.float*float(np*nc))
      rankSum(pl)
    toc("plaq sum")
  result = pl
  toc("plaq end", flops=lo.nSites.float*float(np*(2*8*nc*nc*nc-1)))

discard """
# s[mu] = a_mu s[mu] + f_mu_nu Unu Vmu Unu^+ + b_mu_nu Unu^+ Vmu Unu
proc staples*[T,A,F,B](staples,uu,vv:openArray[T]; aa:openArray[A];
                       ff:openArray[F]; bb:openArray[B]):auto =
  let s = cast[ptr cArray[T]](unsafeAddr(ss[0]))
  let u = cast[ptr cArray[T]](unsafeAddr(uu[0]))
  let v = cast[ptr cArray[T]](unsafeAddr(vv[0]))
  let a = cast[ptr cArray[A]](unsafeAddr(aa[0]))
  let f = cast[ptr cArray[F]](unsafeAddr(ff[0]))
  let b = cast[ptr cArray[B]](unsafeAddr(bb[0]))
  let n = staples.len
  # sbv[mu,nu].startSB(v[nu][ix])
  # sbu[mu,nu].startSB(u[nu][ix])
  # sbv[mu,nu].localSB(ir,vt[mu][nu][ir]:=u[mu][ir]*it,v[nu][ix])
  # sbu[mu,nu].localSB(ir,ut[mu][nu][ir]:=v[mu][ir]*it,u[nu][ix])
  # sbv[mu,nu].boundarySB(ir,vt[mu][nu][ir]:=u[mu][ir]*it)
  # sbu[mu,nu].boundarySB(ir,ut[mu][nu][ir]:=v[mu][ir]*it)

  # s[mu] := a[mu] * s[mu]

  # sbv[nu,mu].startSB(ut[mu][nu][ix])
  # sbv[nu,mu].localSB(ir,s[mu][ir]+=b[mu][nu]*(u[nu][ir].adj*it),ut[mu][nu][ix])
  # sbv[nu,mu].boundarySB(ir,s[mu][ir]+=b[mu][nu]*(u[nu][ir].adj*it))

  # sbu[mu,nu].localSB(ir, s[mu][ir] += f[mu][nu]*(vt[nu][mu]*it.adj), u[nu][ix])
  # sbu[mu,nu].boundarySB(ir, s[mu][ir] += f[mu][nu]*(vt[nu][mu]*it.adj))

  #threads:
  #  for mu in 0..<n:
"""

proc plaq2*[T](gg:openArray[T]):auto =
  mixin adj
  tic()
  let g = cast[ptr cArray[T]](unsafeAddr(gg[0]))
  let lo = g[0].l
  let nd = lo.nDim
  let nc = g[0][0].ncols
  var m = lo.ColorMatrix()
  var s0 = lo.ColorMatrix()
  #var t0 = lo.ColorMatrix()
  var s1 = lo.ColorMatrix()
  #var t1 = lo.ColorMatrix()
  var tr:type(trace(m))
  toc("plaq2 setup")
  threads:
    #tic()
    m := 0
    #toc("plaq2 zero")
    for mu in 1..<nd:
      for nu in 0..<mu:
        #tic()
        shift(s0, mu,1, g[nu])
        #toc("plaq2 shift1")
        shift(s1, nu,1, g[mu])
        #toc("plaq2 shift2")
        #echo "s0: ", trace(s0)
        #echo "s1: ", trace(s1)
        m += (g[mu]*s0) * (g[nu]*s1).adj
        #m += (g[mu]*s0) * (g[nu]*s1)
        #echo mu, " ", nu, " ", trace(m)/nc
        #toc("plaq2 mul")
    #toc("plaq2 work")
    tr = trace(m)
    #toc("plaq2 trace")
  toc("plaq2 threads")
  result = tr/(lo.physVol.float*0.5*float(nd*(nd-1)*nc))

proc plaq3*[T](g: seq[T]): auto =
  mixin adj, newTransporters
  tic()
  let lo = g[0].l
  let nd = lo.nDim
  let nc = g[0][0].ncols
  let t = newTransporters(g, g[0], 1)
  var m = lo.ColorMatrix()
  var tr: type(trace(m))
  toc("plaq3 setup")
  threads:
    tic()
    m := 0
    toc("plaq3 zero")
    for mu in 1..<nd:
      for nu in 0..<mu:
        tic()
        m += (t[mu]^*g[nu]) * (t[nu]^*g[mu]).adj
        #echo mu, " ", nu, " ", trace(m)/nc
        toc("plaq3 mul")
    toc("plaq3 work")
    tr = trace(m)
    toc("plaq3 trace")
  toc("plaq3 threads")
  result = tr/(lo.physVol.float*0.5*float(nd*(nd-1)*nc))

proc echoPlaq*(g: auto) =
  var pl = plaq(g)
  let nt = g.len - 1
  let ns = pl.len - nt
  var ps, pt: float
  for i in 0..<ns: ps += pl[i]
  for i in ns..<pl.len: pt += pl[i]
  echo "plaqS: ", 2*ps, "  plaqT: ", 2*pt, "  plaq: ", ps+pt

type
  Link[F:ref] = object
    field:F
    forward:bool
    slen:seq[int]
proc wlineReduce(line:openarray[Link]):seq[Link] =
  ## Compute any (possibly gapped) Wilson line in pairs of links,
  ## and return a reduced line.
  ## `line` is a list of fields, directions and shift lengths,
  ## where `forward` denotes whether taking an adjoint of `field`,
  ## `slen[i]` is the shift length in `i` direction,
  ## and `abs(line[i])>1` means gaps in the line.
  tic()
  type
    F = Link.F
    S = typeof(line[0].field[0])
    OPs = tuple
      u1:F
      f1:bool
      u2:F
      f2:bool
      slen:seq[int]
      res:F
  let
    lo = line[0].field.l
    nd = lo.nDim
    nl = (line.len+1) div 2
  var
    lnew = newseq[Link]()
    fs = newseq[OPs]()
    sf = newseq[Shifter[F,S]](nd)
    sb = newseq[Shifter[F,S]](nd)
    pd = newseq[int](nd)
    pf = true
  for i in 0..<nl:
    let
      a = i+i
      b = a+1
    if b < line.len:
      let
        (u1, f1, p) = block:
          let la = line[a]
          (la.field, la.forward, la.slen)
        (u2, f2, s) = block:
          let lb = line[b]
          (lb.field, lb.forward, lb.slen)
      var
        res:F = nil
        forward = true
      for k in 0..<fs.len:
        if f1 == fs[k].f1 and f2 == fs[k].f2 and
            u1 == fs[k].u1 and u2 == fs[k].u2 and s == fs[k].slen:
          res = fs[k].res
          break
        elif f1 != fs[k].f2 and f2 != fs[k].f1 and
            u1 == fs[k].u2 and u2 == fs[k].u1 and s == fs[k].slen:
          res = fs[k].res
          forward = false
          break
      if res.isnil:
        res = newOneOf(u1)
        if f1 or f2:
          fs.add (u1, f1, u2, f2, s, res)
        else:
          let ms = s.mapIt(-it)
          fs.add (u2, not f2, u1, not f1, ms, res)
          forward = false
      # echo s, p
      var dis = p
      if pf:
        for k in 0..<p.len:
          dis[k] += pd[k]
      if forward:
        pd = s
      else:
        for k in 0..<p.len:
          dis[k] += s[k]
        pd = s.mapIt(-it)
      pf = forward
      # echo forward," ",dis
      lnew.add Link(field:res, forward:forward, slen:dis)
    else:
      lnew.add line[a]
  toc("wlineReduce prepare")

  # Only use single shift in the following until we resolve the issue with multi-shifts.
  block:
    var f,b = newseq[bool](nd)
    for i in 0..<fs.len:
      let sl = fs[i].slen
      for d in 0..<sl.len:
        if (not f[d]) and sl[d] > 0:
          f[d] = true
          sf[d] = newShifter(line[0].field, d, 1)
        elif (not b[d]) and sl[d] < 0:
          b[d] = true
          sb[d] = newShifter(line[0].field, d, -1)
  toc("wlineReduce init shifts")

  threads:
    tic()
    for (u1, f1, u2, f2, slen, res) in fs:
      res := u2
      # echo f1," ",f2," ",slen
      for mu in 0..<nd:
        let n = slen[mu]
        # echo n
        for i in 1..abs(n):
          if n > 0:
            res := sf[mu] ^* res
          elif n < 0:
            res := sb[mu] ^* res
      toc("wlineReduce Thr shifts")
      if f1 and f2:
        res := u1 * res
      elif not (f1 or f2):
        qexError("Internal logic error.")
      elif f1:
        res := u1 * res.adj
      else:
        res := u1.adj * res
      toc("wlineReduce Thr mul")
    toc("wlineReduce Thr")
  toc("wlineReduce compute")
  lnew

proc wline0*(g:auto, line:openarray[int]):auto =
  ## Compute the trace of ordered product of gauge links, the Wilson Line.
  ## The line is given as a list of integers +/- 1..nd, where the sign
  ## denotes forward/backward and the number denotes the dimension.
  tic()
  # echo line
  type L = Link[typeof(g[0])]
  const nc = g[0][0].ncols
  let
    vol = g[0].l.physVol
    nd = g[0].l.nDim
  var
    link = newseq[L]()
    pd = -1
    pf = true
  for d in line:
    let
      mu = abs(d) - 1
      forward = d > 0
    if mu<0 or mu>=nd:
      qexError("wline: Unable to parse line: ", line)
    var s = newseq[int](nd)
    if pd>=0 and pf:
      s[pd] += 1
    if not forward:
      s[mu] -= 1
    link.add L(field:g[mu], forward:forward, slen:s)
    pd = mu
    pf = forward
  while link.len > 1:
    # echo "Call wlineReduce: ",link.mapIt((it.forward,it.slen))
    link = link.wlineReduce
  toc("wline wlineReduce")
  var r: typeof(trace(g[0]))
  if link.len == 1:
    threads:
      r = link[0].field.trace
  result = r / float(vol*nc)
  toc("wline trace")

type
  Coord = object
    x: seq[int]

proc `$`(x:Coord):string =
  return "Coord" & $x.x

proc `+`(x,y:Coord):Coord =
  let
    nx = x.x.len
    ny = y.x.len
    n = max(nx,ny)
    m = min(nx,ny)
  result.x.newseq n
  for i in 0..<m:
    result.x[i] = x.x[i] + y.x[i]
  if nx>ny:
    for i in m..<n:
      result.x[i] = x.x[i]
  else:
    for i in m..<n:
      result.x[i] = y.x[i]

proc `-`(x,y:Coord):Coord =
  let
    nx = x.x.len
    ny = y.x.len
    n = max(nx,ny)
    m = min(nx,ny)
  result.x.newseq n
  for i in 0..<m:
    result.x[i] = x.x[i] - y.x[i]
  if nx>ny:
    for i in m..<n:
      result.x[i] = x.x[i]
  else:
    for i in m..<n:
      result.x[i] = -y.x[i]

proc `-`(x:Coord):Coord =
  let
    n = x.x.len
  result.x.newseq n
  for i in 0..<n:
    result.x[i] = -x.x[i]

proc `+=`(x:var Coord,y:Coord) =
  let
    nx = x.x.len
    ny = y.x.len
  if nx<ny:
    x.x.setLen ny
  for i in 0..<ny:
    x.x[i] += y.x[i]

proc `-=`(x:var Coord,y:Coord) =
  let
    nx = x.x.len
    ny = y.x.len
  if nx<ny:
    x.x.setLen ny
  for i in 0..<ny:
    x.x[i] -= y.x[i]

proc `[]`(x:Coord, i:SomeInteger):int =
  if i>=x.x.len:
    result = 0
  else:
    result = x.x[i]

proc `[]`(x:var Coord, i:SomeInteger):var int =
  if i>=x.x.len:
    x.x.setLen(i+1)
  return x.x[i]

proc `[]=`(x:var Coord, i:SomeInteger, d:int) =
  if i>=x.x.len:
    x.x.setLen(i+1)
  x.x[i] = d

type
  OrdPathKind = enum
    opInt, opPair, opList, opAdj
  OrdPath = ref object
    case k: OrdPathKind
    of opInt: d:int
    of opPair: l,r:OrdPath
    of opList: s:seq[OrdPath]
    of opAdj: p:OrdPath
  OrdPathTree* = ref object
    paths: seq[OrdPath]
    segments: seq[OrdPath]
    counts: seq[int]

proc `$`(path:OrdPath):string =
  return
    case path.k:
    of opInt: $path.d
    of opPair: "(" & $path.l & " " & $path.r & ")"
    of opList: $path.s
    of opAdj: "Adj" & $path.p

proc newOrdPath(path:openarray[int]):OrdPath =
  let n = path.len
  var op = OrdPath(k:opList, s:newseq[OrdPath](n))
  for i in 0..<n:
    let d = path[i]
    if d==0:
      qexError("newOrdPath: unable to parse path: ", path)
    op.s[i] = OrdPath(k:opInt, d:d)
  return op

proc adjoint(p:OrdPath):OrdPath =
  case p.k:
  of opInt:
    return OrdPath(k:opInt, d: -p.d)
  of opPair:
    return OrdPath(k:opPair, l:p.r.adjoint, r:p.l.adjoint)
  of opList:
    var
      l = p.s.len
      l1 = l-1
      n = newseq[OrdPath](l)
    for i in 0..<l:
      n[i] = p.s[l1-i].adjoint
    return OrdPath(k:opList, s:n)
  of opAdj:
    return p.p

proc adjointOp(p:OrdPath):OrdPath =
  ## Like adjoint(), but create new OrdPath(k:opAdj) if it's not already an opInt or opAdj.
  case p.k:
  of opInt:
    return p.adjoint
  of opAdj:
    return p.p
  else:
    return OrdPath(k:opAdj, p:p)

proc flatten(path:OrdPath):seq[int] =
  case path.k:
  of opInt:
    return @[path.d]
  of opPair:
    var r = path.l.flatten
    r.add path.r.flatten
    return r
  of opList:
    var r = newseq[int]()
    for l in path.s:
      r.add l.flatten
    return r
  of opAdj:
    return path.p.adjoint.flatten

proc deltaX(path:OrdPath):Coord =
  # return the Coord difference between the end of the path to the beginning.
  case path.k:
  of opInt:
    let d = abs(path.d)-1
    if path.d>0:
      result[d] = 1
    else:
      result[d] = -1
  of opPair:
    result = path.l.deltaX + path.r.deltaX
  of opList:
    for l in path.s:
      result += l.deltaX
  of opAdj:
    result = -path.p.deltaX

proc position(path:OrdPath):Coord =
  # return the Coord relative to the left most starting point.
  case path.k:
  of opInt:
    if path.d<0:
      result[-1-path.d] = -1
  of opPair:
    result = path.l.position
  of opList:
    if path.s.len>0:
      result = path.s[0].position
  of opAdj:
    result = path.p.position - path.p.deltaX

proc `$`*(path:OrdPathTree):string =
  result = "OrdPathTree:\n  paths:"
  for p in path.paths:
    result.add "\n    " & $p.flatten & " \t " & $p & " \t shift " & $p.position.x
  result.add "\n  segments:"
  for i in 0..<path.segments.len:
    let
      l = path.segments[i].l
      r = path.segments[i].r
    result.add "\n    " & $path.counts[i] & "  " & $l & " \t " & $r & " \t shift " & $(l.position - l.deltaX - r.position).x

proc hash(path:OrdPath):Hash =
  return
    case path.k:
    of opInt: !$(opInt.hash !& path.d.hash)
    of opPair: !$(opPair.hash !& path.l.hash !& path.r.hash)
    of opList: !$(opList.hash !& path.s.hash)
    of opAdj: !$(opAdj.hash !& path.p.hash)

proc `==`(x,y:OrdPath):bool {.noSideEffect.} =  # noSideEffect needed for Nim devel
  # strictly follows hash(path:OrdPath)
  return
    if x.k != y.k:
      false
    else:
      case x.k:
      of opInt: x.d == y.d
      of opPair: x.l == y.l and x.r == y.r
      of opList: x.s == y.s
      of opAdj: x.p == y.p

proc mostSharedPair(paths:openarray[OrdPath]):(OrdPath,int) =
  ## Receives paths and search each element of OrdPath(k:opList).
  ## Return an OrdPath(k:opPair) occured most frequently among all the paths,
  ## and its count.
  ## If there are no pairs, return OrdPath(k:opList) with len 0, and 0.
  var pc = initTable[OrdPath,int]()
  var p,pa: OrdPath
  var fp,fpa: seq[int]
  var c,ca: int
  while true:
    c = 0
    ca = 0
    for ps in paths:
      if not (ps.k==opList and ps.s.len>1):
        continue
      var i = 0
      while i<ps.s.len-1:  # FIXME test reversed
        i.inc
        let t = OrdPath(k:opPair, l:ps.s[i-1], r:ps.s[i])
        let ft = t.flatten
        if c==0:
          p = t
          pa = OrdPath(k:opPair, l:t.l.adjointOp, r:t.r.adjointOp)
          if p in pc or pa in pc:
            continue
          fp = p.flatten()
          fpa.newseq(fp.len)
          for i in 0..<fp.len:
            fpa[i] = -fp[fp.len-i-1]
          c = 1
          i.inc
        elif ft == fp:
          c.inc
          i.inc
        elif ft == fpa:
          ca.inc
          i.inc
    if c==0:
      # no new pairs
      break
    else:
      let ct = c+ca
      if c>=ca:
        pc[p] = ct
      else:
        pc[pa] = ct
  c = 0
  for k,v in pc.pairs:
    # If there are multiple paris with the max count,
    # which one we return depends on implementation of Table.
    if v>c:
      c = v
      p = k
  if c==0:
    return (OrdPath(k:opList, s: @[]), 0)
  else:
    return (p, c)

proc groupByPair(paths:openarray[OrdPath], pair:OrdPath):seq[OrdPath] =
  ## Receives paths and group each element of OrdPath(k:opList) by the pair OrdPath(k:opPair).
  let
    fp = pair.flatten
    pairAdj = OrdPath(k:opAdj, p:pair)
  var afp = newseq[int](fp.len)
  for i in 0..<fp.len:
    afp[i] = -fp[fp.len-i-1]
  var newpaths = newseq[OrdPath]()
  for ps in paths:
    if ps.k!=opList:
      newpaths.add ps
      continue
    var
      p = newseq[OrdPath]()
      i = 0
      n = ps.s.len
    while i<n-1:
      var t = ps.s[i].flatten
      t.add ps.s[i+1].flatten
      if t == fp:
        p.add pair
        # echo "groupByPair: add ",p[^1]
        i.inc 2
      elif t==afp:
        p.add pairAdj
        # echo "groupByPair: add ",p[^1]
        i.inc 2
      else:
        p.add ps.s[i]
        i.inc
    if i==n-1:
      p.add ps.s[i]
    newpaths.add(
      if p.len==1:
        p[0]
      else:
        OrdPath(k:opList, s:p)
    )
  return newpaths

proc optimalPairs(paths:openarray[OrdPath]):OrdPathTree =
  ## Runtime O(N^3) with N=paths.len
  tic()
  var
    pgroup = @paths
    segments = newseq[OrdPath]()
    counts = newseq[int]()
    (p,c) = pgroup.mostSharedPair
  while p.k==opPair:
    # echo "optimalPairs: add ",p
    segments.add p
    counts.add c
    pgroup = pgroup.groupByPair p
    (p,c) = pgroup.mostSharedPair
  toc("optimalPairs")
  return OrdPathTree(paths:pgroup, segments:segments, counts:counts)

proc optimalPairs*(paths:openarray[seq[int]]):OrdPathTree =
  ## Runtime O(N^3) with N=paths.len
  let n = paths.len
  var ps = newseq[OrdPath](n)
  for i,p in paths.pairs:
    ps[i] = p.newOrdPath
  ps.optimalPairs

proc singleshift[F,S](f:F, sh:openarray[int], sf,sb:openarray[Shifter[F,S]]):F =
  # Result will be in one of the shift buffers.
  result = f
  for mu,n in sh.pairs:
    if n>0:
      result = sb[mu] ^* result
    elif n<0:
      result = sf[mu] ^* result

proc multishifts[F,S](f:F, sh:openarray[int], sf,sb:openarray[Shifter[F,S]]):F =
  # f will be used as temporary and overwritten.
  # Result will be in one of the shift buffers.
  result = f
  for mu,n in sh.pairs:
    if n>0:
      result = sb[mu] ^* result
      for i in 1..<n:
        f := result
        result = sb[mu] ^* f
    elif n<0:
      result = sf[mu] ^* result
      for i in 1..<(-n):
        f := result
        result = sf[mu] ^* f

proc lrmul(res:auto, l:auto, r:auto, la,ra:bool) =
  tic("lrmul")
  if la and ra:
    res := l.adj * r.adj
    toc("la ra")
  elif la:
    res := l.adj * r
    toc("la")
  elif ra:
    res := l * r.adj
    toc("ra")
  else:
    res := l * r
    toc("direct")

proc gaugeProd*(g:auto, ptree:OrdPathTree, origin=true):auto =
  tic("gaugeProd")
  type
    F = typeof(g[0])
    S = typeof(g[0][0])
  let nd = g[0].l.nDim
  var
    sf = newseq[Shifter[F,S]](nd)
    sb = newseq[Shifter[F,S]](nd)
    sl = newseq[seq[int]](ptree.segments.len)
    gp = initTable[seq[int],F]()
    sfi = newseq[bool](nd)
    sbi = newseq[bool](nd)
  for i,s in ptree.segments.pairs:
    if s.k==opPair:
      let sh = (s.l.position - s.l.deltaX - s.r.position).x
      sl[i] = sh
      for mu,n in sh.pairs:
        if n>0:
          sbi[mu] = true
        elif n<0:
          sfi[mu] = true
    else:
      qexError("Internal logic error.")
  for i in 0..<nd:
    if sfi[i]:
      sf[i] = newShifter(g[0], i, 1)
    if sbi[i]:
      sb[i] = newShifter(g[0], i, -1)
  toc("init shifts")

  proc fetch(s:OrdPath):auto =
    case s.k:
    of opInt:
      if s.d>0:
        return (g[s.d-1],false)
      else:
        return (g[-s.d-1],true)
    of opPair:
      return (gp[s.flatten],false)
    of opList:
      qexError("Internal logic error.")
    of opAdj:
      let (f,a) = fetch s.p
      return (f, a xor true)

  for i,s in ptree.segments.pairs:
    if s.k==opPair:
      let
        (l,la) = fetch s.l
        (r,ra) = fetch s.r
        sh = sl[i]
      # echo "gaugeProd:\n","  l: ",la," ",s.l.flatten,"\n    ",s.l.position,"  ",s.l.deltaX,"\n  r: ",ra," ",s.r.flatten,"\n    ",s.r.position,"  ",s.r.deltaX,"\n  sh: ",sh
      var needs = 0
      for x in sh:
        let ax = abs x
        if needs < ax:
          needs = ax
      var res = newOneOf(r)
      threads:
        let rr =
          if needs==0:
            r
          elif needs==1:
            r.singleshift(sh, sf, sb)
          else:
            res := r
            res.multishifts(sh, sf, sb)
        res.lrmul(l,rr,la,ra)
      # echo "gaugeProd: add ",s.flatten," ",res.trace/float(g[0].l.physVol*g[0][0].ncols)
      gp[s.flatten] = res
    else:
      qexError("Internal logic error.")
  toc("gaugeProd prod")
  let n = ptree.paths.len
  var res = newseq[F](n)
  var resAlloc = newseq[bool](n)  # TODO: implement tracing ref counting
  for i,p in ptree.paths.pairs:
    if p.k==opAdj:
      let t = gp[p.p.flatten]
      # echo "gaugeProd result ",i," adj path ",p.flatten," ",$p
      res[i] = newOneOf t
      resAlloc[i] = true
      threads:
        res[i] := t.adj
    else:
      res[i] = gp[p.flatten]
  if origin:
    var
      s = newseq[seq[int]](n)
      fi = newseq[bool](nd)
      bi = newseq[bool](nd)
    for i,p in ptree.paths.pairs:
      s[i] = (-p.position).x
      for mu,l in s[i].pairs:
        if l>0:
          bi[mu] = true
        elif l<0:
          fi[mu] = true
    for i in 0..<nd:
      if fi[i] and not sfi[i]:
        sf[i] = newShifter(g[0], i, 1)
      if bi[i] and not sbi[i]:
        sb[i] = newShifter(g[0], i, -1)
    for i,p in ptree.paths.pairs:
      # echo "gaugeProd result ",i," shifts ",s[i]," path ",p.flatten," ",$p
      var needs = 0
      for x in s[i]:
        let ax = abs x
        if needs < ax:
          needs = ax
      if needs>0:
        let r = res[i]
        if not resAlloc[i]:
          res[i] = newOneOf r
          resAlloc[i] = true
        threads:
          let t =
            if needs==1:
              r.singleshift(s[i], sf, sb)
            else:
              res[i] := r
              res[i].multishifts(s[i], sf, sb)
          res[i] := t
  result = res
  toc("done")

proc wilsonLines*(g:auto, lines:openarray[seq[int]]):auto =
  ## Compute the trace of ordered product of gauge links, the Wilson Line.
  ## Each line is given as a list of integers +/- 1..nd, where the sign
  ## denotes forward/backward and the number denotes the dimension.
  ## Note that the simplification has a runtime complexity of O(N^3) for N=lines.len.
  tic()
  let n = lines.len
  let rs = g.gaugeProd(lines.optimalPairs, false)
  toc("wilsonLines prod")
  type R = typeof(trace(g[0]))
  var
    r: R
    ts = newseq[R](n)
  const nc = g[0][0].ncols
  let
    vol = g[0].l.physVol
    fac = 1.0/float(vol*nc)
  for i,s in rs.pairs:
    threads:
      r = s.trace
    ts[i] := r*fac
  toc("wilsonLines trace")
  return ts

proc wline*(g:auto, line:openarray[int]):auto =
  return g.wilsonLines([@line])[0]

proc allCorners(path:openarray[int]):seq[seq[int]] =
  let np = path.len
  result = newseq[seq[int]]()
  var old = 0
  for i,x in path.pairs:
    if x!=old:
      var p = newseq[int](np)
      for j in 0..<np:
        p[j] = path[(j+i) mod np]
      result.add p
      old = x
  # echo "allCorners: ",result

proc fmunuCoeffs_fun(loop:int):auto =
  var k:array[5,float]
  case loop
  of 1:
    return [1.0,0,0,0,0]
  of 3:
    k[4] = 1.0/90.0
  of 4:
    k[4] = 0.0
  of 5:
    k[4] = 1.0/180.0
  else:
    discard
  k[0] = 19.0/9.0 - 55.0*k[4]
  k[1] = 1.0/36.0 - 16.0*k[4]
  k[2] = 64.0*k[4] - 32.0/45.0
  k[3] = 1.0/15.0 - 6.0*k[4]
  return k

proc fmunuCoeffs(loop:int):auto =
  const FmunuCoeffs =
    [ fmunuCoeffs_fun(1)
    , fmunuCoeffs_fun(3)
    , fmunuCoeffs_fun(4)
    , fmunuCoeffs_fun(5)
    ]
  case loop
  of 1:
    result = FmunuCoeffs[0]
  of 3:
    result = FmunuCoeffs[1]
  of 4:
    result = FmunuCoeffs[2]
  of 5:
    {.linearScanEnd.}
    result = FmunuCoeffs[3]
  else:
    qexError("fmunuCoeffs uses loop in [1,3,4,5], but got ",loop)

proc fmunuPTree_fun(mu,nu,loop:int):auto =
  tic()
  var lp = newseq[seq[int]]()
  block:
    let
      mu = mu+1
      nu = nu+1
    lp.add allCorners @[-mu,-nu,mu,nu]  # 1x1
    if loop>=3:
      lp.add allCorners @[-mu,-mu,-nu,-nu,mu,mu,nu,nu]  # 2x2
    if loop>=4:
      lp.add allCorners @[-mu,-mu,-nu,mu,mu,nu]  # 2x1
      lp.add allCorners @[-mu,-nu,-nu,mu,nu,nu]  # 1x2
      lp.add allCorners @[-mu,-mu,-mu,-nu,mu,mu,mu,nu]  # 3x1
      lp.add allCorners @[-mu,-nu,-nu,-nu,mu,nu,nu,nu]  # 1x3
    if loop==3 or loop==5:
      lp.add allCorners @[-mu,-mu,-mu,-nu,-nu,-nu,mu,mu,mu,nu,nu,nu]  # 3x3
  result = lp.optimalPairs
  toc("compute PTree")

proc fmunuPTree(mu,nu,loop:int):auto =
  const lps = [-1,0,-1,1,2,3]
  let lp = lps[loop]
  memoize(lp,mu,nu):
    fmunuPTree_fun(mu,nu,loop)

proc fmunu*(g:auto, mu,nu:int, loop=1):auto =
  ## mu,nu: 0..<nd
  ## loop: 1,3,4,5
  ## returns the traceless antihermitian F_munu from the clover
  tic("fmunu")
  if loop notin [1,3,4,5]:
    qexError("fmunu uses loop in [1,3,4,5], but got ",loop)
  const lpc = [4,4,8,8,4]  # counts of distinct loops in order of 1x1, 2x2, 1x2, 1x3, 3x3
  let cs = fmunuCoeffs loop
  let ptree = fmunuPTree(mu,nu,loop)
  toc("fmunuPTree")
  let fs = g.gaugeProd ptree
  toc("gaugeProd")
  let f = newOneOf(fs[0])
  threads:
    f := 0
    var i0 = 0
    for jc in 0..<loop:
      var j = jc
      if loop==3 and jc==2:
        j = 4
      let n = lpc[j]
      let ni = cs[j]/n.float
      for ir in f:
        var t:type(load1(f[0]))
        for i in 0..<n:
          t += fs[i0+i][ir]
        f[ir] += ni*t
      i0 += n
    f.projectTAH
  toc("single fmunu")
  return f

proc fmunu*(g:auto, loop=1):auto =
  ## loop: 1,3,4,5
  ## returns the traceless antihermitian F[mu][nu] where mu>nu
  tic("fmunu")
  type F = type(g[0])
  let nd = g[0].l.nDim
  var f = newseq[seq[F]](nd)  # only partially initialized
  for mu in 1..<nd:
    f[mu].newseq(mu)
    for nu in 0..<mu:
      f[mu][nu] = g.fmunu(mu,nu,loop)
  toc("fmunu tensor")
  return f

proc reTrMul(x,y:auto):auto =
  var d: type(eval(toDouble(redot(x[0],y[0]))))
  for ir in x:
    d += redot(x[ir].adj, y[ir])
  result = simdSum(d)
  x.l.threadRankSum(result)

proc densityE*(f:auto):auto =
  ## construct E from F_munu
  ## returns (E_s, E_t)
  tic("densityE")
  let nd = f[1][0].l.nDim
  var es,et:float
  for mu in 1..<nd:
    for nu in 0..<mu:
      threads:
        let t = reTrMul(f[mu][nu], f[mu][nu])
        threadSingle:
          if mu<nd-1:
            es += t
          else:
            et += t
  let vi = -1.0/f[1][0].l.physVol.float
  toc("end")
  return (vi*es,vi*et)

proc topoQ*(f:auto):auto =
  ## construct Q from F_munu
  ## returns Q
  tic("topoQ")
  var q:float
  threads:
    let
      a = reTrMul(f[1][0], f[3][2])
      b = reTrMul(f[2][0], f[3][1])
      c = reTrMul(f[2][1], f[3][0])
    threadSingle:
      q = -1.0/(4.0*PI*PI)*(a-b+c)
  toc("end")
  return q

template defaultSetup*:untyped {.dirty.} =
  bind paramCount, paramStr, isInteger, parseInt, fileExists, getFileLattice
  echo "rank ", myRank, "/", nRanks
  threads:
    echo "thread ", threadNum, "/", numThreads
  var fn:string
  var lat:seq[int]
  when declared(defaultGaugeFile):
    if fileExists(defaultGaugeFile):
      fn = defaultGaugeFile
  if paramCount()>0:
    if (not isInteger(paramStr(1))) and paramStr(1)[0]!='-':
      fn = paramStr(1)
  if fn != "":
    lat = getFileLattice(fn)
  else:
    if paramCount()>0 and isInteger(paramStr(1)):
      lat.newSeq(0)
      var pc = paramCount()
      for i in 1..pc:
        if not isInteger(paramStr(i)): break
        lat.add parseInt(paramStr(i))
    else:
      when declared(defaultLat):
        when defaultLat is array:
          lat = @defaultLat
        else:
          lat = defaultLat
      else:
        lat = @[8,8,8,8]
  var rankGeom = intSeqParam("rg", newSeq[int](0))
  var lo = newLayout(lat, rankGeom)
  var g = newSeq[type(lo.ColorMatrix())](lat.len)
  echo "Gauge field float type: ", $(g[0].numberType)
  for i in 0..<lat.len:
    g[i] = lo.ColorMatrix()
  if fn != "":
    let status = g.loadGauge(fn)
    if status!=0:
      echo "ERROR: couldn't load gauge file ", fn
      qexFinalize()
      quit(-1)
  else:
    for i in 0..<lat.len:
      g[i] := 1

proc projectU*(x:Field) =
  for i in x: x[i].projectU
proc projectU*(x:Field, y:Field) =
  for i in x: x[i].projectU y[i]

proc projectU*[F:Field](x: openArray[F]) =
  for i in x.low..x.high: x[i].projectU
proc projectU*[F:Field](x: openArray[F], y: openArray[F]) =
  for i in x.low..x.high: x[i].projectU y[i]

proc projectSU*(x:Field) =
  for i in x: x[i].projectSU
proc projectSU*(x:Field, y:Field) =
  for i in x: x[i].projectSU y[i]

proc projectSU*[F:Field](x: openArray[F]) =
  for i in x.low..x.high: x[i].projectSU
proc projectSU*[F:Field](x: openArray[F], y: openArray[F]) =
  for i in x.low..x.high: x[i].projectSU y[i]

proc projectTAH*(x:Field) =
  for i in x: x[i].projectTAH
proc projectTAH*(x:Field, y:Field) =
  for i in x: x[i].projectTAH y[i]

proc projectTAH*[F:Field](x: openArray[F]) =
  for i in x.low..x.high: x[i].projectTAH
proc projectTAH*[F:Field](x: openArray[F], y: openArray[F]) =
  for i in x.low..x.high: x[i].projectTAH y[i]

proc randomU*(x: Field, r: var RNGField) =
  x.gaussian r
  #threadBarrier()
  #echo "x.gaus: ", x.norm2
  #threadBarrier()
  x.projectU
  #threadBarrier()
  #echo "x.proj: ", x.norm2

proc randomSU*(x: Field, r: var RNGField) =
  x.gaussian r
  x.projectSU

proc randTah3(m: var auto, s: var auto) =
  let s2 = 0.70710678118654752440;  # sqrt(1/2)
  let s3 = 0.57735026918962576450;  # sqrt(1/3)
  let r3 = s2 * gaussian(s)
  let r8 = s2 * s3 * gaussian(s)
  m[0,0].set 0, r8+r3
  m[1,1].set 0, r8-r3
  m[2,2].set 0, -2*r8
  let r01 = s2 * gaussian(s)
  let r02 = s2 * gaussian(s)
  let r12 = s2 * gaussian(s)
  let i01 = s2 * gaussian(s)
  let i02 = s2 * gaussian(s)
  let i12 = s2 * gaussian(s)
  m[0,1].set  r01, i01
  m[1,0].set -r01, i01
  m[0,2].set  r02, i02
  m[2,0].set -r02, i02
  m[1,2].set  r12, i12
  m[2,1].set -r12, i12

proc randomTAH*(x: Field, r: var RNGField) =
  when x[0].nrows == 3:
    for i in x.sites:
      randTah3(x{i}, r[i])
  else:
    x.gaussian r
    x.projectTAH

proc checkU*[F:Field](x: openArray[F]): tuple[avg,max:float] {.noinit.} =
  var a,b:float
  for mu in x.low..x.high:
    for s in x[mu]:
      let d = x[mu][s].checkU
      a += d.simdSum
      let m = d.simdMax
      if b < m: b = m
  threadRankSum a
  threadRankMax b
  const nc = x[0][0].nrows
  let vol = x[0].l.physVol
  let c = float(2*(nc*nc+1))
  a = sqrt( a / (c*float(x.len*vol)) )
  b = sqrt( b / c )
  return (a, b)

proc checkSU*[F:Field](x: openArray[F]): tuple[avg,max:float] {.noinit.} =
  var a,b:float
  for mu in x.low..x.high:
    for s in x[mu]:
      let d = x[mu][s].checkSU
      a += d.simdSum
      let m = d.simdMax
      if b < m: b = m
  threadRankSum a
  threadRankMax b
  const nc = x[0][0].nrows
  let vol = x[0].l.physVol
  let c = float(2*(nc*nc+1))
  a = sqrt( a / (c*float(x.len*vol)) )
  b = sqrt( b / c )
  return (a, b)

proc random*[F:Field](g: openArray[F], r: var RNGField) =
  for mu in g.low..g.high:
    when g[mu][0].nrows==1:
      randomU(g[mu], r)
    else:
      randomSU(g[mu], r)

proc warm*[F:Field](g: openArray[F], s: float, r: var RNGField) =
  for mu in g.low..g.high:
    when g[mu][0].nrows==1:
      g[mu].gaussian r
      g[mu] := (1-s) + s*g[mu]
      g[mu].projectU
    else:
      g[mu].gaussian r
      g[mu] := (1-s) + s*g[mu]
      g[mu].projectSU

proc random*(g: array or seq) =
  var r = newRNGField(RngMilc6, g[0].l)
  threads:
    g.random r

proc unit*(g: array or seq) =
  for i in 0..<g.len:
    g[i] := 1

proc randomTAH*[F:Field](g: openArray[F], r: var RNGField) =
  for mu in g.low..g.high:
    randomTAH(g[mu], r)

proc setupLattice*(lat:openarray[int]):auto =
  var
    lat:seq[int] = @lat
    fn = ""
  let pc = paramCount()
  if pc > 0:
    if paramStr(1).isInteger:
      lat = @[]
      for i in 1..pc:
        if not paramStr(i).isInteger: break
        lat.add paramStr(i).parseInt
    elif paramStr(1)[0] != '-':
      fn = paramStr(1)
      lat = fn.getFileLattice
      if lat.len == 0:
        echo "ERROR: getFileLattice failed on '", fn, "'"
        qexExit 1
  var
    lo = lat.newLayout
    g = lo.newGauge
    r = newRNGField(RngMilc6, lo, intParam("seed", 823543).uint64)
  if fn.len > 0:
    if 0 != g.loadGauge(fn):
      echo "ERROR: loadGauge failed on '", fn, "'"
      qexExit 1
  else:
    threads: g.random r
  return (lo, g, r)


when isMainModule:
  import qex
  import physics/qcdTypes
  qexInit()
  #let defaultLat = @[2,2,2,2]
  let defaultLat = @[8,8,8,8]
  defaultSetup()
  g.random

  var pl2 = plaq2(g)
  echo pl2

  var pl3 = plaq3(g)
  echo pl3

  var pl = plaq(g)
  echo pl
  echo pl.sum
  echo pl.mapIt(it*float(lo.nDim*(lo.nDim-1) div 2))

  echo g.wilsonLines [
    @[1,2,-1,-2],
    @[1,3,-1,-3],
    @[2,3,-2,-3],
    @[1,4,-1,-4],
    @[2,4,-2,-4],
    @[3,4,-3,-4],
  ]

  echoTimers()
  resetTimers()

  let wl = g.wilsonLines [
    @[2,-1,-2,1],
    @[3,-1,-3,1],
    @[3,-2,-3,2],
    @[4,-1,-4,1],
    @[4,-2,-4,2],
    @[4,-3,-4,3],
    @[-1,-2,1,2],
    @[-1,-3,1,3],
    @[-2,-3,2,3],
    @[-1,-4,1,4],
    @[-2,-4,2,4],
    @[-3,-4,3,4],
    @[-2,1,2,-1],
    @[-3,1,3,-1],
    @[-3,2,3,-2],
    @[-4,1,4,-1],
    @[-4,2,4,-2],
    @[-4,3,4,-3],
  ]
  for i in 0..<3:
    var r = newseq[typeof wl[0]]()
    for j in 0..<6:
      r.add wl[j+6*i]
    echo r

  #[
  var st = lo.newGauge()
  for gg in g: gg := 0.5
  for s in st: s := 0
  var c = [[1,1,1,1],[1,1,1,1],[1,1,1,1],[1,1,1,1]]
  for i in 0..<lo.nDim:
    echo st[i].norm2
  staples(st, g, g, c, c)
  for i in 0..<lo.nDim:
    echo st[i].norm2
  ]#

  echoTimers()
  qexFinalize()
