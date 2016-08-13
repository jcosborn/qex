import qex
import qcdTypes
import times
import os
import strUtils
import stdUtils
import metaUtils
import random
import profile

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
  return 0

template makeShifts(f:untyped):untyped {.dirty.} =
  template f(mu:int; b:expr):expr =
    tc = (tc+1) mod nf
    shiftFwd(tf[tc], mu, b)
    tf[tc]
  template `f "d"`(mu:int; b:expr):expr {.dirty.} =
    tc = (tc+1) mod nf
    shiftBck(tf[tc], mu, b)
    tf[tc]
template makeTransporters(f,g:untyped):untyped {.dirty.} =
  template f(mu:int; b:expr):expr =
    tc = (tc+1) mod nf
    transportFwd(tf[tc], g[mu], mu, b)
    tf[tc]
  template `f "d"`(mu:int; b:expr):expr {.dirty.} =
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


proc plaq*[T](uu: openArray[T]): auto =
  tic()
  template getIp(mu,nu: int): int = ((mu*(mu-1)) div 2) + nu
  let u = cast[ptr cArray[T]](unsafeAddr(uu[0]))
  let lo = u[0].l
  let nd = lo.nDim
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
    var plt = newSeq[float64](np)
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
      for i,v in pairs(plt):
        pl[i] = v/(lo.physVol.float*0.5*float(nd*(nd-1)*nc))
      rankSum(pl)
    toc("plaq sum")
  result = pl
  toc("plaq end", flops=lo.nSites.float*6*(2*66+36))

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
  tic()
  let g = cast[ptr cArray[T]](unsafeAddr(gg[0]))
  let lo = g[0].l
  let nd = lo.nDim
  var m = lo.ColorMatrix()
  var s0 = lo.ColorMatrix()
  var t0 = lo.ColorMatrix()
  var s1 = lo.ColorMatrix()
  var t1 = lo.ColorMatrix()
  var tr:type(trace(m))
  toc("plaq2 setup")
  threads:
    tic()
    m := 0
    toc("plaq2 zero")
    for mu in 1..<nd:
      for nu in 0..<mu:
        tic()
        shift(s0, mu,1, g[nu])
        toc("plaq2 shift1")
        shift(s1, nu,1, g[mu])
        toc("plaq2 shift2")
        #echo "s0: ", trace(s0)
        #echo "s1: ", trace(s1)
        m += (g[mu]*s0) * (g[nu]*s1).adj
        #m += (g[mu]*s0) * (g[nu]*s1)
        #echo mu, " ", nu, " ", trace(m)/nc
        toc("plaq2 mul")
    toc("plaq2 work")
    tr = trace(m)
    toc("plaq2 trace")
  toc("plaq2 threads")
  result = tr/(lo.physVol.float*0.5*float(nd*(nd-1)*nc))

proc newGauge*(l:Layout):auto =
  result = newSeq[type(l.ColorMatrix())](l.nDim)
  for i in 0..<l.nDim:
    result[i] = l.ColorMatrix()

template defaultSetup*:untyped {.dirty.} =
  bind paramCount, paramStr, isDigit, parseInt, fileExists
  echo "rank ", myRank, "/", nRanks
  threads:
    echo "thread ", threadNum, "/", numThreads
  var fn:string
  var lat:seq[int]
  when declared(defaultGaugeFile):
    if fileExists(defaultGaugeFile):
      fn = defaultGaugeFile
  if paramCount()>0:
    if not isDigit(paramStr(1)):
      fn = paramStr(1)
  if fn != nil:
    lat = getFileLattice(fn)
  else:
    if paramCount()>0 and isDigit(paramStr(1)):
      var pc = paramCount()
      lat.newSeq(pc)
      for i in 0..<pc:
        lat[i] = parseInt(paramStr(i+1))
    else:
      when declared(defaultLat):
        when defaultLat is array:
          lat = @defaultLat
        else:
          lat = defaultLat
      else:
        lat = @[8,8,8,8]
  var lo = newLayout(lat)
  var g = newSeq[type(lo.ColorMatrix())](lat.len)
  for i in 0..<lat.len:
    g[i] = lo.ColorMatrix()
  if fn != nil:
    let status = g.loadGauge(fn)
    if status!=0:
      echo "ERROR: couldn't load gauge file ", fn
      qexFinalize()
      quit(-1)
  else:
    for i in 0..<lat.len:
      g[i] := 1

proc random*(g:Field) =
  let l = g.l
  let nrm1 = g[0].nrows - 1
  let ncm1 = g[0].ncols - 1
  tfor s, 0..<l.nSites:
    for ic in countup(0,nrm1):
      for jc in countup(0,ncm1):
        g{s}[ic,jc].re = random(2.0) - 1.0
        g{s}[ic,jc].im = random(2.0) - 1.0
proc random*[T](g:openArray[T]) =
  for mu in 0..<g.len:
    g[mu].random

when isMainModule:
  qexInit()
  #let defaultLat = @[2,2,2,2]
  let defaultLat = @[8,8,8,8]
  defaultSetup()
  g.random

  var pl2 = plaq2(g)
  echo pl2

  var pl = plaq(g)
  echo pl
  echo pl.sum

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

  qexFinalize()
