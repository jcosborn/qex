import qex
import gauge, physics/qcdTypes
import os, strutils, times

const QEXDEBUG {.booldefine.} = false
template debugPut(n:untyped) =
  when QEXDEBUG:
    echo asttostr(n),": ",n

proc sumEnergy(fr,fi:auto, J, h:auto, g:auto, sf,sb:auto) =
  fr := 0
  fi := 0
  for nu in 0..<g.l.nDim:
    discard sf[nu] ^* g
    discard sb[nu] ^* g
    threadBarrier()
    for i in fr:
      fr[i] += cos(sf[nu].field[i]) + cos(sb[nu].field[i])
    for i in fi:
      fi[i] += sin(sf[nu].field[i]) + sin(sb[nu].field[i])
  fr := J*fr + h
  fi := J*fi

const NMax = 32

type GlobalP = object
  n:int
  N,N2,N3:float
  sn,cn,pidn:float

proc init(globalP:var GlobalP, N:float) =
  if N<=0:
    qexError "N = " & $N
  globalP.N = N
  globalP.N2 = N*N
  globalP.N3 = N*N*N
  globalP.n = N.ceil.int
  let ck = globalP.n mod 2
  globalP.sn = if ck == 0: -1.0 else: 1.0
  globalP.cn = float(ck-globalP.n)
  globalP.pidn = PI/N
  if 2*globalP.n+3>NMax:
    let n = floor((NMax-3).float/2.0)
    qexError "N = " & $N & " exceeds maximum " & $n

var globalP {.align(64).}:GlobalP

type
  ExpDist[N:static[int],D] = object
    ## exp(bi+ai*x) for x in [ci,ci+di], i in 0..<n, n<N
    ## d[i] = c[i+1]-c[i]
    ## ei = exp(ai*di)-1
    ## zi = integrate(f(x), x, c0, ci)
    a,b,c,e,z:array[N,D]
    n:int

proc linearParam[D](lambda:D):auto =
  var
    a {.noinit.}:D
    b {.noinit.}:D
  if lambda>1.904538388056459:
    # optimal in the limit of large lambda
    a = sqrt(lambda)
    b = a*arcsin(a/lambda)+sqrt(lambda*(lambda-1.0));
  else:
    # keep the line above cos(x), b-a*%pi > cos(%pi) = -1
    a = 0.7246113537767085*lambda
    b = 1.276433705732662*lambda
  (-a,b)  # convention: exp(b+a*x)

proc prepareCosCosN[D](dist:var ExpDist, lambda, phi, sigma:D) =
  let
    N = globalP.N
    (a1,b1) = linearParam(lambda)
  var (an,bn) = linearParam(sigma)
  an *= N
  # linear regions for: sigma*cos(N*x)
  # -pi, .. , -3pi/N, -2pi/N, -pi/N, 0, pi/N, 2pi/N, 3pi/N, .. , pi
  # peaks: .., -2pi/N, 0, 2pi/N, ..
  # linear regions for: lambda*cos(x-phi)
  # if phi>0: -pi \ phi-pi / phi \ pi
  # else: -pi / phi \ phi+pi / pi
  var
    c1,phia:float
    s1:float
  if phi<0:
    c1 = phi
    phia = phi
    s1 = -1.0
  else:
    c1 = phi-2.0*PI
    phia = phi-PI
    s1 = 1.0
  let
    n = globalP.n
    pidn = globalP.pidn
    phib = phia+PI
    ka = floor(phia/pidn).int
    kb = floor(phib/pidn).int
  var
    sn = globalP.sn
    cn = globalP.cn
  debugput (a1,b1,an,bn)
  debugput (phia,phib)
  debugput (ka,kb)
  template ak:untyped = s1*a1+sn*an
  template bk:untyped = b1 - s1*a1*c1 + bn - sn*an*cn*pidn
  # populate c
  for k in 0..n+ka:
    sn = -sn
    if k-n>cn.int: cn += 2
    dist.a[k] = ak()
    dist.b[k] = bk()
    dist.c[k] = float(k-n)*pidn
  s1 = -s1
  if phi>=0: c1 = phi
  dist.a[n+ka+1] = ak()
  dist.b[n+ka+1] = bk()
  dist.c[n+ka+1] = phia
  for k in n+ka+2..n+kb+1:
    sn = -sn
    if k-(n+1)>cn.int: cn += 2
    dist.a[k] = ak()
    dist.b[k] = bk()
    dist.c[k] = float(k-(n+1))*pidn
  s1 = -s1
  if phi<0: c1 = phi+2.0*PI
  dist.a[n+kb+2] = ak()
  dist.b[n+kb+2] = bk()
  dist.c[n+kb+2] = phib
  for k in n+kb+3..2*n+1:
    sn = -sn
    if k-(n+2)>cn.int: cn += 2
    dist.a[k] = ak()
    dist.b[k] = bk()
    dist.c[k] = float(k-(n+2))*pidn
  # fix boundary c
  dist.c[0] = -PI
  dist.c[2*n+2] = PI
  let nn = 2*n+2
  dist.n = nn
  # integrate
  var z = 0.0
  dist.z[0] = 0.0
  for k in 0..<nn:
    let
      a = dist.a[k]
      b = dist.b[k]
      c0 = dist.c[k]
      c1 = dist.c[k+1]
      e = expm1(a*(c1-c0))
    z += e*exp(b+a*c0)/a
    dist.e[k] = e
    dist.z[k+1] = z
  z = 1.0/z
  for k in 1..<nn:
    dist.z[k] = dist.z[k]*z
  dist.z[nn] = 1.0

proc bisect[D](v:openarray[D], n:int, x:D):int =
  var
    lo = 0
    hi = n
  while hi-lo>1:
    let mid = (hi+lo) div 2
    if v[mid]<=x:
      lo = mid
    else:
      hi = mid
  lo

proc draw[D](dist:ExpDist[NMax,D], r:D):auto =
  ## r in [0,1)
  let
    n = bisect(dist.z, dist.n, r)
    zn = dist.z[n]
    zn1 = dist.z[n+1]
    u = (r-zn)/(zn1-zn)
    a = dist.a[n]
    e = dist.e[n]  # exp(ai*di)-1
    x = dist.c[n] + log1p(u*e)/a
    logpdf = dist.b[n]+a*x
  (x, logpdf)

proc expCosPlusCosN[D](rng:var RNG, lambda, phi, sigma:D):auto =
  ## sample x ~ exp(lambda*cos(x-phi) + sigma*cos(N*x))
  ## using piece-wise exponential distribution for rejection sampling
  let N = globalP.N
  var dist {.noinit,align(64).}:ExpDist[NMax,D]
  dist.prepareCosCosN(lambda, phi, sigma)
  while true:
    let
      r = rng.uniform
      u = rng.uniform
      (x, logpdf) = dist.draw r
    if u < exp(lambda*cos(x-phi) + sigma*cos(N*x) - logpdf):
      return x

template findRoot(x0,x1:SomeNumber, fexpr,dfexpr:untyped):untyped =
  ## Newton-Raphson with bounds x0 and x1.
  ## Requires f(x0)<=0 and f(x1)>=0.
  ## f and df should use x, which is injected.
  block:
    var
      lo = x0
      hi = x1
      x {.inject.} = 0.5*(lo+hi)
      dxp = abs(hi-lo)
      dx = dxp
      f = fexpr
      df = dfexpr
    # debugput ("findRoot",lo,hi)
    while true:
      # debugput (lo,x,hi)
      if ((x-hi)*df-f)*((x-lo)*df-f) > 0.0 or
          abs(2.0*f) > abs(dxp*df):
        dxp = dx
        dx = 0.5*(hi-lo)
        x = lo+dx
        # debugput (lo,x,dx)
        if lo == x or hi == x: break  # precision loss
      else:
        dxp = dx
        dx = f/df
        let t = x
        x -= dx
        # debugput (t,x,dx)
        if t == x: break  # precision loss
      if abs(dx) == 0.0: break  # or use a small value?
      f = fexpr
      df = dfexpr
      if f<0.0:
        lo=x
      else:
        hi=x
    x
template findRoot(x0,fx0,x1,fx1:SomeNumber, fexpr,dfexpr:untyped):untyped =
  ## Newton-Raphson with bounds x0 and x1,
  ## given fx0 = f(x0) and fx1 = f(x1).
  ## fx0 and fx1 should bracket 0.
  ## f and df should use x, which is injected.
  block:
    var
      lo = x0
      hi = x1
      f0 = fx0
      f1 = fx1
    if f0>f1:
      let t = hi
      hi = lo
      lo = t
    findRoot(lo,hi,fexpr,dfexpr)

proc inrange[D](x,a,b:D):bool =
  (x>=a and x<b) or (x<=a and x>b)
proc inrangeExcl[D](x,a,b:D):bool =
  (x>a and x<b) or (x<a and x>b)

type CosCosN[D] = object
  xs: array[NMax,D]
  nx: int
  phi,sigma: D
  x,fx: D

proc f[D](ccn:var CosCosN[D], x:D):D =
  cos(x-ccn.phi) + ccn.sigma*cos(globalP.N*x)
proc df[D](ccn:var CosCosN[D], x:D):D =
  sin(ccn.phi-x) - globalP.N*ccn.sigma*sin(globalP.N*x)
proc ddf[D](ccn:var CosCosN[D], x:D):D =
  -cos(ccn.phi-x) - globalP.N2*ccn.sigma*cos(globalP.N*x)
proc dddf[D](ccn:var CosCosN[D], x:D):D =
  sin(x-ccn.phi) + globalP.N3*ccn.sigma*sin(globalP.N*x)

proc findEqual[D](ccn:var CosCosN[D], lo,flo,hi,fhi:D) =
  let
    x = ccn.x
    fx = ccn.fx
  debugput (lo,flo,hi,fhi,x,fx)
  if fx.inrange(flo,fhi) and (x<lo or x>hi):
    let nx = ccn.nx
    ccn.xs[nx] = findRoot(lo,flo-fx,hi,fhi-fx,ccn.f(x)-fx,ccn.df(x))
    ccn.nx = nx+1

proc searchMixed[D](ccn:var CosCosN[D], lo,flo,hi,fhi:D) =
  let
    dflo = ccn.df(lo)
    dfhi = ccn.df(hi)
  if inrangeExcl(0.0,dflo,dfhi):
    # find the extremum of f(x)
    let
      xf = findRoot(lo,dflo,hi,dfhi,ccn.df(x),ccn.ddf(x))
      fxf = ccn.f(xf)
    ccn.findEqual(lo,flo,xf,fxf)
    ccn.findEqual(xf,fxf,hi,fhi)
  else:  # dflo and dfhi has the same sign or dflo is zero.
    # find the extremum of df
    let
      xdfm = findRoot(lo,ccn.ddf(lo),hi,ccn.ddf(hi),ccn.ddf(x),ccn.dddf(x))
      dfm = ccn.df(xdfm)
    if inrangeExcl(0.0,dfm,dflo):
      # There are two zeros of df.
      # The maximum number of zeros we can have.
      # find the two extrema of f(x).
      let
        xf0 = findRoot(xdfm,dfm,lo,dflo,ccn.df(x),ccn.ddf(x))
        xf1 = findRoot(xdfm,dfm,hi,dfhi,ccn.df(x),ccn.ddf(x))
        f0 = ccn.f(xf0)
        f1 = ccn.f(xf1)
      ccn.findEqual(lo,flo,xf0,f0)
      ccn.findEqual(xf0,f0,xf1,f1)
      ccn.findEqual(xf1,f1,hi,fhi)
    else:  # otherwise monotonic
      ccn.findEqual(lo,flo,hi,fhi)

proc solveCosCosN[D](ccn:var CosCosN[D]) =
  ## find other values of x that keeps the same value of
  ## cos(x-phi) + sigma*cos(N*x)
  ## for x in [-pi,pi)
  ## Return number of roots saved in xs.
  # extrema at sin(φ - x) = Nσ sin(N x)
  # Is it possible to have two extrema within [k*π/N, (k+1)*π/N]?
  # Or two solutions in sin(φ-y/N)=Nσ sin(y) for y in [0,π)?
  # Yes, but at most two.
  let
    phi = ccn.phi
  var
    phia:float
    s1:bool
  if phi<0:
    phia = phi
    s1 = true
  else:
    phia = phi-PI
    s1 = false
  let
    n = globalP.n
    pidn = globalP.pidn
    phib = phia+PI
    ka = floor(phia/pidn).int
    kb = floor(phib/pidn).int
  var sn = globalP.sn > 0.0
  debugput (phia,phib)
  debugput (ka,kb)
  var
    lo = -PI
    flo = ccn.f lo
  template go(z:untyped):untyped =
    let
      hi = z
      fhi = ccn.f hi
    debugput (s1,sn)
    if s1 xor sn:
      ccn.searchMixed(lo,flo,hi,fhi)
    else:
      ccn.findEqual(lo,flo,hi,fhi)
    lo = hi
    flo = fhi
  # going through regions with 1st derivatives of definite signs
  for k in 1..n+ka:
    go(float(k-n)*pidn)
    sn = not sn
  go phia
  s1 = not s1
  for k in n+ka+2..n+kb+1:
    go(float(k-(n+1))*pidn)
    sn = not sn
  go phib
  s1 = not s1
  for k in n+kb+3..2*n+1:
    go(float(k-(n+2))*pidn)
    sn = not sn
  go PI

proc pickRootCosCosN[D](rng:var RNG, x, phi, sigma:D): auto =
  # We will find solutions from -π to π excluding π.
  # Our action is not always periodic, but
  # the action at -π is always the same as the action at π.
  # We exclude π as the input value to avoid overweighting ±π.
  let x = if x == PI: -PI else: x
  var ccn {.noinit,align(64).}:CosCosN[D]
  ccn.nx = 0
  ccn.phi = phi
  ccn.sigma = sigma
  ccn.x = x
  ccn.fx = ccn.f(x)
  ccn.solveCosCosN
  debugput ccn
  let n = ccn.nx
  if n == 0:
    return x
  elif n == 1:
    return ccn.xs[0]
  else:
    let r = int(floor(n * rng.uniform))
    return ccn.xs[r]

proc probZN(ss:var openarray[float], g:auto, betah:float) =
  let
    n = globalP.n
    pidn = globalP.pidn
  var s:float
  for k in 0..<n:
    let d = 2.0*pidn*k
    threads:
      var t:typeof(g[0])
      for i in g:
        t += cos(g[i]+d)
      var v = t.simdSum
      v.threadRankSum
      threadSingle: s = v
    ss[k] = s
  s = ss[0]
  ss[0] = 1.0
  for k in 1..<n:
    ss[k] = exp(betah*(ss[k]-s))

type HeatBath[F,E] = object
  fr,fi: F
  sf,sb: array[2,seq[Shifter[F,E]]]
  subs: array[2,Subset]

proc newHeatBath(lo:auto):auto =
  let
    fr = lo.Real
    fi = lo.Real
  type
    F = typeof(fr)
    E = typeof(fr[0])
  var r = HeatBath[F,E](fr:fr, fi:fi)
  const p = ["even","odd"]
  for j in 0..1:
    r.sf[j] = newseq[Shifter[F,E]](lo.nDim)
    r.sb[j] = newseq[Shifter[F,E]](lo.nDim)
    for i in 0..<lo.nDim:
      r.sf[j][i] = newShifter(fr, i, 1, p[j])
      r.sb[j][i] = newShifter(fr, i, -1, p[j])
  r.subs[0].layoutSubset(lo,"e")
  r.subs[1].layoutSubset(lo,"o")
  r

proc evolve(H:HeatBath, g:auto, gc:auto, r:auto, R:var RngMilc6, sample = true, jump = true, jumpZN = true) =
  tic("heatbath")
  let
    lo = g.l
    nd = lo.nDim
    (beta, J, h, hn) = gc
    sigma = beta*J*hn
  if H.subs.len != 2:
    qexError "HeatBath only works with even-odd subsets for now."
  threads:
    tic("threads")
    if sample:
      # sample
      for j in 0..<H.subs.len:
        let
          s = H.subs[j]
          so = H.subs[(j+1) mod 2]
        sumEnergy(H.fr[s], H.fi[s], J, h, g, H.sf[j], H.sb[j])
        threadBarrier()
        for i in g[s].sites:
          let
            yr = eval H.fr{i}
            yi = eval H.fi{i}
            lambda = beta*hypot(yi, yr)
            phi = arctan2(yi, yr)
          g{i} := expCosPlusCosN(r{i}, lambda, phi, sigma)
      threadBarrier()
      toc("sample")
    if jump:
      # over relaxation: flip
      for j in 0..<H.subs.len:
        let
          s = H.subs[j]
          so = H.subs[(j+1) mod 2]
        sumEnergy(H.fr[s], H.fi[s], J, h, g, H.sf[j], H.sb[j])
        threadBarrier()
        for i in g[s].sites:
          let
            yr = eval H.fr{i}
            yi = eval H.fi{i}
            lambda = beta*hypot(yi, yr)
            phi = arctan2(yi, yr)
          g{i} := pickRootCosCosN(r{i}, eval g{i}, phi, sigma/lambda)
      toc("flip")
  # end of the threads block above
  if jumpZN:
    let
      n = globalP.n
      pidn = globalP.pidn
    var ss {.noinit.}:array[NMax,float]
    ss.probZN(g, beta*h)
    for k in 1..<n: ss[k] += ss[k-1]
    # echo "# Cumulative Z[N] probability: ",ss
    let u = ss[n-1] * R.uniform
    var k = 0
    while u>=ss[k]: inc k
    if k>0:
      # echo "# Flip to k = ",k
      let d = 2.0*pidn*k
      threads:
        for i in g:
          g[i] += d
    toc("tryZN")
  toc("end")

proc magnet(g:auto):auto =
  tic("magnet")
  var mr,mi = 0.0
  threads:
    var t,s:typeof(g[0])
    for i in g:
      t += cos(g[i])
      s += sin(g[i])
    var
      v = t.simdSum
      u = s.simdSum
    v.threadRankSum
    u.threadRankSum
    threadSingle:
      mr = v
      mi = u
  toc("done")
  (mr, mi)

proc magnet(g:auto,k:int):auto =
  tic("magnetK")
  let d = 2.0*globalP.pidn*k
  var mr,mi = 0.0
  threads:
    var t,s,x:typeof(g[0])
    for i in g:
      x = g[i]+d
      t += cos(x)
      s += sin(x)
    var
      v = t.simdSum
      u = s.simdSum
    v.threadRankSum
    u.threadRankSum
    threadSingle:
      mr = v
      mi = u
  toc("done")
  (mr, mi)

type PhaseDiff[F,E] = object
  cosd,sind:seq[float]
  f: seq[Shifter[F,E]]

proc newPhaseDiff(g:auto):auto =
  type
    F = typeof(g)
    E = typeof(g[0])
  var r {.noinit.}: PhaseDiff[F,E]
  let nd = g.l.nDim
  r.cosd = newseq[float](nd)
  r.sind = newseq[float](nd)
  r.f = newseq[Shifter[F,E]](nd)
  for i in 0..<nd:
    r.f[i] = newShifter(g, i, 1)
  r

proc phaseDiff(del:var PhaseDiff,g:auto):auto =
  let
    # del cannot be captured by nim in threads
    f = del.f
    cosd = cast[ptr UnCheckedArray[float]](del.cosd[0].addr)
    sind = cast[ptr UnCheckedArray[float]](del.sind[0].addr)
  threads:
    for nu in 0..<g.l.nDim:
      var d,t,s:typeof(g[0])
      discard f[nu] ^* g
      threadBarrier()
      for i in g:
        d := f[nu].field[i] - g[i]
        t += cos(d)
        s += sin(d)
      var
        v = t.simdSum
        u = s.simdSum
      v.threadRankSum
      u.threadRankSum
      threadSingle:
        cosd[nu] = v
        sind[nu] = u

proc showMeasure[F,E](del:var PhaseDiff[F,E],g:F,gc:auto,label="") =
  let
    (beta, J, h, hn) = gc
    (mr,mi) = g.magnet
    v = 1.0/g.l.physVol.float
    s = (mr*mr+mi*mi)*v
  var
    ps {.noinit.}: array[NMax, float]
    mrzn = mr
    mizn = mi
    p:float = 1.0
  ps.probZN(g, beta*h)
  for k in 1..<globalP.n:
    let (mrk,mik) = g.magnet k  # mrk²+mik² = mr²+mi²
    # echo "# ",label,"magnet",k,": ",mrk," ",mik," prob: ",ps[k]
    mrzn += mrk*ps[k]
    mizn += mik*ps[k]
    p += ps[k]
  mrzn /= p
  mizn /= p
  del.phaseDiff g
  echo label,"magnet: ",mr," ",mi," ",s
  echo label,"magnetZN: ",mrzn," ",mizn
  var diff = ""
  for i in 0..<g.l.nDim:
    diff &= "CosSinDel" & $i & ": " & $(del.cosd[i]*v) & " " & $(del.sind[i]*v) & "\t"
  diff.setlen(diff.len-1)
  echo label,diff

proc hitFreq(num, freq:int):bool = freq>0 and 0==num mod freq


qexinit()
tic()

letParam:
  #lat = @[8,8,8,8]
  #lat = @[8,8,8]
  lat = @[32,32]
  #lat = @[1024,1024]
  beta = 1.0
  J = 1.0
  h = 0.0
  N = 5.0
  hn = 0.0
  sweeps = 10
  sampleFreq = 1
  jumpFreq = 1
  jumpZNFreq = 1
  measureFreq = 1
  seed:uint64 = int(1000*epochTime())

echoParams()
echo "rank ", myRank, "/", nRanks
threads: echo "thread ", threadNum, "/", numThreads

globalP.init(N)

when QEXDEBUG:
  let
    lambda = 5.0
    phi = -0.05
    sigma = 1.9
  proc testSample(lambda,phi,sigma:float) =
    var dist {.noinit,align(64).}:ExpDist[NMax,float]
    dist.prepareCosCosN(lambda, phi, sigma)
    echo dist
    for i in 0..<dist.n:
      echo i,": ",dist.a[i]," ",dist.b[i]," ",dist.c[i]
      echo "    ",dist.e[i]," ",dist.z[i]
  proc testSearch(x,lambda,phi,sigma:float) =
    var ccn {.noinit,align(64).}:CosCosN[float]
    ccn.nx = 0
    ccn.phi = phi
    ccn.sigma = sigma/lambda
    ccn.x = x
    ccn.fx = ccn.f(x)
    ccn.solveCosCosN
    debugput ccn
    for i in 0..<ccn.nx:
      let x = ccn.xs[i]
      echo i,": ",x," ",cos(x-phi)+sigma/lambda*cos(N*x)
  testSample(lambda,phi,sigma)
  testSearch(2.25,lambda,phi,sigma)
  testSearch(-2.3,lambda,phi,sigma)
  testSearch(-PI,lambda,phi,sigma)  # we always exclude π from the solutions
  testSearch(PI,lambda,phi,sigma)  # This will find -π
  testSearch(-2.855993321445266,lambda,phi,sigma)
  testSearch(-1.427996660722633,lambda,phi,sigma)
  testSearch(-0.05,lambda,phi,sigma)
  testSearch(-0.01761145434218007,lambda,phi,sigma)
  testSearch(2.106314709054424,lambda,phi,sigma)
  qexExit()

let
  gc = (beta:beta, J:J, h:h, hn:hn)
  lo = lat.newLayout
  vol = lo.physVol

var
  r = lo.newRNGField(RngMilc6, seed)
  R:RngMilc6  # global RNG for global update
  g = lo.Real
  del = newPhaseDiff g
R.seed(seed,987654321)

threads:
  for i in g.sites:
    let u = uniform r{i}
    g{i} := PI*(2.0*u-1.0)

del.showMeasure(g, gc, "Initial: ")

var H = lo.newHeatBath

toc("init")

for n in 1..sweeps:
  tic("sweep")
  echo "Begin sweep: ",n

  H.evolve(g,gc,r,R,
    hitFreq(n,sampleFreq),
    hitFreq(n,jumpFreq),
    hitFreq(n,jumpZNFreq))
  toc("evolve")

  if hitFreq(n,measureFreq): del.showMeasure(g, gc)
  toc("measure")

toc("done")
echoTimers()
qexfinalize()
