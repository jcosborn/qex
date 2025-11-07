import qex
import sequtils, json, strutils, io, sets
import gauge
import gauge/gaugeUtils

proc reunit*(g:auto) =
  tic()
  threads:
    let d = g.checkSU
    threadBarrier()
    echo "unitary deviation avg: ",d.avg," max: ",d.max
    g.projectSU
    threadBarrier()
    let dd = g.checkSU
    echo "new unitary deviation avg: ",dd.avg," max: ",dd.max
  toc("reunit")

proc mplaq*(g:auto) =
  tic()
  let
    pl = g.plaq
    nl = pl.len div 2
    ps = pl[0..<nl].sum * 2.0
    pt = pl[nl..^1].sum * 2.0
  echo "MEASplaq ss: ",ps,"  st: ",pt,"  tot: ",0.5*(ps+pt)
  toc("plaq")

proc ploop*(g:auto) =
  tic()
  let pg = g[0].l.physGeom
  var pl = newseq[typeof(g.wline @[1])](pg.len)
  for i in 0..<pg.len:
    pl[i] = g.wline repeat(i+1, pg[i])
  let
    pls = pl[0..^2].sum / float(pl.len-1)
    plt = pl[^1]
  echo "MEASploop spatial: ",pls.re," ",pls.im," temporal: ",plt.re," ",plt.im
  toc("ploop")

proc reTrMul(x,y:auto):auto =
  var d: type(eval(toDouble(redot(x[0],y[0]))))
  for ir in x:
    d += redot(x[ir].adj, y[ir])
  result = simdSum(d)
  x.l.threadRankSum(result)

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

proc rmunu*(g:auto, mu,nu:int, loop=1):auto =
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
  let 
    f = newOneOf(fs[0])
    r = newOneOf(fs[0])
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
    #f.projectTAH
    threadBarrier()
    r := 0.5*(f-f.adj)
  toc("single fmunu")
  return r

proc rmunu*(g:auto, loop=1):auto =
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

proc leviCivita(a,n,r,s: int): int =
  var
    idx = @[a,n,r,s]
    swapCount = 0
  let idxSet = idx.toSet

  if idxSet.len != 4: return 0

  for idx1 in 0..<idx.len:
    for idx2 in idx1+1..<idx.len:
      if idx[idx1] > idx[idx2]:
        let tmp = idx[idx1]
        idx[idx1] = idx[idx2]
        idx[idx2] = tmp
        swapCount.inc

  result = case (swapCount mod 2 == 0):
    of true: 1
    of false: -1

proc topologicalAction*(gc: GaugeActionCoeffs; u: auto): float =
  let 
    f = u.rmunu
    betaq = gc.plaq
  var action: float
  threads:
    let
      a = reTrMul(f[1][0], f[3][2])
      b = reTrMul(f[2][0], f[3][1])
      c = reTrMul(f[2][1], f[3][0])
    threadMaster: action = -betaq*(a-b+c)
  result = action

proc topologicalDerivative*(
    gc: GaugeActionCoeffs; 
    u,f: auto; 
    project: bool = false
  ) =
  # This code is ugly as can be, but it does the job
  const
    FWD =  1
    BWD = -1

  let
    rrs = u.rmunu

  var
    uShift: Shifter[type(u[0]),type(u[0][0])]
    rShift: Shifter[type(rrs[0][0]),type(rrs[0][0][0])]

  var
    up = newSeq[seq[type(uShift)]](u.len) # U_alpha/nu(n+nu/alpha)
    un = newSeq[seq[type(uShift)]](u.len) # U_alpha/nu(n-nu/alpha)
    upn = newSeq[seq[type(uShift)]](u.len) # U_alpha/nu(n+nu/alpha-alpha/nu)

    rp = newSeq[seq[seq[type(rShift)]]](rrs.len) # R_{rho,sigma}(n+nu/alpha)
    rn = newSeq[seq[seq[type(rShift)]]](rrs.len) # R_{rho,sigma}(n-nu/alpha)

    rpp = newSeq[seq[seq[seq[type(rShift)]]]](rrs.len) # R_{rho,sigma}(n+alpha+nu)
    rpn = newSeq[seq[seq[seq[type(rShift)]]]](rrs.len) # R_{rho,sigma}(n+alpha-nu)

  # Single shifts
  for a in 0..<u.len:
    # Gauge field shifts
    up[a] = newSeq[type(uShift)](u.len)
    un[a] = newSeq[type(uShift)](u.len)
    for n in 0..<u.len:
      up[a][n] = newShifter(u[a],n,FWD)
      un[a][n] = newShifter(u[a],n,BWD)
      threads:
        discard up[a][n] ^* u[a] # U_alpha/nu(n+nu/alpha)
        discard un[a][n] ^* u[a] # U_alpha/nu(n-nu/alpha)

    # Clover shifts
    rp[a] = newSeq[seq[type(rShift)]](u.len-1)
    rn[a] = newSeq[seq[type(rShift)]](u.len-1)
    for r in 1..<rrs.len:
      rp[a][r-1] = newSeq[type(rShift)](r)
      rn[a][r-1] = newSeq[type(rShift)](r)
      for s in 0..<r:
        rp[a][r-1][s] = newShifter(rrs[r][s],a,FWD)
        rn[a][r-1][s] = newShifter(rrs[r][s],a,BWD)
        threads:
          discard rp[a][r-1][s] ^* rrs[r][s] # R_{rho,sigma}(n+nu/alpha)
          discard rn[a][r-1][s] ^* rrs[r][s] # R_{rho,sigma}(n-nu/alpha)

  # Double shifts
  for a in 0..<u.len:
    upn[a] = newSeq[type(uShift)](u.len)
    for n in 0..<u.len:
      upn[a][n] = newShifter(up[a][n].field,a,BWD)
      threads:
        discard upn[a][n] ^* up[a][n].field # U_alpha/nu(n+nu/alpha-alpha/nu)
    rpp[a] = newSeq[seq[seq[type(rShift)]]](u.len)
    rpn[a] = newSeq[seq[seq[type(rShift)]]](u.len)
    for n in 0..<u.len:
      rpp[a][n] = newSeq[seq[type(rShift)]](u.len-1)
      rpn[a][n] = newSeq[seq[type(rShift)]](u.len-1)
      if n != a:
        for r in 1..<rrs.len:
          rpp[a][n][r-1] = newSeq[type(rShift)](r)
          rpn[a][n][r-1] = newSeq[type(rShift)](r)
          for s in 0..<r:
            rpp[a][n][r-1][s] = newShifter(rp[a][r-1][s].field,n,FWD)
            rpn[a][n][r-1][s] = newShifter(rp[a][r-1][s].field,n,BWD)
            threads:
              discard rpp[a][n][r-1][s] ^* rp[a][r-1][s].field # R_{r,s}(n+a+n)
              discard rpn[a][n][r-1][s] ^* rp[a][r-1][s].field # R_{r,s}(n+a-n)

  # Calculate A_{alpha,nu,rho,sigma}
  threads:
    for a in 0..<u.len:
      for n in 0..<u.len:
        if n != a:
          for r in 1..<u.len:
            for s in 0..<r:
              if (r != a) and (r != n) and (s != a) and (s != n):
                let eps = gc.plaq*leviCivita(a,n,r,s).float/16.0
                for st in u[0]:
                  var 
                    f1 {.noinit.}:  type(u[0][0])
                    f11 {.noinit.}: type(u[0][0])
                    f12 {.noinit.}: type(u[0][0])

                    f2 {.noinit.}: type(u[0][0])
                    f21 {.noinit.}: type(u[0][0])
                    f22 {.noinit.}: type(u[0][0])

                    f3 {.noinit.}: type(u[0][0])
                    f31 {.noinit.}: type(u[0][0])

                    fmat {.noinit.}: type(u[0][0])

                  fmat := 0

                  case project:
                    of true: f21 := u[a][st] * upn[a][n].field[st].adj
                    of false: f21 := upn[a][n].field[st].adj
                  f11 :=  u[a][st]*up[n][a].field[st]
                  f12 := (u[n][st]*up[a][n].field[st]).adj
                  f1 := f12 * rrs[r][s][st]
                  f1 += up[a][n].field[st].adj * rp[n][r-1][s].field[st] * u[n][st].adj
                  f1 += rpp[a][n][r-1][s].field[st] * f12
                  fmat += f11*f1

                  case project:
                    of true: f21 := u[a][st] * upn[a][n].field[st].adj
                    of false: f21 := upn[a][n].field[st].adj
                  f22 := un[a][n].field[st].adj * un[n][n].field[st]
                  f2 := f22 * rrs[r][s][st]
                  f2 += un[a][n].field[st].adj * rn[n][r-1][s].field[st] * un[n][n].field[st]
                  f2 += rpn[a][n][r-1][s].field[st] * f22
                  fmat -= f21*f2

                  case project:
                    of true: f31 := u[a][st] * rp[a][r-1][s].field[st]
                    of false: f31 := rp[a][r-1][s].field[st]
                  f3 := up[n][a].field[st] * f12
                  f3 -= upn[a][n].field[st].adj * f22
                  fmat += f31*f3

                  threadBarrier()

                  case project:
                    of true: f[a][st].projectTAH(eps*fmat)
                    of false: f[a][st] := eps*fmat

proc setGauge*(u: auto; v: auto) =
  threads:
    for mu in 0..<u.len:
      for s in u[mu]: u[mu][s] := v[mu][s]

template gradientFlow(
    gc: GaugeActionCoeffs; 
    g: array|seq; 
    steps: int;
    eps: float; 
    measure: untyped
  ): untyped =
  #[ Gradient flow w/ Wilson or rectangle action
     Originally written by James Osborn & Xiaoyong Jin.
     d/dt Vt = Z(Vt) Vt
     Runge-Kutta:
     W0 <- Vt
     W1 <- exp(1/4 Z0) W0
     W2 <- exp(8/9 Z1 - 17/36 Z0) W1
     V(t+eps) <- exp(3/4 Z2 - 8/9 Z1 + 17/36 Z0) W2
     where
     Zi = eps Z(Wi)
  ]#
  proc flowProc {.gensym.} =
    tic("flowProc")
    const nc = g[0][0].nrows.float
    var
      p = g[0].l.newGauge  # mom
      f = g[0].l.newGauge  # force
      n = 1
    while true:
      let t = n * eps
      let epsnc = eps * nc  # compensate force normalization
      gc.gaugeForce(g,f)
      threads:
        for mu in 0..<f.len:
          for e in g[mu]:
            var v {.noinit.}:type(load1(f[0][0]))
            v := (-1.0/4.0)*epsnc*f[mu][e]
            let t = exp(v)*g[mu][e]
            p[mu][e] := v
            g[mu][e] := t
      gc.gaugeForce(g,f)
      threads:
        for mu in 0..<f.len:
          for e in g[mu]:
            var v {.noinit.}:type(load1(f[0][0]))
            v := (-8.0/9.0)*epsnc*f[mu][e] + (-17.0/9.0)*p[mu][e]
            let t = exp(v)*g[mu][e]
            p[mu][e] := v
            g[mu][e] := t
      gc.gaugeForce(g,f)
      threads:
        for mu in 0..<f.len:
          for e in g[mu]:
            var v {.noinit.}:type(load1(f[0][0]))
            v := (-3.0/4.0)*epsnc*f[mu][e] - p[mu][e]
            let t = exp(v)*g[mu][e]
            g[mu][e] := t
      let wflowT {.inject.} = t
      measure
      inc n
      if steps>0 and n>steps: break
    toc("end")
  flowProc()

template gradientFlow(
    gc: GaugeActionCoeffs;
    g: array|seq; 
    eps: float; 
    measure: untyped
  ): untyped = 
  gc.gradientFlow(g,0,eps):
    let flowTime {.inject.} = wflowT
    measure

proc flowMeasurements(u: auto; loop: int; tau: float): JsonNode =
  var
    pls,plt: ComplexProxy[ComplexObj[float64,float64]]
    poly: seq[ComplexProxy[ComplexObj[float64,float64]]]
    t2Ess,t2Est,t2Ees,t2Eet: float
  let
    f = u.fmunu(loop)
    (es, et) = f.densityE
    q = f.topoQ
    pl = u.plaq
    nl = pl.len div 2
    ss = 6.0*pl[0..<nl].sum
    st = 6.0*pl[nl..^1].sum
    pg = u[0].l.physGeom
  poly = newSeq[ComplexProxy[ComplexObj[float64,float64]]](pg.len)
  for i in 0..<pg.len: poly[i] = u.wline repeat(i+1,pg[i])
  pls = poly[0..^2].sum/float(poly.len-1)
  plt = poly[^1]
  (t2Ess,t2Est) = (6.0*tau*tau*(3.0-ss),6.0*tau*tau*(3.0-st))
  (t2Ees,t2Eet) = (tau*tau*es,tau*tau*et)
  result = %* {
    "flow-time": tau,
    "plaquette":0.5*ss+0.5*st,
    "clover":es+et,
    "t2E-plaquette":t2Ess+t2Est,
    "t2E-clover":t2Ees+t2Eet,
    "t2E-spacelike-plaquette":t2Ess,
    "t2E-timelike-plaquette":t2Est,
    "t2E-spacelike-clover":t2Ees,
    "t2E-timelike-clover":t2Eet,
    "topological-charge":q,
    "Re(spacelike-Polyakov-loop)":3.0*pls.re,
    "Im(spacelike-Polyakov-loop)":3.0*pls.im,
    "Re(timelike-Polyakov-loop)":3.0*plt.re,
    "Im(timelike-Polyakov-loop)":3.0*plt.im,
  }

proc get(info:JsonNode;key:string): seq[float] = 
  result = newSeq[float]()
  for el in info[key].getElems(): result.add getFloat(el)

template gradientFlow*(u: auto; info: JsonNode; body: untyped) =
  var 
    v = u[0].l.newGauge
    f {.inject,used.}: File
    tau {.inject.}: float
    measurements {.inject.}: JsonNode
  for flow,flowInfo in info:
    let
      loops = case info[flow].hasKey("loops")
        of true: info[flow]["loops"].getInt()
        of false: 1
      beta = case info[flow].hasKey("beta")
        of true: info[flow]["beta"].getFloat()
        of false: 1.0
      cr = case info[flow]["action"].getStr()
        of "Rectangle": info[flow]["cr"].getFloat()
        else: 0.0
      gc = case info[flow]["action"].getStr()
        of "Rectangle": gaugeActRect(beta,cr)
        else: GaugeActionCoeffs(plaq:beta)
      fn = info[flow]["path"].getStr() & info[flow]["filename"].getStr()
      dts = info[flow].get("time-increments")
      maxFlts = info[flow].get("maximum-flow-times")
    var lastMaxFlt = 0.0
    threads: v := u
    v.reunit
    f = fn.open(fmWrite)
    for (dt,maxFlt) in zip(dts,maxFlts):
      gc.gradientFlow(v,dt):
        tau = flowTime + lastMaxFlt
        measurements = v.flowMeasurements(loops,tau)
        body
        if flowTime>=maxFlt: break
      lastMaxFlt = maxFlt
    f.close()
      
proc formatMeasurements*(
    measurements: JsonNode;
    style: string = "default"
  ): string =
  result = ""
  let tau = measurements["flow-time"].getFloat()
  case style:
    of "default":
      var concatenation = @["FLOW",tau.formatFloat(ffDecimal,3)]
      result = concatenation.join(" ")
      for key,measurement in measurements:
        let O = measurement.getFloat()
        concatenation = @[result,O.formatFloat(ffDecimal,13)] 
        result = concatenation.join(" ")
    of "KS_nHYP_FA": # Style of https://github.com/daschaich/KS_nHYP_FA
      let
        standIn = %* {"stand-in": 0.0}
        observables = @[
          measurements["plaquette"],
          measurements["clover"],
          measurements["t2E-clover"],
          standIn["stand-in"],
          measurements["t2E-plaquette"],
          measurements["topological-charge"],
          measurements["t2E-spacelike-clover"],
          measurements["t2E-timelike-clover"],
          measurements["Re(timelike-Polyakov-loop)"],
          measurements["Im(timelike-Polyakov-loop)"],
          measurements["Re(spacelike-Polyakov-loop)"],
          measurements["Im(spacelike-Polyakov-loop)"]
        ]
      var concatenation = @["FLOW",tau.formatFloat(ffDecimal,2)]
      result = concatenation.join(" ")
      for observable in observables:
        let O = observable.getFloat()
        concatenation = @[result,O.formatFloat(ffDecimal,13)] 
        result = concatenation.join(" ")
    else: discard

proc setFilename*(info: auto; fn: string) =
  info["filename"] = %* fn 

if isMainModule:
  qexInit()
  let
    lo = newLayout(@[8,8,8,8])
    gc = GaugeActionCoeffs(plaq:6.0)
  var 
    u = lo.newGauge()
    f = lo.newGauge()
  u.random()
  echo "Levi-Civita 1234: ", leviCivita(1,2,3,4)
  echo "Levi-Civita 2134: ", leviCivita(2,1,3,4)
  echo "Levi-Civita 2143: ", leviCivita(2,1,4,3)
  echo "Levi-Civita 4123: ", leviCivita(4,2,3,1)
  echo "Levi-Civita 1123: ", leviCivita(1,1,2,3)
  echo "action: ", gc.topologicalAction(u)
  gc.topologicalDerivative(f,u)
  for mu in 0..<f.len: echo "|f[",mu,"]|^2: ",f[mu].norm2
  qexFinalize()