import base
#import stdUtils
#import profile
import layout
import layout/shifts
import gaugeUtils
import staples

type
  GaugeActionCoeffs* = object
    plaq*: float
    rect*: float
    pgm*: float
    adjplaq*: float

# plaq: 6 types
# rect: 12 types
# pgm: 32=4*2*4=4*3*2+4*2 types
# shift corners: u[mu],nu mu != nu (12)
# make staples: s[mu][nu] mu != nu (12)
# plaq traces:
#  plaq: U[mu]^+ * sum_{nu!=mu} s[mu][nu] (6)
#  rect: shift(s[mu][nu], nu) (12 shifts)
#  pgm: shift(s[mu][nu], sig) (24 shifts)
proc gaugeAction1*[T](uu: openarray[T]): auto =
  mixin mul, redot, load1
  tic()
  let u = cast[ptr cArray[T]](unsafeAddr(uu[0]))
  let lo = u[0].l
  let nd = lo.nDim
  let np = (nd*(nd-1)) div 2
  let nc = u[0][0].ncols
  var cs = startCornerShifts(uu)
  toc("gaugeAction startCornerShifts")
  var (stf,stu,ss) = makeStaples(uu, cs)
  toc("gaugeAction makeStaples")
  #var ss = startStapleShifts(st)
  #toc("gaugeAction startStapleShifts")
  let maxThreads = getMaxThreads()
  var nth = 0
  var act = newSeq[float](3*maxThreads)
  toc("gaugeAction setup")
  threads:
    tic()
    var plaq = 0.0
    var rect = 0.0
    var pgm = 0.0
    for ir in u[0]:
      for mu in 1..<nd:
        for nu in 0..<mu:
          # plaq
          let p1 = redot(u[mu][ir], stf[mu,nu][ir])
          plaq += simdSum(p1)
          if isLocal(ss[mu][nu],ir) and isLocal(ss[nu][mu],ir):
            var bmu,bnu: type(load1(u[0][0]))
            localSB(ss[mu][nu], ir, assign(bmu,it), stu[mu,nu][ix])
            localSB(ss[nu][mu], ir, assign(bnu,it), stu[nu,mu][ix])
            # rect
            let r1 = redot(bmu, stf[mu,nu][ir])
            rect += simdSum(r1)
            let r2 = redot(bnu, stf[nu,mu][ir])
            rect += simdSum(r2)
    toc("gaugeAction local")
    for mu in 1..<nd:
      for nu in 0..<mu:
        var needBoundary = false
        boundaryWaitSB(ss[mu][nu]): needBoundary = true
        boundaryWaitSB(ss[nu][mu]): needBoundary = true
        if needBoundary:
          for ir in lo:
            if not isLocal(ss[mu][nu],ir) or not isLocal(ss[nu][mu],ir):
              var bmu,bnu: type(load1(u[0][0]))
              getSB(ss[mu][nu], ir, assign(bmu,it), stu[mu,nu][ix])
              getSB(ss[nu][mu], ir, assign(bnu,it), stu[nu,mu][ix])
              # rect
              let r1 = redot(bmu, stf[mu,nu][ir])
              rect += simdSum(r1)
              let r2 = redot(bnu, stf[nu,mu][ir])
              rect += simdSum(r2)
    act[threadNum*3]   = plaq
    act[threadNum*3+1] = rect
    act[threadNum*3+2] = pgm
    if threadNum==0: nth = numThreads
    toc("gaugeAction boundary")
  toc("gaugeAction threads")
  var a = [0.0, 0.0, 0.0]
  for i in 0..<nth:
    a[0] += act[i*3]
    a[1] += act[i*3+1]
    a[2] += act[i*3+2]
  #for i in 0..<3:
  #  a[i] = a[i]/(lo.physVol.float*float(np*nc))
  rankSum(a)
  #echo "plaq: ", a[0]
  #echo "rect: ", a[1]
  #echo "pgm: ", a[2]
  result = (-1.0/nc.float) * a[0]
  toc("gaugeAction end")

proc gaugeForce*[T](uu: openArray[T]): auto =
  mixin load1, adj
  tic()
  let u = cast[ptr cArray[T]](unsafeAddr(uu[0]))
  let lo = u[0].l
  let nd = lo.nDim
  let np = (nd*(nd-1)) div 2
  let nc = u[0][0].ncols
  var cs = startCornerShifts(uu)
  toc("gaugeForce startCornerShifts")
  var (stf,stu,ss) = makeStaples(uu, cs)
  toc("gaugeForce makeStaples")
  var f = newSeq[type(stf[0])](nd)
  for i in 0..<f.len:
    f[i].new(lo)
    f[i] := 0
  threads:
    tic()
    for ir in u[0]:
      for mu in 1..<nd:
        for nu in 0..<mu:
          # plaq
          f[mu][ir] += stf[mu,nu][ir]
          f[nu][ir] += stf[nu,mu][ir]
          if isLocal(ss[mu][nu],ir):
            var bmu: type(load1(u[0][0]))
            localSB(ss[mu][nu], ir, assign(bmu,it), stu[mu,nu][ix])
            f[mu][ir] += bmu
          if isLocal(ss[nu][mu],ir):
            var bnu: type(load1(u[0][0]))
            localSB(ss[nu][mu], ir, assign(bnu,it), stu[nu,mu][ix])
            f[nu][ir] += bnu
    toc("gaugeForce local")
    for mu in 1..<nd:
      for nu in 0..<mu:
        var needBoundary = false
        boundaryWaitSB(ss[mu][nu]): needBoundary = true
        boundaryWaitSB(ss[nu][mu]): needBoundary = true
        if needBoundary:
          boundarySyncSB()
          for ir in lo:
            if not isLocal(ss[mu][nu],ir):
              var bmu: type(load1(u[0][0]))
              getSB(ss[mu][nu], ir, assign(bmu,it), stu[mu,nu][ix])
              f[mu][ir] += bmu
            if not isLocal(ss[nu][mu],ir):
              var bnu: type(load1(u[0][0]))
              localSB(ss[nu][mu], ir, assign(bnu,it), stu[nu,mu][ix])
              f[nu][ir] += bnu
    toc("gaugeForce boundary")
  toc("gaugeForce threads")
  #toc("gaugeAction end")
  for mu in 0..<f.len:
    for e in f[mu]:
      mixin trace
      let s = u[mu][e]*f[mu][e].adj
      f[mu][e].projectTAH s
  toc("gaugeForce end")
  return f

proc gaugeAction2*(c: GaugeActionCoeffs, g: array|seq): auto =
  mixin redot
  tic()
  const nc = g[0][0].nrows
  let lo = g[0].l
  let nd = lo.nDim
  let t = newTransporters(g, g[0], 1)
  let t2 = newTransporters(g, g[0], 1)
  let td = newTransporters(g, g[0], -1)
  var pl = 0.0
  var rt = 0.0
  var pg = 0.0
  toc("gaugeAction2 setup")
  threads:
    tic()
    toc("gaugeAction2 zero")
    var ip = 0
    for mu in 1..<nd:
      for nu in 0..<mu:
        tic()
        var tpl = redot(t[mu]^*g[nu], t[nu]^*g[mu])
        if threadNum==0:
          pl += tpl
        #echo mu, " ", nu, " ", trace(m)/nc
        toc("gaugeAction2 pl")
        if c.rect != 0:
          var tr1 = redot(t[mu]^*t[nu]^*g[nu], t2[nu]^*t[nu]^*g[mu])
          var tr2 = redot(t2[mu]^*t[mu]^*g[nu], t[nu]^*t[mu]^*g[mu])
          if threadNum==0:
            rt += tr1 + tr2
          toc("gaugeAction2 rt")
        if c.pgm != 0:
          for sg in 0..<nu:
            var ts1 = redot(t[mu]^*t[nu]^*g[sg], t[sg]^*t[nu]^*g[mu])
            var ts2 = redot(t[mu]^*t[sg]^*g[nu], t[nu]^*t[sg]^*g[mu])
            var ts3 = redot(t[nu]^*t[mu]^*g[sg], t[sg]^*t[mu]^*g[nu])
            var ts4 = redot(t[nu]^*t[sg]^*g[mu], t[mu]^*t[sg]^*g[nu])
            var ts5 = redot(t[sg]^*t[mu]^*g[nu], t[nu]^*t[mu]^*g[sg])
            var ts6 = redot(t[sg]^*t[nu]^*g[mu], t[mu]^*t[nu]^*g[sg])
            #var ts7 = redot(td[sg]^*t[mu]^*td[nu]^*g[sg], td[nu]^*g[mu])
            #var ts8 = redot(td[sg]^*t[nu]^*td[mu]^*g[sg], td[mu]^*g[nu])
            var ts7 = redot(t[mu]^*td[nu]^*g[sg], t[sg]^*td[nu]^*g[mu])
            var ts8 = redot(t[mu]^*td[sg]^*g[nu], t[nu]^*td[sg]^*g[mu])
            if threadNum==0:
              pg += ts1 + ts2 + ts3 + ts4 + ts5 + ts6 + ts7 + ts8
          toc("gaugeAction2 pg")
    toc("gaugeAction2 work")
  toc("gaugeAction2 threads")
  #echo "plaq: ", pl, "  rect: ", rt, "  pgm: ", pg
  #result = (pl,rt,pg)
  result = (-1.0/nc.float) * (c.plaq*pl + c.rect*rt + c.pgm*pg)
template gaugeAction2*(g: array|seq, c: GaugeActionCoeffs): untyped =
  gaugeAction2(c, g)
proc gaugeAction2*(g: array|seq): auto =
  var c = GaugeActionCoeffs(plaq:1.0)
  gaugeAction2(c, g)

proc gaugeForce2*(f,g: array|seq, c: GaugeActionCoeffs) =
  mixin adj,projectTAH
  tic()
  let lo = g[0].l
  let nd = lo.nDim
  const nc = g[0][0].nrows
  let t = newTransporters(g, g[0], 1)
  let td = newTransporters(g, g[0], -1)
  toc("gaugeForce2 setup")
  threads:
    for mu in 0..<nd:
      #let mu = (mux + 1) mod nd
      f[mu] := 0
      for nu in 0..<nd:
        if nu==mu: continue
        discard t[nu] ^* g[mu]
        #shiftExpr(t[mu].sb, f[mu][ir] += t[nu].field[ir] * it, g[nu][ix].adj)
        shiftExpr(t[mu].sb, f[mu][ir] += t[nu].field[ir]*adj(it), g[nu][ix])
        f[mu] += td[nu] ^* t[mu] ^* g[nu]
    for mu in 0..<nd:
      for e in f[mu]:
        mixin trace
        let s = (c.plaq/nc.float) * g[mu][e] * f[mu][e].adj
        f[mu][e].projectTAH s
proc gaugeForce2*(f,g: array|seq): auto =
  var c = GaugeActionCoeffs(plaq:1.0)
  gaugeForce2(f,g,c)

proc actionA*(c: GaugeActionCoeffs, g: auto): auto =
  ## Specialized gauge action for plaq + adjplaq
  mixin mul, load1, createShiftBufs, re
  tic()
  let lo = g[0].l
  let nd = lo.nDim
  let nc = g[0][0].ncols
  var sf = newSeq[type(createShiftBufs(g[0],1,"all"))](nd)
  for i in 0..<nd-1:
    sf[i] = createShiftBufs(g[0], 1, "all")
  sf[nd-1].newSeq(nd)
  for i in 0..<nd-1: sf[nd-1][i] = sf[i][i]
  var pl = [0.0, 0.0]
  toc("plaq setup")
  threads:
    tic()
    var plt = [0.0, 0.0]
    var umunu,unumu: type(load1(g[0][0]))
    for mu in 0..<nd:
      for nu in 0..<nd:
        if mu != nu:
          startSB(sf[mu][nu], g[mu][ix])
    toc("plaq start shifts")
    for ir in g[0]:
      for mu in 1..<nd:
        for nu in 0..<mu:
          if isLocal(sf[mu][nu],ir) and isLocal(sf[nu][mu],ir):
            localSB(sf[mu][nu], ir, mul(unumu,g[nu][ir],it), g[mu][ix])
            localSB(sf[nu][mu], ir, mul(umunu,g[mu][ir],it), g[nu][ix])
            let dt = dot(umunu,unumu)
            plt[0] += simdSum(dt.re)
            plt[1] += simdSum(dt.norm2)
    toc("plaq local")
    var needBoundary = false
    for mu in 0..<nd:
      for nu in 0..<nd:
        if mu != nu:
          boundaryWaitSB(sf[mu][nu]): needBoundary = true
    toc("plaq wait")
    if needBoundary:
      boundarySyncSB()
      for ir in g[0]:
        for mu in 1..<nd:
          for nu in 0..<mu:
            if not isLocal(sf[mu][nu],ir) or not isLocal(sf[nu][mu],ir):
              if isLocal(sf[mu][nu], ir):
                localSB(sf[mu][nu], ir, mul(unumu,g[nu][ir],it), g[mu][ix])
              else:
                boundaryGetSB(sf[mu][nu], ir):
                  mul(unumu, g[nu][ir], it)
              if isLocal(sf[nu][mu], ir):
                localSB(sf[nu][mu], ir, mul(umunu,g[mu][ir],it), g[nu][ix])
              else:
                boundaryGetSB(sf[nu][mu], ir):
                  mul(umunu, g[mu][ir], it)
              let dt = dot(umunu,unumu)
              plt[0] += simdSum(dt.re)
              plt[1] += simdSum(dt.norm2)
    toc("plaq boundary")
    threadSum(plt)
    if threadNum == 0:
      pl[0] = plt[0] / float(nc)
      pl[1] = plt[1] / float(nc*nc)
      rankSum(pl)
    toc("plaq sum")
  let a0 = 0.5 * float(nd*(nd-1)*lo.physVol)
  result = c.plaq*(a0-pl[0]) + c.adjplaq*(a0-pl[1])
  toc("plaq end", flops=lo.nSites.float*float(2*8*nc*nc*nc-1))

proc forceA*(c: GaugeActionCoeffs, g,f: auto) =
  ## Specialized gauge force for plaq + adjplaq
  mixin load1, adj
  tic()
  let lo = g[0].l
  let nd = lo.nDim
  let nc = g[0][0].ncols
  let cp = c.plaq / float(nc)
  let ca = 2.0 * c.adjplaq / float(nc*nc)
  var cs = startCornerShifts(g)
  toc("gaugeForce startCornerShifts")
  var (stf,stu,ss) = makeStaples(g, cs)
  toc("gaugeForce makeStaples")
  for i in 0..<nd:
    f[i] := 0
  threads:
    tic()
    for ir in g[0]:
      for mu in 1..<nd:
        for nu in 0..<mu:
          let tmn = dot(stf[mu,nu][ir], g[mu][ir])
          f[mu][ir] += (cp+ca*tmn) * stf[mu,nu][ir]
          let tnm = dot(stf[nu,mu][ir], g[nu][ir])
          f[nu][ir] += (cp+ca*tnm) * stf[nu,mu][ir]
          if isLocal(ss[mu][nu],ir):
            var bmu: type(load1(g[0][0]))
            localSB(ss[mu][nu], ir, assign(bmu,it), stu[mu,nu][ix])
            let tmu = dot(bmu, g[mu][ir])
            f[mu][ir] += (cp+ca*tmu) * bmu
          if isLocal(ss[nu][mu],ir):
            var bnu: type(load1(g[0][0]))
            localSB(ss[nu][mu], ir, assign(bnu,it), stu[nu,mu][ix])
            let tnu = dot(bnu, g[nu][ir])
            f[nu][ir] += (cp+ca*tnu) * bnu
    toc("gaugeForce local")
    for mu in 1..<nd:
      for nu in 0..<mu:
        var needBoundary = false
        boundaryWaitSB(ss[mu][nu]): needBoundary = true
        boundaryWaitSB(ss[nu][mu]): needBoundary = true
        if needBoundary:
          boundarySyncSB()
          for ir in lo:
            if not isLocal(ss[mu][nu],ir):
              var bmu: type(load1(g[0][0]))
              getSB(ss[mu][nu], ir, assign(bmu,it), stu[mu,nu][ix])
              let tmu = dot(bmu, g[mu][ir])
              f[mu][ir] += (cp+ca*tmu) * bmu
            if not isLocal(ss[nu][mu],ir):
              var bnu: type(load1(g[0][0]))
              getSB(ss[nu][mu], ir, assign(bnu,it), stu[nu,mu][ix])
              let tnu = dot(bnu, g[nu][ir])
              f[nu][ir] += (cp+ca*tnu) * bnu
    #toc("gaugeForce boundary")
  toc("gaugeForce threads")
  #toc("gaugeAction end")
  for mu in 0..<f.len:
    for e in f[mu]:
      mixin trace
      let s = g[mu][e]*f[mu][e].adj
      f[mu][e].projectTAH s
  toc("gaugeForce end")

when isMainModule:
  import qex
  import physics/qcdTypes
  #import matrixFunctions
  qexInit()
  var defaultGaugeFile = "l88.scidac"
  #let defaultLat = @[2,2,2,2]
  let defaultLat = @[8,8,8,8]
  #let defaultLat = @[8,8,8]
  #let defaultLat = @[8,8]
  defaultSetup()
  #for mu in 0..<g.len: g[mu] := 1
  g.random

  proc test(g:auto) =
    var pl = plaq(g)
    echo "plaq:"
    echo pl
    echo pl.sum
    var gc = GaugeActionCoeffs(plaq:1.0)
    #var ga = gaugeAction(g)
    var ga = gaugeAction.gaugeAction1(g)
    var ga2 = gaugeAction.gaugeAction2(gc,g)
    echo "ga: ", ga, "\t", ga2
    var f = gaugeAction.gaugeForce(g)
    var f2 = g[0].l.newGauge
    gaugeAction.gaugeForce2(f2,g)
    for i in 0..<f.len:
      echo "f[", i, "]: ", f[i].norm2, "\t", f2[i].norm2

  test(g)
  echoTimers()
  resetTimers()
  test(g)

  proc updateX(g,p,eps:auto) =
    mixin exp
    for mu in 0..<g.len:
      for e in g[mu]:
        let t = exp(eps*p[mu][e])*g[mu][e]
        g[mu][e] := t
      #echo "g[", mu, "]: ", g[mu].norm2

  proc updateP(g,p,eps:auto) =
    #let f = gaugeForce(g)
    var f = g[0].l.newGauge
    gaugeAction.gaugeForce2(f, g)
    #gaugeForce2(f, g)
    #gaugeForce2(f, g)
    for mu in 0..<f.len:
      #echo "f[", mu, "]: ", f[mu].norm2
      p[mu] += (-eps)*f[mu]

  var g0 = g[0].l.newGauge
  for mu in 0..<g.len: g0[mu] := g[mu]
  proc test2(steps:int) =
    const t = 0.02
    let eps = t/steps.float
    var p = newSeq[type(g[0])](g.len)
    for mu in 0..<p.len:
      g[mu] := g0[mu]
      p[mu].new(g[0].l)
      for e in p[mu]:
        when p[mu][e].nrows==1:
          #p[mu][e] := asImag(1)
          p[mu][e] := newComplex(0,1)
        else:
          p[mu][e] := 0
          let t = (2*(e mod 2)-1).float
          #let t = 1.0
          p[mu][e][0,1] := t
          p[mu][e][1,0] := -t
    var gc = GaugeActionCoeffs(plaq:1.0)
    let ga = gaugeAction.gaugeAction2(gc,g)
    var p2 = 0.0
    for mu in 0..<p.len: p2 += p[mu].norm2
    let s0 = ga + 0.5*p2
    echo "ACT: ", ga, "\t", 0.5*p2, "\t", s0

    for n in 1..steps:
      #echo "pdiff: ", (p[0]-p[1]).norm2
      #echo "gdiff: ", (g[0]-g[1]).norm2
      #echo "ga: ", gaugeAction.gaugeAction2(gc,g)
      updateX(g,p,0.5*eps)
      #echo "pdiff: ", (p[0]-p[1]).norm2
      #echo "gdiff: ", (g[0]-g[1]).norm2
      #echo "ga: ", gaugeAction.gaugeAction2(gc,g)
      updateP(g,p,eps)
      #echo "pdiff: ", (p[0]-p[1]).norm2
      #echo "gdiff: ", (g[0]-g[1]).norm2
      #echo "ga: ", gaugeAction.gaugeAction2(gc,g)
      updateX(g,p,0.5*eps)

    let ga2 = gaugeAction.gaugeAction2(gc,g)
    p2 = 0.0
    for mu in 0..<p.len: p2 += p[mu].norm2
    let s2 = ga2 + 0.5*p2
    echo "ACT2: ", ga2, "\t", 0.5*p2, "\t", s2
    echo "dH: ", s2 - s0
    let sr = (s2-s0)/(eps*eps)
    echo "error rate: ", sr

  #test2(2000)
  test2(200)
  test2(20)
  test2(2)

  qexFinalize()
