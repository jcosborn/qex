import qex
import stdUtils
import profile
import gaugeUtils
import staples

type
  GaugeActionCoeffs* = object
    plaq*: float
    rect*: float
    pgm*: float

# plaq: 6 types
# rect: 12 types
# pgm: 32=4*2*4=4*3*2+4*2 types
# shift corners: u[mu],nu mu != nu (12)
# make staples: s[mu][nu] mu != nu (12)
# plaq traces:
#  plaq: U[mu]^+ * sum_{nu!=mu} s[mu][nu] (6)
#  rect: shift(s[mu][nu], nu) (12 shifts)
#  pgm: shift(s[mu][nu], sig) (24 shifts)
proc gaugeAction*[T](uu: openArray[T]): auto =
  mixin mul
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
          let p1 = redot(u[mu][ir], stf[mu][nu][ir])
          plaq += simdSum(p1)
          if isLocal(ss[mu][nu],ir) and isLocal(ss[nu][mu],ir):
            var bmu,bnu: type(load1(u[0][0]))
            localSB(ss[mu][nu], ir, assign(bmu,it), stu[mu][nu][ix])
            localSB(ss[nu][mu], ir, assign(bnu,it), stu[nu][mu][ix])
            # rect
            let r1 = redot(bmu, stf[mu][nu][ir])
            rect += simdSum(r1)
            let r2 = redot(bnu, stf[nu][mu][ir])
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
              getSB(ss[mu][nu], ir, assign(bmu,it), stu[mu][nu][ix])
              getSB(ss[nu][mu], ir, assign(bnu,it), stu[nu][mu][ix])
              # rect
              let r1 = redot(bmu, stf[mu][nu][ir])
              rect += simdSum(r1)
              let r2 = redot(bnu, stf[nu][mu][ir])
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
  result = a[0]
  toc("gaugeAction end")

proc gaugeForce*[T](uu: openArray[T]): auto =
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
  var f:type(stf[0])
  f.newSeq(nd)
  for i in 0..<f.len:
    f[i].new(lo)
    f[i] := 0
  threads:
    tic()
    for ir in u[0]:
      for mu in 1..<nd:
        for nu in 0..<mu:
          # plaq
          f[mu][ir] += stf[mu][nu][ir]
          f[nu][ir] += stf[nu][mu][ir]
          if isLocal(ss[mu][nu],ir):
            var bmu: type(load1(u[0][0]))
            localSB(ss[mu][nu], ir, assign(bmu,it), stu[mu][nu][ix])
            f[mu][ir] += bmu
          if isLocal(ss[nu][mu],ir):
            var bnu: type(load1(u[0][0]))
            localSB(ss[nu][mu], ir, assign(bnu,it), stu[nu][mu][ix])
            f[nu][ir] += bnu
    toc("gaugeForce local")
    for mu in 1..<nd:
      for nu in 0..<mu:
        var needBoundary = false
        boundaryWaitSB(ss[mu][nu]): needBoundary = true
        boundaryWaitSB(ss[nu][mu]): needBoundary = true
        if needBoundary:
          for ir in lo:
            if not isLocal(ss[mu][nu],ir):
              var bmu: type(load1(u[0][0]))
              getSB(ss[mu][nu], ir, assign(bmu,it), stu[mu][nu][ix])
              f[mu][ir] += bmu
            if not isLocal(ss[nu][mu],ir):
              var bnu: type(load1(u[0][0]))
              localSB(ss[nu][mu], ir, assign(bnu,it), stu[nu][mu][ix])
              f[nu][ir] += bnu
    toc("gaugeForce boundary")
  toc("gaugeForce threads")
  #toc("gaugeAction end")
  for mu in 0..<f.len:
    for e in f[mu]:
      mixin trace
      let s = u[mu][e]*f[mu][e].adj
      let t = -0.5*(s-s.adj)
      f[mu][e] := t - (trace(t)*(1.0/t.nrows.float))
  toc("gaugeForce end")
  return f

proc gaugeAction2*(c: GaugeActionCoeffs, g: array|seq): auto =
  tic()
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
        var tr1 = redot(t[mu]^*t[nu]^*g[nu], t2[nu]^*t[nu]^*g[mu])
        var tr2 = redot(t2[mu]^*t[mu]^*g[nu], t[nu]^*t[mu]^*g[mu])
        if threadNum==0:
          rt += tr1 + tr2
        toc("gaugeAction2 rt")
        for sg in 0..<nu:
          var ts1 = redot(t[mu]^*t[nu]^*g[sg], t[sg]^*t[nu]^*g[mu])
          var ts2 = redot(t[mu]^*t[sg]^*g[nu], t[nu]^*t[sg]^*g[mu])
          var ts3 = redot(t[nu]^*t[mu]^*g[sg], t[sg]^*t[mu]^*g[nu])
          var ts4 = redot(t[nu]^*t[sg]^*g[mu], t[mu]^*t[sg]^*g[nu])
          var ts5 = redot(t[sg]^*t[mu]^*g[nu], t[nu]^*t[mu]^*g[sg])
          var ts6 = redot(t[sg]^*t[nu]^*g[mu], t[mu]^*t[nu]^*g[sg])
          var ts7 = redot(td[sg]^*t[mu]^*td[nu]^*g[sg], td[nu]^*g[mu])
          var ts8 = redot(td[sg]^*t[nu]^*td[mu]^*g[sg], td[mu]^*g[nu])
          if threadNum==0:
            pg += ts1 + ts2 + ts3 + ts4 + ts5 + ts6 + ts7 + ts8
        toc("gaugeAction2 pg")
    toc("gaugeAction2 work")
  toc("gaugeAction2 threads")
  echo "plaq: ", pl, "  rect: ", rt, "  pgm: ", pg
  #result = (pl,rt,pg)
  result = c.plaq*pl + c.rect*rt + c.pgm*pg
proc gaugeAction2*(g: array|seq): auto =
  var c: GaugeActionCoeffs
  gaugeAction2(c, g)

proc gaugeForce2*(f,g: array|seq) =
  tic()
  let lo = g[0].l
  let nd = lo.nDim
  let nc = g[0][0].nrows
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
        shiftExpr(t[mu].sb, f[mu][ir] += t[nu].field[ir] * it, g[nu][ix].adj)
        f[mu] += td[nu] ^* t[mu] ^* g[nu]
    for mu in 0..<nd:
      for e in f[mu]:
        mixin trace
        let s = g[mu][e] * f[mu][e].adj
        let t = 0.5*(s-s.adj)
        f[mu][e] := t - (trace(t)*(1.0/nc.float))

when isMainModule:
  import qcdTypes
  import matrixFunctions
  qexInit()
  var defaultGaugeFile = "l88.scidac"
  #let defaultLat = @[2,2,2,2]
  let defaultLat = @[8,8,8,8]
  defaultSetup()
  for mu in 0..<g.len: g[mu] := 1
  #g.random

  proc test(g:any) =
    var pl = plaq(g)
    echo "plaq:"
    echo pl
    echo pl.sum
    var gc: GaugeActionCoeffs
    var ga = gaugeAction(g)
    var ga2 = gc.gaugeAction2(g)
    echo "ga: ", ga, "\t", ga2
    var f = gaugeForce(g)
    var f2 = g[0].l.newGauge
    gaugeForce2(f2,g)
    for i in 0..<f.len:
      echo "f[", i, "]: ", f[i].norm2, "\t", f2[i].norm2

  test(g)
  echoTimers()
  resetTimers()
  test(g)

  proc updateX(g,p,eps:any) =
    for mu in 0..<g.len:
      for e in g[mu]:
        let t = exp(eps*p[mu][e])*g[mu][e]
        g[mu][e] := t
      echo "g[", mu, "]: ", g[mu].norm2

  proc updateP(g,p,eps:any) =
    #let f = gaugeForce(g)
    var f = g[0].l.newGauge
    gaugeForce2(f, g)
    #gaugeForce2(f, g)
    #gaugeForce2(f, g)
    for mu in 0..<f.len:
      echo "f[", mu, "]: ", f[mu].norm2
      p[mu] += eps*f[mu]

  proc test2(g,eps:any) =
    var p = newSeq[type(g[0])](g.len)
    for mu in 0..<p.len:
      g[mu] := 1
      p[mu].new(g[0].l)
      for e in p[mu]:
        p[mu][e] := 0
        let t = (2*(e mod 2)-1).float
        #let t = 1.0
        p[mu][e][0,1] := t
        p[mu][e][1,0] := -t
    var gc: GaugeActionCoeffs
    let ga = gc.gaugeAction2(g)
    var p2 = 0.0
    for mu in 0..<p.len: p2 += p[mu].norm2
    let s0 = ga + 0.5*p2
    echo "ACT: ", ga, "\t", 0.5*p2, "\t", s0

    echo "pdiff: ", (p[0]-p[1]).norm2
    echo "gdiff: ", (g[0]-g[1]).norm2
    echo "ga: ", gc.gaugeAction2(g)
    updateX(g,p,0.5*eps)
    echo "pdiff: ", (p[0]-p[1]).norm2
    echo "gdiff: ", (g[0]-g[1]).norm2
    echo "ga: ", gc.gaugeAction2(g)
    updateP(g,p,eps)
    echo "pdiff: ", (p[0]-p[1]).norm2
    echo "gdiff: ", (g[0]-g[1]).norm2
    echo "ga: ", gc.gaugeAction2(g)
    updateX(g,p,0.5*eps)

    let ga2 = gc.gaugeAction2(g)
    p2 = 0.0
    for mu in 0..<p.len: p2 += p[mu].norm2
    let s2 = ga2 + 0.5*p2
    echo ga2, "\t", 0.5*p2, "\t", s2
    echo s2 - s0
    let sr = (s2-s0)/(eps*eps)
    echo sr

  test2(g, 1e-5)
  test2(g, 1e-4)
  test2(g, 1e-3)
  test2(g, 1e-2)

  qexFinalize()
