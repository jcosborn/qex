import base
#import stdUtils
#import profile
import maths
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

proc DBW2*(beta:float, c1:float = -1.4088):GaugeActionCoeffs =
  result.plaq = (1.0-8.0*c1)*beta
  result.rect = c1*beta

# plaq: 6 types
# rect: 12 types
# pgm: 32=4*2*4=4*3*2+4*2 types
# shift corners: u[mu],nu mu != nu (12)
# make staples: s[mu][nu] mu != nu (12)
# plaq traces:
#  plaq: U[mu]^+ * sum_{nu!=mu} s[mu][nu] (6)
#  rect: shift(s[mu][nu], nu) (12 shifts)
#  pgm: shift(s[mu][nu], sig) (24 shifts)
proc gaugeAction1*[T](c: GaugeActionCoeffs, uu: openarray[T]): auto =
  mixin mul, redot, load1
  tic("gaugeAction1")
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
          if c.rect!=0:
            if isLocal(ss[mu][nu],ir):
              var bmu: type(load1(u[0][0]))
              localSB(ss[mu][nu], ir, assign(bmu,it), stu[mu,nu][ix])
              # rect
              let r = redot(bmu, stf[mu,nu][ir])
              rect += simdSum(r)
            if isLocal(ss[nu][mu],ir):
              var bnu: type(load1(u[0][0]))
              localSB(ss[nu][mu], ir, assign(bnu,it), stu[nu,mu][ix])
              # rect
              let r = redot(bnu, stf[nu,mu][ir])
              rect += simdSum(r)
    toc("gaugeAction local")
    if c.rect!=0:
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
                # rect
                let r = redot(bmu, stf[mu,nu][ir])
                rect += simdSum(r)
              if not isLocal(ss[nu][mu],ir):
                var bnu: type(load1(u[0][0]))
                getSB(ss[nu][mu], ir, assign(bnu,it), stu[nu,mu][ix])
                # rect
                let r = redot(bnu, stf[nu,mu][ir])
                rect += simdSum(r)
    act[threadNum*3]   = plaq
    act[threadNum*3+1] = rect
    act[threadNum*3+2] = pgm
    if threadNum==0: nth = numThreads
    # toc("gaugeAction boundary")
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
  result = (-1.0/nc.float) * (c.plaq*a[0] + c.rect*a[1] + c.pgm*a[2])
  toc("gaugeAction end")

proc gaugeAction1*[T](uu: openarray[T]): auto =
  let gc = GaugeActionCoeffs(plaq:1.0)
  return gc.gaugeAction1(uu)

proc gaugeForce*[T](c: GaugeActionCoeffs, uu: openArray[T], f: array|seq) =
  mixin load1, adj
  tic("gaugeForce")
  let u = cast[ptr cArray[T]](unsafeAddr(uu[0]))
  let lo = u[0].l
  let nd = lo.nDim
  let np = (nd*(nd-1)) div 2
  let nc = u[0][0].ncols
  let cp = c.plaq / float(nc)
  let cr = c.rect / float(nc)
  var cs = startCornerShifts(uu)
  var ru:FieldArray[type(u[0]).V,type(u[0]).T]  # the rect parts of 3
  var sb:seq[seq[ShiftB[type(u[0][0])]]]  # backward ru
  var sf:seq[seq[ShiftB[type(u[0][0])]]]  # forward stf
  if cr!=0:
    ru = newFieldArray2(lo,type(u[0]),[nd,nd],mu!=nu)
    sb.newseq(nd)
    for mu in 0..<nd:
      sb[mu].newseq(nd)
      for nu in 0..<nd:
        if mu==nu: continue
        sb[mu][nu].initShiftB(ru[mu,nu], nu, -1, "all")
    sf.newseq(nd)
    for mu in 0..<nd:
      sf[mu].newseq(nd)
      for nu in 0..<nd:
        if mu==nu: continue
        sf[mu][nu].initShiftB(u[mu], nu, 1, "all")
  toc("gaugeForce init")
  var (stf,stu,ss) = makeStaples(uu, cs)
  toc("gaugeForce makeStaples")
  threads:
    tic()
    if cr!=0:
      for mu in 1..<nd:
        for nu in 0..<mu:
          sf[mu][nu].startSB(stf[mu,nu][ix])
          sf[nu][mu].startSB(stf[nu,mu][ix])
    for mu in 0..<nd:
      f[mu] := 0
      if cr!=0:
        for nu in 0..<nd:
          if mu!=nu:
            ru[mu,nu] := 0
    for ir in u[0]:
      for mu in 1..<nd:
        for nu in 0..<mu:
          # plaq
          f[mu][ir] += cp * stf[mu,nu][ir]
          f[nu][ir] += cp * stf[nu,mu][ir]
          if isLocal(ss[mu][nu],ir):
            var bmu: type(load1(u[0][0]))
            localSB(ss[mu][nu], ir, assign(bmu,it), stu[mu,nu][ix])
            f[mu][ir] += cp * bmu
            if cr!=0:
              var umu,unu,bmunu: type(load1(u[0][0]))
              getSB(cs[nu][mu], ir, assign(unu,it), u[nu][ix])
              getSB(cs[mu][nu], ir, assign(umu,it), u[mu][ix])
              bmunu := bmu * unu
              f[nu][ir] += cr * bmunu * umu.adj
              ru[nu,mu][ir] += bmu.adj * u[nu][ir] * umu
              ru[mu,nu][ir] += u[nu][ir].adj * bmunu
          if isLocal(ss[nu][mu],ir):
            var bnu: type(load1(u[0][0]))
            localSB(ss[nu][mu], ir, assign(bnu,it), stu[nu,mu][ix])
            f[nu][ir] += cp * bnu
            if cr!=0:
              var unu,umu,bnumu: type(load1(u[0][0]))
              getSB(cs[mu][nu], ir, assign(umu,it), u[mu][ix])
              getSB(cs[nu][mu], ir, assign(unu,it), u[nu][ix])
              bnumu := bnu * umu
              f[mu][ir] += cr * bnumu * unu.adj
              ru[mu,nu][ir] += bnu.adj * u[mu][ir] * unu
              ru[nu,mu][ir] += u[mu][ir].adj * bnumu
          if cr!=0:
            if isLocal(sf[mu][nu],ir):
              var smu,unu,smunu: type(load1(u[0][0]))
              localSB(sf[mu][nu], ir, assign(smu,it), stf[mu,nu][ix])
              getSB(cs[nu][mu], ir, assign(unu,it), u[nu][ix])
              smunu := smu * unu.adj
              f[mu][ir] += cr * u[nu][ir] * smunu
              f[nu][ir] += cr * u[mu][ir] * smunu.adj
              ru[nu,mu][ir] += u[mu][ir].adj * u[nu][ir] * smu
            if isLocal(sf[nu][mu],ir):
              var snu,umu,snumu: type(load1(u[0][0]))
              localSB(sf[nu][mu], ir, assign(snu,it), stf[nu,mu][ix])
              getSB(cs[mu][nu], ir, assign(umu,it), u[mu][ix])
              snumu := snu * umu.adj
              f[nu][ir] += cr * u[mu][ir] * snumu
              f[mu][ir] += cr * u[nu][ir] * snumu.adj
              ru[mu,nu][ir] += u[nu][ir].adj * u[mu][ir] * snu
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
              f[mu][ir] += cp * bmu
              if cr!=0:
                var umu,unu,bmunu: type(load1(u[0][0]))
                getSB(cs[nu][mu], ir, assign(unu,it), u[nu][ix])
                getSB(cs[mu][nu], ir, assign(umu,it), u[mu][ix])
                bmunu := bmu * unu
                f[nu][ir] += cr * bmunu * umu.adj
                ru[nu,mu][ir] += bmu.adj * u[nu][ir] * umu
                ru[mu,nu][ir] += u[nu][ir].adj * bmunu
            if not isLocal(ss[nu][mu],ir):
              var bnu: type(load1(u[0][0]))
              getSB(ss[nu][mu], ir, assign(bnu,it), stu[nu,mu][ix])
              f[nu][ir] += cp * bnu
              if cr!=0:
                var unu,umu,bnumu: type(load1(u[0][0]))
                getSB(cs[mu][nu], ir, assign(umu,it), u[mu][ix])
                getSB(cs[nu][mu], ir, assign(unu,it), u[nu][ix])
                bnumu := bnu * umu
                f[mu][ir] += cr * bnumu * unu.adj
                ru[mu,nu][ir] += bnu.adj * u[mu][ir] * unu
                ru[nu,mu][ir] += u[mu][ir].adj * bnumu
    if cr!=0:
      for mu in 1..<nd:
        for nu in 0..<mu:
          var needBoundary = false
          boundaryWaitSB(sf[mu][nu]): needBoundary = true
          boundaryWaitSB(sf[nu][mu]): needBoundary = true
          if needBoundary:
            boundarySyncSB()
            for ir in lo:
              if not isLocal(sf[mu][nu],ir):
                var smu,unu,smunu: type(load1(u[0][0]))
                getSB(sf[mu][nu], ir, assign(smu,it), stf[mu,nu][ix])
                getSB(cs[nu][mu], ir, assign(unu,it), u[nu][ix])
                smunu := smu * unu.adj
                f[mu][ir] += cr * u[nu][ir] * smunu
                f[nu][ir] += cr * u[mu][ir] * smunu.adj
                ru[nu,mu][ir] += u[mu][ir].adj * u[nu][ir] * smu
              if not isLocal(sf[nu][mu],ir):
                var snu,umu,snumu: type(load1(u[0][0]))
                getSB(sf[nu][mu], ir, assign(snu,it), stf[nu,mu][ix])
                getSB(cs[mu][nu], ir, assign(umu,it), u[mu][ix])
                snumu := snu * umu.adj
                f[nu][ir] += cr * u[mu][ir] * snumu
                f[mu][ir] += cr * u[nu][ir] * snumu.adj
                ru[mu,nu][ir] += u[nu][ir].adj * u[mu][ir] * snu
          threadBarrier()
          sb[mu][nu].startSB(ru[mu,nu][ix])
          sb[nu][mu].startSB(ru[nu,mu][ix])
      toc("gaugeForce staple boundary")
      for ir in u[0]:
        for mu in 1..<nd:
          for nu in 0..<mu:
            if isLocal(sb[mu][nu],ir):
              var b: type(load1(u[0][0]))
              localSB(sb[mu][nu], ir, assign(b,it), ru[mu,nu][ix])
              f[mu][ir] += cr * b
            if isLocal(sb[nu][mu],ir):
              var b: type(load1(u[0][0]))
              localSB(sb[nu][mu], ir, assign(b,it), ru[nu,mu][ix])
              f[nu][ir] += cr * b
      toc("gaugeForce back rect local")
      for mu in 1..<nd:
        for nu in 0..<mu:
          var needBoundary = false
          boundaryWaitSB(sb[mu][nu]): needBoundary = true
          boundaryWaitSB(sb[nu][mu]): needBoundary = true
          if needBoundary:
            boundarySyncSB()
            for ir in lo:
              if not isLocal(sb[mu][nu],ir):
                var b: type(load1(u[0][0]))
                getSB(sb[mu][nu], ir, assign(b,it), ru[mu,nu][ix])
                f[mu][ir] += cr * b
              if not isLocal(sb[nu][mu],ir):
                var b: type(load1(u[0][0]))
                getSB(sb[nu][mu], ir, assign(b,it), ru[nu,mu][ix])
                f[nu][ir] += cr * b
    for mu in 0..<f.len:
      for e in f[mu]:
        mixin trace
        let s = u[mu][e]*f[mu][e].adj
        f[mu][e].projectTAH s
  toc("gaugeForce end")

proc gaugeForce*[T](uu: openArray[T]): auto =
  let lo = uu[0].l
  var f = newOneOf @uu
  let gc = GaugeActionCoeffs(plaq:1.0)
  gc.gaugeForce(uu,f)
  return f

proc gaugeForce*(f,g: array|seq) =
  var c = GaugeActionCoeffs(plaq:1.0)
  gaugeForce(c,g,f)

proc gaugeAction2*(c: GaugeActionCoeffs, g: array|seq): auto =
  mixin redot
  tic("gaugeAction2")
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

proc gaugeForce2*(c: GaugeActionCoeffs, g,f: array|seq) =
  mixin adj,projectTAH
  tic("gaugeForce2")
  let lo = g[0].l
  let nd = lo.nDim
  const nc = g[0][0].nrows
  let cp = c.plaq / float(nc)
  let cr = c.rect / float(nc)
  let t = newTransporters(g, g[0], 1)
  let t2 = newTransporters(g, g[0], 1)
  let tg = newTransporters(g, g[0], 1)
  let td = newTransporters(g, g[0], -1)
  let td2 = newTransporters(g, g[0], -1)
  toc("gaugeForce2 setup")
  threads:
    for mu in 0..<nd:
      #let mu = (mux + 1) mod nd
      f[mu] := 0
      for nu in 0..<nd:
        if nu==mu: continue
        discard t[nu] ^* g[mu]
        shiftExpr(t[mu].sb, f[mu][ir] += cp * t[nu].field[ir]*adj(it), g[nu][ix])
        f[mu] += cp * td[nu] ^* t[mu] ^* g[nu]
        if cr != 0:
          discard t2[nu] ^* t[nu] ^* g[mu]
          discard tg[nu] ^* g[nu]
          shiftExpr(t[mu].sb, f[mu][ir] += cr * t2[nu].field[ir]*adj(it), tg[nu].field[ix])
          f[mu] += cr * td2[nu] ^* td[nu] ^* tg[mu] ^* t[nu] ^* g[nu]
          f[mu] += cr * td2[mu] ^* td[nu] ^* tg[mu] ^* t[mu] ^* g[nu]
          discard td[nu] ^* tg[mu] ^* t[mu] ^* g[nu]
          shiftExpr(t2[mu].sb, f[mu][ir] += cr * td[nu].field[ir]*adj(it), g[mu][ix])
          discard td[mu] ^* t[nu] ^* tg[mu] ^* g[mu]
          shiftExpr(t2[mu].sb, f[mu][ir] += cr * td[mu].field[ir]*adj(it), g[nu][ix])
          shiftExpr(t2[mu].sb, f[mu][ir] += cr * t[nu].field[ir]*adj(it), t[mu].field[ix])
    for mu in 0..<nd:
      for e in f[mu]:
        let s = g[mu][e] * f[mu][e].adj
        f[mu][e].projectTAH s
  toc("end")
proc gaugeForce2*(f,g: array|seq) =
  var c = GaugeActionCoeffs(plaq:1.0)
  gaugeForce2(c,g,f)

proc gaugeAction3*(c: GaugeActionCoeffs, g: array|seq): auto =
  tic("gaugeAction3")
  const nc = g[0][0].nrows
  let lo = g[0].l
  var pl = 0.0
  var rt = 0.0
  for mu in 1..<lo.nDim:
    for nu in 0..<mu:
      var ls = newseq[seq[int]]()
      block:
        let
          mu = mu+1
          nu = nu+1
        if c.plaq!=0:
          ls.add @[mu, nu, -mu, -nu]
        if c.rect!=0:
          ls.add [@[mu, nu, nu, -mu, -nu, -nu], @[mu, mu, nu, -mu, -mu, -nu]]
      let ws = g.wilsonLines ls
      var i = 0
      if c.plaq!=0:
        pl += ws[0].re
        i = 1
      if c.rect!=0:
        rt += ws[i].re + ws[i+1].re
  result = (-lo.physVol.float) * (c.plaq*pl + c.rect*rt)
  toc("end")
proc gaugeAction3*(g: array|seq): auto =
  var c = GaugeActionCoeffs(plaq:1.0)
  gaugeAction3(c, g)

proc plaqRectPath_fun(c:GaugeActionCoeffs, mu,nu:int):auto =
  let
    mu = mu+1
    nu = nu+1
  var ls = newseq[seq[int]]()
  if c.plaq!=0:
    ls.add [ @[nu, mu, -nu], @[-nu, mu, nu]
           , @[mu, nu, -mu], @[-mu, nu, mu]
           ]
  if c.rect!=0:
    ls.add [ @[nu, nu, mu, -nu, -nu], @[-nu, -nu, mu, nu, nu]
           , @[nu, mu, mu, -nu, -mu], @[-mu, nu, mu, mu, -nu]
           , @[-nu, mu, mu, nu, -mu], @[-mu, -nu, mu, mu, nu]
           , @[mu, mu, nu, -mu, -mu], @[-mu, -mu, nu, mu, mu]
           , @[mu, nu, nu, -mu, -nu], @[-nu, mu, nu, nu, -mu]
           , @[-mu, nu, nu, mu, -nu], @[-nu, -mu, nu, nu, mu]
           ]
  ls.optimalPairs

proc plaqRectPath(c:GaugeActionCoeffs, mu,nu:int):auto =
  var j = 0
  if c.plaq!=0:
    inc j
  if c.rect!=0:
    j += 2
  memoize(j,mu,nu):
    c.plaqRectPath_fun(mu,nu)

proc gaugeForce3*(c: GaugeActionCoeffs, g,f: auto) =
  tic("gaugeForce3")
  const nc = g[0][0].nrows
  let nd = g[0].l.nDim
  let cp = c.plaq / nc.float
  let cr = c.rect / nc.float
  threads:
    for mu in 0..<nd:
      f[mu] := 0
  for mu in 1..<nd:
    for nu in 0..<mu:
      let ptree = c.plaqRectPath(mu,nu)
      let ws = g.gaugeProd ptree
      threads:
        for ir in f[mu]:
          var pmu,rmu,pnu,rnu: type(load1(f[0][0]))
          var i = 0
          if c.plaq!=0:
            for j in 0..<2:
              pmu += ws[j][ir]
              pnu += ws[2+j][ir]
            i = 4
            f[mu][ir] += cp * pmu
            f[nu][ir] += cp * pnu
          if c.rect!=0:
            for j in 0..<6:
              rmu += ws[i+j][ir]
              rnu += ws[6+i+j][ir]
            f[mu][ir] += cr * rmu
            f[nu][ir] += cr * rnu
  threads:
    for mu in 0..<nd:
      for e in f[mu]:
        let s = g[mu][e] * f[mu][e].adj
        f[mu][e].projectTAH s
  toc("end")
proc gaugeForce3*(f,g: array|seq) =
  var c = GaugeActionCoeffs(plaq:1.0)
  gaugeForce3(c,g,f)

proc actionA*(c: GaugeActionCoeffs, g: auto): auto =
  ## Specialized gauge action for plaq + adjplaq
  mixin mul, load1, createShiftBufs, re
  tic("actionA")
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
  tic("forceA")
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
    tic("test")
    echo "Test C_plaq = 1"
    var pl = plaq(g)
    echo "plaq:"
    echo pl
    echo pl.sum
    var gc = GaugeActionCoeffs(plaq:1.0)
    #var ga = gaugeAction(g)
    toc("plaq")
    var ga = gaugeAction.gaugeAction1(g)
    toc("ga1")
    var ga2 = gaugeAction.gaugeAction2(gc,g)
    toc("ga2")
    var ga3 = gaugeAction.gaugeAction3(gc,g)
    toc("ga3")
    var gaA = gaugeAction.actionA(gc,g)
    toc("aA")
    echo "ga: ", ga, "\t", ga2, "\t", ga3, "\t", gaA
    var f = gaugeAction.gaugeForce(g)
    toc("gf")
    var f2 = newOneOf g
    var f3 = newOneOf g
    var fa = newOneOf g
    gaugeAction.gaugeForce2(f2,g)
    toc("gf2")
    gaugeAction.gaugeForce3(f3,g)
    toc("gf3")
    gaugeAction.forceA(gc,g,fa)
    toc("fA")
    for i in 0..<f.len:
      echo "f[", i, "]: ", f[i].norm2, "\t", f2[i].norm2, "\t", f3[i].norm2, "\t", fa[i].norm2
    toc("end")

  proc testR(g:auto) =
    tic("testR")
    echo "Test Rectangle"
    var gc = gaugeAction.DBW2(0.7796)
    echo gc
    var ga = gaugeAction.gaugeAction1(gc,g)
    toc("ga1")
    var ga2 = gaugeAction.gaugeAction2(gc,g)
    toc("ga2")
    var ga3 = gaugeAction.gaugeAction3(gc,g)
    toc("ga3")
    echo "ga: ", ga, "\t", ga2, "\t", ga3
    var f = newOneOf g
    var f2 = newOneOf g
    var f3 = newOneOf g
    toc("init f")
    gaugeAction.gaugeForce(gc,g,f)
    toc("gf")
    gaugeAction.gaugeForce2(gc,g,f2)
    toc("gf2")
    gaugeAction.gaugeForce3(gc,g,f3)
    toc("gf3")
    echo "gf: \t",  f[0].norm2, "\t",  f[1].norm2, "\t",  f[2].norm2, "\t",  f[3].norm2
    echo "gf2:\t", f2[0].norm2, "\t", f2[1].norm2, "\t", f2[2].norm2, "\t", f2[3].norm2
    echo "gf3:\t", f3[0].norm2, "\t", f3[1].norm2, "\t", f3[2].norm2, "\t", f3[3].norm2
    toc("end")

  test(g)
  testR(g)

  proc updateX(g,p,eps:auto) =
    mixin exp
    for mu in 0..<g.len:
      for e in g[mu]:
        let t = exp(eps*p[mu][e])*g[mu][e]
        g[mu][e] := t
      #echo "g[", mu, "]: ", g[mu].norm2

  proc updateP(c:GaugeActionCoeffs, g,p,eps:auto) =
    var f = newOneOf g
    gaugeAction.gaugeForce(c, g, f)
    for mu in 0..<f.len:
      #echo "f[", mu, "]: ", f[mu].norm2
      p[mu] += (-eps)*f[mu]

  var g0 = g[0].l.newGauge
  for mu in 0..<g.len: g0[mu] := g[mu]
  proc test2(steps:int, lambda=0.1931833):auto {.discardable.} =
    const t = 0.02
    let eps = t/steps.float
    echo "eps: ",eps
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
    var gc = gaugeAction.DBW2(0.7796)
    let ga = gaugeAction.gaugeAction1(gc,g)
    var p2 = 0.0
    for mu in 0..<p.len: p2 += p[mu].norm2
    let s0 = ga + 0.5*p2
    echo "ACT: ", ga, "\t", 0.5*p2, "\t", s0

    for n in 1..steps:
      #echo "pdiff: ", (p[0]-p[1]).norm2
      #echo "gdiff: ", (g[0]-g[1]).norm2
      #echo "ga: ", gaugeAction.gaugeAction1(gc,g)
      updateX(g,p,lambda*eps)
      #echo "pdiff: ", (p[0]-p[1]).norm2
      #echo "gdiff: ", (g[0]-g[1]).norm2
      #echo "ga: ", gaugeAction.gaugeAction1(gc,g)
      gc.updateP(g,p,0.5*eps)
      #echo "pdiff: ", (p[0]-p[1]).norm2
      #echo "gdiff: ", (g[0]-g[1]).norm2
      #echo "ga: ", gaugeAction.gaugeAction1(gc,g)
      updateX(g,p,(1.0-2.0*lambda)*eps)
      #echo "pdiff: ", (p[0]-p[1]).norm2
      #echo "gdiff: ", (g[0]-g[1]).norm2
      #echo "ga: ", gaugeAction.gaugeAction1(gc,g)
      gc.updateP(g,p,0.5*eps)
      #echo "pdiff: ", (p[0]-p[1]).norm2
      #echo "gdiff: ", (g[0]-g[1]).norm2
      #echo "ga: ", gaugeAction.gaugeAction1(gc,g)
      updateX(g,p,lambda*eps)

    let ga2 = gaugeAction.gaugeAction1(gc,g)
    p2 = 0.0
    for mu in 0..<p.len: p2 += p[mu].norm2
    let s2 = ga2 + 0.5*p2
    echo "ACT2: ", ga2, "\t", 0.5*p2, "\t", s2
    echo "dH: ", s2 - s0
    let sr = (s2-s0)/(eps*eps)
    echo "error rate: ", sr
    return (s2-s0)/s0

  proc testE4(lambda:float, steps=4):auto =
    # e_n = a*t^2/n^2 + b*t^4/n^4 + ...
    let
      e1 = test2(steps, lambda)
      e4 = test2(4*steps, lambda)
    return e1-16.0*e4

  var lambda = 0.23748

  when false:
    # Search for the lambda that cancels higher order terms.
    let tol = 1e-14
    var
      xlo = 0.15
      xhi = 0.25
      elo = testE4(xlo)
      ehi = testE4(xhi)
      x,e:float
    while elo>0 xor ehi>0:
      x = (ehi*xlo-elo*xhi)/(ehi-elo)
      e = testE4(x)
      echo "lambda: ",x," err_4/s: ",e
      if abs(e)<tol:
        break
      if e>0 xor ehi>0:
        xlo = x
        elo = e
      else:
        xhi = x
        ehi = e
    lambda = x

  test2(200, lambda)
  test2(20, lambda)
  test2(2, lambda)

  let dev = testE4(lambda,10)
  echo "Relative deviation from dt^2 scaling: ",dev
  if abs(dev)>1e-14:
    qexError "Large deviation."

  # echoTimers()
  qexFinalize()
