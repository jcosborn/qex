import base
import layout
import layout/shifts
#import profile
#import stdUtils
#import gaugeUtils

proc startCornerShifts*[T](u0: openArray[T]): auto =
  var u = @u0
  var s:seq[seq[ShiftB[type(u[0][0])]]]
  let nd = u.len
  s.newSeq(nd)
  for mu in 0..<nd:
    s[mu].newSeq(nd)
    for nu in 0..<nd:
      if mu!=nu:
        s[mu][nu].initShiftB(u[mu], nu, 1, "all")
  threads:
    for mu in 0..<nd:
      for nu in 0..<nd:
        if mu!=nu:
          s[mu][nu].startSB(u[mu][ix])
  return s

proc startPlaqShifts*[T](u0: openArray[T]): auto =
  var u = @u0
  var sf:seq[seq[ShiftB[type(u[0][0])]]]
  var sb:seq[seq[ShiftB[type(u[0][0])]]]
  let nd = u.len
  sf.newSeq(nd)
  sb.newSeq(nd)
  for mu in 0..<nd:
    sf[mu].newSeq(nd)
    sb[mu].newSeq(nd)
    for nu in 0..<nd:
      sb[mu][nu].initShiftB(u[mu], nu, -1, "all")
      if mu!=nu:
        sf[mu][nu].initShiftB(u[mu], nu, 1, "all")
  threads:
    for mu in 0..<nd:
      sb[mu][mu].startSB(u[mu][ix])
    for mu in 0..<nd:
      for nu in 0..<nd:
        if mu!=nu:
          sf[mu][nu].startSB(u[mu][ix])
          sb[mu][nu].startSB(u[mu][ix])
  return (sf,sb)

proc finishPlaqShifts*[T](u0,b0: openArray[T], sf,sb: seq): auto =
  var u = @u0
  var b = @b0
  var sc:seq[seq[ShiftB[type(u[0][0])]]]
  let nd = u.len
  sc.newSeq(nd)
  for mu in 0..<nd:
    sc[mu].newSeq(nd)
    for nu in 0..<nd:
      if mu!=nu:
        sc[mu][nu].initShiftB(u[mu], nu, 1, "all")
  threads:
    for mu in 0..<nd:
      for ir in b[mu]:
        localSB(sb[mu][mu], ir, assign(b[mu][ir], it), u[mu][ix])
    for mu in 0..<nd:
      boundarySB(sb[mu][mu], assign(b[mu][ir], it))
    threadBarrier()
    for mu in 0..<nd:
      for nu in 0..<nd:
        if mu!=nu:
          sc[mu][nu].startSB(b[mu][ix])
  return (sc)

proc startStapleShifts*[T](u: openArray[T]): auto =
  var s:seq[seq[seq[ShiftB[type(u[0][0][0])]]]]
  let nd = u.len
  s.newSeq(nd)
  for mu in 0..<nd:
    s[mu].newSeq(nd)
    for nu in 0..<nd:
      if nu!=mu:
        s[mu][nu].newSeq(nd)
        for sig in 0..<nd:
          if sig!=mu:
            s[mu][nu][sig].initShiftB(u[mu][nu], sig, 1, "all")
            s[mu][nu][sig].startSB(u[mu][nu][ix])
  return s

proc makeFwdStaples*[T](uu: openArray[T], s: auto): auto =
  mixin mul
  tic()
  let u = cast[ptr cArray[T]](unsafeAddr(uu[0]))
  let lo = u[0].l
  let nd = lo.nDim
  let nc = u[0][0].ncols
  let flops = lo.nSites.float*float(nd*(nd-1)*3*(4*nc-1)*nc*nc)
  var st: seq[seq[type(uu[0])]]
  st.newSeq(nd)
  for mu in 0..<nd:
    st[mu].newSeq(nd)
    for nu in 0..<nd:
      if mu!=nu:
        st[mu][nu].new(lo)
  toc("makeStaples setup")
  threads:
    tic()
    var umu,unu,umunu: type(load1(u[0][0]))
    for ir in lo:
      for mu in 1..<nd:
        for nu in 0..<mu:
          if isLocal(s[mu][nu],ir) and isLocal(s[nu][mu],ir):
            localSB(s[mu][nu], ir, assign(umu,it), u[mu][ix])
            localSB(s[nu][mu], ir, assign(unu,it), u[nu][ix])
            mul(umunu, umu, unu.adj)
            mul(st[mu][nu][ir], u[nu][ir], umunu)
            mul(st[nu][mu][ir], u[mu][ir], umunu.adj)
    toc("makeStaples local")
    #[
    var needBoundary = false
    for mu in 0..<nd:
      for nu in 0..<nd:
        if mu != nu:
          boundaryWaitSB(s[mu][nu]): needBoundary = true
    toc("makeStaples wait")
    if needBoundary:
      boundarySyncSB()
      for ir in lo:
        for mu in 1..<nd:
          for nu in 0..<mu:
            if not isLocal(s[mu][nu],ir) or not isLocal(s[nu][mu],ir):
              getSB(s[mu][nu], ir, assign(umu,it), u[mu][ix])
              getSB(s[nu][mu], ir, assign(unu,it), u[nu][ix])
              mul(umunu, umu, unu.adj)
              mul(st[mu][nu][ir], u[nu][ir], umunu)
              mul(st[nu][mu][ir], u[mu][ir], umunu.adj)
    ]#
    for mu in 1..<nd:
      for nu in 0..<mu:
        var needBoundary = false
        boundaryWaitSB(s[mu][nu]): needBoundary = true
        boundaryWaitSB(s[nu][mu]): needBoundary = true
        if needBoundary:
          for ir in lo:
            if not isLocal(s[mu][nu],ir) or not isLocal(s[nu][mu],ir):
              getSB(s[mu][nu], ir, assign(umu,it), u[mu][ix])
              getSB(s[nu][mu], ir, assign(unu,it), u[nu][ix])
              mul(umunu, umu, unu.adj)
              mul(st[mu][nu][ir], u[nu][ir], umunu)
              mul(st[nu][mu][ir], u[mu][ir], umunu.adj)
    toc("makeStaples boundary")
  toc("makeStaples threads", flops=flops)
  return st

proc makeStaples*[T](uu: openArray[T], s: auto): auto =
  ## sft: fwd staples
  ## stu: bck staples, offset up
  ## ss: shifts for stu
  mixin mul, load1, adj
  tic()
  let u = cast[ptr cArray[T]](unsafeAddr(uu[0]))
  let lo = u[0].l
  let nd = lo.nDim
  let nc = u[0][0].ncols
  let flops = lo.nSites.float*float(nd*(nd-1)*6*(4*nc-1)*nc*nc)
  var
    stf = newFieldArray2(lo,type(uu[0]),[nd,nd],mu!=nu)
    stu = newFieldArray2(lo,type(uu[0]),[nd,nd],mu!=nu)
  var ss: seq[seq[ShiftB[type(uu[0][0])]]]
  ss.newSeq(nd)
  for mu in 0..<nd:
    ss[mu].newSeq(nd)
    for nu in 0..<nd:
      if mu!=nu:
        ss[mu][nu].initShiftB(stu[mu,nu], nu, -1, "all")
  toc("makeStaples setup")
  threads:
    tic()
    var umu,unu,umunu,unumu {.noInit.}: evalType(u[0][0])
    for ir in lo:
      for mu in 1..<nd:
        for nu in 0..<mu:
          if isLocal(s[mu][nu],ir) and isLocal(s[nu][mu],ir):
            localSB(s[mu][nu], ir, assign(umu,it), u[mu][ix])
            localSB(s[nu][mu], ir, assign(unu,it), u[nu][ix])
            mul(umunu, umu, unu.adj)
            mul(stf[mu,nu][ir], u[nu][ir], umunu)
            mul(stf[nu,mu][ir], u[mu][ir], umunu.adj)
            mul(unumu, u[nu][ir].adj, u[mu][ir])
            mul(stu[mu,nu][ir], unumu, unu)
            mul(stu[nu,mu][ir], unumu.adj, umu)
    toc("makeStaples local")
    var needBoundary = false
    for mu in 1..<nd:
      for nu in 0..<mu:
        var needBoundaryU = false
        boundaryWaitSB(s[mu][nu]): needBoundaryU = true
        boundaryWaitSB(s[nu][mu]): needBoundaryU = true
        needBoundary = needBoundary or needBoundaryU
        if needBoundaryU:
          boundarySyncSB()
          for ir in lo:
            if not isLocal(s[mu][nu],ir) or not isLocal(s[nu][mu],ir):
              getSB(s[mu][nu], ir, assign(umu,it), u[mu][ix])
              getSB(s[nu][mu], ir, assign(unu,it), u[nu][ix])
              mul(unumu, u[nu][ir].adj, u[mu][ir])
              mul(stu[mu,nu][ir], unumu, unu)
              mul(stu[nu,mu][ir], unumu.adj, umu)
        threadBarrier()
        ss[mu][nu].startSB(stu[mu,nu][ix])
        ss[nu][mu].startSB(stu[nu,mu][ix])
    toc("makeStaplesU boundary")
    if needBoundary:
      boundarySyncSB()
      for ir in lo:
        for mu in 1..<nd:
          for nu in 0..<mu:
            if not isLocal(s[mu][nu],ir) or not isLocal(s[nu][mu],ir):
              getSB(s[mu][nu], ir, assign(umu,it), u[mu][ix])
              getSB(s[nu][mu], ir, assign(unu,it), u[nu][ix])
              mul(umunu, umu, unu.adj)
              mul(stf[mu,nu][ir], u[nu][ir], umunu)
              mul(stf[nu,mu][ir], u[mu][ir], umunu.adj)
    toc("makeStaplesF boundary")
  toc("makeStaples threads", flops=flops)
  return (stf,stu,ss)


when isMainModule:
  import qex
  import physics/qcdTypes

  qexInit()
  #var defaultGaugeFile = "l88.scidac"
  #let defaultLat = @[2,2,2,2]
  let defaultLat = @[8,8,8,8]
  defaultSetup()
  g.unit
  #g.random

  proc test(g: auto) =
    tic("test")
    var cs = startCornerShifts(g)
    toc("startCornerShifts")
    var (stf,stu,ss) = makeStaples(g, cs)
    toc("makeStaples")

  test(g)
  resetTimers()
  test(g)
  echoTimers()

  qexFinalize()
