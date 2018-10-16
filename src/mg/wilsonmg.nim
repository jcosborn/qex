import base
import layout
import field
import solvers/gcr
import solvers/cgls
import physics/qcdTypes
import physics/wilsonD
import mg/mgblocks
import mg/mgargs
import mg/wmgutils
import rng
import mg/wmgsetup

type WilsonMg*[W] = object
  w: W
  sp: SolverParams
  #solver: S
  #op: O

proc newWilsonMg[W: WilsonD](w: W): WilsonMg[W] =
  result.w = w
  result.sp.init

proc apply*(o: tuple; r: any; x: any) =
  mixin D
  threadBarrier()
  o.w.D(r, x, o.m)
  threadBarrier()
proc applyAdj*(o: tuple; r: any; x: any) =
  mixin Ddag
  threadBarrier()
  o.w.Ddag(r, x, o.m)
  threadBarrier()

type opArgs[T] = object
  a: T
proc apply*(oa: opArgs; r: any; x: any) =
  oa.a.apply(r, x)
proc applyAdj*(oa: opArgs; r: any; x: any) =
  oa.a.applyAdj(r, x)
proc preconditioner*(oa: opArgs; z: any; gs: GcrState) =
  z := gs.r

proc solveCgls(w: Wilson, x: Field, b: Field2, m: float,
               sp: var SolverParams) =
  sp.subset.layoutSubset(x.l, "all")
  let t = (w:w, m:m)
  var oa = opArgs[type(t)](a: t)
  var zero = x.newOneOf()
  zero := 0
  var cgls = newCglsState(x=x, b1=zero, b2=b)
  cgls.solve(oa, sp)
  echo sp.finalIterations

proc solveGcr*(w: Wilson; x,b: Field; m: SomeNumber; sp: var SolverParams) =
  sp.subset.layoutSubset(x.l, "all")
  let t = (w:w, m:m)
  var oa = opArgs[type(t)](a: t)
  var gcr = newGcrState(x=x, b=b)
  gcr.solve(oa, sp)
  echo sp.finalIterations

type Weo[T] = object
  w: T
template D*(ww: Weo, r, x: Field; m: SomeNumber) =
  let m4 = 4.0 + m
  let m2 = m4*m4
  threadBarrier()
  wilsonD2ee(ww.w.se, ww.w.so, r, ww.w.g, x, m2)
  threadBarrier()
proc solveGcrEo*(w: Wilson; x,b: Field; m: SomeNumber; sp: var SolverParams) =
  sp.subset.layoutSubset(x.l, "even")
  let wx = Weo[type(w)](w: w)
  let t = (w:wx, m:m)
  var oa = opArgs[type(t)](a: t)
  let r = b.newOneOf
  threads:
    w.eoReduce(r, b, m)
  var gcr = newGcrState(x=x, b=r)
  gcr.solve(oa, sp)
  threads:
    w.eoReconstruct(x, b, m)
  echo sp.finalIterations

proc normalize(x: Field) =
  let t = x.norm2
  let s = 1/sqrt(t)
  x := s*x

proc projectOut(x: Field, y: Field2) =
  let yx = dot(y,x)
  let yy = y.norm2
  let s = yx / yy
  x -= s * y

proc Anorm(op: any, x: Field, y: Field2): auto =
  mixin apply
  op.apply(y, x)
  let x2 = x.norm2
  let y2 = y.norm2
  y2/x2

proc mgsetup(r,p: var MgTransfer, op: any, x: Field) =
  let nmgv1 = p.v[0].len
  var b = newOneOf(x)
  var t = newOneOf(x)
  var sp: SolverParams
  sp.maxits = 100
  sp.verbosity = 0
  sp.r2req = 0.0
  sp.subset.layoutSubset(x.l, "all")
  var rs = newRNGField(RngMilc6, x.l, 987654321)
  r.v.wmgzero
  p.v.wmgzero
  x.gaussian rs
  for i in 0..<nmgv1:
    echo "setup: ", i
    for k in 1..4:
      b := x
      #p.wmgBlockProject(b, t, op.cb)
      b.normalize
      b.wmgProject(p)
      b.normalize
      b.wmgProject(p)
      b.normalize
      #x := 0
      #op.a.w.solveGcr(x, b, op.a.m, sp)
      t := gamma5 * b
      x := 0
      op.a.w.solveGcr(x, t, op.a.m, sp)
      t := gamma5 * x
      x := 0
      op.a.w.solveGcr(x, t, op.a.m, sp)
      #op.a.w.solveGcrEo(x, b, op.a.m, sp)
      echo "Anorm:  ", op.Anorm(x, t)
      x.normalize
      x.wmgProject(p)
      x.normalize
      x.wmgProject(p)
      x.normalize
      echo "Anormp: ", op.Anorm(x, t)

    #p.wmgBlockNormalizeInsert(x, i, t, op.cb)
    p.v.wmgInsert(x, i)

    mixin apply
    #op.apply(t, x)
    t := x
    t.wmgProject(r)
    t.normalize
    #r.wmgBlockNormalizeInsert(t, i, x, op.cb)
    r.v.wmgInsert(t, i)

    x := b

type opArgs2[T] = object
  a: T
proc apply*(oa: opArgs2; r: any; x: any) =
  threadBarrier()
  oa.a.p.prolong(oa.a.f, x)
  threadBarrier()
  oa.a.apply(r, oa.a.f)
  threadBarrier()
proc applyDag*(oa: opArgs2; r: any; x: any) =
  mixin applyDag
  threadBarrier()
  oa.a.applyDag(oa.a.f, x)
  threadBarrier()
  oa.a.p.restrict(r, oa.a.f)
  threadBarrier()
proc preconditioner*(oa: opArgs2; z: any; gs: GcrState) =
  #z := gs.r
  threadBarrier()
  oa.applyDag(z, gs.r)
  threadBarrier()

type opArgs3[T] = object
  a: T
proc apply*(oa: opArgs3; r: any; x: any) =
  threadBarrier()
  oa.a.p.prolong(oa.a.f, x)
  threadBarrier()
  oa.a.apply(oa.a.f2, oa.a.f)
  threadBarrier()
  oa.a.r.restrict(r, oa.a.f2)
  threadBarrier()
proc applyDag*(oa: opArgs3; r: any; x: any) =
  mixin applyDag
  threadBarrier()
  oa.a.r.prolong(oa.a.f, x)
  threadBarrier()
  oa.a.applyDag(oa.a.f2, oa.a.f)
  threadBarrier()
  oa.a.p.restrict(r, oa.a.f2)
  threadBarrier()
proc preconditioner*(oa: opArgs3; z: any; gs: GcrState) =
  z := gs.r
  #threadBarrier()
  #oa.applyDag(z, gs.r)
  #threadBarrier()

type opArgsVc[T,C,F,P,R] = object
  a: T
  cb,cx: C
  f,f2: F
  p: P
  r: R
  lp,nv,csType,eo: int
proc apply*(oa: opArgsVc; r: any; x: any) =
  mixin D
  if oa.eo == 0:
    threadBarrier()
    oa.a.w.D(r, x, oa.a.m)
    threadBarrier()
  else:
    let m4 = 4.0 + oa.a.m
    let m2 = m4*m4
    threadBarrier()
    wilsonD2ee(oa.a.w.se, oa.a.w.so, r, oa.a.w.g, x, m2)
    threadBarrier()
proc applyDag*(oa: opArgsVc; r: any; x: any) =
  mixin Ddag
  threadBarrier()
  oa.a.w.Ddag(r, x, oa.a.m)
  threadBarrier()
proc smoother*(oa: opArgsVc; x: any; b: any) =
  var sp: SolverParams
  sp.maxits = 6
  sp.verbosity = 2
  sp.r2req = 0.01
  if oa.eo == 0:
    sp.subset.layoutSubset(x.l, "all")
    var gcr = newGcrState(x=x, b=b)
    echo "starting gcr"
    var oa2 = opArgs[type(oa.a)](a: oa.a)
    gcr.solve(oa2, sp)
    echo sp.finalIterations
    echo gcr.r2
  else:
    sp.subset.layoutSubset(x.l, "even")
    var gcr = newGcrState(x=x, b=b)
    echo "starting gcr"
    var oa2 = opArgs[type(oa)](a: oa)
    gcr.solve(oa2, sp)
    echo sp.finalIterations
    echo gcr.r2
proc csolve*(oa: opArgsVc; x: any; b: any) =
  var sp: SolverParams
  sp.maxits = 100
  sp.verbosity = 2
  sp.r2req = 0.01
  sp.subset.layoutSubset(x.l, "all")
  var gcr = newGcrState(x=x, b=b)
  echo "starting gcr"
  var oa2 = opArgs2[type(oa)](a: oa)
  x := 0
  gcr.solve(oa2, sp)
  echo sp.finalIterations
  echo gcr.r2
proc csolve1*(oa: opArgsVc; x: any; b: any) =
  var sp: SolverParams
  sp.maxits = 100
  sp.verbosity = 2
  sp.r2req = 0.01
  sp.subset.layoutSubset(x.l, "all")
  var gcr = newGcrState(x=x, b=b)
  echo "starting gcr"
  var oa3 = opArgs3[type(oa)](a: oa)
  x := 0
  gcr.solve(oa3, sp)
  echo sp.finalIterations
  echo gcr.r2
proc csolve2*(oa: opArgsVc; x: any; b: any) =
  var sp: SolverParams
  sp.maxits = 100
  sp.verbosity = 2
  sp.r2req = 0.01
  sp.subset.layoutSubset(x.l, "all")
  var gcr = newGcrState(x=x, b=b)
  echo "starting gcr"
  var oa3 = opArgs3[type(oa)](a: oa)
  x := 0
  gcr.solve(oa3, sp)
  echo sp.finalIterations
  echo gcr.r2
proc preconditioner*(oa: var opArgsVc; z: var any; gs: GcrState) =
  if oa.eo == 0: # all sites
    # z = r (r'A'r)/(r'A'Ar)
    oa.apply(oa.f, gs.r)
    let rar = dot(oa.f, gs.r)
    let raar = oa.f.norm2
    let c = rar/raar
    z := c * gs.r
    oa.apply(oa.f, z)
    oa.f -= gs.r
    echo "r: ", oa.f.norm2

    z := gs.r
    z.wmgProject(oa.r)
    echo "wmgProj: ", z.norm2

    oa.applyDag(z, gs.r)
    z.wmgBlockProject(oa.p, oa.f, oa.cb)
    echo "wmgBlockProj: ", z.norm2

    #oa.p.restrict(oa.cb, gs.r)
    #oa.cx := oa.cb
    #[
    var p = oa.p
    z := gs.r
    #p.v.wmgzero(0)
    let lp = oa.lp
    echo "lp: ", lp
    p.wmgBlockNormalizeInsert(z, lp, oa.f, oa.cb)
    oa.lp = (lp+1) mod oa.nv
    oa.p.restrict(oa.cx, gs.r)
    oa.p.prolong(z, oa.cx)
    z -= gs.r
    echo "z: ", z.norm2
    ]#

    case oa.csType:
    of 0:
      oa.csolve(oa.cx, gs.r)
      #oa.gcrc.solve(oa.c, oa.sp)
      oa.p.prolong(z, oa.cx)
    else:
      oa.r.restrict(oa.cb, gs.r)
      oa.csolve1(oa.cx, oa.cb)
      #oa.gcrc.solve(oa.c, oa.sp)
      oa.p.prolong(z, oa.cx)

    oa.apply(oa.f, z)
    oa.f -= gs.r
    echo "csolve: ", oa.f.norm2
    oa.smoother(z, gs.r)
    oa.apply(oa.f, z)
    oa.f -= gs.r
    echo "prec: ", oa.f.norm2
    #z := gs.r
  else:  # even only
    #z := gs.r
    oa.cb := 0
    oa.r.restrict(oa.cb, gs.r)
    oa.cx := 0
    oa.csolve2(oa.cx, oa.cb)
    z := 0
    oa.p.prolong(z, oa.cx)
    z.odd := 0

    oa.f := 0
    oa.apply(oa.f, z)
    oa.f -= gs.r
    echo "csolve: ", oa.f.norm2
    oa.smoother(z, gs.r)
    oa.f := 0
    oa.apply(oa.f, z)
    oa.f -= gs.r
    echo "prec: ", oa.f.norm2
proc solveGcrVc*(w: Wilson; x: var Field; b: Field; m: SomeNumber;
                 sp: var SolverParams, cst=0) =
  sp.subset.layoutSubset(x.l, "all")
  let t = (w:w, m:m)
  #let rs = b.toSingle.newOneOf
  #let zs = b.toSingle.newOneOf
  let loF = x.l
  let latC = [4,4,4,4]
  let loC = newLayout(latC, loF.V, loF.rankGeom, loF.innerGeom)
  let mgb1 = newMgBlock(loF, loC)
  const nmgv1 {.intDefine.} = 48
  var rv,pv: LatticeWmgVectorsV[nmgv1]
  rv.new(loF)
  pv.new(loF)
  var r = newMgTransfer(mgb1, rv)
  var p = newMgTransfer(mgb1, pv)

  var f = x.newOneOf()
  var f2 = x.newOneOf()
  var cb: LatticeMgColorVectorV[nmgv1]
  cb.new(loC)
  var cx: LatticeMgColorVectorV[nmgv1]
  cx.new(loC)
  var oa = opArgsVc[type(t),type(cb),type(f),type(p),type(r)](
                    a: t, cb: cb, cx: cx, f: f, f2: f2, p: p, r: r)
  oa.nv = nmgv1
  oa.csType = cst
  oa.eo = 0
  let setupKind = intParam("s", 1)
  case setupKind:
    of 0: mgsetup(r, p, oa, x)
    else: mgsetupSvd(r, p, oa, x)
  var gcr = newGcrState(x=x, b=b)
  x := 0
  gcr.solve(oa, sp)
  echo sp.finalIterations

proc solveGcrEoVc*(w: Wilson; x: Field; b: Field; m: SomeNumber;
                   sp: var SolverParams, cst=0) =
  sp.subset.layoutSubset(x.l, "even")
  let t = (w:w, m:m)
  #let rs = b.toSingle.newOneOf
  #let zs = b.toSingle.newOneOf
  let loF = x.l
  let latC = [4,4,4,4]
  let loC = newLayout(latC, loF.V, loF.rankGeom, loF.innerGeom)
  let mgb1 = newMgBlock(loF, loC)
  const nmgv1 {.intDefine.} = 48
  var rv,pv: LatticeWmgVectorsV[nmgv1]
  rv.new(loF)
  pv.new(loF)
  var r = newMgTransfer(mgb1, rv)
  var p = newMgTransfer(mgb1, pv)

  var f = x.newOneOf()
  var f2 = x.newOneOf()
  var cb: LatticeMgColorVectorV[nmgv1]
  cb.new(loC)
  var cx: LatticeMgColorVectorV[nmgv1]
  cx.new(loC)
  var oa = opArgsVc[type(t),type(cb),type(f),type(p),type(r)](
                    a: t, cb: cb, cx: cx, f: f, f2: f2, p: p, r: r)
  oa.nv = nmgv1
  oa.csType = cst
  oa.eo = 0  # no EO for setup
  let setupKind = intParam("s", 1)
  case setupKind:
    of 0: mgsetup(r, p, oa, x)
    else: mgsetupSvd(r, p, oa, x)
  oa.eo = 1
  var be = b.newOneOf()
  be := 0
  threads:
    w.eoReduce(be, b, m)
  var gcr = newGcrState(x=x, b=be)
  x := 0
  gcr.solve(oa, sp)
  threads:
    w.eoReconstruct(x, b, m)
  echo sp.finalIterations




when isMainModule:
  import qex
  qexInit()
  var defaultLat = @[16,16,16,16]
  defaultSetup()
  var rs = newRNGField(RngMilc6, lo, intParam("seed", 987654321).uint64)
  var warmT = floatParam("warmt", 0.1)
  threads:
    if fn == "":
      g.warm warmT, rs
    g.setBC

  var src1 = lo.DiracFermion()
  var src2 = lo.DiracFermion()
  var soln1 = lo.DiracFermion()
  var soln2 = lo.DiracFermion()
  var v1 = lo.DiracFermion()
  var v2 = lo.DiracFermion()
  var r = lo.DiracFermion()
  var tv = lo.DiracFermion()

  var w = newWilson(g, v1)
  echo "done newWilson"
  var m = floatParam("mass", 0.1)
  echo "mass: ", m

  threads:
    src1 := 0
    echo "Init Wilson D:"
    w.D(soln1, src1, m)
    echo "... done"

  var sp = initSolverParams()
  sp.maxits = int(1e9/lo.physVol.float)
  sp.verbosity = intParam("spverb", 0)
  sp.r2req = 1e-16
  echo sp
  if myRank==0:
    src1{0}[0][0] := 1
  w.solveEO(soln1, src1, m, sp)
  threads:
    soln2 := src1
    threadBarrier()
    w.D(src1, soln1, m)
    w.D(src2, soln2, m)

  proc checkResid(src,soln,x: any) =
    threads:
      #echo "v2: ", v2.norm2
      #echo "v2.even: ", v2.even.norm2
      #echo "v2.odd: ", v2.odd.norm2
      w.D(r, x, m)
      threadBarrier()
      r := src - r
      let r2 = r.norm2
      tv := soln - x
      let e2 = tv.norm2
      echo "r2: ", r2, "  e2: ", e2

      #var fr20 = relResid(r, v2, 1e-20)
      #var fr30 = relResid(r, v2, 1e-30)
      #echo "fnalResid: ", fr20, "  ", fr30

  proc test1(src,soln: any) =
    #[
    echo "EO biCGStab:"
    threads:
      v2 := 0
    #w.solveEO(v2, src, m, sp)
    #resetTimers()
    w.solve(v2, src, m, sp)
    #echoTimers()
    checkResid(src, soln, v2)
    #var fr = relResid(r,v2,1e-30)
    #echo "fnalResid: ", fr
    ]#

    echo "EO CG:"
    threads:
      v2 := 0
    #resetTimers()
    w.solveEO(v2, src, m, sp)
    #echoTimers()
    checkResid(src, soln, v2)
    #fr = relResid(r,v2,1e-30)
    #echo "fnalResid: ", fr

    echo "CGLS:"
    threads:
      v2 := 0
    #resetTimers()
    w.solveCgls(v2, src, m, sp)
    #echoTimers()
    checkResid(src, soln, v2)
    #fr = relResid(r,v2,1e-30)
    #echo "fnalResid: ", fr

    echo "GCR:"
    threads:
      v2 := 0
    #resetTimers()
    w.solveGcr(v2, src, m, sp)
    #echoTimers()
    checkResid(src, soln, v2)
    #fr = relResid(r,v2,1e-30)
    #echo "fnalResid: ", fr

    echo "EO GCR:"
    threads:
      v2 := 0
    #resetTimers()
    w.solveGcrEo(v2, src, m, sp)
    #echoTimers()
    checkResid(src, soln, v2)
    #fr = relResid(r,v2,1e-30)
    #echo "fnalResid: ", fr

    #[
    echo "GCR VC:"
    threads:
      v2 := 0
    #resetTimers()
    sp.verbosity = 3
    w.solveGcrVc(v2, src, m, sp, 1)
    sp.verbosity = 1
    #echoTimers()
    checkResid(src, soln, v2)
    #fr = relResid(r,v2,1e-30)
    #echo "fnalResid: ", fr
    ]#

    echo "GCR EO VC:"
    threads:
      v2 := 0
    #resetTimers()
    sp.verbosity = 3
    w.solveGcrEoVc(v2, src, m, sp, 1)
    sp.verbosity = 1
    #echoTimers()
    checkResid(src, soln, v2)
    #fr = relResid(r,v2,1e-30)
    #echo "fnalResid: ", fr

  sp.verbosity = 1

  #sp.r2req = 1e-12
  #test1(src1, soln1)

  sp.r2req = 1e-16
  test1(src1, soln1)

  #echoTimers()
  qexFinalize()


#[
  import qex
  import physics/qcdTypes
  qexInit()
  echo "rank ", myRank, "/", nRanks
  #var lat = [8,8,8,8]
  var lat = [4,4,4,4]
  var lo = newLayout(lat)
  var m = lo.ColorMatrixD()
  var v1 = lo.ColorVectorD()
  var v2 = lo.ColorVectorD()
  var v3 = lo.ColorVectorD()
  type opArgs = object
    m: type(m)
  var oa = opArgs(m: m)
  proc apply*(oa: opArgs; r: type(v1); x: type(v1)) =
    r := oa.m*x
  #proc apply2*(oa: opArgs; r: type(v1); x: type(v1)) =
  #  r := oa.m*oa.m*x
  proc preconditioner*(oa: opArgs; r: type(v1); gs: GcrState) =
    #r := x
    #r := oa.m*gs.r
    #for i in r:
    #  let t = oa.m[i].norm2
    #  let s = 1.0/sqrt(t)
    #  r[i] := asReal(s)*x[i]
    #let a = 1.0/sqrt(gs.r2)
    for i in r:
    #  #let x2 = gs.x[i].norm2 + 1e-30
      let r2 = gs.r[i].norm2
      #let s = 1.0
      #let s = sqrt(r2)
      let s = r2
      #let s = a*r2
      #let s = r2*sqrt(r2)
      r[i] := asReal(s)*gs.r[i]
      #r[i] := asReal(s)*(oa.m[i]*gs.r[i])
  template resid(r,b,x,oa: untyped) =
    oa.apply(r, x)
    r := b - r
  var sp: SolverParams
  sp.r2req = 1e-30
  sp.maxits = 200
  sp.verbosity = 3
  sp.subset.layoutSubset(lo, "all")
  threads:
    m.even := 1
    m.odd := 10
    threadBarrier()
    tfor i, 0, lo.nSites-1:
      m{i} := i+1
      #m{i} := sqrt(sqrt(i+0.01))
    threadBarrier()
    v1.even := 1
    v1.odd := 2
    v2 := 0
    threadBarrier()
    echo v1.norm2
    echo m.norm2

  var gcr = newGcrState(x=v2, b=v1)
  gcr.solve(oa, sp)
  echo sp.finalIterations
  v3.resid(v1,v2,oa)
  echo "rsq: ", v3.norm2/gcr.b2
]#
