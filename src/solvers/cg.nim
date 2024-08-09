import base
import layout
import field
import solverBase
export solverBase

type
  CgPrecon* = enum
    cpNone,
    cpHerm,
    cpLeftRight,
    #cpRightNonHerm  # x = Ry, r = b-ARy, min r'A^-1r, p'r = 0
  CgState*[T] = object
    r*,Ap,b*: T
    p,x*: T
    z,q,LAp: T
    b2,r2,r2old,r2stop,rz,rzold,alpha: float
    iterations: int
    precon: CgPrecon

proc reset*(cgs: var CgState) =
  cgs.b2 = -1
  cgs.iterations = 0
  cgs.r2 = 1.0
  cgs.r2old = 1.0
  cgs.rzold = 1.0
  cgs.r2stop = 0.0

proc initPrecon*(state: var CgState) =
  case state.precon
  of cpNone:
    state.z = state.r
    state.q = state.p
    state.LAp = state.Ap
  of cpHerm:
    state.z = newOneof state.r
    state.q = state.p
    state.LAp = state.Ap
  of cpLeftRight:
    state.z = newOneof state.r
    state.q = newOneof state.p
    state.LAp = newOneOf state.z

proc newCgState*[T](x,b: T): CgState[T] =
  result.r = newOneOf b
  result.Ap = newOneOf b
  result.b = b
  result.p = newOneOf x
  result.x = x
  result.precon = cpNone
  result.initPrecon
  result.reset

# solves: A x = b
proc solve*(state: var CgState; op: auto; sp: var SolverParams) =
  mixin apply, applyPrecon
  tic("solve")
  let vrb = sp.verbosity
  template verb(n:int; body:untyped) =
    if vrb>=n: body
  let sub = sp.subset
  template subset(body:untyped) =
    onNoSync(sub):
      body
  template mythreads(body:untyped) =
    threads:
      onNoSync(sub):
        body

  let precon = op.precon
  if precon != state.precon:
    state.precon = precon
    state.initPrecon
    state.reset
  let
    r = state.r
    p = state.p
    Ap = state.Ap
    x = state.x
    b = state.b
    z = state.z
    q = state.q
    LAp = state.LAp
  var
    b2 = state.b2
    r2 = state.r2
    rz = state.rz
    #qLAp = state.qLAp

  if precon == cpHerm:
    when not compiles(op.applyPrecon(z, r)):
      qexError("cg.solve: precon == cpHerm but op.applyPrecon not found")
  if precon == cpLeftRight:
    when not compiles(op.applyPreconL(z, r)):
      qexError("cg.solve: precon == cpLeftRight but op.applyPreconL not found")
    when not compiles(op.applyPreconR(p, q)):
      qexError("cg.solve: precon == cpLeftRight but op.applyPreconR not found")

  template getRz =
    case precon
    of cpNone:
      rz = r2
    of cpHerm:
      subset:
        rz = r.redot z
    of cpLeftRight:
      subset:
        rz = z.norm2  # convenient to use rz here for z2
    #of cpRightNonHerm:
    #  subset:
    #    rz = Ap.dot z  # convenient to use rz here
  template preconL(z, r) =
    case precon
    of cpNone:
      discard
    of cpHerm:
      when compiles(op.applyPrecon(z, r)):
        op.applyPrecon(z, r)
    of cpLeftRight:
      when compiles(op.applyPreconL(z, r)):
        op.applyPreconL(z, r)
  template preconR(p, q) =
    case precon
    of cpNone:
      discard
    of cpHerm:
      discard
    of cpLeftRight:
      when compiles(op.applyPreconR(p, q)):
        op.applyPreconR(p, q)

  if b2<0:  # first call
    mythreads:
      b2 = b.norm2
    state.b2 = b2
    verb(1):
      echo("input norm2: ", b2)
    if b2 == 0.0:
      mythreads:
        x := 0
        r := 0
      r2 = 0.0
      rz = 0.0
    else:
      threads:
        op.apply(Ap, x)
        subset:
          r := b - Ap
          r2 = r.norm2
          p := 0
          verb(3):
            echo("p2: ", p.norm2)
            echo("r2: ", r2)

  let r2stop = sp.r2req * b2
  state.r2stop = r2stop
  let maxits = sp.maxits
  var itn0 = state.iterations
  var r2o0 = state.r2old
  var rzo0 = state.rzold
  var alpha0 = state.alpha
  toc("cg setup")

  if r2 > r2stop:
    threads:
      var itn = itn0
      var r2o = r2o0
      var rzo = rzo0
      var alpha = alpha0
      var qlap = 0.0
      verb(1):
        echo("CG iteration: ", itn, "  r2/b2: ", r2/b2)

      while itn<maxits and r2>r2stop:
        tic("cg loop")
        if itn == 0 or precon != cpLeftRight:
          preconL(z, r)  # z = L r  or  z = R r for RightNonHerm
        else:
          subset:
            z -= alpha * LAp
        getRz()  # r.z or z.z or Ap.z
        var beta = 0.0
        #if precon == cpRightNonHerm:
        #  beta = -rz / qLAp
        #else:
        beta = rz/rzo
        r2o = r2
        rzo = rz
        subset:
          if itn == 0:
            q := z
          else:
            q := z + beta*q
        toc("q update", flops=2*numNumbers(q[0])*sub.lenOuter)
        verb(3):
          echo "beta: ", beta
        preconR(p, q)  # p = R q
        toc("preconR")
        inc itn
        op.apply(Ap, p)
        toc("Ap")
        if precon == cpLeftRight:
          preconL(LAp, Ap)  # LAp = L Ap
        subset:
          #let pAp = p.redot(Ap)
          qLAp = q.redot(LAp)
          toc("qLAp", flops=2*numNumbers(p[0])*sub.lenOuter)
          alpha = rz/qLAp
          x += alpha*p
          toc("x", flops=2*numNumbers(p[0])*sub.lenOuter)
          r -= alpha*Ap
          toc("r", flops=2*numNumbers(r[0])*sub.lenOuter)
          r2 = r.norm2
          toc("r2", flops=2*numNumbers(r[0])*sub.lenOuter)
        verb(2):
          #echo(itn, " ", r2)
          echo("CG iteration: ", itn, "   r2/b2: ", r2/b2)
        verb(3):
          subset:
            #qLAp = q.redot(LAp)
            #echo "p2: ", p.norm2
            #echo "Ap2: ", Ap.norm2
            echo "rz: ", rz
            echo "qLAp: ", qLAp
            echo "alpha: ", alpha
            echo "x2: ", x.norm2
            echo "r2: ", r2
            echo "z2: ", z.norm2
            echo "q2: ", q.norm2
            echo "Ap2: ", Ap.norm2
            echo "LAp2: ", LAp.norm2
          op.apply(Ap, x)
          var fr2: float
          subset:
            threadBarrier()
            Ap -= b
            threadBarrier()
            fr2 = Ap.norm2
          echo "fr2: ", fr2, "  fr2/b2: ", fr2/b2
        if itn mod 64 == 0: aggregateTimers()
      toc("cg iterations")
      if threadNum==0:
        itn0 = itn
        r2o0 = r2o
        rzo0 = rzo
        alpha0 = alpha
      #var fr2: float
      #op.apply(Ap, x)
      #subset:
      #  r := b - Ap
      #  fr2 = r.norm2
      #verb(1):
      #  echo iterations, " acc r2: ", r2/b2
      #  echo iterations, " tru r2: ", fr2/b2

  state.iterations = itn0
  state.r2old = r2o0
  state.r2 = r2
  state.rzold = rzo0
  state.rz = rz
  state.alpha = alpha0
  #state.qLAp = qLAp
  verb(1):
    echo "CG final iterations: ", state.iterations, "  r2/b2: ", r2/b2
    #threads:
    #  op.apply(Ap, x)
    #  var fr2: float
    #  subset:
    #    fr2 = (b - Ap).norm2
    #  echo "   ", fr2/b2
  sp.finalIterations = state.iterations
  toc("cg final")

proc solve*(state: var CgState; x: Field; b: Field2; op: auto;
            sp: var SolverParams) =
  state.x = x
  state.b = b
  state.reset
  state.solve(op, sp)

proc cgSolve*(x: Field; b: Field2; op: auto; sp: var SolverParams) =
  var cg = newCgState(x, b)
  cg.solve x, b, op, sp

when isMainModule:
  import qex
  import physics/qcdTypes
  qexInit()
  echo "rank ", myRank, "/", nRanks
  #var lat = [8,8,8,8]
  var lat = [4,4,4,4]
  var lo = newLayout(lat)
  var m = lo.ColorMatrix()
  var v1 = lo.ColorVector()
  var v2 = lo.ColorVector()
  var v3 = lo.ColorVector()

  type opArgs = object
    m: type(m)
    precon: CgPrecon
  var oa = opArgs(m: m, precon: cpNone)
  proc apply*(oa: opArgs; r: type(v1); x: type(v1)) =
    r := oa.m*x
    #mul(r, m, x)

  type opArgsP = object
    m: type(m)
    precon: CgPrecon
  var oap = opArgsP(m: m, precon: cpHerm)
  proc apply*(oa: opArgsP; r: type(v1); x: type(v1)) =
    r := oa.m*x
    #mul(r, m, x)
  proc applyPrecon*(oa: opArgsP; r: type(v1); x: type(v1)) =
    for e in r:
      let t = sqrt(1.0 / m[e][0,0])
      r[e] := t * x[e]
    #mul(r, m, x)
  var precL = true
  var precR = true
  proc applyPreconL*(oa: opArgsP; r: type(v1); x: type(v1)) =
    if precL:
      applyPrecon(oa, r, x)
    else:
      r := x
  proc applyPreconR*(oa: opArgsP; r: type(v1); x: type(v1)) =
    if precR:
      applyPrecon(oa, r, x)
    else:
      r := x

  var sp:SolverParams
  sp.r2req = 1e-20
  sp.maxits = 200
  sp.verbosity = 2
  sp.subset.layoutSubset(lo, "all")
  threads:
    m.even := 1
    m.odd := 10
    threadBarrier()
    tfor i, 0..<lo.nSites:
      m{i} := i+1
    threadBarrier()
    v1.even := 1
    v1.odd := 2
    v2 := 0
    echo v1.norm2
    echo m.norm2
  template resid(r,b,x,oa: untyped) =
    oa.apply(r, x)
    r := b - r
  template checkResid(r, b, x, oa: auto) =
    resid(r, b, x, oa)
    let r2 = r.norm2
    let b2 = b.norm2
    echo "true r2/b2: ", r2/b2

  #cgSolve(v2, v1, oa, sp)
  var cg = newCgState(x=v2, b=v1)
  echo "starting cg.solve"
  cg.solve(oa, sp)
  checkResid(v3, v1, v2, oa)
  echo "end cg.solve iterations: ", sp.finalIterations

  echo "starting cg.solve restart test"
  v2 := 0
  cg.reset
  sp.maxits = 0
  sp.verbosity = 0
  while cg.r2 > cg.r2stop:
    sp.maxits += 10
    cg.solve(oa, sp)
    v3.resid(v1,v2,oa)
    let tr2 = v3.norm2
    echo cg.iterations, " ", cg.r2, "/", cg.r2stop, " ", tr2
    #cg.r := v3
    #cg.r2 = tr2
  echo sp.finalIterations, " ", cg.r2, "/", cg.r2stop
  checkResid(v3, v1, v2, oa)
  echo "end cg.solve restart test"

  echo "starting cg.solve restart test 2"
  v2 := 0
  cg.reset
  sp.maxits = 0
  while cg.r2 > cg.r2stop:
    sp.maxits += 10
    cg.solve(oa, sp)
    let c = cg.x.norm2
    echo cg.iterations, ": ", c, " ", cg.r2
  checkResid(v3, v1, v2, oa)
  echo "end cg.solve restart test 2"

  echo "starting cg.solve cpHerm restart test"
  v2 := 0
  cg.reset
  sp.maxits = 0
  sp.verbosity = 0
  while cg.r2 > cg.r2stop:
    sp.maxits += 10
    cg.solve(oap, sp)
    let c = cg.x.norm2
    echo cg.iterations, ": ", c, " ", cg.r2
  checkResid(v3, v1, v2, oa)
  echo "end cg.solve cpHerm restart test"

  echo "starting cg.solve cpLeftRight restart test"
  oap.precon = cpLeftRight
  v2 := 0
  cg.reset
  sp.maxits = 0
  sp.verbosity = 1
  while cg.r2 > cg.r2stop:
    sp.maxits += 10
    cg.solve(oap, sp)
    let c = cg.x.norm2
    echo cg.iterations, ": ", c, " ", cg.r2
  checkResid(v3, v1, v2, oa)
  echo "end cg.solve cpLeftRight restart test"

  echo "starting cg.solve cpLeftRight R restart test"
  oap.precon = cpLeftRight
  precL = false
  precR = true
  v2 := 0
  cg.reset
  sp.maxits = 0
  sp.verbosity = 0
  while cg.r2 > cg.r2stop:
    sp.maxits += 10
    cg.solve(oap, sp)
    let c = cg.x.norm2
    echo cg.iterations, ": ", c, " ", cg.r2
  checkResid(v3, v1, v2, oa)
  echo "end cg.solve cpLeftRight R restart test"

  echo "starting cg.solve cpLeftRight L restart test"
  oap.precon = cpLeftRight
  precL = true
  precR = false
  v2 := 0
  cg.reset
  sp.maxits = 0
  sp.verbosity = 0
  while cg.r2 > cg.r2stop:
    sp.maxits += 10
    cg.solve(oap, sp)
    let c = cg.x.norm2
    echo cg.iterations, ": ", c, " ", cg.r2
  checkResid(v3, v1, v2, oa)
  echo "end cg.solve cpLeftRight L restart test"
