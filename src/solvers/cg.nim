import base
import layout
import field
import solverBase
export solverBase

type
  CgState*[T] = object
    r,Ap,b: T
    p,x: T
    b2,r2,r2old,r2stop: float
    iterations: int

proc reset*(cgs: var CgState) =
  cgs.b2 = -1
  cgs.iterations = 0
  cgs.r2old = 1.0
  cgs.r2stop = 0.0

proc newCgState*[T](x,b: T): CgState[T] =
  result.r = newOneOf(b)
  result.Ap = newOneOf(b)
  result.b = b
  result.p = newOneOf(x)
  result.x = x
  result.reset

# solves: A x = b
proc solve*(state: var CgState; op: auto; sp: var SolverParams) =
  mixin apply
  tic()
  let vrb = sp.verbosity
  template verb(n:int; body:untyped):untyped =
    if vrb>=n: body
  let sub = sp.subset
  template subset(body:untyped):untyped =
    onNoSync(sub):
      body
  template mythreads(body:untyped):untyped =
    threads:
      onNoSync(sub):
        body

  let
    r = state.r
    p = state.p
    Ap = state.Ap
    x = state.x
    b = state.b
  var
    b2 = state.b2
    r2 = state.r2

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
    else:
      threads:
        op.apply(Ap, x)
        subset:
          r := b - Ap
          p := 0
          r2 = r.norm2
          verb(3):
            echo("p2: ", p.norm2)
            echo("r2: ", r2)

  let r2stop = sp.r2req * b2
  state.r2stop = r2stop
  let maxits = sp.maxits
  var itn0 = state.iterations
  var r2o0 = state.r2old

  toc("cg setup")
  if r2 > r2stop:
    threads:
      var itn = itn0
      var r2o = r2o0
      verb(1):
        #echo(-1, " ", r2)
        echo(itn, " ", r2/b2)

      while itn<maxits and r2>r2stop:
        tic()
        let beta = r2/r2o
        r2o = r2
        subset:
          p := r + beta*p
        toc("p update", flops=2*numNumbers(r[0])*sub.lenOuter)
        verb(3):
          echo "beta: ", beta
        inc itn
        op.apply(Ap, p)
        toc("Ap")
        subset:
          let pAp = p.redot(Ap)
          toc("pAp", flops=2*numNumbers(p[0])*sub.lenOuter)
          let alpha = r2/pAp
          x += alpha*p
          toc("x", flops=2*numNumbers(p[0])*sub.lenOuter)
          r -= alpha*Ap
          toc("r", flops=2*numNumbers(r[0])*sub.lenOuter)
          r2 = r.norm2
          toc("r2", flops=2*numNumbers(r[0])*sub.lenOuter)
        verb(2):
          #echo(itn, " ", r2)
          echo(itn, " ", r2/b2)
        verb(3):
          subset:
            let pAp = p.redot(Ap)
            echo "p2: ", p.norm2
            echo "Ap2: ", Ap.norm2
            echo "pAp: ", pAp
            echo "alpha: ", r2o/pAp
            echo "x2: ", x.norm2
            echo "r2: ", r2
          op.apply(Ap, x)
          var fr2: float
          subset:
            fr2 = (b - Ap).norm2
          echo "    ", fr2, "    ", fr2/b2
        if itn mod 64 == 0: aggregateTimers()
      toc("cg iterations")
      if threadNum==0:
        itn0 = itn
        r2o0 = r2o
      #var fr2: float
      #op.apply(Ap, x)
      #subset:
      #  r := b - Ap
      #  fr2 = r.norm2
      #verb(1):
      #  echo iterations, " acc r2:", r2/b2
      #  echo iterations, " tru r2:", fr2/b2

  state.iterations = itn0
  state.r2old = r2o0
  state.r2 = r2
  verb(1):
    echo state.iterations, " acc r2:", r2/b2
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
  var oa = opArgs(m: m)
  proc apply*(oa: opArgs; r: type(v1); x: type(v1)) =
    r := oa.m*x
    #mul(r, m, x)
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

  #cgSolve(v2, v1, oa, sp)
  var cg = newCgState(x=v2, b=v1)
  cg.solve(oa, sp)
  echo sp.finalIterations

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

  v2 := 0
  cg.reset
  sp.maxits = 0
  while cg.r2 > cg.r2stop:
    sp.maxits += 10
    cg.solve(oa, sp)
    let c = cg.x.norm2
    echo cg.iterations, ": ", c
