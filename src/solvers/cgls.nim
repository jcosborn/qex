import base
import layout
import field
import solverBase
export solverBase

type
  CglsState*[T1,T2] = object
    r1,p,x,b1: T1
    r2,Ap,b2: T2
    shift*,b1sq,b2sq,r1sqold,r1sq,r1sqstop,r2sq,r2sqstop: float
    iterations: int

proc reset*(cs: var CglsState) =
  cs.b1sq = -1
  cs.b2sq = -1
  cs.r1sqold = 1.0
  cs.r1sq = 99.0
  cs.r1sqstop = 0.0
  cs.r2sq = 99.0
  cs.r2sqstop = 0.0
  cs.iterations = 0

proc newCglsState*[T1,T2](x,b1: T1; b2: T2): CglsState[T1,T2] =
  result.r1 = newOneOf(b1)
  result.p = newOneOf(x)
  result.x = x
  result.b1 = b1
  result.r2 = newOneOf(b2)
  result.Ap = newOneOf(b2)
  result.b2 = b2
  result.shift = 0.0
  result.reset

#[
template subset(body: auto) {.dirty.} =
  onNoSync(sp.subset):
    body
template mythreads(body: auto) {.dirty.} =
  threads:
    onNoSync(sp.subset):
      body
]#

proc setup*(cs0: var CglsState, op: auto) =
  var cs = cs0
  threads:
    cs.b1sq = cs.b1.norm2
    cs.b2sq = cs.b2.norm2

  threads:
    cs.x := 0
    cs.r2 := cs.b2
    threadBarrier()
    op.applyAdj(cs.r1, cs.r2)
    threadBarrier()
    cs.r1 += cs.b1
    cs.r1sq = cs.r1.norm2
    cs.p := 0

  cs.r1sqold = 1.0
  cs.r2sq = cs.b2sq
  cs0 = cs
  #[
  if b1sq == 0.0:
    mythreads:
      x := 0
      r := 0
    r1sq = 0.0
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
  ]#


# solves:  (A' A + s) x = b1 + A' b2
# with r2 = b2 - A x
# and  r1 = b1 + A' r2 - s x
#guess: b1 <- b1 - s x0;  b2 <- b2 - A x0
proc solve*(state: var CglsState; op: auto; sp: var SolverParams) =
  mixin apply
  tic()
  template verb(n: int; body: auto) =
    if sp.verbosity>=n: body

  if state.b1sq<0:  # first call
    state.setup(op)
  verb(1):
    echo("CGLS b1 norm2: ", state.b1sq)
    echo("CGLS b2 norm2: ", state.b2sq)

  var
    r1 = state.r1
    r2 = state.r2
    p = state.p
    Ap = state.Ap
    x = state.x
    b1 = state.b1
    b1sq = state.b1sq
    b2sq = state.b2sq
    r1sq = state.r1sq
    r2sq = state.r2sq
    itn = state.iterations
    r1sqo = state.r1sqold
    shift = state.shift

  let r1sqstop = sp.r2req * b1sq
  let r2sqstop = sp.r2req * b2sq
  state.r1sqstop = r1sqstop
  state.r2sqstop = r2sqstop
  let maxits = sp.maxits

  toc("cgls setup")
  while (r1sq>=r1sqstop and r2sq>=r2sqstop) and itn<maxits:
    verb(2):
      #echo(itn, ": ", r1sq/b1sq, "  ", r2sq/b2sq)
      echo(itn, ": ", r1sq, "\t", r2sq)
    inc itn
    let beta = r1sq / r1sqo
    r1sqo = r1sq
    threads:
      p := r1 + beta*p
      threadBarrier()
      op.apply(Ap, p)
      threadBarrier()
      var psq = 0.0
      if shift != 0.0:
        psq = p.norm2
      #sub=sub2;
      let Apsq = Ap.norm2
      #sub=sub1;
      let alpha = r1sq/(Apsq + shift*psq)
      x += alpha * p
      #sub=sub2
      r2 -= alpha * Ap
      #sub=sub1
      threadBarrier()
      op.applyAdj(r1, r2)
      threadBarrier()
      if shift != 0.0:
        r1 -= shift * x
      if b1sq > 0.0:
        r1 += b1

      let r1sq0 = r1.norm2
      if threadNum==0:
        r1sq = r1sq0
      if r2sqstop!=0.0:
        #sub=sub2
        let r2sq0 = r2.norm2
        #sub=sub1
        if threadNum==0:
          r2sq = r2sq0

      #if(cgls->verbose>1) printf0("%*s%-3i r1sq = %-12g  r2sq = %g\n", cgls->indent, "", itn, r1sq/b1sq, r2sq/b2sq);

      #toc("cgls iterations")
      #var fr2: float
      #op.apply(Ap, x)
      #subset:
      #  r := b - Ap
      #  fr2 = r.norm2
      #verb(1):
      #  echo iterations, " acc r2:", r2/b2
      #  echo iterations, " tru r2:", fr2/b2

  state.iterations = itn
  state.r1sqold = r1sqo
  state.r1sq = r1sq
  state.r2sq = r2sq
  verb(1):
    echo state.iterations, " acc r2: ", r2sq/b2sq
    #threads:
    #  op.apply(Ap, x)
    #  var fr2: float
    #  subset:
    #    fr2 = (b - Ap).norm2
    #  echo "   ", fr2/b2
  sp.iterations = state.iterations
  toc("cgls final")

proc solve*(state: var CglsState; x: Field; b: Field2; op: auto;
            sp: var SolverParams) =
  state.x = x
  state.b = b
  state.reset
  state.solve(op, sp)

proc cglsSolve*(x: Field; b: Field2; op: auto; sp: var SolverParams) =
  var cgls = newCglsState(x, b)
  cgls.solve x, b, op, sp

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
  var zero = lo.ColorVector()
  type opArgs = object
    m: type(m)
  var oa = opArgs(m: m)
  proc apply*(oa: opArgs; r: type(v1); x: type(v1)) =
    r := oa.m * x
    #mul(r, m, x)
  proc applyAdj*(oa: opArgs; r: type(v1); x: type(v1)) =
    r := oa.m.adj * x
    #mul(r, m, x)
  var sp: SolverParams
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
    zero := 0
    echo v1.norm2
    echo m.norm2
  template resid(r,b,x,oa: auto) =
    oa.apply(r, x)
    r := b - r

  #cglsSolve(v2, v1, oa, sp)
  var cgls = newCglsState(x=v2, b1=zero, b2=v1)
  cgls.solve(oa, sp)
  echo sp.finalIterations

  v2 := 0
  cgls.reset
  sp.maxits = 0
  sp.verbosity = 0
  var stop = 0.0
  while cgls.r2sq > stop:
    sp.maxits += 10
    cgls.solve(oa, sp)
    v3.resid(v1,v2,oa)
    let tr2 = v3.norm2
    echo cgls.iterations, " ", cgls.r1sq, "/", cgls.r1sqstop, " ", tr2
    echo cgls.iterations, " ", cgls.r2sq, "/", cgls.r2sqstop
    #cgls.r := v3
    #cgls.r2 = tr2
    stop = 1e-6*cgls.b2sq
  echo sp.finalIterations, " ", cgls.r1sq, "/", cgls.r1sqstop

  v2 := 0
  cgls.reset
  sp.maxits = 0
  while cgls.r2sq > stop:
    sp.maxits += 10
    cgls.solve(oa, sp)
    let c = cgls.x.norm2
    echo cgls.iterations, ": ", c
