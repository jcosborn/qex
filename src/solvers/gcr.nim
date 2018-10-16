import base
import layout
import field
import solverBase
export solverBase
import physics/qcdTypes

type
  Cmplx = DComplex
  GcrVec*[T,U] = object
    level*: int
    vec*: T
    Avec*: U
    Avn*: float
    alpha*: Cmplx
    beta*: seq[Cmplx]
  GcrState*[T,U] = object
    nv*: int
    vecs*: seq[GcrVec[T,U]]
    p*,x*: T
    r*,Ap*,b*: U
    b2*,r2*,r2old*,r2stop*: float
    iterations*: int

proc reset*(gs: var GcrState) =
  ## reset GcrState, forgetting previous vectors
  gs.b2 = -1
  gs.iterations = 0
  gs.r2old = 1.0
  gs.r2stop = 0.0
  gs.nv = 0

#[
proc restart*(gs: var GcrState) =
  ## restart solve keeping previous vectors
  gs.b2 = -1
  gs.iterations = 0
  gs.r2old = 1.0
  gs.r2stop = 0.0
  gs.nv = 0
]#

proc newGcrState*[T,U](x: T; b: U): GcrState[T,U] =
  result.p = newOneOf(x)
  result.x = x
  result.r = newOneOf(b)
  result.Ap = newOneOf(b)
  result.b = b
  result.vecs.newSeq(0)
  result.reset

template `[]`(gs: GcrState, i: int): untyped = gs.vecs[i]
#template level(gs: GcrState, i: int): untyped = gs[i].level

proc combine*(gs: var GcrState, n: int) =
  var c = gs[n-1].alpha / gs[n].alpha
  gs[n].Avec += c * gs[n-1].Avec
  swap(gs[n-1].Avec, gs[n].Avec)
  gs[n-1].Avn = gs[n].Avn + c.norm2 * gs[n-1].Avn
  c += gs[n].beta[n-1]
  gs[n].vec += c * gs[n-1].vec
  swap(gs[n-1].vec, gs[n].vec)
  gs[n-1].alpha = gs[n].alpha
  for i in 0..<n:
    gs[n].beta[i] += c * gs[n-1].beta[i]
    gs[n-1].beta[i] = gs[n].beta[i]

proc addvec*(gs: var GcrState) =
  var nv = gs.nv
  #echo("[GCR] Starting addvec; nv: ", nv)
  # #[
  if nv > 31:
    var n = nv - 1
    while n > 0 and gs[n-1].level <= gs[n].level + 0:
      gs.combine(n)
      dec n
      inc gs[n].level
    nv = n + 1
    if n==0: gs[n].level = 1
  # ]#
  let nvalloc = gs.vecs.len
  if nv >= nvalloc:
    gs.vecs.setLen(nv+1)
    for i in nvalloc .. nv:
      gs[i].vec = newOneOf(gs.x)
      gs[i].Avec = newOneOf(gs.b)
      gs[i].beta.newSeq(i+1)
  #echo "nva: ", nvalloc, "  nv: ", nv
  gs[nv].level = 1
  gs.nv = nv+1

proc orth*(gs: var GcrState) =
  var n = gs.nv - 1
  for i in 0 ..< n:
    let d = dot(gs[i].Avec, gs[n].Avec)
    let z = d/(-gs[i].Avn)
    gs[n].beta[i] = z
    gs[n].Avec += z * gs[i].Avec

proc getx*(gs: GcrState) =
  let nv = gs.nv
  var b = newSeq[Cmplx](nv)
  gs.x += gs[nv-1].alpha * gs[nv-1].vec
  for i in 0 .. (nv-2):
    b[i] := gs[nv-1].alpha * gs[nv-1].beta[i]
  for i in countdown(nv-2,0):
    let c = gs[i].alpha + b[i]
    gs.x += c * gs[i].vec
    for j in 0..<i:
      b[j] += c * gs[i].beta[j]

#proc solve*(linop: ptr linop_t; inv_arg: ptr QOP_invert_arg_t;
#         res_arg: ptr QOP_resid_arg_t; `out`: ptr Vector;
#         `in`: ptr Vector; r: ptr Vector; subset: QDP_Subset): QOP_status_t =
# solves: A x = b
proc solve*(gs: var GcrState; opx: var any; sp: var SolverParams) =
  mixin apply
  tic()
  let vrb = sp.verbosity
  template verb(n:int; body:untyped):untyped =
    if vrb>=n: body
  let sub = sp.subset
  template subset(body:untyped):untyped =
    #onNoSync(sub):
      body
  template mythreads(body:untyped):untyped =
    threads:
      #onNoSync(sub):
        body

  let
    r = gs.r
    #p = gs.p
    Ap = gs.Ap
    x = gs.x
    b = gs.b
  var
    insq = gs.b2
    rsq = gs.r2
    op = opx

  if insq<0:  # first call
    mythreads:
      #echo b.norm2
      insq = b.norm2
    gs.b2 = insq
    verb(1):
      echo("input norm2: ", insq)
    if insq == 0.0:
      mythreads:
        x := 0
        r := 0
      rsq = 0.0
    else:
      threads:
        op.apply(Ap, x)
        subset:
          r := b - Ap
          #p := 0
          rsq = r.norm2
          verb(3):
            #echo("p2: ", p.norm2)
            echo("r2: ", rsq)

  var
    #rsq: float
    #insq: float
    rsqstop: float
    alpha: float
    iteration = 0
    total_iterations = 0
    max_iterations = sp.maxits
    max_restarts = 1 #sp.max_restarts
  #if max_restarts < 0: max_restarts = 5
  # Default output values unless reassigned
  #res_arg.final_rsq = 0
  #res_arg.final_rel = 0
  #res_arg.final_iter = 0
  #res_arg.final_restart = 0
  #insq = gs.b2
  rsqstop = sp.r2req * insq
  #VERB(LOW, "GCR2: rsqstop = %g\x0A", rsqstop)

  #gs.init()
  #V_eq_V(r, `in`, subset)
  #r := b
  #r_eq_norm2_V(addr(rsq), r, subset)
  #rsq = r.norm2
  #gs.r2 = rsq
  verb(1): echo("[GCR] Starting iterations")
  while rsq > rsqstop and total_iterations < max_iterations:
    inc(iteration)
    inc(total_iterations)
    #echo "begin addvec"
    gs.addvec
    #echo "end addvec"
    let nv = gs.nv - 1
    if nv == 1:
      gs[0].vec := gs[0].alpha * gs[0].vec
      op.apply(gs[0].Avec, gs[0].vec)
      gs[0].Avn = gs[0].Avec.norm2
      op.apply(Ap, x)
      r := b - Ap
      let ctmp = dot(gs[0].Avec, r)
      gs[0].alpha = ctmp / gs[0].Avn
      r -= gs[0].alpha * gs[0].Avec
      rsq = r.norm2
      verb(2):
        echo iteration, "  rsq: ", gs.r2, " -> ", rsq
      gs.r2 = rsq
    op.preconditioner(gs[nv].vec, gs)
    op.apply(gs[nv].Avec, gs[nv].vec)
    gs.orth
    gs[nv].Avn = gs[nv].Avec.norm2
    let ctmp = dot(gs[nv].Avec, r)
    gs[nv].alpha = ctmp / gs[nv].Avn
    #echo "pAAp: ", gs[nv].Avn
    #echo "pAr: ", ctmp
    #echo "alpha: ", gs[nv].alpha
    r -= gs[nv].alpha * gs[nv].Avec
    #rsq -= gs[nv].alpha.norm2 * gs[nv].Avn
    rsq = r.norm2
    gs.r2 = rsq
    #rsq = r.norm2
    #VERB(HI, "GCR2: iter %i rsq = %g rel = %g\x0A", total_iterations, rsq,
    #     relnorm2)
    verb(3):
      echo iteration, ": ", rsq/gs.b2
    when false:
      getx(r, subset)
      linop(`out`, r, subset)
      V_eq_V_minus_V(r, `in`, `out`, subset)
      r_eq_norm2_V(addr(rsq), r, subset)
      #VERB(HI, "GCR2: iter %i rsq = %g rel = %g\x0A", total_iterations, rsq,
      #     relnorm2)

  #echo "[GCR] end iterations"
  #VERB(LOW, "GCR2: done: iter %i rsq = %g rel = %g\x0A", total_iterations,
  #     rsq, relnorm2)
  gs.getx()
  #gs.fini()
  opx = op
  #res_arg.final_rsq = rsq div insq
  #res_arg.final_rel = relnorm2
  sp.finalIterations = iteration
  verb(1):
    echo "GCR: its: ", iteration, "  rsq: ", rsq
  #return QOP_SUCCESS

when isMainModule:
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
