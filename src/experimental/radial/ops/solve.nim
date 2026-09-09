## CG and Jegerlehner shifted multishift CG (hep-lat/9612014).
##
## Generic over the vector type V.  A type qualifies when it provides, with plain
## Euclidean products,
##   newLike(x), sameShape(x, y), zero(x), x := y, axpy(x, a: float, y),
##   axpby(x, a, y, b), redot(x, y), norm2(x).
## `Spin` (core/spinor), `Gauge` (ops/gaugeact) and `seq[float]` (ops/gaugeact,
## the site scalars of the gauge projection) qualify; this is the one CG in the
## tree.  The operator is a closure `op(dst, src)` applying the UNSHIFTED
## Hermitian positive-definite matrix A.  `cgmSolve` solves (A + shifts[j]) x_j = b
## for every j out of one Krylov space, rebased on the smallest shift.
##
## Reported `r2` is the RELATIVE squared residual |b - (A + s) x|^2 / |b|^2,
## always recomputed from x (never the recursive estimate), so it can be
## compared to `r2req` directly.  Acceptance uses `r2 <= 1.001*r2req`: the
## recursive and recomputed residuals differ at the roundoff floor and an exact
## comparison rejects solves that are converged for every practical purpose.
## The 0.1% margin is not enough when r2req approaches the operator's own
## roundoff floor (~(eps*cond)^2) -- see doc/06 WP-D.
##
## Both solvers always start from x = 0.  No chronological guess: MD
## reversibility depends on the solve being a pure function of (b, op).
##
## Scratch: every solver takes a bank `work: var seq[V]` that it grows on demand
## (`ensureWork`) and reuses across calls, so a caller that keeps the bank makes
## the solve allocation-free after the first call.  The overloads without a bank
## allocate a temporary one.  A bank belongs to one lattice shape.

type
  CgInfo* = object
    iters*: int
    r2*: float               ## |b - A x|^2 / |b|^2, recomputed
    converged*: bool
  MultiCgInfo* = object
    iters*: int              ## multishift recurrence iterations
    r2*: seq[float]          ## per shift, |b - (A+s_j) x_j|^2 / |b|^2, recomputed
    r2pre*: seq[float]       ## same, measured before refinement -- how far the recurrence lied
    converged*: bool
    refined*: int            ## shifts that needed a single-shift refinement pass
    refits*: int             ## CG iterations spent in those refinement passes

proc ensureWork*[V](w: var seq[V], n: int, like: V) =
  ## Grow the bank to at least n vectors shaped like `like`.
  mixin newLike
  let m = w.len
  if m < n:
    w.setLen n
    for i in m..<n: w[i] = newLike(like)

proc conform*[V](x: var V, like: V) =
  ## Give x the shape of `like` (contents are not preserved).
  mixin newLike, sameShape
  if not sameShape(x, like): x = newLike(like)

proc resid2[V](b: V, ax: var V, x: V, s: float): float =
  ## |b - (A + s) x|^2 given ax = A x; `ax` is consumed as scratch.
  mixin axpy, axpby, norm2
  axpby(ax, 1.0, b, -1.0)
  if s != 0.0: axpy(ax, -s, x)
  norm2(ax)

proc cgRun[V](x: var V, b: V, s, r2stop: float, maxits: int,
              op: proc(dst: var V, src: V), r, p, ap: var V): tuple[iters: int, r2: float] =
  ## CG for (A + s) x = b from x = 0.  All scratch is supplied by the caller.
  ## Returns the recursive residual; the caller recomputes the true one.
  ## A non-positive curvature p.Ap stops the iteration: it can only come from a
  ## semidefinite A fed a right-hand side with a kernel component.
  mixin zero, `:=`, axpy, axpby, redot, norm2
  x.zero
  r := b
  p := b
  var r2 = norm2(r)
  var its = 0
  while r2 > r2stop and its < maxits:
    op(ap, p)
    if s != 0.0: axpy(ap, s, p)
    let pap = redot(p, ap)
    if pap <= 0.0: break
    let a = r2/pap
    axpy(x, a, p)
    axpy(r, -a, ap)
    let r2n = norm2(r)
    axpby(p, 1.0, r, r2n/r2)
    r2 = r2n
    inc its
  (its, r2)

proc cgSolve*[V](x: var V, b: V, r2req: float, maxits: int,
                 op: proc(dst: var V, src: V), work: var seq[V]): CgInfo =
  ## A x = b from x = 0; `work` supplies three scratch vectors.
  mixin zero, norm2
  conform(x, b)
  let b2 = norm2(b)
  if b2 == 0.0:
    x.zero
    return CgInfo(iters: 0, r2: 0.0, converged: true)
  ensureWork(work, 3, b)
  let it = cgRun(x, b, 0.0, r2req*b2, maxits, op, work[0], work[1], work[2])
  op(work[0], x)
  let r2 = resid2(b, work[0], x, 0.0)/b2
  CgInfo(iters: it.iters, r2: r2, converged: r2 <= 1.001*r2req)

proc cgSolve*[V](x: var V, b: V, r2req: float, maxits: int,
                 op: proc(dst: var V, src: V)): CgInfo =
  var w: seq[V]
  cgSolve(x, b, r2req, maxits, op, w)

proc cgmSolve*[V](xs: var seq[V], b: V, shifts: openArray[float],
                  r2req: float, maxits: int,
                  op: proc(dst: var V, src: V), work: var seq[V]): MultiCgInfo =
  ## (A + shifts[j]) x_j = b for every j; `work` supplies 3 + shifts.len vectors.
  mixin zero, `:=`, axpy, axpby, redot, norm2
  let ns = shifts.len
  doAssert ns > 0
  doAssert shifts[0] > 0.0, "shifts must be positive"
  for j in 1..<ns:
    doAssert shifts[j] > shifts[j-1], "shifts must be ascending"

  if xs.len != ns: xs.setLen ns
  for j in 0..<ns:
    conform(xs[j], b)
    xs[j].zero

  result.r2 = newSeq[float](ns)
  result.r2pre = newSeq[float](ns)
  let b2 = norm2(b)
  if b2 == 0.0:
    result.converged = true
    return

  ensureWork(work, 3 + ns, b)
  template r: untyped = work[0]
  template ap: untyped = work[1]
  template w: untyped = work[2]
  template ps(j: int): untyped = work[3 + j]
  var zi = newSeq[float](ns)
  var zim1 = newSeq[float](ns)
  var live = newSeq[bool](ns)

  let s0 = shifts[0]
  let r2stop = r2req*b2
  for j in 0..<ns:
    ps(j) := b
    zi[j] = 1.0
    zim1[j] = 1.0
    live[j] = true
  r := b

  # Base system is (A + s0); shift j enters the auxiliary recurrence as
  # sg = shifts[j] - s0 >= 0, so |z_j| <= 1 and the base residual bounds all.
  var r2 = b2
  var aim1 = 1.0        # alpha_{i-1}; with beta_{-1} = 0 this gives z_1 = 1/(1+sg*alpha_0)
  var bim1 = 0.0        # beta_{i-1}
  var its = 0
  while r2 > r2stop and its < maxits:
    op(ap, ps(0))
    axpy(ap, s0, ps(0))
    let a = r2/redot(ps(0), ap)
    axpy(xs[0], a, ps(0))
    axpy(r, -a, ap)
    let r2n = norm2(r)
    let bt = r2n/r2
    for j in 1..<ns:
      if not live[j]: continue
      let sg = shifts[j] - s0
      let d = a*bim1*(zim1[j] - zi[j]) + zim1[j]*aim1*(1.0 + sg*a)
      let zn = zi[j]*zim1[j]*aim1/d
      let zr = zn/zi[j]
      axpy(xs[j], a*zr, ps(j))
      zim1[j] = zi[j]
      zi[j] = zn
      # Shifted residual is z_j r, so freeze once it meets the target.  Not
      # just an optimization: z_j falls off like prod 1/(1+sg*alpha) and
      # underflows to 0, after which the recurrence for z is 0/0.
      if zn*zn*r2n <= r2stop: live[j] = false
      else: axpby(ps(j), zn, r, bt*zr*zr)
    axpby(ps(0), 1.0, r, bt)
    aim1 = a
    bim1 = bt
    r2 = r2n
    inc its
  result.iters = its

  # The recurrence residual z_j^2 r2 is only an estimate for the shifted
  # systems.  Recompute the true one and re-solve whatever missed.
  let acc = 1.001*r2req
  result.converged = true
  for j in 0..<ns:
    op(w, xs[j])
    result.r2pre[j] = resid2(b, w, xs[j], shifts[j])/b2
    result.r2[j] = result.r2pre[j]
    # j = 0 is the seed system: the recurrence above IS plain CG on (A + s0),
    # step for step, so refining it would reproduce xs[0] bit for bit.  If it
    # missed, it is at the roundoff floor of the operator and nothing here helps.
    if j > 0 and result.r2[j] > acc:
      let it = cgRun(xs[j], b, shifts[j], r2stop, maxits, op, r, w, ap)
      op(w, xs[j])
      result.r2[j] = resid2(b, w, xs[j], shifts[j])/b2
      inc result.refined
      result.refits += it.iters
    if result.r2[j] > acc: result.converged = false

proc cgmSolve*[V](xs: var seq[V], b: V, shifts: openArray[float],
                  r2req: float, maxits: int,
                  op: proc(dst: var V, src: V)): MultiCgInfo =
  var w: seq[V]
  cgmSolve(xs, b, shifts, r2req, maxits, op, w)
