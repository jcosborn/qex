## Overlap operator on S^2 x R  (WP-F).
##
## Normative reference: doc/02-formulation.md section 4, doc/04-interfaces.md section 10.
##
##   X = D_lat - M          the RAW Wilson operator of (IV.1), plain matrix adjoints,
##                          M in raw D_lat units (settled: doc/04 section 10; M = 1 default)
##   D_ov = 1 + X (X^dag X)^{-1/2}                                   (IV.9)
##   D(mass) = D_ov + mass  (additive mass convention)
##
## (X^dag X)^{-1/2} is the frozen Zolotarev rational of ops/zolotarev.nim,
##   R(x) = cst + sum_j res_j/(x + pole_j),  x in [smin^2, smax^2],  poles in sigma^2
## units, applied through one multishift CG (ops/solve.nim cgmSolve) per operator
## apply.  The window is frozen at construction; `kernelWindow` only monitors it --
## out of window is the caller's hard stop, never a rebuild.
##
## Solve accounting is accumulated in `Ov.stats`; `stats.ok` goes false if any inner
## solve missed its tolerance (keep r2inner at least ~2 decades above the operator's
## roundoff floor (eps*cond)^2, doc/06 WP-D).  Solver initial guesses are always zero.
##
## Workspace: everything the apply path needs is preallocated in `newOv` (`work`
## slots, `xs`/`xt` multishift banks); this module allocates nothing per apply.
## `cgmSolve`/`cgSolve` allocate their own scratch per call (nshift+3 fields; a
## documented WP-D property of the stateless section-9 signatures).  `dst` must not
## alias `src` in any apply.  Fixed `work` slots:
##   0 applyH intermediate, 1 applyOv/applyOvAdj, 2 applyNormal,
##   3..5 ovGradient (z, X^dag left, X s_j / X t_j).

import std/[math, complex]
import ../core/lattice
import ../core/spinor
import eigens/linalgFuncs
import wilson
import zolotarev
import solve

export wilson, zolotarev, solve

type
  SolveStats* = object
    ## Cumulative counters; clear with `clearStats`.
    nx*: int          ## D_W-level applies: X and X^dag count 1 each, H counts 2
    nmulti*: int      ## multishift solves
    miters*: int      ## multishift recurrence iterations
    mrefits*: int     ## multishift refinement iterations
    ncg*: int         ## plain CG solves (solveNormal outer, kernelWindow inner)
    cgiters*: int     ## their iterations
    ok*: bool         ## AND of every inner-solve `converged` since the last clear

  Ov* = ref object
    l*: Lat
    m*: float                  ## domain-wall height M, raw D_lat units; paper: 1
    rat*: Rat                  ## frozen rational; window [rat.smin, rat.smax]
    work*: seq[Spin]           ## persistent scratch, fixed slots (module header)
    xs*: seq[Spin]             ## multishift solutions (H + pole_j)^{-1} b
    xt*: seq[Spin]             ## second bank: ovGradient's t_j
    r2inner*, r2outer*: float  ## multishift / outer-CG relative residual targets
    maxits*: int
    stats*: SolveStats
    cu: ptr Gauge              ## gauge field for the solver closures, set per call
    cmass: float               ## mass for the solveNormal closure
    hop: proc(dst: var Spin, src: Spin)   ## H = X^dag X at cu[]
    nop: proc(dst: var Spin, src: Spin)   ## D(cmass)^dag D(cmass) at cu[]

const nwork = 6

proc clearStats*(o: Ov) =
  o.stats = SolveStats(ok: true)

proc applyX*(o: Ov, dst: var Spin, src: Spin, u: Gauge) =
  ## dst = X src = (D_W - M) src, the raw operator of (IV.1).
  applyDw(o.l, dst, src, u, o.m)
  inc o.stats.nx

proc applyXAdj*(o: Ov, dst: var Spin, src: Spin, u: Gauge) =
  ## dst = X^dag src.
  applyDwAdj(o.l, dst, src, u, o.m)
  inc o.stats.nx

proc applyH*(o: Ov, dst: var Spin, src: Spin, u: Gauge) =
  ## dst = X^dag X src.  Uses work[0]; neither argument may be work[0].
  applyX(o, o.work[0], src, u)
  applyXAdj(o, dst, o.work[0], u)

proc msolve(o: Ov, xs: var seq[Spin], b: Spin) =
  ## (H + pole_j) xs_j = b for every Zolotarev pole out of one Krylov space.
  ## `o.cu` must already point at the gauge field.
  let mi = cgmSolve(xs, b, o.rat.pole, o.r2inner, o.maxits, o.hop)
  inc o.stats.nmulti
  o.stats.miters += mi.iters
  o.stats.mrefits += mi.refits
  if not mi.converged: o.stats.ok = false

proc applyOv*(o: Ov, dst: var Spin, src: Spin, u: Gauge, mass = 0.0) =
  ## dst = (D_ov + mass) src = (1 + mass) src + X R(H) src.  One multishift solve.
  o.cu = addr u
  msolve(o, o.xs, src)
  o.work[1] := src                    # z = R(H) src
  scale(o.work[1], o.rat.cst)
  for j in 0..<o.rat.npole: axpy(o.work[1], o.rat.res[j], o.xs[j])
  applyX(o, dst, o.work[1], u)
  axpy(dst, 1.0 + mass, src)

proc applyOvAdj*(o: Ov, dst: var Spin, src: Spin, u: Gauge, mass = 0.0) =
  ## dst = (D_ov + mass)^dag src = (1 + mass) src + R(H) X^dag src.
  ## One multishift solve; the exact adjoint of `applyOv` up to solve residuals,
  ## because R(H) is Hermitian and (X R(H))^dag = R(H) X^dag.
  o.cu = addr u
  applyXAdj(o, o.work[1], src, u)
  msolve(o, o.xs, o.work[1])
  dst := src
  scale(dst, 1.0 + mass)
  axpy(dst, o.rat.cst, o.work[1])
  for j in 0..<o.rat.npole: axpy(dst, o.rat.res[j], o.xs[j])

proc applyNormal*(o: Ov, dst: var Spin, src: Spin, u: Gauge, mass = 0.0) =
  ## dst = (D_ov + mass)^dag (D_ov + mass) src.  Two multishift solves.
  applyOv(o, o.work[2], src, u, mass)
  applyOvAdj(o, dst, o.work[2], u, mass)

proc solveNormal*(o: Ov, x: var Spin, b: Spin, u: Gauge, mass = 0.0): CgInfo =
  ## x = [(D_ov + mass)^dag (D_ov + mass)]^{-1} b, CG from x = 0 at r2outer.
  ## Every CG iteration costs two multishift solves (plus one extra applyNormal
  ## for cgSolve's recomputed true residual).
  o.cu = addr u
  o.cmass = mass
  result = cgSolve(x, b, o.r2outer, o.maxits, o.nop)
  inc o.stats.ncg
  o.stats.cgiters += result.iters
  if not result.converged: o.stats.ok = false

proc ovGradient*(o: Ov, f: var Gauge, left, right: Spin, u: Gauge,
                 scale = 1.0, add = false) =
  ## f_link (+)= scale * d[ 2 Re <left, D_ov right> ] / d theta_link.
  ##
  ## THE single pullback: the HMC force, the Hasenbusch frames, the conserved
  ## current and the Ward test all go through here -- never derive a second one.
  ##
  ## With R(H) = c0 + sum_j r_j G_j, G_j = (H + q_j)^{-1}:
  ##   delta D_ov = (delta X) R(H) + X delta R,
  ##   delta R = -sum_j r_j G_j (delta H) G_j,  delta H = (delta X)^dag X + X^dag delta X,
  ## and with z = R(H) right, s_j = G_j right, t_j = G_j X^dag left,
  ##   dF = 2 Re[ <left, dX z> - sum_j r_j ( <X s_j, dX t_j> + <X t_j, dX s_j> ) ],
  ## each term one dwPullback (delta X = delta D_W, M constant).  The mass term of
  ## D(m) = D_ov + m carries no link, so there is no mass parameter here.
  ## Cost: two multishift solves (s_j and t_j), 2 npole + 1 X applies and pullbacks.
  o.cu = addr u
  msolve(o, o.xs, right)                       # s_j
  o.work[3].zero                               # z = R(H) right
  axpy(o.work[3], o.rat.cst, right)
  for j in 0..<o.rat.npole: axpy(o.work[3], o.rat.res[j], o.xs[j])
  applyXAdj(o, o.work[4], left, u)
  msolve(o, o.xt, o.work[4])                   # t_j
  dwPullback(o.l, f, left, o.work[3], u, scale, add)
  for j in 0..<o.rat.npole:
    let rj = scale*o.rat.res[j]
    applyX(o, o.work[5], o.xs[j], u)           # X s_j
    dwPullback(o.l, f, o.work[5], o.xt[j], u, -rj, add = true)
    applyX(o, o.work[5], o.xt[j], u)           # X t_j
    dwPullback(o.l, f, o.work[5], o.xs[j], u, -rj, add = true)

proc kernelWindow*(o: Ov, u: Gauge, iters = 32):
    tuple[smin, smax, lo, hi: float; inside: bool] =
  ## Monitor of spec(X^dag X) against the frozen rational window.  sigma_max by
  ## power iteration on H, sigma_min by inverse iteration (a CG solve of H x = v at
  ## shift 0 and tolerance r2inner per step); both from a fixed seeded start and at
  ## most `iters` steps, stopping early once the eigenpair residual
  ## rho = |H v - rq v|/|v| falls below 1e-8 rq.  `smin`/`smax` are the Rayleigh
  ## quotients (sigma units); `lo`/`hi` are expanded by rho, so [lo, hi] brackets
  ## the true extremes once the iteration is anywhere near converged
  ## (rq <= lambda_max and rq >= lambda_min hold unconditionally).
  ## `inside` false means the frozen window is violated: the caller must STOP --
  ## never rebuild the rational mid-ensemble.  Allocates its own vectors: this is a
  ## monitor, not an inner-loop routine.
  o.cu = addr u
  let n = o.l.nsite
  var
    v = newSpin(n)
    hv = newSpin(n)
    x = newSpin(n)
    r: Threefry4x64
  r.seedIndep(20260821, 0)
  proc eigpair(v, hv: var Spin): tuple[rq, rho: float] =
    ## Normalizes v, forms hv = H v, returns the Rayleigh quotient and residual.
    scale(v, 1.0/sqrt(norm2(v)))
    applyH(o, hv, v, o.cu[])
    result.rq = redot(v, hv)
    var d2 = 0.0
    for i in 0..<n:
      for c in 0..1:
        let
          dr = hv[i][c].re - result.rq*v[i][c].re
          di = hv[i][c].im - result.rq*v[i][c].im
        d2 += dr*dr + di*di
    result.rho = sqrt(d2)
  # power iteration for lambda_max
  v.gaussian r
  var top = eigpair(v, hv)
  for k in 1..<iters:
    if top.rho <= 1e-8*top.rq: break
    v := hv
    top = eigpair(v, hv)
  # inverse iteration for lambda_min
  v.gaussian r
  var bot = eigpair(v, hv)
  for k in 1..<iters:
    if bot.rho <= 1e-8*bot.rq: break
    let ci = cgSolve(x, v, o.r2inner, o.maxits, o.hop)
    inc o.stats.ncg
    o.stats.cgiters += ci.iters
    if not ci.converged: o.stats.ok = false
    v := x
    bot = eigpair(v, hv)
  result.smin = sqrt(bot.rq)
  result.smax = sqrt(top.rq)
  result.lo = sqrt(max(0.0, bot.rq - bot.rho))
  result.hi = sqrt(top.rq + top.rho)
  result.inside = result.lo >= o.rat.smin and result.hi <= o.rat.smax

proc denseOv*(o: Ov, u: Gauge): seq[Complex64] =
  ## Exact dense D_ov = 1 + X (X^dag X)^{-1/2}: X from denseDw, X^dag X
  ## eigendecomposed with zheev (zeigs), the inverse square root formed exactly.
  ## Column-major, dimension 2*nsite, row index 2*sIdx + spin.  Tests only.
  let
    nd = 2*o.l.nsite
    x = denseDw(o.l, u, o.m)
  var h = newSeq[Complex64](nd*nd)       # h = X^dag X, exactly Hermitian
  for j in 0..<nd:
    for i in j..<nd:
      var s = complex64(0.0, 0.0)
      for k in 0..<nd: s += conjugate(x[k + nd*i])*x[k + nd*j]
      h[i + nd*j] = s
      h[j + nd*i] = conjugate(s)
  var ev = newSeq[float](nd)
  zeigs(cast[ptr float64](addr h[0]), addr ev[0], nd)   # h <- eigenvectors V
  # g = V diag(ev^{-1/2}) V^dag
  var g = newSeq[Complex64](nd*nd)
  for j in 0..<nd:
    for i in 0..<nd:
      var s = complex64(0.0, 0.0)
      for k in 0..<nd:
        s += (1.0/sqrt(ev[k]))*h[i + nd*k]*conjugate(h[j + nd*k])
      g[i + nd*j] = s
  # result = 1 + X g
  result = newSeq[Complex64](nd*nd)
  for j in 0..<nd:
    for i in 0..<nd:
      var s = complex64(0.0, 0.0)
      for k in 0..<nd: s += x[i + nd*k]*g[k + nd*j]
      result[i + nd*j] = s
    result[j + nd*j] += complex64(1.0, 0.0)

proc newOv*(l: Lat, m: float, rat: Rat, r2inner, r2outer: float, maxits: int): Ov =
  ## All workspace is allocated here; the apply path allocates nothing.
  let o = Ov(l: l, m: m, rat: rat, r2inner: r2inner, r2outer: r2outer,
             maxits: maxits)
  o.work = newSeq[Spin](nwork)
  for i in 0..<nwork: o.work[i] = newSpin(l.nsite)
  o.xs = newSeq[Spin](rat.npole)
  o.xt = newSeq[Spin](rat.npole)
  for j in 0..<rat.npole:
    o.xs[j] = newSpin(l.nsite)
    o.xt[j] = newSpin(l.nsite)
  o.clearStats
  o.hop = proc(dst: var Spin, src: Spin) = applyH(o, dst, src, o.cu[])
  o.nop = proc(dst: var Spin, src: Spin) = applyNormal(o, dst, src, o.cu[], o.cmass)
  o
