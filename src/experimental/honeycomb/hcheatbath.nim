## Cabibbo-Marinari SU(N) heatbath and SU(2)-subgroup microcanonical
## overrelaxation for the triangle gauge action on the 16-cell honeycomb
## (task **M2**) -- the algorithm the paper used.
##
## Local action (pinned numerically in tests/thmc.nim against `hcActionDeriv`
## and against single-link finite differences of `hcAction`):
##
##   S(U)  =  (beta/2) sum_x sum_{i=1}^{32} (1 - Re Tr P_i / N)
##         =  const  -  (beta/2N) Re Tr( U_l Sigma_l^dag )  + (U_l-independent)
##
## per link, where `Sigma_l = sum_{k=1}^{8} V_k` is the raw sum of the link's
## 8 triangle staples in QEX's V-convention (P_k = U_l V_k^dag) -- exactly the
## output of `hcActionDeriv` called with `beta = 2N` (its prefactor is
## `beta/2N`).  Hence each link is distributed as
##
##   P(U_l)  ~  exp( c Re Tr(U_l Sigma_l^dag) ) ,     c = beta/(2N).
##
## Update schedule (the structural facts are verified programmatically in
## tests/thmc.nim by scanning `hcgeom.triPath`):
##
##   (a) every triangle has exactly ONE axis edge, so the staple of any axis
##       link (uA or uB) is a product of diagonal links only  ->  all 8 axis
##       fields can be updated in one pass with their staples held fixed;
##   (b) the two diagonal edges of a triangle always carry DIFFERENT delta
##       indices, so no two links of the same uD[delta] field share a
##       triangle  ->  each uD[delta] field can be updated in one pass.
##
##   one sweep  =  axis staples -> update uA[0..3], uB[0..3]
##                 -> refresh the uA shift trees
##                 -> for delta in 0..15: staple of uD[delta] -> update it.
##
## Per link the update is the standard Cabibbo-Marinari sweep over the
## N(N-1)/2 SU(2) subgroups (rows/cols (k,l), k<l).  With R = U Sigma^dag and
## W the 2x2 submatrix of R, only the "real quaternion" part of W enters:
##
##   Re Tr(a W) = 2 (a0 w0 - a1 w1 - a2 w2 - a3 w3)    for a in SU(2),
##   w0 = (Re W00 + Re W11)/2,   w1 = (Im W01 + Im W10)/2,
##   w2 = (Re W01 - Re W10)/2,   w3 = (Im W00 - Im W11)/2,
##
## so with kq = |w|, V = quat(w)/kq in SU(2), and b = a V:
## Re Tr(a W) = kq Re Tr(a V) = 2 kq b0.  The heatbath draws b from
##
##   P(b) db ~ sqrt(1-b0^2) exp(alpha b0) db0 dOmega ,   alpha = 2 c kq,
##
## (`sampleA0` below: Kennedy-Pendleton for alpha > 1, plain rejection for
## alpha <= 1; both sample the same exact distribution, the split is only for
## efficiency), sets a = b V^dag and updates U <- E(a) U, R <- E(a) R where
## E(a) embeds the SU(2) matrix at rows/cols (k,l).
##
## The overrelaxation move is the microcanonical reflection a = (V^dag)^2,
## which preserves Re Tr(a W) exactly (Re Tr(V^dag^2 W) = kq Re Tr(V^dag)
## = kq Re Tr(V) = Re Tr(W)) and is an involution, hence detailed-balanced
## at constant action.
##
## All staple algebra reuses the `HcActionWork` shifters of task A; the
## per-site SU(2) work is done in plain scalar arrays via the lane-proxy
## accessors (`f{i}`), with one independent per-site RNG stream (`r{i}`), so
## results are independent of the thread count.  Nothing allocates inside
## `threads:` blocks.

import std/math
import base, layout, field, maths, rng
import physics/qcdTypes
import gauge
import hcgeom, hclayout, hcgauge, hcaction

export hcaction

# ---------------------------------------------------------------------------
# SU(2) sampling
# ---------------------------------------------------------------------------

proc sampleA0*[R](r: var R, alpha: float): float =
  ## Draw a0 from P(a0) ~ sqrt(1-a0^2) exp(alpha*a0) on [-1,1], alpha >= 0.
  ## alpha <= 1: plain rejection (envelope uniform, accept prob
  ##   sqrt(1-x^2) e^{alpha(x-1)} <= 1);
  ## alpha > 1: Kennedy-Pendleton (Phys.Lett.B 156 (1985) 393):
  ##   a0 = 1-2 lam^2, lam^2 = -(ln r1 + cos^2(2 pi r2) ln r3)/(2 alpha),
  ##   accepted with probability sqrt(1-lam^2).
  ## Both branches are exact; the split is efficiency only.
  const tiny = 1e-300
  if alpha <= 1.0:
    while true:
      let x = 2.0*float(uniform(r)) - 1.0
      let u = float(uniform(r))
      if u*u <= (1.0 - x*x)*exp(2.0*alpha*(x - 1.0)):
        return x
  else:
    while true:
      let r1 = float(uniform(r))
      let r2 = float(uniform(r))
      let r3 = float(uniform(r))
      let cc = cos(2.0*PI*r2)
      let lam2 = -(ln(r1 + tiny) + cc*cc*ln(r3 + tiny))/(2.0*alpha)
      if lam2 > 1.0: continue
      let u = float(uniform(r))
      if u*u <= 1.0 - lam2:
        return 1.0 - 2.0*lam2

proc sampleSphere*[R](r: var R): array[3, float] =
  ## uniform direction on S^2
  let ct = 2.0*float(uniform(r)) - 1.0
  let st = sqrt(max(0.0, 1.0 - ct*ct))
  let phi = 2.0*PI*float(uniform(r))
  [st*cos(phi), st*sin(phi), ct]

type Quat* = array[4, float]  ## q0 + i (q1 s1 + q2 s2 + q3 s3), s_j = Pauli

proc qmul*(p, q: Quat): Quat =
  ## quaternion product (matrix product of the associated 2x2 matrices)
  [p[0]*q[0] - p[1]*q[1] - p[2]*q[2] - p[3]*q[3],
   p[0]*q[1] + q[0]*p[1] - (p[2]*q[3] - p[3]*q[2]),
   p[0]*q[2] + q[0]*p[2] - (p[3]*q[1] - p[1]*q[3]),
   p[0]*q[3] + q[0]*p[3] - (p[1]*q[2] - p[2]*q[1])]

proc qconj*(q: Quat): Quat = [q[0], -q[1], -q[2], -q[3]]

proc su2HeatbathQuat*[R](r: var R, w: Quat, c: float): Quat =
  ## One SU(2) heatbath draw: given the quaternion part `w` of the 2x2
  ## submatrix W (see module docs) and the coupling c, return `a` distributed
  ## as P(a) ~ exp(c Re Tr(a W)) on SU(2).
  let kq = sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2] + w[3]*w[3])
  var v: Quat = [1.0, 0.0, 0.0, 0.0]
  if kq > 1e-300:
    v = [w[0]/kq, w[1]/kq, w[2]/kq, w[3]/kq]
  let alpha = 2.0*c*kq
  let b0 = sampleA0(r, alpha)
  let rho = sqrt(max(0.0, 1.0 - b0*b0))
  let n = sampleSphere(r)
  let b: Quat = [b0, rho*n[0], rho*n[1], rho*n[2]]
  qmul(b, qconj(v))

proc su2OverrelaxQuat*(w: Quat): Quat =
  ## Microcanonical reflection a = (V^dag)^2 for the same local action.
  let k2 = w[0]*w[0] + w[1]*w[1] + w[2]*w[2] + w[3]*w[3]
  if k2 < 1e-300:
    return [1.0, 0.0, 0.0, 0.0]
  # (V^dag)^2 = (v0^2 - |vv|^2, -2 v0 vv) with v = w/kq
  [(w[0]*w[0] - (w[1]*w[1] + w[2]*w[2] + w[3]*w[3]))/k2,
   -2.0*w[0]*w[1]/k2, -2.0*w[0]*w[2]/k2, -2.0*w[0]*w[3]/k2]

# ---------------------------------------------------------------------------
# per-site Cabibbo-Marinari update on plain scalar matrices
# ---------------------------------------------------------------------------

template toF(x: untyped): float =
  ## QEX single-site element accessors return lane proxies, not raw floats
  block:
    var v: float
    v := x
    v

# The site kernels work on Mat = array[nc][nc][2] (re,im), filled from the
# lane proxies once per site and written back once.

template siteMatT(nc: untyped): untyped = array[nc, array[nc, array[2, float]]]

template loadMat(m, f, i, nc: untyped) =
  for a in 0..<nc:
    for b in 0..<nc:
      m[a][b][0] = toF f{i}[a, b].re
      m[a][b][1] = toF f{i}[a, b].im

template storeMat(f, i, m, nc: untyped) =
  for a in 0..<nc:
    for b in 0..<nc:
      f{i}[a, b].re := m[a][b][0]
      f{i}[a, b].im := m[a][b][1]

template mulAdjInto(rm, um, sm, nc: untyped) =
  ## rm = um * sm^dag
  for a in 0..<nc:
    for b in 0..<nc:
      var re = 0.0
      var im = 0.0
      for k in 0..<nc:
        re += um[a][k][0]*sm[b][k][0] + um[a][k][1]*sm[b][k][1]
        im += um[a][k][1]*sm[b][k][0] - um[a][k][0]*sm[b][k][1]
      rm[a][b][0] = re
      rm[a][b][1] = im

template quatOf(rm, k, l: untyped): Quat =
  ## real quaternion components of the (k,l) 2x2 submatrix of rm
  [0.5*(rm[k][k][0] + rm[l][l][0]),
   0.5*(rm[k][l][1] + rm[l][k][1]),
   0.5*(rm[k][l][0] - rm[l][k][0]),
   0.5*(rm[k][k][1] - rm[l][l][1])]

template applySu2(m, a, k, l, nc: untyped) =
  ## m <- E(a) m, E(a) the SU(2) matrix a embedded at rows (k,l):
  ##   A = [[a0+i a3, a2+i a1], [-a2+i a1, a0-i a3]]
  for b in 0..<nc:
    let xr = m[k][b][0]
    let xi = m[k][b][1]
    let yr = m[l][b][0]
    let yi = m[l][b][1]
    m[k][b][0] = a[0]*xr - a[3]*xi + a[2]*yr - a[1]*yi
    m[k][b][1] = a[0]*xi + a[3]*xr + a[2]*yi + a[1]*yr
    m[l][b][0] = -a[2]*xr - a[1]*xi + a[0]*yr + a[3]*yi
    m[l][b][1] = -a[2]*xi + a[1]*xr + a[0]*yi - a[3]*yr

# ---------------------------------------------------------------------------
# field-level sweeps
# ---------------------------------------------------------------------------

type
  HcHeatbath*[V: static[int], F, W] = ref object
    ## Staple work space + parameters.  `w` is a task-A `HcActionWork` (shift
    ## trees etc.); `stA`/`stB` hold the 8 axis staple sums, `stD` the staple
    ## sum of the uD field currently being updated.  Create once with
    ## `newHcHeatbath`; sweeps are then allocation free.
    beta*: float
    w*: W
    stA*, stB*: array[nDim, F]
    stD*: F

proc newHcHeatbath*[V: static[int], F](g: HcGauge[V, F], beta: float): auto =
  var w = newHcActionWork(g)
  var h = HcHeatbath[V, F, type(w)](beta: beta, w: w)
  for mu in 0..<nDim:
    h.stA[mu] = g.uA[0].newOneOf
    h.stB[mu] = g.uA[0].newOneOf
  h.stD = g.uA[0].newOneOf
  h

proc rebind[V: static[int], F, W](h: HcHeatbath[V, F, W], g: HcGauge[V, F]) =
  ## point the uA shift trees at this g (ref assignments; outside threads:)
  for mu in 0..<nDim:
    h.w.shA[mu].setSrc g.uA[mu]

proc axisStaples*[V: static[int], F, W](h: HcHeatbath[V, F, W],
                                        g: HcGauge[V, F]) =
  ## stA[mu], stB[mu] := raw staple sums Sigma of uA[mu], uB[mu].
  ## These involve only uD links (fact (a) of the module docs), so they stay
  ## valid while all 8 axis fields are updated.
  let w = h.w
  let stA = h.stA
  let stB = h.stB
  threads:
    for mu in 0..<nDim:
      stA[mu] := 0
      stB[mu] := 0
    for t in apexTris:
      let
        mu = t.mu
        delta = t.delta
        deltaP = t.deltaP
        db = delta xor 15
        dbp = deltaP xor 15
      # V(uA[mu])(y) = [uD[delta]^dag uD[deltaP]](y-delta)
      w.t := g.uD[delta].adj * g.uD[deltaP]
      block:
        var cur = w.t
        var lev = 0
        for b in 0..<nDim:
          if ((delta shr b) and 1) != 0:
            cur = w.sB[lev][b] ^* cur
            inc lev
        stA[mu] += cur
      # V(uB[mu])(z) = uD[db](z) . uD[dbp](z+e_mu)^dag
      let s2 = w.sF[mu] ^* g.uD[dbp]
      stB[mu] += g.uD[db] * s2.adj

proc dStaple*[V: static[int], F, W](h: HcHeatbath[V, F, W], g: HcGauge[V, F],
                                    delta0: int) =
  ## stD := raw staple sum Sigma of uD[delta0].  Requires the uA shift trees
  ## to be current (call `refreshTrees` after any uA change).  The staple
  ## involves uA, uB and uD[delta0 xor 2^mu] only -- never uD[delta0] itself
  ## (fact (b)).
  let w = h.w
  let stD = h.stD
  threads:
    stD := 0
    for mu in 0..<nDim:
      let dx = delta0 xor (1 shl mu)
      if ((delta0 shr mu) and 1) == 0:
        # as uD[delta] of apex-B (mu,delta0):  V = uD[deltaP](z) uA[mu](z+delta0)^dag
        stD += g.uD[dx] * w.shA[mu].f[delta0].adj
        # as uD[dbp] of apex-A (mu, delta0 xor 15 xor 2^mu ...):
        #   V(y) = uB[mu](y-e_mu)^dag uD[db](y-e_mu),  db = delta0 xor 2^mu
        w.t := g.uB[mu].adj * g.uD[dx]
        let rs = w.sB[0][mu] ^* w.t
        stD += rs
      else:
        # as uD[deltaP] of apex-B (mu, delta0 xor 2^mu):  V = uD[delta](z) uA[mu](z+delta)
        stD += g.uD[dx] * w.shA[mu].f[dx]
        # as uD[db] of apex-A:  V = uB[mu](z) uD[dbp](z+e_mu),  dbp = delta0 xor 2^mu
        let s2 = w.sF[mu] ^* g.uD[dx]
        stD += g.uB[mu] * s2

proc refreshTrees*[V: static[int], F, W](h: HcHeatbath[V, F, W]) =
  ## refresh the 16-way shift trees of the (possibly just updated) uA fields;
  ## call inside `threads:`
  for mu in 0..<nDim:
    h.w.shA[mu].run

# The two per-field site loops.  `u` is the link field being updated, `st`
# its staple-sum field, `c = beta/(2N)`, `r` the per-site RNG field.
# Inside `threads:`; the `sites` iterator gives every thread whole SIMD
# blocks, so the lane writes below do not race.

proc hbField[F, RF](u, st: F, c: float, r: RF) =
  const nc = u[0].nrows
  var um, sm, rm: siteMatT(nc)
  for i in u.sites:
    loadMat(um, u, i, nc)
    loadMat(sm, st, i, nc)
    mulAdjInto(rm, um, sm, nc)          # R = U Sigma^dag
    for k in 0..<(nc-1):
      for l in (k+1)..<nc:
        let w = quatOf(rm, k, l)
        let a = su2HeatbathQuat(r{i}, w, c)
        applySu2(um, a, k, l, nc)
        applySu2(rm, a, k, l, nc)
    storeMat(u, i, um, nc)

proc orField[F](u, st: F) =
  const nc = u[0].nrows
  var um, sm, rm: siteMatT(nc)
  for i in u.sites:
    loadMat(um, u, i, nc)
    loadMat(sm, st, i, nc)
    mulAdjInto(rm, um, sm, nc)
    for k in 0..<(nc-1):
      for l in (k+1)..<nc:
        let w = quatOf(rm, k, l)
        let a = su2OverrelaxQuat(w)
        applySu2(um, a, k, l, nc)
        applySu2(rm, a, k, l, nc)
    storeMat(u, i, um, nc)

# ---------------------------------------------------------------------------
# sweeps
# ---------------------------------------------------------------------------

template sweepImpl(h, g, updateCall: untyped) {.dirty.} =
  ## one full lattice sweep: axis phase, tree refresh, 16 uD phases.
  ## `updateCall` is expanded inside `threads:` once per link field, with
  ## `u` bound to the field and `st` to its staple sums.
  const nc = g.uA[0][0].nrows
  let c = h.beta/(2.0*float(nc))
  rebind(h, g)
  block:
    let stA = h.stA
    let stB = h.stB
    let stD = h.stD
    axisStaples(h, g)
    threads:
      threadBarrier()
      for hcMu in 0..<nDim:
        block:
          template u: untyped = g.uA[hcMu]
          template st: untyped = stA[hcMu]
          updateCall
        block:
          template u: untyped = g.uB[hcMu]
          template st: untyped = stB[hcMu]
          updateCall
      threadBarrier()
      refreshTrees(h)
    for hcDelta in 0..<nDiag:
      dStaple(h, g, hcDelta)
      threads:
        threadBarrier()
        block:
          template u: untyped = g.uD[hcDelta]
          template st: untyped = stD
          updateCall

proc hbSweep*[V: static[int], F, W, RF](h: HcHeatbath[V, F, W],
                                        g: HcGauge[V, F], r: RF) =
  ## one Cabibbo-Marinari heatbath sweep over all 24 link fields
  sweepImpl(h, g):
    hbField(u, st, c, r)

proc orSweep*[V: static[int], F, W](h: HcHeatbath[V, F, W],
                                    g: HcGauge[V, F]) =
  ## one microcanonical overrelaxation sweep over all 24 link fields
  sweepImpl(h, g):
    orField(u, st)

proc update*[V: static[int], F, W, RF](h: HcHeatbath[V, F, W],
                                       g: HcGauge[V, F], r: RF,
                                       norSweeps = 3) =
  ## one "update" = 1 heatbath sweep + norSweeps overrelaxation sweeps
  hbSweep(h, g, r)
  for k in 0..<norSweeps:
    orSweep(h, g)

when isMainModule:
  import std/[strformat, monotimes, times]
  qexInit()
  echo "hcheatbath smoke test + timing"
  let hl = newHcLayout([4, 4, 4, 4])
  var r = hl.lo.newRNGField(MRG32k3a, 246813579'u64)
  var g = newHcGauge(hl)
  var hb = newHcHeatbath(g, 6.0)
  var w = hb.w
  echo &"unit start: S = {hcAction(w, hb.beta, g):.6f}"
  for n in 1..10:
    hb.update(g, r, 2)
    let s = hcAction(w, hb.beta, g)
    let ts = 1.0 - 2.0*s/(hb.beta*float(nTriPerSite*hl.nSites))
    echo &"update {n}: S = {s:.4f}  triangleSum = {ts:.6f}"
  let d = g.checkSU
  echo &"checkSU avg {d.avg:.3e} max {d.max:.3e}"
  # timing on 8^4
  block:
    let hl8 = newHcLayout([8, 8, 8, 8])
    var r8 = hl8.lo.newRNGField(MRG32k3a, 111'u64)
    var g8 = newHcGauge(hl8)
    var hb8 = newHcHeatbath(g8, 6.0)
    hb8.update(g8, r8, 1)   # warm buffers
    let t0 = getMonoTime()
    for n in 1..5:
      hb8.hbSweep(g8, r8)
    let t1 = getMonoTime()
    for n in 1..5:
      hb8.orSweep(g8)
    let t2 = getMonoTime()
    echo &"8^4 cells: hbSweep {(t1-t0).inMicroseconds.float/5e6:.4f} s  ",
         &"orSweep {(t2-t1).inMicroseconds.float/5e6:.4f} s"
  qexFinalize()
