## Cabibbo-Marinari heatbath and SU(2) microcanonical overrelaxation for the
## triangle action.
##
## Per link, with Sigma_l the raw sum of its 8 triangle staples
## (actionDeriv at beta = 2N):
##   P(U_l) ~ exp(c Re Tr(U_l Sigma_l^dag)),   c = beta/(2N)
## For the SU(2) subgroup (k,l) of R = U Sigma^dag, with w the real-quaternion
## part of the 2x2 block (w0 = (Re R_kk + Re R_ll)/2, w1 = (Im R_kl + Im R_lk)/2,
## w2 = (Re R_kl - Re R_lk)/2, w3 = (Im R_kk - Im R_ll)/2), kq = |w|, V = w/kq:
##   Re Tr(a W) = kq Re Tr(a V) = 2 kq b0,   b = a V
## so b is drawn from P(b) ~ sqrt(1-b0^2) exp(alpha b0), alpha = 2 c kq
## (Kennedy-Pendleton for alpha > 1, plain rejection below), a = b V^dag, and
## U <- E(a) U, R <- E(a) R.  Overrelaxation: a = (V^dag)^2, which preserves
## Re Tr(a W) and is an involution.
##
## Sweep order: every triangle has one axis edge and two diagonals of
## different delta, so the 8 axis fields update together with fixed staples
## (built from uD only), then each uD[delta] field updates whole.  Per-site
## RNG streams keep results independent of the thread count.

import std/math
import base, layout, field, maths, rng
import physics/qcdTypes
import gauge
import hcgeom, hcgauge, hcaction

export hcaction

# ---------------------------------------------------------------------------
# SU(2) sampling
# ---------------------------------------------------------------------------

proc sampleA0*[R](r: var R, alpha: float): float =
  ## a0 ~ sqrt(1-a0^2) exp(alpha a0) on [-1,1], alpha >= 0
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
      let lam2 = -(ln(r1) + cc*cc*ln(r3))/(2.0*alpha)
      if lam2 > 1.0: continue
      let u = float(uniform(r))
      if u*u <= 1.0 - lam2:
        return 1.0 - 2.0*lam2

proc sampleSphere[R](r: var R): array[3, float] =
  let ct = 2.0*float(uniform(r)) - 1.0
  let st = sqrt(max(0.0, 1.0 - ct*ct))
  let phi = 2.0*PI*float(uniform(r))
  [st*cos(phi), st*sin(phi), ct]

type Quat* = array[4, float]  ## q0 + i (q1 s1 + q2 s2 + q3 s3), s_j Pauli

func qmul(p, q: Quat): Quat =
  [p[0]*q[0] - p[1]*q[1] - p[2]*q[2] - p[3]*q[3],
   p[0]*q[1] + q[0]*p[1] - (p[2]*q[3] - p[3]*q[2]),
   p[0]*q[2] + q[0]*p[2] - (p[3]*q[1] - p[1]*q[3]),
   p[0]*q[3] + q[0]*p[3] - (p[1]*q[2] - p[2]*q[1])]

proc su2Heatbath*[R](r: var R, w: Quat, c: float): Quat =
  ## a ~ exp(c Re Tr(a W)) on SU(2), w the quaternion part of W
  let kq = sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2] + w[3]*w[3])
  var v: Quat = [1.0, 0.0, 0.0, 0.0]
  if kq > 0.0:
    v = [w[0]/kq, w[1]/kq, w[2]/kq, w[3]/kq]
  let b0 = sampleA0(r, 2.0*c*kq)
  let rho = sqrt(max(0.0, 1.0 - b0*b0))
  let n = sampleSphere(r)
  qmul([b0, rho*n[0], rho*n[1], rho*n[2]], [v[0], -v[1], -v[2], -v[3]])

func su2Overrelax*(w: Quat): Quat =
  ## a = (V^dag)^2, V = w/|w|
  let k2 = w[0]*w[0] + w[1]*w[1] + w[2]*w[2] + w[3]*w[3]
  if k2 == 0.0:
    return [1.0, 0.0, 0.0, 0.0]
  [(w[0]*w[0] - (w[1]*w[1] + w[2]*w[2] + w[3]*w[3]))/k2,
   -2.0*w[0]*w[1]/k2, -2.0*w[0]*w[2]/k2, -2.0*w[0]*w[3]/k2]

# ---------------------------------------------------------------------------
# per-site Cabibbo-Marinari update on scalar matrices
# ---------------------------------------------------------------------------

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
  ## rm = um sm^dag
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
  [0.5*(rm[k][k][0] + rm[l][l][0]),
   0.5*(rm[k][l][1] + rm[l][k][1]),
   0.5*(rm[k][l][0] - rm[l][k][0]),
   0.5*(rm[k][k][1] - rm[l][l][1])]

template applySu2(m, a, k, l, nc: untyped) =
  ## m <- E(a) m,  A = [[a0+i a3, a2+i a1], [-a2+i a1, a0-i a3]] at rows (k,l)
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
  HcHeatbath*[F, W] = ref object
    beta*: float
    w*: W                     ## ActionWork (shift trees, scratch)
    stA*, stB*: array[nDim, F]  ## axis staple sums
    stD*: F                   ## staple sum of the uD field being updated

proc newHcHeatbath*[F, W](g: HcGauge[F], w: W, beta: float): auto =
  ## allocates the 9 staple fields; outside `threads:`
  var h = HcHeatbath[F, W](beta: beta, w: w)
  for mu in 0..<nDim:
    h.stA[mu] = g.uA[0].newOneOf
    h.stB[mu] = g.uA[0].newOneOf
  h.stD = g.uA[0].newOneOf
  h

proc rebind[F, W](h: HcHeatbath[F, W], g: HcGauge[F]) =
  for mu in 0..<nDim:
    h.w.shA[mu].setSrc g.uA[mu]

proc axisStaples*[F, W](h: HcHeatbath[F, W], g: HcGauge[F]) =
  ## stA[mu], stB[mu] := Sigma of uA[mu], uB[mu] (uD links only)
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
      w.t := g.uD[delta].adj * g.uD[deltaP]
      block:
        var cur = w.t
        var lev = 0
        for b in 0..<nDim:
          if ((delta shr b) and 1) != 0:
            cur = w.sB[lev][b] ^* cur
            inc lev
        stA[mu] += cur
      let s2 = w.sF[mu] ^* g.uD[dbp]
      stB[mu] += g.uD[db] * s2.adj

proc dStaple*[F, W](h: HcHeatbath[F, W], g: HcGauge[F], delta0: int) =
  ## stD := Sigma of uD[delta0] (uA, uB and uD[delta0 xor 2^mu] only);
  ## the uA shift trees must be current
  let w = h.w
  let stD = h.stD
  threads:
    stD := 0
    for mu in 0..<nDim:
      let dx = delta0 xor (1 shl mu)
      if ((delta0 shr mu) and 1) == 0:
        stD += g.uD[dx] * w.shA[mu].f[delta0].adj
        w.t := g.uB[mu].adj * g.uD[dx]
        let rs = w.sB[0][mu] ^* w.t
        stD += rs
      else:
        stD += g.uD[dx] * w.shA[mu].f[dx]
        let s2 = w.sF[mu] ^* g.uD[dx]
        stD += g.uB[mu] * s2

proc refreshTrees*[F, W](h: HcHeatbath[F, W]) =
  ## inside `threads:`
  for mu in 0..<nDim:
    h.w.shA[mu].run

proc hbField[F, RF](u, st: F, c: float, r: RF) =
  const nc = u[0].nrows
  var um, sm, rm: siteMatT(nc)
  for i in u.sites:
    loadMat(um, u, i, nc)
    loadMat(sm, st, i, nc)
    mulAdjInto(rm, um, sm, nc)
    for k in 0..<(nc-1):
      for l in (k+1)..<nc:
        let a = su2Heatbath(r{i}, quatOf(rm, k, l), c)
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
        let a = su2Overrelax(quatOf(rm, k, l))
        applySu2(um, a, k, l, nc)
        applySu2(rm, a, k, l, nc)
    storeMat(u, i, um, nc)

template sweepImpl(h, g, updateCall: untyped) {.dirty.} =
  ## axis phase, tree refresh, 16 uD phases; `u` and `st` bound per field
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

proc hbSweep*[F, W, RF](h: HcHeatbath[F, W], g: HcGauge[F], r: RF) =
  ## one heatbath sweep over all 24 link fields
  sweepImpl(h, g):
    hbField(u, st, c, r)

proc orSweep*[F, W](h: HcHeatbath[F, W], g: HcGauge[F]) =
  ## one overrelaxation sweep over all 24 link fields
  sweepImpl(h, g):
    orField(u, st)

proc update*[F, W, RF](h: HcHeatbath[F, W], g: HcGauge[F], r: RF,
                       norSweeps = 3) =
  ## 1 heatbath sweep + norSweeps overrelaxation sweeps
  hbSweep(h, g, r)
  for k in 0..<norSweeps:
    orSweep(h, g)
