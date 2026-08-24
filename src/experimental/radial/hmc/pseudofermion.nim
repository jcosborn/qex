## Hasenbusch pseudofermions for the dynamical overlap (WP-H).
##
## Normative reference: doc/02-formulation.md section 4.2, doc/04-interfaces.md
## section 13.  N_f even; one pseudofermion pair per two flavors (copies = nf/2),
## each copy carrying a Hasenbusch ladder over the strictly increasing masses
## m_0 < m_1 < ... < m_K, m_0 = physical, with the additive convention
## D(m) = D_ov + m and Q_i = D(m_i)^dag D(m_i).  Frames per copy:
##
##   ratio frame i (i = 0..K-1):  S_i = phi^dag D_{i+1} Q_i^{-1} D_{i+1}^dag phi
##     heatbath  phi = D_{i+1} Q_{i+1}^{-1} D_i^dag xi   =>  S_i = |xi|^2 exactly
##     (covariance D_{i+1}^{-dag} Q_i D_{i+1}^{-1}, so the frame weight is
##      det Q_i / det Q_{i+1}; no non-Hermitian solve anywhere)
##   heaviest frame K:            S_K = phi^dag Q_K^{-1} phi
##     heatbath  phi = D_K^dag xi                        =>  S_K = |xi|^2 exactly
##
## and the product of the frame weights telescopes to det Q_0 = |det D(m_0)|^2,
## one flavor pair.  masses.len == 1 is the no-Hasenbusch case (heaviest only).
##
## Forces differentiate ALL theta dependence and go through `ovGradient` only
## (doc/04 section 10: THE single pullback; the additive mass carries no link):
##   heaviest:  dS = -2 Re<y, dD_ov eta>,          eta = Q_K^{-1} phi, y = D_K eta
##   ratio:     dS = +2 Re<phi, dD_ov eta> - 2 Re<y, dD_ov eta>,
##              chi = D_{i+1}^dag phi, eta = Q_i^{-1} chi, y = D_i eta.
##
## Two overlap operators share the lattice and M but not the rational order:
## `actOp` (order 31) does every heatbath and every Hamiltonian; `frcOp`
## (order 11) does MD forces only.  Orders are never mixed inside one gradient.
##
## Solve failures do not raise (the apply path is silent by WP-F contract);
## they land in actOp.stats.ok / frcOp.stats.ok -- callers check.

import ../core/lattice
import ../core/spinor
import ../ops/overlap

export overlap

type Pf* = ref object
  l*: Lat
  actOp*, frcOp*: Ov         ## same lattice and M, rational order 31 vs 11
  nf*: int                   ## even, >= 2; copies = nf div 2
  masses*: seq[float]        ## strictly increasing; masses[0] = physical
  phi*: seq[seq[Spin]]       ## [copy][frame]
  xi2*: seq[seq[float]]      ## |xi|^2 recorded by the last refresh, [copy][frame]
  wa, wb, wc: Spin           ## frame scratch (disjoint from Ov.work)

func ncopy*(p: Pf): int = p.nf div 2
func nframe*(p: Pf): int = p.masses.len

proc newPf*(l: Lat, actOp, frcOp: Ov, nf: int, masses: seq[float]): Pf =
  doAssert nf >= 2 and nf mod 2 == 0, "nf must be even and >= 2"
  doAssert masses.len >= 1, "need at least one mass"
  for i in 1..<masses.len:
    doAssert masses[i] > masses[i-1], "masses must be strictly increasing"
  doAssert actOp.l == l and frcOp.l == l, "operators must share the lattice"
  doAssert actOp.m == frcOp.m, "operators must share the kernel height M"
  let p = Pf(l: l, actOp: actOp, frcOp: frcOp, nf: nf, masses: masses)
  p.phi = newSeq[seq[Spin]](nf div 2)
  p.xi2 = newSeq[seq[float]](nf div 2)
  for c in 0..<nf div 2:
    p.phi[c] = newSeq[Spin](masses.len)
    p.xi2[c] = newSeq[float](masses.len)
    for i in 0..<masses.len: p.phi[c][i] = newSpin(l.nsite)
  p.wa = newSpin(l.nsite)
  p.wb = newSpin(l.nsite)
  p.wc = newSpin(l.nsite)
  p

proc refreshFrame*(p: Pf, u: Gauge, c, i: int, r: var Threefry4x64) =
  ## Heatbath for frame i of copy c from the stream `r`, with actOp.
  ## Records |xi|^2 in p.xi2[c][i]; by construction S_i(phi) == |xi|^2 up to
  ## solve residuals (the rational cancels between heatbath and action).
  let k = p.masses.len - 1
  p.wa.gaussian r
  p.xi2[c][i] = norm2(p.wa)
  if i == k:
    applyOvAdj(p.actOp, p.phi[c][i], p.wa, u, p.masses[k])
  else:
    applyOvAdj(p.actOp, p.wb, p.wa, u, p.masses[i])            # D_i^dag xi
    discard solveNormal(p.actOp, p.wc, p.wb, u, p.masses[i+1]) # Q_{i+1}^{-1} .
    applyOv(p.actOp, p.phi[c][i], p.wc, u, p.masses[i+1])      # D_{i+1} .

proc refresh*(p: Pf, u: Gauge, r: var Threefry4x64) =
  ## All frames of all copies drawn sequentially from one stream.  The HMC
  ## driver keys a fresh stream per (copy, frame) instead -- see trajectory.nim.
  for c in 0..<p.ncopy:
    for i in 0..<p.nframe:
      refreshFrame(p, u, c, i, r)

proc frameAction*(p: Pf, o: Ov, u: Gauge, c, i: int): float =
  ## S_i of copy c evaluated with operator `o` (actOp for the Hamiltonian;
  ## frcOp only for the finite-difference force test).
  let k = p.masses.len - 1
  if i == k:
    discard solveNormal(o, p.wa, p.phi[c][i], u, p.masses[k])
    redot(p.phi[c][i], p.wa)
  else:
    applyOvAdj(o, p.wb, p.phi[c][i], u, p.masses[i+1])         # chi
    discard solveNormal(o, p.wa, p.wb, u, p.masses[i])         # eta = Q_i^{-1} chi
    redot(p.wb, p.wa)

proc pfAction*(p: Pf, u: Gauge): float =
  ## Sum of every frame action, actOp.  Enters the Hamiltonian and the accept
  ## test; never evaluate it with the low-order rational.
  for c in 0..<p.ncopy:
    for i in 0..<p.nframe:
      result += frameAction(p, p.actOp, u, c, i)

proc frameForce*(p: Pf, o: Ov, f: var Gauge, u: Gauge, c, i: int, add = false) =
  ## f (+)= dS_i/dtheta for copy c with operator `o`, via ovGradient only.
  ## Unprojected: the kernel projection is the trajectory driver's job.
  let k = p.masses.len - 1
  if i == k:
    discard solveNormal(o, p.wa, p.phi[c][i], u, p.masses[k])  # eta
    applyOv(o, p.wb, p.wa, u, p.masses[k])                     # y = D_K eta
    ovGradient(o, f, p.wb, p.wa, u, -1.0, add)
  else:
    applyOvAdj(o, p.wb, p.phi[c][i], u, p.masses[i+1])         # chi
    discard solveNormal(o, p.wa, p.wb, u, p.masses[i])         # eta
    applyOv(o, p.wc, p.wa, u, p.masses[i])                     # y = D_i eta
    ovGradient(o, f, p.phi[c][i], p.wa, u, 1.0, add)
    ovGradient(o, f, p.wc, p.wa, u, -1.0, add = true)

proc pfForce*(p: Pf, f: var Gauge, u: Gauge, level: int) =
  ## MD force of Hasenbusch level `level` (1 = heaviest frame ... nframe =
  ## lightest; level 0 is the gauge force and lives in the trajectory driver),
  ## summed over copies, with frcOp.
  let i = p.masses.len - level
  for c in 0..<p.ncopy:
    frameForce(p, p.frcOp, f, u, c, i, add = c > 0)
