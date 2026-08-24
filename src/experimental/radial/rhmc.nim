## HMC for QED3 on S^2 x R: dynamical overlap fermions (Hasenbusch ladder) +
## non-compact Gaussian gauge action (WP-H production app).
##
## Theory (doc/02 sections 4.2, 5): S = S_g(theta) + pseudofermion frames, N_f
## even, one pair per two flavors, additive mass D(m) = D_ov + m, m_0 physical.
## actOp (Zolotarev order 31) does heatbath and accept/reject; frcOp (order 11)
## does MD forces only.  Nested Omelyan 2MN: gauge innermost with
## steps*innerSteps, each Hasenbusch frame with `steps`.  Gauge zero modes
## (gauge orbit + Polyakov flat mode) are projected out of the committed field,
## the refreshed momentum and every MD force (hmc/trajectory.nim).
##
## Couplings: parameter g2R = g^2 R with R = 1, so g2 = g2R; the deck labels
## runs by g^2 a with a = 1/L, reported here as g2a = g2R/lev.
##
## Randomness is trajectory addressed (baseSeed, trajectory number, purpose,
## copy, frame); restart depends only on the committed trajectory counter.
## `ckpt` is written every ckptFreq trajectories and, when it already exists at
## startup, the run RESUMES from it (every manifest field is validated).
## Configurations for the measurement app are the same format, written to
## `<cfg>.t<traj>` every measEvery.
##
## The kernel window is monitored every windowEvery trajectories;
## `inside == false` is a hard stop (the frozen rationals are never rebuilt).
##
## Mass convention in every manifest: additive, D(m) = D_ov + m.

import base
import std/[math, os, strformat, times]
import hmc/trajectory

qexInit()
freezeTimers()

letParam:
  lev = 1              ## icosahedral refinement L
  nt = 60              ## time slices
  at = 0.2             ## temporal spacing (units of R)
  g2R = 1.0            ## g^2 R; R = 1 so g2 = g2R (deck label: g2a = g2R/lev)
  gconv = 1            ## beta_l area convention: 0 geodesic, 1 exact (paper), 2 flat
  M = 1.0              ## domain-wall height, raw D_lat units
  ratLo = 0.3          ## frozen rational window, sigma units.  The FREE window
                       ## at L=1, nt=60, at=0.2 is [1.1234, 10.74] (WP-F), but
                       ## the interacting ensemble moves sigma_min: measured
                       ## ~0.58 on a g2R=1 heatbath start, hence the margin.
  ratHi = 12.5
  actOrder = 31        ## Zolotarev order, action/heatbath operator
  frcOrder = 11        ## Zolotarev order, MD force operator
  nf = 2               ## flavors, even
  masses = @[0.0, 0.5] ## Hasenbusch ladder, strictly increasing, masses[0] physical
  tau = 1.0            ## trajectory length
  steps = 4            ## 2MN steps per trajectory, every fermion frame
  innerSteps = 5       ## gauge steps = steps*innerSteps (innermost level)
  ntraj = 30           ## trajectories this run
  warmup = 5           ## forced-accept trajectories (absolute trajectory index)
  seed = 20260821      ## base seed of the trajectory-addressed RNG
  ckpt = ""            ## checkpoint path; resumes from it if it exists
  ckptFreq = 10        ## checkpoint every this many trajectories (0 = never)
  cfg = ""             ## configuration path stem for the measurement app
  measEvery = 0        ## save <cfg>.t<traj> every this many trajectories
  windowEvery = 10     ## kernelWindow monitor cadence (0 = never)
  r2inner = 1e-22      ## multishift relative residual target
  r2outer = 1e-18      ## outer normal-CG relative residual target
  maxits = 100000
  hotStart = 1         ## 1 = exact gauge heatbath start, 0 = cold (free field)

installStandardParams()
echoParams()
processHelpParam()

let
  sph = newSphere(lev)
  lat = newLat(sph, nt, at)
  bt = newBeta(lat, g2R, GeomConv(gconv))
  ratAct = newRat(ratLo, ratHi, actOrder)
  ratFrc = newRat(ratLo, ratHi, frcOrder)
  actOp = newOv(lat, M, ratAct, r2inner, r2outer, maxits)
  frcOp = newOv(lat, M, ratFrc, r2inner, r2outer, maxits)
  pf = newPf(lat, actOp, frcOp, nf, masses)
  sgDof = float(sph.ne*nt)          ## rank M; free theory: <S_g>/dof = 0.5

var lsteps = @[steps*innerSteps]
for i in 0..<masses.len: lsteps.add steps

var m = newRadialHmc(lat, bt, pf, tau, lsteps, uint64(seed))

echo &"lattice: L = {lev} (nv {sph.nv} ne {sph.ne} nf {sph.nf})  nt = {nt}" &
     &"  at = {at}  abar = {sph.abar:.6f}  abar/at = {lat.asOverAt:.3f}"
echo &"couplings: g2 = {g2R} (g2a = {g2R/float(lev)})  conv = {bt.conv}" &
     &"  maxM = {lat.maxM:.4f}"
echo &"overlap: M = {M}  window [{ratLo}, {ratHi}]" &
     &"  maxRelErr({actOrder}) = {ratAct.maxRelErr:.3e}" &
     &"  maxRelErr({frcOrder}) = {ratFrc.maxRelErr:.3e}"
echo &"rational hashes: act {ratAct.hash:#x}  frc {ratFrc.hash:#x}"
echo &"pseudofermions: nf = {nf} ({pf.ncopy} pair(s))  masses = {masses}" &
     &"  (additive: D(m) = D_ov + m)"
echo &"MD: tau = {tau}  steps/level = {lsteps} (gauge first, then heaviest->lightest)"
echo &"dof: links {nlink(lat)}  ker M {sph.nv*nt}  momentum dof {nlink(lat) - sph.nv*nt}"

if ckpt != "" and fileExists(ckpt):
  loadCheckpoint(m, ckpt)
  echo &"resumed from {ckpt} at trajectory {m.traj}"
elif hotStart != 0:
  var r = keyedRng(m.seed, 0, rkMomentum)
  discard heatbath(lat, m.u, bt, r)
  discard projectKernel(lat, m.u)
  echo "hot start: exact gauge heatbath"
else:
  echo "cold start"

var
  win = windowCheck(m)
  nAcc = 0
  nMeas = 0
  sumE = 0.0
  sumE2 = 0.0
  sumT = 0.0
echo &"window at start: [{win.lo:.4f}, {win.hi:.4f}] inside {win.inside}"
if not win.inside:
  echo "STOP: kernel window violated before the first trajectory"
  qexExit(1)

for k in 1..ntraj:
  m.forceAccept = m.traj + 1 <= warmup
  let t0 = epochTime()
  m.update
  let secs = epochTime() - t0
  sumT += secs
  if not m.forceAccept:
    inc nMeas
    if m.accepted: inc nAcc
    sumE += m.expmDeltaH
    sumE2 += m.expmDeltaH*m.expmDeltaH
  if windowEvery > 0 and m.traj mod windowEvery == 0:
    win = windowCheck(m)
  let sg = gaugeAction(lat, m.u, bt)
  echo &"traj: {m.traj} dH: {m.deltaH:.6g} acc: {int(m.accepted)}" &
       &" pacc: {m.pAccept:.4f} plaq-like: {sg/sgDof:.6f}" &
       &" window: [{win.lo:.4f},{win.hi:.4f}] secs: {secs:.2f}"
  if not (actOp.stats.ok and frcOp.stats.ok):
    echo "WARNING: a solve missed its tolerance this run (stats.ok false)"
  if not win.inside:
    if ckpt != "": saveCheckpoint(m, ckpt)
    echo "STOP: kernel window violated -- widen the frozen window and restart"
    qexExit(1)
  if ckpt != "" and ckptFreq > 0 and m.traj mod ckptFreq == 0:
    saveCheckpoint(m, ckpt)
  if cfg != "" and measEvery > 0 and m.traj mod measEvery == 0:
    saveCheckpoint(m, cfg & ".t" & $m.traj)

if ckpt != "": saveCheckpoint(m, ckpt)

echo &"trajectories: {ntraj} (committed counter {m.traj})"
echo &"seconds/trajectory: mean {sumT/float(ntraj):.3f}"
if nMeas > 1:
  let
    acc = float(nAcc)/float(nMeas)
    eMean = sumE/float(nMeas)
    eErr = sqrt((sumE2/float(nMeas) - eMean*eMean)/float(nMeas - 1))
  echo &"acceptance (post-warmup): {acc:.4f} ({nAcc}/{nMeas})"
  echo &"<exp(-dH)> = {eMean:.6f} +- {eErr:.6f}"
echo &"solves: act nmulti {actOp.stats.nmulti} miters {actOp.stats.miters}" &
     &" ncg {actOp.stats.ncg} cgiters {actOp.stats.cgiters} ok {actOp.stats.ok}"
echo &"solves: frc nmulti {frcOp.stats.nmulti} miters {frcOp.stats.miters}" &
     &" ncg {frcOp.stats.ncg} cgiters {frcOp.stats.cgiters} ok {frcOp.stats.ok}"
echo &"force counts (gauge, heaviest..lightest): {m.fcount}"

processSaveParams()
writeParamFile()
qexFinalize()
