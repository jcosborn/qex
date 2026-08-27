#RUNCMD $RUN1
## WP-H acceptance tests: doc/04-interfaces.md section 15 items 7-12.
##  1. heatbath identity S == |xi|^2 per frame
##  2. dense pseudofermion-action oracle + Hasenbusch telescoping (1e-11)
##  3. per-frame force vs centered FD of the frcOp frame action (best < 1e-5),
##     plus the systematic actOp-vs-frcOp force discrepancy (reported + pinned)
##  4. Gaussian momentum measure: transverse, <|p|^2> = dof = rank M
##  5. reversibility: bitwise field restore on forced reject; round-trip drift
##  6. exact integer force counts per trajectory for a pinned schedule
##  7. |dH| ~ dt^2 for the nested 2MN (ratios in (2.5, 6.0) for steps 1, 2, 4)
##  8. checkpoint: bitwise round trip, exact restart, corruption/mismatch raise
##  9. short end-to-end run: acceptance, <exp(-dH)>, seconds/trajectory
## plus the allocation regression (live memory across GC_fullCollect).

import std/[math, complex, os, streams, strformat, times, unittest]
import base/alignedMem
import eigens/linalgFuncs
import ../hmc/trajectory

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))

# --- helpers -------------------------------------------------------------------

proc randGauge(l: Lat, sed: int, amp = 1.0): Gauge =
  result = newGauge(l)
  var r: Threefry4x64
  r.seedIndep(sed, 0)
  for i in 0..<result.s.len: result.s[i] = amp*r.gaussian
  for i in 0..<result.t.len: result.t[i] = amp*r.gaussian

proc denseHBounds(x: seq[Complex64], nd: int): tuple[smin, smax: float] =
  ## sigma bounds of X from the eigenvalues of the dense X^dag X.
  var h = newSeq[Complex64](nd*nd)
  for j in 0..<nd:
    for i in j..<nd:
      var s = complex64(0.0, 0.0)
      for k in 0..<nd: s += conjugate(x[k + nd*i])*x[k + nd*j]
      h[i + nd*j] = s
      h[j + nd*i] = conjugate(s)
  var ev = newSeq[float](nd)
  zeigs(cast[ptr float64](addr h[0]), addr ev[0], nd)
  (sqrt(ev[0]), sqrt(ev[nd-1]))

proc denseAdjApply(a: seq[Complex64], nd: int, dst: var Spin, src: Spin) =
  ## dst = A^dag src, column-major a, spinor flat index i = 2*site + comp.
  for j in 0..<nd:
    var s = complex64(0.0, 0.0)
    for i in 0..<nd: s += conjugate(a[i + nd*j])*src[i shr 1][i and 1]
    dst[j shr 1][j and 1] = s

proc qEig(dov: seq[Complex64], nd: int, mass: float):
    tuple[v: seq[Complex64], ev: seq[float]] =
  ## Eigendecomposition of Q(m) = D(m)^dag D(m),
  ## D(m) = (1-m/2)dov + m (rho=1).
  var d = dov
  let alpha = ovMassAlpha(mass)
  for i in 0..<d.len: d[i] = alpha*d[i]
  for j in 0..<nd: d[j + nd*j] += complex64(mass, 0.0)
  result.v = newSeq[Complex64](nd*nd)
  for j in 0..<nd:
    for i in j..<nd:
      var s = complex64(0.0, 0.0)
      for k in 0..<nd: s += conjugate(d[k + nd*i])*d[k + nd*j]
      result.v[i + nd*j] = s
      result.v[j + nd*i] = conjugate(s)
  result.ev = newSeq[float](nd)
  zeigs(cast[ptr float64](addr result.v[0]), addr result.ev[0], nd)

proc oracleS(v: seq[Complex64], ev: seq[float], nd: int, x: Spin): float =
  ## sum_k |q_k^dag x|^2 / lambda_k
  for k in 0..<nd:
    var s = complex64(0.0, 0.0)
    for i in 0..<nd: s += conjugate(v[i + nd*k])*x[i shr 1][i and 1]
    result += abs2(s)/ev[k]

func logdet(ev: seq[float]): float =
  for x in ev: result += ln(x)

proc bitwiseEq(a, b: Gauge): bool =
  result = true
  for i in 0..<a.s.len:
    if a.s[i] != b.s[i]: return false
  for i in 0..<a.t.len:
    if a.t[i] != b.t[i]: return false

func checkpointHash(s: string): uint64 =
  ## Test-side copy of the checkpoint FNV-1a, used to construct checksum-valid
  ## files with deliberately incompatible semantic headers.
  result = 0xcbf29ce484222325'u64
  for ch in s:
    result = (result xor uint64(ord(ch))) * 0x100000001b3'u64

proc rewriteCheckpointI32(src, dst: string, offset: int, value: int32) =
  var payload = readFile(src)
  payload.setLen(payload.len - 8)       # drop the stored checksum
  var body = newStringStream(payload)
  body.setPosition(offset)
  body.write value
  var packed = newStringStream()
  packed.write body.data
  packed.write checkpointHash(body.data)
  writeFile(dst, packed.data)

# --- fixtures ------------------------------------------------------------------
# L = 1: nv 12, ne 30, nf 20.  Operator fixture nt = 6, at = 0.4 (nsite 72,
# dense dim 144, abar/at = 2.77, M = 1 < maxM = 2.08); momentum fixture nt = 4,
# at = 0.35 to reuse WP-G's measured dim ker M = 48, rank M = 120.

const
  ## Solver targets sit ~2 decades above the operators' roundoff floors
  ## (WP-D guidance).  Measured here: the multishift true-residual floor is
  ## ~8e-27 and the outer normal-CG floor ~5e-23, so the WP-F defaults
  ## (1e-26, 1e-22) sit ON the floor for these fixtures and trip stats.ok.
  r2in = 1e-24
  r2out = 1e-20
  mxit = 20000
  mass2 = @[0.0, 0.5]
  mass3 = @[0.0, 0.3, 0.8]

let
  sph = newSphere(1)
  lat = newLat(sph, 6, 0.4)
  nd = 2*lat.nsite
  u0 = newGauge(lat)
  ur = randGauge(lat, 20260821, 0.25)

# tight window: dense sigma envelope of the two fixed test fields, padded 5%
let
  b0 = denseHBounds(denseDw(lat, u0, 1.0), nd)
  br = denseHBounds(denseDw(lat, ur, 1.0), nd)
  sminFix = min(b0.smin, br.smin)
  smaxFix = max(b0.smax, br.smax)
  ratA = newRat(0.95*sminFix, 1.05*smaxFix, 31)
  ratF = newRat(0.95*sminFix, 1.05*smaxFix, 11)

echo &"fixture window: sigma envelope [{sminFix:.6f}, {smaxFix:.6f}]" &
     &"  maxRelErr(31) = {ratA.maxRelErr:.3e}  maxRelErr(11) = {ratF.maxRelErr:.3e}"

proc newOps(): tuple[a, f: Ov] =
  (newOv(lat, 1.0, ratA, r2in, r2out, mxit),
   newOv(lat, 1.0, ratF, r2in, r2out, mxit))

# MD window: envelope over free field and a g2 = 1 heatbath sample, padded wide
# because the trajectories move the field.
let bt = newBeta(lat, 1.0, gcExactArea)
var uh = newGauge(lat)
block:
  var r: Threefry4x64
  r.seedIndep(775577, 0)
  discard heatbath(lat, uh, bt, r)
let
  bh = denseHBounds(denseDw(lat, uh, 1.0), nd)
  sminMd = 0.75*min(b0.smin, bh.smin)
  smaxMd = 1.30*max(b0.smax, bh.smax)
  ratAmd = newRat(sminMd, smaxMd, 31)
  ratFmd = newRat(sminMd, smaxMd, 11)

echo &"MD window: [{sminMd:.6f}, {smaxMd:.6f}] (heatbath smin {bh.smin:.6f})" &
     &"  maxRelErr(31) = {ratAmd.maxRelErr:.3e}  maxRelErr(11) = {ratFmd.maxRelErr:.3e}"

proc newMdHmc(masses: seq[float], steps: seq[int], tau: float,
              seed: uint64): RadialHmc =
  let
    a = newOv(lat, 1.0, ratAmd, r2in, r2out, mxit)
    f = newOv(lat, 1.0, ratFmd, r2in, r2out, mxit)
  newRadialHmc(lat, bt, newPf(lat, a, f, 2, masses), tau, steps, seed)

# --- 1. heatbath identity --------------------------------------------------------

suite "heatbath identity S == |xi|^2 (ladder item 7)":

  test "pseudofermion ladders reject masses outside the standard interval":
    let (a, f) = newOps()
    expect ValueError:
      discard newPf(lat, a, f, 2, @[-0.1, 0.5])
    expect ValueError:
      discard newPf(lat, a, f, 2, @[0.0, 2.0])

  test "every frame, both mass ladders, free and random field":
    var worst = 0.0
    for masses in [mass2, @[0.0], mass3]:
      let (a, f) = newOps()
      let p = newPf(lat, a, f, 2, masses)
      for (nm, u) in [("free", u0), ("random", ur)]:
        var r: Threefry4x64
        r.seedIndep(31001, 0)
        refresh(p, u, r)
        for c in 0..<p.ncopy:
          for i in 0..<p.nframe:
            let
              s = frameAction(p, p.actOp, u, c, i)
              e = abs(s - p.xi2[c][i])/p.xi2[c][i]
            worst = max(worst, e)
        check p.actOp.stats.ok
      echo &"  masses {masses}: worst |S - |xi|^2| / |xi|^2 = {worst:.3e}"
    check worst < 1e-8

# --- 2. dense oracle ------------------------------------------------------------

suite "dense action oracle and Hasenbusch telescoping (ladder item 8)":

  test "pf action == sum_k |q_k^dag phi|^2/lambda_k from denseOv, 3-mass ladder":
    let (a, f) = newOps()
    let p = newPf(lat, a, f, 2, mass3)
    var r: Threefry4x64
    r.seedIndep(32001, 0)
    refresh(p, ur, r)
    let dov = denseOv(a, ur)
    var qs: seq[tuple[v: seq[Complex64], ev: seq[float]]]
    for m in mass3: qs.add qEig(dov, nd, m)
    let k = mass3.len - 1
    var worst = 0.0
    var chi = newSpin(lat.nsite)
    var stot = 0.0
    var otot = 0.0
    for c in 0..<p.ncopy:
      for i in 0..<p.nframe:
        let s = frameAction(p, p.actOp, ur, c, i)
        var o: float
        if i == k:
          o = oracleS(qs[k].v, qs[k].ev, nd, p.phi[c][i])
        else:
          var d = dov
          let alpha = ovMassAlpha(mass3[i+1])
          for j in 0..<d.len: d[j] = alpha*d[j]
          for j in 0..<nd: d[j + nd*j] += complex64(mass3[i+1], 0.0)
          denseAdjApply(d, nd, chi, p.phi[c][i])
          o = oracleS(qs[i].v, qs[i].ev, nd, chi)
        worst = max(worst, abs(s - o)/o)
        stot += s
        otot += o
    echo &"  worst per-frame |S - S_dense|/S_dense = {worst:.3e}" &
         &"   total S = {stot:.12e} vs dense {otot:.12e}"
    check worst < 1e-8

    # telescoping: sum_i (logdet Q_i - logdet Q_{i+1}) + logdet Q_K == logdet Q_0
    var lhs = 0.0
    for i in 0..<k: lhs += logdet(qs[i].ev) - logdet(qs[i+1].ev)
    lhs += logdet(qs[k].ev)
    let rhs = logdet(qs[0].ev)
    echo &"  telescoping |lhs - logdet Q_0| = {abs(lhs - rhs):.3e}" &
         &"  (logdet Q_0 = {rhs:.12e})"
    check abs(lhs - rhs) < 1e-11

# --- 3. per-frame force vs finite differences -------------------------------------

suite "per-frame force vs centered FD (ladder item 9)":

  test "FD of the frcOp frame action; actOp-frcOp discrepancy":
    let (a, f) = newOps()
    let p = newPf(lat, a, f, 2, mass2)
    var r: Threefry4x64
    r.seedIndep(33001, 0)
    refresh(p, ur, r)
    let du = randGauge(lat, 33002)
    var fr = newGauge(lat)
    var fa = newGauge(lat)
    var up = newGauge(lat)
    var worstFd = 0.0
    var worstDisc = 0.0
    for i in 0..<p.nframe:
      frameForce(p, p.frcOp, fr, ur, 0, i)
      let pred = dot(fr, du)
      var best = 1e30
      var bh = 0.0
      for h in [1e-2, 1e-3, 1e-4, 1e-5]:
        up := ur
        axpy(up, h, du)
        let sp = frameAction(p, p.frcOp, up, 0, i)
        up := ur
        axpy(up, -h, du)
        let sm = frameAction(p, p.frcOp, up, 0, i)
        let e = abs((sp - sm)/(2.0*h) - pred)/abs(pred)
        if e < best:
          best = e
          bh = h
      frameForce(p, p.actOp, fa, ur, 0, i)
      var d2 = 0.0
      for j in 0..<fa.s.len:
        let d = fa.s[j] - fr.s[j]
        d2 += d*d
      for j in 0..<fa.t.len:
        let d = fa.t[j] - fr.t[j]
        d2 += d*d
      let disc = sqrt(d2/norm2(fa))
      echo &"  frame {i}: best FD rel err = {best:.3e} at h = {bh:.0e}" &
           &"   |f11 - f31|/|f31| = {disc:.3e}"
      worstFd = max(worstFd, best)
      worstDisc = max(worstDisc, disc)
    check p.frcOp.stats.ok and p.actOp.stats.ok
    check worstFd < 1e-5
    # The systematic order-11 vs order-31 force difference is O(maxRelErr(11))
    # times an amplification set by the frame's conditioning (worst for the
    # m = 0 ratio frame, measured ~1e2 here).  The pinned factor guards against
    # order-mixing regressions (those show up at ~1e-3, four decades away).
    echo &"  worst discrepancy = {worstDisc:.3e} = {worstDisc/ratF.maxRelErr:.1f} x" &
         &" maxRelErr(11) = {ratF.maxRelErr:.3e}"
    check worstDisc < 2e3*ratF.maxRelErr

# --- 4. momentum measure ----------------------------------------------------------

suite "Gaussian momentum measure (dof = rank M)":

  test "projected momentum transverse; <|p|^2> = (ne+nv)nt - (nv nt - 1) - 1":
    let lat4 = newLat(sph, 4, 0.35)
    let
      nl = nlink(lat4)
      kerdim = lat4.sph.nv*lat4.nt        # gauge orbit (nv nt - 1) + Polyakov 1
      dof = nl - kerdim
    check nl == 168 and kerdim == 48 and dof == 120   # WP-G measured values
    var p = newGauge(lat4)
    var s = 0.0
    var s2 = 0.0
    const n = 1000
    for k in 0..<n:
      var r = keyedRng(424242'u64, k+1, rkMomentum)
      gaussian(p, r)
      discard projectKernel(lat4, p)
      let p2 = norm2(p)
      s += p2
      s2 += p2*p2
      if k == 0:
        let t = transversality(lat4, p)
        echo &"  first draw: |div p|^2/|p|^2 = {t.divP:.3e}  flat = {t.flatP:.3e}"
        check t.divP < 1e-20
        check t.flatP < 1e-12
    let
      mean = s/float(n)
      err = sqrt((s2/float(n) - mean*mean)/float(n-1))
    echo &"  <|p|^2> = {mean:.4f} +- {err:.4f}  vs dof = {dof}" &
         &"  (pull {abs(mean - float(dof))/err:.2f})"
    check abs(mean - float(dof)) < 5.0*err

  test "pure-gauge sampler (single level, pf = nil)":
    let lat4 = newLat(sph, 4, 0.35)
    var m = newRadialHmc(lat4, newBeta(lat4, 1.7, gcExactArea), nil, 0.5,
                         @[16], 111222'u64)
    for k in 0..<4: m.update
    let t = transversality(lat4, m.p)
    echo &"  4 traj: dH = {m.deltaH: .3e}  acc {m.nAccepts}/{m.nUpdates}" &
         &"  fcount {m.fcount}  divP {t.divP:.3e}"
    check m.fcount == @[4*2*16]
    check abs(m.deltaH) < 1e-3       # quadratic action, 16-step 2MN
    check t.divP < 1e-18

# --- 5-6. reversibility, bitwise reject, force counts ------------------------------

suite "reversibility and force counts (ladder items 10, 11)":

  test "forced reject restores the field bitwise; round-trip drift; counts":
    var m = newMdHmc(mass2, @[4, 2, 2], 0.6, 550001'u64)
    # move off the free field: two forced-accept trajectories
    m.forceAccept = true
    m.update
    m.update
    m.forceAccept = false
    echo &"  thermalized 2 traj: dH = {m.deltaH:.3e}  S_g = {gaugeAction(lat, m.u, bt):.4f}"

    # 6. exact integer force counts for the schedule [4, 2, 2] (2MN: 2 per step)
    m.clearForceCounts
    m.forceAccept = true
    m.update
    m.forceAccept = false
    echo &"  force counts per trajectory, steps [4,2,2]: {m.fcount}"
    check m.fcount == @[8, 4, 4]

    # 5a. forced reject restores bitwise
    var uSnap = newGauge(lat)
    uSnap := m.u
    m.forceReject = true
    m.update
    m.forceReject = false
    check not m.accepted
    check bitwiseEq(m.u, uSnap)
    echo "  forced reject: field restored bitwise"

    # 5b. round-trip reversibility on the current (u, p, phi)
    let rv = reversibilityCheck(m)
    echo &"  reverse: |du| rms {rv.duRms:.3e} max {rv.duMax:.3e}" &
         &"  |dp| rms {rv.dpRms:.3e} max {rv.dpMax:.3e}  dH {rv.dh:.3e}"
    echo &"  round-trip momentum: |div p|^2/|p|^2 = {rv.divP:.3e}  flat = {rv.flatP:.3e}"
    check rv.duRms < 1e-6 and rv.dpRms < 1e-6
    check abs(rv.dh) < 1e-6
    check rv.divP < 1e-18
    check rv.flatP < 1e-10

    # final momentum of a normal trajectory is still transverse
    m.update
    let t = transversality(lat, m.p)
    echo &"  post-MD momentum: |div p|^2/|p|^2 = {t.divP:.3e}  flat = {t.flatP:.3e}"
    check t.divP < 1e-18
    check t.flatP < 1e-10
    check m.pf.actOp.stats.ok and m.pf.frcOp.stats.ok

# --- 7. dH scaling -----------------------------------------------------------------

suite "|dH| ~ dt^2 for the nested 2MN (ladder item 12)":

  test "ratios between step counts 1, 2, 4 in (2.5, 6.0)":
    ## The MD integrates the exact gradient of the ORDER-11 Hamiltonian, so
    ## dH11 scales as dt^2.  The accept/reject dH31 equals dH11 plus the
    ## dt-independent action mismatch [S31-S11](u_f) - [S31-S11](u_i), whose
    ## scale is ~2 maxRelErr(11) S_pf; both series are measured, the ratio
    ## assertion is on dH11 and dH31's saturation is reported against that
    ## floor (this is why the accept dH stops shrinking with the step count,
    ## and why Metropolis with the order-31 action is what corrects the
    ## order-11 flow).
    ## Base schedule [4, 2, 2]: the coarser [2, 1, 1] has gauge dt = 0.4 with
    ## dt sqrt(lambda_max(M)) ~ 1.2, outside the asymptotic regime (same
    ## effect as WP-G's flow-order test at n = 2), and its dH is non-monotone.
    var m = newMdHmc(mass2, @[4, 2, 2], 0.8, 660001'u64)
    proc h11(m: RadialHmc): float =
      result = 0.5*norm2(m.p) + gaugeAction(lat, m.u, bt)
      for c in 0..<m.pf.ncopy:
        for i in 0..<m.pf.nframe:
          result += frameAction(m.pf, m.pf.frcOp, m.u, c, i)
    m.forceAccept = true
    m.update
    m.update
    m.forceAccept = false
    var uSnap = newGauge(lat)
    uSnap := m.u
    var dh11: array[3, float]
    var dh31: array[3, float]
    var spf = 0.0
    for k in 0..<3:
      let sc = 1 shl k
      m.setSteps(@[4*sc, 2*sc, 2*sc])
      m.u := uSnap
      refreshMomentum(m, 999)
      refreshPseudo(m, 999)
      let
        ha = h11(m)
        hb = hmcH(m)
      if k == 0: spf = pfAction(m.pf, m.u)
      mdEvolve(m)
      dh11[k] = h11(m) - ha
      dh31[k] = hmcH(m) - hb
      echo &"  steps x{sc}: dH11 = {dh11[k]: .6e}   dH31 = {dh31[k]: .6e}"
    m.u := uSnap
    let
      r01 = abs(dh11[0]/dh11[1])
      r12 = abs(dh11[1]/dh11[2])
      floor31 = 2.0*ratFmd.maxRelErr*spf
    echo &"  dH11 ratios: {r01:.3f}, {r12:.3f}  (2MN: expect ~4)"
    echo &"  dH31 - dH11 = {dh31[0]-dh11[0]: .3e} {dh31[1]-dh11[1]: .3e}" &
         &" {dh31[2]-dh11[2]: .3e}  vs 2 maxRelErr(11) S_pf = {floor31:.3e}"
    check r01 > 2.5 and r01 < 6.0
    check r12 > 2.5 and r12 < 6.0
    check abs(dh31[2] - dh11[2]) < 10.0*floor31
    check m.pf.actOp.stats.ok and m.pf.frcOp.stats.ok

# --- 8. checkpoint -----------------------------------------------------------------

suite "checkpoint round trip and exact restart":

  test "bitwise field, identical restart draws, corruption/mismatch raise":
    let
      dir = getTempDir()
      path = dir / "thmc_ckpt.bin"
      pathBad = dir / "thmc_ckpt_bad.bin"
      pathV1 = dir / "thmc_ckpt_v1.bin"
      pathConv = dir / "thmc_ckpt_convention.bin"
    var ma = newMdHmc(mass2, @[4, 2, 2], 0.6, 770001'u64)
    ma.update                    # trajectory 1
    saveCheckpoint(ma, path)
    var uSave = newGauge(lat)
    uSave := ma.u
    let trajSave = ma.traj
    # chain A continues: trajectories 2 and 3
    var recRnd: seq[float]
    var recDh: seq[float]
    var recAcc: seq[bool]
    for k in 0..<2:
      ma.update
      recRnd.add ma.rnd
      recDh.add ma.deltaH
      recAcc.add ma.accepted
    # chain B: fresh object, same parameters, restart from the checkpoint
    var mb = newMdHmc(mass2, @[4, 2, 2], 0.6, 770001'u64)
    loadCheckpoint(mb, path)
    check mb.traj == trajSave
    check bitwiseEq(mb.u, uSave)
    var worstDh = 0.0
    for k in 0..<2:
      mb.update
      check mb.rnd == recRnd[k]
      check mb.accepted == recAcc[k]
      worstDh = max(worstDh, abs(mb.deltaH - recDh[k]))
    echo &"  restart: draws bitwise, worst |dH_restart - dH_original| = {worstDh:.3e}"
    check worstDh < 1e-10
    check bitwiseEq(mb.u, ma.u)
    echo "  restarted chain's committed field is bitwise identical"

    # corruption: flip one byte in the field payload
    var data = readFile(path)
    let mid = data.len div 2
    data[mid] = chr(ord(data[mid]) xor 1)
    writeFile(pathBad, data)
    expect ValueError:
      loadCheckpoint(mb, pathBad)
    # Semantic compatibility: both files have valid checksums, but one carries
    # the retired v1/additive header and one has the wrong v2 convention id.
    rewriteCheckpointI32(path, pathV1, 8, 1'i32)
    var mismatch = ""
    try:
      loadCheckpoint(mb, pathV1)
    except ValueError as e:
      mismatch = e.msg
    check mismatch == "checkpoint mismatch: version"
    rewriteCheckpointI32(path, pathConv, 12, 99'i32)
    mismatch = ""
    try:
      loadCheckpoint(mb, pathConv)
    except ValueError as e:
      mismatch = e.msg
    check mismatch == "checkpoint mismatch: mass convention"
    # mismatch: same file against a different-parameter sampler
    var mc = newMdHmc(@[0.0, 0.6], @[4, 2, 2], 0.6, 770001'u64)
    expect ValueError:
      loadCheckpoint(mc, path)
    var md = newMdHmc(mass2, @[4, 2, 2], 0.7, 770001'u64)
    expect ValueError:
      loadCheckpoint(md, path)
    echo "  corrupted, legacy-version, convention, and parameter mismatches raise"
    removeFile(path)
    removeFile(pathBad)
    removeFile(pathV1)
    removeFile(pathConv)

# --- 9. end-to-end -----------------------------------------------------------------

suite "end-to-end run (L=1, nt=6, Nf=2)":

  test "acceptance > 0.5 and <exp(-dH)> ~ 1 over 20 trajectories":
    var m = newMdHmc(mass2, @[8, 2, 2], 0.6, 880001'u64)
    const
      warm = 3
      ntraj = 20
    var
      sE = 0.0
      sE2 = 0.0
      nAcc = 0
      nMeas = 0
      tsum = 0.0
      tmax = 0.0
    for k in 1..ntraj:
      m.forceAccept = k <= warm
      let t0 = epochTime()
      m.update
      let dt = epochTime() - t0
      tsum += dt
      tmax = max(tmax, dt)
      if k > warm:
        inc nMeas
        if m.accepted: inc nAcc
        sE += m.expmDeltaH
        sE2 += m.expmDeltaH*m.expmDeltaH
      if k mod 5 == 0:
        let w = windowCheck(m)
        check w.inside
        echo &"  traj {k}: dH {m.deltaH: .4e}  acc {m.accepted}" &
             &"  window [{w.lo:.4f}, {w.hi:.4f}] inside {w.inside}"
    let
      acc = float(nAcc)/float(nMeas)
      eMean = sE/float(nMeas)
      eErr = sqrt((sE2/float(nMeas) - eMean*eMean)/float(nMeas - 1))
    echo &"  acceptance = {acc:.3f} ({nAcc}/{nMeas})"
    echo &"  <exp(-dH)> = {eMean:.6f} +- {eErr:.6f}  (pull {abs(eMean-1.0)/eErr:.2f})"
    echo &"  seconds/trajectory: mean {tsum/float(ntraj):.3f}  max {tmax:.3f}"
    check acc > 0.5
    check abs(eMean - 1.0) < 5.0*eErr
    check m.pf.actOp.stats.ok and m.pf.frcOp.stats.ok

# --- allocation regression ----------------------------------------------------------

suite "allocation regression":

  test "live memory unchanged across GC_fullCollect over 16 force+action rounds":
    var m = newMdHmc(mass2, @[4, 2, 2], 0.6, 990001'u64)
    m.forceAccept = true
    m.update
    m.forceAccept = false
    var f = newGauge(lat)
    GC_fullCollect()
    let
      occ0 = getOccupiedMem()
      raw0 = getRawMemAllocated()
    for k in 0..<16:
      pfForce(m.pf, f, m.u, 1)
      pfForce(m.pf, f, m.u, 2)
      discard projectKernel(lat, f)
      discard pfAction(m.pf, m.u)
    GC_fullCollect()
    let
      occ1 = getOccupiedMem()
      raw1 = getRawMemAllocated()
    echo &"  occupied {occ0} -> {occ1}  raw {raw0} -> {raw1}"
    check occ1 == occ0
    check raw1 == raw0
