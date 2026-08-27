## HMC trajectory driver for dynamical overlap + non-compact gauge (WP-H).
##
## Normative reference: doc/02-formulation.md sections 4.2 and 5,
## doc/04-interfaces.md section 13.  Metropolis via QEX hmc/metropolis.nim
## (`RadialHmc` subclasses MetropolisRoot; hooks start/getH/generate/globalRand/
## accept/reject below); MD via the mdevolve nimble package, one Omelyan 2MN
## per level sharing the position update:
##   level 0            gauge force (innermost, most steps)
##   level 1..nframe    Hasenbusch frames heaviest -> lightest (outermost)
##
## H = |p|^2/2 + S_g(theta) + S_pf(theta; phi), one real momentum per link
## (spatial and temporal), S_pf and the heatbath with actOp (order 31), MD
## forces with frcOp (order 11).
##
## Gauge zero modes (WP-G machinery; ker M = gauge orbit + the uniform temporal
## Polyakov mode, dim n_V L_t): `projectKernel` is applied to
##   (a) the committed field -- at construction/load and again at every commit,
##       which keeps the field entering each trajectory transverse while `reject`
##       stays a bitwise restore of the trajectory's start field;
##   (b) the refreshed momentum;
##   (c) every MD force at every level.
## All three are required: at fixed phi the extended action has a longitudinal
## force even though the integrated determinant is gauge invariant.
##
## Randomness is trajectory addressed: every draw comes from a Threefry stream
## keyed by splitmix64-mixing (baseSeed, trajectory number, purpose, copy,
## frame).  No generator state is serialized; a restart depends only on the
## committed trajectory counter.

import std/[math, streams]
import ../core/lattice
import ../core/spinor
import ../ops/gaugeact
import pseudofermion
import hmc/metropolis
import mdevolve

export gaugeact, pseudofermion, metropolis

# --- trajectory-addressed randomness -----------------------------------------

const
  rkMomentum* = 1            ## purpose keys for `keyedRng`
  rkAccept* = 2
  rkPseudo* = 3

func splitmix64(x: uint64): uint64 =
  var z = x + 0x9E3779B97F4A7C15'u64
  z = (z xor (z shr 30)) * 0xBF58476D1CE4E5B9'u64
  z = (z xor (z shr 27)) * 0x94D049BB133111EB'u64
  z xor (z shr 31)

func mixKey*(seed: uint64, traj, purpose: int, copy = 0, frame = 0): uint64 =
  ## Splitmix64 chain over (seed, traj, purpose, copy, frame).
  result = splitmix64(seed)
  result = splitmix64(result xor uint64(traj))
  result = splitmix64(result xor uint64(purpose))
  result = splitmix64(result xor uint64(copy))
  result = splitmix64(result xor uint64(frame))

proc keyedRng*(seed: uint64, traj, purpose: int, copy = 0, frame = 0): Threefry4x64 =
  result.seedIndep(mixKey(seed, traj, purpose, copy, frame),
                   mixKey(seed, traj, purpose + 0x100, copy, frame))

proc gaussian*(x: var Gauge, r: var Threefry4x64) =
  ## One unit normal per link, so kinetic = |p|^2/2 gives <|p|^2> = dof.
  for i in 0..<x.s.len: x.s[i] = r.gaussian
  for i in 0..<x.t.len: x.t[i] = r.gaussian

# --- the sampler --------------------------------------------------------------

type
  RadialHmc* = ref object of MetropolisRoot
    l*: Lat
    bt*: Beta
    pf*: Pf                  ## nil = pure gauge
    u*, p*: Gauge            ## committed field and momentum
    uOld: Gauge              ## saved by `start`; `reject` restores it bitwise
    tau*: float
    seed*: uint64
    traj*: int               ## trajectory counter; incremented by `start`
    fcount*: seq[int]        ## force evaluations per level since last clear
    projR2*: float           ## kernel-projection tolerance (relative)
    forceAccept*, forceReject*: bool
    levels: seq[Integrator]  ## one 2MN per level, [0] = gauge
    evo: ParIntegrator       ## nil when there is a single level
    f: Gauge                 ## force scratch

proc applyForce(m: RadialHmc, lv: int, t: float) =
  ## p -= t * P f_lv, P = projectKernel: rule (c) of the module header.
  if lv == 0: gaugeForce(m.l, m.f, m.u, m.bt)
  else: pfForce(m.pf, m.f, m.u, lv)
  discard projectKernel(m.l, m.f, m.projR2)
  axpy(m.p, -t, m.f)
  inc m.fcount[lv]

proc newRadialHmc*(l: Lat, bt: Beta, pf: Pf, tau: float, steps: seq[int],
                   seed: uint64, projR2 = 1e-24): RadialHmc =
  ## `steps` is per level: [gauge, heaviest frame, ..., lightest frame], so
  ## steps.len == 1 + (pf.isNil ? 0 : pf.nframe).  The committed field starts
  ## at zero; assign to `m.u` (then `projectKernel`) or `loadCheckpoint` after.
  doAssert steps.len == 1 + (if pf.isNil: 0 else: pf.nframe)
  for s in steps: doAssert s >= 1
  var m = RadialHmc(l: l, bt: bt, pf: pf, tau: tau, seed: seed, projR2: projR2)
  init(m)
  m.u = newGauge(l)
  m.p = newGauge(l)
  m.uOld = newGauge(l)
  m.f = newGauge(l)
  m.fcount = newSeq[int](steps.len)
  let mm = m
  proc mdt(t: float) = axpy(mm.u, t, mm.p)
  proc mdv(ts: openarray[float]) =
    for lv in 0..<ts.len:
      if ts[lv] != 0.0: applyForce(mm, lv, ts[lv])
  let (vv, tt) = newIntegratorPair(mdv, mdt)
  for lv in 0..<steps.len:
    m.levels.add mkOmelyan2MN(steps = steps[lv], V = vv[lv], T = tt)
  if steps.len > 1: m.evo = newParallelEvolution(m.levels)
  m

proc setSteps*(m: RadialHmc, steps: openArray[int]) =
  ## Per-level step counts, same order as the constructor.
  doAssert steps.len == m.levels.len
  for lv in 0..<steps.len: m.levels[lv].steps = steps[lv]

proc clearForceCounts*(m: RadialHmc) =
  for lv in 0..<m.fcount.len: m.fcount[lv] = 0

proc mdEvolve*(m: RadialHmc) =
  ## One MD trajectory of length tau on the current (u, p, phi).
  if m.evo.isNil:
    m.levels[0].evolve m.tau
    m.levels[0].finish
  else:
    m.evo.evolve m.tau
    m.evo.finish

proc hmcH*(m: RadialHmc): float =
  ## H = |p|^2/2 + S_g + S_pf (actOp).
  result = 0.5*norm2(m.p) + gaugeAction(m.l, m.u, m.bt)
  if not m.pf.isNil: result += pfAction(m.pf, m.u)

proc refreshMomentum*(m: RadialHmc, traj: int) =
  var r = keyedRng(m.seed, traj, rkMomentum)
  gaussian(m.p, r)
  discard projectKernel(m.l, m.p, m.projR2)

proc refreshPseudo*(m: RadialHmc, traj: int) =
  for c in 0..<m.pf.ncopy:
    for i in 0..<m.pf.nframe:
      var r = keyedRng(m.seed, traj, rkPseudo, c, i)
      refreshFrame(m.pf, m.u, c, i, r)

# --- metropolis hooks (required: start getH generate globalRand accept reject)

proc start*(m: var RadialHmc) =
  inc m.traj
  m.uOld := m.u
  refreshMomentum(m, m.traj)
  if not m.pf.isNil: refreshPseudo(m, m.traj)

proc getH*(m: RadialHmc): float = hmcH(m)

proc generate*(m: var RadialHmc) = mdEvolve(m)

proc globalRand*(m: RadialHmc): float =
  if m.forceAccept: return 0.0     # below any pAccept; no stream is consumed
  if m.forceReject: return 2.0     # above any pAccept
  var r = keyedRng(m.seed, m.traj, rkAccept)
  r.uniform

proc accept*(m: var RadialHmc) =
  discard projectKernel(m.l, m.u, m.projR2)   # rule (a): commit transverse

proc reject*(m: var RadialHmc) =
  m.u := m.uOld                               # bitwise restore

# --- diagnostics ---------------------------------------------------------------

type RevInfo* = object
  ## Reversibility drifts: forward tau, p -> -p, forward tau, p -> -p.
  dh*: float                 ## H(round trip) - H(start)
  duRms*, duMax*: float      ## per-link |u - u0|
  dpRms*, dpMax*: float      ## per-link |p - p0|
  divP*: float               ## |d^dag p|^2 / |p|^2 of the round-trip momentum
  flatP*: float              ## |mean temporal p| * sqrt(nlink) / |p|

proc transversality*(l: Lat, p: Gauge): tuple[divP, flatP: float] =
  var d: seq[float]
  divergence(l, d, p)
  var d2 = 0.0
  for x in d: d2 += x*x
  var s = 0.0
  for x in p.t: s += x
  s /= float(p.t.len)
  let p2 = norm2(p)
  (d2/p2, abs(s)*sqrt(float(nlink(l))/p2))

proc reversibilityCheck*(m: RadialHmc): RevInfo =
  ## Round-trip integration on the current (u, p, phi); restores u and p.
  var u0 = newGauge(m.l)
  var p0 = newGauge(m.l)
  u0 := m.u
  p0 := m.p
  let h0 = hmcH(m)
  mdEvolve(m)
  scale(m.p, -1.0)
  mdEvolve(m)
  scale(m.p, -1.0)
  result.dh = hmcH(m) - h0
  let t = transversality(m.l, m.p)
  result.divP = t.divP
  result.flatP = t.flatP
  var s = 0.0
  for i in 0..<m.u.s.len:
    let d = m.u.s[i] - u0.s[i]
    s += d*d
    result.duMax = max(result.duMax, abs(d))
  for i in 0..<m.u.t.len:
    let d = m.u.t[i] - u0.t[i]
    s += d*d
    result.duMax = max(result.duMax, abs(d))
  result.duRms = sqrt(s/float(nlink(m.l)))
  s = 0.0
  for i in 0..<m.p.s.len:
    let d = m.p.s[i] - p0.s[i]
    s += d*d
    result.dpMax = max(result.dpMax, abs(d))
  for i in 0..<m.p.t.len:
    let d = m.p.t[i] - p0.t[i]
    s += d*d
    result.dpMax = max(result.dpMax, abs(d))
  result.dpRms = sqrt(s/float(nlink(m.l)))
  m.u := u0
  m.p := p0

proc windowCheck*(m: RadialHmc): tuple[smin, smax, lo, hi: float, inside: bool] =
  ## kernelWindow on the committed field, against BOTH frozen rationals.
  ## `inside == false` is the caller's hard stop (never rebuild mid-ensemble).
  result = kernelWindow(m.pf.actOp, m.u)
  if not (result.lo >= m.pf.frcOp.rat.smin and result.hi <= m.pf.frcOp.rat.smax):
    result.inside = false

# --- checkpoint ----------------------------------------------------------------

const
  ckptMagic = 0x51455852484D4331'u64   ## "QEXRHMC1"
  ckptVersion = 2'i32                  ## v2 adds the mass-convention identifier

func fnv(s: string): uint64 =
  result = 0xcbf29ce484222325'u64
  for ch in s:
    result = (result xor uint64(ord(ch))) * 0x100000001b3'u64

proc saveCheckpoint*(m: RadialHmc, path: string) =
  ## Versioned binary manifest + gauge field + trailing FNV-1a of everything.
  ## Version 1 checkpoints used the legacy additive mass and are intentionally
  ## incompatible.  Version 2 records the standard-overlap convention explicitly.
  var st = newStringStream()
  st.write ckptMagic
  st.write ckptVersion
  st.write ovMassConventionId
  st.write int32(m.l.sph.lev)
  st.write int32(m.l.nt)
  st.write m.l.at
  st.write m.bt.g2
  st.write int32(ord(m.bt.conv))
  if m.pf.isNil:
    st.write int32(0)
  else:
    st.write int32(m.pf.nf)
    st.write m.pf.actOp.m
    st.write m.pf.actOp.rat.hash
    st.write m.pf.frcOp.rat.hash
    st.write int32(m.pf.masses.len)
    for x in m.pf.masses: st.write x
  st.write m.tau
  st.write int32(m.levels.len)
  for lv in m.levels: st.write int32(lv.steps)
  st.write m.seed
  st.write int64(m.traj)
  st.write int32(m.u.s.len)
  st.write int32(m.u.t.len)
  for x in m.u.s: st.write x
  for x in m.u.t: st.write x
  let body = st.data
  st.write fnv(body)
  writeFile(path, st.data)

template ckWant(cond: bool, what: string) =
  if not cond:
    raise newException(ValueError, "checkpoint mismatch: " & what)

proc loadCheckpoint*(m: RadialHmc, path: string) =
  ## Validates every manifest field against the live configuration and refuses
  ## mismatches; a payload-hash failure means corruption.  On success sets the
  ## gauge field and the trajectory counter, bitwise -- the stored field was
  ## projected when it was committed, so rule (a) holds without re-projecting
  ## (which would break the bitwise restart guarantee).
  let data = readFile(path)
  ckWant data.len > 8, "file too short"
  var st = newStringStream(data)
  var body = data
  body.setLen(data.len - 8)
  var h: uint64
  st.setPosition(data.len - 8)
  st.read h
  ckWant h == fnv(body), "payload hash (corrupted file)"
  st.setPosition(0)
  var u64: uint64
  var i32: int32
  var i64: int64
  var f64: float
  st.read u64; ckWant u64 == ckptMagic, "magic"
  st.read i32; ckWant i32 == ckptVersion, "version"
  st.read i32; ckWant i32 == ovMassConventionId, "mass convention"
  st.read i32; ckWant int(i32) == m.l.sph.lev, "lev"
  st.read i32; ckWant int(i32) == m.l.nt, "nt"
  st.read f64; ckWant f64 == m.l.at, "at"
  st.read f64; ckWant f64 == m.bt.g2, "g2"
  st.read i32; ckWant int(i32) == ord(m.bt.conv), "geometry convention"
  st.read i32
  if m.pf.isNil:
    ckWant i32 == 0, "nf (checkpoint has fermions)"
  else:
    ckWant int(i32) == m.pf.nf, "nf"
    st.read f64; ckWant f64 == m.pf.actOp.m, "M"
    st.read u64; ckWant u64 == m.pf.actOp.rat.hash, "action rational hash"
    st.read u64; ckWant u64 == m.pf.frcOp.rat.hash, "force rational hash"
    st.read i32; ckWant int(i32) == m.pf.masses.len, "mass count"
    for x in m.pf.masses:
      st.read f64; ckWant f64 == x, "mass value"
  st.read f64; ckWant f64 == m.tau, "tau"
  st.read i32; ckWant int(i32) == m.levels.len, "level count"
  for lv in m.levels:
    st.read i32; ckWant int(i32) == lv.steps, "step count"
  st.read u64; ckWant u64 == m.seed, "seed"
  st.read i64
  st.read i32; ckWant int(i32) == m.u.s.len, "spatial link count"
  st.read i32; ckWant int(i32) == m.u.t.len, "temporal link count"
  for i in 0..<m.u.s.len: st.read m.u.s[i]
  for i in 0..<m.u.t.len: st.read m.u.t[i]
  m.traj = int(i64)
