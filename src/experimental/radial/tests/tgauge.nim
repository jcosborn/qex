#RUNCMD $RUN1
## WP-G: non-compact U(1) gauge action, zero-mode projection, exact heatbath and the
## pseudo-inverse propagator.  doc/03-targets.md T1.5g, doc/02-formulation.md section 5.
##
## The oracle everywhere is the dense incidence matrix: S = theta^T C^T W C theta / 2
## assembled column by column, and its eigendecomposition for the pseudo-inverse and
## for the kernel count.  Nothing in the oracle shares code with `ops/gaugeact.nim`
## beyond the geometry tables.

import std/[math, strformat, unittest]
import eigens/linalgFuncs
import ../ops/gaugeact

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))

# --- dense incidence oracle -------------------------------------------------

type Dense = object
  n: int                     ## number of links
  m: seq[float]              ## column-major M = C^T W C
  nplaq: int

proc plaqRow(l: Lat, kind, idx, t: int, row: var seq[float]) =
  ## Dense incidence row of a plaquette.  kind 0 = spatial triangle `idx`,
  ## kind 1 = temporal plaquette of spatial edge `idx`.  Written literally from
  ## doc/02 section 5, independently of gaugeact's loops.
  for i in 0..<row.len: row[i] = 0.0
  let ns = l.sph.ne*l.nt
  if kind == 0:
    let fc = l.sph.faces[idx]
    for i in 0..2:
      row[eIdx(l, fc.e[i], t)] += float(fc.s[i])
  else:
    let ed = l.sph.edges[idx]
    row[eIdx(l, idx, t)] += 1.0
    row[ns + tIdx(l, ed.b, t)] += 1.0
    row[eIdx(l, idx, t+1)] -= 1.0
    row[ns + tIdx(l, ed.a, t)] -= 1.0

proc newDense(l: Lat, b: Beta): Dense =
  result.n = nlink(l)
  result.m = newSeq[float](result.n*result.n)
  var row = newSeq[float](result.n)
  var nz: seq[int]
  for kind in 0..1:
    let nk = if kind == 0: l.sph.nf else: l.sph.ne
    for idx in 0..<nk:
      let w = if kind == 0: b.face[idx] else: b.edge[idx]
      for t in 0..<l.nt:
        plaqRow(l, kind, idx, t, row)
        nz.setLen 0
        for i in 0..<result.n:
          if row[i] != 0.0: nz.add i
        for i in nz:
          for j in nz:
            result.m[i + result.n*j] += w*row[i]*row[j]
        inc result.nplaq

proc act(d: Dense, x: openArray[float]): float =
  for j in 0..<d.n:
    var s = 0.0
    for i in 0..<d.n: s += d.m[i + d.n*j]*x[i]
    result += 0.5*s*x[j]

proc mul(d: Dense, x: openArray[float]): seq[float] =
  result = newSeq[float](d.n)
  for j in 0..<d.n:
    let xj = x[j]
    if xj == 0.0: continue
    for i in 0..<d.n: result[i] += d.m[i + d.n*j]*xj

proc eigen(d: Dense): tuple[v: seq[Complex64], w: seq[float]] =
  ## Real symmetric eigenproblem through the Hermitian LAPACK path (zheev).
  var h = newSeq[Complex64](d.n*d.n)
  for i in 0..<d.n*d.n: h[i] = complex64(d.m[i], 0.0)
  var ew = newSeq[float](d.n)
  zeigs(cast[ptr float64](addr h[0]), addr ew[0], d.n)
  (h, ew)

proc pinv(d: Dense, ev: tuple[v: seq[Complex64], w: seq[float]], tol: float): seq[float] =
  ## Moore-Penrose pseudo-inverse, kernel dropped.
  result = newSeq[float](d.n*d.n)
  for k in 0..<d.n:
    if abs(ev.w[k]) <= tol: continue
    let iw = 1.0/ev.w[k]
    for j in 0..<d.n:
      let c = iw*ev.v[j + d.n*k].re
      if c == 0.0: continue
      for i in 0..<d.n: result[i + d.n*j] += c*ev.v[i + d.n*k].re

proc expm(d: Dense, ev: tuple[v: seq[Complex64], w: seq[float]], s: float): seq[float] =
  ## exp(-M s), column major.
  result = newSeq[float](d.n*d.n)
  for k in 0..<d.n:
    let e = exp(-s*ev.w[k])
    for j in 0..<d.n:
      let c = e*ev.v[j + d.n*k].re
      for i in 0..<d.n: result[i + d.n*j] += c*ev.v[i + d.n*k].re

# --- helpers ----------------------------------------------------------------

proc randGauge(l: Lat, r: var Threefry4x64, sc = 1.0): Gauge =
  result = newGauge(l)
  for i in 0..<result.s.len: result.s[i] = sc*r.gaussian
  for i in 0..<result.t.len: result.t[i] = sc*r.gaussian

proc randSeqF(n: int, r: var Threefry4x64, sc = 1.0): seq[float] =
  result = newSeq[float](n)
  for i in 0..<n: result[i] = sc*r.gaussian

proc maxAbs(x: openArray[float]): float =
  for v in x: result = max(result, abs v)

proc rng(sed: int): Threefry4x64 =
  result.seedIndep(sed, 0)

const
  lev = 1
  nt = 4
  at = 0.35
  g2 = 1.7

let
  sph = newSphere(lev)
  lat = newLat(sph, nt, at)
  bet = newBeta(lat, g2)
  den = newDense(lat, bet)
  ev = den.eigen
  evtol = 1e-9*maxAbs(ev.w)

suite "gauge action: structure":

  test "counts and the dense oracle agree on dimensions":
    check den.n == (sph.ne + sph.nv)*nt
    check den.nplaq == (sph.nf + sph.ne)*nt
    echo &"  L={lev} nt={nt} nlink={den.n} nplaq={den.nplaq} abar/at={lat.asOverAt:.4f}"

  test "sum_faces Theta_tri = 0 on every time slice":
    var r = rng(11)
    let u = randGauge(lat, r)
    var sc = 0.0
    for f in 0..<sph.nf:
      for t in 0..<nt: sc = max(sc, abs(plaqSpatial(lat, u, f, t)))
    for t in 0..<nt:
      var s = 0.0
      for f in 0..<sph.nf: s += plaqSpatial(lat, u, f, t)
      check abs(s) < 1e-12*sc

  test "S = 0 on a pure-gauge configuration theta = d alpha":
    var r = rng(12)
    let al = randSeqF(lat.nsite, r)
    var p = newGauge(lat)
    gradient(lat, p, al)
    let s = gaugeAction(lat, p, bet)
    echo &"  S(d alpha) = {s:.6e}"
    check abs(s) < 1e-25

  test "gauge invariance of the action":
    var r = rng(13)
    let
      u = randGauge(lat, r)
      al = randSeqF(lat.nsite, r, 3.0)
    var p = newGauge(lat)
    gradient(lat, p, al)
    var v = newGauge(lat)
    v := u
    axpy(v, 1.0, p)
    let
      s0 = gaugeAction(lat, u, bet)
      s1 = gaugeAction(lat, v, bet)
    echo &"  S = {s0:.15e}  dS/S = {abs(s1-s0)/abs(s0):.3e}"
    check abs(s1 - s0) < 1e-12*abs(s0)

  test "d^dagger is the adjoint of d":
    var r = rng(14)
    let
      p = randGauge(lat, r)
      al = randSeqF(lat.nsite, r)
    var dal = newGauge(lat)
    gradient(lat, dal, al)
    var dp = newSeq[float](lat.nsite)
    divergence(lat, dp, p)
    var a = 0.0
    for i in 0..<al.len: a += al[i]*dp[i]
    let b = dot(p, dal)
    check abs(a - b) <= 1e-13*(abs(a) + abs(b))

  test "laplace == divergence . gradient":
    var r = rng(15)
    let al = randSeqF(lat.nsite, r)
    var p = newGauge(lat)
    gradient(lat, p, al)
    var d0 = newSeq[float](lat.nsite)
    divergence(lat, d0, p)
    var d1 = newSeq[float](lat.nsite)
    laplace(lat, d1, al)
    var e = 0.0
    for i in 0..<d0.len: e = max(e, abs(d0[i] - d1[i]))
    check e < 1e-13*maxAbs(d0)

suite "gauge action: action and force":

  test "action matches the dense incidence oracle":
    var r = rng(21)
    let u = randGauge(lat, r)
    let
      s0 = gaugeAction(lat, u, bet)
      s1 = den.act(u.toSeq)
    echo &"  S = {s0:.15e}   dense = {s1:.15e}   rel = {abs(s0-s1)/abs(s0):.3e}"
    check abs(s0 - s1) < 1e-12*abs(s0)

  test "action parts sum to the total and are separately positive":
    var r = rng(22)
    let u = randGauge(lat, r)
    let p = gaugeActionParts(lat, u, bet)
    check p.sp > 0.0
    check p.tp > 0.0
    check abs(p.sp + p.tp - gaugeAction(lat, u, bet)) < 1e-13*(p.sp + p.tp)

  test "force matches the dense incidence oracle":
    var r = rng(23)
    let u = randGauge(lat, r)
    var f = newGauge(lat)
    gaugeForce(lat, f, u, bet)
    let
      fa = f.toSeq
      fd = den.mul(u.toSeq)
    var e = 0.0
    for i in 0..<fa.len: e = max(e, abs(fa[i] - fd[i]))
    echo &"  max |f - M theta| = {e:.3e}  (|f|_max = {maxAbs(fa):.4f})"
    check e < 1e-12*maxAbs(fd)

  test "force vs centered finite differences, four step sizes":
    var r = rng(24)
    var u = randGauge(lat, r)
    var f = newGauge(lat)
    gaugeForce(lat, f, u, bet)
    let fa = f.toSeq
    var best = 1e300
    for h in [1e-1, 1e-2, 1e-3, 1e-4]:
      var e = 0.0
      for k in countup(0, den.n-1, 7):     # every 7th link: 24 of them, all types
        var v = newGauge(lat)
        v := u
        if k < v.s.len: v.s[k] += h else: v.t[k - v.s.len] += h
        let sp = gaugeAction(lat, v, bet)
        v := u
        if k < v.s.len: v.s[k] -= h else: v.t[k - v.s.len] -= h
        let sm = gaugeAction(lat, v, bet)
        e = max(e, abs((sp - sm)/(2.0*h) - fa[k]))
      echo &"  h = {h:.0e}   max |fd - analytic| = {e:.3e}"
      best = min(best, e)
    check best < 1e-8

  test "force is orthogonal to every pure-gauge direction":
    var r = rng(25)
    let u = randGauge(lat, r)
    var f = newGauge(lat)
    gaugeForce(lat, f, u, bet)
    var e = 0.0
    for k in 0..3:
      let al = randSeqF(lat.nsite, r)
      var p = newGauge(lat)
      gradient(lat, p, al)
      e = max(e, abs(dot(f, p))/sqrt(norm2(f)*norm2(p)))
    echo &"  max |<f, d alpha>|/(|f||d alpha|) = {e:.3e}"
    check e < 1e-14

  test "force is also orthogonal to the flat Polyakov direction":
    var r = rng(26)
    let u = randGauge(lat, r)
    var f = newGauge(lat)
    gaugeForce(lat, f, u, bet)
    var s = 0.0
    for x in f.t: s += x
    check abs(s) < 1e-12*sqrt(norm2(f)*float(f.t.len))

  test "gaugeForce allocates nothing":
    var r = rng(27)
    let u = randGauge(lat, r)
    var f = newGauge(lat)
    gaugeForce(lat, f, u, bet)
    let m0 = getOccupiedMem()
    for i in 0..<64: gaugeForce(lat, f, u, bet)
    check getOccupiedMem() == m0

suite "gauge action: zero modes":

  test "kernel of M is the gauge orbit plus one flat direction":
    var nz = 0
    for w in ev.w:
      if abs(w) <= evtol: inc nz
    echo &"  dim ker M = {nz}  (nv*nt = {sph.nv*nt}), rank = {den.n - nz} " &
         &"(ne*nt = {sph.ne*nt}), |w| in [{maxAbs(ev.w):.3e} max]"
    check nz == sph.nv*nt
    check den.n - nz == sph.ne*nt

  test "every kernel vector is killed by projectKernel":
    var e = 0.0
    for k in 0..<den.n:
      if abs(ev.w[k]) > evtol: continue
      var p = newGauge(lat)
      for i in 0..<den.n:
        if i < p.s.len: p.s[i] = ev.v[i + den.n*k].re
        else: p.t[i - p.s.len] = ev.v[i + den.n*k].re
      discard projectKernel(lat, p)
      e = max(e, sqrt(norm2 p))
    echo &"  max |P_perp v_kernel| = {e:.3e}"
    check e < 1e-10

  test "projectGauge kills a pure-gauge field":
    var r = rng(31)
    let al = randSeqF(lat.nsite, r)
    var p = newGauge(lat)
    gradient(lat, p, al)
    let n0 = sqrt(norm2 p)
    let info = projectGauge(lat, p)
    echo &"  |d alpha| {n0:.4f} -> {sqrt(norm2 p):.3e}  ({info.iters} its)"
    check sqrt(norm2 p) < 1e-10*n0

  test "projectGauge drives the divergence to zero and is idempotent":
    var r = rng(32)
    var p = randGauge(lat, r)
    var d0 = newSeq[float](lat.nsite)
    divergence(lat, d0, p)
    var n0 = 0.0
    for v in d0: n0 += v*v
    let info = projectGauge(lat, p)
    var d1 = newSeq[float](lat.nsite)
    divergence(lat, d1, p)
    var n1 = 0.0
    for v in d1: n1 += v*v
    echo &"  |div p|^2: {n0:.6e} -> {n1:.3e}  ratio {n1/n0:.3e}  ({info.iters} its, r2 {info.r2:.2e})"
    check info.converged
    check n1 < 1e-20*n0
    var q = newGauge(lat)
    q := p
    discard projectGauge(lat, q)
    var e = 0.0
    for i in 0..<q.s.len: e = max(e, abs(q.s[i] - p.s[i]))
    for i in 0..<q.t.len: e = max(e, abs(q.t[i] - p.t[i]))
    echo &"  idempotency: max change on the second pass = {e:.3e}"
    check e < 1e-12*sqrt(norm2 p)

  test "projectGauge leaves the action invariant":
    var r = rng(33)
    var p = randGauge(lat, r)
    let s0 = gaugeAction(lat, p, bet)
    discard projectGauge(lat, p)
    let s1 = gaugeAction(lat, p, bet)
    check abs(s1 - s0) < 1e-12*abs(s0)

suite "gauge action: pseudo-inverse (T1.5g)":

  let mp = den.pinv(ev, evtol)

  test "double-CG pseudo-inverse matches the dense one":
    var e = 0.0
    var sc = 0.0
    for k in countup(0, den.n-1, 13):
      let col = gaugePropagator(lat, g2, k, 1e-24, 40000)
      for i in 0..<den.n:
        e = max(e, abs(col[i] - mp[i + den.n*k]))
        sc = max(sc, abs(mp[i + den.n*k]))
    echo &"  max |Mtilde^-1 - dense pinv| = {e:.3e}  (scale {sc:.4f})"
    check e < 1e-10

  test "pseudo-inverse is symmetric and annihilates the kernel":
    let
      c0 = gaugePropagator(lat, g2, 3, 1e-24, 40000)
      c1 = gaugePropagator(lat, g2, den.n - 5, 1e-24, 40000)
    check abs(c0[den.n-5] - c1[3]) < 1e-12
    var p = fromSeq(lat, c0)
    let n0 = sqrt(norm2 p)
    discard projectKernel(lat, p)
    check abs(sqrt(norm2 p) - n0) < 1e-10*n0     # already transverse

  test "the solve is invariant under adding kernel vectors to the source":
    var r = rng(41)
    let bt = newBeta(lat, g2)
    var b = newGauge(lat)
    triSource(lat, b, 0, 0)
    var x0 = newGauge(lat)
    discard pseudoSolve(lat, x0, b, bt, 1e-24, 40000)
    let al = randSeqF(lat.nsite, r, 2.0)
    var k = newGauge(lat)
    gradient(lat, k, al)
    for i in 0..<k.t.len: k.t[i] += 0.75      # plus the flat Polyakov direction
    var b2 = newGauge(lat)
    b2 := b
    axpy(b2, 1.0, k)
    var x1 = newGauge(lat)
    discard pseudoSolve(lat, x1, b2, bt, 1e-24, 40000)
    var e = 0.0
    for i in 0..<x0.s.len: e = max(e, abs(x0.s[i] - x1.s[i]))
    for i in 0..<x0.t.len: e = max(e, abs(x0.t[i] - x1.t[i]))
    echo &"  max |x(b) - x(b + ker)| = {e:.3e}  (|x| = {sqrt(norm2 x0):.4f})"
    check e < 1e-10

  test "kernel-regularized solve agrees with the double CG, and goes further":
    ## regSolve is what the correlator app uses: the double CG has no control over
    ## the kernel component roundoff injects into the Krylov space, and past
    ## r2 ~ 1e-26 that component diverges (measured on L=1, nt=120).  A is
    ## positive definite, so regSolve just converges.
    let bt = newBeta(lat, g2)
    var reg = newRegOp(lat, bt)
    var b = newGauge(lat)
    triSource(lat, b, 5, 2)
    var x0 = newGauge(lat)
    discard pseudoSolve(lat, x0, b, bt, 1e-24, 40000)
    var x1 = newGauge(lat)
    let i1 = regSolve(lat, x1, b, reg, 1e-24, 40000)
    var e = 0.0
    for i in 0..<x0.s.len: e = max(e, abs(x0.s[i] - x1.s[i]))
    for i in 0..<x0.t.len: e = max(e, abs(x0.t[i] - x1.t[i]))
    echo &"  max |regSolve - pseudoSolve| = {e:.3e}  ({i1.iters} its, r2 {i1.r2:.2e})"
    check i1.converged
    check e < 1e-12
    # the regularized operator is symmetric positive definite on the whole space
    var r = rng(42)
    let
      p = randGauge(lat, r)
      q = randGauge(lat, r)
    var ap = newGauge(lat)
    var aq = newGauge(lat)
    applyReg(lat, reg, ap, p)
    applyReg(lat, reg, aq, q)
    check abs(dot(q, ap) - dot(p, aq)) < 1e-12*abs(dot(q, ap))
    check dot(p, ap) > 0.0
    var kv = newGauge(lat)
    gradient(lat, kv, randSeqF(lat.nsite, r))
    applyReg(lat, reg, ap, kv)
    check dot(kv, ap) > 0.0                # positive definite on the kernel too
    # ... and it agrees with M on the transverse subspace
    var pt = newGauge(lat)
    pt := p
    discard projectKernel(lat, pt)
    applyReg(lat, reg, ap, pt)
    var mp = newGauge(lat)
    gaugeForce(lat, mp, pt, bt)
    var e2 = 0.0
    for i in 0..<mp.s.len: e2 = max(e2, abs(mp.s[i] - ap.s[i]))
    for i in 0..<mp.t.len: e2 = max(e2, abs(mp.t[i] - ap.t[i]))
    check e2 < 1e-10*sqrt(norm2 mp)

  test "<Theta_tri Theta_tri> is gauge invariant and matches the oracle":
    let bt = newBeta(lat, g2)
    var b = newGauge(lat)
    triSource(lat, b, 0, 0)
    var x = newGauge(lat)
    discard pseudoSolve(lat, x, b, bt, 1e-24, 40000)
    let bs = b.toSeq
    var e = 0.0
    for f in 0..<sph.nf:
      for t in 0..<nt:
        var c = newGauge(lat)
        triSource(lat, c, f, t)
        let cs = c.toSeq
        var oracle = 0.0
        for i in 0..<den.n:
          for j in 0..<den.n: oracle += cs[i]*mp[i + den.n*j]*bs[j]
        e = max(e, abs(dot(c, x) - oracle))
    echo &"  max |<Theta Theta>_cg - dense| over all (f,t) = {e:.3e}"
    check e < 1e-10

suite "gauge action: exact heatbath":

  test "<S> = rank(M)/2 and <Theta^2> matches the propagator":
    let
      bt = newBeta(lat, g2)
      nsamp = 600
    var rank = 0
    for w in ev.w:
      if abs(w) > evtol: inc rank
    let mp = den.pinv(ev, evtol)
    # exact <Theta_f(t)^2> for face 0 at t = 0
    var c = newGauge(lat)
    triSource(lat, c, 0, 0)
    let cs = c.toSeq
    var exact = 0.0
    for i in 0..<den.n:
      for j in 0..<den.n: exact += cs[i]*mp[i + den.n*j]*cs[j]
    var
      r = rng(51)
      u = newGauge(lat)
      sm, s2 = 0.0
      tm, t2 = 0.0
      its = 0
    for n in 0..<nsamp:
      let info = heatbath(lat, u, bt, r, 1e-22, 20000)
      check info.converged
      its += info.iters
      let s = gaugeAction(lat, u, bt)
      sm += s
      s2 += s*s
      let th = plaqSpatial(lat, u, 0, 0)
      tm += th*th
      t2 += th*th*th*th
    let
      fn = float(nsamp)
      sbar = sm/fn
      serr = sqrt(max(0.0, s2/fn - sbar*sbar)/(fn - 1.0))
      tbar = tm/fn
      terr = sqrt(max(0.0, t2/fn - tbar*tbar)/(fn - 1.0))
    echo &"  <S> = {sbar:.4f} +- {serr:.4f}   rank/2 = {0.5*float(rank):.4f}" &
         &"   pull {(sbar - 0.5*float(rank))/serr:+.2f}   ({its div nsamp} CG its/sample)"
    echo &"  <Theta^2> = {tbar:.6f} +- {terr:.6f}   exact = {exact:.6f}" &
         &"   pull {(tbar - exact)/terr:+.2f}"
    check abs(sbar - 0.5*float(rank)) < 5.0*serr
    check abs(tbar - exact) < 5.0*terr

  test "the heatbath field is transverse and the Polyakov mode is zero":
    var
      r = rng(52)
      u = newGauge(lat)
    discard heatbath(lat, u, bet, r, 1e-22, 20000)
    var d = newSeq[float](lat.nsite)
    divergence(lat, d, u)
    var n = 0.0
    for v in d: n += v*v
    var s = 0.0
    for x in u.t: s += x
    echo &"  |div u|^2 = {n:.3e}   sum theta^t = {s:.3e}   |u|^2 = {norm2(u):.4f}"
    check n < 1e-16*norm2(u)
    check abs(s) < 1e-8*sqrt(norm2 u)

suite "gauge action: geometry conventions":

  test "the three conventions differ at O(abar^2) only":
    for lv in [1, 2, 4]:
      let
        s2 = newSphere(lv)
        l2 = newLat(s2, 4, 0.35)
        bs = newBeta(l2, 1.0, gcGeodesic)   # explicit: the default is now gcExactArea
        bx = newBeta(l2, 1.0, gcExactArea)
        bf = newBeta(l2, 1.0, gcFlat)
      var ef, ee, ex = 0.0
      for f in 0..<s2.nf: ef = max(ef, abs(bf.face[f]/bs.face[f] - 1.0))
      for e in 0..<s2.ne:
        ee = max(ee, abs(bf.edge[e]/bs.edge[e] - 1.0))
        ex = max(ex, abs(bx.edge[e]/bs.edge[e] - 1.0))
      echo &"  L={lv} abar^2={s2.abar*s2.abar:.5f}  max|beta_tri^flat/beta_tri-1|={ef:.4e}" &
           &"  edge flat {ee:.4e}  edge exact-area {ex:.4e}"
      check ef < 0.6*s2.abar*s2.abar
      check ee < 0.6*s2.abar*s2.abar
      check ex < 0.6*s2.abar*s2.abar
      # gcExactArea only changes beta_l
      for f in 0..<s2.nf: check bx.face[f] == bs.face[f]

  test "the exact spherical diamond areas tile the sphere":
    ## This is why gcExactArea exists: sum_e A_e^exact = sum_f A_tri = 4 pi to
    ## roundoff, whereas the flat form sum_e l(l*_1+l*_2)/2 misses by O(abar^2)
    ## (3.6 % at L=1).  Recovered here from beta_l = 2 A_e/(g2 l^2 at) with
    ## g2 = 1, at = 0.35.  Identifying THIS as the paper's A_l is what reproduces
    ## the published Delta_0(L=1) = 1.33242; see the WP-G entry in doc/06.
    for lv in [1, 2, 4]:
      let
        s2 = newSphere(lv)
        l2 = newLat(s2, 4, 0.35)
        bs = newBeta(l2, 1.0, gcGeodesic)   # explicit: the default is now gcExactArea
        bx = newBeta(l2, 1.0, gcExactArea)
      var sx, ss, sf = 0.0
      for e in 0..<s2.ne:
        let w = 0.5*0.35*s2.edges[e].len*s2.edges[e].len
        sx += w*bx.edge[e]
        ss += w*bs.edge[e]
      for f in 0..<s2.nf: sf += s2.faces[f].area
      echo &"  L={lv}: sum A_e exact = {sx:.12f}  flat-form = {ss:.12f}  " &
           &"sum A_tri = {sf:.12f}  4pi = {4.0*PI:.12f}"
      check abs(sx - 4.0*PI) < 1e-12
      check abs(sf - 4.0*PI) < 1e-12
      check abs(ss - 4.0*PI) > 1e-3*s2.abar*s2.abar    # the flat form does NOT

  test "the flat action is still gauge invariant and its force is exact":
    let bf = newBeta(lat, g2, gcFlat)
    var r = rng(61)
    let
      u = randGauge(lat, r)
      al = randSeqF(lat.nsite, r)
    var p = newGauge(lat)
    gradient(lat, p, al)
    var v = newGauge(lat)
    v := u
    axpy(v, 1.0, p)
    check abs(gaugeAction(lat, v, bf) - gaugeAction(lat, u, bf)) <
          1e-12*gaugeAction(lat, u, bf)
    let df = newDense(lat, bf)
    var f = newGauge(lat)
    gaugeForce(lat, f, u, bf)
    let
      fa = f.toSeq
      fd = df.mul(u.toSeq)
    var e = 0.0
    for i in 0..<fa.len: e = max(e, abs(fa[i] - fd[i]))
    check e < 1e-12*maxAbs(fd)
