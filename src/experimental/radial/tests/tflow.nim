#RUNCMD $RUN1
## WP-G: gradient flow of the non-compact U(1) gauge field.  doc/04-interfaces.md
## section 12, doc/03-targets.md T2.2.
##
## The action is quadratic, so the flow is LINEAR: theta(s) = exp(-M s) theta(0).
## The oracle is therefore exact -- a dense matrix exponential through the
## eigendecomposition of the same incidence matrix that `tgauge.nim` uses.

import std/[math, strformat, unittest]
import eigens/linalgFuncs
import ../ops/flow

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))

# --- dense oracle -----------------------------------------------------------

type Dense = object
  n: int
  m: seq[float]

proc newDense(l: Lat, b: Beta): Dense =
  let
    n = nlink(l)
    ns = l.sph.ne*l.nt
  var m = newSeq[float](n*n)
  var row = newSeq[float](n)
  var nz: seq[int]
  for kind in 0..1:
    for idx in 0..<(if kind == 0: l.sph.nf else: l.sph.ne):
      let w = if kind == 0: b.face[idx] else: b.edge[idx]
      for t in 0..<l.nt:
        for i in 0..<n: row[i] = 0.0
        if kind == 0:
          for i in 0..2:
            row[eIdx(l, l.sph.faces[idx].e[i], t)] += float(l.sph.faces[idx].s[i])
        else:
          row[eIdx(l, idx, t)] += 1.0
          row[ns + tIdx(l, l.sph.edges[idx].b, t)] += 1.0
          row[eIdx(l, idx, t+1)] -= 1.0
          row[ns + tIdx(l, l.sph.edges[idx].a, t)] -= 1.0
        nz.setLen 0
        for i in 0..<n:
          if row[i] != 0.0: nz.add i
        for i in nz:
          for j in nz: m[i + n*j] += w*row[i]*row[j]
  Dense(n: n, m: m)

proc eigen(d: Dense): tuple[v: seq[Complex64], w: seq[float]] =
  var h = newSeq[Complex64](d.n*d.n)
  for i in 0..<d.n*d.n: h[i] = complex64(d.m[i], 0.0)
  var ew = newSeq[float](d.n)
  zeigs(cast[ptr float64](addr h[0]), addr ew[0], d.n)
  (h, ew)

proc expmv(d: Dense, ev: tuple[v: seq[Complex64], w: seq[float]],
           s: float, x: openArray[float]): seq[float] =
  ## exp(-M s) x, formed mode by mode so no dense matrix is ever squared.
  result = newSeq[float](d.n)
  for k in 0..<d.n:
    var c = 0.0
    for i in 0..<d.n: c += ev.v[i + d.n*k].re*x[i]
    c *= exp(-s*ev.w[k])
    for i in 0..<d.n: result[i] += c*ev.v[i + d.n*k].re

# --- helpers ----------------------------------------------------------------

proc rng(sed: int): Threefry4x64 =
  result.seedIndep(sed, 0)

proc randGauge(l: Lat, r: var Threefry4x64): Gauge =
  result = newGauge(l)
  for i in 0..<result.s.len: result.s[i] = r.gaussian
  for i in 0..<result.t.len: result.t[i] = r.gaussian

proc randSeqF(n: int, r: var Threefry4x64): seq[float] =
  result = newSeq[float](n)
  for i in 0..<n: result[i] = r.gaussian

proc maxDiff(a, b: openArray[float]): float =
  for i in 0..<a.len: result = max(result, abs(a[i] - b[i]))

proc maxAbs(x: openArray[float]): float =
  for v in x: result = max(result, abs v)

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

suite "gradient flow":

  test "spectrum of M and the stability window":
    echo &"  lambda(M) in [{maxAbs(ev.w):.4f} max];  ",
         &"nlink = {den.n}, at = {at}, g2 = {g2}"
    check maxAbs(ev.w) > 0.0

  test "RK flow reproduces exp(-M s) theta_0":
    var r = rng(71)
    let u0 = randGauge(lat, r)
    for s in [0.02, 0.1, 0.5]:
      var u = newGauge(lat)
      u := u0
      var got: seq[float]
      flowRun(lat, u, bet, [s], s/512.0, RK4CK_2N,
              proc(t: float, v: Gauge) = got = v.toSeq)
      let want = den.expmv(ev, s, u0.toSeq)
      let e = maxDiff(got, want)
      echo &"  s = {s:.2f}: max |RK - exp(-Ms)theta| = {e:.3e}  (|theta(s)| = {maxAbs(want):.4f})"
      check e < 1e-10

  test "RK order from halving the step":
    ## The asymptotic regime needs h*lambda_max well below 1: lambda_max is 12.4
    ## here, so the n = 2 step (h*lambda = 2.5) is outside it and reports a
    ## spuriously high order.  Start at n = 16, h*lambda = 0.31.
    var r = rng(72)
    let
      u0 = randGauge(lat, r)
      s = 0.4
      want = den.expmv(ev, s, u0.toSeq)
    for (name, order) in [("RK3W6", 3.0), ("RK4CK", 4.0)]:
      var es: seq[float]
      for n in [16, 32, 64, 128]:
        var u = newGauge(lat)
        u := u0
        var got: seq[float]
        let cb = proc(t: float, v: Gauge) = got = v.toSeq
        if name == "RK3W6": flowRun(lat, u, bet, [s], s/float(n) + 1e-12, RK3W6_2N, cb)
        else: flowRun(lat, u, bet, [s], s/float(n) + 1e-12, RK4CK_2N, cb)
        es.add maxDiff(got, want)
      var p: seq[float]
      for i in 1..<es.len: p.add log2(es[i-1]/es[i])
      echo &"  {name}: errs " &
           &"{es[0]:.3e} {es[1]:.3e} {es[2]:.3e} {es[3]:.3e}  " &
           &"orders {p[0]:.2f} {p[1]:.2f} {p[2]:.2f}"
      for x in p: check abs(x - order) < 0.25

  test "the flow is linear":
    var r = rng(73)
    let u0 = randGauge(lat, r)
    var a = newGauge(lat)
    a := u0
    scale(a, 3.0)
    var b = newGauge(lat)
    b := u0
    flowRun(lat, a, bet, [0.3], 0.01, RK4CK_2N, nil)
    flowRun(lat, b, bet, [0.3], 0.01, RK4CK_2N, nil)
    scale(b, 3.0)
    let e = maxDiff(a.toSeq, b.toSeq)
    echo &"  max |flow(3u) - 3 flow(u)| = {e:.3e}"
    check e < 1e-13*maxAbs(a.toSeq)

  test "the flow is gauge covariant":
    var r = rng(74)
    let
      u0 = randGauge(lat, r)
      al = randSeqF(lat.nsite, r)
    var k = newGauge(lat)
    gradient(lat, k, al)
    var a = newGauge(lat)
    a := u0
    var b = newGauge(lat)
    b := u0
    axpy(b, 1.0, k)
    flowRun(lat, a, bet, [0.25], 0.005, RK4CK_2N, nil)
    flowRun(lat, b, bet, [0.25], 0.005, RK4CK_2N, nil)
    # exp(-Ms) leaves ker M pointwise fixed, so the flowed fields differ by d alpha
    axpy(b, -1.0, a)
    let e = maxDiff(b.toSeq, k.toSeq)
    echo &"  max |flow(u + d a) - flow(u) - d a| = {e:.3e}"
    check e < 1e-12*maxAbs(k.toSeq)
    # and every gauge-invariant observable is untouched
    var eo = 0.0
    for f in 0..<sph.nf:
      for t in 0..<nt:
        var c = newGauge(lat)
        c := a
        axpy(c, 1.0, k)
        eo = max(eo, abs(plaqSpatial(lat, c, f, t) - plaqSpatial(lat, a, f, t)))
    check eo < 1e-12

  test "energyDensity decreases monotonically along the flow":
    var
      r = rng(75)
      u = newGauge(lat)
    discard heatbath(lat, u, bet, r, 1e-22, 20000)
    var
      prevS = 1e300
      prevE = 1e300
      prevT = 1e300
      okS = true
      okE = true
      okT = true
      rows: seq[string]
    flowRun(lat, u, bet, [0.0, 0.02, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6], 0.005, RK4CK_2N,
      proc(t: float, v: Gauge) =
        let
          s = gaugeAction(lat, v, bet)
          e = energyDensity(lat, v, bet)
          w = energyDensityT(lat, v, bet)
        rows.add &"    s={t:5.2f}  S={s:12.6f}  E_s={e:12.6f}  E_t={w:12.6f}"
        if s > prevS: okS = false
        if e > prevE: okE = false
        if w > prevT: okT = false
        prevS = s
        prevE = e
        prevT = w)
    for l in rows: echo l
    check okS
    check okE
    check okT

  test "flowStep matches one rk2nStep of flowRun":
    var r = rng(76)
    let u0 = randGauge(lat, r)
    var a = newGauge(lat)
    a := u0
    flowStep(lat, a, bet, 0.01, RK4CK_2N)
    var b = newGauge(lat)
    b := u0
    flowRun(lat, b, bet, [0.01], 0.01, RK4CK_2N, nil)
    check maxDiff(a.toSeq, b.toSeq) == 0.0

  test "free-Maxwell short-flow-time law: E_s ~ s^{-3/2}":
    ## <E_s(s)> = (1/2V) sum_i c_i exp(-2 lambda_i s) exactly, so the ensemble
    ## average is a pure trace -- no Monte Carlo needed.  On a coarse lattice the
    ## continuum t^{-3/2} only shows up in a window; this test just pins the exact
    ## trace formula against the flowed heatbath ensemble.
    let nsamp = 200
    var
      r = rng(77)
      u = newGauge(lat)
      acc = 0.0
    for n in 0..<nsamp:
      discard heatbath(lat, u, bet, r, 1e-22, 20000)
      flowRun(lat, u, bet, [0.05], 0.002, RK4CK_2N, nil)
      acc += energyDensity(lat, u, bet)
    acc /= float(nsamp)
    # exact: E_s = tr(M_s exp(-2 M s) M^+)/(2 V),  M_s = the spatial part of M
    let
      bs = Beta(face: bet.face, edge: newSeq[float](sph.ne), afac: bet.afac,
                g2: bet.g2, conv: bet.conv)
      ds = newDense(lat, bs)
      tolw = 1e-9*maxAbs(ev.w)
    var tr = 0.0
    for k in 0..<den.n:
      if abs(ev.w[k]) <= tolw: continue
      var col = newSeq[float](den.n)
      for i in 0..<den.n: col[i] = ev.v[i + den.n*k].re
      var q = 0.0
      for i in 0..<den.n:
        var s = 0.0
        for j in 0..<den.n: s += ds.m[i + den.n*j]*col[j]
        q += col[i]*s
      tr += q*exp(-2.0*0.05*ev.w[k])/ev.w[k]
    let exact = 0.5*tr/(4.0*PI*at*float(nt))
    echo &"  <E_s(0.05)> = {acc:.6f} over {nsamp} heatbath samples,  exact trace = {exact:.6f}"
    check abs(acc - exact) < 0.05*exact
