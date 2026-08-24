#RUNCMD $RUN1
## WP-C: analytic continuum correlators (doc/03-targets T1.6, doc/02-formulation section 7).
## Also writes the Fig. 13 / Fig. 14 data of arXiv:2510.03085 to output/radial/analytic/.

import std/[math, os, strformat, unittest]
import ../core/analytic

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))

const
  outDir = currentSourcePath().parentDir.parentDir.parentDir.parentDir.parentDir /
           "output" / "radial" / "analytic"
  sqrt2 = sqrt(2.0)

# --- independent oracles ----------------------------------------------------

func fermionThermal(t, T: float, nmax: int): float =
  ## Mode-by-mode antiperiodic finite-T sum, independent of the image construction:
  ##   (1/4pi) sum_n (n+1) [e^{-(n+1)t} + e^{-(n+1)(T-t)}] / (1 + e^{-(n+1)T}),  0 <= t <= T.
  var
    s = 0.0
    n = 0
  while true:
    let
      lam = (n + 1).float
      term = lam*(exp(-lam*t) + exp(-lam*(T - t)))/(1.0 + exp(-lam*T))
    s += term
    if (if nmax >= 0: n >= nmax else: n > 2 and term <= 1e-18*s): break
    inc n
  s/(4.0*PI)

func gaugeThermal(t, T: float, nmax: int): float =
  ## Same for the periodic gauge tower, 1/(1 - e^{-lam T}).
  var
    s = 0.0
    n = 1
  while true:
    let
      lam = sqrt(n.float*(n + 1).float)
      term = lam*(2*n + 1).float*(exp(-lam*t) + exp(-lam*(T - t)))/(1.0 - exp(-lam*T))
    s += term
    if (if nmax >= 0: n >= nmax else: n > 2 and term <= 1e-18*s): break
    inc n
  s/(8.0*PI)

proc gaussLegendre(m: int): tuple[x, w: seq[float]] =
  ## Nodes and weights on [-1,1], exact for polynomials up to degree 2m-1.
  ## Newton on P_m with P_m' = ((m+1)/2) P_{m-1}^(1,1); w = 2/((1-x^2) P_m'^2).
  result.x = newSeq[float](m)
  result.w = newSeq[float](m)
  for i in 0..<m:
    var z = cos(PI*(i.float + 0.75)/(m.float + 0.5))
    for it in 0..99:
      let
        p = jacobiP(m, 0.0, 0.0, z)
        dp = 0.5*(m + 1).float*jacobiP(m - 1, 1.0, 1.0, z)
        dz = p/dp
      z -= dz
      if abs(dz) < 1e-16: break
    let dp = 0.5*(m + 1).float*jacobiP(m - 1, 1.0, 1.0, z)
    result.x[i] = z
    result.w[i] = 2.0/((1.0 - z*z)*dp*dp)

let gl = gaussLegendre(96)

proc glInt(f: proc(z: float): float): float =
  for i in 0..<gl.x.len: result += gl.w[i]*f(gl.x[i])

# --- tests ------------------------------------------------------------------

suite "jacobiP":

  test "closed forms for small j":
    ## P_0 = 1, P_1 = (a+1) + (a+b+2)(z-1)/2, and P_2 from the terminating series of (C.11).
    for (a, b) in [(0.0, 0.0), (0.0, -1.0), (1.0, -1.0), (2.0, 0.0), (0.5, 1.5), (-0.5, 0.25)]:
      for z in [-0.97, -0.4, 0.0, 0.3, 0.88]:
        let
          x = 0.5*(1.0 - z)
          p1 = a + 1.0 + 0.5*(a + b + 2.0)*(z - 1.0)
          p2 = 0.5*(a + 1.0)*(a + 2.0)*
               (1.0 - 2.0*(3.0 + a + b)/(a + 1.0)*x +
                (3.0 + a + b)*(4.0 + a + b)/((a + 1.0)*(a + 2.0))*x*x)
        check abs(jacobiP(0, a, b, z) - 1.0) < 1e-14
        check abs(jacobiP(1, a, b, z) - p1) < 1e-14
        check abs(jacobiP(2, a, b, z) - p2) < 1e-13

  test "Legendre special case and Rodrigues values at z = 1":
    for l in 0..12:
      # P_l(z) = jacobiP(l,0,0,z); check against the Bonnet recurrence
      var p0 = 1.0
      var p1 = 0.37
      for n in 2..l:
        let p = ((2*n - 1).float*0.37*p1 - (n - 1).float*p0)/n.float
        p0 = p1
        p1 = p
      let want = if l == 0: 1.0 else: p1
      check abs(jacobiP(l, 0.0, 0.0, 0.37) - want) < 1e-13
      # P_j^(a,b)(1) = binom(j+a, j)
      for (a, b) in [(0.0, -1.0), (1.0, -1.0), (2.0, 0.0)]:
        var bn = 1.0
        for k in 1..l: bn *= (a + k.float)/k.float
        check abs(jacobiP(l, a, b, 1.0) - bn) < 1e-11*max(1.0, bn)

  test "degenerate beta = -1 reductions used by (C.10) and (C.43)":
    for n in 0..40:
      for z in [-0.999, -0.6, -0.1, 0.45, 0.93, 1.0]:
        # P_{n+1}^(0,-1)(z) = ((1+z)/2) P_n^(0,1)(z)   -> xi_{1/2,n} = (1/2) sqrt(1+z) P_n^(0,1)
        let
          lhs = jacobiP(n + 1, 0.0, -1.0, z)
          rhs = 0.5*(1.0 + z)*jacobiP(n, 0.0, 1.0, z)
        check abs(lhs - rhs) <= 1e-11*max(1.0, abs(rhs))
    for n in 1..60:
      for z in [-0.999, -0.6, -0.1, 0.45, 0.93]:
        # f_{0,n} = (1-z) P_n^(1,-1) = (1-z^2) P_n'/n, so f'_{0,n} = -(n+1) P_n by Legendre
        let
          f = (1.0 - z)*jacobiP(n, 1.0, -1.0, z)
          g = (1.0 - z*z)*0.5*(n + 1).float*jacobiP(n - 1, 1.0, 1.0, z)/n.float
        check abs(f - g) <= 1e-10*max(1.0, abs(g))

suite "fermion propagator (V.3)":

  test "closed form == truncated sum as nmax -> infinity":
    for t in [0.5, 1.0, 2.0, 4.0]:
      let ex = fermionG(t)
      var prev = 1.0
      for nmax in [20, 60, 200]:
        let rel = abs(fermionG(t, nmax) - ex)/ex
        if nmax == 200:
          check rel < 1e-14
        check rel <= prev + 1e-18
        prev = rel
      echo &"  t={t:5.2f} G={ex:.16e} rel(nmax=200)={abs(fermionG(t,200)-ex)/ex:.3e}"

  test "short-distance normalization G(t) 4 pi t^2 -> 1":
    var prev = 1.0
    for t in [1e-1, 1e-2, 1e-3, 1e-4]:
      let d = abs(fermionG(t)*4.0*PI*t*t - 1.0)
      echo &"  t={t:8.1e} |G 4pi t^2 - 1| = {d:.4e}   (expect t^2/12 = {t*t/12.0:.4e})"
      check d < 0.11*t*t          # leading correction is t^2/12
      check d < prev
      prev = d
    # sign(t) and the 1/(16 pi sinh^2) form
    check abs(fermionG(-1.3) + fermionG(1.3)) < 1e-18
    check abs(fermionG(1.3) - 1.0/(16.0*PI*sinh(0.65)^2)) < 1e-18

  test "antiperiodic images: symmetry, antiperiodicity, mode-sum equality":
    ## `T - t` cancels to ~eps*T, and h ~ 1/u^2 doubles that, so the identities hold to
    ## ~eps*T/t, not to eps.  1e-13 covers t = 0.35 at T = 12.
    const T = 12.0
    for t in [0.35, 1.0, 3.0, 5.5, 6.0]:
      let g = fermionGPeriodic(t, T)
      check abs(fermionGPeriodic(T - t, T) - g) <= 1e-13*abs(g)      # fold about T/2
      check abs(fermionGPeriodic(t + T, T) + g) <= 1e-13*abs(g)      # antiperiodic
      check abs(fermionGPeriodic(-t, T) + g) <= 1e-13*abs(g)         # odd in t
      let m = fermionThermal(t, T, -1)
      check abs(g - m) <= 1e-14*abs(m)
    for nmax in [5, 30, 120]:
      for t in [0.35, 2.0, 5.0]:
        let
          g = fermionGPeriodic(t, T, nmax)
          m = fermionThermal(t, T, nmax)
        check abs(g - m) <= 1e-14*abs(m)

suite "gauge current correlator (V.14)":

  test "truncation and periodic images":
    check abs(gaugeG(2.0, -1) - gaugeG(2.0, 400)) <= 1e-15*gaugeG(2.0, -1)
    check abs(gaugeG(-1.7) - gaugeG(1.7)) < 1e-18
    const T = 12.0
    for t in [0.4, 1.0, 3.0, 5.0, 6.0]:
      let g = gaugeGPeriodic(t, T)
      check abs(gaugeGPeriodic(T - t, T) - g) <= 1e-14*g
      check abs(gaugeGPeriodic(t + T, T) - g) <= 1e-14*g
      let m = gaugeThermal(t, T, 200)
      check abs(g - m) <= 1e-13*m
    # leading exponent is sqrt(2), not an integer; the n=2 contamination at t=21 is ~3e-9
    let r = ln(gaugeG(20.0, -1)/gaugeG(22.0, -1))/2.0
    echo &"  -d ln G_g/dt at t=21: {r:.12f}  (sqrt2 = {sqrt2:.12f})"
    check abs(r - sqrt2) < 1e-8

suite "effective dimension (V.4)-(V.5)":

  test "exact cosh recovers Delta":
    const
      T = 16.0
      at = 0.1
      dl = 1.234567
    let nt = int(round(T/at))
    var g = newSeq[float](nt)
    for i in 0..<nt: g[i] = cosh(dl*(0.5*T - i.float*at))
    let d = effDim(g, at, T)
    var worst = 0.0
    for i in 0..<int(round((0.5*T - 1.0)/at)): worst = max(worst, abs(d[i] - dl))
    echo &"  max |Delta_eff - Delta| for t < T/2-1 : {worst:.3e}"
    check worst < 1e-12

  test "fermion plateau -> 1 and gauge plateau -> sqrt(2)":
    const
      T = 40.0
      at = 0.1
    let nt = int(round(T/at))
    var
      gf = newSeq[float](nt)
      gg = newSeq[float](nt)
    for i in 1..<nt:
      let t = i.float*at
      gf[i] = fermionGPeriodic(t, T)
      gg[i] = gaugeGPeriodic(t, T)
    gf[0] = gf[1]                      # t = 0 is the contact singularity; never used below
    gg[0] = gg[1]
    let
      df = effDim(gf, at, T)
      dg = effDim(gg, at, T)
    var
      wf = 0.0
      wg = 0.0
    for i in int(round(16.0/at))..int(round(19.0/at)):
      wf = max(wf, abs(df[i] - 1.0))
      wg = max(wg, abs(dg[i] - sqrt2))
    for tt in [8.0, 12.0, 16.0, 19.0]:
      let i = int(round(tt/at))
      echo &"  t={tt:5.1f}  Delta_eff(fermion)={df[i]:.12f}  Delta_eff(gauge)={dg[i]:.12f}"
    echo &"  plateau 16 <= t <= 19: fermion dev {wf:.3e}, gauge dev {wg:.3e}"
    check wf < 1e-6
    check wg < 1e-6

suite "flat Wilson spectrum (IV.8), Fig. 5":

  test "range, doubler position and conjugate pairing":
    let kap = 1.0/sqrt(3.0)
    for aOverAt in [2.0, 1.0]:                     # 2 kap' above / below 4 kap
      let
        kapT = 0.5*sqrt(3.0)*aOverAt
        ev = flatSpectrum(kap, kapT, 12, 12, 12)
      check ev.len == 2*12*12*12
      var
        remin = ev[0].re
        remax = ev[0].re
        dmin = 4.5*kap + 2.0*kapT
      for i in countup(0, ev.len - 2, 2):
        check abs(ev[i].re - ev[i + 1].re) < 1e-18       # conjugate pair
        check abs(ev[i].im + ev[i + 1].im) < 1e-18
        check ev[i].im >= 0.0
        remin = min(remin, ev[i].re)
        remax = max(remax, ev[i].re)
        if ev[i].im < 1e-10 and ev[i].re > 1e-10:
          dmin = min(dmin, ev[i].re)
      let
        want = min(4.0*kap, 2.0*kapT)
        top = 4.5*kap + 2.0*kapT       # K point (k1,k2) = (2pi/3, 2pi/3), not 4 kap + 2 kap'
      echo &"  abar/at={aOverAt}: Re in [{remin:.12f}, {remax:.12f}], top={top:.12f}, " &
           &"first doubler {dmin:.12f} (min(4kap,2kap')={want:.12f})"
      check remin > -1e-15
      check abs(remin) < 1e-15
      check remax <= top + 1e-12
      check abs(remax - top) < 1e-12
      check abs(dmin - want) < 1e-12

suite "T1.6a: S^2 fermion propagator (C.54)-(C.55)":

  test "xi_{1/2,n} sum reproduces the CFT answer in envelope":
    ## The series is conditionally convergent (Abel/Cesaro), so partial sums oscillate with an
    ## O(nmax^-1/2) envelope.  What must decrease is the envelope, taken over one oscillation
    ## period 2 pi / theta in n, and the bound |err| sin(theta/2) sqrt(nmax).
    for theta in [0.5, 1.0, 2.0]:
      let
        ex = 1.0/(4.0*PI*sin(0.5*theta))
        per = int(ceil(2.0*PI/theta))
      var prev = 1e30
      for nmax in [50, 100, 200]:
        var env = 0.0
        for n in nmax..nmax + per: env = max(env, abs(s2FermionProp(theta, n) - ex))
        echo &"  theta={theta:4.2f} nmax={nmax:4d} exact={ex:.8f} " &
             &"pointwise={s2FermionProp(theta, nmax):.8f} envelope-rel={env/ex:.4e} " &
             &"C=|err| sin(th/2) sqrt(n)={env*sin(0.5*theta)*sqrt(nmax.float):.4f}"
        check env < prev                                        # envelope strictly decreasing
        check env*sin(0.5*theta)*sqrt(nmax.float) < 0.15        # uniform bound
        prev = env

  test "ultraviolet is restored down to theta ~ pi/nmax":
    for theta in [0.05, 0.02]:
      let
        ex = 1.0/(4.0*PI*sin(0.5*theta))
        per = int(ceil(2.0*PI/theta))
      var env: array[2, float]
      for k, nmax in [50, 200]:
        var e = 0.0
        for n in nmax..nmax + per: e = max(e, abs(s2FermionProp(theta, n) - ex))
        env[k] = e/ex
      echo &"  theta={theta:5.3f} envelope-rel: nmax=50 {env[0]:.4e}, nmax=200 {env[1]:.4e}"
      check env[1] < env[0]
    # far below pi/nmax the truncated sum is far short of the exact 1/(4 pi sin(theta/2))
    check s2FermionProp(1e-3, 50) < 0.02*1.0/(4.0*PI*sin(0.5e-3))

  test "literal (C.10) equals the |m| = 1/2 reduction":
    for theta in [0.05, 0.5, 1.0, 2.0, 3.0]:
      let z = cos(theta)
      for nmax in [3, 17, 60]:
        var s = 0.0
        for n in 0..nmax:
          let x = 0.5*sqrt(1.0 - z)*jacobiP(n, 0.0, 1.0, -z)   # (1/2) sqrt(1+w) P_n^(0,1)(w)
          s += (if (n and 1) == 0: x else: -x)
        s /= PI*sqrt2
        check abs(s2FermionProp(theta, nmax) - s) <= 1e-10*max(1e-3, abs(s))

suite "T1.6b: S^2 current correlator (C.56)-(C.57)":

  test "literal (C.57) collapses to the Legendre completeness sum":
    for nmax in [1, 5, 20, 40]:
      for z in [-0.999, -0.5, 0.0, 0.5, 0.9, 0.999, 1.0]:
        var s = 0.0
        for n in 1..nmax: s += (2*n + 1).float*jacobiP(n, 0.0, 0.0, z)
        s /= 4.0*PI
        check abs(s2CurrentCorr(z, nmax) - s) <= 1e-11*max(1.0, abs(s))

  test "exact polynomial moments: delta normalization and the -1/(4 pi) offset":
    ## For a polynomial phi of degree <= nmax the truncated sum integrates exactly to
    ## (1/2pi)[phi(1) - (1/2) int phi], so these pin both pieces of (C.57) with no truncation.
    for nmax in [5, 20, 40]:
      let
        i0 = glInt(proc(z: float): float = s2CurrentCorr(z, nmax))
        i1 = glInt(proc(z: float): float = z*s2CurrentCorr(z, nmax))
        i2 = glInt(proc(z: float): float = z*z*s2CurrentCorr(z, nmax))
        io = glInt(proc(z: float): float = (1.0 - z)*s2CurrentCorr(z, nmax))
      echo &"  nmax={nmax:3d}  <1>={i0:.3e}  <z>={i1:.15f}  <z^2>={i2:.15f}  <1-z>={io:.15f}"
      check abs(i0) < 1e-13                       # (1/2pi)[1 - (1/2)(2)]
      check abs(i1 - 1.0/(2.0*PI)) < 1e-13        # (1/2pi)[1 - 0]
      check abs(i2 - 1.0/(3.0*PI)) < 1e-13        # (1/2pi)[1 - 1/3]
      # phi(1) = 0 kills the delta, so <1-z> is purely the offset: -(1/4pi) int (1-z) dz
      check abs(io + 1.0/(2.0*PI)) < 1e-13
      check abs(io - (-1.0/(4.0*PI))*2.0) < 1e-13

  test "delta function develops against smooth test functions, nmax 20 -> 40":
    ## A smeared delta of width 0.15 and a function with a pole just outside [-1,1] both
    ## need modes up to n ~ 1/width; exp(z) is entire and is already exact at nmax = 20.
    proc bump(z: float): float =
      let d = (z - 1.0)/0.15
      exp(-0.5*d*d)
    proc pole(z: float): float = 1.0/(1.05 - z)
    proc ent(z: float): float = exp(z)
    for (name, fn, tol) in [("gauss(z-1, 0.15)", bump, 1e-7),
                            ("1/(1.05 - z)", pole, 1e-4),
                            ("exp(z)", ent, 1e-13)]:
      let ex = (fn(1.0) - 0.5*glInt(fn))/(2.0*PI)
      var e: array[2, float]
      for k, nmax in [20, 40]:
        let v = glInt(proc(z: float): float = fn(z)*s2CurrentCorr(z, nmax))
        e[k] = abs(v - ex)
      echo &"  {name:18s} exact={ex:.10f}  err(20)={e[0]:.3e}  err(40)={e[1]:.3e}"
      check e[1] < tol
      check e[1] <= max(1e-2*e[0], 1e-13)   # 100x gain, or already at the roundoff floor

  test "Gauss-Legendre oracle is exact":
    for k in 0..40:
      let
        v = glInt(proc(z: float): float = z^k)
        w = if (k and 1) == 1: 0.0 else: 2.0/(k + 1).float
      check abs(v - w) < 1e-13

suite "figure data":

  test "write Fig. 13 and Fig. 14 tables":
    createDir(outDir)
    block:
      let f = open(outDir/"fig13_s2_fermion.tsv", fmWrite)
      f.writeLine("# arXiv:2510.03085 Fig. 13: free fermion propagator on S^2")
      f.writeLine("# columns: theta  exact_C54  nmax50_C55  nmax200_C55")
      const np = 400
      for i in 0..np:
        let theta = exp(ln(1e-3) + (ln(3.0) - ln(1e-3))*i.float/np.float)
        f.writeLine(&"{theta:.10e}\t{1.0/(4.0*PI*sin(0.5*theta)):.10e}\t" &
                    &"{s2FermionProp(theta, 50):.10e}\t{s2FermionProp(theta, 200):.10e}")
      f.close()
    block:
      let f = open(outDir/"fig14_s2_current.tsv", fmWrite)
      f.writeLine("# arXiv:2510.03085 Fig. 14: J^t J^t correlator on S^2, (1/g^2)<J^t J^t>")
      f.writeLine("# columns: z  nmax20_C57  nmax40_C57")
      const np = 2000
      for i in 0..np:
        let z = -1.0 + 2.0*i.float/np.float
        f.writeLine(&"{z:.10e}\t{s2CurrentCorr(z, 20):.10e}\t{s2CurrentCorr(z, 40):.10e}")
      f.close()
    check fileExists(outDir/"fig13_s2_fermion.tsv")
    check fileExists(outDir/"fig14_s2_current.tsv")
    echo "  wrote ", outDir
