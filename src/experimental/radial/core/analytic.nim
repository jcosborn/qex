## Analytic continuum correlators for QED3 on S^2 x R (free limit) and the
## two-dimensional S^2 checks.  Equation tags are those of arXiv:2510.03085.
##
## Sphere radius R = 1 throughout, so t is in units of R and z = cos(theta).
##
## Truncation convention, used by every mode sum here: `nmax >= 0` keeps exactly the
## modes up to n = nmax -- that is what the n_max fits of the paper compare against.
## `nmax < 0` means the untruncated correlator: a closed form where one exists,
## otherwise the sum carried until the terms underflow.

import std/math
import types

const
  sumEps = 1e-17        ## relative term size at which a mode / image sum stops
  maxMode = 1_000_000   ## mode cap, only reached if the sum does not converge
  maxImage = 1024       ## image cap, likewise

# --- Jacobi polynomials -----------------------------------------------------

func jacobiP*(j: int, alpha, beta, z: float): float =
  ## Jacobi polynomial with generalized (real, possibly negative) parameters, (C.11):
  ##   P_j^(a,b)(z) = [Gamma(a+1+j)/(Gamma(a+1) j!)] F(-j, 1+a+b+j; a+1; (1-z)/2),  j >= 0.
  ## Evaluated by the three-term recurrence (DLMF 18.9.2)
  ##   2n(n+a+b)(2n+a+b-2) P_n = (2n+a+b-1)[(2n+a+b)(2n+a+b-2) z + a^2 - b^2] P_{n-1}
  ##                             - 2(n+a-1)(n+b-1)(2n+a+b) P_{n-2}
  ## seeded with P_0 = 1, P_1 = (a+1) + (a+b+2)(z-1)/2.  The terminating series of (C.11)
  ## is deliberately not used: at j ~ 200 its terms cancel over tens of digits.
  ##
  ## The leading coefficient vanishes only for a+b = -n or a+b = 2-2n.  For the parameter
  ## sets needed here -- (0,-1) for xi_{1/2,n}, (1,-1) and (2,0) for f_{0,n}, (0,0) for
  ## Legendre -- those roots sit at n = 0, 1 only, which the explicit seeds cover.
  if j == 0: return 1.0
  let ab = alpha + beta
  var
    p0 = 1.0
    p1 = alpha + 1.0 + 0.5*(ab + 2.0)*(z - 1.0)
  for n in 2..j:
    let
      fn = n.float
      d = 2.0*fn + ab
      c0 = 2.0*fn*(fn + ab)*(d - 2.0)
      c1 = (d - 1.0)*(d*(d - 2.0)*z + alpha*alpha - beta*beta)
      c2 = 2.0*(fn + alpha - 1.0)*(fn + beta - 1.0)*d
      p = (c1*p1 - c2*p0)/c0
    p0 = p1
    p1 = p
  p1

# --- Fermion propagator on S^2 x R ------------------------------------------

func fermionMag(u: float, nmax: int): float =
  ## |G^(1,1)| at separation u > 0.  From (V.3),
  ##   (1/4pi) sum_{n>=0} (n+1) e^{-(n+1)u}, and sum_{n>=0}(n+1)x^{n+1} = x/(1-x)^2 at
  ##   x = e^-u, so the closed form is 1/(16 pi sinh^2(u/2)).
  if nmax < 0:
    let s = sinh(0.5*u)
    result = 1.0/(16.0*PI*s*s)
  else:
    var s = 0.0
    for n in countdown(nmax, 0):        # ascending magnitude
      s += (n + 1).float*exp(-(n + 1).float*u)
    result = s/(4.0*PI)

func fermionG*(t: float, nmax = -1): float =
  ## (V.3), (1,1) component:
  ##   G(t) = sign(t)/(16 pi sinh^2(|t|/2)) -> sign(t)/(4 pi t^2) as t -> 0,
  ## the 3D flat-space normalization.  t = 0 is the short-distance singularity.
  sgn(t).float*fermionMag(abs(t), nmax)

func fermionGPeriodic*(t, T: float, nmax = -1): float =
  ## Antiperiodic in t, so the images alternate: G_T(t) = sum_k (-1)^k G(t + kT).
  ## Write t = t_r + kT with t_r in [0,T).  For the images at t_r + jT the argument is
  ## positive and for those at t_r - (j+1)T it is negative, and in the latter the sign(t)
  ## of (V.3) cancels the extra (-1) of the antiperiodic image, leaving
  ##   G_T(t) = (-1)^k sum_{j>=0} (-1)^j [ h(t_r + jT) + h((j+1)T - t_r) ],   h = |G|.
  ## Both terms swap under t_r -> T - t_r, so the shape is symmetric about T/2 -- this is
  ## the cancellation noted in the paper, and it is what makes cosh fits legitimate.
  ##
  ## Term by term this equals the thermal mode sum
  ##   (1/4pi) sum_n (n+1) [e^{-(n+1)t_r} + e^{-(n+1)(T-t_r)}] / (1 + e^{-(n+1)T}),
  ## but the image form reuses the closed form and needs only a couple of terms.
  let
    k = floor(t/T)
    tr = t - k*T
    sg = if (k.int and 1) == 0: 1.0 else: -1.0
  var
    s = 0.0
    j = 0
  while true:
    let d = fermionMag(tr + j.float*T, nmax) + fermionMag((j + 1).float*T - tr, nmax)
    s += (if (j and 1) == 0: d else: -d)
    if d <= sumEps*abs(s) or j >= maxImage: break
    inc j
  sg*s

# --- Gauge current correlator on S^2 x R ------------------------------------

func gaugeG*(t: float, nmax = 200): float =
  ## (V.14): G_g(t) = (1/8pi) sum_{n>=1} sqrt(n(n+1)) (2n+1) e^{-sqrt(n(n+1))|t|}.
  ## The tower Delta_l = sqrt(l(l+1)) is irrational, so there is no closed form;
  ## nmax < 0 sums until the terms underflow.
  let u = abs(t)
  var
    s = 0.0
    n = 1
  while true:
    let
      lam = sqrt(n.float*(n + 1).float)
      term = lam*(2*n + 1).float*exp(-lam*u)
    s += term
    if (if nmax >= 0: n >= nmax else: (n > 2 and term <= sumEps*s) or n >= maxMode): break
    inc n
  s/(8.0*PI)

func gaugeGPeriodic*(t, T: float, nmax = 200): float =
  ## The gauge field is periodic in t, so the images carry no sign:
  ##   G_T(t) = sum_k G_g(t + kT) = sum_{j>=0} [ G_g(t_r + jT) + G_g((j+1)T - t_r) ],
  ## t_r = t mod T.  Symmetric about T/2 and, mode by mode, equal to the truncated
  ## thermal sum with 1/(1 - e^{-lam T}).
  let tr = t - floor(t/T)*T
  var
    s = 0.0
    j = 0
  while true:
    let d = gaugeG(tr + j.float*T, nmax) + gaugeG((j + 1).float*T - tr, nmax)
    s += d
    if d <= sumEps*s or j >= maxImage: break
    inc j
  s

# --- Flat-limit Wilson spectrum ---------------------------------------------

func flatSpectrum*(kap, kapT: float, n1, n2, nt: int): seq[Complex64] =
  ## (IV.8), the Wilson spectrum on the equilateral triangular lattice with radial extension:
  ##   lam = kap (3 - cos k1 - cos k2 - cos(k1+k2)) + kapT (1 - cos kt)
  ##         +- i sqrt( kap^2 |e12 sin k1 + e23 sin k2 + (e12+e23) sin(k1+k2)|^2
  ##                    + kapT^2 sin^2 kt )
  ## with e12 = (1,0), e23 = (-1/2, sqrt3/2).  Free-limit couplings are kap = 1/sqrt3,
  ## kapT = (sqrt3/2)(abar/at).
  ##
  ## Momenta k = 2 pi i / n, i = 0 ..< n.  Both branches are emitted consecutively, so the
  ## result holds 2*n1*n2*nt entries and is closed under conjugation by construction.
  const
    e12 = [1.0, 0.0]
    e23 = [-0.5, 0.5*sqrt(3.0)]
    esum = [e12[0] + e23[0], e12[1] + e23[1]]
  result = newSeqOfCap[Complex64](2*n1*n2*nt)
  for i1 in 0..<n1:
    let k1 = 2.0*PI*i1.float/n1.float
    for i2 in 0..<n2:
      let
        k2 = 2.0*PI*i2.float/n2.float
        s1 = sin(k1)
        s2 = sin(k2)
        s12 = sin(k1 + k2)
        vx = e12[0]*s1 + e23[0]*s2 + esum[0]*s12
        vy = e12[1]*s1 + e23[1]*s2 + esum[1]*s12
        v2 = kap*kap*(vx*vx + vy*vy)
        re0 = kap*(3.0 - cos(k1) - cos(k2) - cos(k1 + k2))
      for it in 0..<nt:
        let
          kt = 2.0*PI*it.float/nt.float
          st = sin(kt)
          re = re0 + kapT*(1.0 - cos(kt))
          im = sqrt(v2 + kapT*kapT*st*st)
        result.add c(re, im)
        result.add c(re, -im)

# --- Two-dimensional checks on S^2 ------------------------------------------

func xiHalf(n: int, w: float): float =
  ## xi_{|m|,n} of (C.10) specialized to |m| = 1/2:
  ##   xi_{|m|,n}(w) = (1-w)^{(|m|-1/2)/2} (1+w)^{-(|m|+1/2)/2} P_{n+|m|+1/2}^{(|m|-1/2,-|m|-1/2)}(w)
  ## At |m| = 1/2 the first exponent is 0 and the second is -1/2, and the Jacobi indices are
  ##   j = n+1,  alpha = 0,  beta = -1,
  ## so j, alpha+j and beta+j are all non-negative integers and Rodrigues' formula (hence the
  ## ordinary recurrence) applies even though beta < 0.
  ##
  ## The beta = -1 family is degenerate, P_{n+1}^{(0,-1)}(w) = ((1+w)/2) P_n^{(0,1)}(w), so
  ##   xi_{1/2,n}(w) = (1/2) sqrt(1+w) P_n^{(0,1)}(w)
  ## with no removable singularity at w = -1.  tanalytic pins the reduction; the literal
  ## (C.10) form is what is evaluated here.
  jacobiP(n + 1, 0.0, -1.0, w)/sqrt(1.0 + w)

func s2FermionProp*(theta: float, nmax: int): float =
  ## (C.55) truncated at n = nmax, as the coefficient of sigma1 (the paper plots the (1,2)
  ## spinor component, and the matrix is sigma1 times this scalar):
  ##   <psi(theta,0) psibar(0,0)> = sigma1 sum_{n>=0} (-1)^n/(pi sqrt2) xi_{1/2,n}(-z),
  ## z = cos(theta).  The limit is the CFT answer (C.54) sigma1/(4 pi sin(theta/2)).
  ##
  ## The series is only conditionally convergent: with the reduction above the sum is
  ##   (sin(theta/2)/2pi) sum_n (-1)^n P_n^{(0,1)}(-z),
  ## whose Jacobi generating function at t = -1 gives 1/(2 sin^2(theta/2)) -- the exact answer,
  ## but on the boundary of the disc of convergence.  Partial sums therefore oscillate with an
  ## O(n^-1/2) envelope; what nmax buys is the ultraviolet, restored down to theta ~ pi/nmax.
  let z = cos(theta)
  var s = 0.0
  for n in 0..nmax:
    let x = xiHalf(n, -z)
    s += (if (n and 1) == 0: x else: -x)
  s/(PI*sqrt(2.0))

func f0Deriv(n: int, z: float): float =
  ## d/dz f_{0,n}(z) with f_{0,n}(z) = (1-z) P_n^{(1,-1)}(z), (C.43), using
  ##   d/dz P_j^(a,b) = ((j+a+b+1)/2) P_{j-1}^(a+1,b+1).
  -jacobiP(n, 1.0, -1.0, z) + (1.0 - z)*0.5*(n + 1).float*jacobiP(n - 1, 2.0, 0.0, z)

func s2CurrentCorr*(z: float, nmax: int): float =
  ## (C.57) truncated at n = nmax:
  ##   (1/g^2)<J^t(theta,0) J^t(0,0)> = -(1/4pi) sum_{n>=1} ((2n+1)/(n+1)) f'_{0,n}(z)
  ##                                  = (1/2pi)[delta(z-1) - 1/2].
  ## The Legendre equation collapses f'_{0,n} = -(n+1) P_n(z), so the partial sum is
  ## (1/4pi) sum_{n=1..nmax} (2n+1) P_n(z): the completeness sum for the covariant delta on
  ## S^2 with the constant mode (the absent gauge zero mode) removed.  Evaluated in the
  ## literal (C.57) form; tanalytic checks the collapse.
  var s = 0.0
  for n in 1..nmax:
    s += (2*n + 1).float/(n + 1).float*f0Deriv(n, z)
  -s/(4.0*PI)

# --- Effective dimension ----------------------------------------------------

func effDim*(g: openArray[float], at, T: float): seq[float] =
  ## (V.4)-(V.5): f(t) = arccosh(G(t)/G(T/2)), Delta_eff(t) = -[f(t+at) - f(t)]/at.
  ## `g[i]` is the correlator at t = i*at and `result[i]` is Delta_eff(i*at), so the result
  ## is one shorter than the input.  G(T/2) is taken at index round(T/(2 at)).
  ##
  ## For a single state with (anti)periodic images G(t) = A cosh(Delta (T/2 - t)) the ratio is
  ## exactly cosh(Delta (T/2 - t)), f is linear and Delta_eff = Delta with no discretization
  ## error; the residual is pure excited-state contamination.  Past T/2 the sign flips.
  let gh = g[int(round(0.5*T/at))]
  result = newSeq[float](g.len - 1)
  var f0 = arccosh(g[0]/gh)
  for i in 0..<g.len - 1:
    let f1 = arccosh(g[i + 1]/gh)
    result[i] = -(f1 - f0)/at
    f0 = f1
