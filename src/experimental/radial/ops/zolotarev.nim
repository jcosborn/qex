## Zolotarev optimal rational approximation to 1/sqrt(x).
##
##   R(x) = cst + sum_j res_j/(x + pole_j) ~ 1/sqrt(x),   x in [smin^2, smax^2]
##
## the minimax (equioscillating) rational of type (m,m), m = (order-1)/2, order odd.
## Chiu, Hsieh, Huang, Huang, PRD 66 (2002) 114502 [hep-lat/0206007];
## van den Eshof, Frommer, Lippert, Schilling, van der Vorst, CPC 146 (2002) 203.
##
## Everything is stored in sigma^2 units: the operator is H = X^dag X with
## spec(H) in [smin^2, smax^2], so the poles that go into the multishift CG are
## shifts of H, not of X.
##
## Needs the complete elliptic integral K and the Jacobi sn/cn; QEX has neither.

import std/math

const
  agmIt = 60          ## AGM converges quadratically; a hard stop for k -> 1
  agmTol = 1e-17

func agm(x, y: float): float =
  ## a_{n+1} = (a_n + b_n)/2,  b_{n+1} = sqrt(a_n b_n)
  var
    a = x
    b = y
  for _ in 1..agmIt:
    if abs(a-b) <= agmTol*a: break
    let t = 0.5*(a+b)
    b = sqrt(a*b)
    a = t
  a

func ellipticK*(k: float): float =
  ## K(k) = int_0^{pi/2} dt/sqrt(1 - k^2 sin^2 t) = pi/(2 agm(1, k')),  k' = sqrt(1-k^2).
  0.5*PI/agm(1.0, sqrt(1.0 - k*k))

func sncn(u, k: float): tuple[sn, cn: float] =
  ## Jacobi sn, cn by the descending Landen / AGM recursion (A&S 16.4):
  ##   a_0 = 1, b_0 = k', c_0 = k;   a_{n+1},b_{n+1} as in agm, c_{n+1} = (a_n - b_n)/2
  ##   phi_N = 2^N a_N u,   phi_{n-1} = (phi_n + arcsin((c_n/a_n) sin phi_n))/2
  ##   sn = sin phi_0,  cn = cos phi_0
  ## cn comes from cos(phi_0), never from sqrt(1-sn^2): the poles need cn^2 where
  ## sn -> 1 and the subtraction there throws away every digit.
  ## k = 0 gives (sin u, cos u); the AGM stalls at k = 1, whose limit is (tanh u, sech u).
  if abs(k) >= 1.0: return (tanh(u), 1.0/cosh(u))
  var
    a, c: array[agmIt+1, float]
    b = sqrt(1.0 - k*k)
    n = 0
  a[0] = 1.0
  c[0] = abs(k)
  while n < agmIt and c[n] > agmTol*a[n]:
    let bn = sqrt(a[n]*b)
    c[n+1] = 0.5*(a[n] - b)
    a[n+1] = 0.5*(a[n] + b)
    b = bn
    inc n
  var f = a[n]*u*float(1'i64 shl n)
  for i in countdown(n, 1):
    f = 0.5*(f + arcsin(c[i]/a[i]*sin(f)))
  (sin(f), cos(f))

func jacobiSn*(u, k: float): float = sncn(u, k)[0]

type Rat* = object
  order*, npole*: int
  smin*, smax*: float              ## sigma bounds (NOT squared)
  cst*: float                      ## constant term
  pole*, res*, zero*: seq[float]   ## sigma^2 units; 0 < pole_0 < zero_0 < pole_1 < ..., res > 0
  maxRelErr*, maxAbsErr*: float    ## max |1 - sqrt(x) R(x)| and max |R(x) - 1/sqrt(x)|
  hash*: uint64

func ratValue*(r: Rat, x: float): float =
  ## R(x) = cst + sum_j res_j/(x + pole_j)
  result = r.cst
  for j in 0..<r.npole: result += r.res[j]/(x + r.pole[j])

const
  fnvBasis = 0xcbf29ce484222325'u64
  fnvPrime = 0x100000001b3'u64

func fnv(h, v: uint64): uint64 =
  ## FNV-1a over the 8 bytes of v, little end first.
  result = h
  var x = v
  for _ in 0..7:
    result = (result xor (x and 0xff'u64))*fnvPrime
    x = x shr 8

func ratHash(r: Rat): uint64 =
  ## Fingerprint of the frozen rational: (order, smin, smax, cst, poles, residues).
  result = fnv(fnvBasis, uint64(r.order))
  for v in [r.smin, r.smax, r.cst]: result = fnv(result, cast[uint64](v))
  for v in r.pole: result = fnv(result, cast[uint64](v))
  for v in r.res: result = fnv(result, cast[uint64](v))

proc newRat*(smin, smax: float, order: int, nsample = 20001): Rat =
  ## Zolotarev rational for 1/sqrt(x) on [smin^2, smax^2].  `order` odd and >= 3;
  ## `nsample` >= 2 log-spaced probes of the error (a bad count would silently report zero).
  if smin <= 0.0 or smax <= smin:
    raise newException(ValueError, "newRat needs 0 < smin < smax, got " & $smin & " " & $smax)
  if order < 3 or order mod 2 == 0:
    raise newException(ValueError, "newRat needs odd order >= 3, got " & $order)
  if nsample < 2:
    raise newException(ValueError, "newRat needs nsample >= 2, got " & $nsample)
  let
    m = (order-1) div 2
    k = smin/smax
    kp = sqrt(1.0 - k*k)                    ## k' = sqrt(1-k^2)
    du = ellipticK(kp)/float(order)         ## K'/n
  # c_i = sn^2(i K'/n; k')/cn^2(i K'/n; k'), i = 1..2m, on the reference window [k^2, 1]:
  # poles at the odd i, zeros at the even i, both scaled by k^2.
  var
    p = newSeq[float](m)
    z = newSeq[float](m)
  for j in 0..<m:
    let
      (sp, cp) = sncn(float(2*j+1)*du, kp)
      (sz, cz) = sncn(float(2*j+2)*du, kp)
    p[j] = k*k*(sp*sp)/(cp*cp)
    z[j] = k*k*(sz*sz)/(cz*cz)
  # P(x) = prod_i (x+z_i)/(x+p_i);  e(x) = 1 - sqrt(x) d0 P(x).
  # d0 from endpoint equioscillation e(k^2) = -e(1):  2 = d0 (k P(k^2) + P(1)).
  var
    pk = 1.0
    p1 = 1.0
  for i in 0..<m:
    pk *= (k*k + z[i])/(k*k + p[i])
    p1 *= (1.0 + z[i])/(1.0 + p[i])
  let d0 = 2.0/(k*pk + p1)
  # exact partial fractions of the monic ratio: N/D = 1 + sum_j N(-p_j)/(D'(-p_j) (x+p_j)),
  #   a_j = d0 prod_i (z_i - p_j) / prod_{l/=j} (p_l - p_j).
  # Accumulated in log space with the sign counted separately: at order 31 the two
  # products individually span tens of decades while their ratio is O(1).  Interlacing
  # makes both signs (-1)^j, so a_j > 0; the test checks that rather than assuming it.
  var a = newSeq[float](m)
  for j in 0..<m:
    var
      lg = 0.0
      neg = 0
    for i in 0..<m:
      let d = z[i] - p[j]
      lg += ln(abs(d))
      if d < 0.0: inc neg
    for l in 0..<m:
      if l != j:
        let d = p[l] - p[j]
        lg -= ln(abs(d))
        if d < 0.0: inc neg
    a[j] = (if (neg and 1) == 0: d0 else: -d0)*exp(lg)
  # [k^2,1] -> [smin^2,smax^2]:  R(x) = Rs(x/smax^2)/smax, so pole,zero *= smax^2,
  # res *= smax, cst /= smax.  Getting this wrong is the classic Zolotarev bug.
  let s2 = smax*smax
  result.order = order
  result.npole = m
  result.smin = smin
  result.smax = smax
  result.cst = d0/smax
  result.pole = newSeq[float](m)
  result.zero = newSeq[float](m)
  result.res = newSeq[float](m)
  for j in 0..<m:
    result.pole[j] = p[j]*s2
    result.zero[j] = z[j]*s2
    result.res[j] = a[j]*smax
  # measured, not assumed: e(x) = 1 - sqrt(x) R(x) at nsample log-spaced x
  let
    lo = 2.0*ln(smin)
    dl = 2.0*(ln(smax) - ln(smin))/float(nsample-1)
  for i in 0..<nsample:
    let
      x = exp(lo + dl*float(i))
      sx = sqrt(x)
      e = abs(1.0 - sx*ratValue(result, x))
    if e > result.maxRelErr: result.maxRelErr = e
    if e/sx > result.maxAbsErr: result.maxAbsErr = e/sx
  result.hash = ratHash(result)
