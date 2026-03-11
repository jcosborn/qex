import base
#import basicOps
import complexNumbers
#import complexType
import matrixConcept
import types
import matinv
export matinv
import matexp
export matexp
import projUderiv
export projUderiv
getOptimPragmas()

proc determinantN*(a: auto): auto =
  mixin simdMinReduce
  const nc = a.nrows
  var c {.noInit.}: type(a)
  var row: array[nc,int]
  var nswaps = 0
  var r: type(a[0,0])
  r := 1
  template C(i,j): untyped = c[row[i],j]

  for i in 0..<nc:
    for j in 0..<nc:
      c[i,j] = a[i,j]
    row[i] = i

  for j in 0..<nc:
    if j>0:
      for i in j..<nc:
        var t2 = C(i,j)
        for k in 0..<j:
          t2 -= C(i,k) * C(k,j)
        C(i,j) := t2

    var rmax = C(j,j).norm2
    var kmax = j
    for k in (j+1)..<nc:
      var rn = C(k,j).norm2
      if rn.simdMinReduce > rmax.simdMinReduce:
        rmax = rn
        kmax = k
    #if rmax.simdMin == 0: # some matrix is singular
    #  r := 0  # FIXME need to only adjust rmax an continue
    #  return r
    if kmax != j:
      swap(row[j], row[kmax])
      inc nswaps

    r *= C(j,j)

    let ri = 1.0/rmax
    var Cjji = ri * C(j,j).adj
    for i in (j+1)..<nc:
      var t2 = C(j,i)
      for k in 0..<j:
        t2 -= C(j,k) * C(k,i)
      C(j,i) := t2 * Cjji

  if (nswaps and 1) != 0:
    r := -r

  r

proc determinant*(x: auto): auto {.alwaysInline.} =
  assert(x.nrows == x.ncols)
  when x.nrows==1:
    result = x[0,0]
  elif x.nrows==2:
    result = x[0,0]*x[1,1] - x[0,1]*x[1,0]
  elif x.nrows==3:
    result = (x[0,0]*x[1,1]-x[0,1]*x[1,0])*x[2,2] +
             (x[0,2]*x[1,0]-x[0,0]*x[1,2])*x[2,1] +
             (x[0,1]*x[1,2]-x[0,2]*x[1,1])*x[2,0]
  else:
    result = determinantN(x)

proc eigs3(e0,e1,e2: var auto; tr,p2,det: auto) {.alwaysInline.} =
  mixin sin,cos,acos
  let tr3 = (1.0/3.0)*tr
  let p23 = (1.0/3.0)*p2
  let tr32 = tr3*tr3
  let q = abs(0.5*(p23-tr32))
  let r = 0.25*tr3*(5*tr32-p2) - 0.5*det
  let sq = sqrt(q)
  let sq3 = q*sq
  #let rsq3 = r/sq3
  #var minv,maxv {.noinit.}:type(rsq3)
  #minv := -1.0
  #maxv := 1.0
  #let rsq3r = min(maxv, max(minv,rsq3))
  let isq3 = 1.0/sq3
  var minv,maxv {.noinit.}: type(isq3)
  maxv := 3e38
  minv := -3e38
  let isq3c = min(maxv, max(minv,isq3))
  let rsq3c = r * isq3c
  maxv := 1
  minv := -1
  let rsq3 = min(maxv, max(minv,rsq3c))
  let t = (1.0/3.0)*acos(rsq3)
  let st = sin(t)
  let ct = cos(t)
  let sqc = sq*ct
  let sqs = 1.73205080756887729352*sq*st  # sqrt(3)
  let ll = tr3 + sqc
  e0 = tr3 - 2*sqc
  e1 = ll + sqs
  e2 = ll - sqs
  # ~27 flops

template rsqrtPHM2(r:typed; x:typed) =
  let x00 = x[0,0].re
  let x11 = x[1,1].re
  let x01r = 0.5*(x[0,1].re+x[1,0].re)
  let x01i = 0.5*(x[0,1].im-x[1,0].im)
  let det = abs(x00*x11 - x01r*x01r - x01i*x01i)
  let tr = x00 + x11
  let sdet = sqrt(det)
  let trsdet = tr + sdet
  let c1 = 1/(sdet*sqrt(trsdet+sdet))
  let c0 = trsdet*c1
  r := c0 - c1*x

proc rsqrtPHM3f(c0,c1,c2:var auto; tr,p2,det:auto) {.alwaysInline.} =
  #[
  mixin sin,cos,acos
  let tr3 = (1.0/3.0)*tr
  let p23 = (1.0/3.0)*p2
  let tr32 = tr3*tr3
  let q = abs(0.5*(p23-tr32))
  let r = 0.25*tr3*(5*tr32-p2) - 0.5*det
  let sq = sqrt(q)
  let sq3 = q*sq
  #let rsq3 = r/sq3
  #var minv,maxv {.noinit.}:type(rsq3)
  #minv := -1.0
  #maxv := 1.0
  #let rsq3r = min(maxv, max(minv,rsq3))
  let isq3 = 1.0/sq3
  var minv,maxv {.noinit.}: type(isq3)
  maxv := 3e38
  minv := -3e38
  let isq3c = min(maxv, max(minv,isq3))
  let rsq3c = r * isq3c
  maxv := 1
  minv := -1
  let rsq3 = min(maxv, max(minv,rsq3c))
  let t = (1.0/3.0)*acos(rsq3)
  let st = sin(t)
  let ct = cos(t)
  let sqc = sq*ct
  let sqs = 1.73205080756887729352*sq*st  # sqrt(3)
  let l0 = tr3 - 2*sqc
  let ll = tr3 + sqc
  let l1 = ll + sqs
  let l2 = ll - sqs
  ]#
  var l0,l1,l2 {.noInit.}: type(tr)
  eigs3(l0,l1,l2, tr,p2,det)
  let sl0 = sqrt(abs(l0))
  let sl1 = sqrt(abs(l1))
  let sl2 = sqrt(abs(l2))
  let u = sl0 + sl1 + sl2
  let w = sl0 * sl1 * sl2
  let d = w*(sl0+sl1)*(sl0+sl2)*(sl1+sl2)
  let di = 1/d
  c0 = (w*u*u+l0*sl0*(l1+l2)+l1*sl1*(l0+l2)+l2*sl2*(l0+l1))*di
  c1 = -(tr*u+w)*di
  c2 = u*di
  # flops: eigs3(27) + 33 = 60

template rsqrtPHM3(r:typed; x:typed) =
  let tr = trace(x).re
  let x2 = x*x
  let p2 = trace(x2).re
  let det = determinant(x).re
  var c0,c1,c2:type(tr)
  rsqrtPHM3f(c0, c1, c2, tr, p2, det)
  r := c0 + c1*x + c2*x2
  # 2*tr(2) + mm(66) + det(64) + f(60) + 56 = 250

template rsqrtPHMN(r:typed; x:typed) =
  mixin simdMax
  let xi = 1/x
  let xi2 = xi.norm2
  let xit = trace(xi).re
  let ds = xit/xi2
  #var ds = x.norm2.simdMax
  #ds = 0.5*sqrt(ds)
  #echo "ds: ", ds

  var e = (0.5*ds)*xi - 0.5
  var s = 1 + e
  #echo "e: ", e
  #echo "s: ", s

  let estop = epsilon(ds.simdMax)^2
  let maxit = 20
  var nit = 0
  while true:
    inc nit
    #let t = (e/s) * e
    #let t = e * (s \ e)
    let si = 1/s
    let t = e * (si * e)
    e := -0.5 * t
    s += e
    let enorm = e.norm2.simdMax
    #echo nit, " enorm: ", enorm
    if nit>=maxit or enorm<estop: break
  let sds = 1/sqrt(ds)
  r := sds*s

#[
# Bini (https://arxiv.org/pdf/1703.02456.pdf)
template rsqrtPHMN2(r:typed; x:typed) =
  let xn = x.norm2
  let ds = sqrt(xn)
  let dsi = 3/ds
  #echo "ds: ", ds

  var a = dsi * x
  var b {.noInit.} :type(r)
  b := 1
  #var b = (1.4/xn)*x

  let estop = epsilon(ds.simdMax)^2
  let maxit = 20
  var nit = 0
  while true:
    let e = 1 - b*a*b
    let enorm = e.norm2.simdMax
    echo nit, " enorm: ", enorm
    if nit>=maxit or enorm<estop: break
    inc nit
    let t = b*e
    #let t2 = t.norm2
    #let c = 0.5/sqrt(t2)
    let c = 0.5
    b += c*t
  let sds = sqrt(dsi)
  r := sds*b
  #r := b

# -0.5[B+((1-BBA)^3-1)/(BA)] = -0.5B[1-3+3BBA-BBBBAA] = 0.5B[2-3BBA+BBBBAA]
# = 0.5B[1-BBA][2-BBA]
template rsqrtPHMN3(r:typed; x:typed) =
  let xn = x.norm2
  let ds = sqrt(xn)
  let dsi = 1/ds
  #echo "ds: ", ds

  var a = dsi * x
  var b {.noInit.} :type(r)
  b := 1
  #var b = (1.4/xn)*x

  let estop = epsilon(ds.simdMax)^2
  let maxit = 20
  var nit = 0
  while true:
    let e = 1 - b*a*b
    let enorm = e.norm2.simdMax
    echo nit, " enorm: ", enorm
    if nit>=maxit or enorm<estop: break
    inc nit
    let t = b*e*(1+e)
    #let t2 = t.norm2
    #let c = 0.5/sqrt(t2)
    let c = 0.5
    b += c*t
  let sds = sqrt(dsi)
  r := sds*b
  #r := b
]#

template rsqrtPHM(r:typed; x:typed) =
  mixin rsqrt, nrows
  assert(r.nrows == x.nrows)
  assert(r.ncols == x.ncols)
  assert(r.nrows == r.ncols)
  when r.nrows==1:
    let t = rsqrt(x[0,0].re)
    r := t
  elif r.nrows==2:
    rsqrtPHM2(r, x)
  elif r.nrows==3:
    rsqrtPHM3(r, x)
  else:
    rsqrtPHMN(r, x)
    #rsqrtPHMN2(r, x)
    #rsqrtPHMN3(r, x)
proc rsqrtPH*(r: var Mat1; x: Mat2) = rsqrtPHM(r, x)
template rsqrtPH*[T:Mat1](x: T): T =
  var r {.noInit.}: T
  rsqrtPH(r, x)
  r

proc projectUrsqrt*(r: var Mat1; x: Mat2, eps = 1e-20) {.alwaysInline.} =
  #let t = x.adj * x   # issues with gcc
  let xa = x.adj
  var t = xa * x
  t += eps
  rsqrtPHM(r, t)

# x (x'x)^{-1/2}
proc projectU*(r: var Mat1; x: Mat2, eps = 1e-20) {.inline.} =
  var t2{.noInit.}: evalType(x)
  projectUrsqrt(t2, x)
  mul(r, x, t2)
  #echo "t: ", t.norm2, "  t2: ", t2.norm2, "  r: ", r.norm2
template projectUflops*(nc: int): int =
  nc*nc*(2*(6*nc+2*(nc-1))) + 250  # only for nc=3

template projectU*(r: var Mat1, eps = 1e-20) =
  var t{.noInit.}: evalType(r)
  t := r
  r.projectU t, eps

# (d/dX') Tr(U'C+C'U) / 2 = (d/dX') Tr(X'CZ+C'XZ) / 2
# = CZ - (1/2) < Z (X'C + C'X) Z (dY/dX') >
# (dY/dX') Y + Y (dY/dX') = 2X
# Z(dY/dX') Y + Y Z(dY/dX') = 2ZX
# S Y + Y S = U' C Z = Z X' C Z
# cz-xz^3(x'c+c'x)/2
# CH: 4528 flops
proc projectUderiv*(r: var Mat1, u: Mat2, x: Mat3, chain: Mat4, eps = 1e-20) =
  # U = X (X'X)^{-1/2} = (XX')^{-1/2} X
  # Y = sqrt(X'X)
  # Z = (X'X)^{-1/2}
  # F = C Z - z (Cd U + Ud C) z (dY/dX)
  var y, z, t1, t2: Mat1
  #y := x.adj * u
  #inverse(z, y)
  projectUrsqrt(z, x)
  inverse(y, z)
  #echo "inverse: ", z
  #QLA_M_eq_M_times_M(d, c, &z);
  r := chain * z
  #QLA_M_eq_Ma_times_M(&t1, p, d);
  t1 := u.adj * r
  #sylsolve_site(NCARG &t2, &y, &y, &t1);
  sylsolve(t2, y, t1)
  #QLA_M_eq_M(&t1, &t2);
  #QLA_M_peq_Ma(&t1, &t2);
  #QLA_M_eq_M_times_M(&t2, m, &t1);
  #QLA_M_meq_M(d, &t2);
  t1 := t2 + t2.adj
  r -= x * t1

proc projectUderiv*(r: var Mat1, x: Mat2, c: Mat3) =
  var u {.noInit.}: type(r)
  projectU(u, x)
  #echo u, x, c
  projectUderiv(r, u, x, c)

proc projectUVJP*(r: var Mat1, u: Mat2, x: Mat3, chain: Mat4, eps = 1e-20) =
  ## Alias for projectUderiv with naming based on the kernel.
  projectUderiv(r, u, x, chain, eps)

proc projectUVJP*(r: var Mat1, x: Mat2, c: Mat3) =
  ## Convenience alias for projectUderiv without providing u.
  projectUderiv(r, x, c)

proc projectUHVPu*(r: var Mat1, u: Mat2, x: Mat3, c: Mat4, dx: Mat3, dc: Mat4, eps = 1e-20) =
  ## Analytical directional second derivative of the unitary projection.
  ##
  ## Given:
  ##   u     = projectU(x)
  ##   r₁(x) = projectUderiv(u, x, chain)
  ##
  ## this returns
  ##   r = d/dε [ r₁(x + ε dx, chain + ε dc) ] |_{ε=0}
  ##
  ## by differentiating the algebra used in `projectUderiv`:
  ##   Z = (X†X)^(-1/2),  Y = Z^(-1)
  ##   R0 = C Z
  ##   T2 solves  Y T2 + T2 Y = U† R0
  ##   S  = T2 + T2†
  ##   r1 = R0 - X S
  ##
  ## Variations:
  ##   dA  = X† dX + dX† X
  ##   dY  solves Y dY + dY Y = dA
  ##   dZ  = - Z dY Z
  ##   dR0 = dC Z + C dZ
  ##   dU  = dX Z + X dZ   ,   dU† = Z dX† + dZ X†
  ##   dB  = dU† R0 + U† dR0          (B = U† R0)
  ##   dT2 solves Y dT2 + dT2 Y = dB - dY T2 - T2 dY
  ##   dS  = dT2 + dT2†
  ##   r   = dR0 - dX S - X dS
  ##
  ## This keeps all Sylvester solves local and avoids another outer finite
  ## difference.
  discard eps  # parameter kept for API compatibility
  var z, y: Mat1
  projectUrsqrt(z, x)   # Z
  inverse(y, z)         # Y = Z^{-1} = sqrt(X†X)

  # Base first-derivative pieces (same as projectUderiv).
  var r0, t1, s: Mat1
  mul(r0, c, z)       # R0 = C Z
  mul(t1, u.adj, r0)      # B  = U† R0
  # Solve directly for S: Y S + S Y = B + B†
  var bSym: Mat1
  bSym := t1
  bSym += t1.adj
  sylsolve(s, y, bSym)

  # dA = X† dX + dX† X
  var da, tmp: Mat1
  mul(da, x.adj, dx)
  mul(tmp, dx.adj, x)
  da += tmp
  # Keep dA Hermitian (numerical symmetrization).
  da = 0.5*(da + da.adj)

  # dY from sqrt(A): Y dY + dY Y = dA
  var dy: Mat1
  sylsolve(dy, y, da)
  dy = 0.5*(dy + dy.adj)

  # dZ via inverse derivative: dZ = - Z dY Z (exact for Z = Y^{-1}).
  var dz: Mat1
  var zdy: Mat1
  mul(zdy, z, dy)
  mul(dz, zdy, z)
  dz *= -1
  dz = 0.5*(dz + dz.adj)

  # dR0 = dC Z + C dZ
  var dr0: Mat1
  mul(dr0, dc, z)
  mul(tmp, c, dz)
  dr0 += tmp

  # dU = dX Z + X dZ
  var du: Mat1
  mul(du, dx, z)
  mul(tmp, x, dz)
  du += tmp

  # dU† = Z dX† + dZ X†
  var duAdj: Mat1
  mul(duAdj, z, dx.adj)
  mul(tmp, dz, x.adj)
  duAdj += tmp               # full dU† including dZ contribution

  # dB pieces
  # Full variation: dB = dU† R0 + U† dR0
  # However, the dR0 contribution through U† dR0 is separate from the direct dR0 term.
  # We include both here as the Sylvester equation for dS needs the full dB.
  var db: Mat1
  mul(db, duAdj, r0)           # dU† R0
  mul(tmp, u.adj, dr0)         # U† dR0
  db += tmp                    # full dB

  # rhs for dT2: dB - dY T2 - T2 dY
  # Solve directly for dS: Y dS + dS Y = d(B+B†) - dY S - S dY
  var rhs, dyS, Sdy: Mat1
  # Build full symmetric RHS for dS:
  # RHS = (dB + dB†) - (dY S + S dY)
  # Note: dY and S are both Hermitian, so dY S + S dY is Hermitian
  var dBsym: Mat1
  # dB = dU† R0 + U† dR0, already in db
  dBsym := db + db.adj
  # dY S + S dY (already Hermitian since dY and S are Hermitian)
  mul(dyS, dy, s)
  mul(Sdy, s, dy)
  rhs := dBsym
  rhs -= dyS
  rhs -= Sdy

  var ds: Mat1
  sylsolve(ds, y, rhs)
  ds = 0.5*(ds + ds.adj)  # Ensure ds is Hermitian

  # Final directional second derivative: dR0 - dX S - X dS
  tmp := 0
  mul(tmp, dx, s)
  r := dr0
  r -= tmp
  tmp := 0
  mul(tmp, x, ds)
  r -= tmp

proc projectUHVPu*(r: var Mat1, u: Mat2, x: Mat3, c: Mat4, dx: Mat3, eps = 1e-20) =
  ## Compatibility wrapper when no upstream variation dc is provided.
  var dc {.noInit.}: type(c)
  dc := 0
  projectUHVPu(r, u, x, c, dx, dc, eps)

proc projectUHVP*(r: var Mat1, x: Mat2, c: Mat3, dx: Mat2, eps = 1e-20) =
  ## Convenience wrapper for projectUHVP when only x is given.
  var u {.noInit.}: type(r)
  var dc {.noInit.}: type(c)
  projectU(u, x)
  dc := 0
  projectUHVPu(r, u, x, c, dx, dc, eps)

proc projectUHVP*(r: var Mat1, x: Mat2, c: Mat3, dx: Mat2, dc: Mat3, eps = 1e-20) =
  ## Convenience wrapper including an upstream variation dc.
  var u {.noInit.}: type(r)
  projectU(u, x)
  projectUHVPu(r, u, x, c, dx, dc, eps)

proc projectUJVP*(du: var Mat1, u: Mat2, x: Mat3, dx: Mat3, eps = 1e-20) =
  ## Forward tangent of unitary projection.
  ##
  ## Computes: du = d/dε[ projectU(x + ε dx) ]|_{ε=0}
  ##
  ## Given U = X Z where Z = (X†X)^{-1/2}, Y = Z^{-1} = sqrt(X†X):
  ##   dU = dX Z + X dZ
  ## where dZ = -Z dY Z and Y dY + dY Y = dA, dA = X† dX + dX† X
  discard eps
  var z, y: Mat1
  projectUrsqrt(z, x)   # Z = (X†X)^{-1/2}
  inverse(y, z)         # Y = Z^{-1}

  # dA = X† dX + dX† X (Hermitian)
  var da, tmp: Mat1
  mul(tmp, x.adj, dx)
  da := tmp + tmp.adj

  # Solve Y dY + dY Y = dA
  var dy: Mat1
  sylsolve(dy, y, da)

  # dZ = -Z dY Z
  var dz: Mat1
  tmp := dy * z
  dz := z * tmp
  dz := -dz

  # dU = dX Z + X dZ
  mul(du, dx, z)
  tmp := x * dz
  du += tmp

proc projectUJVP*(du: var Mat1, x: Mat2, dx: Mat2, eps = 1e-20) =
  ## Convenience wrapper when u is not provided.
  var u {.noInit.}: type(du)
  projectU(u, x)
  projectUJVP(du, u, x, dx, eps)

proc projectUVJPChain*(chainbar: var Mat1, u: Mat2, x: Mat3, rbar: Mat4, eps = 1e-20) =
  ## Adjoint of projectUderiv with respect to the chain input.
  ##
  ## Given: r = projectUderiv(u, x, chain) which computes
  ##   Z = (X†X)^{-1/2}, Y = Z^{-1}
  ##   R0 = chain * Z
  ##   T1 = U† * R0
  ##   sylsolve(T2, Y, T1)  // Y*T2 + T2*Y = T1
  ##   S = T2 + T2†
  ##   r = R0 - X*S
  ##
  ## This computes chainbar such that ⟨rbar, r⟩ = ⟨chainbar, chain⟩
  ## for any chain (the map chain → r is linear).
  ##
  ## Derivation:
  ##   r = chain*Z - X*(T2 + T2†) where Y*T2 + T2*Y = U†*chain*Z
  ##
  ##   Term 1: chain*Z contributes Z*rbar to chainbar (Z is Hermitian)
  ##
  ##   Term 2: -X*S where S = T2 + T2†
  ##   ⟨-X*S, rbar⟩ = ⟨S, -X†*rbar⟩
  ##   Let Q = -X†*rbar, then ⟨S, Q⟩ = ⟨T2+T2†, Q⟩ = ⟨T2, Q+Q†⟩
  ##   Adjoint of Sylvester Y*T2 + T2*Y = T1 gives: Y*T2bar + T2bar*Y = Q+Q†
  ##   Then ⟨T1, T2bar⟩ = ⟨U†*chain*Z, T2bar⟩ = ⟨chain, U*T2bar*Z⟩
  ##
  ##   Total: chainbar = Z*rbar + U*T2bar*Z
  discard eps
  var z, y: Mat1
  projectUrsqrt(z, x)   # Z = (X†X)^{-1/2}
  inverse(y, z)         # Y = Z^{-1}

  # Q = -X†*rbar, Q_sym = Q + Q†
  var q, qsym, tmp: Mat1
  mul(q, x.adj, rbar)
  q := -q
  qsym := q + q.adj

  # Solve Y*T2bar + T2bar*Y = Q_sym
  var t2bar: Mat1
  sylsolve(t2bar, y, qsym)

  # chainbar = rbar*Z + U*T2bar*Z
  mul(chainbar, rbar, z)       # rbar*Z (Z is Hermitian, so Z† = Z)
  mul(tmp, u, t2bar)
  tmp := tmp * z
  chainbar += tmp

proc projectUVJPChain*(chainbar: var Mat1, x: Mat2, rbar: Mat3, eps = 1e-20) =
  ## Convenience wrapper when u is not provided.
  var u {.noInit.}: type(chainbar)
  projectU(u, x)
  projectUVJPChain(chainbar, u, x, rbar, eps)

proc projectUHVPVJP_dx*(dxbar: var Mat1, u: Mat2, x: Mat3, c: Mat4, rbar: auto, eps = 1e-20) =
  ## Adjoint of projectUHVP with respect to dx.
  ##
  ## Given: r = projectUHVPu(u, x, c, dx) which computes (with dc=0):
  ##   Z = (X†X)^{-1/2}, Y = Z^{-1}
  ##   R0 = C Z
  ##   B = U† R0
  ##   S solves Y S + S Y = B + B†
  ##   dA = X† dx + dx† X
  ##   dY solves Y dY + dY Y = dA
  ##   dZ = -Z dY Z
  ##   dR0 = C dZ  (since dc=0)
  ##   dU = dx Z + X dZ
  ##   dU† = Z dx† + dZ X†
  ##   dB = dU† R0 + U† dR0
  ##   dS solves Y dS + dS Y = (dB + dB†) - dY S - S dY
  ##   r = dR0 - dx S - X dS
  ##
  ## This computes dxbar such that ⟨rbar, r⟩ = ⟨dxbar, dx⟩.
  ##
  ## dx appears in:
  ##   1. Direct: -dx S  →  dxbar += -rbar S†  (S is Hermitian, so S† = S)
  ##   2. dU = dx Z  →  contribution through dU (NOT dU†!)
  ##   3. dU† = Z dx†  →  contribution through dU† to dB
  ##   4. dA = X† dx + dx† X  →  flows through dY, dZ
  ##
  ## Backward pass traces through the chain of dependencies.
  discard eps
  var z, y: Mat1
  projectUrsqrt(z, x)   # Z = (X†X)^{-1/2}
  inverse(y, z)         # Y = Z^{-1}

  # Forward quantities needed for adjoint
  var r0, b, s, tmp: Mat1
  mul(r0, c, z)         # R0 = C Z
  mul(b, u.adj, r0)     # B = U† R0
  var bSym: Mat1
  bSym := b + b.adj
  sylsolve(s, y, bSym)  # Y S + S Y = B + B†

  # === Backward pass ===

  # Step 1: Direct term -dx S
  # ⟨rbar, -dx S⟩ = ⟨-rbar S, dx⟩ (S is Hermitian)
  mul(dxbar, rbar, s)
  dxbar := -dxbar

  # Step 2: Trace through dS
  # r -= X dS, so rbar gives contribution to dSbar
  # ⟨rbar, -X dS⟩ = ⟨-X† rbar, dS⟩
  var dSbar: Mat1
  mul(dSbar, x.adj, rbar)
  dSbar := -dSbar

  # dS solves Y dS + dS Y = RHS where RHS = (dB+dB†) - dY S - S dY
  # RHS is Hermitian (sum of Hermitian terms), so dS is already Hermitian.
  # The 0.5 symmetrization in forward is just numerical cleanup.
  # Sylvester adjoint: Y dRHSbar + dRHSbar Y = dSbar
  var dRHSbar: Mat1
  sylsolve(dRHSbar, y, dSbar)

  # RHS = (dB + dB†) - dY S - S dY
  # From (dB + dB†): ⟨dRHSbar, dB + dB†⟩ = ⟨dRHSbar + dRHSbar†, dB⟩
  var dBbar: Mat1
  dBbar := dRHSbar + dRHSbar.adj

  # From -dY S - S dY: contribution to dYbar
  # ⟨dRHSbar, -dY S - S dY⟩ = -⟨dRHSbar S + S dRHSbar, dY⟩
  var dYbar: Mat1
  var rhsS, Srhs: Mat1
  mul(rhsS, dRHSbar, s)
  mul(Srhs, s, dRHSbar)
  dYbar := -(rhsS + Srhs)

  # Step 3: Trace through dR0
  # r += dR0 (positive sign), so ⟨rbar, dR0⟩
  # dR0 = C dZ (since dc=0)
  # ⟨rbar, C dZ⟩ = ⟨C† rbar, dZ⟩
  var dZbar: Mat1
  mul(dZbar, c.adj, rbar)

  # Step 4: Trace through dB = dU† R0 + U† dR0
  #
  # From dU† R0: ⟨dBbar, dU† R0⟩
  # dU† = Z dx† + dZ X†
  #
  # For Z dx†: ⟨dBbar, Z dx† R0⟩ = Re(tr(dBbar† Z dx† R0))
  #   = Re(tr(R0 dBbar† Z dx†)) = Re(tr((Z† dBbar R0†)† dx†))
  #   = Re(tr(dx (Z† dBbar R0†))) = redot((Z† dBbar R0†)†, dx)
  #   = redot(R0 dBbar† Z, dx)
  # So: dxbar += R0 dBbar† Z
  var tmp2: Mat1
  mul(tmp, r0, dBbar.adj)
  mul(tmp2, tmp, z)
  dxbar += tmp2

  # For dZ X† in dU†: ⟨dBbar, dZ X† R0⟩ = ⟨R0† X dBbar, dZ⟩
  # Actually: Re(tr(dBbar† dZ X† R0)) = Re(tr(R0 dBbar† dZ X†))
  #         = Re(tr(X† R0 dBbar† dZ)) = redot((R0 dBbar†)† X, dZ)
  #         = redot(dBbar R0† X, dZ)? Let me be more careful.
  # ⟨dBbar, dZ X† R0⟩ = Re(tr(dBbar† dZ X† R0))
  # Let A = dBbar†, then tr(A dZ X† R0) = tr(X† R0 A dZ)
  # = redot((X† R0 A)†, dZ) = redot(A† R0† X, dZ) = redot(dBbar R0† X, dZ)
  # Hmm, let me use the identity: redot(A, B) = Re(tr(A† B))
  # redot(dBbar, dZ X† R0) = Re(tr(dBbar† dZ X† R0))
  # We want ⟨?, dZ⟩ form: Re(tr(?† dZ))
  # tr(dBbar† dZ X† R0) = tr((X† R0 dBbar†)† dZ)†... no, that's wrong
  # Actually: tr(A B) = tr(B A), so tr(dBbar† dZ X† R0) = tr(X† R0 dBbar† dZ)
  # = tr((X† R0 dBbar†)† dZ)† = [tr(dZ† (X† R0 dBbar†))]†
  # For real part: Re(tr(X† R0 dBbar† dZ)) = redot((X† R0 dBbar†)†, dZ) = redot(dBbar R0† X, dZ)
  # So: dZbar += dBbar R0† X
  mul(tmp, dBbar, r0.adj)
  mul(tmp2, tmp, x)
  dZbar += tmp2

  # From U† dR0: ⟨dBbar, U† dR0⟩ = ⟨U dBbar, dR0⟩
  # dR0 = C dZ, so ⟨U dBbar, C dZ⟩ = ⟨C† U dBbar, dZ⟩
  mul(tmp, u, dBbar)
  mul(tmp2, c.adj, tmp)
  dZbar += tmp2

  # Step 5: Trace through dZ = -Z dY Z
  # dZ is then numerically symmetrized: dZ_final = 0.5*(dZ + dZ†)
  # But since Z and dY are Hermitian, dZ = -Z dY Z is already Hermitian.
  # So the 0.5 factor is just numerical cleanup and doesn't affect the adjoint.
  # ⟨dZbar, -Z dY Z⟩ = -⟨Z dZbar Z, dY⟩
  var zDzbarZ: Mat1
  mul(tmp, z, dZbar)
  mul(zDzbarZ, tmp, z)
  dYbar -= zDzbarZ

  # Step 6: Trace through dY
  # Y dY + dY Y = dA, then dY is symmetrized: dY_final = 0.5*(dY + dY†)
  # Since dA is Hermitian and Y is Hermitian, dY is already Hermitian.
  # So the 0.5 factor is just numerical cleanup.
  # Sylvester adjoint: Y dAbar + dAbar Y = dYbar
  var dAbar: Mat1
  sylsolve(dAbar, y, dYbar)

  # Step 7: Trace through dA = X† dx + dx† X
  # From X† dx: ⟨dAbar, X† dx⟩ = ⟨X dAbar, dx⟩
  # From dx† X: ⟨dAbar, dx† X⟩ = ⟨X dAbar†, dx⟩ (verified by test!)
  # Combined: dxbar += X dAbar + X dAbar†
  mul(tmp, x, dAbar)
  dxbar += tmp
  mul(tmp, x, dAbar.adj)
  dxbar += tmp

proc projectUHVPVJP_dx*(dxbar: var Mat1, x: Mat2, c: Mat3, rbar: Mat4, eps = 1e-20) =
  ## Convenience wrapper when u is not provided.
  var u {.noInit.}: type(dxbar)
  projectU(u, x)
  projectUHVPVJP_dx(dxbar, u, x, c, rbar, eps)

proc projectUHVPVJP_dc*(dcbar: var Mat1, u: Mat2, x: Mat3, c: Mat4, rbar: auto, eps = 1e-20) =
  ## Adjoint of projectUHVP with respect to dc (chain tangent).
  ## This is identical to the VJP w.r.t. the chain input of projectUVJP.
  projectUVJPChain(dcbar, u, x, rbar, eps)

proc projectUHVPVJP_dc*(dcbar: var Mat1, x: Mat2, c: Mat3, rbar: Mat4, eps = 1e-20) =
  ## Convenience wrapper when u is not provided.
  projectUVJPChain(dcbar, x, rbar, eps)

proc projectSU*(r: var Mat1; x: Mat2) =
  const nc = r.nrows
  var m{.noinit.}: type(r)
  #echo "x: ", x
  m.projectU x
  #echo "m: ", m
  var d = m.determinant    # already unitary: 1=|d
  let p = (1.0/float(-nc)) * atan2(d.im, d.re)
  d.re = cos p
  d.im = sin p
  #echo "d: ", d
  r := d * m

template projectSU*(r: var Mat1) =
  r.projectSU r

proc projectTAH*(r: var Mat1; x: Mat2) =
  r := 0.5*(x-x.adj)
  const nc = x.nrows
  when nc > 1:
    let d = r.trace / nc.float
    r -= d

template projectTAH*(r: var Mat1) =
  var t{.noInit.}: evalType(r)
  t := r
  r.projectTAH t

proc checkU*(x: Mat1): auto {.inline, noinit.} =
  ## Returns the sum of deviations of x^dag x and det(x) from unitarity.
  var d = norm2(-1.0 + x.adj * x)
  return d

proc checkSU*(x: Mat1): auto {.inline, noinit.} =
  ## Returns the sum of deviations of x^dag x and det(x) from unitarity.
  var d = norm2(-1.0 + x.adj * x)
  d += norm2(-1.0 + x.determinant)
  return d

#[
template rsqrtM2(r:typed; x:typed) =
  load(x00, x[0,0].re)
  load(x01, x[0,1])
  #load(x10, x[1,0])
  load(x11, x[1,1].re)
  let det := a00*a11 -
  QLA_r_eq_Re_c_times_c (det, a00, a11);
  QLA_r_meq_Re_c_times_c(det, a01, a10);
  tr = QLA_real(a00) + QLA_real(a11);
  sdet = sqrtP(fabsP(det));
  // c0 = (l2/sl1-l1/sl2)/(l2-l1) = (l2+sl1*sl2+l1)/(sl1*sl2*(sl1+sl2))
  // c1 = (1/sl2-1/sl1)/(l2-l1) = -1/(sl1*sl2*(sl1+sl2))
  c1 = 1/(sdet*sqrtP(fabsP(tr+2*sdet)));
  c0 = (tr+sdet)*c1;
  c1 = -c1;
  // c0 + c1*a
  QLA_c_eq_c_times_r_plus_r(QLA_elem_M(*r,0,0), a00, c1, c0);
  QLA_c_eq_c_times_r(QLA_elem_M(*r,0,1), a01, c1);
  QLA_c_eq_c_times_r(QLA_elem_M(*r,1,0), a10, c1);
  QLA_c_eq_c_times_r_plus_r(QLA_elem_M(*r,1,1), a11, c1, c0);

template rsqrtM(r:typed; x:typed) =
  assert(r.nrows == x.nrows)
  assert(r.ncols == x.ncols)
  assert(r.nrows == r.ncols)
  if r.nrows==1:
    rsqrt(r[0,0], x[0,0])
  elif r.nrows==2:
    rsqrtM2(r, x)
  elif r.nrows==3:
    rsqrtM3(r, x)
  else:
    echo "unimplemented"
    quit(1)
proc rsqrt(r:var Mat1; x:Mat2) = rsqrt(r, x)
]#

proc exp*(m: Mat1): auto {.noInit.} =
  var r{.noInit.}: MatrixArray[m.nrows,m.ncols,type(m[0,0])]
  when m.nrows == 1:
    r := exp(m[0,0])
  else:
    #r := expPoly12(m)
    #r := expPade4(m)
    #r := expPade8(m)
    #r := expPade9(m)
    var p: ExpParam
    p.scale = 20
    p.kind = ekPoly
    p.order = 4
    r := p.exp(m)
    #[
    type ft = numberType(m)
    template term(n,x: typed): untyped =
      when x.type is nil.type: 1 + ft(n)*m
      else: 1 + ft(n)*m*x
    #template r3:untyped = nil
    let r12 = term(1.0/12.0, nil)
    let r11 = term(1.0/11.0, r12)
    let r10 = term(1.0/10.0, r11)
    let r9 = term(1.0/9.0, r10)
    let r8 = term(1.0/8.0, r9)
    let r7 = term(1.0/7.0, r8)
    let r6 = term(1.0/6.0, r7)
    let r5 = term(1.0/5.0, r6)
    let r4 = term(1.0/4.0, r5)
    let r3 = term(1.0/3.0, r4)
    let r2 = term(1.0/2.0, r3)
    r := 1 + m*r2
    ]#
  r

proc expDeriv*(m: Mat1, c:Mat2): auto {.noInit.} =
  var r{.noInit.}: MatrixArray[m.nrows,m.ncols,type(m[0,0])]
  when m.nrows == 1:
    r := exp(m[0,0]) * c[0,0]
  else:
    var p: ExpParam
    p.scale = 20
    p.kind = ekPoly
    p.order = 4
    r := p.expDeriv(m, c)
  r

proc ln*(m: Mat1): auto {.noInit.} =
  var r{.noInit.}: MatrixArray[m.nrows,m.ncols,type(m[0,0])]
  when m.nrows == 1:
    r := ln(m[0,0])
  else:
    static: error("ln of matrix not implimented.")
  r

proc re*(m: Mat1): auto {.noInit.} =
  var r{.noInit.}: MatrixArray[m.nrows,m.ncols,type(m[0,0])]
  for i in 0..<m.nrows:
    for j in 0..<m.ncols:
      r[i,j] := re(m[i,j])
  r
proc im*(m: Mat1): auto {.noInit.} =
  var r{.noInit.}: MatrixArray[m.nrows,m.ncols,type(m[0,0])]
  for i in 0..<m.nrows:
    for j in 0..<m.ncols:
      r[i,j] := im(m[i,j])
  r

when isMainModule:
  import macros
  import simd
  template `+`(x: SimdS4): untyped = x
  template `+`(x: SimdS4, y: ComplexType): untyped =
    asReal(x) + y
  template `-`(x: SimdS4, y: ComplexType): untyped =
    asReal(x) - y
  template `*`(x: SimdS4, y: ComplexType): untyped =
    asReal(x) * y
  template `*`(x: ComplexType, y: SimdS4): untyped =
    x * asReal(y)
  template `/`(x: SomeFloat, y: SomeInteger): untyped = x/(type(x))(y)
  template add(r: ComplexType, x: SimdS4, y: ComplexType): untyped =
    add(r, asReal(x), y)
  template sub(r: ComplexType, x: SimdS4, y: ComplexType): untyped =
    sub(r, asReal(x), y)
  template mul(r: ComplexType, x: SimdS4, y: ComplexType): untyped =
    mul(r, asReal(x), y)
  template check(x:untyped, n:SomeNumber):untyped =
    let r0 = x
    let r = simdSum(r0)/simdLength(r0)
    echo "error/eps: ", r/epsilon(r)
    doAssert(abs(r)<n*epsilon(r))
  proc test(T: typedesc) =
    var m1,m2,m3,m4: T
    let N = m1.nrows
    for i in 0..<N:
      for j in 0..<N:
        let fi = i.float
        let fj = j.float
        m1[i,j].re := 0.5 + 0.7/(0.9+1.3*fi-fj)
        m1[i,j].im := 0.1 + 0.3/(0.4+fi-1.1*fj)
    echo "test " & $N & " " & $T
    #echo "m1: ", m1
    m2 := m1.adj * m1
    #echo m2
    rsqrtPH(m3, m2)
    #echo m3
    m4 := m3*m2*m3
    let err = sqrt((1-m4).norm2)/N
    #echo m4
    echo " rsqrtPH err: ", err
    check(err, 50)

    projectU(m2, m1)
    m3 := m2.adj*m2
    var err2 = sqrt((1-m3).norm2)/N
    echo " projectU err: ", err2
    check(err, 50)

    #m2 := 0.1*(m2 - (trace(m2)/N))
    let m2n = 1/(10*sqrt(m2.norm2))
    m3 := exp(m2n*m2)
    m4 := exp(-m2n*m2)
    m2 := m3*m4
    #echo "exp ",m2,"\n\t= ",m3
    err2 = sqrt((1-m2).norm2/(N*N))
    echo " exp err: ", err2
    check(err2, 5)

    inverse(m4, m1)
    m3 := m1*m4
    #echo m3
    err2 = sqrt((1-m3).norm2/(N*N))
    echo " inverse err: ", err2
    check(err2, 5)

    if N<5:
      projectU(m2, m1)
      let r1 = trace(m4.adj*m2).re
      let seps = sqrt(epsilon(simdSum(r1)))
      m3 := m1 + 3*seps*m1*m1
      #m3 := m1 + 1e-3'f32*m1*m1
      projectU(m2, m3)
      let r2 = trace(m4.adj*m2).re
      projectUderiv(m2, m1, m4)
      let dr = r2 - r1
      let dm = trace((m3-m1).adj * m2).re
      #echo " r1: ", r1
      #echo " r2: ", r2
      #echo " dr: ", dr
      #echo " dm: ", dm
      #echo "m1: ", m1
      #echo "m3: ", m3
      let dd = abs(dr - dm)
      echo " projectUderiv err: ", dd
      #doAssert(simdSum(dd)<simdLength(dd)*N*eps*40)
      check(dd, 20*N)


  type
    Cmplx[T] = ComplexType[T]
    CM[N:static[int],T] = MatrixArray[N,N,Cmplx[T]]
  template doTest(t:untyped) =
    when declared(t):
      test(CM[1,t])
      test(CM[2,t])
      test(CM[3,t])
      test(CM[4,t])
  doTest(float32)
  doTest(float64)
  doTest(SimdS4)
  doTest(SimdD4)
  doTest(SimdS8)
  doTest(SimdD8)
