import base
import complexNumbers
import matrixConcept
import matinv
import matrixFunctions
import types
import strformat
getOptimPragmas()

proc adjugate*(r: var Mat1, x: Mat2) {.alwaysInline.} =
  const nc = r.nrows
  when nc==1:
    r := 1
  elif nc==2:
    r[0,0] :=  x[1,1]
    r[0,1] := -x[0,1]
    r[1,0] := -x[1,0]
    r[1,1] :=  x[0,0]
  elif nc==3:
    let x00 = x[0,0]
    let x01 = x[0,1]
    let x02 = x[0,2]
    let x10 = x[1,0]
    let x11 = x[1,1]
    let x12 = x[1,2]
    let x20 = x[2,0]
    let x21 = x[2,1]
    let x22 = x[2,2]
    r[0,0] := x11*x22 - x12*x21
    r[0,1] := x21*x02 - x22*x01
    r[0,2] := x01*x12 - x02*x11
    r[1,0] := x12*x20 - x10*x22
    r[1,1] := x22*x00 - x20*x02
    r[1,2] := x02*x10 - x00*x12
    r[2,0] := x10*x21 - x11*x20
    r[2,1] := x20*x01 - x21*x00
    r[2,2] := x00*x11 - x01*x10
  else:
    echo &"adjugate n({nc})>3 not supported"
    doAssert(false)

proc sylsolveN*(x: var Mat1, a0: Mat2, c0: Mat3) =
  mixin simdMax
  let nc = x.nrows
  let a2 = a0.norm2
  let ia = rsqrt(a2)

  x := (0.5*ia)*c0
  var a = ia*a0

  let rstop = epsilon(x.norm2.simdMax)
  let maxit = 20
  var nit = 0
  while true:
    #echo nit, ": ", x
    inc nit

    var v = 3 - a*a
    var w = a*x - x*a
    var d = x*v
    x := 0.5*(d-a*w)
    a := 0.5*a*v

    let r = c0-(a0*x+x*a0)
    let rnorm = r.norm2.simdMax
    #echo nit, " r2: ", rnorm
    if nit>=maxit or rnorm<rstop:
      if rnorm>rstop:
        echo "WARNING sylsolveN failed to converge: ", nit, " r2: ", rnorm
      break

proc sylsolveN2*(x: var Mat1, a: Mat2, c: Mat3) =
  mixin simdMax
  let nc = x.nrows
  let aa = a.adj
  let t2 = a.norm2
  let kappa = 0.5/t2

  x := c
  #x := kappa * (aa*c + c*aa)

  let rstop = epsilon(x.norm2.simdMax)
  let maxit = 50
  var nit = 0
  while true:
    inc nit

    let r = c - a*x - x*a
    x += kappa * (aa*r + r*aa)

    let rnorm = r.norm2.simdMax
    echo nit, " r2: ", rnorm
    if nit>=maxit or rnorm<rstop:
      if rnorm>rstop:
        echo "WARNING sylsolveN failed to converge: ", nit, " r2: ", rnorm
      break

proc sylsolve*(x: var Mat1, a: Mat2, c: Mat3) =
  ## solves A X + X A = C for X
  const nc = x.nrows
  when nc==1:
    x[0,0] := c[0,0] / (2*a[0,0])
  elif nc==2:
    # x = (C + |A| A^-1 C A^-1)/2Tr(A)
    let a00 = a[0,0]
    let a01 = a[0,1]
    let a10 = a[1,0]
    let a11 = a[1,1]
    let c00 = c[0,0]
    let c01 = c[0,1]
    let c10 = c[1,0]
    let c11 = c[1,1]
    let idet = 1/(a00*a11 - a01*a10)
    let itr = 0.5/(a00 + a11)
    # ai = [[a11,-a01][-a10,a00]]
    let aic00 = a11 * c00 - a01 * c10
    let aic01 = a11 * c01 - a01 * c11
    let aic10 = a00 * c10 - a10 * c00
    let aic11 = a00 * c11 - a10 * c01
    x[0,0] := itr * (c00 + idet * (aic00*a11-aic01*a10))
    x[0,1] := itr * (c01 + idet * (aic01*a00-aic00*a01))
    x[1,0] := itr * (c10 + idet * (aic10*a11-aic11*a10))
    x[1,1] := itr * (c11 + idet * (aic11*a00-aic10*a01))
  elif nc==3:
    var ad {.noInit.}: type(a)
    adjugate(ad, a)
    let t = a[0,0] + a[1,1] + a[2,2]
    let s = ad[0,0] + ad[1,1] + ad[2,2]
    let r = a[0,0]*ad[0,0] + a[0,1]*ad[1,0] + a[0,2]*ad[2,0]
    var ac {.noInit.}: type(a)
    var ca {.noInit.}: type(a)
    var aca {.noInit.}: type(a)
    var adc {.noInit.}: type(a)
    var cad {.noInit.}: type(a)
    var adcad {.noInit.}: type(a)
    mul(ac, a, c)
    mul(ca, c, a)
    mul(aca, ac, a)
    mul(adc, ad, c)
    mul(cad, c, ad)
    mul(adcad, adc, ad)
    let c2 = 1/(2*(s*t-r))
    let c0 = c2*(s+t*t)
    let c1 = c2*(t/r)
    let c4 = c2*(t)
    for i in 0..2:
      for j in 0..2:
        x[i,j] := c0*c[i,j] + c1*adcad[i,j] + c2*(aca[i,j]-adc[i,j]-cad[i,j]) -
                  c4*(ac[i,j]+ca[i,j])
  else:
    sylsolveN(x, a, c)

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
  mul(r0, c, z)            # R0 = C Z
  mul(t1, u.adj, r0)       # B = U† R0
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
  duAdj += tmp  # full dU† including dZ contribution

  # dB pieces
  # Full variation: dB = dU† R0 + U† dR0
  # However, the dR0 contribution through U† dR0 is separate from the direct dR0 term.
  # We include both here as the Sylvester equation for dS needs the full dB.
  var db: Mat1
  mul(db, duAdj, r0)     # dU† R0
  mul(tmp, u.adj, dr0)   # U† dR0
  db += tmp              # full dB

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
  mul(chainbar, rbar, z)  # rbar*Z (Z is Hermitian, so Z† = Z)
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
