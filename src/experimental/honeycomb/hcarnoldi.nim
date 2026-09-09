## Krylov-Schur (thick-restart) Arnoldi for general complex operators.
##
##   A V_m = V_m G + beta_m v_{m+1} e_m^H          (CGS2 + DGKS, orthonormal to ~1e-14)
## Restart: Schur form of G (zgees), ztrsen moves the kKeep best Ritz values
## (per `which`) to the leading block, the basis is rotated to those Schur
## vectors and the old residual becomes basis vector kKeep; the projected
## matrix keeps a spike row below the retained triangle.  Converged
## directions lock softly (their spike entries are ~0).  Convergence:
## beta_m |y_{m-1}| <= tol max(|theta|, eps^{2/3} |G|).
##
## Vectors are handled through the ArnoldiOp closures and six mixin
## operations in scope at the instantiation site (hcspec, tests/tarnoldi):
##   vcopy(dst, src)  vzero(v)  vscale(v, s)  vaxpy(y, a, x)
##   vdot(x, y) = x^H y (conjugate-linear in the first argument)  vnorm2(x)
## Deterministic: seeded start vectors, LAPACK, no threading here.
## LAPACK link options: -d:hcLapackLib="..." (macOS default Accelerate).

import std/[complex, math, algorithm, strformat]

# ---------------------------------------------------------------------------
# LAPACK bindings (self-contained; see module docs for the link options)
# ---------------------------------------------------------------------------

when defined(macosx):
  const hcLapackLib* {.strdefine.} = "-framework Accelerate"
else:
  const hcLapackLib* {.strdefine.} = "-llapack -lblas"
{.passL: hcLapackLib.}

type lpint* = cint  ## LAPACK integer (32-bit in Accelerate's default LAPACK)

proc zgeev(jobvl, jobvr: cstring; n: ptr lpint; a: ptr Complex64;
           lda: ptr lpint; w: ptr Complex64; vl: ptr Complex64;
           ldvl: ptr lpint; vr: ptr Complex64; ldvr: ptr lpint;
           work: ptr Complex64; lwork: ptr lpint; rwork: ptr float64;
           info: ptr lpint) {.importc: "zgeev_".}

proc zgees(jobvs, srt: cstring; select: pointer; n: ptr lpint;
           a: ptr Complex64; lda: ptr lpint; sdim: ptr lpint;
           w: ptr Complex64; vs: ptr Complex64; ldvs: ptr lpint;
           work: ptr Complex64; lwork: ptr lpint; rwork: ptr float64;
           bwork: pointer; info: ptr lpint) {.importc: "zgees_".}

proc ztrsen(job, compq: cstring; select: ptr lpint; n: ptr lpint;
            t: ptr Complex64; ldt: ptr lpint; q: ptr Complex64;
            ldq: ptr lpint; w: ptr Complex64; m: ptr lpint; s: ptr float64;
            sep: ptr float64; work: ptr Complex64; lwork: ptr lpint;
            info: ptr lpint) {.importc: "ztrsen_".}

# ---------------------------------------------------------------------------
# small dense complex matrices (column-major, LAPACK-compatible)
# ---------------------------------------------------------------------------

type ZMat* = object
  nr*, nc*: int
  d*: seq[Complex64]        ## column-major: [i + j*nr]

proc newZMat*(nr, nc: int): ZMat =
  ZMat(nr: nr, nc: nc, d: newSeq[Complex64](nr*nc))

template `[]`*(a: ZMat; i, j: int): Complex64 = a.d[i + j*a.nr]
template `[]=`*(a: ZMat; i, j: int; x: Complex64) = a.d[i + j*a.nr] = x

proc zeroit*(a: var ZMat) =
  for e in mitems(a.d): e = complex64(0.0, 0.0)

proc frobNorm*(a: ZMat): float =
  var s = 0.0
  for e in a.d: s += abs2(e)
  sqrt(s)

proc zeig*(a: ZMat; wantV: bool = false): tuple[w: seq[Complex64]; v: ZMat] =
  ## Eigenvalues (and, if wantV, unit-normalized right eigenvectors) of a
  ## general complex matrix via LAPACK zgeev.  `a` is left untouched.
  doAssert a.nr == a.nc and a.nr >= 1
  var n = lpint(a.nr)
  var ac = a                       # zgeev destroys its input
  result.w = newSeq[Complex64](a.nr)
  var vp: ptr Complex64 = nil
  if wantV:
    result.v = newZMat(a.nr, a.nr)
    vp = addr result.v.d[0]
  var lwork = lpint(8*a.nr + 16)
  var work = newSeq[Complex64](lwork)
  var rwork = newSeq[float64](2*a.nr)
  var info: lpint
  let jobvr: cstring = if wantV: cstring"V" else: cstring"N"
  zgeev("N", jobvr, addr n, addr ac.d[0], addr n,
        addr result.w[0], nil, addr n, vp, addr n,
        addr work[0], addr lwork, addr rwork[0], addr info)
  doAssert info == 0, "zgeev failed, info = " & $info

proc zschur*(a: ZMat): tuple[t, q: ZMat; w: seq[Complex64]] =
  ## Complex Schur decomposition a = q t q^H (zgees, no sorting).
  doAssert a.nr == a.nc and a.nr >= 1
  var n = lpint(a.nr)
  result.t = a
  result.q = newZMat(a.nr, a.nr)
  result.w = newSeq[Complex64](a.nr)
  var sdim: lpint
  var lwork = lpint(8*a.nr + 16)
  var work = newSeq[Complex64](lwork)
  var rwork = newSeq[float64](a.nr)
  var info: lpint
  zgees("V", "N", nil, addr n, addr result.t.d[0], addr n, addr sdim,
        addr result.w[0], addr result.q.d[0], addr n,
        addr work[0], addr lwork, addr rwork[0], nil, addr info)
  doAssert info == 0, "zgees failed, info = " & $info

proc ztrsenReorder*(t, q: var ZMat; keep: openArray[int]): seq[Complex64] =
  ## Reorder the Schur form so the eigenvalues at positions `keep` lead;
  ## returns the reordered eigenvalue list.  t, q updated in place.
  doAssert t.nr == t.nc and q.nr == t.nr
  var n = lpint(t.nr)
  result = newSeq[Complex64](t.nr)
  var m: lpint
  var s, sep: float64
  var lwork = lpint(t.nr*t.nr + 16)
  var work = newSeq[Complex64](lwork)
  var sel = newSeq[lpint](t.nr)
  for i in keep: sel[i] = 1
  var info: lpint
  ztrsen("N", "V", addr sel[0], addr n, addr t.d[0], addr n,
         addr q.d[0], addr n, addr result[0], addr m, addr s, addr sep,
         addr work[0], addr lwork, addr info)
  doAssert info == 0, "ztrsen failed, info = " & $info

proc solveShifted(g: ZMat; lam: Complex64; rhs: seq[Complex64]): seq[Complex64] =
  ## Solve (g - lam I) x = rhs by LU with partial pivoting.  Near-zero pivots
  ## are floored at eps*|g| (inverse-iteration convention), so a (near-)
  ## singular shift returns a highly amplified null-ish direction.
  let n = g.nr
  var a = g
  for i in 0..<n: a[i, i] = a[i, i] - lam
  let tiny = 2.3e-16*frobNorm(a)
  result = rhs
  for k in 0..<n:
    var p = k
    var mx = abs2(a[k, k])
    for i in k+1..<n:
      let v = abs2(a[i, k])
      if v > mx:
        mx = v
        p = i
    if p != k:
      for j in 0..<n: swap(a.d[k + j*n], a.d[p + j*n])
      swap(result[k], result[p])
    var akk = a[k, k]
    if abs(akk) < tiny: akk = complex64(tiny, 0.0)
    a[k, k] = akk
    for i in k+1..<n:
      let l = a[i, k]/akk
      for j in k+1..<n: a[i, j] = a[i, j] - l*a[k, j]
      result[i] = result[i] - l*result[k]
  for k in countdown(n-1, 0):
    var s = result[k]
    for j in k+1..<n: s -= a[k, j]*result[j]
    result[k] = s/a[k, k]

proc eigvecResid*(g: ZMat; lam: Complex64; y: seq[Complex64]): float =
  ## |g y - lam y|
  var s = 0.0
  for i in 0..<g.nr:
    var t = -lam*y[i]
    for j in 0..<g.nr: t += g[i, j]*y[j]
    s += abs2(t)
  sqrt(s)

proc refineEigvec*(g: ZMat; lam: Complex64; y0: seq[Complex64]): seq[Complex64] =
  ## Improve a computed eigenvector of g by (up to two steps of) inverse
  ## iteration with the computed eigenvalue.  zgeev/ztrevc eigenvectors can
  ## have residuals ~eps/gap when eigenvalues are clustered — which happens
  ## systematically here, since soft-locked converged Ritz values reappear
  ## in the active subspace (roundoff re-seeds them; for truly degenerate
  ## operators the copies are real and wanted).  Inverse iteration recovers
  ## an eps-residual vector in the cluster.  Returns the best vector found
  ## (unit norm), never worse than y0.
  result = y0
  var best = eigvecResid(g, lam, y0)
  var cur = y0
  for it in 0..1:
    var z = solveShifted(g, lam, cur)
    var zn = 0.0
    for e in z: zn += abs2(e)
    if zn <= 0.0: break
    let sc = 1.0/sqrt(zn)
    for e in mitems(z): e = sc*e
    let r = eigvecResid(g, lam, z)
    cur = z
    if r < best:
      best = r
      result = z
    else:
      break

# ---------------------------------------------------------------------------
# eigenvalue ordering
# ---------------------------------------------------------------------------

proc eigOrder*(w: openArray[Complex64]; which: string): seq[int] =
  ## Indices of w sorted best-first according to `which`:
  ## "LM"/"SM" largest/smallest |w|, "LR"/"SR" largest/smallest Re w.
  ## Deterministic: ties broken by (re, im).
  result = newSeq[int](w.len)
  for i in 0..<w.len: result[i] = i
  proc cmpIdx(w: openArray[Complex64]; which: string; a, b: int): int =
    var d = case which
      of "LM": cmp(abs2(w[b]), abs2(w[a]))
      of "SM": cmp(abs2(w[a]), abs2(w[b]))
      of "LR": cmp(w[b].re, w[a].re)
      of "SR": cmp(w[a].re, w[b].re)
      else: 0
    if d == 0: d = cmp(w[a].re, w[b].re)
    if d == 0: d = cmp(w[a].im, w[b].im)
    d
  let wl = @w
  result.sort(proc(a, b: int): int = cmpIdx(wl, which, a, b))

# ---------------------------------------------------------------------------
# the eigensolver
# ---------------------------------------------------------------------------

type ArnoldiOp*[V] = object
  apply*: proc (r: var V; x: V)    ## r = A x
  newVec*: proc (): V              ## fresh zeroed vector
  start*: proc (v: var V)          ## fill v with a (seeded) start vector;
                                   ## also used to replace a breakdown vector

proc arnoldi*[V](op: ArnoldiOp[V]; nev, ncv: int; tol: float;
                 maxRestarts: int; which: string; verb: int = 0):
    tuple[vals: seq[Complex64]; vecs: seq[V]; resids: seq[float]; nApply: int] =
  ## Compute the `nev` eigenvalues of A selected by `which`
  ## ("LM", "SM", "LR", "SR"), with a Krylov subspace of dimension `ncv`
  ## and at most `maxRestarts` restart cycles.
  ##
  ## Returns the `nev` best Ritz pairs found (best-first in the `which`
  ## ordering).  `resids[i] = |A v_i - vals_i v_i| / max(|vals_i|, floor)`
  ## is computed DIRECTLY with an extra operator application per pair
  ## (floor = macheps^(2/3) * |G|_F guards the lambda ~ 0 case);
  ## `nApply` counts every call of `op.apply`.
  ## Convergence of pair i inside the iteration uses the Arnoldi estimate
  ## `beta_m |e_m^H y_i| <= tol * max(|lambda_i|, floor)`.
  mixin vcopy, vzero, vscale, vaxpy, vdot, vnorm2
  doAssert which in ["LM", "SM", "LR", "SR"], "arnoldi: bad which=" & which
  doAssert nev >= 1 and ncv >= nev + 2, "arnoldi: need ncv >= nev+2"
  doAssert op.apply != nil and op.newVec != nil and op.start != nil
  let m = ncv
  let kKeep = nev + (ncv - nev) div 2      # thick-restart size, nev < kKeep < ncv
  let eps23 = pow(2.220446049250313e-16, 2.0/3.0)

  var vb = newSeq[V](m + 1)                # Krylov basis, one seq, no copies
  for i in 0..m: vb[i] = op.newVec()
  var w = op.newVec()                      # work vector
  var scr: seq[V]                          # restart rotation scratch (lazy)
  var g = newZMat(m, m)                    # projected matrix
  var betaM = 0.0                          # residual coefficient
  var nApply = 0

  op.start(vb[0])
  let n0 = vnorm2(vb[0])
  doAssert n0 > 0.0, "arnoldi: start vector is zero"
  vscale(vb[0], 1.0/sqrt(n0))

  # -- Arnoldi expansion of columns j0..m-1 with CGS2 --
  template expandFrom(j0: int) =
    for j {.inject.} in j0..<m:
      op.apply(w, vb[j])
      inc nApply
      let wn0 = sqrt(vnorm2(w))
      # classical Gram-Schmidt at least twice (CGS2); if the norm keeps
      # collapsing (severe cancellation, e.g. w almost in span(V)), repeat
      # up to twice more (DGKS criterion) so orthogonality stays at eps.
      var bta = 0.0
      var wprev = wn0
      for pass in 0..3:
        for i in 0..j:
          let c = vdot(vb[i], w)
          vaxpy(w, -c, vb[i])
          g[i, j] = g[i, j] + c
        bta = sqrt(vnorm2(w))
        if pass >= 1 and (bta > 0.5*wprev or bta <= 1e-14*wn0): break
        wprev = bta
      if bta <= 1e-13*wn0:
        # (near-)invariant subspace: true residual is zero; continue the
        # basis with a fresh orthonormalized direction and a 0 coupling.
        if verb > 0: echo &"arnoldi: breakdown at column {j}, deflating"
        op.start(w)
        for pass in 0..1:
          for i in 0..j:
            let c = vdot(vb[i], w)
            vaxpy(w, -c, vb[i])
        bta = sqrt(vnorm2(w))
        doAssert bta > 1e-12, "arnoldi: Krylov space exhausted; reduce ncv"
        vscale(w, 1.0/bta)
        vcopy(vb[j+1], w)
        if j < m-1: g[j+1, j] = complex64(0.0, 0.0) else: betaM = 0.0
      else:
        vscale(w, 1.0/bta)
        vcopy(vb[j+1], w)
        if j < m-1: g[j+1, j] = complex64(bta, 0.0) else: betaM = bta

  # -- Ritz extraction and convergence count --
  var rw: seq[Complex64]                   # Ritz values of g
  var ry: ZMat                             # right eigenvectors of g
  var rord: seq[int]                       # sorted (best-first) indices
  var floorv = 0.0
  var nconv = 0
  template ritzEval() =
    (rw, ry) = zeig(g, true)
    rord = eigOrder(rw, which)
    floorv = eps23*frobNorm(g)
    nconv = 0
    for i in 0..<nev:
      let idx = rord[i]
      let est = betaM * abs(ry[m-1, idx])  # zgeev vectors have |y| = 1
      if est <= tol * max(abs(rw[idx]), floorv): inc nconv

  # -- Krylov-Schur restart: keep the kKeep best Ritz values --
  template restartStep() =
    var (t, q, sw) = zschur(g)
    let sord = eigOrder(sw, which)
    discard ztrsenReorder(t, q, sord[0..<kKeep])
    var spike = newSeq[Complex64](kKeep)
    for j in 0..<kKeep: spike[j] = betaM * q[m-1, j]
    if scr.len < kKeep:
      let o = scr.len
      scr.setLen(kKeep)
      for i in o..<kKeep: scr[i] = op.newVec()
    for j in 0..<kKeep:                    # basis <- basis * Q(:, 0..kKeep-1)
      vzero(scr[j])
      for i in 0..<m:
        vaxpy(scr[j], q[i, j], vb[i])
    for j in 0..<kKeep: swap(vb[j], scr[j])
    swap(vb[kKeep], vb[m])                 # old residual vector -> v_kKeep
    g.zeroit
    for j in 0..<kKeep:
      for i in 0..j: g[i, j] = t[i, j]     # retained triangular block
      g[kKeep, j] = spike[j]               # the spike row

  expandFrom(0)
  ritzEval()
  var cycles = 0
  while nconv < nev and cycles < maxRestarts:
    inc cycles
    restartStep()
    expandFrom(kKeep)
    ritzEval()
    if verb > 1 or (verb == 1 and cycles mod 50 == 0):
      let i0 = rord[0]
      echo &"arnoldi cycle {cycles}: nconv {nconv}/{nev}, nApply {nApply}, " &
           &"best ({rw[i0].re:.6g},{rw[i0].im:.6g}) " &
           &"est {betaM*abs(ry[m-1,i0]):.3e}"

  if verb > 0:
    echo &"arnoldi({which}): nev {nev}, ncv {ncv}, kKeep {kKeep}: " &
         &"{nconv}/{nev} converged after {cycles} cycles, {nApply} applies"

  if verb >= 3:  # expensive diagnostics: basis orthonormality + factorization
    var oerr = 0.0
    for i in 0..m:
      for j in i..m:
        var d = vdot(vb[i], vb[j])
        if i == j: d = d - complex64(1.0, 0.0)
        oerr = max(oerr, abs(d))
    var ferr = 0.0
    for j in 0..<m:
      op.apply(w, vb[j])
      inc nApply
      for i in 0..<m:
        vaxpy(w, -g[i, j], vb[i])
      if j == m-1: vaxpy(w, complex64(-betaM, 0.0), vb[m])
      ferr = max(ferr, sqrt(vnorm2(w)))
    echo &"arnoldi diag: orthonormality err {oerr:.3e}, " &
         &"factorization err {ferr:.3e}"

  # -- final: Ritz vectors and DIRECT residuals (one extra apply each) --
  result.vals = newSeq[Complex64](nev)
  result.vecs = newSeq[V](nev)
  result.resids = newSeq[float](nev)
  for i in 0..<nev:
    let idx = rord[i]
    let lam = rw[idx]
    var y = newSeq[Complex64](m)
    for j in 0..<m: y[j] = ry[j, idx]
    y = refineEigvec(g, lam, y)
    var x = op.newVec()
    vzero(x)
    for j in 0..<m:
      vaxpy(x, y[j], vb[j])
    let xn = vnorm2(x)
    doAssert xn > 0.0
    vscale(x, 1.0/sqrt(xn))
    op.apply(w, x)
    inc nApply
    vaxpy(w, -lam, x)
    result.vals[i] = lam
    result.vecs[i] = x
    result.resids[i] = sqrt(vnorm2(w)) / max(abs(lam), floorv)
    if verb >= 3:
      echo &"  pair {i}: lam ({lam.re:.9g},{lam.im:.9g}) " &
           &"est {betaM*abs(ry[m-1, idx]):.3e} direct {result.resids[i]:.3e} " &
           &"|gy-ly| {eigvecResid(g, lam, y):.3e} |y_m| {abs(y[m-1]):.3e}"
  result.nApply = nApply

