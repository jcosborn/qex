## Dense complex linear algebra on column-major square matrices, the LAPACK
## layout: element (i, j) of an n x n matrix `a` is a[i + n*j], and Complex64
## has the memory layout of dcomplex.  Shared by the dense oracles of the tests,
## the free-limit Matsubara pipelines (rfree, rspec) and the dense diagnostics
## of the operators; nothing here is used inside a solver loop.

import std/[complex, math]
import eigens/lapack
import eigens/linalgFuncs

export complex

proc zmm*(ta, tb: cstring, n: int, a, b: seq[Complex64], c: var seq[Complex64]) =
  ## c = op(a) op(b) with op = "N" (as is) or "C" (conjugate transpose), BLAS zgemm.
  if c.len != n*n: c = newSeq[Complex64](n*n)
  var
    nn = fint(n)
    one = dcomplex(re: 1.0, im: 0.0)
    zero = dcomplex(re: 0.0, im: 0.0)
  zgemm(ta, tb, addr nn, addr nn, addr nn, addr one,
        cast[ptr dcomplex](unsafeAddr a[0]), addr nn,
        cast[ptr dcomplex](unsafeAddr b[0]), addr nn,
        addr zero, cast[ptr dcomplex](addr c[0]), addr nn)

proc zmm*(a, b: seq[Complex64], n: int): seq[Complex64] =
  ## a b
  zmm("N", "N", n, a, b, result)

proc zmmAdjL*(a, b: seq[Complex64], n: int): seq[Complex64] =
  ## a^dag b
  zmm("C", "N", n, a, b, result)

proc zsolve*(a: var seq[Complex64], b: var seq[Complex64], n, nrhs: int) =
  ## Solve a x = b with LAPACK zgesv.  `a` is overwritten by its LU factors,
  ## `b` (n x nrhs) by the solution.
  var
    nn = fint(n)
    nr = fint(nrhs)
    ipiv = newSeq[fint](n)
    info = fint(0)
  zgesv(addr nn, addr nr, cast[ptr dcomplex](addr a[0]), addr nn,
        addr ipiv[0], cast[ptr dcomplex](addr b[0]), addr nn, addr info)
  doAssert info == 0, "zgesv: info = " & $info

proc zinv*(a: seq[Complex64], n: int): seq[Complex64] =
  ## a^{-1} by LU (zgesv against the identity).
  var lu = a
  result = newSeq[Complex64](n*n)
  for i in 0..<n: result[i + n*i] = complex64(1.0, 0.0)
  zsolve(lu, result, n, n)

proc eigvals*(a: seq[Complex64], n: int): seq[Complex64] =
  ## Eigenvalues of a general matrix (zgeev on a copy).
  var m = a
  result = newSeq[Complex64](n)
  zgeigs(cast[ptr float64](addr m[0]), cast[ptr float64](addr result[0]), n)

proc heig*(a: var seq[Complex64], n: int): seq[float] =
  ## Hermitian eigenproblem (zheev): `a` is replaced by its eigenvectors, one per
  ## column, and the ascending eigenvalues are returned.
  result = newSeq[float](n)
  zeigs(cast[ptr float64](addr a[0]), addr result[0], n)

proc sigmaBounds*(x: seq[Complex64], n: int): tuple[smin, smax: float] =
  ## Extreme singular values of x from the eigenvalues of x^dag x.
  var h = zmmAdjL(x, x, n)
  let ev = heig(h, n)
  (sqrt(ev[0]), sqrt(ev[n-1]))

proc ovFromXInto*(x: var seq[Complex64], n: int,
                  h, w, g: var seq[Complex64], ev: var seq[float]) =
  ## x <- 1 + x (x^dag x)^{-1/2}, the overlap operator of the kernel x, through
  ## the exact eigendecomposition of x^dag x.  h, w, g, ev are scratch, grown on
  ## demand and reusable across calls of the same dimension.
  zmm("C", "N", n, x, x, h)                            # h = X^dag X
  if ev.len != n: ev = newSeq[float](n)
  zeigs(cast[ptr float64](addr h[0]), addr ev[0], n)   # h <- eigenvectors V
  if w.len != n*n: w = newSeq[Complex64](n*n)
  for j in 0..<n:
    let s = 1.0/sqrt(ev[j])
    for i in 0..<n: w[i + n*j] = s*h[i + n*j]          # w = V E^{-1/2}
  zmm("N", "C", n, w, h, g)                            # g = V E^{-1/2} V^dag
  zmm("N", "N", n, x, g, w)                            # w = X g
  for i in 0..<n*n: x[i] = w[i]
  for i in 0..<n: x[i + n*i] += complex64(1.0, 0.0)

proc ovFromX*(x: seq[Complex64], n: int): seq[Complex64] =
  ## 1 + x (x^dag x)^{-1/2} (fresh copy).
  result = x
  var h, w, g: seq[Complex64]
  var ev: seq[float]
  ovFromXInto(result, n, h, w, g, ev)
