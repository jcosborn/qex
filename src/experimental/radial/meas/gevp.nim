## Generalized eigenvalue problem for correlator matrices  (WP-I).
##
##   C(t) v = lambda(t, t0) C(t0) v
##
## solved with QEX's committed `eigens/linalgFuncs.zeigsgv` (zhegv underneath,
## Hermitian generalized, with its automatic diagonal-regularization retry).
## Effective dimensions: Delta_n(t) = -ln(lambda_n(t+at)/lambda_n(t))/at.
##
## Normative reference: doc/07-observables.md section 4.3,
## doc/04-interfaces.md section 14.
##
## Correlator layout: `c[t][i][j]`, real, one matrix per time separation.
## Every matrix is symmetrized as (C + C^T)/2 before use -- the measured
## estimators are only symmetric in expectation.  `gevpCheck` reports the
## spectrum of C(t0): the GEVP is only meaningful where C(t0) is positive
## SEMIdefinite, and the condition number bounds the accuracy of the small
## eigenvalues.
##
## Rank truncation (`cut`).  A correlator matrix is rank deficient whenever
## the channel holds fewer states than operators -- the NORMAL situation for
## a small basis on a coarse lattice (measured: the L = 1 l = 1 gauge channel
## has exactly ONE state, so ANY C(t0) there is rank 1).  Directions of C(t0)
## below cut*evmax are projected out and the remaining problem is whitened and
## solved as an ordinary Hermitian eigenproblem; only when every direction
## passes the cut is the committed zeigsgv (zhegv + regularized retry) used
## directly.  Without the truncation the near-null directions produce garbage
## generalized eigenvalues that can overtake the ground state.  For noisy
## (Monte Carlo) matrices set `cut` at the relative noise level of C(t0).

import std/[math, complex]
import eigens/linalgFuncs

proc toHerm(m: seq[seq[float]]): seq[Complex64] =
  ## Column-major complex copy of the symmetrized real matrix.
  let n = m.len
  result = newSeq[Complex64](n*n)
  for i in 0..<n:
    for j in 0..<n:
      result[i + n*j] = complex64(0.5*(m[i][j] + m[j][i]), 0.0)

proc gevp*(c: seq[seq[seq[float]]], t0, t: int, cut = 1e-10): seq[float] =
  ## The generalized eigenvalues lambda_n(t, t0) of C(t) v = lambda C(t0) v,
  ## DESCENDING (n = 0 is the ground state), on the subspace where C(t0) has
  ## eigenvalues above cut*evmax.  Result length = that rank.
  let n = c[t].len
  doAssert c[t0].len == n
  var
    a = toHerm(c[t])
    b = toHerm(c[t0])
    ev = newSeq[float](n)
  zeigs(cast[ptr float64](addr b[0]), addr ev[0], n)   # b <- eigenvectors, ascending
  doAssert ev[n-1] > 0.0, "gevp: C(t0) has no positive direction"
  var m = 0
  for k in 0..<n:
    if ev[k] > cut*ev[n-1]: inc m
  if m == n:
    var
      a2 = toHerm(c[t])
      b2 = toHerm(c[t0])
      e = newSeq[float](n)
    zeigsgv(cast[ptr float64](addr a2[0]), cast[ptr float64](addr b2[0]),
            addr e[0], n)
    result = newSeq[float](n)
    for k in 0..<n: result[k] = e[n - 1 - k]
  else:
    # whiten: At = W^dag A W, W = V_k/sqrt(ev_k) for the kept (top) directions
    var at = newSeq[Complex64](m*m)
    var w = newSeq[Complex64](n)
    for kj in 0..<m:
      let j = n - m + kj
      for p in 0..<n:
        var s = complex64(0.0, 0.0)
        for q in 0..<n: s += a[p + n*q]*b[q + n*j]
        w[p] = s
      for ki in 0..<m:
        let i = n - m + ki
        var s = complex64(0.0, 0.0)
        for p in 0..<n: s += conjugate(b[p + n*i])*w[p]
        at[ki + m*kj] = s/sqrt(ev[i]*ev[j])
    var e = newSeq[float](m)
    zeigs(cast[ptr float64](addr at[0]), addr e[0], m)
    result = newSeq[float](m)
    for k in 0..<m: result[k] = e[m - 1 - k]

proc gevpDims*(c: seq[seq[seq[float]]], t0: int, at: float,
               cut = 1e-10): seq[seq[float]] =
  ## Delta_n(t) = -ln(lambda_n(t+at)/lambda_n(t))/at from successive GEVPs at
  ## fixed t0.  result[t][n] pairs times t and t+1 (so len = c.len - 1); entries
  ## are NaN wherever a lambda is non-positive (noise), left to propagate.
  result = newSeq[seq[float]](c.len - 1)
  var lam = newSeq[seq[float]](c.len)
  for t in 0..<c.len:
    lam[t] = gevp(c, t0, t, cut)
  let n = lam[0].len
  for t in 0..<c.len - 1:
    result[t] = newSeq[float](n)
    for k in 0..<n:
      result[t][k] = -ln(lam[t+1][k]/lam[t][k])/at

proc gevpCheck*(c: seq[seq[seq[float]]], t0: int):
    tuple[evmin, evmax, cond, asym: float] =
  ## Health of the reference matrix: extreme eigenvalues of the symmetrized
  ## C(t0), their ratio (reported as `cond`; negative `evmin` makes it
  ## meaningless and flags an indefinite C(t0)), and the largest antisymmetric
  ## part |C - C^T| of the raw input relative to the largest entry.
  let n = c[t0].len
  var asym = 0.0
  var scale = 0.0
  for i in 0..<n:
    for j in 0..<n:
      asym = max(asym, abs(c[t0][i][j] - c[t0][j][i]))
      scale = max(scale, abs(c[t0][i][j]))
  var a = toHerm(c[t0])
  var e = newSeq[float](n)
  zeigs(cast[ptr float64](addr a[0]), addr e[0], n)   # ascending
  result.evmin = e[0]
  result.evmax = e[n-1]
  result.cond = e[n-1]/e[0]
  result.asym = if scale > 0.0: asym/scale else: 0.0
