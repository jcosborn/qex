## Site-local matrix-field kernels. Call inside a threads region; layouts must match.
## Each kernel makes one pass over its output's site partition, like field assignment.
## Entries remain SIMD values. Site contractions return 1x1 matrix fields.

import layout, physics/qcdTypes
import maths/groupOps

template siteKernel(r, body: untyped) =
  for e {.inject.} in r:
    var t {.noinit, inject.}: evalType(r[e])
    body
    r[e] := t

proc transpose*[V:static[int],T](r, x: Field[V,T]) =
  static: doAssert x[0].nrows == x[0].ncols, "same-type transpose requires square site matrices"
  for e in r:
    var t {.noinit.}: evalType(x[e])
    for i in 0..<t.nrows:
      for j in 0..<t.ncols:
        t[i,j] := x[e][j,i]
    r[e] := t

proc siteTrace*[V:static[int],R,T](r: Field[V,R], x: Field[V,T]) =
  static: doAssert r[0].nrows == 1 and r[0].ncols == 1
  for e in r:
    r[e][0,0] := x[e].trace

proc siteNorm2*[V:static[int],R,T](r: Field[V,R], x: Field[V,T]) =
  static: doAssert r[0].nrows == 1 and r[0].ncols == 1
  for e in r:
    r[e][0,0] := x[e].norm2

proc siteRedot*[V:static[int],R,T](r: Field[V,R], x, y: Field[V,T]) =
  static: doAssert r[0].nrows == 1 and r[0].ncols == 1
  for e in r:
    r[e][0,0] := redot(x[e], y[e])

proc siteDot*[V:static[int],R,T](r: Field[V,R], x, y: Field[V,T]) =
  static: doAssert r[0].nrows == 1 and r[0].ncols == 1
  for e in r:
    r[e][0,0] := dot(x[e], y[e])

proc scale*[V:static[int],S,T](r: Field[V,T], c: Field[V,S], x: Field[V,T]) =
  ## r = c x for a real 1x1 field c and real or complex x.
  static: doAssert c[0].nrows == 1 and c[0].ncols == 1
  for e in r:
    for i in 0..<r[e].nrows:
      for j in 0..<r[e].ncols:
        r[e][i,j] := asReal(c[e][0,0]) * x[e][i,j]

proc scale*[V,n:static[int],T](r: Field[V,ColorMatrixN[n,ComplexType[T]]],
    c: Field[V,ColorMatrixN[1,ComplexType[T]]], x: Field[V,ColorMatrixN[n,ComplexType[T]]]) =
  ## Complex 1x1 scalar times a complex matrix field.
  for e in r:
    for i in 0..<n:
      for j in 0..<n:
        r[e][i,j] := c[e][0,0] * x[e][i,j]

proc solve*[V:static[int],T](r, a, b: Field[V,T]) =
  ## r = a^-1 b, using unpivoted LU at each site.
  siteKernel(r):
    solve(t[], a[e][], b[e][])

proc inverse*[V:static[int],T](r, a: Field[V,T]) =
  siteKernel(r):
    inverseNoPivot(t[], a[e][])

proc logDet*[V:static[int],R,T](r: Field[V,R], a: Field[V,T]) =
  ## det a > 0, with nonzero leading LU pivots; negative pivots are allowed.
  static: doAssert r[0].nrows == 1 and r[0].ncols == 1
  for e in r:
    r[e][0,0] := logDet(a[e][])

proc sum*[V:static[int],T:SomeNumber|Simd](x: Field[V,ColorMatrixN[1,T]]): auto =
  ## Sum over all physical sites and ranks, unnormalized.
  ## A thread reduction: every thread calls it inside threads.
  x.trace

template realPart(op: untyped) =
  proc op*[V,n:static[int],T](r: Field[V,ColorMatrixN[n,T]], x: Field[V,ColorMatrixN[n,ComplexType[T]]]) =
    for e in r:
      for i in 0..<n:
        for j in 0..<n:
          r[e][i,j] := x[e][i,j].op

realPart(re)
realPart(im)

template complexPart(op, real, imag: untyped) =
  proc op*[V,n:static[int],T](r: Field[V,ColorMatrixN[n,ComplexType[T]]], x: Field[V,ColorMatrixN[n,T]]) =
    for e in r:
      for i in 0..<n:
        for j in 0..<n:
          let a {.inject.} = x[e][i,j]
          r[e][i,j].re := real
          r[e][i,j].im := imag

complexPart(complex, a, 0)
complexPart(imaginary, 0, a)

proc toMatrix*[V:static[int],R,T](r: Field[V,ColorMatrixN[1,R]], x: Field[V,T]) =
  ## Scalar sites into 1x1 matrix sites; the precisions may differ.
  for e in r:
    r[e][0,0] := x[e]

proc toScalar*[V:static[int],R,T](r: Field[V,R], x: Field[V,ColorMatrixN[1,T]]) =
  ## 1x1 matrix sites into scalar sites; the precisions may differ.
  for e in r:
    r[e] := x[e][0,0]

proc blendSubset*[V:static[int],T](r, cand, x: Field[V,T], sub: Subset) =
  ## r = cand on sub's contiguous outer-site interval, x elsewhere.
  ## Whole SIMD sites are selected; either input may alias r.
  for e in r:
    if e >= sub.lowOuter and e < sub.highOuter: r[e] := cand[e]
    else: r[e] := x[e]

proc maskSubset*[V:static[int],T](r, x: Field[V,T], sub: Subset) =
  ## r = x on sub's contiguous outer-site interval, 0 elsewhere.
  ## Whole SIMD sites are selected; x may alias r.
  for e in r:
    if e >= sub.lowOuter and e < sub.highOuter: r[e] := x[e]
    else: r[e] := 0

template scalarFunction(op: untyped) =
  proc op*[V:static[int],T](r, x: Field[V,ColorMatrixN[1,T]]) =
    mixin op
    for e in r:
      r[e][0,0] := op(x[e][0,0])

scalarFunction(exp)
scalarFunction(ln)
scalarFunction(sqrt)
scalarFunction(sin)
scalarFunction(cos)

proc divide*[V:static[int],T](r, x, y: Field[V,ColorMatrixN[1,T]]) =
  for e in r:
    r[e][0,0] := x[e][0,0] / y[e][0,0]

proc expi*[V:static[int],T](r: Field[V,ColorMatrixN[1,ComplexType[T]]],
                           x: Field[V,ColorMatrixN[1,T]]) =
  for e in r:
    r[e][0,0].re := cos(x[e][0,0])
    r[e][0,0].im := sin(x[e][0,0])

proc expi*[V:static[int],T](r: Field[V,ColorMatrixN[1,ComplexType[T]]], x: Field[V,T]) =
  ## Unit phases exp(i x) of a real scalar field.
  for e in r:
    r[e][0,0].re := cos(x[e])
    r[e][0,0].im := sin(x[e])

proc arg*[V:static[int],T](r: Field[V,ColorMatrixN[1,T]],
                          x: Field[V,ColorMatrixN[1,ComplexType[T]]]) =
  ## Principal angle in [-pi,pi]; derivatives exclude zero and the branch cut.
  for e in r:
    r[e][0,0] := atan2(x[e][0,0].im, x[e][0,0].re)

# Layout-independent overloads keep the element kernels available to raw users.
proc su3AdNeg*[V:static[int],R,T](r: Field[V,R], x: Field[V,T]) =
  siteKernel(r):
    su3AdNeg(t[], x[e][])

proc su3ProjectDeriv*[V:static[int],R,T](r: Field[V,R], x: Field[V,T]) =
  siteKernel(r):
    su3ProjectDeriv(t[], x[e][])

proc su3AdNegAdj*[V:static[int],R,T](r: Field[V,R], x: Field[V,T]) =
  siteKernel(r):
    su3AdNegAdj(t[], x[e][])

proc su3ProjectDerivAdj*[V:static[int],R,T](r: Field[V,R], x: Field[V,T]) =
  siteKernel(r):
    su3ProjectDerivAdj(t[], x[e][])

# Fused field kernels. Each runs inside `threads` on the calling thread's site
# partition, on one field or on every direction of a gauge. The destination
# may alias an operand.

proc contractProjectTAH*[F:Field](z: F; x, y: F) =
  ## z = projectTAH(x y†) site by site.
  for e in z:
    let s = x[e] * y[e].adj
    z[e].projectTAH s

proc contractProjectTAH*[F:Field](z: F; x, y: F; sub: Subset) =
  ## The same on sub's outer sites only.
  for e in sub:
    let s = x[e] * y[e].adj
    z[e].projectTAH s

proc contractProjectTAH*[F:Field](z: openArray[F]; x, y: openArray[F]) =
  for mu in 0..<z.len: contractProjectTAH(z[mu], x[mu], y[mu])

template axexpSites(z, ex, a, x, y, sites: typed; kern: untyped; writeE, writeZ: static bool) =
  for e in sites:
    var t {.noinit.}: evalType(x[e])
    t[] := kern(a * x[e][])
    when writeE: ex[e] := t
    when writeZ: z[e] := t * y[e]

template axexpKernels(axexp, axexpmuly, kern, what: untyped) =
  proc axexp*[F:Field](z: F; a: float; x: F) =
    ## z = `what`(a x) site by site.
    axexpSites(z, z, a, x, x, z, kern, true, false)
  proc axexp*[F:Field](z: F; a: float; x: F; sub: Subset) =
    ## z = `what`(a x) on sub's outer sites only; the complement is untouched.
    axexpSites(z, z, a, x, x, sub, kern, true, false)
  proc axexp*[F:Field](z: openArray[F]; a: float; x: openArray[F]) =
    ## z[mu] = `what`(a x[mu]) site by site.
    for mu in 0..<z.len: axexp(z[mu], a, x[mu])
  proc axexpmuly*[F:Field](z: F; a: float; x, y: F; expax: F = nil) =
    ## z = `what`(a x) y site by site; expax, when given, receives `what`(a x)
    ## and must be distinct from z and y.
    if expax == nil: axexpSites(z, expax, a, x, y, z, kern, false, true)
    else: axexpSites(z, expax, a, x, y, z, kern, true, true)
  proc axexpmuly*[F:Field](z: F; a: float; x, y: F; sub: Subset; expax: F = nil) =
    ## The same on sub's outer sites only; the complement is untouched.
    if expax == nil: axexpSites(z, expax, a, x, y, sub, kern, false, true)
    else: axexpSites(z, expax, a, x, y, sub, kern, true, true)
  proc axexpmuly*[F:Field](z: openArray[F]; a: float; x, y: openArray[F]) =
    ## z[mu] = `what`(a x[mu]) y[mu] site by site.
    for mu in 0..<z.len: axexpmuly(z[mu], a, x[mu], y[mu])
  proc axexpmuly*[F:Field](z: openArray[F]; a: float; x, y, expax: openArray[F]) =
    ## The same with expax[mu] receiving `what`(a x[mu]).
    for mu in 0..<z.len: axexpmuly(z[mu], a, x[mu], y[mu], expax[mu])

# axexp/axexpmuly: any x, with the general `exp` of maths/matrixFunctions.
axexpKernels(axexp, axexpmuly, exp, "exp")
# axexpah/axexpahmuly: x anti-Hermitian (traceless for SU(N)), with expAH.
axexpKernels(axexpah, axexpahmuly, expAH, "expAH")
