import field
import simd
import base, maths

template noSimd*[V,T](x: typedesc[Field[V,T]]): untyped =
  mixin noSimd
  Field[V,noSimd(type T)]

#template noSimd[T](x: typedesc[Color[T]]): untyped =
#  mixin noSimd
#  Color[noSimd(type T)]
#template noSimd[T](x: typedesc[AsVector[T]]): untyped =
#  mixin noSimd
#  AsVector[noSimd(type T)]
#template noSimd[N,T](x: typedesc[VectorArrayObj[N,T]]): untyped =
#  mixin noSimd
#  VectorArrayObj[N,noSimd(type T)]
#template noSimd[T](x: typedesc[ComplexType[T]]): untyped =
#  mixin noSimd
#  ComplexType[noSimd(type T)]

proc checkeq*(a: Field, b: Field) =
  let nmax = 10
  var n = 0
  var ax,bx: a[0].type.noSimd
  var crd = newSeq[cint](a.l.nDim)
  for i in a.l.singleSites:
    ax := a{i}
    bx := b{i}
    let d = norm2(ax-bx)
    if d>1e-10:
      if n<nmax:
        inc n
        a.l.coord(crd,i)
        echo i, " ", crd, " ", d
        echo "  ", ax
        echo "  ", bx

type
  Norm2Stats* = tuple[s2, n2min, n2max: float]
    ## Sum, minimum and maximum of |x|^2 over sites and SIMD lanes.
  Norm2Acc* = object
    ## Norm2Stats of one thread's sites; |x|^2 >= 0 starts the maximum at zero.
    s2, mn, mx: float
    n: int

proc add*(s: var Norm2Acc, n2: auto) {.inline.} =
  ## Add |x|^2 of one site, a SIMD vector or a number.
  let lo = n2.simdMin
  let hi = n2.simdMax
  s.s2 += n2.simdSum
  if s.n == 0 or lo < s.mn: s.mn = lo
  if s.mx < hi: s.mx = hi
  inc s.n

proc total*(s: var Norm2Acc): Norm2Stats =
  ## Combine over threads and ranks; every thread calls it.
  s.s2.threadRankSum
  s.mx.threadRankMax
  if s.n == 0: s.mn = s.mx  # a thread without sites cannot lower the minimum
  var m = -s.mn
  m.threadRankMax  # min(x) = -max(-x)
  (s.s2, -m, s.mx)
