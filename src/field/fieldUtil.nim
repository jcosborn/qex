import field
import simd

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
