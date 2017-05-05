import base
import layout
import fieldProxy

type
  FieldObj*[V: static[int], T] = object
    s*: alignedMem[T]
    l*: Layout[V]
    elemSize*: int
  FieldObjRef*[V: static[int], T] = ref FieldObj[V, T]
  Field*[V: static[int], T] = FieldProxy[FieldObjRef[V, T]]

template `[]`*(x: FieldObjRef, y: untyped): untyped = x[].s[y]
template `[]=`*(x: FieldObjRef, y: untyped, z: untyped): untyped =
  x[].s[y] = z

proc initFieldImpl*[V:static[int],T](x: var FieldObjRef[V,T]; l: Layout[V]) =
  x.new()
  x.l = l
  x.s.new(l.nSitesOuter)
  #fence()
  x.elemSize = sizeOf(T)

proc newFieldImpl*[V:static[int],T](x: FieldObjRef[V,T];
                                    T2: typedesc): Field[V,T2] =
  result[].initFieldImpl(x[].l)

proc newField*[V:static[int],T](l: Layout[V]; t: typedesc[T]): Field[V,T] =
  result[].initFieldImpl(l)

template indices*(x: FieldObjRef): untyped = 0..(x[].l.nSitesOuter-1)



#[
proc new*[V:static[int],T](x: var FieldObjRef[V,T]; l: Layout[V]) =
  x.l = l
  x.s.new(l.nSitesOuter)
  #fence()
  x.elemSize = sizeOf(T)

proc new*[V:static[int],T](x: var Field[V,T]; l: Layout[V]) =
  x.new()
  new(x[], l)

proc new*[V:static[int],T](x:var FieldObj[V,T]; y:Field) = x.new(y.l)

proc new*[V:static[int],T](x:var Field[V,T]; y:Field) = x.new(y.l)

proc newOneOf*(x:Field):auto =
  var r:type(x)
  r.new(x.l)
  r
]#

when isMainModule:
  #import basicOps

  var lat = [4,4,4,4]
  var lo = newLayout(lat)

  var x = lo.newField(float)
  var y = lo.newField(float)
  var z = lo.newField(float)

  x := y + z


# subsets
# reductions
