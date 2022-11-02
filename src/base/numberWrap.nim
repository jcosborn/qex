type
  AsFloat32*[T] = distinct T
  Float32* = float32 | AsFloat32
  AsFloat64*[T] = distinct T
  Float64* = float64 | AsFloat64
  AsFloat* = AsFloat32 | AsFloat64
  Float* = SomeFloat | AsFloat
  AsNumber* = AsFloat
  Number* = SomeNumber | AsNumber

template asFloat32*[T](x: T): untyped = AsFloat32[type T](x)
template `[]`*[T](x: AsFloat32[T]): untyped = (type T)(x)
template asWrapper*(x: typedesc[float32], y: typed): untyped = asFloat32(y)
template eval*(x: typedesc[AsFloat32]): typedesc = float32
template eval*(x: AsFloat32): untyped =
  mixin `:=`
  var r: float32
  r := x[]
  r
#converter toFloat32*(x: AsFloat32): float32 {.inline.} = eval(x)

template asFloat64*[T](x: T): untyped = AsFloat64[type T](x)
template `[]`*[T](x: AsFloat64[T]): untyped = (type T)(x)
template asWrapper*(x: typedesc[float64], y: typed): untyped = asFloat64(y)
template eval*(x: typedesc[AsFloat64]): typedesc = float64
template eval*(x: AsFloat64): untyped =
  mixin `:=`
  var r: float64
  r := x[]
  r

template liftUnary(fn: untyped) =
  template fn*(x: AsNumber): untyped =
    mixin fn
    fn(eval(x))

liftUnary(exp)
liftUnary(ln)
liftUnary(norm2)

template liftBinary(fn: untyped) =
  template fn*(x: AsNumber, y: SomeNumber): untyped =
    mixin fn
    fn(eval(x), y)
  template fn*(x: SomeNumber, y: AsNumber): untyped =
    mixin fn
    fn(x, eval(y))
  template fn*[X,Y:AsNumber](x: X, y: Y): untyped =
    mixin fn
    fn(eval(x), eval(y))

liftBinary(`+`)
liftBinary(`*`)

template liftBinaryInplace(fn: untyped) =
  template fn*(x: AsNumber, y: SomeNumber) =
    mixin fn
    fn(x[], y)
  template fn*(x: SomeNumber, y: AsNumber) =
    mixin fn
    fn(x, eval(y))
  template fn*[X,Y:AsNumber](x: X, y: Y) =
    mixin fn
    fn(x[], eval(y))

liftBinaryInplace(`:=`)
liftBinaryInplace(assign)
liftBinaryInplace(`+=`)
liftBinaryInplace(`*=`)
