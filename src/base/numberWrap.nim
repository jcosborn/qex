type
  AsFloat32*[T] = distinct T
  Float32* = float32 | AsFloat32
  AsFloat64*[T] = distinct T
  Float64* = float64 | AsFloat64
  AsFloat* = AsFloat32 | AsFloat64
  Float* = SomeFloat | AsFloat
  AsNumber* = AsFloat
  Number* = SomeNumber | AsNumber

template asFloat32*[T](x: T): auto = AsFloat32[type T](x)
template `[]`*[T](x: AsFloat32[T]): auto = (type T)(x)
template asWrapper*(x: typedesc[float32], y: typed): auto = asFloat32(y)
template eval*(x: typedesc[AsFloat32]): typedesc = float32
template eval*(x: AsFloat32): auto =
  mixin `:=`
  var r {.noInit.}: float32
  r := x[]
  r

template asFloat64*[T](x: T): auto = AsFloat64[type T](x)
template `[]`*[T](x: AsFloat64[T]): auto = (type T)(x)
template asWrapper*(x: typedesc[float64], y: typed): auto = asFloat64(y)
template eval*(x: typedesc[AsFloat64]): typedesc = float64
template eval*(x: AsFloat64): auto =
  mixin `:=`
  var r {.noInit.}: float64
  r := x[]
  r

template liftUnary(fn: untyped) =
  template fn*(x: AsNumber): auto =
    mixin fn
    fn(eval(x))
  template fn*(x: typedesc[AsNumber]): typedesc =
    mixin fn
    fn(eval(x))

liftUnary(numNumbers)
liftUnary(exp)
liftUnary(ln)
liftUnary(norm2)

template liftBinary(fn: untyped) =
  template fn*(x: AsNumber, y: SomeNumber): auto =
    mixin fn
    fn(eval(x), y)
  template fn*(x: SomeNumber, y: AsNumber): auto =
    mixin fn
    fn(x, eval(y))
  template fn*[X,Y:AsNumber](x: X, y: Y): auto =
    mixin fn
    fn(eval(x), eval(y))
  template fn*(x: typedesc[AsNumber], y: typedesc[SomeNumber]): typedesc =
    mixin fn
    fn(eval(x), y)
  template fn*(x: typedesc[SomeNumber], y: typedesc[AsNumber]): typedesc =
    mixin fn
    fn(x, eval(y))
  template fn*[X,Y:AsNumber](x: typedesc[X], y: typedesc[Y]): typedesc =
    mixin fn
    fn(eval(x), eval(y))

liftBinary(`+`)
liftBinary(`-`)
liftBinary(`*`)
liftBinary(`/`)

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
