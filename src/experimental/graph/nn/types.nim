import qex
import ../../../nn
import ../core
import ../field/types
import std/sequtils
export types

type
  Greal*[T: SomeFloat] = GfieldOf[seq[RealField[T]]]
  Gmask* = GfieldOf[RealField[float32]]
  Garray*[T: SomeFloat] = ref object of Gvalue
    ## Spatially replicated parameters; axes describe channels and taps.
    data*: seq[T]
    shape*: seq[int]

template realValues(T: typedesc) =
  proc realNodeLike*(x: Greal[T]): Greal[T] = x.fieldNodeLike
  proc gradSeeded*(dep, x, seed: Greal[T]): Greal[T] =
    Greal[T](core.gradSeeded(Gvalue(dep), Gvalue(x), Gvalue(seed)))

realValues(float32)
realValues(float64)

proc toGvalue*[F](rt: GraphRuntime, fs: seq[F]): auto =
  type T = numberType(F)
  static: doAssert F is RealField[T]
  toGfield(rt, fs)

proc toGarray*[T: SomeFloat](rt: GraphRuntime, data: openArray[T], shape: openArray[int]): Garray[T] =
  var n = 1
  for d in shape:
    if d <= 0: raiseValueError("parameter array extents must be positive")
    n *= d
  if shape.len == 0 or data.len != n:
    raiseValueError("parameter array data length differs from shape")
  result = Garray[T](runtime: rt, data: @data, shape: @shape).assignStableNodeId
  result.updated

proc update*[T: SomeFloat](x: Garray[T], data: openArray[T]) =
  if data.len != x.data.len: raiseValueError("parameter update changes shape")
  for i in 0..<data.len: x.data[i] = data[i]
  x.updated

proc arrayNodeLike*[T: SomeFloat](x: Garray[T]): Garray[T] =
  Garray[T](runtime: x.runtime, data: newSeq[T](x.data.len), shape: x.shape.mapIt(it)).assignStableNodeId

template arrayMethods(T: typedesc) =
  method newOneOf*(x: Garray[T]): Gvalue = x.arrayNodeLike
  method valueLike*(x: Garray[T]): Gvalue = x.arrayNodeLike
  method zeroLike*(x: Garray[T]): Gvalue =
    result = x.arrayNodeLike
    result.staticZeroLeaf = true
  method oneLike*(x: Garray[T]): Gvalue =
    let z = x.arrayNodeLike
    for i in 0..<z.data.len: z.data[i] = 1
    z.updated
    z
  method isZero*(x: Garray[T]): bool = x.staticZeroLeaf
  method copyCompatible*(x: Garray[T], y: Gvalue): bool =
    y of Garray[T] and x.shape == Garray[T](y).shape
  method valCopy*(z: Garray[T], x: Gvalue) =
    if not z.copyCompatible(x): raiseValueError("parameter copy shape or precision differs")
    let a = Garray[T](x)
    for i in 0..<z.data.len: z.data[i] = a.data[i]

fieldMethods(Greal[float32], "RealFields32")
fieldMethods(Greal[float64], "RealFields64")
fieldMethods(Gmask, "RealField32")
arrayMethods(float32)
arrayMethods(float64)
