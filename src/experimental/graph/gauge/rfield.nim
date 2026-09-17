## Precision/channel adapters for the gauge layer's real matrix fields.
import ../[core, scalar]
import ../nn as neural
import ../support/op
import types, matrix
import qex
import field/matrixFields
import ../../../nn as numeric

proc channels[T: SomeFloat](x: Grfield; precision: typedesc[T]): Greal[T]
proc matrixField[T: SomeFloat](x: Greal[T]): Grfield

proc channelsF[T: SomeFloat](v: Gvalue) =
  let x = Grfield(v.inputs[0])
  let z = Greal[T](v)
  threads:
    toScalar(z.fval[0], x.fval)

proc channelsB[T: SomeFloat](zb,z: Gvalue; i: int; input: Gvalue): Gvalue =
  matrixField[T](requireUpstream(zb,"real channels backward",Greal[T]))

proc channels[T: SomeFloat](x: Grfield; precision: typedesc[T]): Greal[T] =
  let f = x.fval.l.newShape(RealField[T].T)
  let op {.global.} = Gfunc(bufferMode:bmFull,forward:channelsF[T],backward:channelsB[T],name:"realChannels")
  graphNode(Greal[T](runtime:x.runtime,fval: @[f]),@[Gvalue(x)],op,"realChannels")

proc matrixF[T: SomeFloat](v: Gvalue) =
  let x = Greal[T](v.inputs[0])
  let z = Grfield(v)
  threads:
    toMatrix(z.fval, x.fval[0])

proc matrixB[T: SomeFloat](zb,z: Gvalue; i: int; input: Gvalue): Gvalue =
  channels(requireUpstream(zb,"real matrix backward",Grfield),T)

proc matrixField[T: SomeFloat](x: Greal[T]): Grfield =
  if x.fval.len != 1:
    raiseValueError("real matrix adapter requires one real channel")
  let f = x.fval[0].l.newShape(DRealMatrixV[1])
  let op {.global.} = Gfunc(bufferMode:bmFull,forward:matrixF[T],backward:matrixB[T],name:"realMatrix")
  graphNode(Grfield(runtime:x.runtime,fval:f),@[Gvalue(x)],op,"realMatrix")

proc real*[T: SomeFloat](x: Gcfield; precision: typedesc[T]): Greal[T] =
  channels(matrix.re(x),T)

proc imag*[T: SomeFloat](x: Gcfield; precision: typedesc[T]): Greal[T] =
  channels(matrix.im(x),T)

proc complexImpl[T: SomeFloat](re,im: Greal[T]): Gcfield =
  matrix.complex(matrixField[T](re)) + matrix.imaginary(matrixField[T](im))

proc scaleImpl[T: SomeFloat](c: Greal[T]; x: Gfield): Gfield =
  matrix.scale(matrixField[T](c),x)

template bridgeOps(T: typedesc) =
  proc complex*(re,im: Greal[T]): Gcfield = complexImpl[T](re,im)
  proc scale*(c: Greal[T]; x: Gfield): Gfield = scaleImpl[T](c,x)

bridgeOps(float32)
bridgeOps(float64)
