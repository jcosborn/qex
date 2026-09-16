## Real matrix fields and real/complex scalar bridges.
## Pairings are sum_x tr(H^T dA) and sum_x Re tr(U^dag dX), without normalization.

import ../[core, scalar]
import ../support/op
import layout, physics/qcdTypes
import field/matrixFields
import types, basic_ops, cfield
export types, basic_ops

proc realNodeLike*[F](x: GfieldOf[F], n: static[int]): Grmat[n] =
  static: doAssert n == 1 or n == 8, "graph real matrix storage supports sizes 1 and 8"
  let f = x.fval.l.newShape(DRealMatrixV[n])
  Grmat[n](runtime: x.runtime, fval: f).assignStableNodeId

proc complexNodeLike(x: Grfield): Gcfield =
  let f = x.fval.l.newShape(ColorMatrixN[1,DComplexV])
  Gcfield(runtime: x.runtime, fval: f).assignStableNodeId

proc scale*[F](c: Grfield, x: GfieldOf[F]): GfieldOf[F]
proc siteRedot*[F](x, y: GfieldOf[F]): Grfield
proc siteNorm2*[F](x: GfieldOf[F]): Grfield
proc trace*[n:static[int]](x: Grmat[n]): Grfield
proc solve*[n:static[int]](a, b: Grmat[n]): Grmat[n]
proc inverse*[n:static[int]](a: Grmat[n]): Grmat[n]
proc logDet*[n:static[int]](a: Grmat[n]): Grfield
proc re*(x: Gcfield): Grfield
proc im*(x: Gcfield): Grfield
proc complex*(x: Grfield): Gcfield
proc imaginary*(x: Grfield): Gcfield
proc `/`*(x, y: Grfield): Grfield
proc exp*(x: Grfield): Grfield
proc ln*(x: Grfield): Grfield
proc sqrt*(x: Grfield): Grfield
proc sin*(x: Grfield): Grfield
proc cos*(x: Grfield): Grfield
proc expi*(x: Grfield): Gcfield
proc arg*(x: Gcfield): Grfield

proc scaleF[F](v: Gvalue) =
  threads:
    scale(GfieldOf[F](v).fval, Grfield(v.inputs[0]).fval, GfieldOf[F](v.inputs[1]).fval)

proc scaleB[F](zb, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let u = requireUpstream(zb, "real scale backward", GfieldOf[F])
  if i == 0:
    return siteRedot(GfieldOf[F](z.inputs[1]), u)
  scale(Grfield(z.inputs[0]), u)

proc scale*[F](c: Grfield, x: GfieldOf[F]): GfieldOf[F] =
  if c.fval.l != x.fval.l:
    raiseValueError("real scale requires matching field shapes")
  let op {.global.} = Gfunc(bufferMode: bmFull, inplace: @[1], forward: scaleF[F], backward: scaleB[F], name: "realScale")
  graphNode(x.fieldNodeLike, @[Gvalue(c), Gvalue(x)], op, "realScale")

proc siteRedotF[F](v: Gvalue) =
  threads:
    siteRedot(Grfield(v).fval, GfieldOf[F](v.inputs[0]).fval, GfieldOf[F](v.inputs[1]).fval)

proc siteRedotB[F](zb, z: Gvalue, i: int, input: Gvalue): Gvalue =
  scale(requireUpstream(zb, "siteRedot backward", Grfield), GfieldOf[F](z.inputs[1-i]))

proc siteRedot*[F](x, y: GfieldOf[F]): Grfield =
  x.requireSameFieldShape(y, "siteRedot")
  let op {.global.} = Gfunc(bufferMode: bmFull, forward: siteRedotF[F], backward: siteRedotB[F], name: "siteRedot")
  graphNode(x.realNodeLike(1), @[Gvalue(x), Gvalue(y)], op, "siteRedot")

proc siteNorm2F[F](v: Gvalue) =
  threads:
    siteNorm2(Grfield(v).fval, GfieldOf[F](v.inputs[0]).fval)

proc siteNorm2B[F](zb, z: Gvalue, i: int, input: Gvalue): Gvalue =
  scale(requireUpstream(zb, "siteNorm2 backward", Grfield), 2.0*GfieldOf[F](z.inputs[0]))

proc siteNorm2*[F](x: GfieldOf[F]): Grfield =
  let op {.global.} = Gfunc(bufferMode: bmFull, forward: siteNorm2F[F], backward: siteNorm2B[F], name: "siteNorm2")
  graphNode(x.realNodeLike(1), @[Gvalue(x)], op, "siteNorm2")

proc traceF[n:static[int]](v: Gvalue) =
  threads:
    siteTrace(Grfield(v).fval, Grmat[n](v.inputs[0]).fval)

proc traceB[n:static[int]](zb, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let x = Grmat[n](z.inputs[0])
  scale(requireUpstream(zb, "real trace backward", Grfield), unitField(x.runtime, x.fval))

proc trace*[n:static[int]](x: Grmat[n]): Grfield =
  let op {.global.} = Gfunc(bufferMode: bmFull, forward: traceF[n], backward: traceB[n], name: "siteTraceRmat")
  graphNode(x.realNodeLike(1), @[Gvalue(x)], op, "siteTraceRmat")

proc sum*(x: Grfield): Gscalar =
  ## Its pullback broadcasts the global upstream to every physical site.
  retr(x)

proc solveF[n:static[int]](v: Gvalue) =
  threads:
    solve(Grmat[n](v).fval, Grmat[n](v.inputs[0]).fval, Grmat[n](v.inputs[1]).fval)

proc solveB[n:static[int]](zb, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let a = Grmat[n](z.inputs[0])
  let u = requireUpstream(zb, "solve backward", Grmat[n])
  let b = solve(a.transpose, u)
  if i == 0:
    return -b*Grmat[n](z).transpose
  b

proc solve*[n:static[int]](a, b: Grmat[n]): Grmat[n] =
  ## a y = b by unpivoted LU; every leading pivot of a must be nonzero.
  a.requireSameFieldShape(b, "solve")
  let op {.global.} = Gfunc(bufferMode: bmFull, inplace: @[0, 1], forward: solveF[n], backward: solveB[n], name: "solveRmat")
  graphNode(a.fieldNodeLike, @[Gvalue(a), Gvalue(b)], op, "solveRmat")

proc inverseF[n:static[int]](v: Gvalue) =
  threads:
    inverse(Grmat[n](v).fval, Grmat[n](v.inputs[0]).fval)

proc inverseB[n:static[int]](zb, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let y = Grmat[n](z).transpose
  -y*requireUpstream(zb, "inverse backward", Grmat[n])*y

proc inverse*[n:static[int]](a: Grmat[n]): Grmat[n] =
  ## Unpivoted LU requires nonzero leading pivots at every site.
  let op {.global.} = Gfunc(bufferMode: bmFull, inplace: @[0], forward: inverseF[n], backward: inverseB[n], name: "inverseRmat")
  graphNode(a.fieldNodeLike, @[Gvalue(a)], op, "inverseRmat")

proc logDetF[n:static[int]](v: Gvalue) =
  threads:
    logDet(Grfield(v).fval, Grmat[n](v.inputs[0]).fval)

proc logDetB[n:static[int]](zb, z: Gvalue, i: int, input: Gvalue): Gvalue =
  scale(requireUpstream(zb, "logDet backward", Grfield), inverse(Grmat[n](z.inputs[0])).transpose)

proc logDet*[n:static[int]](a: Grmat[n]): Grfield =
  ## det a > 0 and nonzero leading LU pivots; individual pivots may be negative.
  let op {.global.} = Gfunc(bufferMode: bmFull, forward: logDetF[n], backward: logDetB[n], name: "logDetRmat")
  graphNode(a.realNodeLike(1), @[Gvalue(a)], op, "logDetRmat")

proc sumLogDet*[n:static[int]](a: Grmat[n]): Gscalar = sum(logDet(a))

template scalarBridge(op, back: untyped, A, B: typedesc) =
  proc op*(x {.inject.}: A): B =
    proc forward(v: Gvalue) =
      threads:
        op(B(v).fval, A(v.inputs[0]).fval)
    proc backward(zb, z: Gvalue, i: int, input: Gvalue): Gvalue =
      back(requireUpstream(zb, "scalar bridge backward", B))
    let fn {.global.} = Gfunc(bufferMode: bmFull, forward: forward, backward: backward, name: astToStr(op))
    when B is Grfield:
      graphNode(x.realNodeLike(1), @[Gvalue(x)], fn, astToStr(op))
    else:
      graphNode(x.complexNodeLike, @[Gvalue(x)], fn, astToStr(op))

scalarBridge(re, complex, Gcfield, Grfield)
scalarBridge(im, imaginary, Gcfield, Grfield)
scalarBridge(complex, re, Grfield, Gcfield)
scalarBridge(imaginary, im, Grfield, Gcfield)

proc divideF(v: Gvalue) =
  threads:
    divide(Grfield(v).fval, Grfield(v.inputs[0]).fval, Grfield(v.inputs[1]).fval)

proc divideB(zb, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let u = requireUpstream(zb, "field division backward", Grfield)
  let y = Grfield(z.inputs[1])
  if i == 0:
    return u/y
  -(u*Grfield(z))/y

let divideG = Gfunc(bufferMode: bmFull, forward: divideF, backward: divideB, name: "divideRfield")

proc `/`*(x, y: Grfield): Grfield =
  ## y must be nonzero at every site.
  graphNode(sameShapeFieldNodeLike(x, y, "real field division"), @[Gvalue(x), Gvalue(y)], divideG, "divideRfield")

proc `/`*(x: Gscalar, y: Grfield): Grfield =
  (x*unitField(y.runtime, y.fval))/y

proc `/`*(x: Grfield, y: Gscalar): Grfield =
  x/(y*unitField(x.runtime, x.fval))

proc `/`*[T:int|float](x: T, y: Grfield): Grfield = toGvalue(y.runtime, float(x))/y
proc `/`*[T:int|float](x: Grfield, y: T): Grfield = x/toGvalue(x.runtime, float(y))

template scalarFunction(op, pullback: untyped) =
  proc op*(x {.inject.}: Grfield): Grfield =
    proc forward(v: Gvalue) =
      threads:
        op(Grfield(v).fval, Grfield(v.inputs[0]).fval)
    proc backward(zb, z0: Gvalue, i: int, input: Gvalue): Gvalue =
      let x {.inject, used.} = Grfield(z0.inputs[0])
      let z {.inject, used.} = Grfield(z0)
      let u {.inject.} = requireUpstream(zb, "real scalar function backward", Grfield)
      pullback
    let fn {.global.} = Gfunc(bufferMode: bmFull, forward: forward, backward: backward, name: astToStr(op) & "Rfield")
    graphNode(x.fieldNodeLike, @[Gvalue(x)], fn, astToStr(op) & "Rfield")

scalarFunction(exp, u*z)
# ln and differentiated sqrt require x > 0 at every site.
scalarFunction(ln, u/x)
scalarFunction(sqrt, u/(2.0*z))
scalarFunction(sin, u*cos(x))
scalarFunction(cos, -u*sin(x))

proc expiF(v: Gvalue) =
  threads:
    expi(Gcfield(v).fval, Grfield(v.inputs[0]).fval)

proc expiB(zb, z: Gvalue, i: int, input: Gvalue): Gvalue =
  im(Gcfield(z).adj*requireUpstream(zb, "expi backward", Gcfield))

let expiG = Gfunc(bufferMode: bmFull, forward: expiF, backward: expiB, name: "expiField")

proc expi*(x: Grfield): Gcfield =
  graphNode(x.complexNodeLike, @[Gvalue(x)], expiG, "expiField")

proc argF(v: Gvalue) =
  threads:
    arg(Grfield(v).fval, Gcfield(v.inputs[0]).fval)

proc argB(zb, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let x = Gcfield(z.inputs[0])
  imaginary(requireUpstream(zb, "arg backward", Grfield)/siteNorm2(x))*x

let argG = Gfunc(bufferMode: bmFull, forward: argF, backward: argB, name: "argField")

proc arg*(x: Gcfield): Grfield =
  ## Principal angle in [-pi,pi]; derivatives exclude zero and the branch cut.
  graphNode(x.realNodeLike(1), @[Gvalue(x)], argG, "argField")

proc blendSubset*[F](parity: int, cand, x: GfieldOf[F]): GfieldOf[F] =
  if parity < 0 or parity > 1:
    raiseValueError("field blendSubset parity must be 0 or 1")
  let sub = x.fval.l.getSubset(if parity == 0: "even" else: "odd")
  proc forward(v: Gvalue) =
    threads:
      blendSubset(GfieldOf[F](v).fval, GfieldOf[F](v.inputs[0]).fval, GfieldOf[F](v.inputs[1]).fval, sub)
  proc backward(zb, z: Gvalue, i: int, input: Gvalue): Gvalue =
    let u = requireUpstream(zb, "field blendSubset backward", GfieldOf[F])
    let zero = GfieldOf[F](u.zeroLike)
    if i == 0:
      return blendSubset(parity, u, zero)
    blendSubset(parity, zero, u)
  graphNode(sameShapeFieldNodeLike(cand, x, "field blendSubset"), @[Gvalue(cand), Gvalue(x)],
    Gfunc(bufferMode: bmFull, forward: forward, backward: backward, name: "blendSubsetField"), "blendSubsetField")

proc maskSubset*[F](parity: int, x: GfieldOf[F]): GfieldOf[F] =
  blendSubset(parity, x, GfieldOf[F](x.zeroLike))

when DColorMatrixV is ColorMatrixN[3,DComplexV]:
  proc su3AdNeg*(x: Gfield): Grmat8
  proc su3ProjectDeriv*(x: Gfield): Grmat8
  proc su3AdNegAdj*(x: Grmat8): Gfield
  proc su3ProjectDerivAdj*(x: Grmat8): Gfield

  proc gaugeFieldLike(x: Grmat8): Gfield =
    let f = x.fval.l.newShape(DColorMatrixV)
    Gfield(runtime: x.runtime, fval: f).assignStableNodeId

  template su3Bridge(op, back: untyped, A, B: typedesc) =
    proc op*(x {.inject.}: A): B =
      proc forward(v: Gvalue) =
        threads:
          op(B(v).fval, A(v.inputs[0]).fval)
      proc backward(zb, z: Gvalue, i: int, input: Gvalue): Gvalue =
        back(requireUpstream(zb, "SU3 bridge backward", B))
      let fn {.global.} = Gfunc(bufferMode: bmFull, forward: forward, backward: backward, name: astToStr(op))
      when B is Grmat8:
        graphNode(x.realNodeLike(8), @[Gvalue(x)], fn, astToStr(op))
      else:
        graphNode(x.gaugeFieldLike, @[Gvalue(x)], fn, astToStr(op))

  su3Bridge(su3AdNeg, su3AdNegAdj, Gfield, Grmat8)
  su3Bridge(su3ProjectDeriv, su3ProjectDerivAdj, Gfield, Grmat8)
  su3Bridge(su3AdNegAdj, su3AdNeg, Grmat8, Gfield)
  su3Bridge(su3ProjectDerivAdj, su3ProjectDeriv, Grmat8, Gfield)
