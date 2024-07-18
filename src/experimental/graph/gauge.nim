#[

- Allocate the field when constructing the graph.
- Call updated after changing the field in the graph after construction.

]#

import core, scalar, multi
import layout, ../../gauge, physics/qcdTypes

type Gauge = seq[DLatticeColorMatrixV]

type Ggauge* {.final.} = ref object of Gvalue
  gval: Gauge

proc getgauge*(x: Gvalue): Gauge = Ggauge(x).gval

proc update*(x: Gvalue, g: Gauge) =
  let u = x.getgauge
  threads:
    for mu in 0..<u.len:
      u[mu] := g[mu]
  x.updated

proc toGvalue*(x: Gauge): Ggauge =
  # proc instead of converter to avoid converting seq
  result = Ggauge(gval: x)
  result.updated

method newOneOf*(x: Ggauge): Gvalue =
  let g = x.gval.newOneOf
  threads:
    for f in g:
      f := 0.0
  Ggauge(gval: g)
method valCopy*(z: Ggauge, x: Ggauge) =
  let u = z.gval
  let v = x.gval
  threads:
    for mu in 0..<u.len:
      u[mu] := v[mu]

method `$`*(x: Ggauge): string = "Gauge"

#
# basic ops
#

method retr*(x: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("retr(" & $x & ")")
method adj*(x: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("adj(" & $x & ")")
method norm2*(x: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("norm2(" & $x & ")")
method redot*(x: Gvalue, y: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("redot(" & $x & "," & $y & ")")
method expDeriv*(b: Gvalue, x: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("expDeriv(" & $b & "," & $x & ")")
method projTAH*(x: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("projTAH(" & $x & ")")

method redot*(x: Gscalar, y: Gscalar): Gvalue = x*y

proc retrgb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  let g = z.inputs[0].getgauge.newOneOf
  threads:
    for f in g:
      f := 1.0
  let one = toGvalue g
  case i
  of 0:
    if zb == nil:
      return one
    else:
      return zb*one
  else:
    raiseValueError("i must be 0, got: " & $i)

proc retrgf(v: Gvalue) =
  let x = Ggauge(v.inputs[0])
  let z = Gscalar(v)
  threads:
    var t = 0.0
    for mu in 0..<x.gval.len:
      t += x.gval[mu].trace.re
    threadMaster: z.getfloat = t

let retrg = newGfunc(forward = retrgf, backward = retrgb, name = "retrg")

method retr*(x: Ggauge): Gvalue = Gscalar(inputs: @[Gvalue(x)], gfunc: retrg)

proc adjgb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  if zb == nil:
    raiseValueError("no backward node")
  case i
  of 0:
    return zb.adj
  else:
    raiseValueError("i must be 0, got: " & $i)

proc adjgf(v: Gvalue) =
  let x = Ggauge(v.inputs[0])
  let z = Ggauge(v)
  threads:
    for mu in 0..<z.gval.len:
      z.gval[mu] := x.gval[mu].adj

let adjg = newGfunc(forward = adjgf, backward = adjgb, name = "adjg")

method adj*(x: Ggauge): Gvalue = Ggauge(gval: x.gval.newOneOf, inputs: @[Gvalue(x)], gfunc: adjg)

proc norm2gb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  case i
  of 0:
    if zb == nil:
      return 2.0 * z.inputs[0]
    else:
      return (2.0 * zb) * z.inputs[0]
  else:
    raiseValueError("i must be 0, got: " & $i)

proc norm2gf(v: Gvalue) =
  let x = Ggauge(v.inputs[0])
  let z = Gscalar(v)
  threads:
    var t = 0.0
    for mu in 0..<x.gval.len:
      t += x.gval[mu].norm2
    threadMaster: z.getfloat = t

let norm2g = newGfunc(forward = norm2gf, backward = norm2gb, name = "norm2g")

method norm2*(x: Ggauge): Gvalue = Gscalar(inputs: @[Gvalue(x)], gfunc: norm2g)

proc neggf(v: Gvalue) =
  let x = Ggauge(v.inputs[0])
  let z = Ggauge(v)
  threads:
    for mu in 0..<z.gval.len:
      z.gval[mu] := -x.gval[mu]

proc neggb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  if zb == nil:
    raiseValueError("no backward node")
  case i
  of 0:
    return -zb
  else:
    raiseValueError("i must be 0, got: " & $i)

let negg = newGfunc(forward = neggf, backward = neggb, name = "-g")

method `-`*(x: Ggauge): Gvalue = Ggauge(gval: x.gval.newOneOf, inputs: @[Gvalue(x)], gfunc: negg)

proc addsgb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  if zb == nil:
    raiseValueError("no backward node")
  case i
  of 0:
    return zb.retr
  of 1:
    return zb
  else:
    raiseValueError("i must be 0 or 1, got: " & $i)

proc addsgf(v: Gvalue) =
  let x = Gscalar(v.inputs[0])
  let y = Ggauge(v.inputs[1])
  let z = Ggauge(v)
  threads:
    for mu in 0..<z.gval.len:
      z.gval[mu] := x.getfloat + y.gval[mu]

let addsg = newGfunc(forward = addsgf, backward = addsgb, name = "s+g")

method `+`*(x: Gscalar, y: Ggauge): Gvalue = Ggauge(gval: y.gval.newOneOf, inputs: @[Gvalue(x), y], gfunc: addsg)

proc addggb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  if zb == nil:
    raiseValueError("no backward node")
  else:
    return zb

proc addggf(v: Gvalue) =
  let x = Ggauge(v.inputs[0])
  let y = Ggauge(v.inputs[1])
  let z = Ggauge(v)
  threads:
    for mu in 0..<z.gval.len:
      z.gval[mu] := x.gval[mu] + y.gval[mu]

let addgg = newGfunc(forward = addggf, backward = addggb, name = "g+g")

method `+`*(x: Ggauge, y: Ggauge): Gvalue = Ggauge(gval: x.gval.newOneOf, inputs: @[Gvalue(x), y], gfunc: addgg)

proc mulsgb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  if zb == nil:
    raiseValueError("no backward node")
  case i
  of 0:
    return redot(zb, z.inputs[1])
  of 1:
    return z.inputs[0]*zb
  else:
    raiseValueError("i must be 0 or 1, got: " & $i)

proc mulsgf(v: Gvalue) =
  let x = Gscalar(v.inputs[0])
  let y = Ggauge(v.inputs[1])
  let z = Ggauge(v)
  threads:
    for mu in 0..<z.gval.len:
      z.gval[mu] := x.getfloat * y.gval[mu]

let mulsg = newGfunc(forward = mulsgf, backward = mulsgb, name = "s*g")

method `*`*(x: Gscalar, y: Ggauge): Gvalue = Ggauge(gval: y.gval.newOneOf, inputs: @[Gvalue(x), y], gfunc: mulsg)

proc mulggb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  if zb == nil:
    raiseValueError("no backward node")
  case i
  of 0:
    return zb*z.inputs[1].adj
  of 1:
    return z.inputs[0].adj*zb
  else:
    raiseValueError("i must be 0 or 1, got: " & $i)

proc mulggf(v: Gvalue) =
  let x = Ggauge(v.inputs[0])
  let y = Ggauge(v.inputs[1])
  let z = Ggauge(v)
  threads:
    for mu in 0..<z.gval.len:
      z.gval[mu] := x.gval[mu] * y.gval[mu]

let mulgg = newGfunc(forward = mulggf, backward = mulggb, name = "g*g")

method `*`*(x: Ggauge, y: Ggauge): Gvalue = Ggauge(gval: x.gval.newOneOf, inputs: @[Gvalue(x), y], gfunc: mulgg)

proc redotggb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  case i
  of 0:
    if zb == nil:
      return z.inputs[1]
    else:
      return zb*z.inputs[1]
  of 1:
    if zb == nil:
      return z.inputs[0]
    else:
      return zb*z.inputs[0]
  else:
    raiseValueError("i must be 0 or 1, got: " & $i)

proc redotggf(v: Gvalue) =
  let x = Ggauge(v.inputs[0])
  let y = Ggauge(v.inputs[1])
  let z = Gscalar(v)
  threads:
    var t = 0.0
    for mu in 0..<x.gval.len:
      t += redot(x.gval[mu], y.gval[mu])
    threadMaster: z.getfloat = t

let redotgg = newGfunc(forward = redotggf, backward = redotggb, name = "redotgg")

method redot*(x: Ggauge, y: Ggauge): Gvalue = Gscalar(inputs: @[Gvalue(x), y], gfunc: redotgg)

proc subgsb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  if zb == nil:
    raiseValueError("no backward node")
  case i
  of 0:
    return zb
  of 1:
    return -zb.retr
  else:
    raiseValueError("i must be 0 or 1, got: " & $i)

proc subgsf(v: Gvalue) =
  let x = Ggauge(v.inputs[0])
  let y = Gscalar(v.inputs[1])
  let z = Ggauge(v)
  threads:
    for mu in 0..<z.gval.len:
      z.gval[mu] := x.gval[mu] - y.getfloat

let subgs = newGfunc(forward = subgsf, backward = subgsb, name = "g-s")

method `-`*(x: Ggauge, y: Gscalar): Gvalue = Ggauge(gval: x.gval.newOneOf, inputs: @[Gvalue(x), y], gfunc: subgs)

proc subggb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  if zb == nil:
    raiseValueError("no backward node")
  case i
  of 0:
    return zb
  of 1:
    return -zb
  else:
    raiseValueError("i must be 0 or 1, got: " & $i)

proc subggf(v: Gvalue) =
  let x = Ggauge(v.inputs[0])
  let y = Ggauge(v.inputs[1])
  let z = Ggauge(v)
  threads:
    for mu in 0..<z.gval.len:
      z.gval[mu] := x.gval[mu] - y.gval[mu]

let subgg = newGfunc(forward = subggf, backward = subggb, name = "g-g")

method `-`*(x: Ggauge, y: Ggauge): Gvalue = Ggauge(gval: x.gval.newOneOf, inputs: @[Gvalue(x), y], gfunc: subgg)

proc expgb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  if zb == nil:
    raiseValueError("no backward node")
  case i
  of 0:
    return expDeriv(zb, z.inputs[0])
  else:
    raiseValueError("i must be 0, got: " & $i)

proc expgf(v: Gvalue) =
  let x = Ggauge(v.inputs[0])
  let z = Ggauge(v)
  threads:
    for mu in 0..<z.gval.len:
      for e in z.gval[mu]:
        z.gval[mu][e] := exp(x.gval[mu][e])

let expg = newGfunc(forward = expgf, backward = expgb, name = "expg")

method exp*(x: Ggauge): Gvalue = Ggauge(gval: x.gval.newOneOf, inputs: @[Gvalue(x)], gfunc: expg)

proc expDerivgb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  if zb == nil:
    raiseValueError("no backward node")
  case i
  of 0:
    raiseValueError("unimplemented")
  else:
    raiseValueError("i must be 0, got: " & $i)

proc expDerivgf(v: Gvalue) =
  let x = Ggauge(v.inputs[0])
  let y = Ggauge(v.inputs[1])
  let z = Ggauge(v)
  threads:
    for mu in 0..<z.gval.len:
      for e in z.gval[mu]:
        z.gval[mu][e] := expDeriv(y.gval[mu][e], x.gval[mu][e])

let expDerivg = newGfunc(forward = expDerivgf, backward = expDerivgb, name = "expDerivg")

method expDeriv*(b: Ggauge, x: Ggauge): Gvalue = Ggauge(gval: x.gval.newOneOf, inputs: @[Gvalue(b), x], gfunc: expDerivg)

proc projTAHb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  if zb == nil:
    raiseValueError("no backward node")
  case i
  of 0:
    return zb.projTAH
  else:
    raiseValueError("i must be 0, got: " & $i)

proc projTAHf(v: Gvalue) =
  let x = Ggauge(v.inputs[0])
  let z = Ggauge(v)
  threads:
    for mu in 0..<z.gval.len:
      for e in z.gval[mu]:
        z.gval[mu][e].projectTAH(x.gval[mu][e])

let projTAHg = newGfunc(forward = projTAHf, backward = projTAHb, name = "projTAH")

method projTAH*(x: Ggauge): Gvalue = Ggauge(gval: x.gval.newOneOf, inputs: @[Gvalue(x)], gfunc: projTAHg)

#
# Fused ops
#

#[ TODO:
- reuse common terms in the backward
  - save & reuse local vars attached to a node
- fuse the backward
]#

method contractProjTAH*(x: Gvalue, y: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("contractProjTAH(" & $x & "," & $y & ")")
method axexp*(a: Gvalue, x: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("axexp(" & $a & "," & $x & ")")
method axexpmuly*(a: Gvalue, x: Gvalue, y: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("axexpmuly(" & $a & "," & $x & "," & $y & ")")

proc contractProjTAHb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  if zb == nil:
    raiseValueError("no backward node")
  let zbTAH = zb.projTAH
  case i
  of 0:
    return zbTAH * z.inputs[1]
  of 1:
    return zbTAH.adj * z.inputs[0]
  else:
    raiseValueError("i must be 0 or 1, got: " & $i)

proc contractProjTAHf(v: Gvalue) =
  let x = Ggauge(v.inputs[0])
  let y = Ggauge(v.inputs[1])
  let z = Ggauge(v)
  threads:
    for mu in 0..<z.gval.len:
      for e in z.gval[mu]:
        let s = x.gval[mu][e]*y.gval[mu][e].adj
        z.gval[mu][e].projectTAH s

let contractProjTAHg = newGfunc(forward = contractProjTAHf, backward = contractProjTAHb, name = "projTAHg")

method contractProjTAH*(x: Ggauge, y: Ggauge): Gvalue = Ggauge(gval: x.gval.newOneOf, inputs: @[Gvalue(x), y], gfunc: contractProjTAHg)

proc axexpb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  # z = exp(ax)
  # D[ax] = exp'(Dz, ax)
  # Da = redot(exp'(Dz, ax), x)
  # Dx = a exp'(Dz, ax)
  if zb == nil:
    raiseValueError("no backward node")
  let dexp = expDeriv(zb, z.inputs[0]*z.inputs[1])
  case i
  of 0:
    return dexp.redot z.inputs[1]
  of 1:
    return z.inputs[0]*dexp
  else:
    raiseValueError("i must be 0 or 1, got: " & $i)

proc axexpf(v: Gvalue) =
  let a = Gscalar(v.inputs[0]).getfloat
  let x = Ggauge(v.inputs[1])
  let z = Ggauge(v)
  threads:
    for mu in 0..<z.gval.len:
      for e in z.gval[mu]:
        z.gval[mu][e] := exp(a*x.gval[mu][e])

let axexpg = newGfunc(forward = axexpf, backward = axexpb, name = "axexp")

method axexp*(a: Gscalar, x: Ggauge): Gvalue = Ggauge(gval: x.gval.newOneOf, inputs: @[Gvalue(a), x], gfunc: axexpg)

proc axexpmulyb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  # z = exp(ax)y
  # D[exp(ax)] = Dz y†
  # D[ax] = exp'(Dz y†, ax)
  # Da = redot(exp'(Dz y†, ax), x)
  # Dx = a exp'(Dz y†, ax)
  # Dy = exp(ax)† Dz
  if zb == nil:
    raiseValueError("no backward node")
  let dexp = expDeriv(zb*z.inputs[2].adj, z.inputs[0]*z.inputs[1])
  case i
  of 0:
    return dexp.redot z.inputs[1]
  of 1:
    return z.inputs[0]*dexp
  of 2:
    return axexp(z.inputs[0], z.inputs[1]).adj*zb
  else:
    raiseValueError("i must be 0, 1, or 2, got: " & $i)

proc axexpmulyf(v: Gvalue) =
  let a = Gscalar(v.inputs[0]).getfloat
  let x = Ggauge(v.inputs[1])
  let y = Ggauge(v.inputs[2])
  let z = Ggauge(v)
  threads:
    for mu in 0..<z.gval.len:
      for e in z.gval[mu]:
        var t{.noinit.}: evalType(x.gval[mu][e])
        t := exp(a*x.gval[mu][e])
        z.gval[mu][e] := t*y.gval[mu][e]

let axexpmulyg = newGfunc(forward = axexpmulyf, backward = axexpmulyb, name = "axexpmuly")

method axexpmuly*(a: Gscalar, x: Ggauge, y: Ggauge): Gvalue = Ggauge(gval: x.gval.newOneOf, inputs: @[Gvalue(a), x, y], gfunc: axexpmulyg)

#
# gauge action
#

type Gactcoeff* {.final.} = ref object of Gvalue
  cval: GaugeActionCoeffs

proc getactcoeff*(x: Gvalue): GaugeActionCoeffs = Gactcoeff(x).cval

converter toGvalue*(x: GaugeActionCoeffs): Gvalue =
  result = Gactcoeff(cval: x)
  result.updated

method newOneOf*(x: Gactcoeff): Gvalue = Gactcoeff()
method valCopy*(z: Gactcoeff, x: Gactcoeff) = z.cval = x.cval

method `$`*(x: Gactcoeff): string = $x.cval

proc mulscb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  if zb == nil:
    raiseValueError("no backward node")
  case i
  of 0:
    return redot(zb, z.inputs[1])
  of 1:
    return z.inputs[0]*zb
  else:
    raiseValueError("i must be 0 or 1, got: " & $i)

proc mulscf(v: Gvalue) =
  let x = Gscalar(v.inputs[0]).getfloat
  let c = Gactcoeff(v.inputs[1])
  let z = Gactcoeff(v)
  z.cval = x * c.cval

let mulsc = newGfunc(forward = mulscf, backward = mulscb, name = "s*c")

method `*`*(x: Gscalar, c: Gactcoeff): Gvalue = Gactcoeff(inputs: @[Gvalue(x), c], gfunc: mulsc)

proc redotccb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  case i
  of 0:
    if zb == nil:
      return z.inputs[1]
    else:
      return zb*z.inputs[1]
  of 1:
    if zb == nil:
      return z.inputs[0]
    else:
      return zb*z.inputs[0]
  else:
    raiseValueError("i must be 0 or 1, got: " & $i)

proc redotccf(v: Gvalue) =
  let x = Gactcoeff(v.inputs[0])
  let y = Gactcoeff(v.inputs[1])
  let z = Gscalar(v)
  var t = 0.0
  for a, b in fields(x.cval, y.cval):
    t += a * b
  z.getfloat = t

let redotcc = newGfunc(forward = redotccf, backward = redotccb, name = "redotcc")

method redot*(x: Gactcoeff, y: Gactcoeff): Gvalue = Gscalar(inputs: @[Gvalue(x), y], gfunc: redotcc)

const
  C1Symanzik = -1.0/12.0  # tree-level
  C1Iwasaki = -0.331
  C1DBW2 = -1.4088
let
  gacWilson = Gactcoeff(cval: GaugeActionCoeffs(plaq: 1.0))
  gacSymanzik = Gactcoeff(cval: GaugeActionCoeffs(plaq: 1.0-8.0*C1Symanzik, rect: C1Symanzik))
  gacIwasaki = Gactcoeff(cval: GaugeActionCoeffs(plaq: 1.0-8.0*C1Iwasaki, rect: C1Iwasaki))
  gacDBW2 = Gactcoeff(cval: GaugeActionCoeffs(plaq: 1.0-8.0*C1DBW2, rect: C1DBW2))

method actWilson*(beta: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("actWilson(" & $beta & ")")
method actSymanzik*(beta: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("actSymanzik(" & $beta & ")")
method actIwasaki*(beta: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("actIwasaki(" & $beta & ")")
method actDBW2*(beta: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("actDBW2(" & $beta & ")")
method actAdj*(beta: Gvalue, adjFac: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("actAdj(" & $beta & "," & $adjFac & ")")

method actWilson*(beta: Gscalar): Gvalue = beta * gacWilson
method actSymanzik*(beta: Gscalar): Gvalue = beta * gacSymanzik
method actIwasaki*(beta: Gscalar): Gvalue = beta * gacIwasaki
method actDBW2*(beta: Gscalar): Gvalue = beta * gacDBW2
method actAdj*(beta: Gscalar, adjFac: Gscalar): Gvalue = beta * updateAt(gacWilson, 3, adjFac)

method gaugeAction*(c: Gvalue, g: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("gaugeAction(" & $c & "," & $g & ")")
method gaugeActionDeriv*(c: Gvalue, g: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("gaugeActionDeriv(" & $c & "," & $g & ")")
method gaugeActionDeriv2*(b: Gvalue, c: Gvalue, g: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("gaugeActionDeriv2(" & $b & "," & $c & "," & $g & ")")

proc gaugeForce*(c: Gvalue, g: Gvalue): Gvalue = contractProjTAH(gaugeActionDeriv(c, g), g)

proc gaugeActionb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  case i
  of 0:
    raiseValueError("unimplemented")
#[ this is only for beta_plaq
    if zb == nil:
      return z/z.inputs[0]
    else:
      return (zb*z)/z.inputs[0]
]#
  of 1:
    if zb == nil:
      return gaugeActionDeriv(z.inputs[0], z.inputs[1])
    else:
      return zb*gaugeActionDeriv(z.inputs[0], z.inputs[1])
  else:
    raiseValueError("i must be 0 or 1, got: " & $i)

proc gaugeActionf(v: Gvalue) =
  let gc = Gactcoeff(v.inputs[0]).cval
  let g = Ggauge(v.inputs[1])
  let z = Gscalar(v)
  if gc.adjplaq == 0:
    z.getfloat = gc.gaugeAction1 g.gval
  elif gc.rect == 0 and gc.pgm == 0:
    z.getfloat = gc.actionA g.gval
  else:
    raiseValueError("Gauge coefficient unsupported: " & $gc)

let gaugeActiong = newGfunc(forward = gaugeActionf, backward = gaugeActionb, name = "gaugeAction")

method gaugeAction*(c: Gactcoeff, g: Ggauge): Gvalue = Gscalar(inputs: @[Gvalue(c), g], gfunc: gaugeActiong)

proc gaugeActionDerivb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  if zb == nil:
    raiseValueError("no backward node")
  case i
  of 0:
    raiseValueError("unimplemented")
#[ this is only for beta_plaq
    return (zb*z)/z.inputs[0]
]#
  of 1:
    return gaugeActionDeriv2(zb, z.inputs[0], z.inputs[1])
  else:
    raiseValueError("i must be 0 or 1, got: " & $i)

proc gaugeActionDerivf(v: Gvalue) =
  var gc = Gactcoeff(v.inputs[0]).cval
  let g = Ggauge(v.inputs[1])
  let z = Ggauge(v)
  for f in gc.fields:  # gaugeActionDeriv/gaugeADeriv has a different sign than we want!
    f = -f
  if gc.adjplaq == 0:
    gc.gaugeActionDeriv(g.gval, z.gval)
  elif gc.rect == 0 and gc.pgm == 0:
    gc.gaugeADeriv(g.gval, z.gval)
  else:
    raiseValueError("Gauge coefficient unsupported: " & $gc)

let gaugeActionDerivg = newGfunc(forward = gaugeActionDerivf, backward = gaugeActionDerivb, name = "gaugeActionDeriv")

method gaugeActionDeriv*(c: Gactcoeff, g: Ggauge): Gvalue = Ggauge(gval: g.gval.newOneOf, inputs: @[Gvalue(c), g], gfunc: gaugeActionDerivg)

proc gaugeActionDeriv2b(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  if zb == nil:
    raiseValueError("no backward node")
  case i
  of 0:
    raiseValueError("unimplemented")
  of 1:
    raiseValueError("unimplemented")
  of 2:
    raiseValueError("unimplemented")
  else:
    raiseValueError("i must be 0, 1, or 2, got: " & $i)

proc gaugeActionDeriv2f(v: Gvalue) =
  let b = Ggauge(v.inputs[0])
  let gc = Gactcoeff(v.inputs[1]).cval
  let g = Ggauge(v.inputs[2])
  let z = Ggauge(v)
  if gc.adjplaq == 0:
    gc.gaugeDerivDeriv2(g.gval, b.gval, z.gval)
  elif gc.rect == 0 and gc.pgm == 0:
    raiseValueError("unimplemented")
  else:
    raiseValueError("Gauge coefficient unsupported: " & $gc)

let gaugeActionDeriv2g = newGfunc(forward = gaugeActionDeriv2f, backward = gaugeActionDeriv2b, name = "gaugeActionDeriv2")

method gaugeActionDeriv2*(b: Ggauge, c: Gactcoeff, g: Ggauge): Gvalue = Ggauge(gval: g.gval.newOneOf, inputs: @[Gvalue(b), c, g], gfunc: gaugeActionDeriv2g)
