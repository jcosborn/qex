## Bridge between matrix fields and the complex scalar field (1x1 matrices):
##   trace(x)(e) = tr x(e)          x-bar = scale(u, 1)
##   scale(c, x)(e) = c(e) x(e)     c-bar = dot(x, u),     x-bar = scale(adj c, u)
##   dot(x, y)(e) = tr(x(e)^dag y(e))  x-bar = scale(adj u, y), y-bar = scale(u, x)
## under the pairing dL = Re tr(G^dag dF) applied per site to the 1x1 values.
## With these, |tr P|^2 summed over sites is norm2(trace(P)).

import ../[core, scalar]
import ../support/op
import layout, gauge, physics/qcdTypes
import types, basic_ops

proc trace*(x: Gfield): Gcfield
proc scale*(c: Gcfield, x: Gfield): Gfield
proc dot*(x, y: Gfield): Gcfield

proc cfieldNodeLike(x: Gfield): Gcfield =
  when cfieldIsGfield:
    x.fieldNodeLike
  else:
    let f = x.fval.l.newField(ColorMatrixN[1, DComplexV])
    f.zeroFieldStorage
    Gcfield(runtime: x.runtime, fval: f).assignStableNodeId

proc traceb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let x = Gfield(z.inputs[0])
  scale(requireUpstream(zb, "trace backward", Gcfield), unitField(x.runtime, x.fval))

proc tracef(v: Gvalue) =
  let x = Gfield(v.inputs[0])
  let z = Gcfield(v)
  threads:
    for e in z.fval:
      z.fval[e][0,0] := x.fval[e].trace

let traceg = Gfunc(forward: tracef, backward: traceb, name: "trace")

proc trace*(x: Gfield): Gcfield =
  ## Per-site trace as a complex scalar field.
  graphNode(x.cfieldNodeLike, @[Gvalue(x)], traceg, "trace")

proc scaleb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let c = Gcfield(z.inputs[0])
  let x = Gfield(z.inputs[1])
  let u = requireUpstream(zb, "scale backward", Gfield)
  if i == 0:
    return dot(x, u)
  scale(c.adj, u)

proc scalef(v: Gvalue) =
  let c = Gcfield(v.inputs[0])
  let x = Gfield(v.inputs[1])
  let z = Gfield(v)
  threads:
    for e in z.fval:
      z.fval[e] := c.fval[e][0,0] * x.fval[e]

let scaleg = Gfunc(forward: scalef, backward: scaleb, name: "scale")

proc scale*(c: Gcfield, x: Gfield): Gfield =
  ## Per-site complex scalar times matrix.
  if c.fval.l != x.fval.l:
    raiseValueError("scale requires matching field shapes")
  graphNode(x.fieldNodeLike, @[Gvalue(c), Gvalue(x)], scaleg, "scale")

proc dotb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let x = Gfield(z.inputs[0])
  let y = Gfield(z.inputs[1])
  let u = requireUpstream(zb, "dot backward", Gcfield)
  if i == 0:
    return scale(u.adj, y)
  scale(u, x)

proc dotf(v: Gvalue) =
  let x = Gfield(v.inputs[0])
  let y = Gfield(v.inputs[1])
  let z = Gcfield(v)
  threads:
    for e in z.fval:
      z.fval[e][0,0] := dot(x.fval[e], y.fval[e])

let dotg = Gfunc(forward: dotf, backward: dotb, name: "dot")

proc dot*(x, y: Gfield): Gcfield =
  ## Per-site tr(x^dag y) as a complex scalar field.
  x.requireSameFieldShape(y, "dot")
  graphNode(x.cfieldNodeLike, @[Gvalue(x), Gvalue(y)], dotg, "dot")
