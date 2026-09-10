## Single-direction gauge-field graph values and grad-complete basic ops.
##
## `Gfield` holds one lattice color-matrix field (one direction of a `Gauge`).
## Applications keep whole-`Ggauge` signatures; `Gfield` exists so that
## cross-direction expressions (Wilson lines, staples) and the per-direction
## cotangents of gauge derivatives can be written as ordinary graph values.
##
## Closure contract: backward graphs use this module, the `Ggauge` algebra in
## `basic_ops` (which supplies `addLike` for gauge-valued `injectLink` adjoints),
## constant leaves, and the scalar layer. Graphs composed from these operations
## can be differentiated to arbitrary order by the core gradient engine.

import ../[core, scalar]
import ../support/op
import layout, gauge, physics/qcdTypes
import shared, basic_ops

type
  Gfield* = ref object of Gvalue
    ## Graph-owned single-direction storage; public writes must mark freshness.
    fval*: DLatticeColorMatrixV

type FieldLiteral = int | float

proc requireSameFieldShape(x: Gfield, y: Gfield, label: string) =
  if x.fval.l != y.fval.l:
    raiseValueError(label & " requires matching field shapes")

proc zeroFieldStorage*(f: DLatticeColorMatrixV) =
  threads:
    f := 0.0

proc requireLinkShape*(g: Ggauge, mu: int, f: DLatticeColorMatrixV, label: string) =
  if mu < 0 or mu >= g.gval.len:
    raiseValueError(label & " direction out of range")
  if f.l != g.gval[mu].l:
    raiseValueError(label & " requires matching field shapes")

proc fieldNodeLike(x: Gfield): Gfield =
  let f = x.fval.newOneOf
  f.zeroFieldStorage
  Gfield(runtime: x.runtime, fval: f).assignStableNodeId

proc sameShapeFieldNodeLike(x: Gfield, y: Gfield, label: string): Gfield =
  x.requireSameFieldShape(y, label)
  x.fieldNodeLike

proc unitField(grt: GraphRuntime, proto: DLatticeColorMatrixV): Gfield =
  ## Constant identity-matrix field leaf.
  let f = proto.newOneOf
  threads:
    f := 1.0
  result = Gfield(runtime: grt, fval: f).assignStableNodeId
  result.updated

proc unitFieldLike*(g: Ggauge): Gfield =
  ## Constant identity-matrix field leaf shaped like one direction of g.
  unitField(g.runtime, g.gval[0])

method newOneOf*(x: Gfield): Gvalue =
  x.fieldNodeLike

method zeroLike*(x: Gfield): Gvalue =
  result = x.fieldNodeLike
  result.staticZeroLeaf = true

method isZero*(x: Gfield): bool =
  ## Field zero leaves are marked when constructed; other fields are not scanned.
  x.staticZeroLeaf

method valCopy*(z: Gfield, x: Gvalue) =
  let src = Gfield(x)
  if z.fval.l != src.fval.l:
    raiseValueError("field copy requires matching field shapes")
  threads:
    z.fval := src.fval

method copyCompatible*(prototype: Gfield, value: Gvalue): bool =
  value of Gfield and prototype.fval.l == Gfield(value).fval.l

method `$`*(x: Gfield): string =
  let v = x.fval[0][0,0]
  result = "GaugeField (" & $v.re[0] & ", " & $v.im[0] & ")"

# --- forward declarations for backward builders -------------------------------

proc retr*(x: Gfield): Gscalar
proc adj*(x: Gfield): Gfield
proc norm2*(x: Gfield): Gscalar
proc redot*(x: Gfield, y: Gfield): Gscalar
proc projTAH*(x: Gfield): Gfield
proc shift*(x: Gfield, dir, len: int): Gfield
proc linkField*(g: Ggauge, mu: int): Gfield
proc injectLink*(x: Gfield, mu: int, like: Ggauge): Ggauge
proc `-`*(x: Gfield): Gfield
proc `+`*(x: Gscalar, y: Gfield): Gfield
proc `+`*(x: Gfield, y: Gfield): Gfield
proc `*`*(x: Gscalar, y: Gfield): Gfield
proc `*`*(x: Gfield, y: Gfield): Gfield
proc `-`*(x: Gfield, y: Gfield): Gfield

proc `+`*[T: FieldLiteral](x: T, y: Gfield): Gfield =
  toGvalue(y.runtime, float(x)) + y

proc `*`*[T: FieldLiteral](x: T, y: Gfield): Gfield =
  toGvalue(y.runtime, float(x)) * y

method addLike*(prototype: Gfield, x: Gvalue, y: Gvalue): Gvalue =
  Gfield(x) + Gfield(y)

method scaleLike*(contribution: Gfield, upstream: Gvalue): Gvalue =
  if upstream of Gscalar:
    return Gscalar(upstream) * contribution
  if upstream of Gfield:
    return Gfield(upstream) * contribution
  raiseValueError("field scale upstream expects scalar or field value, got:\n" & upstream.nodeRepr)

# --- site-local algebra --------------------------------------------------------

proc retrfb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let x = Gfield(z.inputs[0])
  scaledUpstreamOr(zb, Gscalar, unitField(x.runtime, x.fval))

proc retrff(v: Gvalue) =
  let x = Gfield(v.inputs[0])
  let z = Gscalar(v)
  threads:
    let t = x.fval.trace.re
    threadMaster: z.sval = t

let retrfg = Gfunc(forward: retrff, backward: retrfb, name: "retrf")

proc retr*(x: Gfield): Gscalar =
  graphNode(scalarNodeLike(x), @[Gvalue(x)], retrfg, "retrf")

proc adjfb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  requireUpstream(zb, "adjf backward", Gfield).adj

proc adjff(v: Gvalue) =
  let x = Gfield(v.inputs[0])
  let z = Gfield(v)
  threads:
    z.fval := x.fval.adj

let adjfg = Gfunc(forward: adjff, backward: adjfb, name: "adjf")

proc adj*(x: Gfield): Gfield =
  graphNode(x.fieldNodeLike, @[Gvalue(x)], adjfg, "adjf")

proc norm2fb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let x = Gfield(z.inputs[0])
  scaledUpstreamOr(zb, Gscalar, 2.0 * x)

proc norm2ff(v: Gvalue) =
  let x = Gfield(v.inputs[0])
  let z = Gscalar(v)
  threads:
    let t = x.fval.norm2
    threadMaster: z.sval = t

let norm2fg = Gfunc(forward: norm2ff, backward: norm2fb, name: "norm2f")

proc norm2*(x: Gfield): Gscalar =
  graphNode(scalarNodeLike(x), @[Gvalue(x)], norm2fg, "norm2f")

proc redotfb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  bilinearBackward(zb, z, i, Gfield)

proc redotff(v: Gvalue) =
  let x = Gfield(v.inputs[0])
  let y = Gfield(v.inputs[1])
  let z = Gscalar(v)
  threads:
    let t = redot(x.fval, y.fval)
    threadMaster: z.sval = t

let redotfg = Gfunc(forward: redotff, backward: redotfb, name: "redotf")

proc redot*(x: Gfield, y: Gfield): Gscalar =
  x.requireSameFieldShape(y, "redotf")
  graphNode(scalarNodeLike(x), @[Gvalue(x), Gvalue(y)], redotfg, "redotf")

proc negfb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  -requireUpstream(zb, "-f backward", Gfield)

proc negff(v: Gvalue) =
  let x = Gfield(v.inputs[0])
  let z = Gfield(v)
  threads:
    z.fval := -x.fval

let negfg = Gfunc(forward: negff, backward: negfb, name: "-f")

proc `-`*(x: Gfield): Gfield =
  graphNode(x.fieldNodeLike, @[Gvalue(x)], negfg, "-f")

proc addsfb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let upstream = requireUpstream(zb, "s+f backward", Gfield)
  if i == 0:
    return retr(upstream)
  upstream

proc addsff(v: Gvalue) =
  let x = Gscalar(v.inputs[0])
  let y = Gfield(v.inputs[1])
  let z = Gfield(v)
  threads:
    z.fval := x.sval + y.fval

let addsfg = Gfunc(forward: addsff, backward: addsfb, name: "s+f")

proc `+`*(x: Gscalar, y: Gfield): Gfield =
  graphNode(y.fieldNodeLike, @[Gvalue(x), Gvalue(y)], addsfg, "s+f")

proc addffb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  requireUpstream(zb, "f+f backward", Gfield)

proc addfff(v: Gvalue) =
  let x = Gfield(v.inputs[0])
  let y = Gfield(v.inputs[1])
  let z = Gfield(v)
  threads:
    z.fval := x.fval + y.fval

let addffg = Gfunc(forward: addfff, backward: addffb, name: "f+f")

proc `+`*(x: Gfield, y: Gfield): Gfield =
  graphNode(sameShapeFieldNodeLike(x, y, "f+f"), @[Gvalue(x), Gvalue(y)], addffg, "f+f")

proc subffb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let upstream = requireUpstream(zb, "f-f backward", Gfield)
  if i == 0:
    return upstream
  -upstream

proc subfff(v: Gvalue) =
  let x = Gfield(v.inputs[0])
  let y = Gfield(v.inputs[1])
  let z = Gfield(v)
  threads:
    z.fval := x.fval - y.fval

let subffg = Gfunc(forward: subfff, backward: subffb, name: "f-f")

proc `-`*(x: Gfield, y: Gfield): Gfield =
  graphNode(sameShapeFieldNodeLike(x, y, "f-f"), @[Gvalue(x), Gvalue(y)], subffg, "f-f")

proc mulsfb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let x = Gscalar(z.inputs[0])
  let y = Gfield(z.inputs[1])
  let upstream = requireUpstream(zb, "s*f backward", Gfield)
  if i == 0:
    return redot(upstream, y)
  x * upstream

proc mulsff(v: Gvalue) =
  let x = Gscalar(v.inputs[0])
  let y = Gfield(v.inputs[1])
  let z = Gfield(v)
  threads:
    z.fval := x.sval * y.fval

let mulsfg = Gfunc(forward: mulsff, backward: mulsfb, name: "s*f")

proc `*`*(x: Gscalar, y: Gfield): Gfield =
  graphNode(y.fieldNodeLike, @[Gvalue(x), Gvalue(y)], mulsfg, "s*f")

proc mulffb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let x = Gfield(z.inputs[0])
  let y = Gfield(z.inputs[1])
  let upstream = requireUpstream(zb, "f*f backward", Gfield)
  if i == 0:
    return upstream * y.adj
  x.adj * upstream

proc mulfff(v: Gvalue) =
  let x = Gfield(v.inputs[0])
  let y = Gfield(v.inputs[1])
  let z = Gfield(v)
  threads:
    z.fval := x.fval * y.fval

let mulffg = Gfunc(forward: mulfff, backward: mulffb, name: "f*f")

proc `*`*(x: Gfield, y: Gfield): Gfield =
  graphNode(sameShapeFieldNodeLike(x, y, "f*f"), @[Gvalue(x), Gvalue(y)], mulffg, "f*f")

proc projTAHfb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  projTAH(requireUpstream(zb, "projTAHf backward", Gfield))

proc projTAHff(v: Gvalue) =
  let x = Gfield(v.inputs[0])
  let z = Gfield(v)
  threads:
    for e in z.fval:
      z.fval[e].projectTAH(x.fval[e])

let projTAHfg = Gfunc(forward: projTAHff, backward: projTAHfb, name: "projTAHf")

proc projTAH*(x: Gfield): Gfield =
  graphNode(x.fieldNodeLike, @[Gvalue(x)], projTAHfg, "projTAHf")

# --- lattice shift --------------------------------------------------------------

type GfieldShift = ref object of Gfield
  ## Owns the shift comm buffers; `newOneOf` clones them so generic graph
  ## cloning preserves the node (DESIGN.md exceptional-node checklist).
  sh: Shifter[DLatticeColorMatrixV, DColorMatrixV]
  dir, len: int

proc shiftNodeLike(x: Gfield, dir, len: int): GfieldShift =
  # Node storage aliases the shifter's receive buffer, so the apply writes
  # the result directly into this node's value; newOneOf clones them as a pair.
  result = GfieldShift(
    runtime: x.runtime,
    sh: newShifter(x.fval, dir, len),
    dir: dir,
    len: len)
  result.fval = result.sh.field
  result.fval.zeroFieldStorage
  result.assignStableNodeId

method newOneOf(x: GfieldShift): Gvalue =
  shiftNodeLike(x, x.dir, x.len)

proc shiftfb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let s = GfieldShift(z)
  shift(requireUpstream(zb, "shiftf backward", Gfield), s.dir, -s.len)

proc shiftff(v: Gvalue) =
  let x = Gfield(v.inputs[0])
  let z = GfieldShift(v)
  threads:
    discard z.sh ^* x.fval

let shiftfg = Gfunc(forward: shiftff, backward: shiftfb, name: "shiftf")

proc shift*(x: Gfield, dir, len: int): Gfield =
  ## x(pos + len*dir); the adjoint of a lattice shift is the opposite shift.
  if dir < 0 or dir >= x.fval.l.nDim:
    raiseValueError("shift direction out of range")
  graphNode(x.shiftNodeLike(dir, len), @[Gvalue(x)], shiftfg, "shiftf")

# --- direction extraction / injection -------------------------------------------

proc linkField*(g: Ggauge, mu: int): Gfield =
  ## One direction of a gauge field as a single-field graph value.
  if mu < 0 or mu >= g.gval.len:
    raiseValueError("linkField direction out of range")
  let f = g.gval[mu].newOneOf
  f.zeroFieldStorage
  proc forward(v: Gvalue) =
    let
      g = Ggauge(v.inputs[0])
      z = Gfield(v)
    threads:
      z.fval := g.gval[mu]
  proc backward(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
    injectLink(
      requireUpstream(zb, "linkField backward", Gfield),
      mu,
      Ggauge(z.inputs[0]))
  graphNode(
    Gfield(runtime: g.runtime, fval: f),
    @[Gvalue(g)],
    Gfunc(forward: forward, backward: backward, name: "linkField"),
    "linkField")

proc injectLink*(x: Gfield, mu: int, like: Ggauge): Ggauge =
  ## Gauge value with x in direction slot mu and zero elsewhere.
  ## `like` supplies the gauge shape only; it is not a graph dependency.
  like.requireLinkShape(mu, x.fval, "injectLink")
  proc forward(v: Gvalue) =
    # Off-direction slots stay zero: gaugeNodeLike zeroes construction and clones.
    let
      x = Gfield(v.inputs[0])
      z = Ggauge(v)
    threads:
      z.gval[mu] := x.fval
  proc backward(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
    linkField(requireUpstream(zb, "injectLink backward", Ggauge), mu)
  graphNode(
    like.gaugeNodeLike,
    @[Gvalue(x)],
    Gfunc(forward: forward, backward: backward, name: "injectLink"),
    "injectLink")
