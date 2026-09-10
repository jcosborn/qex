## Lattice moves of a single-direction field: shift, and the direction
## extraction/injection that connects Gfield with the gauge bundle. The site
## algebra of both types is stamped in basic_ops; this module closes with it.

import ../[core, scalar]
import ../support/op
import layout, gauge, physics/qcdTypes
import types, basic_ops

proc shift*(x: Gfield, dir, len: int): Gfield
proc linkField*(g: Ggauge, mu: int): Gfield
proc injectLink*(x: Gfield, mu: int, like: Ggauge): Ggauge

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
