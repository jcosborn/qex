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

method ensureStorage*(x: GfieldShift) =
  procCall Gfield(x).ensureStorage
  if x.sh.sb.sb.isNil:
    x.sh = newShifter(x.fval, x.dir, x.len, dest = x.fval)
  x.sh.field = x.fval

method bindBuffer*(x: GfieldShift, buffer: Gvalue) =
  procCall Gfield(x).bindBuffer(buffer)
  x.sh.field = x.fval

method releaseWork*(x: GfieldShift) =
  x.sh = default(typeof(x.sh))

method releaseStorage*(x: GfieldShift) =
  x.releaseWork
  procCall Gfield(x).releaseStorage

proc shiftNodeLike(x: Gfield, dir, len: int): GfieldShift =
  result = GfieldShift(runtime: x.runtime, fval: x.fval.newShape, dir: dir, len: len)
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

let shiftfg = Gfunc(bufferMode: bmFull, forward: shiftff, backward: shiftfb, name: "shiftf")

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
  let f = g.gval[mu].newShape
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
    Gfunc(bufferMode: bmFull, forward: forward, backward: backward, name: "linkField"),
    "linkField")

proc injectLink*(x: Gfield, mu: int, like: Ggauge): Ggauge =
  ## Gauge value with x in direction slot mu and zero elsewhere.
  ## `like` supplies the gauge shape only; it is not a graph dependency.
  like.requireLinkShape(mu, x.fval, "injectLink")
  proc forward(v: Gvalue) =
    # bmZero clears the other directions before this output buffer is reused.
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
    Gfunc(bufferMode: bmZero, forward: forward, backward: backward, name: "injectLink"),
    "injectLink")
