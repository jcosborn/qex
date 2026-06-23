import ../core

type
  Gscalar* {.final.} = ref object of Gvalue
    ## Float scalar node. `update` also marks dependent graphs stale.
    sval*: float
  Gint* {.final.} = ref object of Gvalue
    ival*: int

proc toGvalue*(grt: GraphRuntime,
               x: float): Gscalar =
  result = Gscalar(
    runtime: grt,
    sval: x).assignStableNodeId
  result.updated

proc scalarNodeLike*(anchor: Gvalue): Gscalar =
  Gscalar(runtime: anchor.runtime).assignStableNodeId

proc toGvalue*(grt: GraphRuntime,
               x: int): Gint =
  result = Gint(
    runtime: grt,
    ival: x).assignStableNodeId
  result.updated

proc intNodeLike*(anchor: Gvalue): Gint =
  Gint(runtime: anchor.runtime).assignStableNodeId

proc numericLeafLike*(anchor: Gvalue,
                      value: int): Gvalue =
  if anchor of Gint:
    return toGvalue(anchor.runtime, value)
  toGvalue(anchor.runtime, float(value))

proc numericLeafLike*(anchor: Gvalue,
                      value: float): Gvalue =
  if anchor of Gint:
    raiseValueError("float literal is incompatible with int graph value")
  toGvalue(anchor.runtime, value)

# Fresh, current scalar/int leaf nodes; used as local variable / lambda-parameter
# placeholders. They are just the zero-valued leaf constructors.
proc localScalar*(grt: GraphRuntime): Gscalar = toGvalue(grt, 0.0)
proc localInt*(grt: GraphRuntime): Gint = toGvalue(grt, 0)

proc update*(x: Gscalar, y: float) =
  x.sval = y
  x.updated

method newOneOf*(x: Gscalar): Gvalue =
  result = scalarNodeLike(x)
method zeroLike*(x: Gscalar): Gvalue =
  result = scalarNodeLike(x)
  result.markStaticZeroLeaf
method oneLike*(x: Gscalar): Gvalue =
  toGvalue(x.runtime, 1.0)
method valCopy*(z: Gscalar, x: Gvalue) =
  z.sval = Gscalar(x).sval
method copyCompatible*(prototype: Gscalar, value: Gvalue): bool =
  value of Gscalar

method `$`*(x: Gscalar): string = $x.sval

method isZero*(x: Gscalar): bool = x.sval == 0.0

proc update*(x: Gint, y: int) =
  x.ival = y
  x.updated

method newOneOf*(x: Gint): Gvalue =
  result = intNodeLike(x)
method zeroLike*(x: Gint): Gvalue =
  result = intNodeLike(x)
  result.markStaticZeroLeaf
method oneLike*(x: Gint): Gvalue =
  toGvalue(x.runtime, 1)
method valCopy*(z: Gint, x: Gvalue) =
  z.ival = Gint(x).ival
method copyCompatible*(prototype: Gint, value: Gvalue): bool =
  value of Gint

method `$`*(x: Gint): string = $x.ival

method isZero*(x: Gint): bool = x.ival == 0
