import ../core
import ../core/base

type
  Gscalar* {.final.} = ref object of Gvalue
    ## Float scalar node. `update` also marks dependent graphs stale.
    sval*: float
  Gint* {.final.} = ref object of Gvalue
    ival*: int

proc scalarNodeIn*(grt: GraphRuntime): Gscalar =
  Gscalar().attachRuntime(grt)

proc scalarLeafIn*(grt: GraphRuntime,
                   value: float): Gscalar =
  result = Gscalar(
    sval: value).attachRuntime(grt)
  result.updated

proc scalarNodeLike*(anchor: Gvalue): Gscalar =
  scalarNodeIn(anchor.runtime)

proc scalarLeafLike*(anchor: Gvalue,
                     value: float): Gscalar =
  scalarLeafIn(anchor.runtime, value)

proc scalarLeafLike*(anchor: Gvalue,
                     value: int): Gscalar =
  scalarLeafLike(anchor, float(value))

proc intNodeIn*(grt: GraphRuntime): Gint =
  Gint().attachRuntime(grt)

proc intLeafIn*(grt: GraphRuntime,
                value: int): Gint =
  result = Gint(
    ival: value).attachRuntime(grt)
  result.updated

proc intNodeLike*(anchor: Gvalue): Gint =
  intNodeIn(anchor.runtime)

proc intLeafLike*(anchor: Gvalue,
                  value: int): Gint =
  intLeafIn(anchor.runtime, value)

proc numericLeafLike*(anchor: Gvalue,
                      value: int): Gvalue =
  if anchor of Gint:
    return intLeafLike(anchor, value)
  scalarLeafLike(anchor, value)

proc numericLeafLike*(anchor: Gvalue,
                      value: float): Gvalue =
  if anchor of Gint:
    raiseValueError("float literal is incompatible with int graph value")
  scalarLeafLike(anchor, value)

proc requireScalar*(value: Gvalue,
                    label: string): Gscalar =
  if not (value of Gscalar):
    raiseValueError(label & " expects scalar value, got:\n" & value.nodeRepr)
  Gscalar(value)

proc requireInt*(value: Gvalue,
                 label: string): Gint =
  if not (value of Gint):
    raiseValueError(label & " expects int value, got:\n" & value.nodeRepr)
  Gint(value)

proc update*(x: Gscalar, y: float) =
  x.sval = y
  x.updated

proc toGvalue*(grt: GraphRuntime,
               x: float): Gscalar =
  scalarLeafIn(grt, x)

method newOneOf*(x: Gscalar): Gvalue =
  result = Gscalar().attachRuntime(x.runtime)
method oneLike*(x: Gscalar): Gvalue =
  toGvalue(x.runtime, 1.0)
method valCopy*(z: Gscalar, x: Gvalue) =
  z.sval = x.requireScalar("scalar copy").sval
method copyCompatible*(prototype: Gscalar, value: Gvalue): bool =
  value of Gscalar

method `$`*(x: Gscalar): string = $x.sval

method isZero*(x: Gscalar): bool = x.sval == 0.0
method supportsCondSelection*(x: Gscalar): bool = true

proc update*(x: Gscalar, y: int) =
  x.sval = float(y)
  x.updated

proc update*(x: Gint, y: int) =
  x.ival = y
  x.updated

proc toGvalue*(grt: GraphRuntime,
               x: int): Gint =
  intLeafIn(grt, x)

method newOneOf*(x: Gint): Gvalue =
  result = Gint().attachRuntime(x.runtime)
method oneLike*(x: Gint): Gvalue =
  toGvalue(x.runtime, 1)
method valCopy*(z: Gint, x: Gvalue) =
  z.ival = x.requireInt("int copy").ival
method copyCompatible*(prototype: Gint, value: Gvalue): bool =
  value of Gint

method `$`*(x: Gint): string = $x.ival

method isZero*(x: Gint): bool = x.ival == 0
method supportsCondSelection*(x: Gint): bool = true
