import ../core

type
  Gscalar* {.final.} = ref object of Gvalue
    ## Float scalar node. `update` also marks dependent graphs stale.
    sval: float
  Gint* {.final.} = ref object of Gvalue
    ival: int

proc scalarNodeIn*(grt: GraphRuntime): Gscalar =
  Gscalar(runtime: grt)

proc scalarLeafIn*(grt: GraphRuntime,
                   value: float): Gscalar =
  result = Gscalar(
    sval: value,
    runtime: grt)
  result.updated

proc scalarNodeLike*(anchor: Gvalue): Gscalar =
  scalarNodeIn(anchor.runtime)

proc scalarLeafLike*(anchor: Gvalue,
                     value: float): Gscalar =
  if anchor == nil:
    raiseValueError("scalar leaf requires non-nil anchor")
  scalarLeafIn(anchor.runtime, value)

proc scalarLeafLike*(anchor: Gvalue,
                     value: int): Gscalar =
  scalarLeafLike(anchor, float(value))

proc intNodeIn*(grt: GraphRuntime): Gint =
  Gint(runtime: grt)

proc intLeafIn*(grt: GraphRuntime,
                value: int): Gint =
  result = Gint(
    ival: value,
    runtime: grt)
  result.updated

proc intNodeLike*(anchor: Gvalue): Gint =
  intNodeIn(anchor.runtime)

proc intLeafLike*(anchor: Gvalue,
                  value: int): Gint =
  if anchor == nil:
    raiseValueError("int leaf requires non-nil anchor")
  intLeafIn(anchor.runtime, value)

proc numericLeafLike*(anchor: Gvalue,
                      value: int): Gvalue =
  if anchor of Gint:
    return intLeafLike(anchor, value)
  scalarLeafLike(anchor, value)

proc numericLeafLike*(anchor: Gvalue,
                      value: float): Gvalue =
  if anchor == nil:
    raiseValueError("numeric leaf requires non-nil anchor")
  if anchor of Gint:
    raiseValueError("float literal is incompatible with int graph value")
  scalarLeafLike(anchor, value)

proc requireScalar*(value: Gvalue,
                    label: string): Gscalar =
  if value == nil:
    raiseValueError(label & " requires non-nil scalar value")
  if not (value of Gscalar):
    raiseValueError(label & " expects scalar value, got:\n" & value.nodeRepr)
  Gscalar(value)

proc requireInt*(value: Gvalue,
                 label: string): Gint =
  if value == nil:
    raiseValueError(label & " requires non-nil int value")
  if not (value of Gint):
    raiseValueError(label & " expects int value, got:\n" & value.nodeRepr)
  Gint(value)

proc getfloat*(x: Gvalue): float = Gscalar(x).sval

proc `getfloat=`*(x: Gvalue, y: float) =
  let xs = Gscalar(x)
  xs.sval = y

method update*(x: Gscalar, y: float) =
  x.getfloat = y
  x.updated

proc toGvalue*(grt: GraphRuntime,
               x: float): Gscalar =
  scalarLeafIn(grt, x)

proc valCopy*(z: Gvalue,
              x: float) =
  z.valCopy(toGvalue(z.runtime, x))

method newOneOf*(x: Gscalar): Gvalue =
  result = Gscalar(runtime: x.runtime)
method valCopy*(z: Gscalar, x: Gscalar) = z.sval = x.sval
method copyCompatible*(prototype: Gscalar, value: Gscalar): bool =
  prototype != nil and value != nil

method `$`*(x: Gscalar): string = $x.sval

method isZero*(x: Gscalar): bool = x.sval == 0.0

proc `getfloat=`*(x: Gvalue, y: int) =
  let xs = Gscalar(x)
  xs.sval = float(y)

method update*(x: Gscalar, y: int) =
  x.getfloat = y
  x.updated

proc getint*(x: Gvalue): int = Gint(x).ival

proc `getint=`*(x: Gvalue, y: int) =
  let xs = Gint(x)
  xs.ival = y

method update*(x: Gint, y: int) =
  x.getint = y
  x.updated

proc toGvalue*(grt: GraphRuntime,
               x: int): Gint =
  intLeafIn(grt, x)

proc valCopy*(z: Gvalue,
              x: int) =
  z.valCopy(toGvalue(z.runtime, x))

method newOneOf*(x: Gint): Gvalue =
  result = Gint(runtime: x.runtime)
method valCopy*(z: Gint, x: Gint) = z.ival = x.ival
method copyCompatible*(prototype: Gint, value: Gint): bool =
  prototype != nil and value != nil

method `$`*(x: Gint): string = $x.ival

method isZero*(x: Gint): bool = x.ival == 0
