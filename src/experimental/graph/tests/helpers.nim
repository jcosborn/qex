import math, unittest

import ../[core, scalar]

type
  SampleScalarValues* = tuple[a, b: float]
  SampleScalarSweep* = tuple[a, b, c, d: float]
  ScalarLeafPair* = tuple[a, b: float, x, y: Gvalue]
  MutableScalarPair* = tuple[a, b, c, d: float, x, y: Gscalar]

proc sampleScalarValues*(): SampleScalarValues =
  (
    a: 0.5 * (sqrt(5.0) - 1.0),
    b: sqrt(2.0) - 1.0)

proc sampleScalarSweep*(): SampleScalarSweep =
  let values = sampleScalarValues()
  (
    a: values.a,
    b: values.b,
    c: 2.0 * values.a - 1.0,
    d: values.a + 3.0 * values.b - 1.0)

proc initScalarLeafPair*(grt: GraphRuntime): ScalarLeafPair =
  let values = sampleScalarValues()
  (
    a: values.a,
    b: values.b,
    x: toGvalue(grt, values.a),
    y: toGvalue(grt, values.b))

proc initMutableScalarPair*(grt: GraphRuntime): MutableScalarPair =
  let values = sampleScalarSweep()
  let x = toGvalue(grt, 0.0)
  let y = toGvalue(grt, 0.0)
  x.update values.a
  y.update values.b
  (
    a: values.a,
    b: values.b,
    c: values.c,
    d: values.d,
    x: x,
    y: y)

template checkeq*(ii: tuple[filename: string, line: int, column: int],
                  sa: string,
                  a: float,
                  sb: string,
                  b: float) =
  if not almostEqual(a, b, unitsInLastPlace = 64):
    checkpoint(ii.filename & ":" & $ii.line & ":" & $ii.column & ": Check failed: " & sa & " :~ " & sb)
    checkpoint("  " & sa & ": " & $a)
    checkpoint("  " & sb & ": " & $b)
    fail()

template checkeq*(ii: tuple[filename: string, line: int, column: int],
                  sa: string,
                  a: int,
                  sb: string,
                  b: int) =
  if a != b:
    checkpoint(ii.filename & ":" & $ii.line & ":" & $ii.column & ": Check failed: " & sa & " :~ " & sb)
    checkpoint("  " & sa & ": " & $a)
    checkpoint("  " & sb & ": " & $b)
    fail()

template `:~`*(a: Gvalue, b: float) =
  checkeq(instantiationInfo(), astToStr a, a.eval.getfloat, astToStr b, b)

template `:~`*(a: Gvalue, b: int) =
  checkeq(instantiationInfo(), astToStr a, a.eval.getint, astToStr b, b)

template `:~`*(a: Gvalue, b: Gvalue) =
  let av = a.eval
  let bv = b.eval
  if (av of Gscalar) and (bv of Gscalar):
    checkeq(instantiationInfo(), astToStr a, av.getfloat, astToStr b, bv.getfloat)
  elif (av of Gint) and (bv of Gint):
    checkeq(instantiationInfo(), astToStr a, av.getint, astToStr b, bv.getint)
  else:
    raise newException(GraphValueError,
      "Gvalue :~ Gvalue only supports scalar or int nodes; use norm-based checks for gauge values")

template `:<`*(a: Gvalue, b: float) =
  let av = abs(a.eval.getfloat)
  if av >= b:
    let ii = instantiationInfo()
    let sa = astToStr a
    let sb = astToStr b
    checkpoint(ii.filename & ":" & $ii.line & ":" & $ii.column & ": Check failed: " & sa & " :< " & sb)
    checkpoint("  " & sa & ": " & $av)
    checkpoint("  " & sb & ": " & $b)
    fail()
