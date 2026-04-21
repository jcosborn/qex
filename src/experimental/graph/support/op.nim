import ../core/base

proc requireUpstream*(zb: Gvalue, label: string): Gvalue =
  if zb == nil:
    raiseValueError(label & " requires non-nil upstream gradient")
  zb

proc scaledUpstreamOr*(zb: Gvalue, scale: Gvalue): Gvalue =
  if zb == nil:
    return scale
  zb * scale

proc raiseInputIndexError*(label: string,
                           expected: string,
                           i: int) {.noreturn.} =
  raiseValueError(label & " input index must be " & expected & ", got: " & $i)

proc raiseUnsupportedPath*(label: string,
                           detail = "") {.noreturn.} =
  var msg = label & " is not implemented"
  if detail.len > 0:
    msg &= ": " & detail
  raiseValueError(msg)

template unaryBackwardCase*(label: string, i: int, body0: untyped) =
  case i
  of 0:
    body0
  else:
    raiseInputIndexError(label, "0", i)

template binaryBackwardCase*(label: string, i: int, body0, body1: untyped) =
  case i
  of 0:
    body0
  of 1:
    body1
  else:
    raiseInputIndexError(label, "0 or 1", i)

template ternaryBackwardCase*(label: string, i: int, body0, body1, body2: untyped) =
  case i
  of 0:
    body0
  of 1:
    body1
  of 2:
    body2
  else:
    raiseInputIndexError(label, "0, 1, or 2", i)
