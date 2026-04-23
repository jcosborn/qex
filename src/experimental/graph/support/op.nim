import ../core/base

proc requireUpstream*(zb: Gvalue, label: string): Gvalue =
  if zb == nil:
    raiseValueError(label & " requires non-nil upstream gradient")
  zb

proc requireUpstream*[T: Gvalue](zb: Gvalue,
                                 upstreamType: typedesc[T],
                                 label: string): T =
  let value = requireUpstream(zb, label)
  if not (value of T):
    raiseValueError(
      label & " expects upstream gradient of type " & $upstreamType &
      ", got:\n" & value.nodeRepr)
  T(value)

template scaledUpstreamOr*(zb: Gvalue, scale: Gvalue): untyped =
  if zb == nil:
    scale
  else:
    scale.scaleLike zb

template scaledUpstreamOr*[U: Gvalue, S: Gvalue](zb: Gvalue,
                                                 upstreamType: typedesc[U],
                                                 scale: S,
                                                 label: string): untyped =
  if zb == nil:
    scale
  else:
    scale.scaleLike requireUpstream(zb, upstreamType, label)

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
