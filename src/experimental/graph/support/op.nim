import ../core/base

template scaledUpstreamOr*[U: Gvalue, S: Gvalue](zb: Gvalue,
                                                 upstreamType: typedesc[U],
                                                 scale: S): untyped =
  if zb == nil:
    scale
  else:
    scale.scaleLike upstreamType(zb)

template bilinearBackward*(zb: Gvalue, z: Gvalue, i: int,
                          T: typedesc): untyped {.dirty.} =
  ## Backward for a symmetric scalar-valued bilinear op `<x, y>`: each operand's
  ## adjoint is the other operand scaled by the upstream cotangent. `T` recovers
  ## the concrete operand type from the erased inputs.
  if i == 0:
    scaledUpstreamOr(zb, Gscalar, T(z.inputs[1]))
  else:
    scaledUpstreamOr(zb, Gscalar, T(z.inputs[0]))

proc raiseUnsupportedPath*(label: string,
                           detail = "") {.noreturn.} =
  var msg = label & " is not implemented"
  if detail.len > 0:
    msg &= ": " & detail
  raiseValueError(msg)

template requireUpstream*(zb: Gvalue,
                          label: string,
                          upstreamType: typedesc): untyped =
  block:
    if zb == nil:
      raiseValueError(label & " requires an explicit upstream gradient")
    upstreamType(zb)
