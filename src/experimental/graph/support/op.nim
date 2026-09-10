import ../core/base
from ../core/grad_engine import gradSeeded
from ../core/slotvar import slotVar

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

proc secondPullback*[T: Gvalue](x: T, seed: Gvalue, upstream: Gvalue,
                                replica: proc(slot: T): Gvalue): Gvalue =
  ## Computes `d/d slot [ (d replica/d slot)^T seed ]^T upstream` at
  ## `slot = x`.
  ##
  ## `replica` must build the primal over `slot`. Every occurrence of the
  ## differentiated argument must be spelled as `slot`. A builder that captures
  ## `x` directly compiles but returns a wrong partial. `seed` and `upstream` are
  ## ordinary graph values and stay live, so the result stays exact under
  ## further differentiation even when they are `x` or depend on `x`. The
  ## engine walks their paths to `x` through their own slots.
  discard sharedGraphRuntime(
    [Gvalue(x), seed, upstream], "secondPullback")
  let slot = slotVar(x)
  gradSeeded(gradSeeded(replica(slot), slot, seed), slot, upstream)
