import ../core
import callable
import closure/normalize

proc lambda*(param: Gvalue, body: Gvalue): Gvalue =
  if param == nil:
    raiseValueError("lambda parameter cannot be nil")
  if body == nil:
    raiseValueError("lambda body cannot be nil")
  let grt = param.runtime
  # Closure captures must stay in one runtime for cloning and cache identity.
  if grt != body.runtime:
    raiseValueError("lambda mixes multiple graph runtimes")
  let fn = Glambda(param: param, body: body).attachRuntime(grt)
  normalizeClosure(fn)
  result = fn
  result.updated
