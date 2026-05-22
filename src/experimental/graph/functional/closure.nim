import ../core
import callable
import closure/normalize

proc lambda*(param: Gvalue, body: Gvalue): Gvalue =
  # Closure captures must stay in one runtime for cloning and cache identity.
  let grt = requireSameGraphRuntime(param, body, "lambda", "parameter", "body")
  let fn = Glambda(param: param, body: body).attachRuntime(grt)
  normalizeClosure(fn)
  result = fn
  result.updated
