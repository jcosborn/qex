# Gauge-test gradient-check machinery, shared by the gauge test programs via
# `include`. The include site must provide `grt`, `qex`, the graph/gauge
# imports, `math`, and `algorithms/numdiff`.

# Directional-derivative check:
# z(t) = f(x + t A) B†, so d/dt z should match redot(dz/dx, A).

proc ndiff(zt: Gscalar, t: Gscalar, dx = 0.125, ordMax: static int = 3; noise: static bool = false): auto =
  when noise:
    var mag: array[ordMax, float]
    var n = 0
  proc z(v:float):float =
    t.update v
    result = zt.eval.sval
    when noise:
      mag[n div 2] += abs(result)
      inc n
  var dzdt,e: float
  ndiff(dzdt, e, z, 0.0, dx, ordMax=ordMax)
  when noise:
    # D = sum_i w_i [S(h_i)-S(-h_i)]/(2 h_i), h_i = dx/2^i.
    # w_i = product_(j != i) h_j^2/(h_j^2-h_i^2).
    # Sample resolution eps |S| gives eps sum_i |w_i| (|S+|+|S-|)/(2 |h_i|).
    var err = 0.0
    for i in 0..<ordMax:
      let hi = dx / float(1 shl i)
      var w = 1.0
      for j in 0..<ordMax:
        if i != j:
          let hj = dx / float(1 shl j)
          w *= hj*hj/(hj*hj-hi*hi)
      err += abs(w)*mag[i]/(2.0*abs(hi))
    (dzdt, e, DBL_EPSILON*err)
  else:
    (dzdt, e)

template check(ii: tuple[filename:string, line:int, column:int], ast: string, dzdt, e, gdota: float; noise: float = 0.0) =
  let tol = noise + DBL_EPSILON*abs(gdota+dzdt)*float(512*1024)
  # Extra roundoff allowance requires the extrapolation to have converged to it.
  if not (almostEqual(gdota, dzdt, unitsInLastPlace = 512*1024) or
          (noise > 0.0 and abs(dzdt-gdota) <= tol and e <= tol)):
    checkpoint(ii.filename & ":" & $ii.line & ":" & $ii.column & ": Check failed: " & ast)
    checkpoint("  ndiff: " & $dzdt & " +/- " & $e)
    checkpoint("  grad: " & $gdota)
    checkpoint("  reldelta: " & $(abs(dzdt-gdota)/abs(dzdt+gdota)))
    if noise > 0.0:
      checkpoint("  roundoff allowance: " & $noise)
    fail()

template ckforce(s: untyped, f: untyped, x: untyped, p: untyped) =
  let t = grt.toGvalue(0.0)
  # S(exp(t*p)*x): wider steps limit cancellation; six levels cancel h^2 through h^10.
  let (dsdt, e, noise) = ndiff(s(exp(t*p)*x), t, 2.0, ordMax=6, noise=true)
  let pdotf = eval(redot(p, f(x))).sval
  check(instantiationInfo(), astTostr(s(x) -> f(x)), dsdt, e, pdotf, noise)

template ckgrad(f: untyped, x: untyped, a: untyped) =
  let t = grt.toGvalue(0.0)
  let (dzdt, e) = ndiff(f(x+t*a), t)
  let ff = f(x)
  let gdota = eval(redot(grad(ff, x), a)).sval
  check(instantiationInfo(), astTostr(f(x)), dzdt, e, gdota)

template ckgrad2(f: untyped, x: untyped, y: untyped, ax: untyped, ay: untyped) =
  let t = grt.toGvalue(0.0)
  let (dzdt, e) = ndiff(f(x+t*ax, y+t*ay), t)
  let ff = f(x, y)
  let gdota = eval(redot(grad(ff, x), ax) + redot(grad(ff, y), ay)).sval
  check(instantiationInfo(), astTostr(f(x,y)), dzdt, e, gdota)

template ckgradm(f: untyped, x: untyped, a: untyped, b: untyped) =
  let t = grt.toGvalue(0.0)
  let (dzdt, e) = ndiff(f(x+t*a).redot b, t)
  let ff = f(x).redot b
  let gdota = eval(redot(grad(ff, x), a)).sval
  check(instantiationInfo(), astTostr(f(x)), dzdt, e, gdota)

template ckgradm2(f: untyped, x: untyped, y: untyped, ax: untyped, ay: untyped, b: untyped) =
  let t = grt.toGvalue(0.0)
  let (dzdt, e) = ndiff(f(x+t*ax, y+t*ay).redot b, t)
  let ff = f(x, y).redot b
  let gdota = eval(redot(grad(ff, x), ax) + redot(grad(ff, y), ay)).sval
  check(instantiationInfo(), astTostr(f(x,y)), dzdt, e, gdota)

template ckgradm3(f: untyped, x: untyped, y: untyped, u: untyped, ax: untyped, ay: untyped, au: untyped, b: untyped) =
  let t = grt.toGvalue(0.0)
  let (dzdt, e) = ndiff(f(x+t*ax, y+t*ay, u+t*au).redot b, t)
  let ff = f(x, y, u).redot b
  let gdota = eval(redot(grad(ff, x), ax) + redot(grad(ff, y), ay) + redot(grad(ff, u), au)).sval
  check(instantiationInfo(), astToStr(f(x,y,u)), dzdt, e, gdota)

template ckbinarynorm2grad(fusedExpr: untyped,
                           refExpr: untyped,
                           x: untyped,
                           y: untyped,
                           tol: float) =
  let fused = fusedExpr
  let refv = refExpr
  norm2(fused - refv) :< tol
  let fusedNorm2 = fused.norm2
  let refNorm2 = refv.norm2
  norm2(grad(fusedNorm2, x) - grad(refNorm2, x)) :< tol
  norm2(grad(fusedNorm2, y) - grad(refNorm2, y)) :< tol
