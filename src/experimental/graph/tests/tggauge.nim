import math, strutils, unittest

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))

# basicOps.epsilon collides with fenv.epsilon
import qex except epsilon
import algorithms/numdiff, gauge/stoutsmear
import helpers
import ../[core, scalar, multi, gauge]
import ../functional
import ../gauge/shared as graphGaugeShared
import ../gauge/basic_ops as graphGaugeBasic
import ../hmcgauge/optimizer, ../hmcgauge/integrator
import ../hmcgauge/trajectory
import ../hmcgauge/training
import ../hmcgauge/config
import ../hmcgauge/gauge_io
import ../hmcgauge/params
import ../hmcgauge/rng

let grt = initGraphRuntime()

# Directional-derivative check:
# z(t) = f(x + t A) B†, so d/dt z should match redot(dz/dx, A).

proc ndiff(zt: Gscalar, t: Gscalar): (float, float) =
  proc z(v:float):float =
    t.update v
    zt.eval.sval
  var dzdt,e: float
  ndiff(dzdt, e, z, 0.0, 0.125, ordMax=3)
  (dzdt, e)

template check(ii: tuple[filename:string, line:int, column:int], ast: string, dzdt, e, gdota: float) =
  if not almostEqual(gdota, dzdt, unitsInLastPlace = 512*1024):
    checkpoint(ii.filename & ":" & $ii.line & ":" & $ii.column & ": Check failed: " & ast)
    checkpoint("  ndiff: " & $dzdt & " +/- " & $e)
    checkpoint("  grad: " & $gdota)
    checkpoint("  reldelta: " & $(abs(dzdt-gdota)/abs(dzdt+gdota)))
    fail()

template ckforce(s: untyped, f: untyped, x: untyped, p: untyped) =
  let t = grt.toGvalue(0.0)
  let (dsdt, e) = ndiff(s(exp(t*p)*x), t)
  let pdotf = eval(redot(p, f(x))).sval
  check(instantiationInfo(), astTostr(s(x) -> f(x)), dsdt, e, pdotf)

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

qexInit()

let
  lat = @[8,8,8,16]
  lo = lat.newLayout
  seed = 1234567891u64
  vol = lo.physVol
var
  r = lo.newRNGField(Philox4x64, seed)
  g = lo.newgauge
  u = lo.newgauge
  p = lo.newgauge
  q = lo.newgauge
  m = lo.newgauge
  ss = lo.newStoutSmear(0.1)
const nc = g[0][0].nrows
threads:
  g.random r
  u.random r
  p.randomTAH r
  q.randomTAH r
  m.randomTAH r
for i in 0..4:
  ss.smear(g, g)
  ss.smear(u, u)
threads:
  for t in m:
    t *= 0.01

let scalarValues = sampleScalarValues()
let a = scalarValues.a
let b = scalarValues.b

proc zeroGaugeLike(source: graphGaugeShared.Gauge): graphGaugeShared.Gauge =
  result = source.newOneOf
  graphGaugeShared.zeroGaugeStorage(result)

include gauge/coeffs
include gauge/basic
include gauge/fused
include gauge/action
include hmcgauge/basic

qexFinalize()
