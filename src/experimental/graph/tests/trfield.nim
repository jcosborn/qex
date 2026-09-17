## Field bridge contracts independent of Python and checkpoint fixtures.
import qex
import std/[unittest, math]
import ../[core, scalar, functional, plan]
import ../nn as neural
import ../gauge/[types, basic_ops, rfield]
import ../../../nn as numeric

proc runBridges*[T: SomeFloat](precision: typedesc[T]) =
  let lo = newLayout(@[8,8])
  let rt = initGraphRuntime()
  let re = @[numeric.realField(lo,T)]
  let im = @[numeric.realField(lo,T)]
  for s in 0..<lo.nSites:
    re[0]{s} := T(0.25)+T(0.125)*T(lo.coords[0][s])
    im[0]{s} := T(-0.5)+T(0.0625)*T(lo.coords[1][s])
  let r = neural.toGvalue(rt,re)
  let i = neural.toGvalue(rt,im)
  let z = rfield.complex(r,i)
  let tol = when T is float32: 2e-6 else: 2e-13

  proc checkFields(a,b: Greal[T]; factor = 1.0) =
    discard a.eval
    discard b.eval
    check a.fval.len == b.fval.len
    for c in 0..<a.fval.len:
      for s in 0..<lo.nSites:
        var x,y: float
        x := a.fval[c]{s}
        y := b.fval[c]{s}
        check abs(x-factor*y) <= tol*max(1.0,abs(factor*y))

  suite "real field bridge " & $T:
    test "real and imaginary components retain channel values":
      check z.runCount == 0
      checkFields(rfield.real(z,T),r)
      checkFields(rfield.imag(z,T),i)

    test "complex embedding and extraction have the expected adjoints":
      let loss = norm2(z)
      checkFields(grad(loss,r),r,2.0)
      checkFields(grad(loss,i),i,2.0)
      let zr = rfield.real(z,T)
      let zi = rfield.imag(z,T)
      checkFields(grad(neural.redot(zr,zr),r),r,2.0)
      checkFields(grad(neural.redot(zi,zi),i),i,2.0)

    test "derivatives preserve aliased and dependent seeds":
      let dr = grad(norm2(z),r)
      checkFields(grad(neural.redot(dr,r),r),r,4.0)
      checkFields(grad(neural.redot(dr,r*r),r),r*r,6.0)

    test "real weighting works for matrix fields":
      let raw = lo.ColorMatrixD()
      threads: raw := 2.0
      let x = toGfield(rt,raw)
      let outg = rfield.scale(r,x)
      let loss = norm2(outg)
      const nc = raw[0].nrows
      checkFields(grad(loss,r),r,8.0*float(nc))
      let dx = grad(loss,x)
      discard dx.eval
      for s in 0..<lo.nSites:
        var a,b: float
        a := dx.fval{s}[0,0].re
        b := r.fval[0]{s}
        check abs(a-4*b*b) < tol

    test "updates and lambda clones preserve input ownership":
      let formal = Greal[T](r.newOneOf)
      let fun = lambda(formal,rfield.complex(formal,i))
      let applied = Gcfield(apply(fun,r))
      checkFields(rfield.real(applied,T),r)
      let before = z.runCount
      r.update(im)
      checkFields(rfield.real(applied,T),r)
      checkFields(rfield.real(z,T),r)
      check z.runCount == before+1
      discard z.eval
      check z.runCount == before+1
      r.update(re)

    test "plans reuse bridge storage after updates and clear":
      let value = rfield.real(z,T)
      let loss = norm2(z)
      let dr = grad(loss,r)
      let run = plan(value,loss,dr)
      for k in 0..1:
        if k == 1: r.update(im)
        let want = loss.eval.sval
        discard run.eval
        checkFields(Greal[T](run[0]),r)
        check abs(Gscalar(run[1]).sval-want) < tol
        checkFields(Greal[T](run[2]),r,2.0)
        let forwards = run.stats.forwards
        discard run.eval
        check run.stats.forwards == forwards
        run.clear
      r.update(re)

    test "construction rejects incompatible channels and layouts":
      let two = neural.toGvalue(rt,@[re[0],im[0]])
      expect GraphValueError: discard rfield.complex(two,i)
      let other = newLayout(@[8,12])
      let raw = other.ColorMatrixD()
      threads: raw := 1.0
      expect GraphValueError: discard rfield.scale(r,toGfield(rt,raw))

when isMainModule:
  qexInit()
  runBridges(float32)
  runBridges(float64)
  qexFinalize()
