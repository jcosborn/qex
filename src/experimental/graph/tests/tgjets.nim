import math, unittest
import qex except epsilon
import algorithms/numdiff
import helpers
import ../[core, scalar, gauge]
from ../gauge/matfun import expJet, expTopReplica

let grt = initGraphRuntime()
include gauge/gaugehelpers

template ckjet(f, x, a: untyped) =
  block:
    let t = grt.toGvalue(0.0)
    let (dv, err) = ndiff(f(x+t*a), t, 0.05, ordMax=5)
    let av = redot(grad(f(x), x), a).eval.sval
    checkpoint("  derivative " & $av & ", finite difference " & $dv &
      ", error estimate " & $err)
    check abs(dv-av) < 1e-9 * max(1.0, abs(dv))

proc runJetTests*() =
  qexInit()
  letParam:
    lat = latticeFromLocalLattice(@[4,4,8,8], nRanks)
  let lo = lat.newLayout
  var
    rng = lo.newRNGField(Philox4x64, 9040911u64)
    y = lo.newGauge
    dirs: array[4, typeof(y)]
    b = lo.newGauge
    refv = lo.newGauge
  for d in dirs.mitems:
    d = lo.newGauge
  threads:
    y.randomTAH rng
    b.random rng
    for d in dirs.mitems:
      d.randomTAH rng
      for f in d:
        f *= 0.1
    for f in y:
      f *= 0.1
  let
    gy = grt.toGvalue(y)
    gb = grt.toGvalue(b)
  var ds: array[4, Ggauge]
  for i in 0..3:
    ds[i] = grt.toGvalue(dirs[i])
  const nc = y[0][0].nrows

  when nc == 1:
    proc cplx(re, im: float): Ggauge =
      let g = lo.newGauge
      threads:
        for mu in 0..<g.len:
          for e in g[mu]:
            g[mu][e][0,0].re := re
            g[mu][e][0,0].im := im
      grt.toGvalue(g)

  suite "exponential jet dispatch boundary":
    test "four directions agree with the numerical jet":
      threads:
        for mu in 0..<y.len:
          for e in y[mu]:
            var yy, outv {.noinit.}: evalType(y[mu][e])
            var dd {.noinit.}: array[4, evalType(y[mu][e])]
            yy := y[mu][e]
            for i in 0..3:
              dd[i] := dirs[i][mu][e]
            outv[] := expTop(yy[], [dd[0][], dd[1][], dd[2][], dd[3][]])
            refv[mu][e] := outv
      check norm2(expJet(gy, ds) - grt.toGvalue(refv)).eval.sval < 1e-19

    test "backward of three directions agrees with finite differences":
      proc loss(x: Ggauge): Gscalar = redot(expJet(x, ds[0..2]), gb)
      ckjet(loss, gy, ds[3])

    when nc == 1:
      test "zero and four directions preserve the scalar identity":
        check norm2(expJet(gy, []) - exp(gy)).eval.sval < 1e-24
        check norm2(expJet(gy, ds) - exp(gy)*ds[0]*ds[1]*ds[2]*ds[3]).eval.sval < 1e-24
        check norm2(expTopReplica(gy, []) - exp(gy)).eval.sval < 1e-24

      test "large imaginary scalar exponent uses exact exp at zero and four directions":
        let
          x = cplx(0.0, 10000.0)
          cs = [newComplex(0.7, -0.4), newComplex(-0.6, 0.5),
                newComplex(0.9, 0.2), newComplex(0.3, -0.8)]
          z = newComplex(cos(10000.0), sin(10000.0))
          w = z * cs[0] * cs[1] * cs[2] * cs[3]
          ref0 = cplx(z.re, z.im)
          ref4 = cplx(w.re, w.im)
        var d: array[4, Ggauge]
        for i in 0..3:
          d[i] = cplx(cs[i].re, cs[i].im)
        # The scale/square polynomial has an O(1e-6) phase defect at 10000i.
        check norm2(expJet(x, []) - ref0).eval.sval < 1e-24
        check norm2(expTopReplica(x, []) - ref0).eval.sval < 1e-24
        check norm2(expJet(x, d) - ref4).eval.sval < 1e-24

      test "fallback pullbacks keep repeated directions and cotangents live":
        let
          x = cplx(0.12, -0.19)
          d = cplx(0.7, 0.4)
          a = cplx(-0.3, 0.6)
          b = cplx(0.5, -0.2)
          c = cplx(-0.4, -0.7)
        proc loss(x, d: Ggauge): Gscalar =
          # x aliases the exponent, d repeats, x*d depends on both slots.
          redot(expJet(x, [x, d, d, x*d]), b + x*c + d*a)
        proc fx(x: Ggauge): Gscalar = loss(x, d)
        proc fd(d: Ggauge): Gscalar = loss(x, d)
        proc fx2(x: Ggauge): Gscalar = redot(grad(fx(x), x), a + x*b)
        proc fd2(d: Ggauge): Gscalar = redot(grad(fd(d), d), c + d*b)
        ckjet(fx, x, a)
        ckjet(fd, d, c)
        ckjet(fx2, x, c)
        ckjet(fd2, d, a)
  qexFinalize()

when isMainModule:
  runJetTests()
