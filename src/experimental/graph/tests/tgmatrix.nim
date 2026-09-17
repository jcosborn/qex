import math, unittest
import qex except epsilon
import algorithms/numdiff
import ../[core, scalar]
import ../gauge/matrix
import helpers

let grt = initGraphRuntime()
include gauge/gaugehelpers

template ckmat(f, x, a: untyped) =
  block:
    let t = grt.toGvalue(0.0)
    let (dv, err) = ndiff(f(x+t*a), t, 0.03, ordMax=5)
    let av = redot(grad(f(x),x),a).eval.sval
    echo "  matrix derivative residual ", abs(dv-av), ", FD estimate ", err, ", derivative ", av
    check abs(dv-av) < 1e-9*max(1.0,abs(dv))

template ckmat2(f, x, y, a, b: untyped) =
  block:
    let t = grt.toGvalue(0.0)
    let (dv, err) = ndiff(f(x+t*a,y+t*b), t, 0.03, ordMax=5)
    let ff = f(x,y)
    let av = (redot(grad(ff,x),a)+redot(grad(ff,y),b)).eval.sval
    echo "  matrix derivative residual ", abs(dv-av), ", FD estimate ", err, ", derivative ", av
    check abs(dv-av) < 1e-9*max(1.0,abs(dv))

qexInit()
letParam:
  lat = latticeFromLocalLattice(@[4,4,8,8], nRanks)
let lo = lat.newLayout
let vol = float(lo.physVol)
let a = lo.RealMatrix(8)
let b = lo.RealMatrix(8)
let d = lo.RealMatrix(8)
let r = lo.RealMatrix(1)
let s = lo.RealMatrix(1)
let w = lo.RealMatrix(1)
let c = lo.ColorMatrix(1)
let dc = lo.ColorMatrix(1)
threads:
  for e in a:
    for i in 0..<8:
      for j in 0..<8:
        a[e][i,j] := (if i == j: (if i < 2: -2.0-float(i) else: 1.5+0.1*float(i)) else: 0.006*float(2*i-j+1))
        b[e][i,j] := 0.013*float(i+3*j-7)
        d[e][i,j] := 0.02*float(3*i-j+2)
  for e in r:
    # Bound finite-difference directions independently of volume and SIMD width.
    let z = float(e mod 4)
    r[e][0,0] := 1.2 + 0.01*z
    s[e][0,0] := 0.2 + 0.03*z
    w[e][0,0] := (if e mod 2 == 0: -0.7 else: 0.4) + 0.02*z
    c[e][0,0].re := 1.1+0.02*z
    c[e][0,0].im := -0.2+0.01*z
    dc[e][0,0].re := 0.1
    dc[e][0,0].im := 0.3

let ga = grt.toGvalue(a)
let gb = grt.toGvalue(b)
let gd = grt.toGvalue(d)
let gr = grt.toGvalue(r)
let gs = grt.toGvalue(s)
let gw = grt.toGvalue(w)
let gc = grt.toGvalue(c)
let gdc = grt.toGvalue(dc)

suite "graph real matrix fields":
  test "allocation copying and updates preserve concrete types":
    static:
      doAssert not compiles(grt.toGvalue(lo.RealMatrix(2)))
    check Gvalue(ga).newOneOf of Grmat8
    check Gvalue(gr).newOneOf of Grfield
    let copy = grt.toGvalue(a)
    check norm2(copy-ga).eval.sval < 1e-28
    copy.update b
    check norm2(copy-gb).eval.sval < 1e-28

  test "local reductions and global sums":
    # Compare densities: extensive reduction roundoff grows with lattice volume.
    check abs((sum(trace(ga))-retr(ga)).eval.sval)/vol < 1e-12
    check abs((sum(siteNorm2(ga))-norm2(ga)).eval.sval)/vol < 1e-12
    check abs((sum(siteRedot(ga,gb))-redot(ga,gb)).eval.sval)/vol < 1e-12
    check norm2(grad(3.0*sum(gr),gr)-3.0*unitField(gr.runtime,gr.fval)).eval.sval < 1e-26
    # The implicit root seed of a matrix field is the site identity, not all ones.
    check norm2(grad(ga,ga)-unitField(ga.runtime,ga.fval)).eval.sval < 1e-26
    check abs((sum(maskSubset(0,unitField(gr.runtime,gr.fval)))-0.5*float(lo.physVol)).eval.sval) < 1e-12

  test "1x1 solve inverse and logdet values and weighted pullbacks":
    check norm2(solve(gr,gs)-gs/gr).eval.sval < 1e-28
    check norm2(inverse(gr)-1.0/gr).eval.sval < 1e-28
    check norm2(logDet(gr)-ln(gr)).eval.sval < 1e-28
    let f = redot(solve(gr,gs),gw)
    check norm2(grad(f,gr)+gw*gs/(gr*gr)).eval.sval < 1e-26
    check norm2(grad(f,gs)-gw/gr).eval.sval < 1e-26
    check norm2(grad(redot(inverse(gr),gw),gr)+gw/(gr*gr)).eval.sval < 1e-26
    check norm2(grad(redot(logDet(gr),gw),gr)-gw/gr).eval.sval < 1e-26

  test "general 8x8 solve and inverse values":
    check norm2(ga*solve(ga,gb)-gb).eval.sval < 1e-24
    check norm2(ga*inverse(ga)-unitField(ga.runtime,ga.fval)).eval.sval < 1e-24
    check norm2(transpose(transpose(ga))-ga).eval.sval < 1e-28

  test "solve pullbacks include both operands and arbitrary upstreams":
    proc f(x,y: Grmat8): Gscalar = redot(solve(x,y),gd)
    ckmat2(f, ga, gb, gd, ga)
    proc f2(x: Grmat8): Gscalar = redot(grad(f(x,gb),x),x.transpose*gd+gb)
    ckmat(f2,ga,gd)

  test "solve repeated operand pullbacks cancel":
    let y = solve(ga,ga)
    check norm2(y-unitField(ga.runtime,ga.fval)).eval.sval < 1e-24
    check norm2(grad(redot(y,gd),ga)).eval.sval < 1e-24

  test "weighted logdet differentiates through dependent cotangents":
    proc f(x: Grmat8): Gscalar = redot(logDet(x),gs)
    ckmat(f,ga,gd)
    proc constant(x: Grmat8): Gscalar = redot(grad(f(x),x),x*gd)
    check abs((constant(ga)-redot(gs,trace(gd))).eval.sval)/vol < 1e-12
    check norm2(grad(constant(ga),ga)).eval.sval < 1e-22
    proc f2(x: Grmat8): Gscalar = redot(grad(f(x),x),x.transpose*gd+gb)
    ckmat(f2,ga,gb)
    proc f3(x: Grmat8): Gscalar = redot(grad(f2(x),x),x*gb+gd)
    ckmat(f3,ga,gd)
    check abs((sumLogDet(ga)-sum(logDet(ga))).eval.sval)/vol < 1e-12

  test "real scalar functions and division":
    proc f(x: Grfield): Gscalar = sum(exp(sin(x))+ln(x)+sqrt(x)+cos(x)+x/gs)
    ckgrad(f,gr,gs)
    proc f2(x: Grfield): Gscalar = redot(grad(f(x),x),x*gs)
    ckgrad(f2,gr,gs)

  test "real and imaginary embeddings carry arbitrary complex cotangents":
    check norm2(re(complex(gr)+imaginary(gs))-gr).eval.sval < 1e-26
    check norm2(im(complex(gr)+imaginary(gs))-gs).eval.sval < 1e-26
    proc f(x: Grfield): Gscalar = redot(complex(x)+imaginary(x*x),gc)
    ckgrad(f,gr,gs)
    proc f2(x: Grfield): Gscalar = redot(grad(f(x),x),x)
    ckgrad(f2,gr,gs)
    proc fc(x: Gcfield): Gscalar = sum(re(x)*gr+im(x)*gs)
    ckgrad(fc,gc,gdc)

  test "real complex scaling and angle functions close under differentiation":
    proc f(x: Grfield): Gscalar = norm2(scale(x,gc)*expi(x))+redot(expi(x),gc)
    ckgrad(f,gr,gs)
    proc f2(x: Grfield): Gscalar = redot(grad(f(x),x),x)
    ckgrad(f2,gr,gs)
    proc fc(x: Gcfield): Gscalar = sum(gs*arg(x))+norm2(scale(gr,x))
    ckgrad(fc,gc,gdc)
    proc fc2(x: Gcfield): Gscalar = redot(grad(fc(x),x),complex(gs)*x)
    ckgrad(fc2,gc,gdc)

  when DColorMatrixV is ColorMatrixN[3,DComplexV]:
    test "SU3 bridge graph adjoints":
      let m = lo.ColorMatrixD
      let dm = lo.ColorMatrixD
      threads:
        for e in m:
          for i in 0..<3:
            for j in 0..<3:
              m[e][i,j].re := 0.02*float(i+2*j-1)
              m[e][i,j].im := 0.03*float(2*i-j+1)
              dm[e][i,j].re := 0.017*float(2*i-3*j+2)
              dm[e][i,j].im := -0.011*float(i+j+1)
      let gm = grt.toGvalue(m)
      let gdm = grt.toGvalue(dm)
      check abs((redot(su3AdNeg(gm),gd)-redot(gm,su3AdNegAdj(gd))).eval.sval) < 1e-12
      check abs((redot(su3ProjectDeriv(gm),gd)-redot(gm,su3ProjectDerivAdj(gd))).eval.sval) < 1e-12
      proc f(x: Gfield): Gscalar = norm2(su3AdNeg(x)*su3ProjectDeriv(x))
      ckgrad(f,gm,gm)
      ckmat(f,gm,gdm)

qexFinalize()
