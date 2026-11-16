import math
import testutils
import qex

# The fused expAH kernels at the build's SIMD width, VLEN lanes, and at 2 VLEN
# lanes when that SIMD type exists. expAH evaluates entries wider than one SIMD
# register (VLEN lanes without intrinsics) in slices of one register with the
# scale of the whole site, which rounds differently from the expression for one
# lane: compare to 8 eps in norm.

qexInit()
letParam:
  lat = latticeFromLocalLattice(@[4,4,4,4], nRanks)
let tol = (8.0*epsilon(float))^2

proc laneRelNorm2(z, x: Field; a: float; y: Field = nil): float =
  ## relNorm2(z, w) for w = expAH(a x) y, y = 1 when nil, with each lane of w
  ## from the scalar matrices of that lane: their own scale and no slices.
  const n = x[0].nrows
  const complex = evalType(x[0][0,0]) is ComplexProxy
  var s = [0.0, 0.0]
  for e in 0..<x.l.nSitesOuter:
    for k in 0..<x.l.V:
      when complex:
        var m, ym {.noinit.}: MatrixArray[n, n, ComplexType[float64]]
      else:
        var m, ym {.noinit.}: MatrixArray[n, n, float64]
      for i in 0..<n:
        for j in 0..<n:
          when complex:
            m[i,j].re := x[e][i,j].re[k]
            m[i,j].im := x[e][i,j].im[k]
            if y != nil:
              ym[i,j].re := y[e][i,j].re[k]
              ym[i,j].im := y[e][i,j].im[k]
          else:
            m[i,j] := x[e][i,j][k]
            if y != nil: ym[i,j] := y[e][i,j][k]
      var w = expAH(a*m)
      if y != nil: w = w*ym
      for i in 0..<n:
        for j in 0..<n:
          when complex:
            s[0] += (z[e][i,j].re[k] - w[i,j].re)^2 + (z[e][i,j].im[k] - w[i,j].im)^2
            s[1] += w[i,j].re^2 + w[i,j].im^2
          else:
            s[0] += (z[e][i,j][k] - w[i,j])^2
            s[1] += w[i,j]^2
  rankSum(s)
  s[0]/(1.0+s[1])

proc adjointViews(x: Field) =
  let e = x.newOneOf
  let re = x.newOneOf
  const kind = if evalType(x[0][0,0]) is ComplexProxy: "complex" else: "real"
  suite "expAH " & kind & " adjoint views at " & $x.l.V & " lanes":
    test "views agree with materialized adjoints":
      for a0 in [0.0, 0.0153, 0.125, -0.75, 3.0]:
        let a = a0
        threads:
          for i in x:
            var m, ma {.noinit.}: evalType(x[i][])
            m := a*x[i][]
            ma := m.adj
            re[i][] := expAH(ma)
            e[i][] := expAH(m.adj)
        check relNorm2(e, re) <= tol
        check laneRelNorm2(e, x, -a) <= tol

proc thresholds[V: static int](lo: Layout[V], n: static int) =
  let ah = lo.ColorMatrix(n)
  let src = lo.ColorMatrix(n)
  let other = lo.ColorMatrix(n)
  let e = lo.ColorMatrix(n)
  let z = lo.ColorMatrix(n)
  let re = lo.ColorMatrix(n)
  let rz = lo.ColorMatrix(n)
  let saved = lo.ColorMatrix(n)
  threads:
    for k in ah:
      ah[k] := 0.0
      var imag, neg: array[V, float64]
      for lane in 0..<imag.len:
        let factor = [0.5, 1.0-2e-14, 1.0+2e-14, 4.0][lane mod 4]
        imag[lane] = factor*0.1767766952966369
        neg[lane] = -imag[lane]
      ah[k][0,0].im := imag
      when n > 1: ah[k][1,1].im := neg
      for i in 0..<n:
        for j in 0..<n:
          other[k][i,j].re := 0.11*float(i-j)+(if i == j: 0.2 else: 0.0)
          other[k][i,j].im := 0.03*float(i+2*j+1)
  suite "native exponential thresholds Nc=" & $n & " at " & $V & " lanes":
    test "mixed SIMD lanes preserve algorithms adjoints and selected outputs":
      for (gen, adj) in [(false,false), (false,true), (true,false), (true,true)]:
        var general = gen
        var flip = adj
        for alpha in [0.5, 1.0, -2.0]:
          var coefficient = alpha
          threads:
            if flip: src := ah.adj
            else: src := ah
            saved := src
            for k in src:
              var arg {.noinit.}: evalType(src[k])
              arg := coefficient*src[k]
              if general: re[k][] := exp(arg[])
              else: re[k][] := expAH(arg[])
            rz := re*other
            e := 9.0
            z := -7.0
            if general: axexp(e, coefficient, src)
            else: axexpah(e, coefficient, src)
          check relNorm2(e, re) <= tol
          if not general: check laneRelNorm2(e, src, coefficient) <= tol
          threads:
            if general: axexpmuly(z, coefficient, src, other)
            else: axexpahmuly(z, coefficient, src, other)
          check relNorm2(z, rz) <= tol
          if not general: check laneRelNorm2(z, src, coefficient, other) <= tol
          threads:
            e := 9.0
            z := -7.0
            if general: axexpmuly(z, coefficient, src, other, e)
            else: axexpahmuly(z, coefficient, src, other, e)
          check relNorm2(e, re) <= tol
          check relNorm2(z, rz) <= tol
          check relNorm2(src, saved) == 0.0

proc realMatrices[V: static int](lo: Layout[V], n: static int) =
  ## A real anti-Hermitian matrix is antisymmetric; expAH of 1x1 is 1.
  ## e and z share one allocation with `next`, which must stay unchanged.
  let x = lo.RealMatrix(n)
  let y = lo.RealMatrix(n)
  let outs = newFieldArray(lo, typeof(x), 3)
  let e = outs[0]
  let z = outs[1]
  let next = outs[2]
  threads:
    for k in x:
      for i in 0..<n:
        for j in 0..<n:
          var v, w: array[V, float64]
          for lane in 0..<V:
            v[lane] = 0.37*sin(float(97*k+13*lane+5*i+3*j+1))
            w[lane] = 0.21*cos(float(31*k+7*lane+i-2*j))
          y[k][i,j] := w
          if i < j: x[k][i,j] := v
      for i in 0..<n:
        x[k][i,i] := 0.0
        for j in 0..<i: x[k][i,j] := -x[k][j,i]
    next := 0.25
  when n > 1: adjointViews(x)
  suite "real matrices Nc=" & $n & " at " & $V & " lanes":
    test "expAH of real antisymmetric sites":
      for a0 in [0.0, 0.0153, 0.125, -0.75, 3.0]:
        let a = a0  # threads cannot capture the lent loop value
        threads:
          axexpah(e, a, x)
          axexpahmuly(z, a, x, y)
        check laneRelNorm2(e, x, a) <= tol
        check laneRelNorm2(z, x, a, y) <= tol
        threads:
          e := 0
          axexpahmuly(z, a, x, y, e)
        check laneRelNorm2(e, x, a) <= tol
        check laneRelNorm2(z, x, a, y) <= tol
        var n2: float
        threads:
          let t = next.norm2
          threadSingle: n2 = t
        check n2 == 0.0625*float(n*lo.physVol)  # next = 1/4 on the diagonal

proc run[V: static int](lo: Layout[V]) =
  let x = lo.ColorMatrix(3)
  let y = lo.ColorMatrix(3)
  let e = lo.ColorMatrix(3)
  let z = lo.ColorMatrix(3)
  let re = lo.ColorMatrix(3)
  let rz = lo.ColorMatrix(3)
  let tmp = lo.ColorMatrix(3)
  var rng = lo.newRNGField(Philox4x64, 20261005'u64)
  threads:
    x.randomTAH rng
    y.gaussian rng
  adjointViews(x)
  suite "expAH kernels at " & $V & " lanes":
    test "whole fields with and without the stored exponential":
      # 0.125 puts |a x|^2 near 1/16, where lanes would choose different scales.
      for a0 in [0.0, 0.0153, 0.125, -0.75, 3.0]:
        let a = a0  # threads cannot capture the lent loop value
        threads:
          for i in x:
            re[i][] := expAH(a * x[i][])
            rz[i] := re[i] * y[i]
          axexpah(e, a, x)
          axexpahmuly(z, a, x, y)
        check relNorm2(e, re) <= tol
        check relNorm2(z, rz) <= tol
        check laneRelNorm2(e, x, a) <= tol
        check laneRelNorm2(z, x, a, y) <= tol
        threads:
          e := 0
          axexpahmuly(z, a, x, y, e)
        check relNorm2(e, re) <= tol
        check relNorm2(z, rz) <= tol

    test "aliased operands and subsets":
      const a = 0.375
      threads:
        for i in x:
          re[i][] := expAH(a * x[i][])
          rz[i] := re[i] * y[i]
        tmp := x
        axexpah(tmp, a, tmp)
      check relNorm2(tmp, re) <= tol
      threads:
        tmp := y
        axexpahmuly(tmp, a, x, tmp)
      check relNorm2(tmp, rz) <= tol
      for par in ["even", "odd"]:
        let sub = lo.getSubset(par)
        threads:
          z := -0.875
          tmp := -0.875
        threads:
          for i in sub: tmp[i] := rz[i]
          axexpahmuly(z, a, x, y, sub)
        check relNorm2(z, tmp) <= tol

  thresholds(lo, 1)
  thresholds(lo, 3)
  realMatrices(lo, 1)
  realMatrices(lo, 3)

run(lat.newLayout)
const wide = 2*VLEN
when wide <= 16:
  run(lat.newLayout(wide))

qexFinalize()
