#RUNCMD $RUN1
import std/[unittest, math, complex, strformat]
import base/alignedMem
import ../core/spinor

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))

proc mkTest(n: int, off: float): Spin =
  ## Exactly representable entries: x_i = (i+off) + i*(2i+off)/4.
  result = newSpin(n)
  for i in 0..<n:
    for c in 0..1:
      result[i][c].re = float(4*i + c) + off
      result[i][c].im = 0.25*float(2*i - c) - off

suite "spinor":

  test "dot symmetry and norms":
    let x = mkTest(37, 0.5)
    let y = mkTest(37, -1.25)
    let dxy = dot(x, y)
    let dyx = dot(y, x)
    # -ffast-math lets each proc reassociate its own reduction, so compare at
    # the roundoff level of the vector norms rather than bit for bit.
    let sc = 1e-13*sqrt(norm2(x)*norm2(y))
    check abs(dxy.re - dyx.re) <= sc
    check abs(dxy.im + dyx.im) <= sc
    check abs(redot(x, y) - dxy.re) <= sc
    let dxx = dot(x, x)
    check dxx.im == 0.0           ## xr*xi - xi*xr cancels exactly, term by term
    check abs(norm2(x) - dxx.re) <= 1e-13*norm2(x)
    check abs(norm2(x) - redot(x, x)) <= 1e-13*norm2(x)

  test "axpy and scale are exact":
    var x = mkTest(11, 0.5)
    let y = mkTest(11, -1.25)
    let x0 = x
    axpy(x, 0.5, y)
    for i in 0..<x.len:
      for c in 0..1:
        check x[i][c].re == x0[i][c].re + 0.5*y[i][c].re
        check x[i][c].im == x0[i][c].im + 0.5*y[i][c].im
    # a*y with a = (0, 1) is an exact rotation of the components
    var z = newSpin(11)
    axpy(z, complex64(0.0, 1.0), y)
    for i in 0..<z.len:
      for c in 0..1:
        check z[i][c].re == -y[i][c].im
        check z[i][c].im == y[i][c].re
    var s = mkTest(11, 0.5)
    scale(s, 0.25)
    for i in 0..<s.len:
      for c in 0..1:
        check s[i][c].re == 0.25*x0[i][c].re
        check s[i][c].im == 0.25*x0[i][c].im
    # axpby: x = a*y + b*x
    var p = mkTest(11, 0.5)
    axpby(p, 2.0, y, 0.5)
    for i in 0..<p.len:
      for c in 0..1:
        check p[i][c].re == 2.0*y[i][c].re + 0.5*x0[i][c].re

  test "zero and :=":
    var x = mkTest(9, 0.5)
    x.zero
    check norm2(x) == 0.0
    let y = mkTest(9, 0.5)
    x := y
    check norm2(x) == norm2(y)
    for i in 0..<x.len:
      check x[i][0] == y[i][0]
      check x[i][1] == y[i][1]

  test "gaussian normalization and stream independence":
    const n = 20000
    var r: Threefry4x64
    r.seedIndep(987654321, 0)
    var x = newSpin(n)
    x.gaussian r
    let m = norm2(x)/float(2*n)   ## mean |x_i|^2 per complex component
    echo "  gaussian <|x|^2> = ", m, "  (1 +- ", 1.0/sqrt(float(2*n)), ")"
    check abs(m - 1.0) < 0.03     ## 6 sigma

    var r1: Threefry4x64
    r1.seedIndep(987654321, 1)
    var y = newSpin(n)
    y.gaussian r1
    let corr = redot(x, y)/sqrt(norm2(x)*norm2(y))
    echo "  stream 0 vs 1 correlation = ", corr, "  (0 +- ", 0.5/sqrt(float(n)), ")"
    check abs(corr) < 0.02        ## ~6 sigma
    var d = x
    axpy(d, -1.0, y)
    check norm2(d) > 0.5*(norm2(x) + norm2(y))

    # same seed and index reproduces the stream bit for bit
    var r2: Threefry4x64
    r2.seedIndep(987654321, 0)
    var x2 = newSpin(n)
    x2.gaussian r2
    var e = x
    axpy(e, -1.0, x2)
    check norm2(e) == 0.0

  test "pointSource":
    let n = 40
    let x = pointSource(n, 17, 1)
    check norm2(x) == 1.0
    var nz = 0
    for i in 0..<n:
      for c in 0..1:
        if x[i][c].re != 0.0 or x[i][c].im != 0.0: inc nz
    check nz == 1
    check x[17][1].re == 1.0

  test "allocation regression":
    const n = 512
    var r: Threefry4x64
    r.seedIndep(11, 0)
    var x = newSpin(n)
    var y = newSpin(n)
    x.gaussian r
    y.gaussian r
    # warm up so any lazily grown internal buffer is already there
    axpy(x, 1e-3, y)
    discard dot(x, y)
    discard norm2(x)
    axpby(x, 1e-3, y, 1.0)
    let raw0 = getRawMemAllocated()
    let occ0 = getOccupiedMem()
    var acc = 0.0
    for k in 0..<64:
      axpy(x, 1e-9, y)
      scale(x, 1.0 - 1e-9)
      axpby(y, 1e-9, x, 1.0)
      acc += dot(x, y).re + norm2(x) + redot(x, y)
    let raw1 = getRawMemAllocated()
    let occ1 = getOccupiedMem()
    # getRawMemAllocated only tracks QEX alignedMem; a Spin is a plain seq, so
    # it stays 0 here and getOccupiedMem is the check with teeth.
    echo &"  rawMem {raw0} -> {raw1}   occupied {occ0} -> {occ1}"
    check raw1 == raw0
    check occ1 == occ0
    check acc != 0.0    ## keep the loop alive
