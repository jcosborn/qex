import base/globals
setVLENmax(4)

import math
import testutils
import qex
import field/matrixFields
import maths/groupOps

qexInit()
letParam:
  lat = latticeFromLocalLattice(@[4,4], nRanks)
let lo = lat.newLayout
let a = lo.RealMatrix(8)
let b = lo.RealMatrix(8)
let y = lo.RealMatrix(8)
let t = lo.RealMatrix(8)
let r = lo.RealMatrix(1)
let s = lo.RealMatrix(1)
let c = lo.ColorMatrix(1)
let d = lo.ColorMatrix(1)
let old = lo.RealD

threads:
  for e in a:
    for i in 0..<8:
      for j in 0..<8:
        a[e][i,j] := (if i == j: (if i < 2: -2.0-float(i) else: 1.5+0.1*float(i)) else: 0.006*float(2*i-j+1))
        b[e][i,j] := 0.013*float(i+3*j-7)
  for e in r:
    r[e][0,0] := 1.1 + 0.01*float(e)
    s[e][0,0] := 0.3 + 0.02*float(e)

suite "numerical matrix fields":
  test "real scalar storage and explicit conversions":
    static:
      doAssert sizeof(DRealMatrixV[1]) * 64 == sizeof(DRealMatrixV[8])
      doAssert DLatticeRealMatrixV[1] isnot DLatticeRealV
    threads:
      toScalar(old, r)
      toMatrix(s, old)
    check diffNorm2(r, s) < 1e-28
    threads:
      complex(c, r)
      imaginary(d, s)
      c += d
      re(r, c)
      im(s, c)
    check diffNorm2(r, s) < 1e-28
    static:
      doAssert sizeof(evalType(re(c[0][])[0,0])) == sizeof(evalType(r[0][0,0]))

  test "LU helpers enforce square factors and right hand sides":
    static:
      doAssert not compiles(block:
        var a: MatrixArray[2,3,float64]
        a.luNoPivot)
      doAssert not compiles(block:
        var a: MatrixArray[2,2,float64]
        var b: MatrixArray[2,3,float64]
        a.solveLNoPivot(b))
      doAssert not compiles(block:
        var a: MatrixArray[2,2,float64]
        var b: MatrixArray[3,2,float64]
        a.solveRNoPivot(b))
    var aa, lu, bb, ll, rr: MatrixArray[2,2,float64]
    aa[0,0] = 3; aa[0,1] = 1; aa[1,0] = 2; aa[1,1] = 4
    bb[0,0] = 1; bb[0,1] = 2; bb[1,0] = 3; bb[1,1] = 5
    lu := aa
    ll := bb
    rr := bb
    lu.luNoPivot
    lu.solveLNoPivot(ll)
    lu.solveRNoPivot(rr)
    check norm2(aa*ll-bb) < 1e-28
    check norm2(rr*aa-bb) < 1e-28

  test "noncommuting solve inverse and transpose":
    threads:
      solve(y, a, b)
      t := a*y
    check diffNorm2(t, b) < 1e-25
    threads:
      inverse(y, a)
      t := a*y - 1.0
      let err = t.norm2
      threadSingle: check err < 1e-25
      transpose(y, b)
      transpose(t, y)
    check diffNorm2(t, b) < 1e-28

  test "positive determinant permits negative pivots":
    threads:
      a := 1.0
      for e in a:
        a[e][0,0] := -2.0
        a[e][1,1] := -3.0
      logDet(r, a)
      let v = sum(r)
      threadSingle: check abs(v - float(lo.physVol)*ln(6.0)) < 1e-12

  test "log determinant avoids products outside the floating point range":
    for v in [1e40, 1e-50]:
      threads:
        a := v
        for e in a:
          a[e][0,0] := -v
          a[e][1,1] := -v
        logDet(r, a)
        s := 8.0*ln(v)
      check diffNorm2(r, s)/float(lo.physVol) < 1e-22
    var m: MatrixArray[4,4,float64]
    m := 0.0
    for i, v in [1e200, 1e200, 1e-200, 1e-200]: m[i,i] = v
    check abs(logDet(m)) < 1e-12

  test "scalar sums have the same numerical type at every layout width":
    forStatic k, 0, 2:
      let l = newLayout(latticeFromLocalLattice(@[4,4], nRanks), 1 shl k)
      let x = l.RealMatrix(1)
      static: doAssert typeof(sum(x)) is float64
      threads:
        x := 2.0
        let v = sum(x)
        threadSingle: check v == 2.0*float(l.physVol)

  test "local contractions and masks retain physical normalization":
    threads:
      siteRedot(r, b, b)
      siteNorm2(s, b)
      let sr = sum(r)
      let sb = norm2(b)
      threadSingle: check abs(sr-sb) < 1e-12
    check diffNorm2(r, s) < 1e-28
    threads:
      r := 2.0
      maskSubset(s, r, lo.getSubset("even"))
      let v = sum(s)
      threadSingle: check abs(v-float(lo.physVol)) < 1e-12

  test "subset blends and masks permit either input to own the output":
    let z = lo.RealMatrix(1)
    for par in ["even", "odd"]:
      let sub = lo.getSubset(par)
      threads:
        r := 2.0
        s := 3.0
        for e in z:
          z[e][0,0] := (if e >= sub.lowOuter and e < sub.highOuter: 2.0 else: 3.0)
        blendSubset(r, r, s, sub)
      check diffNorm2(r, z) == 0.0
      threads:
        r := 2.0
        blendSubset(s, r, s, sub)
      check diffNorm2(s, z) == 0.0
      threads:
        for e in z:
          z[e][0,0] := (if e >= sub.lowOuter and e < sub.highOuter: 2.0 else: 0.0)
        maskSubset(r, r, sub)
      check diffNorm2(r, z) == 0.0

  test "scalar functions use SIMD entries":
    threads:
      r := 1.2
      exp(s, r)
      ln(r, s)
      s := 1.2
    check diffNorm2(r, s) < 1e-26
    threads:
      expi(c, r)
      arg(s, c)
    check diffNorm2(r, s) < 1e-26

  test "real scalar scaling of real and complex matrices":
    let m = lo.ColorMatrix(3)
    let n = lo.ColorMatrix(3)
    threads:
      r := 0.7
      scale(t, r, b)
      y := 0.7*b
      for e in m:
        for i in 0..<3:
          for j in 0..<3:
            m[e][i,j].re := 0.1*float(i-j)
            m[e][i,j].im := 0.2*float(i+j)
      scale(n, r, m)
      m := 0.7*m
    check diffNorm2(t, y) < 1e-28
    check diffNorm2(n, m) < 1e-28

  test "SU3 bridge adjoints accept general real matrices":
    let m = lo.ColorMatrix(3)
    let g = lo.ColorMatrix(3)
    threads:
      for e in m:
        for i in 0..<3:
          for j in 0..<3:
            m[e][i,j].re := 0.03*float(i+2*j-1)
            m[e][i,j].im := 0.04*float(2*i-j+1)
      su3AdNeg(a, m)
      su3AdNegAdj(g, b)
      let xa = redot(a,b)
      let xb = redot(m,g)
      threadSingle: check abs(xa-xb) < 1e-12
      su3ProjectDeriv(a, m)
      su3ProjectDerivAdj(g, b)
      let da = redot(a,b)
      let db = redot(m,g)
      threadSingle: check abs(da-db) < 1e-12


suite "complex scalar matrix scaling":
  test "multiplication by i preserves norm and rotates the inner product":
    let u=lo.ColorMatrix(3)
    let v=lo.ColorMatrix(3)
    let c=lo.ColorMatrix(1)
    var nr,dr,di:float
    threads:
      for e in u:
        c[e][0,0] := newComplex(0.0,1.0)
        for i in 0..<3:
          for j in 0..<3: u[e][i,j] := newComplex(float(i+2*j+1),float(3*i-j))
      scale(v,c,u)
      let n=u.norm2
      let d=dot(u,v)
      threadSingle:
        nr=n
        dr=d.re
        di=d.im
    check abs(dr)<1e-14*(1+nr)
    check abs(di-nr)<1e-14*(1+nr)
    threads:
      let n=v.norm2
      threadSingle: dr=n
    check abs(dr-nr)<1e-14*(1+nr)

qexFinalize()
