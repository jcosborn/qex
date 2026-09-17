import qex, nn
import std/math
import common

proc testPointwise*[T: SomeFloat]() =
  suite "NN pointwise " & $T:
    let lo = newLayout(@[8,12])

    test "GELU value, zero derivative, and reflection identity":
      let x = fields[T](lo,2)
      let neg = fields[T](lo,2)
      let a = fields[T](lo,2)
      let b = fields[T](lo,2)
      let seed = fields[T](lo,2,T(1))
      let dx = fields[T](lo,2)
      for c in 0..1:
        for s in 0..<lo.nSites:
          let v = T(c) + T(lo.coords[1][s]-6)/T(4)
          x[c]{s} := v
          neg[c]{s} := -v
      gelu(a,x)
      gelu(b,neg)
      geluVjp(dx,x,seed)
      for c in 0..1:
        for s in 0..<lo.nSites:
          let v = x[c].sample(s)
          let expected = T(0.5)*v*(T(1)+erf(v/sqrt(T(2))))
          check close(a[c].sample(s),expected)
          check close(a[c].sample(s)-b[c].sample(s),v)
          if v == 0:
            check a[c].sample(s) == 0
            check dx[c].sample(s) == T(0.5)
      gelu(x,x)
      check sameFields(x,a)

    test "exponential, logarithm, erfc, clipping and division kernels":
      let x = fields[T](lo,2)
      let y = fields[T](lo,2)
      let z = fields[T](lo,2)
      for c in 0..1:
        for s in 0..<lo.nSites:
          x[c]{s} := T(0.25)+T(c)/T(3)+T(lo.coords[0][s])/T(16)
          y[c]{s} := T(1.5)-T(lo.coords[1][s])/T(32)
      exp(z,x)
      for c in 0..1:
        for s in 0..<lo.nSites: check close(z[c].sample(s),exp(x[c].sample(s)))
      ln(z,z)
      check sameFields(z,x)
      erfc(z,x)
      for c in 0..1:
        for s in 0..<lo.nSites: check close(z[c].sample(s),erfc(x[c].sample(s)))
      divide(z,x,y)
      for c in 0..1:
        for s in 0..<lo.nSites: check close(z[c].sample(s),x[c].sample(s)/y[c].sample(s))
      let floor = T(0.5)
      clipMin(z,x,floor)
      clipSlope(y,x,floor)
      for c in 0..1:
        for s in 0..<lo.nSites:
          let v = x[c].sample(s)
          check z[c].sample(s) == max(v,floor)
          check y[c].sample(s) == (if v < floor: T(0) elif v == floor: T(0.5) else: T(1))
      check abs(redot(x,z)-dot(x,z)) <= 1e-12*abs(dot(x,z))

    test "GELU and arctan pullbacks agree with independent finite differences":
      let x = fields[T](lo,1,T(0.37))
      let dy = fields[T](lo,1,T(0.71))
      let dx = fields[T](lo,1)
      let v = T(0.37)
      let h = when T is float32: T(0.002) else: T(0.00001)
      let tol = when T is float32: T(8e-5) else: T(2e-10)
      proc g(t: T): T = T(0.5)*t*(T(1)+erf(t/sqrt(T(2))))
      geluVjp(dx,x,dy)
      let gd = T(0.71)*(g(v+h)-g(v-h))/(T(2)*h)
      for s in 0..<lo.nSites: check abs(dx[0].sample(s)-gd) < tol
      arctanVjp(dx,x,dy)
      let ad = T(0.71)*(math.arctan(v+h)-math.arctan(v-h))/(T(2)*h)
      for s in 0..<lo.nSites: check abs(dx[0].sample(s)-ad) < tol
      arctanVjp(dy,x,dy)
      check sameFields(dy,dx)

    test "arctan at zero and one, with in-place values and pullbacks":
      let x = fields[T](lo,2)
      x[1] := T(1)
      let dy = fields[T](lo,2,T(2))
      let dx = fields[T](lo,2)
      arctanVjp(dx,x,dy)
      arctan(x,x)
      for s in 0..<lo.nSites:
        check x[0].sample(s) == 0
        check close(x[1].sample(s),T(PI)/T(4))
        check dx[0].sample(s) == T(2)
        check dx[1].sample(s) == T(1)

    test "channel bias and scale preserve isolation and support in-place updates":
      let x = fields[T](lo,3)
      for c in 0..2:
        for s in 0..<lo.nSites: x[c]{s} := T(c+1)
      bias(x,x,@[T(3),T(-7),T(0.5)])
      scale(x,x,@[T(2),T(-3),T(0)])
      for s in 0..<lo.nSites:
        check x[0].sample(s) == T(8)
        check x[1].sample(s) == T(15)
        check x[2].sample(s) == T(0)
      let dy = fields[T](lo,3,T(1))
      scaleVjp(dy,dy,@[T(2),T(-3),T(0)])
      for s in 0..<lo.nSites:
        check dy[0].sample(s) == T(2)
        check dy[1].sample(s) == T(-3)
        check dy[2].sample(s) == T(0)

    test "division retains each precision boundary and the reverse order":
      let x = fields[T](lo,1)
      let a = fields[T](lo,1)
      let dy = fields[T](lo,1)
      for s in 0..<lo.nSites:
        x[0]{s} := T(0.19)+T(lo.coords[1][s])/T(7)
        dy[0]{s} := T(0.31)+T(lo.coords[0][s])/T(11)
      divide(a,x,T(7.3))
      divide(a,a,T(5.9))
      let dx = fields[T](lo,1)
      divideVjp(dx,dy,T(5.9))
      divideVjp(dx,dx,T(7.3))
      # Fast-math may divide through a reciprocal on either side: a few ulps of T.
      const eps = when T is float64: 2.220446049250313e-16 else: 1.1920928955078125e-7
      for s in 0..<lo.nSites:
        let q = float(x[0].sample(s)/T(7.3)/T(5.9))
        let r = float(dy[0].sample(s)/T(5.9)/T(7.3))
        check abs(float(a[0].sample(s))-q) <= 4*eps*abs(q)
        check abs(float(dx[0].sample(s))-r) <= 4*eps*abs(r)

    test "channel reduction and broadcast are global adjoints and overwrite outputs":
      for dims in [@[8,12],@[16,12]]:
        let layout = newLayout(dims)
        let x = fields[T](layout,3)
        let y = fields[T](layout,3,T(-73))
        for c in 0..<3:
          for s in layout.sites:
            x[c]{s} := T(c+1)+T(layout.coords[0][s])/T(4)+T(layout.coords[1][s])/T(8)
        var sums = newSeq[T](3)
        channelSum(sums,x)
        let values = @[T(0.5),T(-0.25),T(2)]
        broadcast(y,values)
        var pair = 0.0
        for c in 0..<3:
          let want = T(layout.physVol)*(T(c+1)+T(dims[0]-1)/T(8)+T(dims[1]-1)/T(16))
          check sums[c] == want
          for s in layout.sites: check y[c].sample(s) == values[c]
          pair += float(sums[c])*float(values[c])
        check pair == dot(x,y)
        x.fill(T(0))
        channelSum(sums,x)
        for v in sums: check v == T(0)
        broadcast(y,sums)
        check dot(y,y) == 0
      var wrong = newSeq[T](2)
      expect ValueError: channelSum(wrong,fields[T](lo,3))
      expect ValueError: broadcast(fields[T](lo,3),wrong)
