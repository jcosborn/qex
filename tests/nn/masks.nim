import qex, nn
import std/math
import common

proc testMasks*[T: SomeFloat]() =
  suite "NN field masks " & $T:
    let lo = newLayout(@[8,16])

    test "global coordinates select SIMD lanes and subset complements stay unchanged":
      let x = fields[T](lo,2)
      let y = fields[T](lo,2,T(-91))
      let ms = realField(lo,float32)
      for s in 0..<lo.nSites:
        ms{s} := (if lo.coords[0][s] mod 2 == 0: -1'f32 else: 0'f32)
        for c in 0..1: x[c]{s} := T(1+c)+T(lo.coords[1][s])/T(4)
      maskedCopy(y,x,ms,"odd")
      for s in 0..<lo.nSites:
        let active = lo.coords[0][s] mod 2 == 0 and (lo.coords[0][s]+lo.coords[1][s]) mod 2 == 1
        for c in 0..1:
          check y[c].sample(s) == (if active: x[c].sample(s) else: T(-91))
      scale(y,y,@[T(2),T(3)],"odd",ms)
      bias(y,y,@[T(7),T(-5)],"odd",ms)
      for s in 0..<lo.nSites:
        let active = lo.coords[0][s] mod 2 == 0 and (lo.coords[0][s]+lo.coords[1][s]) mod 2 == 1
        for c in 0..1:
          let expected = if c == 0: T(2)*x[c].sample(s)+T(7) else: T(3)*x[c].sample(s)-T(5)
          check close(y[c].sample(s),(if active: expected else: T(-91)))
      for s in 0..<lo.nSites: ms{s} := 0'f32
      gelu(y,x,"all",ms)
      for s in 0..<lo.nSites:
        let active = lo.coords[0][s] mod 2 == 0 and (lo.coords[0][s]+lo.coords[1][s]) mod 2 == 1
        for c in 0..1:
          let expected = if c == 0: T(2)*x[c].sample(s)+T(7) else: T(3)*x[c].sample(s)-T(5)
          check close(y[c].sample(s),(if active: expected else: T(-91)))

    test "mask pullbacks partition the seed and preserve the requested passthrough":
      let x = fields[T](lo,1,T(0))
      let dy = fields[T](lo,1,T(2))
      let dx = fields[T](lo,1,T(91))
      let df = fields[T](lo,1,T(-37))
      let ms = realField(lo,float32)
      for s in 0..<lo.nSites:
        ms{s} := (if lo.coords[1][s] mod 3 == 1: 1'f32 else: 0'f32)
      maskVjp(dx,dy,"even",ms)
      maskVjp(df,dy,"even",ms,complement=true)
      for s in 0..<lo.nSites:
        let active = lo.coords[1][s] mod 3 == 1 and (lo.coords[0][s]+lo.coords[1][s]) mod 2 == 0
        check dx[0].sample(s) == (if active: T(2) else: T(0))
        check df[0].sample(s) == (if active: T(0) else: T(2))
        check dx[0].sample(s)+df[0].sample(s) == dy[0].sample(s)
      geluVjp(dx,x,dy,"even",ms,passthrough=true)
      for s in 0..<lo.nSites:
        let active = lo.coords[1][s] mod 3 == 1 and (lo.coords[0][s]+lo.coords[1][s]) mod 2 == 0
        check dx[0].sample(s) == (if active: T(1) else: T(2))

    test "masked convolution adjoints compose the kernel with maskVjp":
      # y = select_A(conv(x)) has the adjoint dx = convVjp(select_A(dy)); A is the
      # odd sites of even rows and the asymmetric taps cross both boundaries.
      proc value(r,c: int): T = T(0.3+0.05*float(r)-0.02*float(c))
      proc seed(r,c: int): T = T(0.7-0.04*float(r)+0.03*float(c))
      let x = fields[T](lo,1)
      let y = fields[T](lo,1)
      let dy = fields[T](lo,1)
      let sdy = fields[T](lo,1)
      let dx = fields[T](lo,1,T(91))
      let ms = realField(lo,float32)
      for s in 0..<lo.nSites:
        let r = lo.coords[0][s].int
        let c = lo.coords[1][s].int
        ms{s} := (if r mod 2 == 0: 1'f32 else: 0'f32)
        x[0]{s} := value(r,c)
        dy[0]{s} := seed(r,c)
      var w = newSeq[T](15)
      for t in 0..<15: w[t] = T(0.1*float(t)-0.6)
      let p = convParams(1,1,[3,5],w)
      let ws = convWorkspace(x[0],p)
      conv(y,x,p,ws)
      maskVjp(y,y,"odd",ms)
      maskVjp(sdy,dy,"odd",ms)
      convVjp(dx,sdy,p,ws)
      let tol = when T is float32: 2e-5 else: 2e-12
      let pair = dot(y,dy)
      check abs(pair-dot(x,dx)) <= tol*max(1.0,abs(pair))
      for s in 0..<lo.nSites:
        # dx(u) = sum_t w_t dy(u-offset_t) over selected output sites u-offset_t.
        var want = 0.0
        for t in 0..<15:
          let r = (lo.coords[0][s].int-(t div 5-1)+8) mod 8
          let c = (lo.coords[1][s].int-(t mod 5-2)+16) mod 16
          if r mod 2 == 0 and (r+c) mod 2 == 1: want += float(w[t])*float(seed(r,c))
        check abs(float(dx[0].sample(s))-want) <= tol*max(1.0,abs(want))
      for s in 0..<lo.nSites: ms{s} := 0'f32
      maskVjp(sdy,dy,"odd",ms)
      convVjp(dx,sdy,p,ws)
      check dot(dx,dx) == 0
