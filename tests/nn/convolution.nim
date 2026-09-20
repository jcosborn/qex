import qex, nn
import std/[math, sequtils]
import common

proc testConvolution*[T: SomeFloat]() =
  suite "NN convolution " & $T:
    let lo = newLayout(@[8,16])

    test "one-tap identity, channel mixing, and bias ordering":
      let x = fields[T](lo,2)
      let y = fields[T](lo,2)
      for s in 0..<lo.nSites:
        x[0]{s} := T(lo.coords[0][s]+1)
        x[1]{s} := T(lo.coords[1][s]-3)/T(2)
      let id = convParams(2,2,[1,1],@[T(1),T(0),T(0),T(1)])
      let ws = convWorkspace(x[0],id)
      conv(y,x,id,ws)
      check sameFields(x,y)
      let p = convParams(2,2,[1,1],@[T(2),T(-3),T(5),T(7)],@[T(0.5),T(-1)])
      conv(y,x,p,ws)
      for s in 0..<lo.nSites:
        let a = x[0].sample(s)
        let b = x[1].sample(s)
        check close(y[0].sample(s),T(2)*a-T(3)*b+T(0.5))
        check close(y[1].sample(s),T(5)*a+T(7)*b-T(1))

    test "asymmetric taps wrap a boundary impulse on a rectangular lattice":
      let x = fields[T](lo,1)
      let y = fields[T](lo,1)
      for s in 0..<lo.nSites:
        if lo.coords[0][s] == 0 and lo.coords[1][s] == 15: x[0]{s} := T(1)
      let p = convParams(1,1,[3,3],@[T(1),T(2),T(4),T(8),T(16),T(32),T(64),T(128),T(256)])
      let ws = convWorkspace(x[0],p)
      conv(y,x,p,ws)
      for s in 0..<lo.nSites:
        let r = lo.coords[0][s]
        let c = lo.coords[1][s]
        var expected = T(0)
        if r == 1:
          if c == 0: expected = T(1)
          elif c == 15: expected = T(2)
          elif c == 14: expected = T(4)
        elif r == 0:
          if c == 0: expected = T(8)
          elif c == 15: expected = T(16)
          elif c == 14: expected = T(32)
        elif r == 7:
          if c == 0: expected = T(64)
          elif c == 15: expected = T(128)
          elif c == 14: expected = T(256)
        check y[0].sample(s) == expected

    test "adjoint identity and repeated reverse calls clear all accumulated storage":
      let x = fields[T](lo,2)
      let y = fields[T](lo,3)
      let dy = fields[T](lo,3)
      let dx = fields[T](lo,2,T(91))
      let saved = fields[T](lo,2)
      for c in 0..<x.len:
        for s in 0..<lo.nSites:
          x[c]{s} := T(0.2+0.03*float(c)+0.07*float(lo.coords[0][s])-0.05*float(lo.coords[1][s]))
      for c in 0..<dy.len:
        for s in 0..<lo.nSites:
          dy[c]{s} := T(0.31*sin(float(c)+0.11*float(lo.coords[0][s])+0.17*float(lo.coords[1][s])))
      let p = convParams(2,3,[3,3],(0..<54).toSeq.mapIt(T(0.03*sin(float(it)+0.2))))
      let ws = convWorkspace(x[0],p)
      conv(y,x,p,ws)
      convVjp(dx,dy,p,ws)
      let lhs = dot(x,dx)
      let rhs = dot(y,dy)
      let tol = when T is float32: 2e-5 else: 2e-12
      check abs(lhs-rhs) <= tol*max(1.0,max(abs(lhs),abs(rhs)))
      threads:
        for c in 0..<dx.len: saved[c] := dx[c]
      dx.fill(T(-73))
      convVjp(dx,dy,p,ws)
      check sameFields(dx,saved)
      conv(y,x,p,ws)
      check abs(dot(y,dy)-rhs) <= tol*max(1.0,abs(rhs))
      dy.fill(T(0))
      convVjp(dx,dy,p,ws)
      check dot(dx,dx) == 0

    test "equivalent workspaces share plans and own their halo buffers":
      let x = fields[T](lo,1,T(2))
      let y = fields[T](lo,1)
      let z = fields[T](lo,1)
      let p = convParams(1,1,[3,3],newSeqWith(9,T(1)))
      let a = convWorkspace(x[0],p)
      let b = convWorkspace(x[0],p)
      check a.layout == b.layout
      check a.map == b.map
      check a.halo[0] != b.halo[0]
      conv(y,x,p,a)
      x.fill(T(3))
      conv(z,x,p,b)
      for s in 0..<lo.nSites:
        check y[0].sample(s) == T(18)
        check z[0].sample(s) == T(27)

    test "wide kernel pullbacks contract periodic taps and clear reused output":
      proc value(c,r,s: int): T = T(0.17+0.03*float(c)+0.07*float(r)-0.02*float(s))
      proc seed(c,r,s: int): T = T(0.09-0.05*float(c)+0.02*float(r)+0.03*float(s))
      let x = fields[T](lo,2)
      let dy = fields[T](lo,3)
      let y = fields[T](lo,3)
      let dx = fields[T](lo,2)
      for c in 0..<x.len:
        for s in lo.sites: x[c]{s} := value(c,lo.coords[0][s].int,lo.coords[1][s].int)
      for c in 0..<dy.len:
        for s in lo.sites: dy[c]{s} := seed(c,lo.coords[0][s].int,lo.coords[1][s].int)
      let p = convParams(2,3,[3,5],(0..<90).toSeq.mapIt(T(0.03*cos(float(it)))))
      let ws = convWorkspace(x[0],p)
      var dw = newSeqWith(90,T(71))
      convWeightVjp(dw,x,dy,p,ws)
      let tol = when T is float32: 8e-6 else: 2e-12
      var pair = 0.0
      for o in 0..<3:
        for i in 0..<2:
          for t in 0..<15:
            var want = 0.0
            for r in 0..<8:
              for s in 0..<16:
                let rr = (r+t div 5-1+8) mod 8
                let ss = (s+t mod 5-2+16) mod 16
                want += float(value(i,rr,ss)*seed(o,r,s))
            let k = (o*2+i)*15+t
            check abs(float(dw[k])-want) <= tol*max(1.0,abs(want))
            pair += float(dw[k])*float(p.weights[k])
      conv(y,x,p,ws)
      check abs(pair-dot(y,dy)) <= tol*max(1.0,abs(pair))
      convVjp(dx,dy,p,ws)
      check abs(pair-dot(x,dx)) <= tol*max(1.0,abs(pair))
      let before = dw.mapIt(it)
      for i in 0..<dw.len: dw[i] = T(-81)
      convWeightVjp(dw,x,dy,p,ws)
      check dw == before
      dy.fill(T(0))
      convWeightVjp(dw,x,dy,p,ws)
      for v in dw: check v == T(0)

    test "replicated weight gradients scale once with the global volume":
      for dims in [@[8,16],@[16,16]]:
        let layout = newLayout(dims)
        let x = fields[T](layout,2,T(2))
        let dy = fields[T](layout,3,T(3))
        let p = convParams(2,3,[3,3],newSeq[T](54))
        let ws = convWorkspace(x[0],p)
        var dw = newSeq[T](54)
        convWeightVjp(dw,x,dy,p,ws)
        for v in dw: check v == T(6*layout.physVol)

    test "weight pullbacks of masked outputs take the masked seed":
      let x = fields[T](lo,1,T(2))
      let dy = fields[T](lo,1,T(3))
      let sdy = fields[T](lo,1)
      let y = fields[T](lo,1)
      let mask = realField(lo,float32)
      for s in lo.sites: mask{s} := (if lo.coords[0][s] mod 4 == 0: 1'f32 else: 0'f32)
      let p = convParams(1,1,[3,3],newSeqWith(9,T(0.25)))
      let ws = convWorkspace(x[0],p)
      var dw = newSeq[T](9)
      maskVjp(sdy,dy,"even",mask)
      convWeightVjp(dw,x,sdy,p,ws)
      conv(y,x,p,ws)
      maskVjp(y,y,"even",mask)
      var pair = 0.0
      for v in dw: pair += 0.25*float(v)
      check pair == dot(y,dy)
      let zero = fields[T](lo,1)
      convWeightVjp(dw,zero,sdy,p,ws)
      for v in dw: check v == T(0)

    test "construction and alias mistakes fail before convolution writes":
      let x = fields[T](lo,1)
      let y = fields[T](lo,1)
      let p = convParams(1,1,[1,1],@[T(1)])
      let ws = convWorkspace(x[0],p)
      expect ValueError: discard convParams(1,1,[2,1],@[T(1),T(1)])
      expect ValueError: discard convParams(1,2,[1,1],@[T(1)])
      expect ValueError: conv(x,x,p,ws)
      expect ValueError: convVjp(x,x,p,ws)
      expect ValueError: conv(fields[T](newLayout(@[8,16]),1),x,p,ws)
      expect ValueError: maskedCopy(y,x,realField(newLayout(@[8,16]),float32))
      var dw = newSeq[T](2)
      expect ValueError: convWeightVjp(dw,x,y,p,ws)
