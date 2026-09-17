## Both precisions share one runtime and retain distinct typed storage.
proc precision(rt: GraphRuntime) =
  suite "NN graph precision":
    test "values, gradients, copies, and caches coexist in one runtime":
      let lo = newLayout(@[8,8])
      proc field[T: SomeFloat](value: T): seq[RealField[T]] =
        result = @[numeric.realField(lo,T)]
        for i in 0..<lo.nSites: result[0]{i} := value
      let x32 = gnn.toGvalue(rt,field(0.37'f32))
      let x64 = gnn.toGvalue(rt,field(0.3700000001))
      let w32 = toGarray(rt,[0.71'f32],[1,1,1,1])
      let w64 = toGarray(rt,[0.7100000000001],[1,1,1,1])
      let b32 = toGarray(rt,[0.11'f32],[1])
      let b64 = toGarray(rt,[0.110000000000003],[1])
      let y32 = gnn.gelu(gnn.conv(x32,w32,b32))
      let y64 = gnn.gelu(gnn.conv(x64,w64,b64))
      let loss = gnn.redot(y32,y32)+gnn.redot(y64,y64)
      let dx32: Greal[float32] = core.grad(loss,x32)
      let dx64: Greal[float64] = core.grad(loss,x64)
      static:
        doAssert x32 is Greal[float32]
        doAssert x64 is Greal[float64]
        doAssert dx32 is Greal[float32]
        doAssert dx64 is Greal[float64]
      proc verify[T: SomeFloat](x: Greal[T]; w,b: Garray[T]; y,dx: Greal[T]) =
        let p = numeric.convParams(1,1,[1,1],w.data,b.data)
        let ws = numeric.convWorkspace(x.fval[0],p)
        let z = field(T(0))
        let outp = field(T(0))
        let seed = field(T(0))
        let gradp = field(T(0))
        numeric.conv(z,x.fval,p,ws)
        numeric.gelu(outp,z)
        numeric.scale(seed,outp,@[T(2)])
        numeric.geluVjp(seed,z,seed)
        numeric.convVjp(gradp,seed,p,ws)
        discard y.eval
        discard dx.eval
        let tol = when T is float64: 2e-12 else: 2e-6
        for i in 0..<lo.nSites:
          var a,b,c,d: float64
          a := y.fval[0]{i}
          b := outp[0]{i}
          c := dx.fval[0]{i}
          d := gradp[0]{i}
          doAssert abs(a-b) < tol and abs(c-d) < tol
      verify[float32](x32,w32,b32,y32,dx32)
      verify[float64](x64,w64,b64,y64,dx64)
      let runs64 = y64.runCount
      w32.update([0.83'f32])
      discard y64.eval
      doAssert y64.runCount == runs64
      verify[float32](x32,w32,b32,y32,dx32)
      verify[float64](x64,w64,b64,y64,dx64)
      let runs32 = y32.runCount
      x64.update(field(0.410000000003))
      w64.update([0.7300000000001])
      discard y32.eval
      doAssert y32.runCount == runs32
      verify[float64](x64,w64,b64,y64,dx64)
      verify[float32](x32,w32,b32,y32,dx32)
      let values = @[Gvalue(x32),Gvalue(x64),Gvalue(y32),Gvalue(y64),Gvalue(dx32),Gvalue(dx64)]
      for value in values:
        let copied = value.newOneOf
        copied.valCopy(value.eval)
        let zero = value.zeroLike
        doAssert value.copyCompatible(copied) and value.copyCompatible(zero)
        if value of Greal[float32]:
          doAssert copied of Greal[float32] and zero of Greal[float32]
          doAssert not copied.copyCompatible(x64)
        else:
          doAssert copied of Greal[float64] and zero of Greal[float64]
          doAssert not copied.copyCompatible(x32)
      doAssert not w32.copyCompatible(w64) and not w64.copyCompatible(w32)
      for value in @[Gvalue(w32),Gvalue(w64)]:
        let copied = value.newOneOf
        copied.valCopy(value)
        let zero = value.zeroLike
        doAssert value.copyCompatible(copied) and value.copyCompatible(zero)
        if value of Garray[float32]:
          doAssert copied of Garray[float32] and zero of Garray[float32]
          doAssert Garray[float32](zero).data == @[0'f32]
        else:
          doAssert copied of Garray[float64] and zero of Garray[float64]
          doAssert Garray[float64](zero).data == @[0.0]
      doAssertRaises(GraphValueError): x32.valCopy(x64)
      doAssertRaises(GraphValueError): x64.valCopy(x32)
      doAssertRaises(GraphValueError): w32.valCopy(w64)
      doAssertRaises(GraphValueError): w64.valCopy(w32)
      doAssertRaises(GraphValueError): discard x32.addLike(x32,x64)
      doAssertRaises(GraphValueError): discard core.gradSeeded(Gvalue(y32),Gvalue(x32),x64.oneLike)
      let independent: Greal[float64] = core.grad(gnn.redot(y32,y32),x64)
      doAssert independent.isStaticZeroLeaf
      doAssert independent of Greal[float64]
