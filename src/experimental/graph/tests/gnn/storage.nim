## Storage and graph runtime contracts using synthetic inputs.
proc storage[S: SomeFloat](rt: GraphRuntime) =
  let tol = when S is float64: 2e-12 else: 2e-6
  proc fields(lo: auto; n: int; value = S(0)): seq[RealField[S]] =
    result = newSeq[RealField[S]](n)
    for i in 0..<n:
      result[i] = numeric.realField(lo,S)
      for s in 0..<lo.nSites: result[i]{s} := value

  proc checkFields(a, b: seq[RealField[S]]; label: string) =
    doAssert a.len == b.len
    let lo = a[0].l
    for c in 0..<a.len:
      for s in 0..<lo.nSites:
        var x, y: float
        x := a[c]{s}
        y := b[c]{s}
        doAssert abs(x-y) <= tol*max(1.0,abs(y)), label & " differs"

  suite "NN graph storage " & $S:
    test "ownership, freshness, shape errors, and lambda clone isolation":
      let lo = newLayout(@[8,8])
      let fs = fields(lo,2)
      for c in 0..1:
        for s in 0..<lo.nSites:
          fs[c]{s} := S(0.21+0.13*float(c)+0.03*float(lo.coords[0][s])-0.017*float(lo.coords[1][s]))
      let x = gnn.toGvalue(rt,fs)
      var old: S
      old := x.fval[0]{0}
      fs[0]{0} := S(99)
      var value: S
      value := x.fval[0]{0}
      doAssert value == old
      fs[0]{0} := old
      let seed = gnn.toGvalue(rt,fields(lo,3,S(0.37)))
      let w = toGarray(rt,(0..<54).toSeq.mapIt(S(0.04*sin(float(it)+0.2))),[3,2,3,3])
      let b = toGarray(rt,[S(0.11),S(-0.07),S(0.23)],[3])
      let sc = toGarray(rt,[S(0.9),S(-0.6),S(1.2)],[3])
      let cv = gnn.conv(x,w,b)
      let y = gnn.scale(gnn.gelu(cv),sc)
      let dx: Greal[S] = gnn.gradSeeded(y,x,seed)
      let dw: Garray[S] = grad(gnn.redot(y,seed),w)
      let db: Garray[S] = grad(gnn.redot(y,seed),b)
      let ds: Garray[S] = grad(gnn.redot(y,seed),sc)
      doAssert y.runCount == 0 and cv.runCount == 0 and dx.runCount == 0
      proc direct() =
        let p = numeric.convParams(2,3,[3,3],w.data,b.data)
        let ws = numeric.convWorkspace(x.fval[0],p)
        let a = fields(lo,3)
        let z = fields(lo,3)
        let act = fields(lo,3)
        let d = fields(lo,3)
        let gx = fields(lo,2)
        numeric.conv(a,x.fval,p,ws)
        numeric.gelu(act,a)
        numeric.scale(z,act,sc.data)
        numeric.scaleVjp(d,seed.fval,sc.data)
        numeric.geluVjp(d,a,d)
        numeric.convVjp(gx,d,p,ws)
        checkFields(y.eval.fval,z,"graph forward")
        checkFields(dx.eval.fval,gx,"graph VJP")
        var dwant = newSeq[S](w.data.len)
        var bwant = newSeq[S](b.data.len)
        var swant = newSeq[S](sc.data.len)
        numeric.convWeightVjp(dwant,x.fval,d,p,ws)
        numeric.channelSum(bwant,d)
        threads:
          for c in 0..<z.len: z[c] := seed.fval[c]*act[c]
        numeric.channelSum(swant,z)
        for pair in [(dw,dwant),(db,bwant),(ds,swant)]:
          let got = pair[0].eval.data
          for i in 0..<got.len:
            doAssert abs(got[i]-pair[1][i]) <= S(tol)*max(S(1),abs(pair[1][i]))
      direct()
      let runs = y.runCount
      discard y.eval
      discard dx.eval
      doAssert y.runCount == runs
      fs[0]{0} := S(0.83)
      x.update(fs)
      direct()
      doAssert y.runCount == runs+1
      var wd = w.data.mapIt(it)
      wd[4] += S(0.19)
      w.update(wd)
      direct()
      var bd = b.data.mapIt(it)
      bd[1] -= S(0.17)
      b.update(bd)
      var sd = sc.data.mapIt(it)
      sd[2] = S(-0.4)
      sc.update(sd)
      direct()
      doAssertRaises(GraphValueError): discard gnn.gradSeeded(y,x,gnn.toGvalue(rt,fields(lo,2)))
      doAssertRaises(GraphValueError): discard gnn.conv(x,toGarray(rt,newSeq[S](81),[3,3,3,3]))
      doAssertRaises(ValueError): discard gnn.conv(x,toGarray(rt,newSeq[S](24),[3,2,2,2]))
      doAssertRaises(GraphValueError): discard toGarray(rt,newSeq[S](7),[2,4])
      doAssertRaises(GraphValueError): w.update(newSeq[S](3))
      doAssertRaises(GraphValueError): x.update(fields(newLayout(@[8,8]),2))
      doAssertRaises(GraphValueError): discard x+gnn.toGvalue(initGraphRuntime(),fs)
      let shared = x+x
      let loss = gnn.redot(shared,shared)
      let total: Greal[S] = grad(loss,x)
      let expected = fields(lo,2)
      threads:
        for c in 0..1: expected[c] := S(8)*x.fval[c]
      checkFields(total.eval.fval,expected,"shared input accumulation")
      let copy = Greal[S](x.newOneOf)
      copy.valCopy(x)
      copy.updated
      doAssert copy.fval[0] != x.fval[0]
      copy.fval[0]{0} := S(27)
      value := x.fval[0]{0}
      doAssert value != S(27)
      let acopy = Garray[S](w.newOneOf)
      acopy.valCopy(w)
      acopy.data[0] = S(31)
      doAssert w.data[0] != S(31)
      doAssert not x.copyCompatible(seed)
      doAssert not w.copyCompatible(sc)
      doAssert Greal[S](x.zeroLike).isStaticZeroLeaf
      doAssert Garray[S](w.zeroLike).isStaticZeroLeaf
      let zeroConv = gnn.conv(Greal[S](x.zeroLike),Garray[S](w.zeroLike),Garray[S](b.zeroLike))
      checkFields(zeroConv.eval.fval,fields(lo,3),"zero-valued graph inputs")
      doAssert zeroConv.runCount == 1
      let one = Greal[S](x.oneLike)
      for f in one.eval.fval:
        for s in lo.sites:
          value := f{s}
          doAssert value == 1
      let proto = gnn.toGvalue(rt,fields(lo,2))
      let fn = lambda(proto,gnn.scale(gnn.gelu(gnn.conv(proto,w,b)),sc))
      let otherFields = fields(lo,2,S(-0.29))
      let other = gnn.toGvalue(rt,otherFields)
      let a = Greal[S](apply(fn,x))
      let z = Greal[S](apply(fn,other))
      let da: Greal[S] = gnn.gradSeeded(a,x,seed)
      let dz: Greal[S] = gnn.gradSeeded(z,other,seed)
      checkFields(a.eval.fval,y.eval.fval,"cloned forward")
      checkFields(da.eval.fval,dx.eval.fval,"cloned VJP")
      doAssert a.fval[0] != z.fval[0]
      let directOther = gnn.scale(gnn.gelu(gnn.conv(other,w,b)),sc)
      checkFields(z.eval.fval,directOther.eval.fval,"second clone forward")
      checkFields(dz.eval.fval,gnn.gradSeeded(directOther,other,seed).eval.fval,"second clone VJP")
      otherFields[0]{0} := S(1.1)
      other.update(otherFields)
      checkFields(z.eval.fval,directOther.eval.fval,"clone input refresh")
      checkFields(a.eval.fval,y.eval.fval,"source isolation after clone refresh")
      wd[5] -= S(0.23)
      w.update(wd)
      checkFields(a.eval.fval,y.eval.fval,"clone parameter refresh")
      checkFields(z.eval.fval,directOther.eval.fval,"second clone parameter refresh")
      checkFields(da.eval.fval,dx.eval.fval,"cloned force parameter refresh")
