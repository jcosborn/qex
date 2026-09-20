## Analytic graph contracts from elementary values and adjoint identities.
proc contracts[T: SomeFloat](rt: GraphRuntime) =
  let lo = newLayout(@[8,16])
  let tol = when T is float64: 3e-12 else: 8e-6
  proc fields(n: int, value: T): seq[RealField[T]] =
    result = newSeq[RealField[T]](n)
    for c in 0..<n:
      result[c] = numeric.realField(lo,T)
      for s in lo.sites: result[c]{s} := value
  proc leaf(value: T, n = 1): Greal[T] = gnn.toGvalue(rt,fields(n,value))
  proc sample(x: Greal[T], c=0, s=0): float =
    result := x.eval.fval[c]{s}
  proc close(a,b: float): bool = abs(a-b) <= tol*max(1.0,abs(b))
  proc checkConstant(x: Greal[T], expected: float) =
    let fs = x.eval.fval
    for f in fs:
      for s in lo.sites:
        var v: float
        v := f{s}
        doAssert close(v,expected), $v & " != " & $expected

  suite "NN graph contracts " & $T:
    test "shared field storage retains shape, type, channels, and owned values":
      let fs = fields(2,T(0.3))
      fs[1]{0} := T(-0.8)
      let x = gnn.toGvalue(rt,fs)
      static: doAssert x is GfieldOf[seq[RealField[T]]]
      fs[0]{0} := T(99)
      check close(sample(x),0.3)
      check close(sample(x,1),-0.8)
      let z = Greal[T](x.newOneOf)
      z.valCopy(x)
      check z.fval[0] != x.fval[0]
      z.fval[0]{0} := T(-7)
      check close(sample(x),0.3)
      expect GraphValueError:
        discard x+leaf(T(0),3)
      let one = Greal[T](x.oneLike)
      checkConstant(one,1)
      check x.copyCompatible(one)
      let zero = Greal[T](x.zeroLike)
      checkConstant(gnn.exp(zero),1)
      checkConstant(gnn.erfc(zero),1)
      checkConstant(gnn.clipMin(zero,T(0.5)),0.5)
      checkConstant(zero+T(0.7),0.7)
      when T is float32:
        # The single-field storage family exists for Gmask only.
        let site = toGfield(rt,fs[0])
        static: doAssert site is GfieldOf[RealField[T]]
        let scopy = GfieldOf[RealField[T]](site.newOneOf)
        scopy.valCopy(site)
        check site.copyCompatible(scopy)
        check scopy.fval != site.fval

    test "activation values and derivatives are analytic through second order":
      for v in [T(-0.7),T(0),T(0.6)]:
        let x = leaf(v)
        let seed = leaf(T(0.23))
        let g = gnn.gelu(x)
        let a = gnn.arctan(x)
        let dg: Greal[T] = grad(gnn.sum(g),x)
        let da: Greal[T] = grad(gnn.sum(a),x)
        let hg: Greal[T] = grad(gnn.redot(dg,seed),x)
        let ha: Greal[T] = grad(gnn.redot(da,seed),x)
        check g.runCount == 0 and dg.runCount == 0 and hg.runCount == 0
        let d = float(v)
        let pdf = exp(-0.5*d*d)/sqrt(2*PI)
        let cdf = 0.5*erfc(-d/sqrt(2.0))
        checkConstant(g,d*cdf)
        checkConstant(dg,cdf+d*pdf)
        checkConstant(hg,0.23*(2-d*d)*pdf)
        checkConstant(da,1/(1+d*d))
        checkConstant(ha,-0.46*d/((1+d*d)*(1+d*d)))
        checkConstant(g-gnn.gelu(-x),d)
        let same = gnn.gradSeeded(g,x,x)
        let live = gnn.gradSeeded(g,x,x*x)
        let dsame: Greal[T] = grad(gnn.sum(same),x)
        let dlive: Greal[T] = grad(gnn.sum(live),x)
        checkConstant(dsame,cdf+d*pdf+d*(2-d*d)*pdf)
        checkConstant(dlive,2*d*(cdf+d*pdf)+d*d*(2-d*d)*pdf)

    test "convolution and its transpose close the fixed weight input Hessian":
      let x = leaf(T(0.31))
      let seed = leaf(T(-0.17))
      let w = toGarray(rt,[T(2)],[1,1,1,1])
      let y = gnn.conv(x,w)
      let dx: Greal[T] = grad(gnn.redot(y,y),x)
      let h: Greal[T] = grad(gnn.redot(dx,seed),x)
      check y.runCount == 0 and dx.runCount == 0 and h.runCount == 0
      checkConstant(dx,8*0.31)
      checkConstant(h,8*(-0.17))
      w.update([T(3)])
      checkConstant(h,18*(-0.17))
      let trans = gnn.convTranspose(seed,w)
      check close(gnn.redot(gnn.conv(x,w),seed).eval.sval,gnn.redot(x,trans).eval.sval)

    test "field masks own samples, refresh, and accumulate shared input slots":
      let x = leaf(T(0.31))
      let fill = leaf(T(-0.42))
      let raw = numeric.realField(lo,float32)
      for s in lo.sites:
        raw{s} := (if (lo.coords[0][s].int+lo.coords[1][s].int) mod 3 == 0: 1'f32 else: 0'f32)
      let mask = toGfield(rt,raw)
      let y = gnn.maskedCopy(x,fill,mask,"even")
      let same = gnn.maskedCopy(x,x,mask,"even")
      let shared: Greal[T] = gnn.gradSeeded(same,x,fill)
      checkConstant(shared,-0.42)
      let dx: Greal[T] = grad(gnn.sum(y),x)
      let df: Greal[T] = grad(gnn.sum(y),fill)
      let sub = lo.getSubset("even")
      discard y.eval
      discard dx.eval
      discard df.eval
      for s in lo.sites:
        let selected = s >= sub.low and s < sub.high and (lo.coords[0][s].int+lo.coords[1][s].int) mod 3 == 0
        check close(sample(y,0,s),if selected: 0.31 else: -0.42)
        check close(sample(dx,0,s),if selected: 1.0 else: 0.0)
        check close(sample(df,0,s),if selected: 0.0 else: 1.0)
      for s in lo.sites: raw{s} := 0'f32
      let runs = y.runCount
      discard y.eval
      check y.runCount == runs
      mask.update(raw)
      checkConstant(y,-0.42)
      check y.runCount == runs+1
      checkConstant(dx,0)
      checkConstant(df,1)

    test "conditionals evaluate only the selected branch and refresh gradients":
      let x = leaf(T(0.23))
      let w = toGarray(rt,[T(2)],[1,1,1,1])
      let a = gnn.conv(x,w)
      let b = x*x
      let sel = scalar.toGvalue(rt,0)
      let y = cond(sel,a,b)
      let dx: Greal[T] = grad(gnn.sum(y),x)
      checkConstant(y,0.23*0.23)
      checkConstant(dx,2*0.23)
      check a.runCount == 0
      sel.update(1)
      checkConstant(y,2*0.23)
      checkConstant(dx,2)
      check a.runCount == 1
      x.update(fields(1,T(-0.4)))
      sel.update(0)
      checkConstant(y,0.16)
      checkConstant(dx,-0.8)
      check a.runCount == 1

    test "channel composition preserves order and shared channel adjoints":
      let x = leaf(T(0.21),2)
      let fs = fields(2,T(0.21))
      for s in lo.sites: fs[1]{s} := T(-0.7)
      x.update(fs)
      let a = gnn.channel(x,1)
      let b = gnn.channel(x,0)
      let y = gnn.concat([a,b,a])
      check y.fval.len == 3
      checkConstant(gnn.channel(y,0),-0.7)
      checkConstant(gnn.channel(y,1),0.21)
      let dx: Greal[T] = grad(gnn.sum(y),x)
      checkConstant(gnn.channel(dx,0),1)
      checkConstant(gnn.channel(dx,1),2)

    test "clipped logarithm uses half derivative at equality":
      for v in [T(0.2),T(0.5),T(0.8)]:
        let x = leaf(v)
        let y = gnn.ln(gnn.clipMin(x,T(0.5)))
        let dx: Greal[T] = grad(gnn.sum(y),x)
        checkConstant(y,ln(max(float(v),0.5)))
        let slope = if v < T(0.5): 0.0 elif v == T(0.5): 1.0 else: 1.0/float(v)
        checkConstant(dx,slope)
    test "typed divisions round each quotient in the field precision":
      # Fast-math may divide through a reciprocal on either side: a few ulps of T.
      const eps = when T is float64: 2.220446049250313e-16 else: 1.1920928955078125e-7
      for v in [T(0.21),T(-1.7),T(2.37)]:
        let x = leaf(v)
        let q = gnn.divide(gnn.divide(x,T(7.3)),T(5.9))
        let b = v/T(7.3)
        let c = b/T(5.9)
        check abs(sample(q)-float(c)) <= 4*eps*abs(float(c))

    test "replicated parameter contractions count each global site once":
      for dims in [@[8,16],@[16,16]]:
        let layout = newLayout(dims)
        let raw = @[numeric.realField(layout,T),numeric.realField(layout,T)]
        threads:
          raw[0] := T(2)
          raw[1] := T(-3)
        let x = gnn.toGvalue(rt,raw)
        let a = toGarray(rt,[T(1.5),T(-0.5)],[2])
        let sums = gnn.channelSum(x)
        let lhs = gnn.redot(sums,a)
        let rhs = gnn.redot(x,gnn.broadcast(a,x))
        let vol = layout.physVol.float
        check lhs.eval.sval == 4.5*vol
        check rhs.eval.sval == lhs.eval.sval
        let dx: Greal[T] = grad(lhs,x)
        let da: Garray[T] = grad(lhs,a)
        check da.eval.data == @[T(2*vol),T(-3*vol)]
        let dfields = dx.eval.fval
        for c in 0..<2:
          for s in layout.sites:
            var v: T
            v := dfields[c]{s}
            check v == a.data[c]
        let both = gnn.redot(gnn.channelSum(gnn.broadcast(a,x)),a)
        let shared: Garray[T] = grad(both,a)
        check shared.eval.data == @[T(3*vol),T(-vol)]
        let repeated = gnn.broadcast(sums,x)
        let sharedField: Greal[T] = grad(gnn.redot(repeated,x),x)
        let fs = sharedField.eval.fval
        for c in 0..<2:
          for s in layout.sites:
            var v: T
            v := fs[c]{s}
            check v == (if c == 0: T(4*vol) else: T(-6*vol))

    test "convolution weight pullbacks remain bilinear graph values":
      let x = leaf(T(0.25))
      let seed = leaf(T(-0.5))
      let w = toGarray(rt,[T(2)],[1,1,1,1])
      let wd = toGarray(rt,[T(3)],[1,1,1,1])
      let loss = gnn.redot(gnn.conv(x,w),seed)
      let reverse = gnn.redot(gnn.convTranspose(seed,w),x)
      let dw: Garray[T] = grad(loss,w)
      let dt: Garray[T] = grad(reverse,w)
      check dw.eval.data == @[T(-0.125*lo.physVol.float)]
      check dt.eval.data == dw.eval.data
      let mixed: Greal[T] = grad(gnn.redot(dw,wd),x)
      let ds: Greal[T] = grad(gnn.redot(dw,wd),seed)
      checkConstant(mixed,-1.5)
      checkConstant(ds,0.75)
      let h: Greal[T] = grad(gnn.sum(mixed),x)
      checkConstant(h,0)
      let shared = gnn.convWeightVjp(x,x,[1,1])
      let twice: Greal[T] = grad(gnn.redot(shared,wd),x)
      checkConstant(twice,1.5)
      x.update(fields(1,T(-0.5)))
      seed.update(fields(1,T(0.25)))
      checkConstant(mixed,0.75)
      checkConstant(ds,-1.5)
      checkConstant(twice,-3)
