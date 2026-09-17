## Field lifetimes and NN evaluation with reusable graph storage.
proc lifetime[T: SomeFloat]() =
  let lo = newLayout(@[8,12])
  let tol = when T is float64: 3e-12 else: 8e-6
  proc fields(n: int, value: T): seq[RealField[T]] =
    result = newSeq[RealField[T]](n)
    for c in 0..<n:
      result[c] = numeric.realField(lo,T)
      for s in lo.sites: result[c]{s} := value
  proc same(a, b: seq[RealField[T]]) =
    check a.len == b.len
    for c in 0..<a.len:
      for s in lo.sites:
        var x,y: float
        x := a[c]{s}
        y := b[c]{s}
        check abs(x-y) <= tol*max(1.0,abs(y))
  proc constant(x: Greal[T], value: T) = same(x.fval,fields(x.fval.len,value))

  suite "NN graph lifetime " & $T:
    test "symbolic values allocate on demand and released aliases retain samples":
      let rt = initGraphRuntime()
      let fs = fields(2,T(0.3))
      let x = gnn.toGvalue(rt,fs)
      let w = toGarray(rt,newSeq[T](36),[2,2,3,3])
      let mem = getRawMemAllocated()
      let z = Greal[T](x.zeroLike)
      let one = Greal[T](x.oneLike)
      let cv = gnn.conv(x,w)
      let y = gnn.gelu(cv)
      let dx: Greal[T] = grad(gnn.sum(y),x)
      check not z.hasStorage and not one.hasStorage
      check not cv.hasStorage and not y.hasStorage and not dx.hasStorage
      check getRawMemAllocated() == mem
      check x.bufferCompatible(z) and z.bufferCompatible(x)
      check z.bufferBytes == x.bufferBytes
      let a = x.fval
      let bytes = x.bufferBytes
      x.releaseStorage
      check not x.hasStorage and x.bufferBytes == bytes
      same(a,fs)
      x.update(fs)
      check x.hasStorage and x.fval[0] != a[0]
      same(x.fval,fs)
      constant(one.eval,T(1))
      constant(z.eval,T(0))
      one.releaseStorage
      z.releaseStorage
      constant(one.eval,T(1))
      constant(z.eval,T(0))
      constant(y.eval,T(0))
      let previous = y.fval
      y.releaseStorage
      check not y.hasStorage
      constant(y.eval,T(0))
      same(y.fval,previous)
      check y.fval[0] != previous[0]
      check y.runCount == 2

    test "plans reuse field and convolution storage through input and parameter updates":
      let rt = initGraphRuntime()
      let fs = fields(2,T(0.19))
      let x = gnn.toGvalue(rt,fs)
      let seed = gnn.toGvalue(rt,fields(2,T(-0.27)))
      let w = toGarray(rt,(0..<36).toSeq.mapIt(T(0.07*cos(float(it)+0.4))),[2,2,3,3])
      let b = toGarray(rt,[T(0.03),T(-0.11)],[2])
      let a = gnn.gelu(gnn.conv(x,w,b))
      let y = gnn.gelu(gnn.conv(a,w,b))
      let dx: Greal[T] = gnn.gradSeeded(y,x,seed)
      let h: Greal[T] = grad(gnn.redot(dx,seed),x)
      let dw: Garray[T] = grad(gnn.redot(y,seed),w)
      let db: Garray[T] = grad(gnn.redot(y,seed),b)
      let p = plan(y,dx,h,dw,db)
      check not y.hasStorage and not dx.hasStorage and not h.hasStorage
      discard p.eval
      check not y.hasStorage and not dx.hasStorage and not h.hasStorage
      check p.stats.buffers > 0 and p.stats.reuses > 0
      check p.stats.workspaces == 1
      proc compare() =
        same(Greal[T](p[0]).fval,y.eval.fval)
        same(Greal[T](p[1]).fval,dx.eval.fval)
        same(Greal[T](p[2]).fval,h.eval.fval)
        for i, want in [dw.eval.data,db.eval.data]:
          let got = Garray[T](p[3+i]).data
          for j in 0..<want.len:
            check abs(got[j]-want[j]) <= T(tol)*max(T(1),abs(want[j]))
      compare()
      let count = p.stats.forwards
      let mem = getRawMemAllocated()
      discard p.eval
      check p.stats.forwards == count
      check getRawMemAllocated() == mem
      let buffers = p.stats.buffers
      for c in 0..<fs.len:
        for s in lo.sites: fs[c]{s} := T(0.13+0.07*float(c)-0.02*float(lo.coords[0][s]))
      x.update(fs)
      var weights = w.data.mapIt(it)
      weights[4] += T(0.17)
      w.update(weights)
      discard p.eval
      check p.stats.buffers == buffers and p.stats.workspaces == 1
      check getRawMemAllocated() == mem
      compare()
      let old = Greal[T](p[0])
      let alias = old.fval
      let saved = alias.newOneOf
      saved.copyFieldStorage(alias)
      p.clear
      check p.stats.buffers == 0 and p.stats.workspaces == 0
      x.update(fields(2,T(-0.31)))
      discard p.eval
      check old.stale
      expect GraphError: discard old.eval
      same(alias,saved)
      compare()

    test "planned channel pullbacks clear untouched channels on reuse":
      let rt = initGraphRuntime()
      let x = gnn.toGvalue(rt,fields(2,T(0.3)))
      let dx: Greal[T] = grad(gnn.sum(gnn.channel(x,0)),x)
      let p = plan(dx)
      let expected = fields(2,T(0))
      threads: expected[0] := 1
      discard p.eval
      same(Greal[T](p[0]).fval,expected)
      Greal[T](p[0]).update(fields(2,T(17)))
      discard p.eval
      same(Greal[T](p[0]).fval,expected)

    test "planned conditionals allocate work only for the selected branch":
      let rt = initGraphRuntime()
      let x = gnn.toGvalue(rt,fields(1,T(0.25)))
      let w = toGarray(rt,[T(2)],[1,1,1,1])
      let sel = scalar.toGvalue(rt,0)
      let y = cond(sel,gnn.conv(x,w),x*x)
      let dx: Greal[T] = grad(gnn.sum(y),x)
      let p = plan(y,dx)
      discard p.eval
      constant(Greal[T](p[0]),T(0.0625))
      constant(Greal[T](p[1]),T(0.5))
      check p.stats.workspaces == 0
      sel.update(1)
      discard p.eval
      constant(Greal[T](p[0]),T(0.5))
      constant(Greal[T](p[1]),T(2))
      check p.stats.workspaces == 1
      x.update(fields(1,T(-0.25)))
      sel.update(0)
      discard p.eval
      constant(Greal[T](p[0]),T(0.0625))
      constant(Greal[T](p[1]),T(-0.5))
