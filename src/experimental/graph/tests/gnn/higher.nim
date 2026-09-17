## Seeded input pullbacks and losses that train through input derivatives.
proc higher[T: SomeFloat]() =
  let lo = newLayout(@[8,12])
  let tol = when T is float64: 2e-11 else: 2e-5
  proc fields(n: int, phase: float): seq[RealField[T]] =
    result = newSeq[RealField[T]](n)
    for c in 0..<n:
      result[c] = numeric.realField(lo,T)
      for s in lo.sites:
        let r = float(lo.coords[0][s])
        let k = float(lo.coords[1][s])
        result[c]{s} := T(0.47*sin(0.31*r+0.19*k+0.73*float(c)+phase)-0.13)
  proc close(a,b: float): bool = abs(a-b) <= tol*max(1.0,abs(b))
  proc derivatives(x: float, atan: bool): array[4,float] =
    if atan:
      let q = 1+x*x
      [1/q,-2*x/(q*q),(6*x*x-2)/(q*q*q),24*x*(1-x*x)/(q*q*q*q)]
    else:
      let pdf = exp(-0.5*x*x)/sqrt(2*PI)
      [0.5*erfc(-x/sqrt(2.0))+x*pdf,(2-x*x)*pdf,
       (x*x*x-4*x)*pdf,(-x*x*x*x+7*x*x-4)*pdf]

  suite "NN higher derivatives " & $T:
    test "nonconstant input pullbacks through fourth order match analytic derivatives":
      for flags in [(false,false),(true,false),(false,true)]:
        let atan = flags[0]
        let live = flags[1]
        let rt = initGraphRuntime()
        let x = gnn.toGvalue(rt,fields(2,0.27))
        var seeds: seq[Greal[T]]
        for k in 0..3: seeds.add gnn.toGvalue(rt,fields(2,0.41+0.37*float(k)))
        let y = if atan: gnn.arctan(x) else: gnn.gelu(x)
        var dx = gnn.gradSeeded(y,x,if live: x*x else: seeds[0])
        var roots: seq[Gvalue]
        for k in 0..3:
          roots.add dx
          if k < 3: dx = gnn.gradSeeded(dx,x,seeds[k+1])
        let clones = cloneValues(roots)
        let p = plan(roots)
        proc checkRoots(xs: seq[Gvalue], label: string) =
          for c in 0..<x.fval.len:
            for s in lo.sites:
              var a: float
              a := x.fval[c]{s}
              var ds = derivatives(a,atan)
              var fac = 1.0
              if live:
                let orig = ds
                ds = [a*a*orig[0],2*a*orig[0]+a*a*orig[1],
                      2*orig[0]+4*a*orig[1]+a*a*orig[2],6*orig[1]+6*a*orig[2]+a*a*orig[3]]
              else: fac := seeds[0].fval[c]{s}
              for k in 0..3:
                if k > 0:
                  var cot: float
                  cot := seeds[k].fval[c]{s}
                  fac *= cot
                var got: float
                got := Greal[T](xs[k]).fval[c]{s}
                doAssert close(got,ds[k]*fac), label & " atan=" & $atan & " live=" & $live & " input order " & $(k+1) & ": " & $got & " != " & $(ds[k]*fac)
        for state in 0..1:
          if state == 1:
            x.update(fields(2,-0.63))
            for k in 0..3: seeds[k].update(fields(2,-0.19-0.53*float(k)))
          for z in roots: discard z.eval
          for z in clones: discard z.eval
          checkRoots(roots,"direct " & $state)
          checkRoots(clones,"clone " & $state)
          checkRoots(p.eval,"plan " & $state)
        p.clear
        checkRoots(p.eval,"rebuild")
        p.clear

    test "mixed parameter and input derivatives match one affine GELU composition":
      let rt = initGraphRuntime()
      let x = gnn.toGvalue(rt,fields(1,0.17))
      let seed = gnn.toGvalue(rt,fields(1,0.61))
      var dirs: seq[Greal[T]]
      for v in [0.75,-0.875,0.625,-1.125]:
        let raw = @[numeric.realField(lo,T)]
        let value = T(v)
        threads: raw[0] := value
        dirs.add gnn.toGvalue(rt,raw)
      let w = toGarray(rt,[T(0.75)],[1,1,1,1])
      let b = toGarray(rt,[T(-0.125)],[1])
      let sc = toGarray(rt,[T(1.5)],[1])
      let y = gnn.scale(gnn.gelu(gnn.bias(gnn.conv(x,w),b)),sc)
      let loss = gnn.redot(y,seed)
      var inputs = @[loss]
      for d in dirs: inputs.add gnn.redot(grad(inputs[^1],x),d)
      var roots: seq[Gvalue]
      for par in [w,b,sc]:
        var first = gnn.redot(grad(loss,par),Garray[T](par.oneLike))
        for n in 0..4:
          roots.add grad(inputs[n],par)
          roots.add first
          if n < 4: first = gnn.redot(grad(first,x),dirs[n])
      # d_x^n [s GELU(wx+b)] = s w^n GELU^(n)(wx+b).
      # Differentiate this expression once in w, b, or s.
      var want = newSeq[float](15)
      let weight = float(w.data[0])
      let scale = float(sc.data[0])
      for site in lo.sites:
        var value, fac: float
        value := x.fval[0]{site}
        fac := seed.fval[0]{site}
        let a = weight*value+float(b.data[0])
        let pdf = exp(-0.5*a*a)/sqrt(2*PI)
        let cdf = 0.5*erfc(-a/sqrt(2.0))
        let ds = [a*cdf,cdf+a*pdf,(2-a*a)*pdf,(a*a*a-4*a)*pdf,
                  (-a*a*a*a+7*a*a-4)*pdf,(a*a*a*a*a-11*a*a*a+18*a)*pdf]
        for n in 0..4:
          let power = pow(weight,float(n))
          let slope = if n == 0: 0.0 else: float(n)*pow(weight,float(n-1))*ds[n]
          want[n] += fac*scale*(slope+power*value*ds[n+1])
          want[5+n] += fac*scale*power*ds[n+1]
          want[10+n] += fac*power*ds[n]
          if n < 4:
            var d: float
            d := dirs[n].fval[0]{site}
            fac *= d
      lo.comm.rankSum(want)
      let p = plan(roots)
      let values = p.eval
      for j in 0..<want.len:
        let direct = float(Garray[T](roots[2*j].eval).data[0])
        let planned = float(Garray[T](values[2*j]).data[0])
        let reversed = Gscalar(values[2*j+1]).sval
        check close(direct,want[j])
        check close(planned,want[j])
        check close(reversed,want[j])
      p.clear
