import base
import layout
import field
import solverBase
export solverBase
import physics/qcdTypes

# BICGSTAB solutions to Mphi = b
#  see https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method
proc bicgstabSolve*(x:Field; b:Field2; A:proc; sp:var SolverParams) =
  tic()
  let vrb = sp.verbosity
  template verb(n:int; body:untyped):untyped =
    if vrb>=n: body
  let sub = sp.subset
  template subset(body:untyped):untyped =
    onNoSync(sub):
      body
  template mythreads(body:untyped):untyped =
    threads:
      onNoSync(sub):
        body
  #type Cmplx = type(x.dot(b))
  type Cmplx = DComplex

  var b2: float
  mythreads:
    x := 0
    b2 = b.norm2
  verb(1):
    echo("input norm2: ", b2)
  if b2 == 0.0:
    sp.finalIterations = 0
    return

  var r = newOneOf(x)
  var r0 = newOneOf(x)
  var v = newOneOf(x)
  var p = newOneOf(x)
  var s = newOneOf(x)
  var t = newOneOf(x)
  let r2stop = sp.r2req * b2;
  let maxits = sp.maxits
  var finalIterations = 0

  threads:
    subset:
      p := 0
      v := 0
      s := 0
      t := 0
      r := b
      r0 := r

    var rho: Cmplx
    rho := 1
    var alpha: Cmplx
    alpha := 1
    var omega: Cmplx
    omega := 1

    var itn = 0
    var r2 = b2
    verb(1):
      #echo(-1, " ", r2)
      echo(itn, " ", r2/b2)
    toc("bicgstab setup")

    while itn<maxits and r2>=r2stop:
      tic()
      inc itn
      # rhoNew = <rhat0, r>
      var rhoNew: Cmplx
      subset:
        rhoNew := r0.dot(r)
      #echo "rhoNew: ", rhoNew

      # beta = (rhoNew/rho)*(alpha/omega)
      let beta = (rhoNew/rho)*(alpha/omega)
      let rho = rhoNew

      # p = r + beta(p - omega v)
      subset:
        p := r + beta*(p - omega*v)
      toc("p update", flops=4*numNumbers(r[0])*sub.lenOuter)
      #echo("p2: ", p.norm2)

      # Apply the matrix.
      threadBarrier()
      A(v, p)
      threadBarrier()
      toc("v")

      # Update s.
      subset:
        # alpha = rho/<rhat0, v>
        let a = r0.dot(v)
        alpha := rho/a
        #alpha := rho/r0.dot(v)
        s := r - alpha*v
        toc("s", flops=2*numNumbers(s[0])*sub.lenOuter)
        r2 := s.norm2
        toc("r2", flops=2*numNumbers(s[0])*sub.lenOuter)

      verb(2):
        echo(itn, " ", r2)

      if r2<=r2stop:
        tic()
        # Update the solution one last time.
        subset:
          x += alpha*p
        toc("x", flops=2*numNumbers(s[0])*sub.lenOuter)
      else:
        tic()
        # t = As
        threadBarrier()
        A(t,s)
        threadBarrier()

        # Update x, r.
        subset:
          # omega = <t,s>/<t,t>
          omega := t.dot(s)/t.norm2
          x += alpha*p + omega*s
          toc("x", flops=4*numNumbers(s[0])*sub.lenOuter)
          r := s - omega*t
          toc("r", flops=3*numNumbers(s[0])*sub.lenOuter)
          r2 := r.norm2
          toc("r2", flops=2*numNumbers(s[0])*sub.lenOuter)

    toc("bicgstab iterations")
    if threadNum==0: finalIterations = itn

    var fr2: float
    A(t, x)
    subset:
      r := b - t
      fr2 = r.norm2
    verb(1):
      echo finalIterations, " acc r2:", r2/b2
      echo finalIterations, " tru r2:", fr2/b2

  sp.finalIterations = finalIterations
  toc("bicgstab final")

when isMainModule:
  import qex
  import physics/qcdTypes
  qexInit()
  echo "rank ", myRank, "/", nRanks
  #var lat = [8,8,8,8]
  var lat = [4,4,4,4]
  var lo = newLayout(lat)
  var m = lo.ColorMatrixD()
  var v1 = lo.ColorVectorD()
  var v2 = lo.ColorVectorD()
  proc op*(r:type(v1); x:type(v1)) =
    r := m*x
    #mul(r, m, x)
  var sp: SolverParams
  sp.r2req = 1e-10
  sp.maxits = 200
  sp.verbosity = 3
  sp.subset.layoutSubset(lo, "all")
  threads:
    m.even := 1
    m.odd := 10
    threadBarrier()
    tfor i, 0, lo.nSites-1:
      m{i} := i+1
    threadBarrier()
    v1.even := 1
    v1.odd := 2
    v2 := 0
    threadBarrier()
    echo v1.norm2
    echo m.norm2

  bicgstabSolve(v2, v1, op, sp)
  echo sp.finalIterations
