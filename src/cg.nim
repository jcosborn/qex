import qex
import qcdTypes
import profile

type
  SolverParams* = object
    r2req*:float
    maxits*:int
    verbosity*:int
    finalIterations*:int
    subset*:Subset

proc cgSolve*(x:Field; b:Field2; A:proc; sp:var SolverParams) =
  tic()
  let vrb = sp.verbosity
  template verb(n:int; body:expr):untyped =
    if vrb>=n: body
  let subset = sp.subset
  template mythreads(body:untyped):untyped =
    threads:
      onNoSync(subset):
        body

  var b2:float
  mythreads:
    x := 0
    b2 = b.norm2
  verb(1):
    echo("input norm2: ", b2)
  if b2 == 0.0: 
    sp.finalIterations = 0
    return

  var r = newOneOf(x)
  var p = newOneOf(x)
  var Ap = newOneOf(x)
  var r2,r2o,r2stop,fr2:float
  var itn = 0

  mythreads:
    p := 0
    r := b
  r2 = b2
  r2o = r2
  r2stop = sp.r2req * b2;
  verb(1):
    #echo(-1, " ", r2)
    echo(itn, " ", r2/b2)

  toc("cg setup")
  while itn<sp.maxits and r2>=r2stop:
    tic()
    inc itn
    let beta = r2/r2o;
    r2o = r2
    threads:
      #toc("begin threads")
      p[subset] := r + beta*p
      #toc("p update", flops=2*numNumbers(r[0])*subset.lenOuter)
      #echo("p2: ", p.norm2)
      A(Ap, p)
      toc("Ap")
      #echo("Ap2: ", Ap.norm2)
      onNoSync(subset):
        #toc("onNoSync")
        #threadBarrier()
        let pAp = p.redot(Ap)
        #threadBarrier()
        #toc("pAp")
        let alpha = r2/pAp
        x += alpha*p
        #toc("x")
        r -= alpha*Ap
        #toc("r")
        #threadBarrier()
        r2 = r.norm2
        #threadBarrier()
        toc("r2")
      verb(2):
        #echo(itn, " ", r2)
        echo(itn, " ", r2/b2)
      verb(3):
        A(Ap, x)
        onNoSync(subset):
          fr2 = (b - Ap).norm2
        echo "   ", fr2/b2
  toc("cg iterations")

  threads:
    A(Ap, x)
    onNoSync(subset):
      r := b - Ap
      fr2 = r.norm2
  verb(1):
    echo itn, " ", fr2/b2
  sp.finalIterations = itn
  toc("cg final")

when isMainModule:
  qexInit()
  echo "rank ", myRank, "/", nRanks
  #var lat = [8,8,8,8]
  var lat = [4,4,4,4]
  var lo = newLayout(lat)
  var m = lo.ColorMatrix()
  var v1 = lo.ColorVector()
  var v2 = lo.ColorVector()
  proc op*(r:type(v1); x:type(v1)) =
    r := m*x
    #mul(r, m, x)
  var sp:SolverParams
  sp.r2req = 1e-14
  sp.maxits = 100
  sp.verbosity = 3
  sp.subset.layoutSubset(lo, "all")
  threads:
    m.even := 1
    m.odd := 10
    threadBarrier()
    for i in m:
      m{i} := i+1
    threadBarrier()
    v1.even := 1
    v1.odd := 2
    v2 := 0
    echo v1.norm2
    echo m.norm2

  cgSolve(v2, v1, op, sp)
  echo sp.finalIterations
