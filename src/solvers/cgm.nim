# Multi-mass (multi-shift) conjugate gradient
# Based on "Krylov space solvers for shifted linear systems" (arXiv:9612014) 
# and "Reduced-Shifted Conjugate-Gradient Method for a Green's Function: 
# Efficient Numerical Approach in a Nano-structured Superconductor"
# (arXiv:1607.03992v2)

import qex
import physics/qcdTypes
import base
import layout
import field
import solverBase

import math

export solverBase

type
  # Multi-shift conjugate gradient type
  CgmState*[T] = object
    r,b,z,q,Ap,LAp: T
    xs,ps,bs: seq[T]
    sigmas*: seq[float]
    b2,r2*,r2stop*: float
    iterations: int
    precon: CgPrecon

# Reset multi-shift CG state
proc reset*(cgs: var CgmState) =
  cgs.iterations = 0
  cgs.b2 = -1
  cgs.r2 = 1.0
  cgs.r2stop = 0.0

# Initialize preconditioning
proc initPrecon*(state: var CgmState) =
  case state.precon:
    of cpNone:
      state.z = state.r
      state.q = state.ps[0]
      state.LAp = state.Ap
    of cpHerm:
      state.z = newOneof(state.r)
      state.q = state.ps[0]
      state.LAp = state.Ap
    of cpLeftRight:
      state.z = newOneof(state.r)
      state.q = newOneof(state.ps[0])
      state.LAp = newOneOf(state.z)

# Construct new multi-shift CG state
proc createCgmState*[T](
    xs: seq[T];
    b: T; 
    sigmas: seq[float];
    precon: CgPrecon = cpNone
  ): CgmState[T] =
  let nmass = sigmas.len
  result = CgmState[T]()
  result.sigmas = sigmas
  result.xs = xs
  result.ps = newSeq[T]()
  result.bs = newSeq[T]()
  for m in 0..<nmass: result.ps.add newOneOf(xs[m])
  result.precon = precon
  result.reset

proc newCgmState*[T](
    xs: seq[T];
    bs: seq[T]; 
    sigmas: seq[float];
    precon: CgPrecon = cpNone
  ): CgmState[T] = 
  result = createCgmState(xs,bs[0],sigmas,precon=precon)
  for m in 0..<bs.len: result.bs.add bs[m]
  result.b = newOneOf(bs[0])
  result.r = newOneOf(bs[0])
  result.Ap = newOneOf(bs[0])
  result.initPrecon

proc newCgmState*[T](
    xs: seq[T];
    b: T; 
    sigmas: seq[float];
    precon: CgPrecon = cpNone
  ): CgmState[T] = 
  result = createCgmState(xs,b,sigmas,precon=precon)
  result.b = b
  result.r = newOneOf(b)
  result.Ap = newOneOf(b)
  result.initPrecon

proc newCgmState*[T](
    x: T;
    b: T; 
    sigmas: seq[float];
    precon: CgPrecon = cpNone
  ): CgmState[T] =
  var xs = newSeq[T](sigmas.len)
  for m in 0..<sigmas.len:
    xs[m] = case m == 0:
      of true: x
      of false: newOneOf(x)
  result = createCgmState(xs,b,sigmas,precon=precon)
  result.b = b
  result.r = newOneOf(b)
  result.Ap = newOneOf(b)
  result.initPrecon

# Multi-shift CG solver (arXiv:9612014,1607.03992v2)
proc solve*(state: var CgmState; op: auto; sp: var SolverParams) = 
  mixin apply,applyPrecon

  var 
    r2,r2stop: float
    b2 = state.b2
    itn0 = 0
    maxits = sp.maxits
    sub = sp.subset
    restart = false
    sg = newSeq[float]() #state.sigmas
  let 
    nmass = state.xs.len
    vrb = sp.verbosity
    r = state.r
    q = state.q
    z = state.z
    ps = state.ps
    Ap = state.Ap
    LAp = state.LAp
    xs = state.xs
    b = state.b
    bs = state.bs
    precon = state.precon
  for m in 0..<xs.len: 
    case m == 0:
      of true: sg.add 0.0
      of false: sg.add state.sigmas[m]

  if precon != state.precon:
    state.precon = precon
    state.initPrecon
    state.reset
  case precon:
    of cpNone: discard
    of cpHerm:
      when not compiles(op.applyPrecon(z,r)):
        qexError "cg.solve: precon == cpHerm but op.applyPrecon not found"
    of cpLeftRight:
      when not compiles(op.applyPreconL(z,r)):
        qexError "cg.solve: precon == cpLeftRight but op.applyPreconL not found"
      when not compiles(op.applyPreconR(ps[0],q)):
        qexError "cg.solve: precon == cpLeftRight but op.applyPreconR not found"

  template forMass(body:untyped) =
    for k in 0..<nmass: 
      let m {.inject.} = k
      body

  template verb(n:int;body:untyped) =
    if vrb>=n: body

  template mythreads(body:untyped) =
    threads:
      onNoSync(sub): body

  template subset(body:untyped) =
    onNoSync(sub): body

  template preconL(z,r) =
    case precon:
      of cpNone: discard
      of cpHerm:
        when compiles(op.applyPrecon(z,r)): 
          op.applyPrecon(z,r)
      of cpLeftRight: 
        when compiles(op.applyPreconL(z,r)):
          op.applyPreconL(z,r)

  template preconR(p,q) =
    case precon:
      of cpNone,cpHerm: discard
      of cpLeftRight: 
        when compiles(op.applyPreconR(p,q)):
          op.applyPreconR(p,q)

  template prod[T](x:T;m,k:int) =
    for j in 0..<xs.len:
      threadBarrier()
      if (j != m) and (j != k):
        Ap := 0
        threadBarrier()
        op.apply(Ap,x,0.25*sg[j])
        threadBarrier()
        Ap *= 1.0/(sg[j]-sg[k])
        threadBarrier()
        x := Ap

  template x: untyped = xs[0]
  template p: untyped = ps[0]

  mythreads:
    var r2t,b2t: float
    case bs.len == 0:
      of true:
        case b2 < 0.0: # new solution
          of true:
            r := b
            b2t = b.norm2
            forMass: xs[m] := 0.0 
            verb(1): echo "input norm2: ", b2t
          of false: # CG restart
            op.apply(Ap,xs[0])
            r := b - Ap
            b2t = b2
      of false: # arXiv:0810.1081
        for m in 0..<xs.len:
          xs[m] := 0.0
          for k in 0..<xs.len:
            threadBarrier()
            if k != m:
              b := bs[m]-bs[k]
              threadBarrier()
              b *= 1.0/(sg[m]-sg[k])
              threadBarrier() 
              b.prod(m,k)
              threadBarrier()
              xs[m] += b
        threadBarrier()
        op.apply(Ap,xs[0])
        threadBarrier()
        r := bs[0] - Ap
        b2t = bs[0].norm2
    threadBarrier()
    r2t = r.norm2
    verb(3): 
      echo "p2: ", ps[0].norm2
      echo "r2: ", r2t
    threadMaster:
      r2 = r2t
      b2 = b2t
  r2stop = sp.r2req * b2
  state.r2stop = r2stop

  if r2 > r2stop:
    threads:
      var
        itn = 0
        r2r,r2i,r2ip1,qLAp: float
        alpha,beta: float
        alphaim1,betaim1: float
        zi = newSeq[float](nmass)
        zim1 = newSeq[float](nmass)
        continuing = true

      # iteration 0
      (alphaim1,betaim1) = (-1.0,0.0)
      forMass: (zim1[m],zi[m]) = (1.0,1.0)
      preconL(z,r)
      subset:
        q := z
        r2i = case precon
          of cpNone: r2
          of cpHerm: redot(r,z)
          of cpLeftRight: z.norm2 
      preconR(p,q)
      forMass: 
        subset:
          case m == 0:
            of true: discard
            of false: ps[m] := p

      verb(1): echo "CG iteration: ",itn,"  r2/b2: ",r2/b2
      
      # iteration 1,2,3,...
      while continuing:
        forMass:
          threadBarrier()
          case m == 0:
            of true: # sigma = 0 solve w/ optional precond.
              subset: op.apply(Ap,p)
              case precon:
                of cpLeftRight: preconL(LAp,Ap)
                of cpNone,cpHerm: discard
              inc itn
              subset: qLAp = redot(q,LAp)
              alpha = case qLAp != 0.0
                of true: r2i/qLAp
                of false: 0.0
              subset: 
                r -= alpha*Ap
                x += alpha*p
              case precon:
                of cpNone,cpHerm: preconL(z,r)
                of cpLeftRight: 
                  subset: z -= alpha*LAp
              subset: 
                r2ip1 = case precon
                  of cpNone: r.norm2
                  of cpHerm: redot(r,z)
                  of cpLeftRight: z.norm2
                r2r = case precon
                  of cpNone: r2ip1
                  of cpHerm,cpLeftRight: r.norm2
              beta = case r2i != 0.0
                of true: r2ip1/r2i
                of false: 0.0
              verb(3): echo "beta: ", beta
              continuing = (itn < maxits) and (r2r > r2stop)
              if continuing:
                subset: q := z + beta*q
                preconR(p,q)
            of false: # shifted solutions
              var zip1,zip1d,zr: float
              zip1d = alpha*betaim1*(zim1[m]-zi[m]) 
              zip1d += zim1[m]*alphaim1*(1.0+sg[m]*alpha)
              zip1 = case zip1d != 0.0 
                of true: zi[m]*zim1[m]*alphaim1/zip1d
                of false: 0.0
              zr = case zi[m] != 0.0
                of true: zip1/zi[m]
                of false: 0.0
              subset: xs[m] += alpha*zr*ps[m]
              if continuing:
                subset: ps[m] := zip1*r + beta*zr*zr*ps[m]
                (zim1[m],zi[m]) = (zi[m],zip1)
        threadBarrier()
        (alphaim1,betaim1,r2i) = (alpha,beta,r2ip1)
        if threadNum == 0: (r2,itn0) = (r2r,itn)
        threadBarrier()
        verb(3):
          var 
            rip12,diff2,diffs2: float
            rip1 = newOneOf(r)
            diff = newOneOf(r)
            diffs = newOneOf(r)
            Ax = newOneOf(x)
            Asx = newOneOf(x)
          subset:
            op.apply(Ax,x)
            threadBarrier()
            diff := b - Ax
            threadBarrier()
            diff2 = diff.norm2
            echo "rz: ", r2i
            echo "qLAp: ", qLAp
            echo "alpha: ", alpha
            echo "x2: ", x.norm2
            echo "r2: ", r2
            echo "z2: ", z.norm2
            echo "q2: ", q.norm2
            echo "Ap2: ", Ap.norm2
            echo "LAp2: ", LAp.norm2
            echo "|b-Ax|^2: ", diff2
          forMass:
            subset: 
              rip1 := zi[m]*r
              op.apply(Asx,xs[m],sg[m])
              threadBarrier()
              diffs := b - Asx
              threadBarrier()
              rip12 = rip1.norm2
              diffs2 = diffs.norm2
            echo "zim1: ", zim1[m]
            echo "zi: ", zi[m]
            echo "rzs2: ", rip12
            echo "|b-(A+sg)x|^2: ", diffs2
        verb(2): echo "CG iteration: ",itn,"  r2/b2: ",r2/b2
        threadBarrier()
  state.r2 = r2
  state.b2 = b2
  state.xs = xs
  state.ps = ps
  state.iterations = itn0
  sp.finalIterations = state.iterations
  # other things need to be saved?

when isMainModule:
  import cg

  qexInit()

  # For testing purposes
  type
    OpArgs[T] = object 
      m: T
      precon*: CgPrecon

  var
    # Lattice vectors/matrices
    lo = newLayout([4,4,4,4])
    m = lo.ColorMatrix()
    x = lo.ColorVector()
    xr = lo.ColorVector()
    xp = lo.ColorVector()
    xhp = lo.ColorVector()
    xlp = lo.ColorVector()
    b = lo.ColorVector()

    # Mass shifts
    sigmas = @[1.0,2.0,3.0,4.0]

    # Solver parameters
    sp: SolverParams
    
    # No preconditioning - "r" is ref. comparison from cg.nim
    oa = OpArgs[type(m)](m:m,precon:cpNone)
    oar = OpArgs[type(m)](m:m,precon:cpNone)
    
    # "Hermitian" preconditioning
    oah = OpArgs[type(m)](m:m,precon:cpHerm)
    oahr = OpArgs[type(m)](m:m,precon:cpHerm)

    # "Left/right" preconditioning
    oal = OpArgs[type(m)](m:m,precon:cpLeftRight)
    oalr = OpArgs[type(m)](m:m,precon:cpLeftRight)

    # Preconditioning options
    precL = true
    precR = true

  # Procs for applying preconditioning
  proc apply*(oa: OpArgs; r,x: auto) = r := oa.m*x
  proc applyPrecon*(oa: OpArgs; r,x: auto) =
    for e in r:
      let t = sqrt(1.0 / m[e][0,0]) #?
      r[e] := t * x[e]
  proc applyPreconL*(oa: OpArgs; r,x: auto) =
    if precL: applyPrecon(oa, r, x)
    else: r := x
  proc applyPreconR*(oa: OpArgs; r,x: auto) =
    if precR: applyPrecon(oa, r, x)
    else: r := x
  template reset(cgm,cgr: auto) =
    threads:
      x := 0
      xr := 0
    cgm.reset
    cgr.reset
  proc runRestart[T](cgm: var T; oa: auto) =
    sp.maxits = 0
    sp.verbosity = 0
    while cgm.r2 > cgm.r2stop:
      sp.maxits += 10
      cgm.solve(oa,sp)
      echo cgm.iterations, " ", cgm.r2, "/", cgm.r2stop, "... r2/b2: ", cgm.r2/cgm.b2

  # Set solver parameters
  sp.r2req = 1e-20
  sp.maxits = 200
  sp.verbosity = 2
  sp.subset.layoutSubset(lo,"all")

  # Initialize test lattice vectors/matrices
  threads:
    m.even := 1
    m.odd := 10
    threadBarrier()
    tfor i, 0..<lo.nSites: m{i} := i+1
    threadBarrier()
    b.even := 1
    b.odd := 2
    x := 0
    xr := 0

  # Create conjugate gradient objects
  var 
    cgm = newCgmState(x,b,sigmas)
    cgr = newCgState(xr,b)
    cgh = newCgmState(x,b,sigmas,precon=cpHerm)
    cghr = newCgState(xr,b)
    cgl = newCgmState(x,b,sigmas,precon=cpLeftRight)
    cglr = newCgState(xr,b)

  echo "note: shifted vectors tested in stagSolve.nim"
  echo "the following only test the non-shifted vector"

  echo "----------------------------------------------"
  echo "starting cgm.solve base solve test w/o restart"
  echo " ----- cg ----- "
  cgr.solve(oar,sp)
  echo " ----- cgm ----- "
  cgm.solve(oa,sp)
  threads: 
    xp := x
    echo "norm2(cg - cgm): ", (x-xr).norm2
  reset(cgr,cgm)
  
  echo "--------------------------------------------------------"
  echo "start cgm.solve hermitian preconditioning test w/o reset"
  echo " ----- cg ----- "
  cghr.solve(oahr,sp)
  echo " ----- cgm ----- "
  cgh.solve(oah,sp)
  threads: 
    xhp := x
    echo "norm2(cg - cgm): ", (x-xr).norm2
  reset(cghr,cgh)
  
  echo "------------------------------------------------------------"
  echo "starting cgm.solve left/right preconditioning test w/o reset"
  echo " ----- cg ----- "
  cglr.solve(oalr,sp)
  echo " ----- cgm ----- "
  cgl.solve(oal,sp)
  threads: 
    xlp := x
    echo "norm2(cg - cgm): ", (x-xr).norm2
  reset(cglr,cgl)
  
  echo "---------------------------------------------------"
  echo "starting cgm.solve restart test w/o preconditioning"
  cgm.runRestart(oa)
  threads:
    echo "norm2(cg - cgp): ", (x-xp).norm2
    x := 0

  echo "----------------------------------------------------"
  echo "starting cmg.solve restart test w/ hermitian precond."
  cgh.runRestart(oah)
  threads:
    echo "norm2(cg - cgp): ", (x-xhp).norm2
    x := 0

  echo "--------------------------------------------------"
  echo "starting cmg.solve restart test w/ left/right precond."
  cgl.runRestart(oal)
  threads:
    echo "norm2(cg - cgp): ", (x-xlp).norm2
    x := 0

  qexFinalize()
  