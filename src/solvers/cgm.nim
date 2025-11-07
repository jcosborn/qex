## Brief: Multi-shift solver
## 
## Author: Curtis Taylor Peterson <curtistaylorpetersonwork@gmail.com>
## 
## Details: Multi-mass (multi-shift) conjugate gradient
## Based on "Krylov space solvers for shifted linear systems" (arXiv:9612014) 
## and "Reduced-Shifted Conjugate-Gradient Method for a Green's Function: 
## Efficient Numerical Approach in a Nano-structured Superconductor"
## (arXiv:1607.03992v2).

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
    r,b,Ap: T
    xs,ps: seq[T]
    sigmas*: seq[float]
    b2,r2*,r2stop*: float
    iterations: int

# Reset multi-shift CG state
proc reset*(cgs: var CgmState; recycle: bool) =
  cgs.iterations = 0
  cgs.b2 = case recycle 
    of true: 1
    of false: -1
  cgs.r2 = 1.0
  cgs.r2stop = 0.0

# Construct new multi-shift CG state
proc createCgmState*[T](
    xs: seq[T];
    b: T; 
    sigmas: seq[float];
    recycle: bool
  ): CgmState[T] =
  let nmass = sigmas.len
  result = CgmState[T]()
  result.sigmas = sigmas
  result.xs = xs
  result.ps = newSeq[T]()
  for m in 0..<nmass: result.ps.add newOneOf(xs[m])
  result.reset(recycle)

# Construct new multi-shift CG state
proc newCgmState*[T](
    xs: seq[T];
    b: T; 
    sigmas: seq[float];
    recycle: bool
  ): CgmState[T] = 
  result = createCgmState(xs,b,sigmas,recycle)
  result.b = b
  result.r = newOneOf(b)
  result.Ap = newOneOf(b)

# Multi-shift CG solver (arXiv:9612014,1607.03992v2)
proc solve*[V,T](
    state: var CgmState; 
    op: proc(a,b: Field[V,T]; shift: float = 0.0); 
    sp: var SolverParams
  ) = 
  ## Brief: multi-shift ("multi-mass") solver
  ## Author: Curtis Taylor Peterson
  ## 
  ## Details:
  ##   Implements multi-shift conjugate gradient a la arXiv:9612014
  ##   and arXiv:1607.03992v2; see for details
  ## 
  ## Input:
  ##   state [CgmState]: conjugate gradient state object
  ##   op    [proc(...)]: proc implementing linear operator
  ##   cp    [SolverParams]: linear solver state object
  var 
    r2,r2stop: float
    b2 = state.b2
    itn0 = 0
    maxits = sp.maxits
    sub = sp.subset
    restart = false
    sg = newSeq[float]()
  let 
    nmass = state.xs.len
    vrb = sp.verbosity
    r = state.r
    ps = state.ps
    Ap = state.Ap
    xs = state.xs
    b = state.b

  # helper templates
  template verb(n: int; body: untyped) = (if vrb >= n: body)

  template mythreads(body: untyped) =
    threads:
      onNoSync(sub): body

  template subset(body: untyped) =
    onNoSync(sub): body

  template x: untyped = xs[0]
  template p: untyped = ps[0]

  # preparation
  for m in 0..<nmass: 
    case m == 0:
      of true: sg.add 0.0
      of false: sg.add state.sigmas[m]

  mythreads:
    var r2t,b2t: float
    case b2 < 0.0: # new solution
      of true:
        r := b
        for m in 0..<nmass: xs[m] := 0.0
        b2t = b.norm2
        r2t = b2t 
        verb(1): echo "input norm2: ", b2t
      of false: # CG restart
        op(Ap,xs[0])
        r := b - Ap
        threadBarrier()
        (b2t,r2t) = (b2,r.norm2)
    verb(3): 
      echo "p2: ", ps[0].norm2
      echo "r2: ", r2t
    threadMaster:
      r2 = r2t
      b2 = b2t
  r2stop = sp.r2req * b2
  state.r2stop = r2stop

  # solve
  if r2 > r2stop:
    threads:
      # extra variable initialization
      var
        r2r,r2i,r2ip1,pAp: float
        alpha,beta: float
        alphaim1,betaim1: float
        (zi,zim1) = (newSeq[float](nmass),newSeq[float](nmass))
        (itn,continuing) = (0,true)

      # iteration 0
      r2i = r2
      (alphaim1,betaim1) = (-1.0,0.0)
      for m in 0..<nmass: 
        (zim1[m],zi[m]) = (1.0,1.0)
        subset: ps[m] := r
      verb(1): echo "CG iteration: ",itn,"  r2/b2: ",r2/b2
      
      # iteration 1,2,3,...
      while continuing:
        # iteration i
        for m in 0..<nmass:
          threadBarrier()
          case m == 0:
            of true: # "sigma = 0"
              subset: op(Ap,p)
              inc itn
              subset: pAp = redot(p,Ap)
              alpha = case pAp != 0.0
                of true: r2i/pAp
                of false: 0.0
              subset: 
                r -= alpha*Ap
                x += alpha*p
              threadBarrier()
              r2ip1 = r.norm2
              r2r = r2ip1
              beta = case r2i != 0.0
                of true: r2ip1/r2i
                of false: 0.0
              verb(3): echo "beta: ", beta
              continuing = (itn < maxits) and (r2r > r2stop)
              if continuing: p := r + beta*p
            of false: # "sigma != 0" (i.e., shifted)
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

        # optional printout
        verb(3):
          var 
            rip12,diff2,diffs2: float
            rip1 = newOneOf(r)
            diff = newOneOf(r)
            diffs = newOneOf(r)
            Ax = newOneOf(x)
            Asx = newOneOf(x)
          subset:
            op(Ax,x)
            threadBarrier()
            diff := b - Ax
            threadBarrier()
            diff2 = diff.norm2
            echo "rz: ", r2i
            echo "pAp: ", pAp
            echo "alpha: ", alpha
            echo "x2: ", x.norm2
            echo "r2: ", r2
            echo "Ap2: ", Ap.norm2
            echo "|b-Ax|^2: ", diff2
          for m in 0..<nmass:
            subset: 
              rip1 := zi[m]*r
              op(Asx,xs[m],sg[m])
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
  state.r2 = r2
  state.b2 = b2
  state.xs = xs
  state.ps = ps
  state.iterations = itn0
  sp.finalIterations = state.iterations

if isMainModule: echo "!!! NO TEST: Run stagSolve for test !!!"
