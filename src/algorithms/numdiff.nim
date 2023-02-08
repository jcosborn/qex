# Numerical derivatives

type UpTri[N:static int,T] = object
  d: array[(N*(N+1) div 2),T]

func upTriIndex[N:static int](i,j:int): int =
  # 0<=i<=j<N
  ((N+N+1-i)*i div 2) + j - i

proc `[]`[N,T](x: UpTri[N,T], i,j: int): T = x.d[upTriIndex[N](i,j)]
proc `[]`[N,T](x: var UpTri[N,T], i,j: int): var T = x.d[upTriIndex[N](i,j)]
proc `[]=`[N,T](x: var UpTri[N,T], i,j: int, y: T) = x.d[upTriIndex[N](i,j)] = y

template ndiffTemplate(derivative: untyped) {.dirty.} =
  #[
    Accurate computation of F'(x) and F'(x) F''(x)
    C.J.F. Ridders
    Advances in Engineering Software (1978)
    Volume 4, Issue 2, April 1982, Pages 75-76

    Note: Ridders refer the method as Romberg's, though Romberg's is for integrals.
    It looks like an extension to Neville's algorithm.
  ]#
  let s2 = scale*scale
  var dx = dx
  var A: UpTri[ordMax, F]
  for i in 0..<ordMax:
    A[0,i] = (f(x+dx) - f(x-dx)) / (2*dx)
    # echo "A[0,",i,"] = ",A[0,i]
    dx /= scale
    var b = s2
    for j in 1..i:
      A[j,i] = (b*A[j-1,i] - A[j-1,i-1]) / (b - 1.0)
      # echo "A[",j,",",i,"] = ",A[j,i]
      b *= s2
  let a = A[ordMax-1,ordMax-1]
  err = max(abs(a-A[ordMax-2,ordMax-1]), abs(a-A[ordMax-2,ordMax-2]))
  r = a

proc ndiff*[X,F](r: var F, err: var F, f: proc(x:X):F, x: X, dx: auto, scale:float =1.618, ordMax:static int=8) =
  ## return r = f'(x), using polynomial approximation with points x+/-dx, with dx/=scale, upto ordMax
  ## err is the difference in the estimate of the current and the previous order.
  ndiffTemplate:
    (f(x+dx) - f(x-dx)) / (2*dx)

proc ndiff2*[X,F](r: var F, err: var F, f: proc(x:X):F, x: X, dx: auto, scale:float=2.0, ordMax:static int=8) =
  ## return r = f''(x), using polynomial approximation with points x+/-dx, with dx/=scale, upto ordMax
  ## err is the difference in the estimate of the current and the previous order.
  ndiffTemplate:
    (f(x+dx) - 2.0*f(x) + f(x-dx)) / (dx*dx)

when isMainModule:
  import math
  # example in Ridder's paper
  proc f(x:float):float = exp(x)/(sin(x)-x*x)
  var r,e:float
  ndiff(r, e, f, 1.0, 0.01, ordMax=5)
  let ex = exp(1.0)
  let d = sin(1.0)-1.0
  let a = ex/d - ex*(cos(1.0)-2.0)/(d*d)
  echo "Analytic: ", a
  echo "Numeric:  ", r, "  error: ",r-a,"  estimated: ",e

  echo "ordMax = 8"
  ndiff(r, e, f, 1.0, 0.01, ordMax=8)
  echo "Numeric:  ", r, "  error: ",r-a,"  estimated: ",e

  # example 2, exp(x)
  echo "exp(1)"
  ndiff2(r, e, exp, 1.0, 1.0, ordMax=5)
  echo "Analytic: ", exp(1.0)
  echo "Numeric:  ", r, "  error: ",r-exp(1.0),"  estimated: ",e
  echo "ordMax = 8"
  ndiff2(r, e, exp, 1.0, 1.0, ordMax=8)
  echo "Numeric:  ", r, "  error: ",r-exp(1.0),"  estimated: ",e

  echo "exp(10)"
  ndiff2(r, e, exp, 10.0, 1.0, ordMax=5)
  echo "Analytic: ", exp(10.0)
  echo "Numeric:  ", r, "  error: ",r-exp(10.0),"  estimated: ",e
  echo "ordMax = 8"
  ndiff2(r, e, exp, 10.0, 1.0, ordMax=8)
  echo "Numeric:  ", r, "  error: ",r-exp(10.0),"  estimated: ",e

#[
  import simd
  var v, vr, ve: SimdD4
  v[0] = 1.0
  v[1] = -0.1
  v[2] = -1.0
  v[3] = 10.0
  proc f(x:SimdD4):SimdD4 = exp(x)/(sin(x)-x*x)
  ndiff(vr, ve, f , v, 0.01)
  let vex = exp(v)
  let vd = sin(v)-v*v
  let va = ex/d - ex*(cos(v)-2.0*v*v)/(d*d)
  echo "Analytic: ", va
  echo "Numeric:  ", vr, "  error: ",vr-va,"  estimated: ",ve
]#
