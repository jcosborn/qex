# Numerical derivatives

proc maxAbsElem(z: var auto, x: auto, y: auto) =
  when compiles(len(z)) and compiles(z[0]):
    # array, try element wise max and abs
    for i in 0..<z.len:
      maxAbsElem(z[i], x[i], y[i])
  elif compiles(nrows(z)) and compiles(ncols(z)) and compiles(z[0,0]):
    # a kind of matrix
    for i in 0..<z.nrows:
      for j in 0..<z.ncols:
        maxAbsElem(z[i,j], x[i,j], y[i,j])
  else:
    # hope for the best
    z = max(abs(x), abs(y))

template ndiffTemplate(derivative: untyped) {.dirty.} =
  #[
    Accurate computation of F'(x) and F'(x) F''(x)
    C.J.F. Ridders
    Advances in Engineering Software (1978)
    Volume 4, Issue 2, April 1982, Pages 75-76

    Note: Ridders refers the method as Romberg's, though Romberg's is for integrals.
    It looks like an extension to Neville's algorithm.
  ]#
  let s2 = scale*scale
  var dx = dx
  var A: array[ordMax, F]
  for i in 0..<ordMax:
    A[i] = derivative
    # echo "A[0,",i,"] = ",A[i]
    dx /= scale
  var b = s2
  var c = 1.0/(b-1.0)
  for j in countdown(ordMax-1,2):
    for i in 0..<j:
      var ai = A[i+1]
      ai *= b
      let ai1 = A[i]
      A[i] = (ai - ai1) * c
      # echo "A[",ordMax-j+1,",",i,"] = ",A[i]
    b *= s2
    c = 1.0/(b-1.0)
  let a1 = A[1]
  let a0 = A[0]
  let a1b = a1*b
  var a = a1b - a0
  a *= c
  maxAbsElem(err, a-a0, a-a1)
  r = a

proc ndiff*[X,F](r: var F, err: var F, f: proc(x:X):F, x: X, dx: auto, scale:float=2.0, ordMax:static int=8) =
  ## return r = f'(x), using polynomial approximation with points x+/-dx, with dx/=scale, upto ordMax
  ## err is the difference in the estimate of the current and the previous order.
  ## F can be multidimensional, X is scalar
  ndiffTemplate:
    (f(x+dx) - f(x-dx)) * (0.5/dx)

proc ndiff2*[X,F](r: var F, err: var F, f: proc(x:X):F, x: X, dx: auto, scale:float=2.0, ordMax:static int=8) =
  ## return r = f''(x), using polynomial approximation with points x+/-dx, with dx/=scale, upto ordMax
  ## err is the difference in the estimate of the current and the previous order.
  ## F can be multidimensional, X is scalar
  ndiffTemplate:
    let idx = 1.0/dx
    (f(x+dx) - 2.0*f(x) + f(x-dx)) * (idx*idx)

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

  echo "dx = 0.1, scale=5.0, ordMax=5"
  ndiff(r, e, f, 1.0, 0.1, scale=5.0, ordMax=5)
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

  import simd
  var v, vr, ve: SimdD4
  v[0] = 1.0
  v[1] = -0.1
  v[2] = -1.0
  v[3] = 10.0
  proc f(x:SimdD4):SimdD4 = exp(x)/(sin(x)-x*x)
  ndiff(vr, ve, f , v, 0.05, scale=5.0, ordMax=5)
  let vex = exp(v)
  let vd = sin(v)-v*v
  let va = vex/vd*(1.0 - (cos(v)-2.0*v)/vd)
  echo "Analytic: ", va
  echo "Numeric:  ", vr, "  error: ",vr-va,"  estimated: ",ve

  ndiff2(vr, ve, f , v, 0.05, scale=3.5, ordMax=6)
  let vc = cos(v)-2.0*v
  let va2 = vex/vd*(1.0 - (2.0*vc - sin(v) - 2.0 - 2.0*vc*vc/vd)/vd)
  echo "2nd derivative:"
  echo "Analytic: ", va2
  echo "Numeric:  ", vr, "  error: ",vr-va2,"  estimated: ",ve
