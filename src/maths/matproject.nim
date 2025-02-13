# some experimental general projection methods

import base
import matrixConcept, matrixFunctions
import types

proc projectU_newton*(x: var Mat1) =
  mixin simdMax
  let estop = (5*x.nrows*epsilon(x.numberType))^2
  let maxit = 20
  var it = 0
  while true:
    var f = x.adj * x
    let fn = f.norm2
    let ft = f.trace.re
    let w = ft/fn
    let ws = sqrt(w)
    x *= ws
    f *= w
    let r = 1 - f
    let rn = r.norm2.simdMax
    echo it, " ", rn
    if rn < estop or it > maxit: break
    if rn < 1e-16:
      x += 0.5*x*r
    else:
      #let q1 = 1 + r
      let r2 = r*r
      #let q2 = 1 + r2
      #let q4 = 1 + r2*r2
      #let q24 = 1 + r2 + r2*r2
      #let q12 = 1 + r + r2
      #let q3 = 1 + r*r2
      #let t = r
      #let t2 = r + r*r
      #let t4 = t2 + r*t2*r
      #let t6 = t2 + r*t4*r
      #let t8 = t2 + r*t6*r
      #let t = r*q1*q2*q4
      #let t = r*q1*q24
      #let t = r*q12*q3
      #let s1 = r + r2
      #let s2 = r + r*r2
      #let s3 = r + r*r2*r2
      #let t = r + r2 + s1*s2*s3
      let s1 = r + r2
      let s2 = 1 + r2
      let t = s1*s2
      x += 0.5*x*t
    inc it

proc projectSU_newton*(x: var Mat1) =
  mixin simdMax
  let n = x.nrows
  let estop = (5*n*epsilon(x.numberType))^2
  let maxit = 50
  var it = 0
  while true:
    #echo x
    var f = x.adj * x
    #let fn = f.norm2
    #let ft = f.trace.re
    #let w = ft/fn
    var dx = determinant(x)
    #if n>1:
    #  # wf = s^2(1-n) f/dx.norm2 -> s^2(1-n) = wdx.norm2
    #  let s = pow(w*dx.norm2, 1.0/(2*(1-n)))
    #  x *= s
    #  dx *= pow(s, float n)
    let w = pow(dx.norm2, -1.0/n)
    f *= w
    let r = 1 - f
    #echo r
    let rn = r.norm2.simdMax
    echo it, " ", rn
    if rn < estop or it > maxit: break
    if rn < 99999.0:
      let dx2 = dx*dx
      let t = ((1-dx2)/dx2)*x + (1/dx2)*x*r
      x += t
    else:
      #let q1 = 1 + r
      let r2 = r*r
      #let q2 = 1 + r2
      #let q4 = 1 + r2*r2
      #let q24 = 1 + r2 + r2*r2
      #let q12 = 1 + r + r2
      #let q3 = 1 + r*r2
      #let t = r
      #let t2 = r + r*r
      #let t4 = t2 + r*t2*r
      #let t6 = t2 + r*t4*r
      #let t8 = t2 + r*t6*r
      #let t = r*q1*q2*q4
      #let t = r*q1*q24
      #let t = r*q12*q3
      #let s1 = r + r2
      #let s2 = r + r*r2
      #let s3 = r + r*r2*r2
      #let t = r + r2 + s1*s2*s3
      let s1 = r + r2
      let s2 = 1 + r2
      let t = s1*s2
      x += 0.5*x*t
    inc it

proc projectSU_exp*(x: var Mat1; a: Mat1) =
  mixin simdMax
  let n = x.nrows
  let m2stop = epsilon(x.numberType)^2 * a.norm2
  echo "m2stop: ", m2stop
  let maxit = 50
  var it = 0
  while true:
    let ax = a * x.adj
    let xa = x * a.adj
    let m = ax - xa
    #let pm = m
    var pm {.noInit.}: Mat1
    projectTAH(pm, m)
    let pm2 = pm.norm2
    echo "projectSU_exp pm2: ", pm2
    if pm2 < m2stop: break
    let s = ax + xa
    let pmm = -redot(pm,m)
    let pmpm = pm * pm
    let pmpms = redot(pmpm,s)
    let eps = 1e-20
    let c = pmm*pmpms/(pmpms*pmpms+eps)
    echo "projectSU_exp c: ", c
    let ecm = exp(c*pm) * x
    x := ecm
    inc it
    if it >= maxit: break
  echo "projectSU_exp it: ", it

#proc projectU_SVD*(x: var Mat1) =
  # tridiag -> t't + eps -> evecs

when isMainModule:
  import complexNumbers, matrixFunctions, simd

  proc setMat(a: var Mat1) =
    let N = a.nrows
    for i in 0..<N:
      let fi = i.float
      for j in 0..<N:
        let fj = j.float
        let tr = 0.5 + 0.7/(0.9+1.3*fi-fj)
        let ti = 0.1 + 0.3/(0.45+fi-1.1*fj)
        a[i,j].re := tr
        a[i,j].im := ti

  proc testprojU(T: typedesc) =
    var m: T
    setMat(m)
    #echo m
    echo "N: ", m.nrows
    var u = m
    projectU_newton(u)
    let um = u.adj*m
    let umh = (um-um.adj).norm2.simdMaxReduce
    let mu = m*u.adj
    let muh = (mu-mu.adj).norm2.simdMaxReduce
    echo ">> ", (1-m.adj*m).norm2.simdMaxReduce, "  ", (1-u.adj*u).norm2.simdMaxReduce,
        "  ", umh, "  ", muh

  proc testprojSU(T: typedesc) =
    var m: T
    setMat(m)
    #echo m
    echo "N: ", m.nrows
    var u = m
    projectSU_newton(u)
    let um = u.adj*m
    let umh = (um-um.adj).norm2.simdMaxReduce
    let mu = m*u.adj
    let muh = (mu-mu.adj).norm2.simdMaxReduce
    let du = determinant(u)
    let due = (du-1).norm2.simdMaxReduce
    echo ">> ", (1-u.adj*u).norm2.simdMaxReduce, "  ", umh, "  ", muh, "  ", due

  proc testprojSU_exp(T: typedesc) =
    var m: T
    setMat(m)
    #echo m
    echo "N: ", m.nrows
    var u = m
    u := 1
    projectSU_exp(u, m)
    let uu = u.adj * u
    let uu1 = (uu-1).norm2.simdMaxReduce
    let du = determinant(u)
    let due = (du-1).norm2.simdMaxReduce
    #let d = (m-u)
    #let d2 = d.norm2.simdMaxReduce
    let mu = m*u.adj
    var pmu {.noInit.}: T
    projectTAH(pmu, mu)
    let muh = 4.0 * pmu.norm2.simdMaxReduce
    echo ">> ", uu1, "  ", due, "  ", muh

  type
    Cmplx[T] = ComplexType[T]
    CM[N:static[int],T] = MatrixArray[N,N,Cmplx[T]]

  #testprojU(CM[1,SimdD8])
  #testprojU(CM[2,SimdD8])
  testprojU(CM[3,SimdD8])
  #testprojU(CM[4,SimdD8])
  #testprojU(CM[5,SimdD8])

  #testprojU(CM[1,float])
  #testprojU(CM[2,float])
  #testprojU(CM[3,float])
  #testprojU(CM[4,float])
  #testprojU(CM[5,float])
  #testprojU(CM[6,float])
  #testprojU(CM[7,float])
  #testprojU(CM[8,float])
  #testprojU(CM[9,float])
  #testprojU(CM[10,float])

  #testprojSU(CM[1,float])
  #testprojSU(CM[2,float])

  #testprojSU_exp(CM[2,float])

# y'y ~ x'x + eps
# y = x + d
# (x+d)'(x+d) ~ x'x + eps = m
# d'x+x'd+d'd ~ eps
#
# y0 = x
# y1 = a y0 - b 
