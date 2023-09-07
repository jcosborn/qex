import base
import matrixConcept
import matrixFunctions
import types

type
  RecipRootParam* = object
    root*: int
    maxit*: int
    tol*: float
    order*: int

# https://arxiv.org/pdf/1703.02456.pdf
proc rsqrtN2(r: var auto; a: auto) =
  mixin sqrt
  let an = a.norm2
  r := 1

  let estop = (20*a.nrows*epsilon(an.simdMax))^2
  let maxit = 50
  var nit = 0
  var done = false
  while true:
    let r2 = r*r
    let f = a*r2
    let fn = f.norm2
    let ft = f.trace.adj
    var w = sqrt(ft/fn)
    let wn = sqrt(w.norm2)
    var c = 1/sqrt(r.norm2*sqrt(an))
    if (f-1).norm2 < 0.01:
      c *= 0.1
    w += c - min(wn,c)
    if (f-1).norm2 > 0.01:
      r *= w
    let g = 1 - (w*w)*f
    let e = 1 - (w*w)*r2*a
    let enorm = e.norm2.simdMax
    #echo nit, " enorm: ", enorm, "   ", (r*a-a*r).norm2
    echo nit, " enorm: ", enorm, "   ", g.norm2, "   ", sqrt(a.norm2*r2.norm2)
    #if nit>=maxit or enorm<estop: break
    if nit>=maxit or done: break
    if enorm<estop: done = true
    inc nit
    if enorm < 0.01:
      #let t = e*r
      #let t = 0.5*(r*e+e*r)
      let t = 0.5*(e*r+r*g)
      r += 0.5*t
    else:
      let ee = e*(1+e+e*(e+e*(e+e*e)))
      let t = ee*r
      #let t = 0.5*(r*ee+ee*r)
      r += 0.5*t
  echo r

when isMainModule:
  import complexNumbers

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
    if a.nrows > 1:
      let ta = a.trace/a.nrows
      a -= ta
    #a *= 100
    a *= 1000
    #a *= 0.01

  proc testrsqrt(T: typedesc) =
    var m, x: T
    setMat(m)
    echo m
    echo m.determinant
    rsqrtN2(x, m)
    let e = (x*m*x-1).norm2
    echo e
    #let s = m*x
    let s = 0.5*(m*x+x*m)
    let d = (s*s-m).norm2/m.norm2
    echo d

  type
    Cmplx[T] = ComplexType[T]
    CM[N:static[int],T] = MatrixArray[N,N,Cmplx[T]]

  testrsqrt(CM[1,float])
  testrsqrt(CM[2,float])
  testrsqrt(CM[3,float])
  testrsqrt(CM[4,float])
  testrsqrt(CM[5,float])
  testrsqrt(CM[6,float])
  testrsqrt(CM[7,float])


    #[
    for i in 1..10:
      let cra = r*a-a*r
      let acra = a.adj*cra-cra*a.adj
      let pa = acra*a-a*acra
      let alp = -dot(pa,cra)/pa.norm2
      r += alp*acra
      echo cra.norm2, " ", (r*a-a*r).norm2, " ", acra.norm2
    if enorm < 1e-12:
      for i in 1..50:
        #let cra = r*a-a*r
        #let acra = a.adj*cra-cra*a.adj
        #let pa = acra*a-a*acra

        let rr = 1 - r*a*r
        let p = a.adj*r.adj*rr + rr*r.adj*a.adj #+ acra
        let xp = r*a*p + p*a*r
        let b0 = redot(xp,rr) #- redot(pa,cra)
        let pp = p*a*p
        let b1 = xp.norm2 - 2*redot(pp,rr) #+ pa.norm2
        let alp = b0/b1
        r += alp*p
        let rr2 = 1 - r*a*r
        echo i, " ", rr.norm2, " ", rr2.norm2
    ]#

