## Brief: Packaged objects/methods for unitary (and special unitary) projection
## 
## Note that only objects for unitary projection have been defined. Objects for
## special unitary projection can be defined if/when the need arises.
## 
## Author: James C. Osborn 
## Author: Curtis Taylor Peterson

import base
import complexNumbers
import maths/[matrixConcept, matrixFunctions]
import maths/[matinv, projUderiv]
import simd

getOptimPragmas()

type
  ProjectionMethod = enum
    CayleyHamilton,
    GoloubKahanReinsch, # not yet in devel branch
    Newton,
    Halley,
    Exponential # only SU

type
  UnitaryProjection = object
    policy: ProjectionMethod
    eps: float
    maxiters: int

proc newUnitaryProjection(
    policy: ProjectionMethod = CayleyHamilton, 
    eps: float = 1e-16, 
    maxiters: int = 100
  ): UnitaryProjection =
  ## Brief: UnitaryProjection constructor
  ## Author: Curtis Taylor Peterson
  ## 
  ## Parameters:
  ##   policy  [ProjectionMethod]: unitary projection method
  ##   eps     [float]: possible epsilon terminating condition
  ##   maxitns [int]: possible max iteration number before terminating
  result = UnitaryProjection(policy: policy)
  case policy:
    of CayleyHamilton: discard
    of GoloubKahanReinsch:
      qexError("\"GoloubKahanReinsch\" not implemented for unitary projection") 
    of Newton: result.maxiters = maxiters
    of Halley: result.maxiters = maxiters
    of Exponential: 
      qexError("\"Exponential\" not implemented for unitary projection") 

proc projectUNewton(x: var Mat1, maxiters: int) =
  ## Brief: Polar-decomposition-based unitary projection (Newton's method)
  ## Author: James C. Osborn
  ## 
  ## Parameters:
  ##   x [MatrixArray[N,N,T]]: matrix to be projected & output of projection
  mixin simdMax
  let estop = (5*x.nrows*epsilon(x.numberType))^2
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
    if rn < estop or it > maxiters: break
    if rn < 1e-16:
      x += 0.5*x*r
    else:
      let r2 = r*r
      let s1 = r + r2
      let s2 = 1 + r2
      let t = s1*s2
      x += 0.5*x*t
    inc it

proc projectSUNewton(x: var Mat1, maxiters: int) =
  ## Brief: Polar-decomposition-based special unitary projection (Newton's method)
  ## 
  ## Parameters:
  ##   x        [MatrixArray[N,N,T]]: matrix to be projected & output of projection
  ##   maxiters [int]: Number of Halley method iterations
  mixin simdMax
  let n = x.nrows
  let estop = (5*n*epsilon(x.numberType))^2
  var it = 0
  while true:
    var f = x.adj * x
    var dx = determinant(x)
    let w = pow(dx.norm2, -1.0/n)
    f *= w
    let r = 1 - f
    let rn = r.norm2.simdMax
    echo it, " ", rn
    if rn < estop or it > maxiters: break
    if rn < 99999.0:
      let dx2 = dx*dx
      let t = ((1-dx2)/dx2)*x + (1/dx2)*x*r
      x += t
    else:
      let r2 = r*r
      let s1 = r + r2
      let s2 = 1 + r2
      let t = s1*s2
      x += 0.5*x*t
    inc it

proc projectUHalley(u: var Mat1; maxiters: int) =
  ## Brief: Polar-decomposition-based unitary projection (Halley's method)
  ## Author: Curtis Taylor Peterson
  ## 
  ## Parameters:
  ##   u        [MatrixArray[N,N,T]]: matrix to be projected & output of projection
  ##   maxiters [int]: Number of Halley method iterations
  var 
    iter = 0
    r {.noinit.}: u.numberType
    udu, ta, tb {.noinit.}: evalType(u)
  let estop = (5*u.nrows*epsilon(u.numberType))^2

  while true:
    if iter > maxiters: break
    udu := u.adj*u
    r := (1 - udu).norm2.simdMax
    if r < estop: break
    tb := 1 + 3.0*udu
    ta.inverse(tb)
    tb := 3 + udu
    u := u*ta*tb
    iter.inc

proc projectSUExponential(x: var Mat1; a: Mat1; maxiters: int) =
  ## Brief: Polar-decomposition-based special unitary projection (exponential)
  ## Author: James C. Osborn
  ## 
  ## Parameters:
  ##   x [MatrixArray[N,N,T]]: output of projection
  ##   a [MatrixArray[N,N,T]]: matrix to be projected    
  mixin simdMax
  let n = x.nrows
  let m2stop = epsilon(x.numberType)^2 * a.norm2
  echo "m2stop: ", m2stop
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
    if it >= maxiters: break
  echo "projectSU_exp it: ", it

proc projectU(self: UnitaryProjection, v: var Mat1; u: Mat1) =
  ## Brief: Performs unitary projection
  ## Author: Curtis Taylor Peterson
  ## 
  ## Parameters:
  ##   self  [UnitaryProjection]: unitary projection object
  ##   v     [MatrixArray[N,N,T]]: output of projection
  ##   u     [MatrixArray[N,N,T]]: matrix to be projected
  v := u
  case self.policy:
    of CayleyHamilton: v.projectU()
    of GoloubKahanReinsch: discard
    of Newton: v.projectUNewton(self.maxiters)
    of Halley: v.projectUHalley(self.maxiters)
    of Exponential: discard

proc projectUderivGeneral(r: var Mat1, u: Mat2, x: Mat3, chain: Mat4) =
  ## Brief: derivative of unitary projection
  ## Author: Curtis Taylor Peterson
  ## 
  ## Parameters:
  ##   v     [MatrixArray[N,N,T]]: result of projection derivative
  ##   u     [MatrixArray[N,N,T]]: unitary projection of x
  ##   x     [MatrixArray[N,N,T]]: matrix that has been unitarily projected
  ##   chain [MatrixArray[N,N,T]]: backprop accumulation
  var y, z, t1, t2 {.noinit.}: Mat1
  t1.inverse(x)
  z := t1*u
  y.inverse(z)
  r := chain*z
  t1 := u.adj*r
  sylsolve(t2, y, t1)
  t1 := t2 + t2.adj
  r -= x*t1

proc projectUderiv(
    self: UnitaryProjection, 
    v: var Mat1, 
    u: Mat2, 
    x: Mat3,
    chain: Mat4
  ) =
  ## Brief: derivative of unitary projection
  ## Author: Curtis Taylor Peterson
  ## 
  ## Parameters:
  ##   self  [UnitaryProjection]: unitary projection object
  ##   v     [MatrixArray[N,N,T]]: result of projection derivative
  ##   u     [MatrixArray[N,N,T]]: unitary projection of x
  ##   x     [MatrixArray[N,N,T]]: matrix that has been unitarily projected
  ##   chain [MatrixArray[N,N,T]]: backprop accumulation
  var chain: typeOf(x)
  chain := 1
  case self.policy:
    of CayleyHamilton: v.projectUderiv(u, x, chain)
    of GoloubKahanReinsch: discard
    of Newton: v.projectUderivGeneral(u, x, chain)
    of Halley: v.projectUderivGeneral(u, x, chain)
    of Exponential: discard

proc projectUderiv(
    self: UnitaryProjection, 
    v: var Mat1, 
    x: Mat3, 
    chain: Mat4
  ) =
  ## Brief: derivative of unitary projection
  ## Author: Curtis Taylor Peterson
  ## 
  ## Parameters:
  ##   self  [UnitaryProjection]: unitary projection object
  ##   v     [MatrixArray[N,N,T]]: result of projection derivative
  ##   x     [MatrixArray[N,N,T]]: matrix that has been unitarily projected
  ##   chain [MatrixArray[N,N,T]]: backprop accumulation
  var u {.noinit.}: typeOf(v)
  self.projectU(u, x)
  self.projectUderiv(v, u, x, chain)

proc projectUderiv(self: UnitaryProjection, v: var Mat1, x: Mat2) =
  ## Brief: derivative of unitary projection
  ## Author: Curtis Taylor Peterson
  ## 
  ## Parameters:
  ##   self  [UnitaryProjection]: unitary projection object
  ##   v     [MatrixArray[N,N,T]]: result of projection derivative
  ##   x     [MatrixArray[N,N,T]]: matrix that has been unitarily projected
  var u, chain {.noinit.}: typeOf(v)
  self.projectU(u, x)
  chain := 1
  self.projectUderiv(v, u, x, chain)

when isMainModule:
  # observation: Newton/Halley methods usually less expensive and more 
  # accurate than Cayley-Hamilton method; they are also both likely to 
  # be far less expensive than Goloub-Kahan-Reinsch (which is not SIMD-
  # friendly anyway)

  import macros
  import times

  type
    Cmplx[T] = ComplexType[T]
    CM[N:static[int],T] = MatrixArray[N,N,Cmplx[T]]

  template `+`(x: SimdS4): untyped = x
  template `+`(x: SimdS4, y: ComplexType): untyped =
    asReal(x) + y
  template `-`(x: SimdS4, y: ComplexType): untyped =
    asReal(x) - y
  template `*`(x: SimdS4, y: ComplexType): untyped =
    asReal(x) * y
  template `*`(x: ComplexType, y: SimdS4): untyped =
    x * asReal(y)
  template `/`(x: SomeFloat, y: SomeInteger): untyped = x/(type(x))(y)
  template add(r: ComplexType, x: SimdS4, y: ComplexType): untyped =
    add(r, asReal(x), y)
  template sub(r: ComplexType, x: SimdS4, y: ComplexType): untyped =
    sub(r, asReal(x), y)
  template mul(r: ComplexType, x: SimdS4, y: ComplexType): untyped =
    mul(r, asReal(x), y)

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

  template check(m:string; x:untyped; n:SomeNumber):untyped =
    let r0 = x
    let r = simdSum(r0)/simdLength(r0)
    echo m & " error/eps: ", r/epsilon(r)

  proc testprojU(T: typedesc) =
    let 
      default = newUnitaryProjection()
      halley = newUnitaryProjection(policy = Halley, maxiters = 100)
      newton = newUnitaryProjection(policy = Newton, maxiters = 100)
      methods = @[default, halley, newton]
    var m: T
    setMat(m)
    let N = m.nrows
    echo "----------------------------------------------------"
    echo "test: " & $N & " " & $T
    for idx in 0..<methods.len:
      echo "\n~~~~~~~~~~~~~~~~~~~~~~~~"
      echo $methods[idx].policy & " projection"
      echo "~~~~~~~~~~~~~~~~~~~~~~~~"
      var u = m
      var m0 = m
      methods[idx].projectU(u,m0)
      let um = u.adj*m
      let umh = (um-um.adj).norm2.simdMaxReduce
      let mu = m*u.adj
      let muh = (mu-mu.adj).norm2.simdMaxReduce
      echo "unitarity error before: ", (1-m.adj*m).norm2.simdMaxReduce 
      echo "unitarity error after: ", (1-u.adj*u).norm2.simdMaxReduce
      echo "hermiticity error: ", umh, "  ", muh
    echo "----------------------------------------------------"

  proc testprojSU(T: typedesc) =
    var m: T
    setMat(m)
    echo "N: ", m.nrows
    var u = m
    projectSUNewton(u, 20)
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
    echo "N: ", m.nrows
    var u = m
    u := 1
    projectSUExponential(u, m, 20)
    let uu = u.adj * u
    let uu1 = (uu-1).norm2.simdMaxReduce
    let du = determinant(u)
    let due = (du-1).norm2.simdMaxReduce
    let mu = m*u.adj
    var pmu {.noInit.}: T
    projectTAH(pmu, mu)
    let muh = 4.0 * pmu.norm2.simdMaxReduce
    echo ">> ", uu1, "  ", due, "  ", muh

  proc test(T: typedesc) =
    var t0: float
    var m,mCH,mPD,mN: T
    var dCH,dPD,dN: T
    let N = m.nrows

    proc ticc() = 
      t0 = cpuTime()
    
    proc tocc(message: string) =
      echo message, " ", cpuTime() - t0, " s"
      t0 = cpuTime()

    setMat(m)

    echo "----------------------------------------------------"
    echo "test: " & $N & " " & $T

    echo "\n~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "projection"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~"
    let eps = epsilon(m.norm2.simdSum)
    let 
      default = newUnitaryProjection()
      halley = newUnitaryProjection(policy = Halley, maxiters = 100)
      newton = newUnitaryProjection(policy = Newton, maxiters = 100)
    ticc()
    default.projectU(mCH,m)
    tocc("Cayley-Hamilton:")
    halley.projectU(mPD,m)
    tocc("Halley:")
    newton.projectU(mN,m)
    tocc("Newton:")
    let err = sqrt((mCH-mPD).norm2)/N
    let err2 = sqrt((mCH-mN).norm2)/N
    let errCH = sqrt((1-mCH*mCH.adj).norm2)/N
    let errPD = sqrt((1-mPD*mPD.adj).norm2)/N
    let errN = sqrt((1-mN*mN.adj).norm2)/N
    echo "Cayley-Hamilton - Halley: ", err
    echo "Cayley-Hamilton - Newton: ", err2
    echo "Haley err / Cayley-Hamilton err: ", errPD/errCH
    echo "Newton err / Cayley-Hamilton err: ", errN/errCH
    check("Cayley-Hamilton:", errCH, 50)
    check("Halley:", errPD, 50)
    check("Newton:", errN, 50)

    echo "\n~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "projection derivative (1)"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~"
    ticc()
    default.projectUderiv(dCH,m)
    tocc("Cayley-Hamilton derivative:")
    halley.projectUderiv(dPD,m)
    tocc("Halley derivative:")
    newton.projectUderiv(dN,m)
    tocc("Newton derivative:")
    var errdH = sqrt((dCH-dPD).norm2)/N
    var errdN = sqrt((dCH-dN).norm2)/N
    echo "Cayley-Hamilton - Halley: ", errdH
    echo "Cayley-Hamilton - Newton: ", errdN

    echo "\n~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "projection derivative (2)"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~"
    ticc()
    var chain = dCH
    default.projectUderiv(dCH,m,chain)
    tocc("Cayley-Hamilton derivative:")
    halley.projectUderiv(dPD,m,chain)
    tocc("Halley derivative:")
    newton.projectUderiv(dN,m,chain)
    tocc("Newton derivative:")
    errdH := sqrt((dCH-dPD).norm2)/N
    errdN := sqrt((dCH-dN).norm2)/N
    echo "Cayley-Hamilton - Halley: ", errdH
    echo "Cayley-Hamilton - Newton: ", errdN

    echo "\n~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "projection derivative (3)"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~"
    ticc()
    default.projectUderiv(dCH,mCH,m,chain)
    tocc("Cayley-Hamilton derivative:")
    halley.projectUderiv(dPD,mPD,m,chain)
    tocc("Halley derivative:")
    newton.projectUderiv(dN,mN,m,chain)
    tocc("Newton derivative:")
    errdH := sqrt((dCH-dPD).norm2)/N
    errdN := sqrt((dCH-dN).norm2)/N
    echo "Cayley-Hamilton - Halley: ", errdH
    echo "Cayley-Hamilton - Newton: ", errdN
    echo "----------------------------------------------------"

  template doTest(t:untyped) =
    when declared(t):
      test(CM[1,t])
      testprojU(CM[1,t])
      test(CM[2,t])
      testprojU(CM[2,t])
      test(CM[3,t])
      testprojU(CM[3,t])
      test(CM[4,t])
      testprojU(CM[4,t])

  doTest(float32)
  doTest(float64)
  doTest(SimdS4)
  doTest(SimdD4)
  doTest(SimdS8)
  doTest(SimdD8)

#[
# y'y ~ x'x + eps
# y = x + d
# (x+d)'(x+d) ~ x'x + eps = m
# d'x+x'd+d'd ~ eps
#
# y0 = x
# y1 = a y0 - b 
]#