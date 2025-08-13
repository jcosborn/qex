import base
import maths/[complexNumbers, matrixConcept, matrixFunctions]
import maths/[types, matinv, matexp, projUderiv]
export matinv
export matexp
getOptimPragmas()

proc eigs3(e0,e1,e2: var auto; tr,p2,det: auto) {.alwaysInline.} =
  mixin sin,cos,acos
  let tr3 = (1.0/3.0)*tr
  let p23 = (1.0/3.0)*p2
  let tr32 = tr3*tr3
  let q = abs(0.5*(p23-tr32))
  let r = 0.25*tr3*(5*tr32-p2) - 0.5*det
  let sq = sqrt(q)
  let sq3 = q*sq
  let isq3 = 1.0/sq3
  var minv,maxv {.noinit.}: type(isq3)
  maxv := 3e38
  minv := -3e38
  let isq3c = min(maxv, max(minv,isq3))
  let rsq3c = r * isq3c
  maxv := 1
  minv := -1
  let rsq3 = min(maxv, max(minv,rsq3c))
  let t = (1.0/3.0)*acos(rsq3)
  let st = sin(t)
  let ct = cos(t)
  let sqc = sq*ct
  let sqs = 1.73205080756887729352*sq*st
  let ll = tr3 + sqc
  e0 = tr3 - 2*sqc
  e1 = ll + sqs
  e2 = ll - sqs

template rsqrtPHM2(r:typed; x:typed) =
  let x00 = x[0,0].re
  let x11 = x[1,1].re
  let x01r = 0.5*(x[0,1].re+x[1,0].re)
  let x01i = 0.5*(x[0,1].im-x[1,0].im)
  let det = abs(x00*x11 - x01r*x01r - x01i*x01i)
  let tr = x00 + x11
  let sdet = sqrt(det)
  let trsdet = tr + sdet
  let c1 = 1/(sdet*sqrt(trsdet+sdet))
  let c0 = trsdet*c1
  r := c0 - c1*x

proc rsqrtPHM3f(c0,c1,c2:var auto; tr,p2,det:auto) {.alwaysInline.} =
  var l0,l1,l2 {.noInit.}: type(tr)
  eigs3(l0,l1,l2, tr,p2,det)
  let sl0 = sqrt(abs(l0))
  let sl1 = sqrt(abs(l1))
  let sl2 = sqrt(abs(l2))
  let u = sl0 + sl1 + sl2
  let w = sl0 * sl1 * sl2
  let d = w*(sl0+sl1)*(sl0+sl2)*(sl1+sl2)
  let di = 1/d
  c0 = (w*u*u+l0*sl0*(l1+l2)+l1*sl1*(l0+l2)+l2*sl2*(l0+l1))*di
  c1 = -(tr*u+w)*di
  c2 = u*di

template rsqrtPHM3(r:typed; x:typed) =
  let tr = trace(x).re
  let x2 = x*x
  let p2 = trace(x2).re
  let det = determinant(x).re
  var c0,c1,c2:type(tr)
  rsqrtPHM3f(c0, c1, c2, tr, p2, det)
  r := c0 + c1*x + c2*x2

template rsqrtPHMN(r:typed; x:typed) =
  mixin simdMax
  let xi = 1/x
  let xi2 = xi.norm2
  let xit = trace(xi).re
  let ds = xit/xi2

  var e = (0.5*ds)*xi - 0.5
  var s = 1 + e

  let estop = epsilon(ds.simdMax)^2
  let maxit = 20
  var nit = 0
  while true:
    inc nit
    let si = 1/s
    let t = e * (si * e)
    e := -0.5 * t
    s += e
    let enorm = e.norm2.simdMax
    if nit>=maxit or enorm<estop: break
  let sds = 1/sqrt(ds)
  r := sds*s

template rsqrtPHM(r:typed; x:typed) =
  mixin rsqrt, nrows
  assert(r.nrows == x.nrows)
  assert(r.ncols == x.ncols)
  assert(r.nrows == r.ncols)
  when r.nrows==1:
    let t = rsqrt(x[0,0].re)
    r := t
  elif r.nrows == 2: rsqrtPHM2(r, x)
  elif r.nrows == 3: rsqrtPHM3(r, x)
  else: rsqrtPHMN(r, x)

proc projectUrsqrtCayleyHamilton(
    r: var Mat1; 
    x: Mat2, 
    eps = 1e-20
  ) {.alwaysInline.} =
  let xa = x.adj
  var t = xa * x
  t += eps
  rsqrtPHM(r, t)

# x (x'x)^{-1/2}
proc projectUCayleyHamilton*(r: var Mat1; x: Mat2, eps = 1e-20) {.inline.} =
  var t2{.noInit.}: evalType(x)
  projectUrsqrtCayleyHamilton(t2, x, eps = eps)
  mul(r, x, t2)

template projectUflops*(nc: int): int =
  nc*nc*(2*(6*nc+2*(nc-1))) + 250 

template projectUCayleyHamilton*(r: var Mat1, eps = 1e-20) =
  var t{.noInit.}: evalType(r)
  t := r
  r.projectUCayleyHamilton t, eps

proc projectUderivCayleyHamilton*(
    r: var Mat1, 
    u: Mat2, 
    x: Mat3, 
    chain: Mat4, 
    eps = 1e-20
  ) =
  var y, z, t1, t2: Mat1
  projectUrsqrtCayleyHamilton(z, x, eps = eps)
  inverse(y, z)
  r := chain * z
  t1 := u.adj * r
  sylsolve(t2, y, t1)
  t1 := t2 + t2.adj
  r -= x * t1

proc projectUderivCayleyHamilton*(r: var Mat1, x: Mat2, c: Mat3) =
  var u {.noInit.}: type(r)
  projectUCayleyHamilton(u, x)
  projectUderivCayleyHamilton(r, u, x, c)