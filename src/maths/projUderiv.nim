import base
import complexNumbers
import matrixConcept
import types
import strformat
getOptimPragmas()

proc adjugate*(r: var Mat1, x: Mat2) {.alwaysInline.} =
  const nc = r.nrows
  when nc==1:
    r := 1
  elif nc==2:
    r[0,0] :=  x[1,1]
    r[0,1] := -x[0,1]
    r[1,0] := -x[1,0]
    r[1,1] :=  x[0,0]
  elif nc==3:
    let x00 = x[0,0]
    let x01 = x[0,1]
    let x02 = x[0,2]
    let x10 = x[1,0]
    let x11 = x[1,1]
    let x12 = x[1,2]
    let x20 = x[2,0]
    let x21 = x[2,1]
    let x22 = x[2,2]
    r[0,0] := x11*x22 - x12*x21
    r[0,1] := x21*x02 - x22*x01
    r[0,2] := x01*x12 - x02*x11
    r[1,0] := x12*x20 - x10*x22
    r[1,1] := x22*x00 - x20*x02
    r[1,2] := x02*x10 - x00*x12
    r[2,0] := x10*x21 - x11*x20
    r[2,1] := x20*x01 - x21*x00
    r[2,2] := x00*x11 - x01*x10
  else:
    echo &"adjugate n({nc})>3 not supported"
    doAssert(false)

proc sylsolveN*(x: var Mat1, a0: Mat2, c0: Mat3) =
  mixin simdMax
  let nc = x.nrows
  let a2 = a0.norm2
  let ia = rsqrt(a2)

  x := (0.5*ia)*c0
  var a = ia*a0

  let rstop = epsilon(x.norm2.simdMax)
  let maxit = 20
  var nit = 0
  while true:
    #echo nit, ": ", x
    inc nit

    var v = 3 - a*a
    var w = a*x - x*a
    var d = x*v
    x := 0.5*(d-a*w)
    a := 0.5*a*v

    let r = c0-(a0*x+x*a0)
    let rnorm = r.norm2.simdMax
    #echo nit, " r2: ", rnorm
    if nit>=maxit or rnorm<rstop:
      if rnorm>rstop:
        echo "WARNING sylsolveN failed to converge: ", nit, " r2: ", rnorm
      break

proc sylsolveN2*(x: var Mat1, a: Mat2, c: Mat3) =
  mixin simdMax
  let nc = x.nrows
  let aa = a.adj
  let t2 = a.norm2
  let kappa = 0.5/t2

  x := c
  #x := kappa * (aa*c + c*aa)

  let rstop = epsilon(x.norm2.simdMax)
  let maxit = 50
  var nit = 0
  while true:
    inc nit

    let r = c - a*x - x*a
    x += kappa * (aa*r + r*aa)

    let rnorm = r.norm2.simdMax
    echo nit, " r2: ", rnorm
    if nit>=maxit or rnorm<rstop:
      if rnorm>rstop:
        echo "WARNING sylsolveN failed to converge: ", nit, " r2: ", rnorm
      break

proc sylsolve*(x: var Mat1, a: Mat2, c: Mat3) =
  ## solves A X + X A = C for X
  const nc = x.nrows
  when nc==1:
    x[0,0] := c[0,0] / (2*a[0,0])
  elif nc==2:
    # x = (C + |A| A^-1 C A^-1)/2Tr(A)
    let a00 = a[0,0]
    let a01 = a[0,1]
    let a10 = a[1,0]
    let a11 = a[1,1]
    let c00 = c[0,0]
    let c01 = c[0,1]
    let c10 = c[1,0]
    let c11 = c[1,1]
    let idet = 1/(a00*a11 - a01*a10)
    let itr = 0.5/(a00 + a11)
    # ai = [[a11,-a01][-a10,a00]]
    let aic00 = a11 * c00 - a01 * c10
    let aic01 = a11 * c01 - a01 * c11
    let aic10 = a00 * c10 - a10 * c00
    let aic11 = a00 * c11 - a10 * c01
    x[0,0] := itr * (c00 + idet * (aic00*a11-aic01*a10))
    x[0,1] := itr * (c01 + idet * (aic01*a00-aic00*a01))
    x[1,0] := itr * (c10 + idet * (aic10*a11-aic11*a10))
    x[1,1] := itr * (c11 + idet * (aic11*a00-aic10*a01))
  elif nc==3:
    var ad {.noInit.}: type(a)
    adjugate(ad, a)
    let t = a[0,0] + a[1,1] + a[2,2]
    let s = ad[0,0] + ad[1,1] + ad[2,2]
    let r = a[0,0]*ad[0,0] + a[0,1]*ad[1,0] + a[0,2]*ad[2,0]
    var ac {.noInit.}: type(a)
    var ca {.noInit.}: type(a)
    var aca {.noInit.}: type(a)
    var adc {.noInit.}: type(a)
    var cad {.noInit.}: type(a)
    var adcad {.noInit.}: type(a)
    mul(ac, a, c)
    mul(ca, c, a)
    mul(aca, ac, a)
    mul(adc, ad, c)
    mul(cad, c, ad)
    mul(adcad, adc, ad)
    let c2 = 1/(2*(s*t-r))
    let c0 = c2*(s+t*t)
    let c1 = c2*(t/r)
    let c4 = c2*(t)
    for i in 0..2:
      for j in 0..2:
        x[i,j] := c0*c[i,j] + c1*adcad[i,j] + c2*(aca[i,j]-adc[i,j]-cad[i,j]) -
                  c4*(ac[i,j]+ca[i,j])
  else:
    sylsolveN(x, a, c)

