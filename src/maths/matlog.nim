import base
import matrixConcept
import types
import matexp

proc log1p*(a: Mat1): auto {.noInit.} =
  mixin simdMax
  var x: MatrixArray[a.nrows,a.ncols,type(a[0,0])]
  var p: ExpParam
  p.scale = 20
  p.kind = ekPoly
  p.order = 4
  x := 0
  let an = (1+a).norm2
  let estop = (10*a.nrows*epsilon(a.numberType))^2
  let maxit = 20
  var it = 0
  while true:
    let ex = p.exp(-0.5*x)
    let f = ex*(1+a)*ex
    let fn = f.norm2
    let ft = f.trace.adj
    var w = ft/fn
    #if (1-f).norm2 < 1e-6: w := 1
    let lw = ln(w)
    x -= lw
    let wf = w*f
    let r = 1 - wf
    let rn = r.norm2.simdMax
    echo it, " ", rn, " ", (x*a-a*x).norm2.simdMax, "   ", sqrt(an*(ex*ex).norm2).simdMax
    if rn < estop or it > maxit: break
    if rn < 1e-66:
      x -= r
    else:
      #let t = r + 0.5*r*r
      #let t = r + 0.5*r*(1+r)*r
      #let t = r + 0.5*r*(1+r+r*r)*r
      #let tt = 1+r+r*r
      var tt = r*r
      tt += r
      tt += 1
      let t = r + 0.5*r*tt*r
      #let t = r + 0.5*r*(1+r+r*(1+r)*r)*r
      x -= t
      #x -= r
    inc it
  x

proc logm*(a: Mat1): auto {.noInit.} =
  log1p(a-1)

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

  proc testlogm(T: typedesc) =
    var m: T
    setMat(m)
    #echo m
    echo "N: ", m.nrows
    var p: ExpParam
    p.scale = 20
    p.kind = ekPoly
    p.order = 4
    let l = logm(m)
    let e = p.exp(l)
    let d = (e-m).norm2.simdMax
    let ei = p.exp(-l)
    let r2 = (1-ei*m).norm2.simdMax
    echo d, " ", r2
    let l2 = logm(e)
    echo (l2-l).norm2.simdMax

  type
    Cmplx[T] = ComplexType[T]
    CM[N:static[int],T] = MatrixArray[N,N,Cmplx[T]]

  #testlogm(CM[1,SimdD8])
  #testlogm(CM[2,SimdD8])

  #testlogm(CM[1,float])
  #testlogm(CM[2,float])
  #testlogm(CM[3,float])
  #testlogm(CM[4,float])
  #testlogm(CM[5,float])
  #testlogm(CM[6,float])
  #testlogm(CM[7,float])
  #testlogm(CM[8,float])
  #testlogm(CM[9,float])
  #testlogm(CM[10,float])

  proc testproj(T: typedesc) =
    var m: T
    setMat(m)
    #echo m
    echo "N: ", m.nrows
    var p: ExpParam
    p.scale = 20
    p.kind = ekPoly
    p.order = 4
    let l = logm(m)
    var t = l
    projectTAH(t)
    for scale in 0..20:
      p.scale = scale
      let e = p.exp(t)
      let d = determinant(e)
      let r = e.adj*e - 1
      echo ">> ", scale, " ", (d-1).norm2.simdMax, " ", r.norm2.simdMax
  testproj(CM[2,SimdD8])
  testproj(CM[3,SimdD8])
  testproj(CM[4,SimdD8])
  testproj(CM[5,SimdD8])

  #testproj(CM[2,float])
  #testproj(CM[3,float])
  #testproj(CM[4,float])
  #testproj(CM[5,float])
