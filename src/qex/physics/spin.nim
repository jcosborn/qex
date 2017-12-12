import qex/base/wrapperTypes
import qex/maths/types
import qex/maths

type
  Spin*[T] = object
    v*: T
  Spin2*[T] = Spin[T]
  Spin3*[T] = Spin[T]
  SpinMatrix[I,J:static[int],T] = Spin[MatrixArray[I,J,T]]

template asSpin*(xx: typed): untyped =
  staticTraceBegin: asSpin
  let x_asSpin = xx
  staticTraceEnd: asSpin
  Spin[type(x_asSpin)](v: x_asSpin)

template isWrapper*(x: Spin): untyped = true
template asWrapper*(x: Spin, y: typed): untyped =
  asSpin(y)
template asVarWrapper*(x: Spin, y: typed): untyped =
  asVar(asSpin(y))
makeDeref(Spin, x.T)
template `[]`*(x: Spin, i: any): untyped = x[][i]
template `[]`*(x: Spin, i: any, j: any): untyped = x[][i,j]
template `[]`*(x: Spin, i: any, j: any, y: any): untyped =
  x[][i,j] = y
forwardFunc(Spin, len)
forwardFunc(Spin, nrows)
forwardFunc(Spin, ncols)
forwardFunc(Spin, numberType)
forwardFunc(Spin, nVectors)
forwardFunc(Spin, simdType)
forwardFunc(Spin, simdLength)
template row*(x: Spin, i: any): untyped =
  mixin row
  asSpin(row(x[],i))
template setRow*(r: Spin; x: Spin2; i: int): untyped =
  setRow(r[], x[], i)

template binDDRet(fn,wr,T1,T2) =
  template fn*(x: T1, y: T2): untyped =
    #echoUntyped: fn
    wr(fn(x[], y[]))

binDDRet(`+`, asSpin, Spin, Spin2)
binDDRet(`-`, asSpin, Spin, Spin2)
binDDRet(`*`, asSpin, Spin, Spin2)
binDDRet(`/`, asSpin, Spin, Spin2)

template numberType*[T](x: typedesc[Spin[T]]): untyped = numberType(type(T))
#template numNuumbers*[T](x: typedesc[Spin[T]]): untyped = numberType(T)
template numNumbers*(x: Spin): untyped = numNumbers(x[])
template load1*(x: Spin): untyped = asSpin(load1(x[]))
template assign*(r: var Spin, x: SomeNumber) =
  assign(r[], x)
template assign*(r: var Spin, x: AsComplex) =
  assign(r[], x)
template assign*(r: var Spin, x: Spin2) =
  assign(r[], x[])
template `:=`*(r: var Spin, x: SomeNumber) =
  `:=`(r[], x)
template `:=`*(r: var Spin, x: Spin2) =
  r[] := x[]
template `+=`*(r: var Spin, x: Spin2) =
  staticTraceBegin: peqSpinSpin
  r[] += x[]
  staticTraceEnd: peqSpinSpin
template `-=`*(r: var Spin, x: Spin2) =
  r[] -= x[]
template `*=`*(r: var Spin, x: SomeNumber) =
  `*=`(r[], x)
template iadd*(r: var Spin, x: AsComplex) =
  iadd(r[], x)
template iadd*(r: var Spin, x: Spin2) =
  iadd(r[], x[])
template isub*(r: var Spin, x: Spin2) =
  isub(r[], x[])
template imul*(r: var Spin, x: SomeNumber) =
  imul(r[], x)
template imadd*(r: var Spin, x: Spin2, y: Spin3) =
  imadd(r[], x[], y[])
template imsub*(r: var Spin, x: Spin2, y: Spin3) =
  imsub(r[], x[], y[])
template `*`*(x: Spin, y: SomeNumber): untyped =
  asSpin(x[] * y)
template `*`*(x: SomeNumber, y: Spin2): untyped =
  asSpin(x * y[])
template `*`*(x: AsComplex, y: Spin2): untyped =
  asSpin(x * y[])
template sub*(r: var Spin, x: Spin2, y: Spin3) =
  sub(r[], x[], y[])
template mul*(r: var Spin, x: Spin2, y: Spin3) =
  mul(r[], x[], y[])
template random*(x: var Spin) =
  gaussian(x[], r)
template gaussian*(x: var Spin, r: var untyped) =
  gaussian(x[], r)
template projectU*(r: var Spin, x: Spin2) =
  projectU(r[], x[])
template norm2*(x: Spin): untyped = norm2(x[])
template inorm2*(r: var any, x: Spin2) = inorm2(r, x[])
template dot*(x: Spin, y: Spin2): untyped =
  dot(x[], y[])
template idot*(r: var any, x: Spin2, y: Spin3) = idot(r, x[], y[])
template redot*(x: Spin, y: Spin2): untyped =
  redot(x[], y[])
template trace*(x: Spin): untyped = trace(x[])

template spinVector*[T](x:static[int], a: untyped): untyped =
  #const
  #  I:int = x
  #type
  #  E = T
  #  VA = VectorArray[I,E]
  #  VAO = VectorArrayObj[I,E]
  #Spin[VA](v: VA(v: a))
  #Spin[VA](v: VA(v: VAO(vec: a)))
  #asSpin(VA(v: a))
  #static: echo "spinVector"
  asSpin(asVectorArray(a))
  #let t1 = asVectorArray(a)
  #static: echo "spinVector1"
  #let t = asSpin(t1)
  #static: echo "spinVector2"
  #t
template spinMatrix*[T](x,y:static[int], a: untyped): untyped =
  const
    I:int = x
    J:int = y
  type
    E = T
    MA = MatrixArray[I,J,E]
    MAO = MatrixArrayObj[I,J,E]
  Spin[MA](v: MA(v: MAO(mat: a)))
#template spinMatrix*[I,J:static[int],T](a: untyped): untyped =
#  Spin[MatrixArray[I,J,T]](v: MatrixArray[I,J,T](v: MatrixArrayObj[I,J,T](mat: a)))
#template spinMatrix*(I,J,T,a: untyped): untyped =
#  Spin[MatrixArray[I,J,T]](v: MatrixArray[I,J,T](v: MatrixArrayObj[I,J,T](mat: a)))

const z0 = ComplexType[float](v: ComplexObj[float](re:  0.0, im:  0.0))
const z1 = ComplexType[float](v: ComplexObj[float](re:  1.0, im:  0.0))
const zi = ComplexType[float](v: ComplexObj[float](re:  0.0, im:  1.0))
const n1 = ComplexType[float](v: ComplexObj[float](re: -1.0, im:  0.0))
const ni = ComplexType[float](v: ComplexObj[float](re:  0.0, im: -1.0))

template s(r,c,x: untyped): untyped =
  spinMatrix[ComplexType[float]](r,c,x)
template g(x: untyped): untyped = s(4,4,x)
template p(x: untyped): untyped = s(2,4,x)
template r(x: untyped): untyped = s(4,2,x)

const
  gamma0* = g([[ z1, z0, z0, z0 ],
               [ z0, z1, z0, z0 ],
               [ z0, z0, z1, z0 ],
               [ z0, z0, z0, z1 ]])
  gamma1* = g([[ z0, z0, z0, zi ],
               [ z0, z0, zi, z0 ],
               [ z0, ni, z0, z0 ],
               [ ni, z0, z0, z0 ]])
  gamma2* = g([[ z0, z0, z0, n1 ],
               [ z0, z0, z1, z0 ],
               [ z0, z1, z0, z0 ],
               [ n1, z0, z0, z0 ]])
  gamma3* = g([[ z0, z0, zi, z0 ],
               [ z0, z0, z0, ni ],
               [ ni, z0, z0, z0 ],
               [ z0, zi, z0, z0 ]])
  gamma4* = g([[ z0, z0, z1, z0 ],
               [ z0, z0, z0, z1 ],
               [ z1, z0, z0, z0 ],
               [ z0, z1, z0, z0 ]])
  gamma5* = g([[ z1, z0, z0, z0 ],
               [ z0, z1, z0, z0 ],
               [ z0, z0, n1, z0 ],
               [ z0, z0, z0, n1 ]])

  spprojmat1p* = p([[ z1, z0, z0, zi ],
                    [ z0, z1, zi, z0 ]])
  spprojmat1m* = p([[ z1, z0, z0, ni ],
                    [ z0, z1, ni, z0 ]])
  spprojmat2p* = p([[ z1, z0, z0, n1 ],
                    [ z0, z1, z1, z0 ]])
  spprojmat2m* = p([[ z1, z0, z0, z1 ],
                    [ z0, z1, n1, z0 ]])
  spprojmat3p* = p([[ z1, z0, zi, z0 ],
                    [ z0, z1, z0, ni ]])
  spprojmat3m* = p([[ z1, z0, ni, z0 ],
                    [ z0, z1, z0, zi ]])
  spprojmat4p* = p([[ z1, z0, z1, z0 ],
                    [ z0, z1, z0, z1 ]])
  spprojmat4m* = p([[ z1, z0, n1, z0 ],
                    [ z0, z1, z0, n1 ]])

  spreconmat1p* = r([[ z1, z0 ],
                     [ z0, z1 ],
                     [ z0, ni ],
                     [ ni, z0 ]])
  spreconmat1m* = r([[ z1, z0 ],
                     [ z0, z1 ],
                     [ z0, zi ],
                     [ zi, z0 ]])
  spreconmat2p* = r([[ z1, z0 ],
                     [ z0, z1 ],
                     [ z0, z1 ],
                     [ n1, z0 ]])
  spreconmat2m* = r([[ z1, z0 ],
                     [ z0, z1 ],
                     [ z0, n1 ],
                     [ z1, z0 ]])
  spreconmat3p* = r([[ z1, z0 ],
                     [ z0, z1 ],
                     [ ni, z0 ],
                     [ z0, zi ]])
  spreconmat3m* = r([[ z1, z0 ],
                     [ z0, z1 ],
                     [ zi, z0 ],
                     [ z0, ni ]])
  spreconmat4p* = r([[ z1, z0 ],
                     [ z0, z1 ],
                     [ z1, z0 ],
                     [ z0, z1 ]])
  spreconmat4m* = r([[ z1, z0 ],
                     [ z0, z1 ],
                     [ n1, z0 ],
                     [ z0, n1 ]])

template I(x: typed): untyped =
  imagType(1)*x

proc spproj1p*(r: var any, x: any) =
  ## r: HalfFermion
  ## x: DiracFermion
  #let nc = x[0].len
  #for i in 0..<nc:
  #  r[0][i] = x[0][i] + x[2][i]
  #  r[1][i] = x[1][i] + x[3][i]
  r := spprojmat1p * x

#[
template spproj1p*(x: any): untyped = spprojmat1p * x
template spproj2p*(x: any): untyped = spprojmat2p * x
template spproj3p*(x: any): untyped = spprojmat3p * x
template spproj4p*(x: any): untyped = spprojmat4p * x
template spproj1m*(x: any): untyped = spprojmat1m * x
template spproj2m*(x: any): untyped = spprojmat2m * x
template spproj3m*(x: any): untyped = spprojmat3m * x
template spproj4m*(x: any): untyped = spprojmat4m * x
]#

template spproj1p*(xx: typed): untyped =
  let x = xx
  let v0 = x[][0] + I(x[][3])
  let v1 = x[][1] + I(x[][2])
  spinVector[type(v0)](2,[v0,v1])
template spproj2p*(xx: typed): untyped =
  let x = xx
  spinVector[type(x[][0])](2,[x[][0]-x[][3],x[][1]+x[][2]])
template spproj3p*(xx: typed): untyped =
  let x = xx
  spinVector[type(x[][0])](2,[x[][0]+I(x[][2]),x[][1]-I(x[][3])])
template spproj4p*(xx: typed): untyped =
  let x = xx
  spinVector[type(x[][0])](2,[x[][0]+x[][2],x[][1]+x[][3]])
template spproj1m*(xx: typed): untyped =
  let x = xx
  let v0 = x[][0] - I(x[][3])
  let v1 = x[][1] - I(x[][2])
  spinVector[type(v0)](2,[v0,v1])
template spproj2m*(xx: typed): untyped =
  let x = xx
  spinVector[type(x[][0])](2,[x[][0]+x[][3],x[][1]-x[][2]])
template spproj3m*(xx: typed): untyped =
  let x = xx
  spinVector[type(x[][0])](2,[x[][0]-I(x[][2]),x[][1]+I(x[][3])])
template spproj4m*(xx: typed): untyped =
  let x = xx
  spinVector[type(x[][0])](2,[x[][0]-x[][2],x[][1]-x[][3]])

proc sprecon1p*(r: var any, x: any) =
  ## r: DiracFermion
  ## x: HalfFermion
  let nc = x[0].len
  for i in 0..<nc:
    r[0][i] = x[0][i]
    r[1][i] = x[1][i]
    r[2][i] = x[0][i]
    r[3][i] = x[1][i]

#[
template sprecon1p*(x: typed): untyped = spreconmat1p * x
template sprecon2p*(x: typed): untyped = spreconmat2p * x
template sprecon3p*(x: typed): untyped = spreconmat3p * x
template sprecon4p*(x: typed): untyped = spreconmat4p * x
template sprecon1m*(x: typed): untyped = spreconmat1m * x
template sprecon2m*(x: typed): untyped = spreconmat2m * x
template sprecon3m*(x: typed): untyped = spreconmat3m * x
template sprecon4m*(x: typed): untyped = spreconmat4m * x
]#

template sprecon1p*(xx: typed): untyped =
  let x = xx
  spinVector[type(x[][0])](4,[x[][0],x[][1],-I(x[][1]),-I(x[][0])])
template sprecon2p*(xx: typed): untyped =
  let x = xx
  spinVector[type(x[][0])](4,[x[][0],x[][1],x[][1],-x[][0]])
template sprecon3p*(xx: typed): untyped =
  let x = xx
  spinVector[type(x[][0])](4,[x[][0],x[][1],-I(x[][0]),I(x[][1])])
template sprecon4p*(xx: typed): untyped =
  let x = xx
  spinVector[type(x[][0])](4,[x[][0],x[][1],x[][0],x[][1]])
template sprecon1m*(xx: typed): untyped =
  let x = xx
  spinVector[type(x[][0])](4,[x[][0],x[][1],I(x[][1]),I(x[][0])])
template sprecon2m*(xx: typed): untyped =
  let x = xx
  spinVector[type(x[][0])](4,[x[][0],x[][1],-x[][1],x[][0]])
template sprecon3m*(xx: typed): untyped =
  let x = xx
  spinVector[type(x[][0])](4,[x[][0],x[][1],I(x[][0]),-I(x[][1])])
template sprecon4m*(xx: typed): untyped =
  let x = xx
  spinVector[type(x[][0])](4,[x[][0],x[][1],-x[][0],-x[][1]])

when isMainModule:
  echo gamma0[0,0]
  let g2 = gamma0 + gamma4
  #echo g2[0,0]
