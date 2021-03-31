{.pragma: syclh, header:"<CL/sycl.hpp>".}

type
  SyclVec*[T;N:static[int]] {.importcpp:"sycl::vec", syclh.} = object

proc newSyclVec*[T;N:static[int]]: SyclVec[T,N] {.noinit,
  importcpp:"'0()", constructor, syclh.}

proc newSyclVec*[T](a: T): SyclVec[T,1] {.
  importcpp:"'0(#)", constructor, syclh.}
proc newSyclVec*[T](a,b: T): SyclVec[T,2] {.
  importcpp:"'0(#,#)", constructor, syclh.}
proc newSyclVec*[T](a,b,c,d: T): SyclVec[T,4] {.
  importcpp:"'0(#,#,#,#)", constructor, syclh.}
proc newSyclVec*[T](a,b,c,d,e,f,g,h: T): SyclVec[T,8] {.
  importcpp:"'0(#,#,#,#,#,#,#,#)", constructor, syclh.}

template `len`*[T;N:static[int]](x: typedesc[SyclVec[T,N]]): untyped = N
template `len`*[T;N:static[int]](x: SyclVec[T,N]): untyped = N
proc `[]`*[T;N:static[int]](v: SyclVec[T,N], i: SomeInteger): T {.
  importcpp:"#[#]", syclh.}

proc `[]=`*[T;N:static[int]](v: var SyclVec[T,N], i: SomeInteger,
                             x: SomeNumber) {.
  importcpp:"#[#]=#", syclh.}
proc `:=`*[T;N:static[int]](v: var SyclVec[T,N], x: SomeNumber) {.
  importcpp:"# = #", syclh.}

proc `$`*[T;N:static[int]](v: SyclVec[T,N]): string =
  result = "vec[" & $v[0]
  for i in 1..<N:
    result &= "," & $v[i]
  result &= "]"

proc `+=`*[T;N:static[int]](v: var SyclVec[T,N], x: SyclVec[T,N]) {.
  importcpp:"# += #", syclh.}

proc `*`*[T;N:static[int]](x: SomeNumber, y: SyclVec[T,N]): SyclVec[T,N] {.
  importcpp:"# * #", syclh.}

proc `+`*[T;N:static[int]](x,y: SyclVec[T,N]): SyclVec[T,N] {.
  importcpp:"# + #", syclh.}
proc `-`*[T;N:static[int]](x,y: SyclVec[T,N]): SyclVec[T,N] {.
  importcpp:"# - #", syclh.}
proc `*`*[T;N:static[int]](x,y: SyclVec[T,N]): SyclVec[T,N] {.
  importcpp:"# * #", syclh.}

proc rmul*(r: var SyclVec, x,y: SyclVec) {.importcpp:"# = #*#".}
template mul*(x,y: SyclVec, r: var SyclVec) = rmul(r, x, y)

when isMainModule:
  var v1 = newSyclVec(1'f32)
  echo v1
  var v2 = newSyclVec(1'f32,2)
  echo v2
  var v4 = newSyclVec(1'f32,2,3,4)
  var w4 = v4
  echo v4
  v4[3] = 2
  echo v4
  echo v4+w4
  var v8 = newSyclVec(1'f32,2,3,4,5,6,7,8)
  echo v8
  #var v = newSyclVec[float32,8]()
  #echo v
  const n = 100
  var a4: array[n,SyclVec[float32,16]]
  var b4: array[n,SyclVec[float32,16]]
  var c4: array[n,SyclVec[float32,16]]
  for i in 0..<n:
    a4[i] := i+1
    b4[i] := i+2
  #echo a4
  for i in 0..<n:
    c4[i] = a4[i] + b4[i]
  echo c4

  proc test(x: var SyclVec) =
    x[0] = 2
  echo v4
  test(v4)
  echo v4
