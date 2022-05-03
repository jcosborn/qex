import unittest
export unittest
import qex

var AT* = -1.0   ## Absolute Tolerance
var CT* = 1e-13  ## (relative) Comparison Tolerance for float
var CT32* = 1e-5  ## (relative) Comparison Tolerance for float32
template getCT*(x: float): untyped = CT
template getCT*(x: float32): untyped = CT32
template getCT*(x: typedesc[float]): untyped = CT
template getCT*(x: typedesc[float32]): untyped = CT32

# passes if relative or absolute tolerance passes
proc `~`*(x:SomeNumber, y:SomeNumber):bool =
  let ct = max(getCT(x),getCT(y))
  if ct > 0:
    if x==0 or y==0: result = x==y
    else: result = abs(x-y)/max(abs(x),abs(y)) <= ct
  if AT > 0:
    result = result or (abs(x-y)<=AT)
  if not result: echo x," !~ ",y
template `~`*(x,y:int):bool = `~`(float(x),float(y))
template `~`*(x:int,y:float):bool = `~`(float(x),y)
template `~`*(x:float,y:int):bool = `~`(x,float(y))

proc `~`*(x:ComplexProxy,y:ComplexProxy2):bool =
  let
    x2 = x.norm2
    y2 = y.norm2
    d2 = norm2(x-y)
    ct = max(getCT(x.numberType),getCT(y.numberType))
  if ct > 0:
    if x2==0 or y2==0: result = x==y
    else: result = d2/max(x2,y2) <= ct*ct
  if AT > 0:
    result = result or d2<=AT*AT
  if not result: echo x," !~ ",y

proc `~`*(x,y:Field):bool =
  result = true
  var
    x2,y2:float
    xt,yt:type(trace x)
    ct = max(getCT(numberType x),getCT(numberType y))
  threads:
    xt = x.trace
    yt = y.trace
    x2 = x.norm2
    y2 = y.norm2
  block:
    let c = xt~yt
    result = result and c
    if not c:
      echo "trace of the fields differ."
  block:
    let c = x2~y2
    result = result and c
    if not c:
      echo "norm2 of the fields differ."
  var
    t = newOneOf x
    t2:float
    tt:type(trace t)
  threads:
    t := x-y
    tt = t.trace
    t2 = t.norm2
  block:
    if ct>0:
      block:
        let c = abs(t2)/max(abs(x2),abs(y2)) <= ct*ct
        result = result and c
        if not c:
          echo "norm2 of the field difference differs (CT)."
      block:
        let c = norm2(tt)/max(norm2(xt),norm2(yt)) <= ct*ct
        result = result and c
        if not c:
          echo "trace of the field difference differs (CT)."
    if AT>0:
      let c = result or (abs(t2) <= AT*AT and norm2(tt) <= AT*AT)
      result = result and c
      if not c:
        echo "norm2 of the field difference differs (AT)."
  if not result:
    echo "x2: ",x2," y2: ",y2," (x-y)^2 ",t2
    echo "TrX: ",xt," TrY: ",yt," TrX-Y: ",tt

proc `~`*[T](x,y:openarray[T]):bool =
  result = x.len == y.len
  for i in 0..<x.len:
    result = result and (x[i] ~ y[i])
  if not result:
    proc show[T](x:openArray[T]):string =
      result = "[ "
      for c in x:
        result &= $c & " "
      result &= "]"
    echo x.show," !~ ",y.show

template `!~`*(x,y:typed):bool =
  bind `~`
  not(x~y)

template subtest*(t: typed): untyped {.dirty.} =
  let prsave = programResult
  programResult = 0
  t
  if programResult != 0:
    fail()
  else:
    programResult = prsave

when isMainModule:
  suite "Test of testutils":
    test "~ and !~":
      check(0~0)
      check(0!~1e-16)
      check(1~(1+1e-16))
      check(1!~(1+1e-12))
    test "CT = 1e-3":
      CT = 1e-3
      check(0~0)
      check(0!~1e-16)
      check(1~(1+1e-4))
      check(1!~(1+1e-2))
