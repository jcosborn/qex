import unittest
export unittest
import qex

var AT* = -1.0   ## Absolute Tolerance
var CT* = 1e-13  ## (relative) Comparison Tolerance

# passes if relative or absolute tolerance passes
proc `~`*(x,y:float):bool =
  if CT > 0:
    if x==0 or y==0: result = x==y
    else: result = abs(x-y)/max(abs(x),abs(y)) <= CT
  if AT > 0:
    result = result or (abs(x-y)<=AT)
  if not result: echo x," !~ ",y

proc `~`*(x:ComplexProxy,y:ComplexProxy2):bool =
  let
    x2 = x.norm2
    y2 = y.norm2
    d2 = norm2(x-y)
  if CT > 0:
    if x2==0 or y2==0: result = x==y
    else: result = d2/max(x2,y2) <= CT*CT
  if AT > 0:
    result = result or d2<=AT*AT
  if not result: echo x," !~ ",y

proc `~`*(x,y:Field):bool =
  var
    x2,y2:float
    xt,yt:type(trace x)
  threads:
    xt = x.trace
    yt = y.trace
    x2 = x.norm2
    y2 = y.norm2
  result = x2~y2 and xt~yt
  if not result:
    echo "norm2 of the fields differ."
  else:
    var
      t = newOneOf x
      t2:float
      tt:type(trace t)
    threads:
      t := x-y
      tt = t.trace
      t2 = t.norm2
    if CT>0:
      result = abs(t2)/max(abs(x2),abs(y2)) <= CT*CT and norm2(tt)/max(norm2(xt),norm2(yt)) <= CT*CT
    if AT>0:
      result = result or (abs(t2) <= AT*AT and norm2(tt) <= AT*AT)
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
