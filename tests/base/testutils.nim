import unittest
export unittest

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
