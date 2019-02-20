import base/metaUtils

echo "* noinit"
proc fr =
  type
    R[K:static[int]] = object
      a:array[K,float]
      s:float
  proc fr[K:static[int]]:R[K] {.noinit.} =
    result.s = 0
    for i in 0..<K:
      result.a[i] = i.float
      result.s += result.a[i]
  var v = fr[5]()
  for x in v.a: echo x
  echo v.s
inlineProcs:
  fr()
