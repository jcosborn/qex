proc lexIndex*[T,U](x: openArray[T], size: openArray[U]): int =
  let n = x.len
  assert(size.len == n)
  for i in countdown(n-1,0):
    result = result*size[i] + x[i]

proc lexIndexN[N:static[int],T,U](x: openArray[T], size: openArray[U]): int =
  assert(x.len == N)
  assert(size.len == N)
  template R(i: int) =
    result = result*size[i] + x[i]
    when i>0: R(i-1)
  R(N-1)

proc lexIndexR[T,U](x: openArray[T], size: openArray[U]): int =
  let n = x.len
  assert(size.len == n)
  for i in 0..<n:
    result = result*size[i] + x[i]

proc lexIndexRN[N:static[int],T,U](x: openArray[T], size: openArray[U]): int =
  assert(x.len == N)
  assert(size.len == N)
  #template R(i: int) =
  #  result = result*size[i] + x[i]
  #  when i<N: R(i+1)
  #R(0)
  for i in 0..<N:
    result = result*size[i] + x[i]

proc lexCoord*[T,U](x: var openArray[T], k: SomeInteger, size: openArray[U]) =
  let n = x.len
  var l = k
  assert(size.len == n)
  for i in 0 ..< n:
    x[i] = T l mod size[i]
    l = l div size[i]

#proc divCoords[T,U,V](r: var openArray[T], x: openArray[U], y: openArray[V]) =
#


when isMainModule:
  echo "testing hyper"
  let c = [1,2,3,4]
  let s = [12,12,12,24]
  var r = [0,0,0,0]
  echo lexIndex(c, s)
  echo lexIndexN[4,int,int](c, s)
  echo lexIndexR(c, s)
  echo lexIndexRN[4,int,int](c, s)
