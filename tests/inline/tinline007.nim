import base/metaUtils

proc g[T;N:static[int]](x:array[N,T]) =
  var s = ""
  for i in 0..<N:
    if i>0: s &= " , "
    s &= $x[i]
  echo "x = [ ",s," ] has size ",N*sizeof(T)
echo "* Types with generic parameters"
proc gt[T] =
  type
    M[N:static[int]] = object
      d:array[N,T]
  var A:M[3]
  proc g[N:static[int]](x:M[N]) = x.d.g
  proc `[]`[N:static[int]](x:M[N],i:int):T = x.d[i]
  proc `[]=`[N:static[int]](x:var M[N],i:int,y:T) = x.d[i] = y
  inlineProcs:
    for i in 0..<A.N:
      A[i] = T(i)
    g(A)
gt[float]()
gt[int]()                     # Note github issue #6126
