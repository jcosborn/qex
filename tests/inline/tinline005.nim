import base/metaUtils

type T {.used.} = array[3,float]
echo "* Generic parameters"
proc g[T;N:static[int]](x:array[N,T]) =
  var s = ""
  for i in 0..<N:
    if i>0: s &= " , "
    s &= $x[i]
  echo "x = [ ",s," ] has size ",N*sizeof(T)
block:
  inlineProcs:
    var v = [0,1,2,3]
    g v
