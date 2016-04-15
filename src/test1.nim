import qex
import qcdTypes

proc test() =
  var lat = [4,4,4,4]
  var lo = newLayout(lat)

  var v1 = lo.ColorVector()
  var v2 = lo.ColorVector()
  var m1 = lo.ColorMatrix()
  var v3 = lo.ColorVector()
  template T0(x:v3.type):expr =
    shift(v3, 3,1, x)
    v3
  
  threads:
    m1 := 2
    v1 := 1
    for e in v1.all:
      var aa:array[lo.V,float32]
      for i in 0..<lo.V:
        aa[i] = ((myRank*100+e)*10+i).float32
      assign(v1[e][0].re, aa)
    threadSingle:
      if myRank==0:
        echo v1[0][0]
    v2 := m1 * v1
    threadSingle:
      if myRank==0:
        echo v2[0][0]
    shift(v1, dir=3, len=1, v2)
    threadSingle:
      if myRank==0:
        echo v1[0][0]
    v2 := m1 * v1.T0
    threadSingle:
      if myRank==0:
        echo v2[0][0]
    var n1 = threadNum.float
    var n2 = 1.0
    threadSum(n1,n2)
    threadSingle:
      echo n1, " ", n2
      rankSum(n1,n2)
      echo n1, " ", n2

  #  var d = sum(v1.adj()*v2)
  #  var l1 = lsum(v1.adj()*v2)
  #  var l2 = lsum((v1.adj()*v1).r)
  #  var s1:type(l1)
  #  var s2:type(l2)
  #  (s1,s2) = gsum((l1,l2))
  #  threadSingle:
  #    echo("s1: %s  s2: %s" % ($s1,$s2))

qexInit()
test()
qexFinalize()
