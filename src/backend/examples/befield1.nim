import qex
import physics/qcdTypes
import backend/[accel,cpugpu,cgfield]
import bench/commonBench
#import strformat
import macros
import base/metaUtils
import parseUtils
import sequtils, strutils

template `*`*[X:SomeNumber,Y:VectorArrayObj](x: typedesc[X], y: typedesc[Y]): typedesc =
  VectorArrayObj[Y.I, X * Y.T]
template `*`*[X:SomeNumber,Y:AsVector](x: typedesc[X], y: typedesc[Y]): typedesc =
  asVector(X * Y[])
template `*`*[X:AsMatrix,Y:AsVector](x: typedesc[X], y: typedesc[Y]): typedesc =
  asVector(X[] * Y[])
template `*`*[N,M:static int;X,Y](x: typedesc[MatrixArrayObj[N,M,X]], y: typedesc[VectorArrayObj[M,Y]]): typedesc =
  VectorArrayObj[N, X * Y]
template `+`*[I:static int,X,Y](x: typedesc[VectorArrayObj[I,X]], y: typedesc[VectorArrayObj[I,Y]]):
  typedesc = VectorArrayObj[I, X + Y]
template `+`*[X:AsVector,Y:AsVector](x: typedesc[X], y: typedesc[Y]): typedesc =
  asVector(X[] + Y[])

template `*`*[X:SomeNumber,Y:Color](x: typedesc[X], y: typedesc[Y]): typedesc =
  asColor(X * Y[])
template `*`*[X:SomeNumber,Y:Color](x: typedesc[X], y: typedesc[ptr Y]): typedesc =
  #static: echo $Y, "  ", $(Y[])
  asColor(X * Y[])
template `+`*[X:Color,Y:Color2](x: typedesc[X], y: typedesc[Y]): typedesc =
  asColor(X[] + Y[])
template `+`*[X:Color,Y:Color2](x: typedesc[X], y: typedesc[ptr Y]): typedesc =
  asColor(X[] + Y[])
template `*`*[X:Color,Y:Color2](x: typedesc[X], y: typedesc[Y]): typedesc =
  asColor(X[] * Y[])
template `*`*[X:Color,Y:Color2](x: typedesc[X], y: typedesc[ptr Y]): typedesc =
  asColor(X[] * Y[])
template `*`*[X:Color,Y:Color2](x: typedesc[ptr X], y: typedesc[ptr Y]): typedesc =
  asColor(X[] * Y[])
proc `:=`*(r: var Color, x: ptr Color2) =
  r := x[]
proc mul*(r: var Color, x: SomeNumber, y: ptr Color2) =
  mul(r[], x, y[][])
proc mul*(r: var Color, x: ptr Color2, y: ptr Color3) =
  mul(r[], x[][], y[][])
proc add*(r: var Color, x: Color, y: ptr Color) =
  add(r[], x[], y[][])

macro exp2string(x:untyped):auto =
  #echo x.treerepr
  #echo x.repr
  #var s = repr fixBracket symToIdent x
  var s = repr symToIdent x
  #echo x.repr
  #echo s
  let n = skipWhitespace(s)
  result = newlit s[n..^1].replace("\n"," ")
  #echo result

template check(r,eqn: untyped) =
  r.cpuReadWrite
  block:
    var s = r.cpu.newOneOf
    threads:
      r := 0
      threadBarrier()
      eqn
      threadBarrier()
      s := r.cpu
      threadBarrier()
      r := 0
    onGpu:
      eqn
    threads:
      s -= r.cpu
      echo s.norm2
  r.cpuUnused

var bs = newSeq[string](0)
var br = newSeq[array[3,string]](0)
var repin = 0
template bench(r,fps,bps,mm,eqn: untyped) =
  check(r, eqn)
  block:
    let N = r.cpu.l.nSitesOuter
    let V = r.cpu.l.V
    let vol = N*V
    var b = newBench()
    case repin
    of 0:
      benchSingle(b):
        let nrep = b.nrep
        onGpu(vol):
          for rep in 0..<nrep:
            eqn
    of 1:
      benchSingle(b):
        let nrep = b.nrep
        for rep in 0..<nrep:
          onGpu(vol):
            eqn
    else:
      r.cpuReadWrite
      benchSingle(b):
        let nrep = b.nrep
        for rep in 0..<nrep:
          onGpu(vol):
            eqn
      r.cpuUnused
    let memMB = vol * mm.float / (1024.0*1024.0)
    let gb = vol * bps.float * b.perNs
    let gf = vol * fps.float * b.perNs
    #echo &"{memMB:8.3f} {gb:8.3f} {gf:8.3f}"
    #echo memMB|(8,-3), " ", gb|(9,-3), " ", gf|(8,-3)
    bs.add exp2string(eqn)
    br.add [memMB|(8,-3), gb|(9,-3), gf|(8,-3)]

proc test(lat:auto, double:static bool=false) =
  let minl = intParam("minl",0)
  if lat[0] < minl: return
  let maxl = intParam("maxl",0)
  if maxl>0:
    if lat[0] > maxl: return
  var lo = newLayout(lat)

  when double:
    type Real = float
    template newCV(): auto = lo.CgColorVectorD()
    template newCM(): auto = lo.CgColorMatrixD()
    #template newHF: untyped = lo.HalfFermion()
    #template newDF: untyped = lo.DiracFermion()
  else:
    type Real = float32
    template newCV(): auto = lo.CgColorVectorS()
    template newCM(): auto = lo.CgColorMatrixS()

  var v1 = newCV()
  var v2 = newCV()
  #var v3 = newCV()
  var m1 = newCM()
  #var m2 = newCM()
  #var m3 = newCM()
  #var m4 = newCM()
  #var m5 = newCM()
  #var h1 = newHF()
  #var d1 = newDF()
  #var d2 = newDF()
  const nc = v1[0].len
  let sf = sizeof(numberType(v1[0]))
  let nc2 = 2*nc
  let vb = nc2*sf
  let mb = nc*vb
  let mvf = (2*nc2-1)*nc2
  echo "Float type: ", $(v1.numberType), "  repin: ", repin
  macro xfer(x: varargs[untyped]): auto =
    result = newNimNode(nnkStmtList)
    for n in x:
      result.add quote do:
        `n`.copyToGpu
        `n`.cpuUnused
  proc reset =
    threads:
      m1 := 2
      v1 := 1
      v2 := 0
      #v3 := 0
    xfer(m1, v1, v2)
  #echo "done setup"
  #var nbench = 0
  template ioGpu:untyped = declared(inOnGpu)

  reset()
  template job1(v2,v1:untyped) =
    when ioGpu:
      #v2 := Real(0.5)*v2 + v1
      for s in gpuRange(v2.n):
        v2.p[s] = Real(0.5)*v2.p[s] + v1.p[s]
      #for s in gpuRange(v2.n):
      #  for ic in 0..<v2.p[0].getNc:
      #    v2.p[s][ic] = Real(0.5)*v2.p[s][ic] + v1.p[s][ic]
    else:
      v2.cpu := Real(0.5)*v2.cpu + v1.cpu
  bench(v2, 2*nc2, 3*vb, 2*vb):
    job1(v2, v1);  ## loop over outer sites

  reset()
  bench(v2, 2*nc2, 3*vb, 2*vb):
    v2 := Real(0.5)*v2 + v1

  reset()
  bench(v2, mvf, mb+2*vb, mb+2*vb):
    v2 := m1 * v1

qexInit()
echo "rank ", myRank, "/", nRanks
threads:
  echo "thread ", threadNum, "/", numThreads

var bss = newSeq[typeof bs](0)
var brs = newSeq[typeof br](0)
template runtests(double:static bool) =
  bs.setLen(0)
  br.setLen(0)
  test([4,4,4,4], double)
  test([8,8,8,8], double)
  test([12,12,12,12], double)
  test([16,16,16,16], double)
  test([24,24,24,24], double)
  test([32,32,32,32], double)
  test([40,40,40,40], double)
  #test([48,48,48,48], double)
  bss.add bs
  brs.add br
template runt(double:static bool) =
  repin = 0
  runtests(double)
  repin = 1
  runtests(double)
  repin = 2
  runtests(double)
runt(false)
let bss32 = bss
let brs32 = brs
bss.setLen(0)
brs.setLen(0)
runt(true)
echoGpuMem()

proc echoResult(p: string, s: seq, r: auto) =
  let sd = s[0].deduplicate
  for t in 0..<sd.len:
    echo "      MB      GB/s     ", p, "  ", sd[t]
    for i in 0..<r[0].len:
      if s[0][i] == sd[t]:
        echo r[0][i][0], " ", r[0][i][1], " ", r[1][i][1], " ", r[2][i][1]
echoResult("float32", bss32, brs32)
echoResult("float64", bss, brs)

echoProf()
qexFinalize()
