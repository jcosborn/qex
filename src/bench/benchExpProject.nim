import qex
import maths/groupOps

const sites = 16
type M[R] = MatrixArray[3,3,ComplexType[R]]

qexInit()
letParam:
  scale = 0
  # pullback, logdet, gradient, or combined
  op = "pullback"
  # scalar or simd
  prec = "scalar"
  # Frobenius norm of the projected input
  norm = 0.03
  reps = 2048
  trials = 5
installHelpParam()
processHelpParam()
doAssert scale in [0, expProjectTAHScale]
doAssert op in ["pullback", "logdet", "gradient", "combined"]
doAssert prec in ["scalar", "simd"]
doAssert norm >= 0 and reps > 0 and trials > 0
doAssert nRanks == 1, "run the site benchmark on one rank"

proc fence(a,b,c,d,e: pointer) {.inline.} =
  # Make input and output batches observable across timed repetitions.
  {.emit: "__asm__ __volatile__(\"\" : : \"g\"(`a`), \"g\"(`b`), \"g\"(`c`), \"g\"(`d`), \"g\"(`e`) : \"memory\");".}

proc sample(k: int): tuple[m,c: M[float64]] =
  var v: VectorArray[8,float64]
  for i in 0..<8:
    v[i] = float(i+1)/13.0 + 1e-4*float((k*(i+3)) mod 11 - 5)
  result.m.suFromVec(v)
  result.m := (norm/sqrt(norm2(result.m)))*result.m
  const
    re = [[0.09,0.013,-0.008],[0.013,0.11,0.017],[-0.008,0.017,0.14]]
    im = [[0.0,0.021,0.011],[-0.021,0.0,-0.014],[-0.011,0.014,0.0]]
  for i in 0..<3:
    for j in 0..<3:
      result.m[i,j].re += re[i][j]
      result.m[i,j].im += im[i][j]
      result.c[i,j].re = float((2*i+5*j+k) mod 7 - 3)/11.0
      result.c[i,j].im = float((3*i+2*j+2*k) mod 9 - 4)/13.0
    result.m[i,i].re += 1e-4*float(k*(i+1))

proc run[R](op: static string, sc: static int) =
  const lanes = simdLength(R)
  static: doAssert sites mod lanes == 0
  var m,c,g,p: array[sites div lanes,M[R]]
  var l: array[sites div lanes,R]
  for k in 0..<m.len:
    var s: array[lanes,tuple[m,c:M[float64]]]
    for lane in 0..<lanes: s[lane] = sample(k*lanes+lane)
    for i in 0..<3:
      for j in 0..<3:
        var mr,mi,cr,ci: array[lanes,float64]
        for lane in 0..<lanes:
          mr[lane] = s[lane].m[i,j].re
          mi[lane] = s[lane].m[i,j].im
          cr[lane] = s[lane].c[i,j].re
          ci[lane] = s[lane].c[i,j].im
        when R is float64:
          m[k][i,j] = newComplex(mr[0],mi[0])
          c[k][i,j] = newComplex(cr[0],ci[0])
        else:
          m[k][i,j].re := mr; m[k][i,j].im := mi
          c[k][i,j].re := cr; c[k][i,j].im := ci
  template step =
    for k in 0..<m.len:
      when op == "pullback": p[k].expProjectTAHPullback(m[k],c[k],scale=sc)
      elif op == "logdet": l[k] := expProjMulLogJac(m[k],scale=sc)
      elif op == "gradient": g[k].expProjMulLogJacGrad(m[k],scale=sc)
      else: g[k].expProjMulLogJacGrad(p[k],m[k],c[k],scale=sc)
    fence(addr m[0],addr c[0],addr g[0],addr p[0],addr l[0])
  for _ in 0..<8: step
  for trial in 0..<trials:
    tic("SU3 " & op)
    for _ in 0..<reps: step
    let secs = getElapsedTime()
    toc("sample")
    echo op, " scale=", sc, " prec=", prec, " norm=", norm,
      " trial=", trial, " ns/site=", 1e9*secs/float(reps*sites)
  var checksum = 0.0
  for k in 0..<m.len:
    checksum += simdSum(norm2(g[k])+norm2(p[k])+l[k]*l[k])
  echo "checksum: ", checksum
  echoProf()

template select(R: typedesc, sc: static int) =
  case op
  of "pullback": run[R]("pullback",sc)
  of "logdet": run[R]("logdet",sc)
  of "gradient": run[R]("gradient",sc)
  of "combined": run[R]("combined",sc)
  else: discard

if prec == "scalar":
  if scale == 0: select(float64,0)
  else: select(float64,expProjectTAHScale)
else:
  type Native = evalType(default(DComplexV).re)
  if scale == 0: select(Native,0)
  else: select(Native,expProjectTAHScale)
qexFinalize()
