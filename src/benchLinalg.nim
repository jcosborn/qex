import qex
import qcdTypes
import times
import macros
import stdUtils
import parseUtils
#import optimize
import metaUtils
import profile

proc fixBracket*(x: NimNode): NimNode =
  if x.kind == nnkCall and eqIdent(x.name,"[]"):
    #echo x.name.treerepr
    result = newNimNode(nnkBracketExpr)
    for i in 1..<x.len:
      result.add fixBracket(x[i])
  else:
    result = x
    for i in 0..<x.len:
      result[i] = fixBracket(x[i])

macro exp2string(x:untyped):auto =
  var s = repr fixBracket symToIdent x
  let n = skipWhitespace(s)
  newlit s[n..^1]

proc checkMem =
  echo("mem: (used+free)/total: (", getOccupiedMem(), "+", getFreeMem(), ")/",
       getTotalMem())
  echo GC_getStatistics()
  GC_fullCollect()
  echo("mem: (used+free)/total: (", getOccupiedMem(), "+", getFreeMem(), ")/",
       getTotalMem())
  echo GC_getStatistics()

template bench(fps,bps:SomeNumber; eqn:untyped) =
  let vol = lo.nSites.float
  let flops = vol * fps.float
  let bytes = vol * bps.float
  let mm = 0
  var nrep = 1
  var dt = 0.0
  while true:
    let t0 = getTics()
    threads:
      for rep in 1..nrep:
        eqn
    let t1 = getTics()
    dt = ticDiffSecs(t1,t0)
    if dt>1: break
    nrep = 1 + int(1.1*nrep.float/dt)
  let mf = int((nrep.float*flops)/(1e6*dt))
  let mb = int((nrep.float*bytes)/(1e6*dt))
  let mem = vol * (mm.float/(1024.0))
  inc nbench
  echo "bench: ",nbench| -6, "secs: ", dt|(6,3), "  mf: ", mf|7, "  mb: ", mb|7, "  mem: ", mem, "  nrep: ", nrep
  echo exp2string(eqn), "\n"

proc test(lat:any) =
  var lo = newLayout(lat)

  var v1 = lo.ColorVector()
  var v2 = lo.ColorVector()
  var v3 = lo.ColorVector()
  var m1 = lo.ColorMatrix()
  var m2 = lo.ColorMatrix()
  var m3 = lo.ColorMatrix()
  const nc = v1[0].len
  let sf = sizeof(numberType(v3[0]))
  threads:
    m1 := 2
    v1 := 1
    v2 := 0
    v3 := 0
  echo "done setup"
  var nbench = 0

  bench(2*2*nc, sf*2*3*nc):
    v2 := 0.5*v2 + v1

  bench((8*nc-2)*nc, sf*2*(nc+2)*nc):
    v2 := m1 * v1

  bench((8*nc-2)*nc, sf*2*(nc+2)*nc):
    for e in v2:
      mulVMV(v2[e], m1[e], v1[e])

  bench((8*nc)*nc, sf*2*(nc+3)*nc):
    v2 += m1 * v1

  bench((8*nc)*nc, sf*2*(nc+3)*nc):
    for e in v2:
      v2[e] += m1[e] * row(m2[e],0)

  bench((8*nc)*nc, sf*2*(nc+3)*nc):
    for e in v2:
      imaddVMV(v2[e], m1[e], v1[e])

  bench((8*nc-2)*nc, sf*2*(nc+2)*nc):
    v2 := m1.adj * v1

  bench((8*nc-2)*nc*nc, sf*2*3*nc*nc):
    m3 := m1 * m2

  bench((8*nc-2)*nc*nc, sf*2*3*nc*nc):
    m3 := m1.adj * m2

  bench(8*nc*nc*nc, sf*2*4*nc*nc):
    m3 += m1 * m2

  bench(8*nc*nc*nc, sf*2*4*nc*nc):
    for e in v2:
      imaddMMM(m3[e], m1[e], m2[e])

  bench((8*nc-2)*nc*nc, sf*2*3*nc*nc):
    for e in v2:
      mulMMM(m3[e], m1[e].adj, m2[e])

  bench((8*nc-2)*nc*nc, sf*2*3*nc*nc):
    for e in m3:
      var vt{.noInit.}:type(v2[e])
      forStatic i, 0, <nc:
        mul(vt, m1[e], row(m2[e],i))
        setRow(m3[e], vt, i)

  bench(8*nc*nc*nc, sf*2*3*nc*nc):
    for e in v2:
      imaddMMM(m3[e], m1[0], m2[e])


qexInit()
#checkMem()
test([4,4,4,4])
test([8,8,8,8])
test([12,12,12,12])
test([16,16,16,16])
test([24,24,24,24])
test([32,32,32,32])
#checkMem()
qexFinalize()
