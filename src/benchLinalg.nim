import qex
import qcdTypes
import times
import macros
import stdUtils
import parseUtils
#import optimize
import metaUtils

proc symToIdent(x: NimNode): NimNode =
  case x.kind:
    of nnkCharLit..nnkUInt64Lit:
      result = newNimNode(x.kind)
      result.intVal = x.intVal
    of nnkFloatLit..nnkFloat64Lit:
      result = newNimNode(x.kind)
      result.floatVal = x.floatVal
    of nnkStrLit..nnkTripleStrLit:
      result = newNimNode(x.kind)
      result.strVal = x.strVal
    of nnkIdent, nnkSym:
      result = newIdentNode($x)
    of nnkOpenSymChoice:
      result = newIdentNode($x[0])
    else:
      result = newNimNode(x.kind)
      for c in x:
        result.add symToIdent(c)

macro exp2string(x:untyped):auto =
  var s = repr symToIdent x
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
  let nrep = int(1e11/flops)
  var t0 = epochTime()
  threads:
    for rep in 1..nrep:
      eqn
  var t1 = epochTime()
  let dt = t1 - t0
  let mf = (nrep.float*flops)/(1e6*dt)
  let mb = (nrep.float*bytes)/(1e6*dt)
  echo "(", exp2string(eqn), ") secs: ", dt|(5,3), "  mf: ", mf.int,
       "  mb: ", mb.int

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

  bench(8*nc*nc*nc, sf*2*4*nc*nc):
    m3 += m1 * m2

  bench(8*nc*nc*nc, sf*2*4*nc*nc):
    for e in v2:
      imaddMMM(m3[e], m1[e], m2[e])

  bench((8*nc-2)*nc*nc, sf*2*3*nc*nc):
    m3 := m1.adj * m2

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
