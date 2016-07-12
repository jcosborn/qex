import qex
import qcdTypes
import times
import macros
import stdUtils
import parseUtils
import optimize

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

macro toString(x:untyped):auto =
  var s = repr symToIdent x
  let n = skipWhitespace(s)
  newlit s[n..^1]

template bench(fps,bps:SomeNumber; eqn:untyped) =
  let vol = lo.nSites.float
  let flops = vol * fps.float
  let bytes = vol * bps.float
  let nrep = int(2e11/flops)
  var t0 = epochTime()
  threads:
    for rep in 1..nrep:
      eqn
  var t1 = epochTime()
  let dt = t1 - t0
  let mf = (nrep.float*flops)/(1e6*dt)
  let mb = (nrep.float*bytes)/(1e6*dt)
  echo "(", toString(eqn), ") secs: ", dt|(5,3), "  mf: ", mf.int,
       "  mb: ", mb.int

proc test(lat:any) =
  var lo = newLayout(lat)

  var v1 = lo.ColorVector()
  var v2 = lo.ColorVector()
  var m1 = lo.ColorMatrix()
  var v3 = lo.ColorVector()
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
      imaddVMV(v2[e], m1[e], v1[e])

proc checkMem =
  echo("mem: (used+free)/total: (", getOccupiedMem(), "+", getFreeMem(), ")/",
       getTotalMem())
  echo GC_getStatistics()
  GC_fullCollect()
  echo("mem: (used+free)/total: (", getOccupiedMem(), "+", getFreeMem(), ")/",
       getTotalMem())
  echo GC_getStatistics()

qexInit()
#checkMem()
test([4,4,4,4])
test([8,8,8,8])
test([12,12,12,12])
test([16,16,16,16])
test([24,24,24,24])
#checkMem()
qexFinalize()
