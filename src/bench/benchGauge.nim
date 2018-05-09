import qex
import physics/qcdTypes
import gauge
import parseUtils
import times
import macros

proc checkMem =
  echo("mem: (used+free)/total: (", getOccupiedMem(), "+", getFreeMem(), ")/",
       getTotalMem())
  echo GC_getStatistics()
  GC_fullCollect()
  echo("mem: (used+free)/total: (", getOccupiedMem(), "+", getFreeMem(), ")/",
       getTotalMem())
  echo GC_getStatistics()
macro exp2string(x:untyped):auto =
  var s = repr symToIdent x
  let n = skipWhitespace(s)
  newlit s[n..^1]
template bench(fps,bps:SomeNumber; eqn:untyped) =
  let vol = lo.nSites.float
  let flops = vol * fps.float
  let bytes = vol * bps.float
  let nrep = 1 + int(1e8/flops)
  #echo nrep
  var t0 = epochTime()
  #threads:
  for rep in 1..nrep:
    #echo rep
    eqn
  var t1 = epochTime()
  let dt = t1 - t0
  let mf = (nrep.float*flops)/(1e6*dt)
  let mb = (nrep.float*bytes)/(1e6*dt)
  echo "(", exp2string(eqn), ") secs: ", dt|(5,3), "  mf: ", mf.int,
       "  mb: ", mb.int

proc test(lat0: any) =
  var scale = 1
  var lat = lat0*scale
  var lo = newLayout(lat)
  let nd = lo.nDim
  let np = (nd*(nd-1)) div 2
  var g = newSeq[type(lo.ColorMatrix())](nd)
  let nc = g[0][0].ncols
  for i in 0..<nd:
    g[i] = lo.ColorMatrix()
  g.random

  var pl2 = plaq2(g)
  echo "plaq2: ", pl2

  var pl = plaq(g)
  echo pl
  echo "plaq: ", pl.sum

  var ga = gaugeAction1(g)
  echo ga

  resetTimers()

  bench(np*(2*8*nc*nc*nc-1), nd*2*nc*nc*sizeof(numberType(g[0][0]))):
    var pl2 = plaq2(g)

  bench(np*(2*8*nc*nc*nc-1), nd*2*nc*nc*sizeof(numberType(g[0][0]))):
    var pl = plaq(g)

  bench(np*(2*8*nc*nc*nc-1), nd*2*nc*nc*sizeof(numberType(g[0][0]))):
    var ga = gaugeAction1(g)

  echoTimers()
  resetTimers()

qexInit()
checkMem()
test([4,4,4,4])
checkMem()
test([8,8,8,8])
checkMem()
test([12,12,12,12])
checkMem()
test([16,16,16,16])
checkMem()
test([24,24,24,24])
checkMem()
qexFinalize()
