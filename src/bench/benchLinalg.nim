import ../base/globals
#setForceInline(false)
setForceInline(true)
setStaticUnroll(false)
#setStaticUnroll(true)
setNoAlias(false)
#setNoAlias(true)

import qex
import physics/qcdTypes
import times
import macros
#import stdUtils
import parseUtils
#import metaUtils
#import profile

proc fixBracket*(x: NimNode): NimNode =
  if x.kind == nnkCall and eqIdent(x[0],"[]"):
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

#[
proc checkMem =
  echo("mem: (used+free)/total: (", getOccupiedMem(), "+", getFreeMem(), ")/",
       getTotalMem())
  echo GC_getStatistics()
  GC_fullCollect()
  echo("mem: (used+free)/total: (", getOccupiedMem(), "+", getFreeMem(), ")/",
       getTotalMem())
  echo GC_getStatistics()
]#

# fps: flops per site
# bps: bytes moved (load+store) per site
# mm: memory footprint (bytes) per site
var minNrep = int.high
template bench(fps,bps,mm,eqn: untyped) {.dirty.} =
  block:
    let vol = lo.nSites.float
    let flops = vol * fps.float
    let bytes = vol * bps.float
    let vmm = vol * mm.float
    var nrep = 1
    var dt = 0.0
    while true:
      threads:
        threadBarrier()
        let t0 = getTics()
        for rep in 1..nrep:
          eqn
        let t1 = getTics()
        var dtt = ticDiffSecs(t1,t0)
        threadSum(dtt)
        threadMaster: dt = dtt/numThreads.float
      if dt>1: break
      let nnrep = 1 + int(1.1*nrep.float/(dt+1e-9))
      nrep = min(10*nrep, nnrep)
    let mf = int((nrep.float*flops)/(1e6*dt))
    #let mb = int((nrep.float*bytes)/(1e6*dt))
    let mb = int((nrep.float*bytes)/(1024.0*1024.0*dt))
    let mem = vmm/1024.0
    inc nbench
    echo "bench: ",nbench|(-6), "secs: ", dt|(6,3), "  mf: ", mf|7, "  mb: ", mb|7, "  mem: ", mem, "  nrep: ", nrep
    echo exp2string(eqn), "\n"
    minNrep = min(minNrep, nrep)
template bench(fps,bps,eqn: untyped) =
  bench(fps,bps,0,eqn)

proc test(lat:auto) =
  let minl = intParam("minl",0)
  if lat[0] < minl: return
  let maxl = intParam("maxl",0)
  if maxl>0:
    if lat[0] > maxl: return
  else:
    if minNrep == 1: return
  minNrep = int.high
  var lo = newLayout(lat)
  template newCV: untyped = lo.ColorVector()
  template newCM: untyped = lo.ColorMatrix()
  template newHF: untyped = lo.HalfFermion()
  template newDF: untyped = lo.DiracFermion()
  #when true:
  #  template newCV: untyped = lo.ColorVectorS()
  #  template newCM: untyped = lo.ColorMatrixS()
  #  template newHF: untyped = lo.HalfFermionS()
  #  template newDF: untyped = lo.DiracFermionS()
  #else:
  #  template newCV: untyped = lo.ColorVectorD()
  #  template newCM: untyped = lo.ColorMatrixD()
  #  template newHF: untyped = lo.HalfFermionD()
  #  template newDF: untyped = lo.DiracFermionD()

  var v1 = newCV()
  var v2 = newCV()
  var v3 = newCV()
  var m1 = newCM()
  var m2 = newCM()
  var m3 = newCM()
  var m4 = newCM()
  var m5 = newCM()
  var h1 = newHF()
  var d1 = newDF()
  var d2 = newDF()
  const nc = v1[0].len
  let sf = sizeof(numberType(v3[0]))
  let nc2 = 2*nc
  let vb = nc2*sf
  let mb = nc*vb
  let mvf = (2*nc2-1)*nc2
  echo "Float type: ", $(v1.numberType)
  threads:
    m1 := 2
    v1 := 1
    v2 := 0
    v3 := 0
  echo "done setup"
  var nbench = 0

  bench(2*nc2, 3*vb, 2*vb):
    v2 := 0.5*v2 + v1

  bench(mvf, mb+2*vb, mb+2*vb):
    v2 := m1 * v1

  bench(mvf, mb+2*vb, mb+2*vb):
    for e in v2:
      mulVMV(v2[e], m1[e], v1[e])

  bench(mvf, mb+2*vb, mb+2*vb):
    v2 := m1.adj * v1

  bench(mvf+nc2, mb+3*vb, mb+2*vb):
    v2 += m1 * v1

  bench(mvf+nc2, mb+3*vb, mb+2*vb):
    for e in v2:
      imaddVMV(v2[e], m1[e], v1[e])

  bench(mvf+nc2, mb+3*vb, mb+2*vb):
    for e in v2:
      v2[e] += m1[e] * row(m2[e],0)

  bench(mvf, mb+2*vb, mb+2*vb):
    for e in v2:
      var t {.noInit.}: evalType(v1[e])
      perm(t, 1, v1[e])
      v2[e] := m1[e] * t

  bench(nc*mvf, 3*mb, 3*mb):
    m3 := m1 * m2

  bench(nc*mvf, 3*mb, 3*mb):
    m3 := m1.adj * m2

  bench(nc*mvf, 3*mb, 3*mb):
    for e in m3:
      mulMMM(m3[e], m1[e].adj, m2[e])
      #mulMMM(m3[e], m1[e], m2[e])

  bench(nc*mvf, 3*mb, 3*mb):
    for e in m3:
      var vt{.noInit.}:type(v2[e])
      #forStatic i, 0, nc.pred:
      for i in 0..<nc:
        mul(vt, m1[e], row(m2[e],i))
        setRow(m3[e], vt, i)

  bench(nc*(mvf+nc2), 4*mb, 3*mb):
    m3 += m1 * m2

  bench(nc*(mvf+nc2), 4*mb, 3*mb):
    for e in m3:
      imaddMMM(m3[e], m1[e], m2[e])

  bench(nc*(mvf+nc2), 3*mb, 2*mb):
    for e in m3:
      imaddMMM(m3[e], m1[0], m2[e])

  bench(nc*(3*mvf+nc2), 6*mb, 5*mb):
    m1 += (m2*m3) * (m4*m5).adj

  bench(4*(mvf+nc2), mb+4*3*vb, mb+4*2*vb):
    d2 += m1 * d1

  bench(2*nc2, 6*vb, 6*vb):
    for e in h1:
      h1[e] := spproj1p(d1[e])

  bench(2*(nc2+mvf), mb+6*vb, mb+6*vb):
    for e in h1:
      h1[e] := m1[e].adj * spproj1m(d1[e])

  bench(2*(nc2+mvf)+4*nc2, mb+12*vb, mb+8*vb):
    for e in d2:
      d2[e] -= sprecon1p( m1[e] * spproj1p(d1[e]) )

  bench(2*(nc2+mvf)+4*nc2, mb+12*vb, mb+8*vb):
    for e in d2:
      d2[e] -= sprecon1m( m1[e].adj * spproj1m(d1[e]) )

qexInit()
echo "rank ", myRank, "/", nRanks
threads:
  echo "thread ", threadNum, "/", numThreads
freezeTimers()

#checkMem()
test([4,4,4,4])
#test([4,4,4,8])
#test([4,4,4,16])
#test([8,8,8,4])
test([8,8,8,8])
#test([8,8,8,16])
#test([12,12,12,8])
test([12,12,12,12])
#test([12,12,12,24])
#test([16,16,16,8])
test([16,16,16,16])
#test([16,16,16,32])
test([20,20,20,20])
#test([24,24,24,12])
test([24,24,24,24])
#test([24,24,24,48])
#test([32,32,32,16])
test([32,32,32,32])
#test([32,32,32,64])

#checkMem()
#echoTimers()
#echoTimersRaw()
qexFinalize()
