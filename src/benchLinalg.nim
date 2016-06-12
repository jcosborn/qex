import qex
import qcdTypes
import times
import macros
import stdUtils
import parseUtils

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

template bench(fps:SomeNumber; eqn:untyped) =
  let vol = lo.nSites.float
  let flops = vol * fps.float
  let nrep = int(1e10/flops)
  let str = toString(eqn)
  var t0,t1:float
  threads:
    if threadNum==0: t0 = epochTime()
    for rep in 1..nrep:
      eqn
    if threadNum==0: t1 = epochTime()
  let dt = t1 - t0
  let mf = (nrep.float*flops)/(1e6*dt)
  echo "(", str, ") secs: ", dt|(5,3), "  mf: ", mf.int

proc test(lat:any) =
  var lo = newLayout(lat)

  var v1 = lo.ColorVector()
  var v2 = lo.ColorVector()
  var m1 = lo.ColorMatrix()
  var v3 = lo.ColorVector()
  threads:
    #m1 := 2
    v1 := 1

  bench((8*nc-2)*nc):
    v2 := m1 * v1

  bench((8*nc-2)*nc):
    for i in v2:
      v2[i] := m1[i] * v1[i]

  bench((8*nc-2)*nc):
    for i in v2:
      mul(v2[i], m1[i], v1[i])

  bench((8*nc-2)*nc):
    for i in v2:
      var vt{.noInit.}:type(perm(v1[i],0))
      assign(vt, v1[i])
      v2[i] := m1[i] * vt

  bench((8*nc-2)*nc):
    for i in v2:
      let vt = perm(v1[i], 0)
      v2[i] := m1[i] * vt

  bench((8*nc-2)*nc):
    for i in v2:
      let vt = perm(v1[i], 1)
      v2[i] := m1[i] * vt

  bench((8*nc-2)*nc):
    for i in v2:
      let vt = perm(v1[i], 1)
      mul(v2[i], m1[i], vt)

  bench((8*nc)*nc):
    v2 += m1 * v1

  bench((8*nc)*nc):
    for i in v2:
      v2[i] += m1[i] * v1[i]

  bench((8*nc)*nc):
    for i in v2:
      imadd(v2[i], m1[i], v1[i])

  bench((8*nc)*nc):
    for i in v2:
      let vt = perm(v1[i], 0)
      v2[i] += m1[i] * vt

  bench((8*nc)*nc):
    for i in v2:
      var vt{.noInit.}:type(load(v1[i]))
      assign(vt, v1[i])
      imadd(v2[i], m1[i], vt)

  bench((8*nc)*nc):
    for i in v2:
      let vt = load(v1[i])
      imadd(v2[i], m1[i], vt)

  bench((8*nc)*nc):
    for i in v2:
      let vt = perm(v1[i], 0)
      imadd(v2[i], m1[i], vt)

  bench((8*nc)*nc):
    for i in v2:
      let vt = perm(v1[i], 1)
      imadd(v2[i], m1[i], vt)

  bench((8*nc)*nc):
    let sch = 0.5
    for i in v2:
      var et{.noInit.}:type(sch*v1[i])
      et := sch*v1[i]
      var vt{.noInit.}:type(perm(et, 1))
      perm(vt, 1, et)
      imadd(v2[i], m1[i], vt)

  bench((8*nc)*nc):
    let sch = 0.5
    for i in v2:
      var et{.noInit.}:type(m1[i].adj*v1[i])
      et := m1[i].adj*v1[i]
      var vt{.noInit.}:type(perm(et, 1))
      perm(vt, 1, et)
      imsub(v2[i], sch, vt)

  bench((8*nc)*nc):
    let sch = 0.5
    for i in v2:
      imadd(v2[i], m1[i], sch*v1[i])

  bench(4*(8*nc)*nc):
    let sch = 0.5
    let prms = [0,1,2,0]
    for i in v2:
      for mu in 0..3:
        let et = sch*v1[i]
        let vt = perm(et, prms[mu])
        imadd(v2[i], m1[i], vt)

  bench(4*(8*nc)*nc):
    let sch = 0.5
    let prms = [1,1,1,1]
    for i in v2:
      for mu in 0..3:
        let et = sch*v1[i]
        let vt = perm(et, prms[mu])
        imadd(v2[i], m1[i], vt)

  bench(4*(8*nc)*nc):
    let sch = 0.5
    for i in v2:
      var vr = load(v2[i])
      for mu in 0..3:
        let et = sch*v1[i]
        let vt = perm(et, 1)
        imadd(vr, m1[i], vt)
      assign(v2[i], vr)

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
#checkMem()
qexFinalize()
