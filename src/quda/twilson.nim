import ../base/globals
setForceInline(false)
#setForceInline(true)
setStaticUnroll(false)
#setStaticUnroll(true)
setNoAlias(false)
#setNoAlias(true)

import qex, physics/qcdTypes
import physics/wilsonD
import physics/wilsonSolve
import qudaWrapper, quda, enum_quda, quda_constants, qudaSet
import strformat
import std/with

proc loadGauge(g: seq) =
  var lo1 = getQudaLayout()
  var qgp = newQudaGaugeParam()
  setLat(lo1.physGeom)
  setWilsonGaugeParam(qgp)

  var g1: array[4,DLatticeColorMatrix]
  for i in 0..<4:
    g1[i].new lo1

  for i in g[0].sites:
    var cv: array[4,cint]
    g[0].l.coord(cv,(g[0].l.myRank,i))
    let ri1 = lo1.rankIndex(cv)
    assert(ri1.rank == r.l.myRank)
    forO mu, 0, 3:
      forO a, 0, 2:
        forO b, 0, 2:
          g1[mu][ri1.index][a,b].re := g[mu]{i}[a,b].re
          g1[mu][ri1.index][a,b].im := g[mu]{i}[a,b].im

  var hGauge: array[4,pointer]
  for i in 0..<4:
    hGauge[i] = cast[pointer](addr g1[i][0])
  loadGaugeQuda(cast[pointer](addr hGauge), addr qgp)

proc convert(r: DLatticeDiracFermionV, lo1: Layout): DLatticeDiracFermion =
  result.new(lo1)
  for i in r.sites:
    var cv: array[4,cint]
    r.l.coord(cv,(r.l.myRank,i))
    let ri1 = lo1.rankIndex(cv)
    assert(ri1.rank == r.l.myRank)
    result[ri1.index] := r{i}

proc convert(r: DLatticeDiracFermionV, x: DLatticeDiracFermion) =
  for i in r.sites:
    var cv: array[4,cint]
    r.l.coord(cv,(r.l.myRank,i))
    let ri1 = x.l.rankIndex(cv)
    assert(ri1.rank == r.l.myRank)
    r{i} := x[ri1.index]

proc checkPlaq(g: seq) =
  var plaq: array[3,cdouble]
  plaqQuda(plaq)
  echo &"QUDA plaq: {plaq[0]:.16g}  spatial: {plaq[1]:.16g}  temporal: {plaq[2]:.16g}"
  let p = plaq(g)
  let ps = p[0] + p[1] + p[2]
  let pt = p[3] + p[4] + p[5]
  echo &"QEX  plaq: {ps+pt:.16g}  spatial: {2*ps:.16g}  temporal: {2*pt:.16g}"

proc getQip(sp: SolverParams, m: float): QudaInvertParam =
  qudaSet.mass = m
  qudaSet.kappa = -1.0
  result = newQudaInvertParam()
  setInvertParam(result)
  result.preserve_source = QUDA_PRESERVE_SOURCE_NO
  result.tune = QUDA_TUNE_YES

proc getQipMg(sp: SolverParams, m: float): QudaInvertParam =
  result = getQip(sp, m)
  var mg_param = newQudaMultigridParam()
  #var mg_inv_param = newQudaInvertParam()
  #var mg_eig_param: array[QUDA_MAX_MG_LEVEL, QudaEigParam]
  #var eig_param = newQudaEigParam()

  setQudaMgSolveTypes()
  setMultigridInvertParam(result)
  mg_param.invert_param = addr result
  for i in 0..<mg_levels:
    #if mg_eig[i]:
    #  mg_eig_param[i] = newQudaEigParam()
    #  setMultigridEigParam(mg_eig_param[i], i)
    #  mg_param.eig_param[i] = addr mg_eig_param[i]
    #else:
      mg_param.eig_param[i] = nil
  setMultigridParam(mg_param)
  var mg_preconditioner = newMultigridQuda(addr mg_param)
  result.preconditioner = mg_preconditioner

proc qudaSolve(w: Wilson; dest,src: Field, m: float, sp: var SolverParams) =
  let lo1 = getQudaLayout()
  var s1 = convert(src, lo1)
  var d1 = convert(dest, lo1)
  var hin = cast[pointer](addr s1[0])
  var hOut = cast[pointer](addr d1[0])
  var qip = getQip(sp, m)
  invertQuda(hout, hin, addr qip)
  convert(dest, d1)
  let s = 1/(4.0+m)
  threads:
    dest *= s

proc qudaSolveMg(w: Wilson; dest,src: Field, m: float, sp: var SolverParams) =
  setQudaDefaultMgTestParams()
  let lo1 = getQudaLayout()
  var s1 = convert(src, lo1)
  var d1 = convert(dest, lo1)
  var hin = cast[pointer](addr s1[0])
  var hOut = cast[pointer](addr d1[0])
  var qip = getQipMg(sp, m)
  invertQuda(hout, hin, addr qip)
  destroyMultigridQuda(qip.preconditioner)
  convert(dest, d1)
  let s = 1/(4.0+m)
  threads:
    dest *= s

proc testW*(w: Wilson; r,t: Field; m: SomeNumber; sp: var SolverParams) =
  let lo1 = r.l.qudaSetup
  loadGauge(w.g)
  checkPlaq(w.g)
  #dslashQuda(hout, hin, addr qip, QUDA_EVEN_PARITY)
  #dslashQuda(hout, hin, addr qip, QUDA_ODD_PARITY)
  #convert(r, r1)
  #threads:
  #  #r := (4.0+m)*(t-2*r)
  #  r := -0.5*r
  qudaSolve(w, r, t, m, sp)

when isMainModule:
  import gauge
  qexInit()
  var defaultLat = @[8,8,8,8]
  defaultSetup()
  echo "rank ", myRank, "/", nRanks
  threads:
    echo "thread ", threadNum, "/", numThreads
  var
    #lo = lat.newLayout
    src = lo.DiracFermion()
    dest = lo.DiracFermion()
    destG = lo.DiracFermion()
    r = lo.DiracFermion()
    #g = lo.newGauge
    rng = RngMilc6.newRNGField lo
  if fn == "":
    g.random rng
  threads:
    g.setBC
    #g.stagPhase
    dest := 0
    destG := 0
    src.z4 rng
    #threadBarrier()
    #src.even = 0
    #threadBarrier()
    #src.odd = 0

  var w = g.newWilson(src)
  let m = floatParam("mass", 0.01)
  let res = floatParam("cg_prec", 1e-9)
  let sloppy = intParam("sloppy", 2).SloppyType
  echo "mass: ", m
  #w.solve(dest, src, m, res, cpuonly = true)
  #w.solve(dest, src, m, res)
  var n2: tuple[a,e,o:float]
  threads:
    echo "src.norm2: ", src.norm2
    echo "src.even: ", src.even.norm2
    echo "src.odd: ", src.odd.norm2
    #n2.a = dest.norm2
    #n2.e = dest.even.norm2
    #n2.o = dest.odd.norm2
    #echo "dest.norm2: ", n2.a
    #echo "dest.even: ", n2.e
    #echo "dest.odd: ", n2.o
    #w.D(r, src, m)
    #w.Ddag(r, src, m)
    #threadBarrier()
    #r := src - r
    #threadBarrier()
    #echo "r.norm2: ", r.norm2
    #echo "r.even: ", r.even.norm2
    #echo "r.odd: ", r.odd.norm2

  var sp = initSolverParams()
  testW(w, dest,src, m, sp)
  n2.a = dest.norm2
  n2.e = dest.even.norm2
  n2.o = dest.odd.norm2
  echo "dest.norm2: ", n2.a
  echo "dest.even: ", n2.e
  echo "dest.odd: ", n2.o
  threads:
    r := 0
  w.D(r, dest, m)
  #w.Ddag(r, dest, m)
  threads:
    r -= src
  echo "diff: ", r.norm2
  echo "diff even: ", r.even.norm2
  echo "diff odd: ", r.odd.norm2
  #echo eval(r{0})
  #echo eval(dest{0})

  threads:
    dest := 0
  qudaSolveMg(w, dest, src, m, sp)
  #qudaSolve(w, dest, src, m, sp)
  n2.a = dest.norm2
  n2.e = dest.even.norm2
  n2.o = dest.odd.norm2
  echo "dest.norm2: ", n2.a
  echo "dest.even: ", n2.e
  echo "dest.odd: ", n2.o
  threads:
    r := 0
  w.D(r, dest, m)
  #w.Ddag(r, dest, m)
  threads:
    r -= src
  echo "diff: ", r.norm2
  echo "diff even: ", r.even.norm2
  echo "diff odd: ", r.odd.norm2

  #[
  s.solve(destG, src, m, res, sloppySolve = sloppy)
  s.solve(destG, src, m, res, sloppySolve = sloppy)
  s.solve(destG, src, m, res, sloppySolve = sloppy)

  threads:
    echo "destG.norm2: ", destG.norm2
    echo "destG.even: ", destG.even.norm2
    echo "destG.odd: ", destG.odd.norm2
    s.D(r, destG, m)
    threadBarrier()
    r := src - r
    threadBarrier()
    echo "r.norm2: ", r.norm2
    echo "r.even: ", r.even.norm2
    echo "r.odd: ", r.odd.norm2
    r := destG - dest
    threadBarrier()
    echo "gpu-cpu: ", r.norm2 / n2.a
    echo "gpu-cpu:even ", r.even.norm2 / n2.e
    echo "gpu-cpu:odd ", r.odd.norm2 / n2.o
  ]#

  #echoTimers()
  qexFinalize()
