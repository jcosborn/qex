import physics/wilsonD

import os
import times
#import strutils
import base
import layout
import field
import qcdTypes
#import stdUtils
import solvers/cg
export cg
import solvers/bicgstab
#import solvers/solverUtils
#import types
#import profile
#import metaUtils
import gauge/gaugeUtils

proc solveEoCg*(s: Wilson; r,x: Field; m: SomeNumber; sp0: var SolverParams) =
  var sp = sp0
  sp.subset.layoutSubset(r.l, sp.subsetName)
  if sp.subsetName=="all": sp.subset.layoutSubset(r.l, "even")
  var t = newOneOf(r)
  var x2 = newOneOf(x)
  var top = 0.0
  let m4 = 4.0 + m
  let m2 = m4*m4
  proc op(a,b:Field) =
    threadBarrier()
    if threadNum==0: top -= epochTime()
    wilsonD2eex(s.se, s.so, t, s.g, b, m2)
    threadBarrier()
    wilsonD2ee(s.se, s.so, a, s.g, t, m2)
    if threadNum==0: top += epochTime()
    threadBarrier()
  threads:
    #echo "x2: ", x.norm2
    s.eoReduce(x2, x, m)
  let t0 = epochTime()
  var oa = (apply: op)
  cgSolve(r, x2, oa, sp)
  threads:
    wilsonD2eex(s.se, s.so, t, s.g, r, m2)
    threadBarrier()
    r[s.se.sub] := t
    threadBarrier()
    s.eoReconstruct(r, x, m)
  let t1 = epochTime()
  let secs = t1-t0
  let flops = (2*2*s.g.len*(12+2*66+24)+2*60)*r.l.nEven*sp.finalIterations
  sp0.finalIterations = sp.finalIterations
  sp0.seconds = secs
  echo "op time: ", top
  echo "solve time: ", secs, "  Gflops: ", 1e-9*flops.float/secs
proc solveEoCg*(s:Wilson; r,x:Field; m:SomeNumber; res:float) =
  var sp = initSolverParams()
  sp.r2req = res
  #sp.maxits = 1000
  sp.verbosity = 1
  solveEO(s, r, x, m, sp)




when isMainModule:
  qexInit()

  var defaultLat = @[8,8,8,8]
  defaultSetup()
  threads:
    g.setBC

  var src1 = lo.DiracFermion()
  var src2 = lo.DiracFermion()
  var soln1 = lo.DiracFermion()
  var soln2 = lo.DiracFermion()
  var v1 = lo.DiracFermion()
  var v2 = lo.DiracFermion()
  var r = lo.DiracFermion()
  var tv = lo.DiracFermion()

  var w = newWilson(g, v1)
  echo "done newWilson"
  var m = floatParam("mass", 0.1)
  echo "mass: ", m

  var sp = initSolverParams()
  sp.maxits = int(1e9/lo.physVol.float)
  threads:
    src1 := 0
  if myRank==0:
    src1{0}[0][0] := 1
  sp.r2req = 1e-16
  sp.verbosity = 1
  w.solve(soln1, src1, m, sp)
  threads:
    soln2 := src1
    threadBarrier()
    w.D(src1, soln1, m)
    w.D(src2, soln2, m)

  proc checkResid(src,soln,x: any) =
    threads:
      #echo "v2: ", v2.norm2
      #echo "v2.even: ", v2.even.norm2
      #echo "v2.odd: ", v2.odd.norm2
      w.D(r, x, m)
      threadBarrier()
      r := src - r
      let r2 = r.norm2
      tv := soln - x
      let e2 = tv.norm2
      echo "r2: ", r2, "  e2: ", e2

      var fr20 = relResid(r, v2, 1e-20)
      var fr30 = relResid(r, v2, 1e-30)
      echo "fnalResid: ", fr20, "  ", fr30

  proc test1(src,soln: any) =
    echo "solve:"
    threads:
      v2 := 0
    #w.solveEO(v2, src, m, sp)
    #resetTimers()
    w.solve(v2, src, m, sp)
    #echoTimers()
    checkResid(src, soln, v2)
    #var fr = relResid(r,v2,1e-30)
    #echo "fnalResid: ", fr

    echo "solveEO:"
    threads:
      v2 := 0
    #resetTimers()
    w.solveEO(v2, src, m, sp)
    #echoTimers()
    checkResid(src, soln, v2)
    #fr = relResid(r,v2,1e-30)
    #echo "fnalResid: ", fr

    echo "solveEoCg:"
    threads:
      v2 := 0
    #resetTimers()
    w.solveEoCg(v2, src, m, sp)
    #echoTimers()
    checkResid(src, soln, v2)
    #fr = relResid(r,v2,1e-30)
    #echo "fnalResid: ", fr

  sp.r2req = 1e-12
  sp.verbosity = 1
  test1(src1, soln1)

  sp.r2req = 1e-16
  test1(src1, soln1)

  #echoTimers()
  qexFinalize()
