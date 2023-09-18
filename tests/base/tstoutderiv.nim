import qex, physics/qcdTypes, algorithms/numdiff, gauge/stoutsmear
import os, sequtils

proc printPlaq(g: auto) =
  let
    p = g.plaq
    sp = 2.0*(p[0]+p[1]+p[2])
    tp = 2.0*(p[3]+p[4]+p[5])
  echo "plaq ",p
  echo "plaq ss: ",sp," st: ",tp," tot: ",p.sum

qexInit()

let
  fn = if paramCount() > 0: paramStr 1 else: ""
  lat = if fn.len == 0: @[8,8,8,8] else: fn.getFileLattice
  lo = lat.newLayout
var g = lo.newGauge
if fn.len == 0:
  g.random
elif 0 != g.loadGauge fn:
  echo "ERROR: couldn't load gauge file: ",fn
  qexFinalize()
  quit(-1)
g.printPlaq

let gc = GaugeActionCoeffs(plaq:6.0)
var ss = lo.newStoutSmear(0.1)

for n in 0..<10:
  ss.smear(g, g)

g.printPlaq
echo "S(g): ",gc.gaugeAction1(g)

var
  f = lo.newGauge
  p = lo.newGauge
  r = lo.newRNGField(MRG32k3a, 4321)

# SU(3) derivative convention, with Tr(T_a T_b) = -1/2
# F'(U) = -2 T_a F'_a(U) = -2 T_a d/dw_a F(exp(w_b T_b) U) |_{w_a=0}
# Tr(F'(U) T_a) = d/dw_a F(exp(w_b T_b) U) |_{w_a=0}

proc addnoise[G](x:float, p:G, g:G, ng:G):auto =
  for mu in 0..<g.len:
    for e in g[mu]:
      let t = x*p[mu][e]
      ng[mu][e] := exp(t)*g[mu][e]
  ng

# Test by computing the derivative of F(exp(t*p)g), dF/dt at t=0
# d/dt F(exp(t*p)g) |_{t=0} = d/dt F(exp(t*p_b T_b)g) |_{t=0}
#   = p_a d/d(t*p_a) F(exp(t*p_b T_b)g) |_{t=0}
#   = p_a Tr(T_a F'(g))
#   = Tr(p F'(g))
template test(action:untyped, deriv:untyped):auto =
  proc testImpl():auto {.gensym.} =
    tic("test")
    var fail = 0
    echo "### Testing ",astToStr(action)," ",astToStr(deriv)
    deriv(g, f)
    let gg = lo.newGauge
    for n in 0..<5:
      p.randomTAH r
      var d,e:float
      proc testAct(x:float):float =
        tic("testAct")
        result = action(addnoise(x, p, g, gg))
        toc("done")
      ndiff(d, e, testAct, 0, 1.0)
      var pf = 0.0
      for mu in 0..<p.len:
        pf += redot(p[mu], f[mu])
      if abs(pf-d)<32*e:    # check for 32 times the estimated error from ndiff
        echo "Test ",n,"  Passed:  p.f: ",pf," \tndiff: ",d," \tdelta: ",pf-d," \terr(ndiff): ",e
      else:
        echo "Test ",n,"  Failed:  p.f: ",pf," \tndiff: ",d," \tdelta: ",pf-d," \terr(ndiff): ",e
        inc fail
    toc("done")
    fail
  qexGC()
  testImpl()

var fail = test(gc.gaugeAction1, gc.gaugeForce)

proc smearTest0[G](ss:var StoutSmear, gf:G, fl:G) =
  let
    alpha = -ss.alpha  # negative from gaugeForce
    f = ss.f
    ds = ss.ds
  ss.gf = gf
  gaugeActionDeriv(GaugeActionCoeffs(plaq:1.0), gf, ds)
  threads:
    for mu in 0..<f.len:
      for e in f[mu]:
        let s = gf[mu][e]*ds[mu][e].adj
        var t{.noinit.}: evalType(f[mu][e])
        t.projectTAH s
        f[mu][e] := t
        fl[mu][e] := exp(alpha*t)

proc smearTest0Deriv[G](ss:StoutSmear, deriv:G, chain:G) =
  let
    alpha = -ss.alpha  # negative from gaugeForce
    f = ss.f

  threads:
    for mu in 0..<f.len:
      for e in deriv[mu]:
        deriv[mu][e] := alpha*expDeriv(alpha*f[mu][e], chain[mu][e])
  gaugeForceDeriv(ss.gf, deriv, deriv, ss.ds, ss.cg)

proc smearedActionTest0(g:auto):auto =
  var sg = lo.newGauge
  ss.smearTest0(g, sg)
  gc.gaugeAction1(sg)
proc smearedForceTest0(g:auto, f:auto) =
  var sg = lo.newGauge
  var ds = lo.newGauge
  ss.smearTest0(g, sg)
  gc.gaugeActionDeriv(sg, ds)
  ss.smearTest0Deriv(f, ds)
  contractProjectTAH(g, f)

fail += test(smearedActionTest0, smearedForceTest0)

proc smearedAction(g:auto):auto =
  var sg = lo.newGauge
  ss.smear(g, sg)
  gc.gaugeAction1(sg)
proc smearedForce(g:auto, f:auto) =
  var sg = lo.newGauge
  var ds = lo.newGauge
  ss.smear(g, sg)
  gc.gaugeActionDeriv(sg, ds)
  ss.smearDeriv(f, ds)
  contractProjectTAH(g, f)

fail += test(smearedAction, smearedForce)

var s2 = lo.newStoutSmear(0.09)

proc smeared2Action(g:auto):auto =
  var sg = lo.newGauge
  var s2g = lo.newGauge
  ss.smear(g, sg)
  s2.smear(sg, s2g)
  gc.gaugeAction1(s2g)
proc smeared2Force(g:auto, f:auto) =
  var sg = lo.newGauge
  var ds = lo.newGauge
  var s2g = lo.newGauge
  var f2 = lo.newGauge
  ss.smear(g, sg)
  s2.smear(sg, s2g)
  gc.gaugeActionDeriv(s2g, ds)
  s2.smearDeriv(f2, ds)
  ss.smearDeriv(f, f2)
  contractProjectTAH(g, f)

fail += test(smeared2Action, smeared2Force)

# echoTimers()

if fail==0:
  qexFinalize()
else:
  qexAbort(fail)
