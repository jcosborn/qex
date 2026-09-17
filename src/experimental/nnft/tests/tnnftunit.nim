## Learned-stage contracts using synthetic parameters, without checkpoint fixtures.
import base/globals
setDefaultNc(1)
setVLENmax(4)

import qex
import ../../../nn
import ../../graph/[core, scalar, gauge, functional, plan]
from ../../graph/nn/types import toGarray
import ../[flow, model, graph]
import std/[math, tables, unittest]

proc parameters[T: SomeFloat](zero = false): NnftParams[T] =
  ## Eight stages of 6 -> 12 -> 12 channels with 3x3 kernels.
  for s in 0..<8:
    var net: NnftNet[T]
    for l, cin in [nnftFeatures, 12]:
      var w = newSeq[T](12*cin*9)
      var b = newSeq[T](12)
      for i in 0..<w.len: w[i] = (if l == 0: T(0.017*sin(float(i+3*s))) else: T(0.023*cos(float(2*i+s))))
      for i in 0..<12: b[i] = (if l == 0: T(0.1*cos(float(i+s))) else: T(0.13*sin(float(i+1))))
      net.layers.add convParams(cin,12,[3,3],w,b)
    net.scale = newSeq[T](nnftCoefs)
    for i in 0..<nnftCoefs: net.scale[i] = (if zero: T(0) else: T(0.4))
    result.add net

proc sample(lo: Layout[VLEN]; shift = 0.0): seq[DLatticeColorMatrixV] =
  result = lo.newGauge
  for d in 0..1:
    for x in 0..<lo.nSites:
      let a = shift+0.43*sin(0.7*float(lo.coords[0][x])+0.4*float(lo.coords[1][x])+float(d))
      result[d]{x}[0,0].re := cos(a)
      result[d]{x}[0,0].im := sin(a)

proc difference(a, b: seq[DLatticeColorMatrixV]): float =
  let lo = a[0].l
  for d in 0..<a.len:
    for x in 0..<lo.nSites:
      var ar, ai, br, bi: float64
      ar := a[d]{x}[0,0].re
      ai := a[d]{x}[0,0].im
      br := b[d]{x}[0,0].re
      bi := b[d]{x}[0,0].im
      result = max(result,max(abs(ar-br),abs(ai-bi)))
  rankMax(result)

proc forwards(rt: GraphRuntime): int =
  for stat in rt.runStatsByNode.values:
    if stat.name == "learnedStageForward": result += stat.count

proc run[T: SomeFloat]() =
  suite "learned stages " & $T:
    let lo = newLayout(@[8,12])
    let g = sample(lo)
    let p = parameters[T]()

    test "zero coefficients give identity, zero logdet, and identity pullback":
      let z = parameters[T](zero=true)
      let dy = sample(lo,0.3)
      let dx = lo.newGauge
      let other = sample(lo,-0.2)
      for s in 0..<8:
        let stage = newNnftStage(g,z[s],s)
        check evalStage(stage,z[s],g) == 0.0
        check difference(stage.output,g) == 0.0
        stageVjp(stage,z[s],dy,2.0,dx)
        check difference(dx,dy) == 0.0
        stageVjp(stage,z[s],other,-0.9,dx)
        check difference(dx,other) == 0.0

    test "inactive links are preserved and active links do not change coefficients":
      var cover = newSeq[int](2*lo.nSites)
      for s in 0..<8:
        let t = newNnftStage(g,p[s],s)
        discard evalStage(t,p[s],g)
        let d = s div 4
        let p0 = (s mod 4) div 2
        let p1 = s mod 2
        for x in 0..<lo.nSites:
          let row = lo.coords[0][x] mod 2
          let col = lo.coords[1][x] mod 2
          let active = row == p0 and col == p1
          check selected(t.active,x) == active
          if active: inc cover[d*lo.nSites+x]
          check selected(t.featureMasks[0],x) == (if d == 0: row != p0 else: col != p1)
          for mu in 0..1:
            check selected(t.featureMasks[mu+1],x) == (mu == 1-d and not active)
        var coef = newSeq[RealField[T]](12)
        for k in 0..<12:
          coef[k] = t.coef[k].newOneOf
          threads: coef[k] := t.coef[k]
        for mu in 0..1:
          for x in 0..<lo.nSites:
            if mu != d or not selected(t.active,x):
              check norm2(t.output[mu]{x}-g[mu]{x}) == 0.0
        let changed = sample(lo)
        for x in 0..<lo.nSites:
          if selected(t.active,x):
            var re, im: float64
            re := changed[d]{x}[0,0].re
            im := changed[d]{x}[0,0].im
            changed[d]{x}[0,0].re := cos(0.27)*re-sin(0.27)*im
            changed[d]{x}[0,0].im := sin(0.27)*re+cos(0.27)*im
        discard evalStage(t,p[s],changed)
        for k in 0..<12:
          for x in 0..<lo.nSites:
            var old, cur: T
            old := coef[k]{x}
            cur := t.coef[k]{x}
            check cur == old
        discard evalStage(t,p[s],g)
        let fresh = newNnftStage(g,p[s],s)
        discard evalStage(fresh,p[s],g)
        check difference(t.output,fresh.output) == 0.0
        check t.logdet == fresh.logdet
        let dy = sample(lo,0.3)
        let dx = lo.newGauge
        for x in 0..<lo.nSites:
          if selected(t.active,x): dy[d]{x} := 0
        stageVjp(t,p[s],dy,0.0,dx)
        check difference(dx,dy) == 0.0
      for n in cover: check n == 1

    test "saturated coefficients use the clipped log determinant and boundary slope":
      let cold = lo.newGauge
      threads:
        for f in cold: f := 1
      var bias, scales = newSeq[T](12)
      for i in 0..<12:
        bias[i] = T(1e30)
        scales[i] = T(1)
      let sat = NnftNet[T](layers: @[
        convParams(6,12,[3,3],newSeq[T](12*6*9),newSeq[T](12)),
        convParams(12,12,[3,3],newSeq[T](12*12*9),bias)],scale:scales)
      let stage = newNnftStage(cold,sat,0)
      let ld = evalStage(stage,sat,cold)
      var jac: float
      jac := stage.m{0}[0,0].re
      check abs(ld-float(lo.physVol div 4)*ln(max(1+jac,1e-8))) < 1e-11
      for j in [-2.0,-1.0,-1.0+1e-9]:
        check learnedLogJac(j) == ln(1e-8)
      for j in [-1.0+1e-8,-0.9999,0.0,0.4]:
        check learnedLogJac(j) == ln(max(1+j,1e-8))

    test "field and logdet share a stage and leaf updates refresh it":
      let rt = initGraphRuntime()
      let v = gauge.toGvalue(rt,g)
      let params = toNnftModel(rt,p)
      let step = learnedStage(v,params[3],3)
      check rt.forwards == 0
      let f = grad(step.lj,v)
      check rt.forwards == 0
      discard step.Wnew.eval
      check rt.forwards == 1
      discard step.lj.eval
      discard f.eval
      check rt.forwards == 1
      let original = step.Wnew.gaugeSnapshot
      discard step.Wnew.eval
      check rt.forwards == 1
      let changed = sample(lo,0.17)
      v.update(changed)
      discard step.Wnew.eval
      discard step.lj.eval
      check rt.forwards == 2
      let numerical = newNnftStage(changed,p[3],3)
      discard evalStage(numerical,p[3],changed)
      check difference(step.Wnew.gval,numerical.output) == 0.0
      check difference(step.Wnew.gval,original) > 1e-3
      params.update(parameters[T](zero=true))
      discard step.Wnew.eval
      check difference(step.Wnew.gval,changed) == 0.0
      check step.lj.eval.sval == 0.0
      check rt.forwards == 3

    test "construction defers tapes and released storage rebuilds from inputs":
      let rt = initGraphRuntime()
      let v = gauge.toGvalue(rt,g)
      let params = toNnftModel(rt,p)
      let before = getRawMemAllocated()
      let stage = learnedStage(v,params[2],2)
      let f = grad(stage.lj,v)
      let run = plan(stage.Wnew,stage.lj,f)
      check getRawMemAllocated() == before
      let base = stage.Wnew.inputs[1]
      check not base.hasStorage
      check not stage.Wnew.hasStorage
      let field = stage.Wnew.eval.gaugeSnapshot
      let ld = stage.lj.eval.sval
      let force = f.eval.gaugeSnapshot
      for k in 0..1:
        if k == 0: base.releaseWork
        else: base.releaseStorage
        check not base.hasStorage
        discard base.eval
        check difference(stage.Wnew.eval.gval,field) == 0.0
        check stage.lj.eval.sval == ld
        check difference(f.eval.gval,force) == 0.0
      check rt.forwards == 3
      discard run.eval
      check difference(Ggauge(run[0]).gval,field) == 0.0
      check Gscalar(run[1]).sval == ld
      check difference(Ggauge(run[2]).gval,force) == 0.0
      run.clear

    test "joint plans refresh fields and forces and rebuild after clear":
      let rt = initGraphRuntime()
      let v = gauge.toGvalue(rt,g)
      let params = toNnftModel(rt,p)
      let stage = learnedStage(v,params[5],5)
      let loss = norm2(stage.Wnew)-0.37*stage.lj
      let f = grad(loss,v)
      let run = plan(stage.Wnew,loss,f)
      for k in 0..2:
        if k == 1: v.update(sample(lo,0.23))
        if k == 2: params.update(parameters[T](zero=true))
        let field = stage.Wnew.eval.gaugeSnapshot
        let value = loss.eval.sval
        let force = f.eval.gaugeSnapshot
        discard run.eval
        check difference(Ggauge(run[0]).gval,field) == 0.0
        check Gscalar(run[1]).sval == value
        check difference(Ggauge(run[2]).gval,force) == 0.0
        let forwards = run.stats.forwards
        discard run.eval
        check run.stats.forwards == forwards
        let old = run[0]
        run.clear
        discard run.eval
        check old.stale
        check old != run[0]
        check difference(Ggauge(run[0]).gval,field) == 0.0
        check difference(Ggauge(run[2]).gval,force) == 0.0
      run.clear

    test "lambda clones retain independent tapes and read current leaves":
      let rt = initGraphRuntime()
      let params = toNnftModel(rt,p)
      let formal = gauge.toGvalue(rt,g)
      let fun = lambda(formal,learnedStage(formal,params[0],0).Wnew)
      let a = gauge.toGvalue(rt,g)
      let gb = sample(lo,0.23)
      let b = gauge.toGvalue(rt,gb)
      let fa = Ggauge(apply(fun,a))
      let fb = Ggauge(apply(fun,b))
      discard fa.eval
      let saved = fa.gaugeSnapshot
      discard fb.eval
      check difference(fa.gval,saved) == 0.0
      check difference(fa.gval,fb.gval) > 1e-3
      let run = plan(fa,fb)
      discard run.eval
      check difference(Ggauge(run[0]).gval,fa.gval) == 0.0
      check difference(Ggauge(run[1]).gval,fb.gval) == 0.0
      formal.update(gb)
      discard fa.eval
      check difference(fa.gval,saved) == 0.0
      a.update(gb)
      discard fa.eval
      check difference(fa.gval,fb.gval) == 0.0
      discard run.eval
      check difference(Ggauge(run[0]).gval,fa.gval) == 0.0
      check difference(Ggauge(run[1]).gval,fb.gval) == 0.0
      params.update(parameters[T](zero=true))
      discard fa.eval
      discard fb.eval
      check difference(fa.gval,gb) == 0.0
      check difference(fb.gval,gb) == 0.0
      run.clear

    test "wrong stage, dimension, extent, and layouts fail at construction":
      let rt = initGraphRuntime()
      let params = toNnftModel(rt,p)
      let v = gauge.toGvalue(rt,g)
      expect GraphValueError: discard learnedStage(v,params[0],-1)
      check stageClass(8) == stageClass(0)
      discard learnedStage(v,params[0],8)
      let three = newLayout(@[4,4,4]).newGauge
      expect ValueError: discard learnedStage(gauge.toGvalue(rt,three),params[0],0)
      when DLatticeColorMatrixV.V == 1:
        let small = newLayout(@[2,8]).newGauge
        expect ValueError: discard learnedStage(gauge.toGvalue(rt,small),params[0],0)
      let other = newLayout(@[8,12]).newGauge
      let mixed = @[g[0],other[1]]
      expect ValueError: discard learnedStage(gauge.toGvalue(rt,mixed),params[0],0)
      check rt.forwards == 0

    test "parameter shapes and runtimes are checked before differentiation":
      let rt = initGraphRuntime()
      let v = gauge.toGvalue(rt,g)
      let params = toNnftModel(rt,p)
      var wrong = params[0]
      wrong.weights[0] = toGarray(rt,p[0].layers[0].weights,[6,12,3,3])
      expect GraphValueError: discard learnedStage(v,wrong,0)
      let other = toNnftModel(initGraphRuntime(),p)
      expect GraphValueError: discard learnedStage(v,other[0],0)
      let stage = learnedStage(v,params[0],0)
      check grad(stage.lj,params[0].biases[1]).copyCompatible(params[0].biases[1])
      check rt.forwards == 0

    test "one graph stage owns one stage of lattice storage":
      let rt = initGraphRuntime()
      let v = gauge.toGvalue(rt,g)
      let params = toNnftModel(rt,p)
      let startStage = getRawMemAllocated()
      let stage = learnedStage(v,params[0],0)
      discard stage.Wnew.eval
      let stageBytes = getRawMemAllocated()-startStage
      let startFlow = getRawMemAllocated()
      let flow = learnedFlow(v,params)
      discard flow.eval
      let flowBytes = getRawMemAllocated()-startFlow
      check stageBytes > 0
      check flowBytes >= 6*stageBytes
      check flowBytes <= 10*stageBytes
      echo "raw allocation bytes: single stage=", stageBytes, ", graph flow=", flowBytes

qexInit()
run[float32]()
run[float64]()
qexFinalize()
