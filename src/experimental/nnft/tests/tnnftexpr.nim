## The composed stage pins the numerical fused map and its joint pullback.
import base/globals
setDefaultNc(1)
setVLENmax(4)
import qex
import std/[unittest, math, tables, sequtils]
import ../../graph/[core, scalar, gauge, functional, plan]
import ../../graph/nn as neural
import ../../graph/gauge/types
import ../[flow, model, expr, graph]
import ../../../nn as numeric
import algorithms/numdiff

proc parameters[T: SomeFloat](): NnftParams[T] =
  for s in 0..<8:
    var net: NnftNet[T]
    for l, cin in [nnftFeatures, 12]:
      var w = newSeq[T](12*cin*9)
      var b = newSeq[T](12)
      for j in 0..<w.len:
        w[j] = (if l == 0: T(0.037*sin(float(j)*0.29+float(s)*0.13)) else: T(0.025*cos(float(j)*0.17+float(s)*0.19)))
      for j in 0..<12:
        b[j] = (if l == 0: T(0.13*cos(float(j+s)*0.3)) else: T(0.11*sin(float(j-s)*0.7)))
      net.layers.add convParams(cin,12,[3,3],w,b)
    net.scale = newSeq[T](nnftCoefs)
    for j in 0..<nnftCoefs: net.scale[j] = T(0.7+0.2*cos(float(j+2*s)*0.4))
    result.add net

proc difference(a,b: seq[DLatticeColorMatrixV]): float =
  for mu in 0..<a.len:
    for s in 0..<a[mu].l.nSites:
      var ar,ai,br,bi: float
      ar := a[mu]{s}[0,0].re
      ai := a[mu]{s}[0,0].im
      br := b[mu]{s}[0,0].re
      bi := b[mu]{s}[0,0].im
      result = max(result,max(abs(ar-br),abs(ai-bi)))
  rankMax(result)

proc run[T: SomeFloat]() =
  let p = parameters[T]()
  let lo = newLayout(@[8,8])
  let g = lo.newGauge
  let seed = lo.newGauge
  for mu in 0..1:
    for s in 0..<lo.nSites:
      let a = 0.19*sin(float(2*lo.coords[0][s]-lo.coords[1][s])+0.41*float(mu))
      g[mu]{s}[0,0].re := cos(a)
      g[mu]{s}[0,0].im := sin(a)
      seed[mu]{s}[0,0].re := 0.13*cos(float(s)*0.31+float(mu))
      seed[mu]{s}[0,0].im := 0.17*sin(float(s)*0.23+float(mu))
  let rt = initGraphRuntime()
  let model = toNnftModel(rt,p)
  let w = gauge.toGvalue(rt,g)
  let u = gauge.toGvalue(rt,seed)
  let tol = when T is float32: 3e-6 else: 3e-12

  suite "composed learned stage " & $T:
    test "all stage fields, log determinants and joint pullbacks agree":
      for s in 0..<8:
        let direct = newNnftStage(g,p[s],s)
        discard evalStage(direct,p[s],g)
        let stage = learnedStageGraph(w,model[s],s)
        let loss = redot(u,stage.Wnew)-0.37*stage.lj
        let grad = grad(loss,w)
        check stage.Wnew.runCount == 0
        discard stage.Wnew.eval
        check difference(stage.Wnew.gaugeSnapshot,direct.output) < tol
        check abs(stage.lj.eval.sval-direct.logdet) < tol
        check abs(logDetJ(stage.Wnew,w).eval.sval-direct.logdet) < tol
        let dx = lo.newGauge
        stageVjp(direct,p[s],seed,-0.37,dx)
        discard grad.eval
        check difference(grad.gaugeSnapshot,dx) < tol

    test "the effective force differentiates the effective action":
      # d/dt S_eff(exp(t p) w) at t = 0 is redot(p, F), F = projectTAH(dS_eff/dw w^dag).
      let beta = actWilson(scalar.toGvalue(rt,3.0))
      proc seff(v: Ggauge): Gscalar =
        let f = learnedFlow(v,model)
        gaugeAction(beta,f)-logDetJ(f,v)
      let pdir = lo.newGauge
      for mu in 0..1:
        for s in 0..<lo.nSites:
          pdir[mu]{s}[0,0].re := 0.0
          pdir[mu]{s}[0,0].im := 0.7*cos(float(s)*0.37+1.3*float(mu))
      let pg = gauge.toGvalue(rt,pdir)
      let pdotf = redot(pg,contractProjTAH(grad(seff(w),w),w)).eval.sval
      let t = scalar.toGvalue(rt,0.0)
      let moved = seff(axexpmuly(t,pg,w))
      proc st(v: float): float =
        t.update v
        moved.eval.sval
      var dsdt, err: float
      ndiff(dsdt,err,st,0.0,0.05,ordMax = 5)
      echo "effective force: ndiff ", dsdt, " +/- ", err, " redot ", pdotf
      check abs(dsdt-pdotf) <= (when T is float32: 2e-5 else: 1e-10)*max(1.0,abs(pdotf))

    test "cloning retains factorization and live model parameters":
      let formal = Ggauge(w.newOneOf)
      let stage = learnedStageGraph(formal,model[0],0)
      let loss = redot(u,stage.Wnew)-stage.lj
      let applied = Gscalar(apply(lambda(formal,loss),w))
      let fieldCall = Ggauge(apply(lambda(formal,stage.Wnew),w))
      discard fieldCall.eval
      let cloned = Ggauge(rt.functional.applyCacheByNode[fieldCall.stableNodeId].instantiated)
      let clonedLd = logDetJ(cloned,w)
      let direct = learnedStageGraph(w,model[0],0)
      let actual = redot(u,direct.Wnew)-direct.lj
      check abs(applied.eval.sval-actual.eval.sval) < tol
      check abs(clonedLd.eval.sval-direct.lj.eval.sval) < tol
      let old = actual.sval
      var bias = model[0].biases[1].data
      bias[0] += T(0.31)
      model[0].biases[1].update(bias)
      check abs(applied.eval.sval-actual.eval.sval) < tol
      check abs(clonedLd.eval.sval-direct.lj.eval.sval) < tol
      check abs(actual.sval-old) > 1e-6
      model.update(p)

    test "factorization rejects parameters depending on the flow input":
      var p = model[0]
      p.biases[1] = neural.`*`(0.001*retr(w),p.biases[1])
      expect GraphValueError:
        discard logDetJ(learnedStage(w,p,0).Wnew,w)
      expect GraphValueError:
        discard logDetJ(learnedStageGraph(w,p,0).Wnew,w)

    test "input Hessian agrees with finite differences of the gradient":
      let stage = learnedStageGraph(w,model[0],0)
      let loss = redot(u,stage.Wnew)-0.37*stage.lj
      let dw = grad(loss,w)
      let hv = grad(redot(u,dw),w)
      discard hv.eval
      let want = hv.gaugeSnapshot
      let gp = lo.newGauge
      let gm = lo.newGauge
      let h = when T is float32: 0.003 else: 1e-5
      let fdTol = when T is float32: 3e-4 else: 2e-7
      threads:
        for mu in 0..1:
          gp[mu] := g[mu]+h*seed[mu]
          gm[mu] := g[mu]-h*seed[mu]
      w.update(gp)
      let plus = dw.eval.gaugeSnapshot
      w.update(gm)
      let minus = dw.eval.gaugeSnapshot
      threads:
        for mu in 0..1: plus[mu] := (0.5/h)*(plus[mu]-minus[mu])
      check difference(want,plus) < fdTol
      w.update(g)

    test "fused derivative fallback preserves independent and live seeds":
      let fused = learnedStage(w,model[0],0)
      let expr = learnedStageGraph(w,model[0],0)
      let fixed = scalar.toGvalue(rt,-0.37)
      for k in 0..2:
        let seed = if k == 0: u elif k == 1: w else: w*w
        let dl = if k == 2: 0.003*retr(w) else: fixed
        let a = redot(seed,fused.Wnew)+dl*fused.lj
        let b = redot(seed,expr.Wnew)+dl*expr.lj
        let da = grad(a,w)
        let db = grad(b,w)
        check difference(da.eval.gaugeSnapshot,db.eval.gaugeSnapshot) < 8*tol
        let ha = grad(redot(u,da),w)
        let hb = grad(redot(u,db),w)
        check difference(ha.eval.gaugeSnapshot,hb.eval.gaugeSnapshot) < 8*tol

    test "plans retain joint values, gradients and derivative replicas":
      let fused = learnedStage(w,model[1],1)
      let expr = learnedStageGraph(w,model[1],1)
      let a = redot(u,fused.Wnew)-0.37*fused.lj
      let b = redot(u,expr.Wnew)-0.37*expr.lj
      let da = grad(a,w)
      let db = grad(b,w)
      let ha = grad(redot(u,da),w)
      let hb = grad(redot(u,db),w)
      let run = plan(fused.Wnew,expr.Wnew,a,b,da,db,ha,hb)
      for k in 0..1:
        discard run.eval
        check difference(Ggauge(run[0]).gval,Ggauge(run[1]).gval) < tol
        check abs(Gscalar(run[2]).sval-Gscalar(run[3]).sval) < tol
        check difference(Ggauge(run[4]).gval,Ggauge(run[5]).gval) < 8*tol
        check difference(Ggauge(run[6]).gval,Ggauge(run[7]).gval) < 8*tol
        let want = ha.eval.gaugeSnapshot
        check difference(Ggauge(run[6]).gval,want) < 8*tol
        run.clear

proc runParameters[T: SomeFloat]() =
  let lo = newLayout(@[4,8])
  let g = lo.newGauge
  let seed = lo.newGauge
  for mu in 0..1:
    for s in 0..<lo.nSites:
      let a = 0.23*sin(float(lo.coords[0][s])+0.39*float(lo.coords[1][s])+float(mu))
      g[mu]{s}[0,0].re := cos(a)
      g[mu]{s}[0,0].im := sin(a)
      seed[mu]{s}[0,0].re := 0.21*cos(float(s)*0.31+float(mu))
      seed[mu]{s}[0,0].im := 0.17*sin(float(s)*0.23+float(mu))
  let rt = initGraphRuntime()
  let w = gauge.toGvalue(rt,g)
  let u = gauge.toGvalue(rt,seed)
  let model = toNnftModel(rt,parameters[T]())
  let tol = when T is float32: 3e-6 else: 3e-11
  let h = when T is float32: 0.004 else: 0.0002
  let fdTol = when T is float32: 2e-6 else: 1e-9
  let fused = learnedStage(w,model[0],0)
  let expr = learnedStageGraph(w,model[0],0)
  let a = redot(u,fused.Wnew)-0.37*fused.lj
  let b = redot(u,expr.Wnew)-0.37*expr.lj
  let da = redot(u,grad(a,w))
  let db = redot(u,grad(b,w))

  proc arrayError(x,y: Garray[T]): float =
    discard x.eval
    discard y.eval
    for k in 0..<x.data.len:
      result = max(result,abs(float(x.data[k])-float(y.data[k])))

  suite "learned parameter pullbacks " & $T:
    test "parameter slots retain mixed input derivatives and finite difference agreement":
      for j,leaf in model[0].inputs:
        let p = Garray[T](leaf)
        let values = p.data.mapIt(it)
        var dir = newSeq[T](values.len)
        for k in 0..<dir.len: dir[k] = T(0.07*cos(float(k)*0.13+float(j)*0.29))
        let v = neural.toGarray(rt,dir,p.shape)
        let ga = grad(a,p)
        let mixed = grad(da,p)
        check arrayError(ga,grad(b,p)) < tol
        check arrayError(mixed,grad(db,p)) < tol
        let first = redot(u,grad(neural.redot(v,ga),w))
        let last = neural.redot(v,mixed)
        check abs(first.eval.sval-last.eval.sval) < tol
        let want = last.eval.sval
        p.update(toSeq(0..<values.len).mapIt(values[it]+T(h)*dir[it]))
        let plus = da.eval.sval
        p.update(toSeq(0..<values.len).mapIt(values[it]-T(h)*dir[it]))
        let minus = da.eval.sval
        p.update(values)
        check abs(want-(plus-minus)/(2*h)) < fdTol

    test "shared parameter slots and dependent seeds preserve mixed pullbacks":
      var p = model[1]
      p.biases[0] = p.biases[1]
      p.scale = p.biases[1]
      let v = neural.toGarray(rt,newSeqWith(12,T(0.1)),[12])
      let weight = neural.redot(v,p.biases[1])
      let bg = weight*w+w*w
      let bl = weight+0.003*retr(w)
      let f = learnedStage(w,p,1)
      let e = learnedStageGraph(w,p,1)
      var x = redot(bg,f.Wnew)+bl*f.lj
      var y = redot(bg,e.Wnew)+bl*e.lj
      for order in 0..1:
        check arrayError(grad(x,p.biases[1]),grad(y,p.biases[1])) < tol
        check difference(grad(x,w).eval.gval,grad(y,w).eval.gval) < tol
        if order == 0:
          x = redot(u,grad(x,w))
          y = redot(u,grad(y,w))

    test "cloned parameter pullbacks retain live leaves and plan storage":
      let p = model[2].weights[1]
      let form = Ggauge(w.newOneOf)
      let f = learnedStage(form,model[2],2)
      let score = redot(u,f.Wnew)-0.37*f.lj
      let mixed = grad(redot(u,grad(score,form)),p)
      let call = Garray[T](apply(lambda(form,mixed),w))
      let direct = learnedStageGraph(w,model[2],2)
      let refLoss = redot(u,direct.Wnew)-0.37*direct.lj
      let want = grad(redot(u,grad(refLoss,w)),p)
      let run = plan(call,want)
      let saved = p.data.mapIt(it)
      for k in 0..1:
        if k == 1: p.update(saved.mapIt(T(0.93)*it))
        check arrayError(call,want) < tol
        discard run.eval
        check arrayError(Garray[T](run[0]),Garray[T](run[1])) < tol
      p.update(saved)
      run.clear

  rt.resetCaches

qexInit()
run[float32]()
run[float64]()
runParameters[float32]()
runParameters[float64]()
qexFinalize()
