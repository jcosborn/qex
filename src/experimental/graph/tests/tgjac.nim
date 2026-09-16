#RUNCMD env OMP_NUM_THREADS=1 $RUN1
## Small fixtures for stout logdet and grouped derivative closure.
import base/globals
setVLENmax(4)

import math, unittest
import helpers
import ../../../../tests/base/scaledexpRef
import qex except epsilon
import algorithms/numdiff
import maths/groupOps
import ../[core, scalar, gauge]
import ../functional
import ../hmcgauge/ftstout

proc runJacTests*(localLat: seq[int]) =
  qexInit()
  defer: qexFinalize()
  letParam:
    expectRanks = nRanks
  check nRanks == expectRanks
  letParam:
    lat = latticeFromLocalLattice(localLat,nRanks)
  let lo = lat.newLayout
  var rng = lo.newRNGField(Philox4x64, 734991'u64)
  let g = lo.newGauge
  let d = lo.newGauge
  let u = lo.newGauge
  let v = lo.newGauge
  let alt = lo.newGauge
  threads:
    g.random rng
    d.random rng
    u.random rng
    v.random rng
    alt.random rng
    for mu in 0..<g.len:
      d[mu] := 0.3*d[mu] + 0.2*alt[mu]
      u[mu] := 0.2*u[mu] + 0.1*d[mu]
      v[mu] := 0.3*v[mu] - 0.15*alt[mu]

  proc same(x, y: Ggauge, tol=2e-18) =
    let err = norm2(x-y).eval.sval
    let refv = norm2(y).eval.sval
    check err < tol*(1.0+refv)

  proc same(x, y: Gscalar, tol=2e-10) =
    let xv = x.eval.sval
    let yv = y.eval.sval
    check abs(xv-yv) < tol*(1.0+abs(yv))

  proc checkT(name: string, f, t: Gscalar, at=0.0, step=0.003) =
    t.update at
    let an = grad(f,t).eval.sval
    proc val(x: float): float =
      t.update x
      f.eval.sval
    var fd, err: float
    ndiff(fd,err,val,at,step,ordMax=4)
    t.update at
    echo "  ", name, " derivative ", an, ", residual ", abs(an-fd), ", FD estimate ", err
    check abs(an-fd) < 3e-7*max(1.0,abs(fd))

  suite "stout Jacobian closure":
    setup:
      let rt = initGraphRuntime()
      defer:
        rt.resetGradCache(false)
        rt.resetApplyCache(false)
        rt.resetLdjCache
      let W = gauge.toGvalue(rt,g)
      let ds = gauge.toGvalue(rt,d)
      let gu = gauge.toGvalue(rt,u)
      let gv = gauge.toGvalue(rt,v)
      let a = scalar.toGvalue(rt,0.09)
      let b = scalar.toGvalue(rt,0.8)
      let t = scalar.toGvalue(rt,0.0)
      let c = actWilson(b)
      const parity = 1
      const dir = 0

    test "finite polynomial replica matches fused value and first gradients":
      let ff = stoutLogDetJ(W,ds,a,parity,dir)
      let fr = stoutLogDetJGraph(W,ds,a,parity,dir)
      same(ff,fr)
      same(grad(ff,W),grad(fr,W))
      same(grad(ff,ds),grad(fr,ds))
      same(grad(ff,a),grad(fr,a))
      same(maskSubset(parity,dir,grad(ff,W)),grad(ff,W))
      same(maskSubset(parity,dir,grad(ff,ds)),grad(ff,ds))
      a.update 0.06
      W.update alt
      ds.update v
      same(ff,fr)
      same(grad(ff,W),grad(fr,W))
      same(grad(ff,ds),grad(fr,ds))
      same(grad(ff,a),grad(fr,a))

    test "logdet second and third derivatives retain dependent cotangents":
      proc second(x,y: Ggauge, alpha, seed: Gscalar, replica=false): Gscalar =
        let f = seed * (if replica: stoutLogDetJGraph(x,y,alpha,parity,dir)
                       else: stoutLogDetJ(x,y,alpha,parity,dir))
        redot(grad(f,x),x.adj*gu+gv) + redot(grad(f,y),y+gu) + alpha*grad(f,alpha)
      let f2 = second(W+t*gu,ds+t*gv,a+0.03*t,b+t)
      let r2 = second(W+t*gu,ds+t*gv,a+0.03*t,b+t,true)
      same(f2,r2)
      same(grad(f2,t),grad(r2,t))
      checkT("logdet second",f2,t)
      let f3 = redot(grad(second(W,ds,a,b),W),W*gv+gu) + a*grad(second(W,ds,a,b),a)
      checkT("logdet third alpha",f3,a,0.09)

    test "aliased field and scalar cotangents retain all input paths":
      proc second(x: Ggauge, alpha: Gscalar, replica=false): Gscalar =
        let f = alpha * (if replica: stoutLogDetJGraph(x,x,alpha,parity,dir)
                        else: stoutLogDetJ(x,x,alpha,parity,dir))
        redot(grad(f,x),x) + alpha*grad(f,alpha)
      let ff = second(W+t*gu,a+t)
      let fr = second(W+t*gu,a+t,true)
      same(ff,fr)
      same(grad(ff,t),grad(fr,t))
      checkT("aliased logdet",ff,t)

    test "grouped fixed-staple update logdet and combined pullbacks":
      for mode in 0..2:
        proc score(x,y: Ggauge, alpha: Gscalar, grouped: bool): Gscalar =
          var z: Ggauge
          var l: Gscalar
          if grouped:
            let st = stoutUpdateLogDetJ(x,y,alpha,parity,dir)
            z = st.Wnew
            l = st.lj
          else:
            z = stoutUpdate(x,y,alpha,parity,dir)
            l = stoutLogDetJ(x,y,alpha,parity,dir)
          case mode
          of 0: redot(z,gu)
          of 1: l
          else:
            let r = redot(z,gu)
            r*r + redot(z,gv) + redot(z,x) - 0.7*l + 0.2*l*l
        proc second(x,y: Ggauge, alpha: Gscalar, grouped: bool): Gscalar =
          let s = score(x,y,alpha,grouped)
          redot(grad(s,x),x.adj*gu+gv) + redot(grad(s,y),y*gv+gu) + alpha*grad(s,alpha)
        let ff = second(W+t*gu,ds+t*gv,a+0.03*t,true)
        let fr = second(W+t*gu,ds+t*gv,a+0.03*t,false)
        same(ff,fr)
        same(grad(ff,t),grad(fr,t))
        checkT("grouped fixed staple " & $mode,ff,t)

    test "grouped action Hessian preserves W coefficient and alpha dependence":
      for mode in 0..2:
        proc score(x: Ggauge, beta, alpha: Gscalar, fused: bool): Gscalar =
          let coeff = actWilson(beta)
          let st = if fused: stoutUpdateLogDetJ(x,coeff,alpha,parity,dir)
                   else: stoutUpdateLogDetJ(x,gaugeActionDeriv(coeff,x,parity,dir),alpha,parity,dir)
          case mode
          of 0: redot(st.Wnew,gu)
          of 1: st.lj
          else:
            let r = redot(st.Wnew,gu)
            r*r + redot(st.Wnew,gv) + redot(st.Wnew,x) - 0.7*st.lj + 0.2*st.lj*st.lj
        proc second(x: Ggauge, beta, alpha: Gscalar, fused: bool): Gscalar =
          let s = score(x,beta,alpha,fused)
          redot(grad(s,x),x.adj*gu+gv) + beta*grad(s,beta) + alpha*grad(s,alpha)
        let ff = second(W+t*gu,b+0.07*t,a+0.03*t,true)
        let fr = second(W+t*gu,b+0.07*t,a+0.03*t,false)
        same(ff,fr)
        same(grad(ff,t),grad(fr,t))
        checkT("grouped action Hessian " & $mode,ff,t)
      let ff = stoutUpdateLogDetJ(W,c,a,parity,dir)
      let gr = grad(ff.lj,b)
      let refg = grad(stoutLogDetJGraph(W,gaugeActionDeriv(c,W,parity,dir),a,parity,dir),b)
      for beta in [0.6,0.0]:
        b.update beta
        same(gr,refg)
        same(grad(gr,b),grad(refg,b))

    test "higher grouped pullbacks clone with their evaluation caches":
      let p = Ggauge(W.newOneOf)
      proc score(x: Ggauge): Gscalar =
        let st = stoutUpdateLogDetJ(x,c,a,parity,dir)
        let s = redot(st.Wnew,gu)-st.lj
        redot(grad(s,x),x*gv+gu) + a*grad(s,a)
      let fun = lambda(p,score(p))
      let cloned = Gscalar(apply(fun,W))
      let direct = score(W)
      p.update d
      same(cloned,direct)
      same(grad(cloned,W),grad(direct,W))
      a.update 0.06
      b.update 0.7
      W.update alt
      p.update g
      same(cloned,direct)
      same(grad(cloned,W),grad(direct,W))

    test "exposed flow rho participates in mixed force derivatives":
      let sa = stoutAction(c,a,1)
      check sa.rho.nodeKey == a.nodeKey
      let f = sa.action(W)
      let forceScore = redot(grad(f,W),gu)
      checkT("flow force rho",forceScore,a,0.09,0.001)
      let rhoScore = grad(sa.action(W+t*gu),a)
      checkT("flow rho field",rhoScore,t,0.0,0.001)

  when g[0][0].nrows == 3:
    let refs = buildRefs()
    let sl = latticeFromLocalLattice(@[4,4],nRanks).newLayout
    type M = MatrixArray[3,3,ComplexType[float64]]

    proc sample(name: string): CaseRef =
      for c in refs.cases:
        if c.name == name: return c
      doAssert false, "missing scaled differential reference " & name

    proc finite(c: CaseRef): FiniteRef =
      for r in c.finite:
        if r.order == 13 and r.scale == expProjectTAHScale: return r
      doAssert false, "missing scaled finite Phi reference"

    proc fill(dst: typeof(g), m: auto) =
      var value: M
      value := m
      threads:
        for f in dst:
          for e in f:
            f[e][] := value

    let base = sample("rotated8")
    let one = sl.newGauge
    let aux = sl.newGauge
    let um = sl.newGauge
    let vm = sl.newGauge
    var id: M
    id := 1.0
    fill(one,id)
    fill(aux,base.m.adj)
    fill(um,base.dm)
    fill(vm,base.c)

    suite "scaled stout Jacobian":
      setup:
        let rt = initGraphRuntime()
        defer:
          rt.resetGradCache(false)
          rt.resetApplyCache(false)
          rt.resetLdjCache
        let W = gauge.toGvalue(rt,one)
        let ds = gauge.toGvalue(rt,aux)
        let gu = gauge.toGvalue(rt,um)
        let gv = gauge.toGvalue(rt,vm)
        let a = scalar.toGvalue(rt,1.0)
        let t = scalar.toGvalue(rt,0.0)
        let s = scalar.toGvalue(rt,0.0)
        const parity = 1
        const dir = 0

      test "bounded oracle spectra match graph value and ambient gradients":
        let ff = stoutLogDetJ(W,ds,a,parity,dir)
        let fr = stoutLogDetJGraph(W,ds,a,parity,dir)
        let fw = grad(ff,W)
        let fd = grad(ff,ds)
        let fa = grad(ff,a)
        let rw = grad(fr,W)
        let rd = grad(fr,ds)
        let ra = grad(fr,a)
        let dv = sl.newGauge
        let gw = sl.newGauge
        let gd = sl.newGauge
        let sub = sl.getSubset("odd")
        let vol = 0.5*float(sl.physVol)
        for c in refs.cases:
          checkpoint c.name
          let r = finite(c)
          let m = c.m
          let mg = r.grad
          let k = r.k
          check determinant(k) > 0.0
          check r.invNorm < 60.0
          echo "  ", c.name, " Fnorm=", c.fNorm,
            " Dnorm=", c.dNorm, " KinvNorm=", r.invNorm
          fill(dv,m.adj)
          ds.update dv
          threads:
            for mu in 0..<gw.len:
              gw[mu] := 0.0
              gd[mu] := 0.0
            threadBarrier()
            for e in sub:
              gw[dir][e][] := mg*m.adj
              gd[dir][e][] := mg.adj
          let ew = gauge.toGvalue(rt,gw)
          let ed = gauge.toGvalue(rt,gd)
          let el = scalar.toGvalue(rt,vol*r.log)
          let ea = scalar.toGvalue(rt,vol*redot(mg,m))
          same(ff,el,2e-11)
          same(fr,el,2e-11)
          same(fw,ew,4e-22)
          same(rw,ew,4e-22)
          same(fd,ed,4e-22)
          same(rd,ed,4e-22)
          same(fa,ea,2e-11)
          same(ra,ea,2e-11)

      test "norm eight mixed fourth logdet derivatives commute":
        # Both ambient directions fail to commute with the base generator.
        let m = base.m
        let dm = base.dm
        let dv = base.c
        var f: M
        f.projectTAH(m)
        check norm2(f*dm-dm*f) > 1.0
        check norm2(f*dv-dv*f) > 1.0
        let x = W+t*gu
        let y = ds+s*gv
        let alpha = a+0.02*t+0.03*s
        let seed = 1.0+0.07*t+0.11*s
        let ff = seed*stoutLogDetJ(x,y,alpha,parity,dir)
        let fr = seed*stoutLogDetJGraph(x,y,alpha,parity,dir)
        let fts = grad(grad(ff,t),s)
        let rts = grad(grad(fr,t),s)
        let ftst = grad(fts,t)
        let ftsts = grad(ftst,s)
        let fstst = grad(grad(grad(grad(ff,s),t),s),t)
        same(fts,rts,2e-11)
        same(ftsts,fstst,2e-9)
        checkT("norm eight mixed fourth",ftst,s,0.0,0.001)

      test "norm eight grouped caches match separate pullbacks after refresh":
        let st = stoutUpdateLogDetJ(W,ds,a,parity,dir)
        let up = stoutUpdate(W,ds,a,parity,dir)
        let lj = stoutLogDetJ(W,ds,a,parity,dir)
        let fg = redot(st.Wnew,gu)+redot(st.Wnew,gv)-0.4*st.lj
        let fr = redot(up,gu)+redot(up,gv)-0.4*lj
        let gw = grad(fg,W)
        let rw = grad(fr,W)
        let gd = grad(fg,ds)
        let rd = grad(fr,ds)
        let ga = grad(fg,a)
        let ra = grad(fr,a)
        let gm = redot(grad(ga,W),gu)
        let rm = redot(grad(ra,W),gu)
        let wv = sl.newGauge
        let dv = sl.newGauge
        let uv = sl.newGauge
        for pass, alpha in [1.0,0.8,1.0]:
          var wm: M
          wm := 1.0
          let angle = 0.13*float(pass)
          wm[0,0].re = cos(angle)
          wm[0,0].im = sin(angle)
          wm[1,1].re = cos(angle)
          wm[1,1].im = -sin(angle)
          # alpha W ds† = M keeps the oracle's positive K throughout refreshes.
          let c = sample(if pass == 1: "extremal8" else: "rotated8")
          fill(wv,wm)
          fill(dv,(1.0/alpha)*(c.m.adj*wm))
          fill(uv,(1.0+0.1*float(pass))*c.dm)
          W.update wv
          ds.update dv
          gu.update uv
          a.update alpha
          if pass > 0:
            rt.resetGradCache(false)
            rt.resetApplyCache(false)
            rt.resetLdjCache
          for repeat in 0..1:
            same(st.Wnew,up,4e-22)
            same(st.lj,lj,2e-11)
            same(fg,fr,2e-11)
            same(gw,rw,4e-22)
            same(gd,rd,4e-22)
            same(ga,ra,2e-11)
            same(gm,rm,2e-10)

      test "norm eight update field pullback matches directional differences":
        var f, df: M
        f.projectTAH(base.m)
        df.projectTAH(base.c*base.m)
        check abs(base.fNorm-8.0) < 1e-12
        check norm2(f*df-df*f) > 1e-12
        let x = W+t*gv
        let st = stoutUpdateLogDetJ(x,ds,a,parity,dir)
        let first = redot(grad(redot(st.Wnew,gu),x),gv)
        # FD evaluates the fused Phi field pullback; grad uses the ordinary
        # exponential replica. Check their directional agreement near norm eight.
        checkT("norm eight update field pullback",first,t,0.0,0.001)
        checkT("norm eight update field pullback half step",first,t,0.0,0.0005)

      test "update pullbacks match differences across norms and inputs":
        let dv = sl.newGauge
        let uv = sl.newGauge
        let vv = sl.newGauge

        proc checkUpdate(name: string, first: Gscalar) =
          let an = grad(first,t).eval.sval
          proc val(x: float): float =
            t.update x
            first.eval.sval
          for step in [0.001,0.0005]:
            var fd, err: float
            ndiff(fd,err,val,0.0,step,ordMax=4)
            t.update 0.0
            let den = max(1.0,abs(fd))
            echo "  ", name, " step ", step, " derivative ", an,
              ", residual ", abs(an-fd), ", FD estimate ", err
            # Keep the approximation bound; require FD uncertainty below a tenth.
            check abs(an-fd) < 3e-7*den
            check err < 3e-8*den

        for name in ["generic012","generic2","rotated8"]:
          checkpoint name
          let c = sample(name)
          let r = finite(c)
          check determinant(r.k) > 0.0
          check r.invNorm < 60.0
          var f, dw, dd: M
          f.projectTAH(c.m)
          dw.projectTAH(c.c*c.m)
          dd.projectTAH(c.dm.adj)
          check norm2(f*dw-dw*f) > 1e-12
          check norm2(f*dd-dd*f) > 1e-12
          fill(dv,c.m.adj)
          fill(uv,c.dm)
          fill(vv,c.c)
          ds.update dv
          gu.update uv
          gv.update vv

          if name != "rotated8":
            let x = W+t*gv
            let st = stoutUpdateLogDetJ(x,ds,a,parity,dir)
            let first = redot(grad(redot(st.Wnew,gu),x),gv)
            checkUpdate(name & " update field pullback",first)
          block:
            let y = ds+t*gu
            let st = stoutUpdateLogDetJ(W,y,a,parity,dir)
            let first = redot(grad(redot(st.Wnew,gv),y),gu)
            checkUpdate(name & " update staple pullback",first)
          block:
            let alpha = a+t
            let st = stoutUpdateLogDetJ(W,ds,alpha,parity,dir)
            let first = redot(grad(redot(st.Wnew,gu),W),gv)
            checkUpdate(name & " update field pullback alpha",first)

when isMainModule:
  runJacTests(@[4,4])
