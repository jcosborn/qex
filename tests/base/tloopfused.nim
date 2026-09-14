import qex, gauge/plaquette
import base/alignedMem
import std/importutils
import testutils

qexInit()
letParam:
  lat = latticeFromLocalLattice(@[4,4,4,4],nRanks)
let lo=lat.newLayout

proc run(nc:static int) =
  let g=lo.newGauge(nc)
  let gp=lo.newGauge(nc)
  let gm=lo.newGauge(nc)
  let f=lo.newGauge(nc)
  let want=lo.newGauge(nc)
  let fp=lo.newGauge(nc)
  let fm=lo.newGauge(nc)
  let hs=[lo.newGauge(nc),lo.newGauge(nc),lo.newGauge(nc),lo.newGauge(nc),lo.newGauge(nc)]
  let work=newLoopWork(g[0])
  var rng=lo.newRNGField(Philox4x64,20260915'u64)
  threads:
    g.random rng
    for h in hs:
      h.gaussian rng
      for mu in 0..<g.len: h[mu] *= 0.07
    for mu in 0..<g.len: g[mu] += hs[0][mu]
  let cases=[GaugeActionCoeffs(plaq:0.8),GaugeActionCoeffs(rect: -0.1),
    GaugeActionCoeffs(pgm:0.07),GaugeActionCoeffs(plaq:0.8,rect: -0.1),
    GaugeActionCoeffs(plaq:0.8,pgm:0.07),GaugeActionCoeffs(rect: -0.1,pgm:0.07),
    GaugeActionCoeffs(plaq:0.8,rect: -0.1,pgm:0.07)]
  privateAccess(typeof(work))
  privateAccess(typeof(work.plan))
  suite "Fused loop kernels Nc=" & $nc:
    test "values, gradients and Hessians match independent production paths":
      for c in cases:
        let a=c.loopAction(g,work=work)
        check abs(a-c.gaugeAction1(g)) < 1e-11*(1+abs(a))
        check abs(a-c.gaugeAction2(g)) < 1e-11*(1+abs(a))
        c.loopDeriv(g,f,work=work)
        (-1.0*c).gaugeActionDeriv(g,want)
        check relativeDiff(f,want) < 2e-12
        c.loopDeriv(g,[hs[0]],f,work=work)
        threads:
          for mu in 0..<g.len: want[mu] := 0
        c.gaugeDerivDeriv2(g,hs[0],want)
        check relativeDiff(f,want) < 2e-12
    test "higher jets and seed permutations":
      let c=cases[^1]
      forStatic k, 2, 5:
        var ds:array[k,type(g)]
        var prev:array[k-1,type(g)]
        for i in 0..<k: ds[i]=hs[i]
        for i in 0..<k-1: prev[i]=hs[i]
        threads:
          for mu in 0..<g.len:
            gp[mu] := g[mu]+1e-4*hs[k-1][mu]
            gm[mu] := g[mu]-1e-4*hs[k-1][mu]
        c.loopDeriv(g,ds,f,work=work)
        c.loopDeriv(gp,prev,fp,work=work)
        c.loopDeriv(gm,prev,fm,work=work)
        threads:
          for mu in 0..<g.len: want[mu] := 5000.0*(fp[mu]-fm[mu])
        check relativeDiff(f,want) < 2e-8
        swap(ds[0],ds[^1])
        c.loopDeriv(g,ds,want,work=work)
        check relativeDiff(f,want) < 2e-12
    test "plaquette specialization agrees through degree three":
      let c=GaugeActionCoeffs(plaq: -float(nc))
      forStatic k, 0, 3:
        var ds:array[k,type(g)]
        for i in 0..<k: ds[i]=hs[i]
        c.loopDeriv(g,ds,f,work=work)
        stapleSum(g,@ds,want)
        check relativeDiff(f,want) < 2e-12
    test "highest jets, repeated seeds and alias contracts":
      let c=cases[3]
      let ds=[hs[0],hs[1],hs[0],hs[1],hs[0]]
      c.loopDeriv(g,ds,want,work=work)
      c.loopDeriv(gp,ds,gp,work=work)
      check relativeDiff(gp,want) == 0
      c.loopDeriv(g,[hs[0],hs[1],hs[0],hs[1],hs[0],hs[1]],f,work=work)
      threads:
        for mu in 0..<g.len: want[mu] := 0
      check relativeDiff(f,want) == 0
      expect(ValueError): c.loopDeriv(g,[hs[0]],hs[0],work=work)
      expect(ValueError): c.loopDeriv(g,g,work=work)
    test "warm workspaces reuse buffers and refresh field values":
      let c=cases[^1]
      discard c.loopAction(g,work=work)
      c.loopDeriv(g,[hs[0],hs[1],hs[2]],f,work=work)
      let raw=getRawMemAllocated()
      discard c.loopAction(gm,work=work)
      c.loopDeriv(gm,[hs[2],hs[1],hs[0]],f,work=work)
      check getRawMemAllocated() == raw
      c.loopDeriv(gm,[hs[2],hs[1],hs[0]],want)
      check relativeDiff(f,want) == 0

    test "shared products reduce work for active families":
      let planes=lo.nDim*(lo.nDim-1) div 2
      let triples=lo.nDim*(lo.nDim-1)*(lo.nDim-2) div 6
      for c in cases:
        let np=if c.plaq!=0: planes else: 0
        let nr=if c.rect!=0: 2*planes else: 0
        let ng=if c.pgm!=0: 4*triples else: 0
        if np+nr+ng == 0: continue
        discard c.loopAction(g,work=work)
        check work.plan.outs.len == np+nr+ng
        check work.plan.steps.len <= 2*np+4*(nr+ng)
        let actionProducts=work.plan.steps.len
        c.loopDeriv(g,f,work=work)
        check work.plan.outs.len == 4*np+6*(nr+ng)
        check work.plan.steps.len < 2*4*np+4*6*(nr+ng)
        let gradientProducts=work.plan.steps.len
        c.loopDeriv(g,[hs[0]],f,work=work)
        check work.plan.steps.len < 5*4*np+11*6*(nr+ng)
        echo "  products per site: ",c," action=",actionProducts,
          " gradient=",gradientProducts," Hessian=",work.plan.steps.len

    test "zero terms and vanishing derivative orders allocate no work":
      let w=newLoopWork(g[0])
      let raw=getRawMemAllocated()
      let z=GaugeActionCoeffs()
      check z.loopAction(g,work=w) == 0
      z.loopDeriv(g,f,work=w)
      check f.norm2 == 0
      cases[0].loopDeriv(g,[hs[0],hs[1],hs[2],hs[3]],f,work=w)
      check f.norm2 == 0
      check w.plan == nil and w.halos.len == 0 and w.products.len == 0
      check getRawMemAllocated() == raw

    test "high jets skip plaquettes and the top jet skips gauge values":
      let c=cases[^1]
      var six=c
      six.plaq=0
      c.loopDeriv(g,[hs[0],hs[1],hs[2],hs[3]],f,work=work)
      let p=work.plan
      six.loopDeriv(g,[hs[0],hs[1],hs[2],hs[3]],want,work=work)
      check work.plan == p
      check relativeDiff(f,want) == 0
      c.loopDeriv(g,hs,f,work=work)
      for load in work.plan.loads:
        privateAccess(typeof(load))
        check load.seed > 0

run(3)
run(1)
qexFinalize()
