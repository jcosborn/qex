import qex
import base/alignedMem, gauge/staples
import testutils

proc run(lo: auto; nc: static int) =
  let
    g = lo.newGauge(nc)
    saved = lo.newGauge(nc)
    f = lo.newGauge(nc)
    cached = lo.newGauge(nc)
    expected = lo.newGauge(nc)
    kept = newStapleWork(g[0])
    reference = newLoopWork(g[0])
    nd = lo.nDim
    c = Iwasaki(5.4)
    wilson = GaugeActionCoeffs(plaq: 5.4)
  var rng = newRNGField(RngMilc6, lo, 992311'u64)
  suite "Temporary rectangle batches Nc=" & $nc & " lattice=" & $lo.physGeom:
    test "changed links and incomplete batches match independent loops":
      for version in 1..3:
        threads:
          g.random rng
          for mu in 0..<nd: saved[mu] := g[mu]
        let tag: LoopTag = (811,version)
        let a = c.gaugeAction1(g)
        let loopValue = c.loopAction(g, work=reference)
        let retainedValue = c.gaugeAction1(g, work=kept, tag=tag)
        check abs(a-loopValue) < 1e-11*(1+abs(loopValue))
        check abs(a-retainedValue) < 1e-11*(1+abs(retainedValue))
        c.gaugeActionDeriv(g, f)
        (-1.0*c).loopDeriv(g, expected, work=reference)
        c.gaugeActionDeriv(g, cached, work=kept, tag=tag)
        check relativeDiff(f, expected) < 2e-12
        check relativeDiff(f, cached) < 2e-14
        let builds = kept.builds
        check c.gaugeAction1(g, work=kept, tag=tag) == retainedValue
        c.gaugeActionDeriv(g, expected, work=kept, tag=tag)
        check relativeDiff(expected, cached) == 0
        check kept.builds == builds
        threads:
          for mu in 0..<nd: expected[mu] := 2.0*f[mu]
        c.gaugeActionDeriv(g, f, accumulate=true)
        check relativeDiff(f, expected) < 2e-14
        check relativeDiff(g, saved) == 0
        wilson.gaugeActionDeriv(g, f)
        (-1.0*wilson).loopDeriv(g, expected, work=reference)
        check relativeDiff(f, expected) < 2e-12
    test "an action without a rect term matches independent loops":
      # This path must not start the stu exchanges: nothing waits for them, and later
      # calls reuse the buffers their receives write into.
      let expected = wilson.loopAction(g)
      for repeat in 0..<6:
        let a = stapleBatches(wilson.plaq, 0.0, g, default(seq[typeof(g[0])]), false)
        check abs(a-expected) < 1e-11*(1+abs(expected))
    test "scratch is bounded and released across repeated calls":
      # Numerical work is warmed before measuring ownership and allocation.
      discard c.gaugeAction1(g)
      c.gaugeActionDeriv(g, f)
      wilson.gaugeActionDeriv(g, f)
      GC_fullCollect()
      let live = getRawMemUsed()
      let slots = min(3,nd*(nd-1) div 2)
      let bytes = lo.nSitesOuter*sizeof(typeof(g[0][0]))
      for repeat in 0..<4:
        let beforeAction = getRawMemAllocated()
        discard c.gaugeAction1(g)
        check getRawMemAllocated()-beforeAction <= 4*slots*bytes+g[0].s.align
        check getRawMemUsed() <= live
        let beforeForce = getRawMemAllocated()
        c.gaugeActionDeriv(g, f)
        check getRawMemAllocated()-beforeForce <= 6*slots*bytes+g[0].s.align
        check getRawMemUsed() <= live
        # Without a rect term a pair needs its forward and backward staples only.
        let beforeWilson = getRawMemAllocated()
        wilson.gaugeActionDeriv(g, f)
        check getRawMemAllocated()-beforeWilson <= 4*slots*bytes+g[0].s.align
        check getRawMemUsed() <= live

qexInit()
# Two dimensions need one pair, three need three, and five need ten: the
# latter exercises a short final batch after three complete batches.
# LayoutX supports at most two SIMD lanes along each direction.
const testVectorLength = VLEN
for local in [@[4,4], @[4,4,4], @[4,4,4,4], @[4,4,4,4,4]]:
  if testVectorLength > (1 shl local.len): continue
  let lo = newLayout(latticeFromLocalLattice(local,nRanks))
  run(lo,3)
  run(lo,1)
qexFinalize()
