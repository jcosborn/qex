import qex
import base/alignedMem
import gauge/symanzik1loopAction
import testutils

proc pair(a, b: auto): float =
  for mu in 0..<a.len:
    result += redot(a[mu], b[mu])

proc explicitAction(c: GaugeActionCoeffs, g: auto): float =
  let nd = g.len
  var paths: seq[seq[int]]
  var coeffs: seq[float]
  for mu in 2..nd:
    for nu in 1..<mu:
      paths.add @[mu,nu,-mu,-nu]
      coeffs.add c.plaq
      paths.add [@[mu,nu,nu,-mu,-nu,-nu], @[mu,mu,nu,-mu,-mu,-nu]]
      coeffs.add [c.rect, c.rect]
      for sg in 1..<nu:
        paths.add [@[mu,nu,sg,-mu,-nu,-sg], @[mu,sg,nu,-mu,-sg,-nu],
                   @[nu,mu,sg,-nu,-mu,-sg], @[mu,-nu,sg,-mu,nu,-sg]]
        coeffs.add [c.pgm, c.pgm, c.pgm, c.pgm]
  for i in 0..<paths.len:
    if coeffs[i] != 0:
      result -= g[0].l.physVol.float * coeffs[i] * g.wline(paths[i]).re

proc fresh(c: GaugeActionCoeffs, g, h, f, work: auto) {.noinline.} =
  var a = newOneOf(g)
  var b = newOneOf(h)
  threads:
    for mu in 0..<g.len:
      a[mu] := g[mu]
      b[mu] := h[mu]
  discard c.gaugeAction1(a, work=work)
  c.gaugeActionDeriv(a, f, work=work)
  c.gaugeDerivDeriv2(a, b, f, work=work)
  c.gaugeDeriv2Subset(a, f, 0, 0, work=work)
  c.gaugeDerivDeriv2Subset(a, b, f, 0, 0, work=work)
  # Drop caller references even if refc still sees a dead stack slot.
  for mu in 0..<g.len:
    a[mu] = nil
    b[mu] = nil

proc transient(c: GaugeActionCoeffs, g, h, f: auto) {.noinline.} =
  var work = newLoopWork(g[0])
  fresh(c, g, h, f, work)
  work = nil

proc runLoops(lo: auto, nc: static int) =
  let
    g = lo.newGauge(nc)
    h = lo.newGauge(nc)
    k = lo.newGauge(nc)
    id = lo.newGauge(nc)
    unit = lo.newGauge(nc)
    gp = lo.newGauge(nc)
    gm = lo.newGauge(nc)
    d = lo.newGauge(nc)
    e = lo.newGauge(nc)
    hp = lo.newGauge(nc)
    hm = lo.newGauge(nc)
    hv = lo.newGauge(nc)
    kv = lo.newGauge(nc)
    tmp = lo.newGauge(nc)
    nd = g.len
    eps = 1.0e-5
    work = newLoopWork(g[0])
  var r = newRNGField(RngMilc6, lo, 738142'u64)
  threads:
    g.random r
    h.gaussian r
    k.gaussian r
    for mu in 0..<nd:
      unit[mu] := g[mu]
      g[mu] += 0.03*h[mu]
      gp[mu] := g[mu]+eps*h[mu]
      gm[mu] := g[mu]-eps*h[mu]
      id[mu] := 1
  let coeffs = [GaugeActionCoeffs(rect: 0.31),
                GaugeActionCoeffs(pgm: -0.23),
                GaugeActionCoeffs(plaq: 0.8, rect: -0.11, pgm: 0.17)]
  suite "Fundamental loops Nc=" & $nc & " lattice=" & $lo.physGeom:
    test "identity normalization and transporter action values":
      for c in coeffs:
        let
          want = -lo.physVol.float * (c.plaq*float(nd*(nd-1) div 2) +
            c.rect*float(nd*(nd-1)) + c.pgm*float(4*nd*(nd-1)*(nd-2) div 6))
          a = c.gaugeAction1(g, work=work)
        check abs(c.gaugeAction1(id, work=work)-want) < 1.0e-10*(1.0+abs(want))
        check abs(a-c.gaugeAction2(g)) < 1.0e-11*(1.0+abs(a))
    test "independent loop products and alternate actions":
      for c in coeffs:
        let a = c.gaugeAction1(g, work=work)
        check abs(a-explicitAction(c, g)) < 1.0e-11*(1.0+abs(a))
        check abs(a-c.gaugeAction3(g, work=work)) < 1.0e-11*(1.0+abs(a))
        if nd == 4:
          var pi: PerfInfo
          # Its shared-staple contractions cancel U*U.adj and assume unitarity.
          let s = c.symanzik1loopAction(unit, pi)
          let off = lo.physVol.float*(6*c.plaq+12*c.rect+16*c.pgm)
          let au = c.gaugeAction1(unit, work=work)
          check abs(au-(s.space+s.time-off)) < 1.0e-11*(1.0+abs(au))

    test "ambient action gradients, accumulation, and coefficient updates":
      for c in coeffs:
        (-1.0*c).gaugeActionDeriv(g, d, work=work)
        threads:
          for mu in 0..<nd:
            e[mu] := 0
        c.gaugeDeriv2(g, e, work=work)
        check relativeDiff(d, e) < 2.0e-12
        let
          num = (c.gaugeAction1(gp, work=work)-c.gaugeAction1(gm, work=work))/(2*eps)
          ana = pair(d, h)
        check abs(num-ana) < 2.0e-8*(1.0+abs(ana))
        threads:
          for mu in 0..<nd:
            e[mu] := k[mu]
            tmp[mu] := k[mu]+d[mu]
        (-1.0*c).gaugeActionDeriv(g, e, accumulate=true, work=work)
        check relativeDiff(tmp, e) < 2.0e-12
      for v in [0.0, 0.19, 0.0, -0.07]:
        let c = GaugeActionCoeffs(plaq: 0.8, rect: v, pgm: 2*v)
        let a = c.gaugeAction1(g, work=work)
        check abs(a-c.gaugeAction2(g)) < 1.0e-11*(1.0+abs(a))

    test "ambient Hessians and symmetry":
      for c in coeffs:
        threads:
          for mu in 0..<nd:
            hv[mu] := 0
            kv[mu] := 0
        c.gaugeDerivDeriv2(g, h, hv, work=work)
        c.gaugeDerivDeriv2(g, k, kv, work=work)
        (-1.0*c).gaugeActionDeriv(gp, hp, work=work)
        (-1.0*c).gaugeActionDeriv(gm, hm, work=work)
        threads:
          for mu in 0..<nd:
            tmp[mu] := (0.5/eps)*(hp[mu]-hm[mu])
        check relativeDiff(hv, tmp) < 2.0e-8
        let a = pair(k, hv)
        check abs(a-pair(h, kv)) < 2.0e-11*(1.0+abs(a))
        threads:
          for mu in 0..<nd:
            e[mu] := k[mu]
            tmp[mu] := k[mu]+hv[mu]
        c.gaugeDerivDeriv2(g, h, e, work=work)
        check relativeDiff(tmp, e) < 2.0e-12

    test "projected force and its full link derivative":
      for c in coeffs:
        c.gaugeForce(g, d, work=work)
        c.gaugeForce2(g, e, work=work)
        check relativeDiff(d, e) < 2.0e-12
        c.gaugeForce3(g, e, work=work)
        check relativeDiff(d, e) < 2.0e-12
        (-1.0*c).gaugeActionDeriv(g, d, work=work)
        threads:
          for mu in 0..<nd:
            hv[mu] := 0
        c.gaugeDerivDeriv2(g, h, hv, work=work)
        c.gaugeForce(gp, hp, work=work)
        c.gaugeForce(gm, hm, work=work)
        threads:
          for mu in 0..<nd:
            tmp[mu] := (0.5/eps)*(hp[mu]-hm[mu])
            for x in lo:
              e[mu][x].projectTAH(-h[mu][x]*d[mu][x].adj-g[mu][x]*hv[mu][x].adj)
        check relativeDiff(tmp, e) < 2.0e-8

    test "subset gradient and Hessian masks, add, base, and sum":
      for c in coeffs:
        (-1.0*c).gaugeActionDeriv(g, d, work=work)
        for parity in 0..1:
          let
            ps = if parity == 0: "even" else: "odd"
            sub = lo.getSubset(ps)
          for dir in 0..<nd:
            threads:
              for mu in 0..<nd:
                hp[mu] := k[mu]
                hm[mu] := k[mu]
                tmp[mu] := 0
                hv[mu] := 0
              threadBarrier()
              hm[dir] := 0
              threadBarrier()
              for x in sub:
                hm[dir][x] := d[dir][x]
                tmp[dir][x] := h[dir][x]
            c.gaugeDeriv2Subset(g, hp, parity, dir, work=work)
            check relativeDiff(hp, hm) < 2.0e-12
            let sd = newShifter(g[0], dir, 1)
            let sf = createShiftBufs(g[0], 1, ps)
            let sb = createShiftBufs(g[0], -1, ps)
            threads:
              hp[dir] := k[dir]
              hm[dir] := k[dir]
              threadBarrier()
              for x in sub: hm[dir][x] := d[dir][x]
            c.gaugeDeriv2SubsetWork(g, hp, sd, sf, sb, parity, dir, false, work=work)
            check relativeDiff(hp, hm) < 2.0e-12
            c.gaugeDerivDeriv2(g, tmp, hv, work=work)
            c.gaugeDerivDeriv2Subset(g, h, e, parity, dir, work=work)
            check relativeDiff(hv, e) < 2.0e-12
            threads:
              for mu in 0..<nd:
                hp[mu] := k[mu]
                hm[mu] := k[mu]+hv[mu]
            c.gaugeDerivDeriv2SubsetAdd(g, h[dir], hp, parity, dir, work=work)
            check relativeDiff(hp, hm) < 2.0e-12
            threads:
              for mu in 0..<nd:
                hp[mu] := 0
              threadBarrier()
              for x in sub:
                hp[dir][x] := h[dir][x]
                hm[dir][x] := h[dir][x]+hv[dir][x]
            c.gaugeDerivDeriv2SubsetAddBase(g, h[dir], k, hp, parity, dir, work=work)
            check relativeDiff(hp, hm) < 2.0e-12
            let w = g[dir].newOneOf
            threads:
              for mu in 0..<nd:
                tmp[mu] := 0
                hv[mu] := 0
              threadBarrier()
              for x in sub:
                tmp[dir][x] := h[dir][x]+k[dir][x]+h[dir][x]
            c.gaugeDerivDeriv2(g, tmp, hv, work=work)
            c.gaugeDerivDeriv2SubsetSum(g, @[h[dir],k[dir],h[dir]], w, hp, parity, dir, work=work)
            check relativeDiff(hp, hv) < 2.0e-12

    test "fundamental and adjoint families remain explicit":
      let c = GaugeActionCoeffs(rect: 0.1, pgm: 0.2, adjplaq: 0.3)
      expect(ValueError): discard c.gaugeAction1(g, work=work)
      expect(ValueError): c.gaugeActionDeriv(g, d, work=work)
      expect(ValueError): c.gaugeDerivDeriv2(g, h, d, work=work)
      expect(ValueError): c.gaugeDeriv2Subset(g, d, 0, 0, work=work)
      expect(ValueError): c.gaugeDerivDeriv2Subset(g, h, d, 0, 0, work=work)
      expect(ValueError): discard c.actionA(g)
      expect(ValueError): c.gaugeADeriv(g, d)

    test "public action and force dispatch by coefficient family":
      for c in [GaugeActionCoeffs(plaq: 0.8), coeffs[0], coeffs[1], coeffs[2],
                GaugeActionCoeffs(plaq: 0.8, adjplaq: 0.13)]:
        let a = c.action(g, work=work)
        check c.action(g, work=nil) == a
        if c.adjplaq != 0:
          check abs(a-c.actionA(g)) < 1.0e-12*(1.0+abs(a))
          c.forceA(g, d)
        else:
          check abs(a-c.gaugeAction2(g)) < 1.0e-11*(1.0+abs(a))
          c.gaugeForce(g, d, work=work)
        c.force(g, e, work=work)
        check relativeDiff(d, e) < 2.0e-12
      let mixed = GaugeActionCoeffs(plaq: 0.8, rect: 0.1, adjplaq: 0.2)
      expect(ValueError): discard mixed.action(g, work=work)
      expect(ValueError): mixed.force(g, e, work=work)

    test "independent workspaces refresh inputs and retain their own storage":
      let other = newLoopWork(g[0])
      proc buffers(w: auto): seq[pointer] =
        for row in w.halos:
          for slot in row.fields:
            for h in slot:
              check h.field == nil
              result.add cast[pointer](h.halo.data)
      for c in coeffs:
        c.gaugeDerivDeriv2(g, h, e, work=other)
        discard c.action(g, work=other)
      let before = buffers(work)
      let separate = buffers(other)
      check before.len > 0
      for p in before: check p notin separate
      for c in coeffs:
        discard c.action(gp, work=work)
        discard c.action(gm, work=other)
        c.force(gm, hp, work=work)
        c.force(gm, hm)
        check relativeDiff(hp, hm) < 2.0e-12
        threads:
          for mu in 0..<nd:
            hp[mu] := 0
            hm[mu] := 0
        c.gaugeDerivDeriv2(gp, k, hp, work=work)
        c.gaugeDerivDeriv2(gp, k, hm, work=other)
        check relativeDiff(hp, hm) < 2.0e-12
      check buffers(work) == before
      check buffers(other) == separate
      let wrong = newLoopWork(newLayout(lo.physGeom).newGauge(nc)[0])
      expect(ValueError): coeffs[0].gaugeDerivDeriv2(g, h, e, work=wrong)

    test "warm workspaces reuse Hessian and subset scratch":
      for c in [GaugeActionCoeffs(plaq: 0.8), coeffs[0], coeffs[1], coeffs[2]]:
        c.gaugeDerivDeriv2(g, h, hp, work=work)
        let raw = getRawMemAllocated()
        c.gaugeDerivDeriv2(gp, k, hp, work=work)
        check getRawMemAllocated() == raw
      for c in coeffs:
        c.gaugeDerivDeriv2Subset(g, h, hp, 0, 0, work=work)
        let raw = getRawMemAllocated()
        c.gaugeDerivDeriv2Subset(gp, k, hp, 1, 1, work=work)
        check getRawMemAllocated() == raw
        # The subset gradient should allocate only the existing action scratch.
        c.gaugeDeriv2Subset(g, hp, 0, 0, work=work)
        let before = getRawMemAllocated()
        (-1.0*c).gaugeActionDeriv(g, e, work=work)
        let full = getRawMemAllocated()-before
        let after = getRawMemAllocated()
        c.gaugeDeriv2Subset(g, hp, 1, 1, work=work)
        check getRawMemAllocated()-after == full

    test "loop workspaces release temporary input fields":
      for c in coeffs:
        discard c.gaugeAction1(g, work=work)
        c.gaugeActionDeriv(g, e, work=work)
        c.gaugeDerivDeriv2(g, h, e, work=work)
        GC_fullCollect()
        GC_fullCollect()
        let before = getRawMemUsed()
        fresh(c, g, h, e, work)
        GC_fullCollect()
        GC_fullCollect()
        check getRawMemUsed() == before
        transient(c, g, h, e)
        GC_fullCollect()
        GC_fullCollect()
        check getRawMemUsed() == before

qexInit()
letParam:
  lat = latticeFromLocalLattice(@[4,4,4,4], nRanks)
let lo = newLayout(lat)
runLoops(lo, 3)
runLoops(lo, 1)
qexFinalize()
