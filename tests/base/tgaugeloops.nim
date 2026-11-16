import qex
import base/alignedMem
import gauge/symanzik1loopAction
import gauge/staples
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
  discard c.gaugeAction1(a, loop=work)
  c.gaugeActionDeriv(a, f, loop=work)
  c.gaugeDerivDeriv2(a, b, f, work=work)
  c.gaugeDeriv2Subset(a, f[0], newSubsetWork(a[0], 0, 0), work=work)
  c.gaugeDerivDeriv2Subset(a, b[0], f, newSubsetWork(a[0], 0, 0), work=work)
  # Drop caller references even if refc still sees a dead stack slot.
  for mu in 0..<g.len:
    a[mu] = nil
    b[mu] = nil

proc transient(c: GaugeActionCoeffs, g, h, f: auto) {.noinline.} =
  # Clear the owning heap slot even if refc sees a dead stack temporary.
  var work = @[newLoopWork(g[0])]
  fresh(c, g, h, f, work[0])
  work[0] = nil

proc transientStaples(c: GaugeActionCoeffs, g, f, loop: auto) {.noinline.} =
  var work = @[newStapleWork(g[0])]
  discard c.gaugeAction1(g, work=work[0], loop=loop)
  c.gaugeActionDeriv(g, f, work=work[0], loop=loop)
  work[0] = nil

proc transientLoop(c: GaugeActionCoeffs, g, f, staples: auto) {.noinline.} =
  var work = @[newLoopWork(g[0])]
  discard c.loopAction(g, work=work[0])
  discard c.gaugeAction1(g, work=staples, loop=work[0])
  c.gaugeActionDeriv(g, f, work=staples, loop=work[0])
  work[0] = nil

const coeffs = [GaugeActionCoeffs(rect: 0.31),
                GaugeActionCoeffs(pgm: -0.23),
                GaugeActionCoeffs(plaq: 0.8, rect: -0.11, pgm: 0.17)]

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
  suite "Fundamental loops Nc=" & $nc & " lattice=" & $lo.physGeom:
    test "identity normalization and transporter action values":
      for c in coeffs:
        let
          want = -lo.physVol.float * (c.plaq*float(nd*(nd-1) div 2) +
            c.rect*float(nd*(nd-1)) + c.pgm*float(4*nd*(nd-1)*(nd-2) div 6))
          a = c.gaugeAction1(g, loop=work)
        check abs(c.gaugeAction1(id, loop=work)-want) < 1.0e-10*(1.0+abs(want))
        check abs(a-c.gaugeAction2(g)) < 1.0e-11*(1.0+abs(a))
    test "independent loop products and alternate actions":
      for c in coeffs:
        let a = c.gaugeAction1(g, loop=work)
        check abs(a-explicitAction(c, g)) < 1.0e-11*(1.0+abs(a))
        check abs(a-c.gaugeAction3(g, work=work)) < 1.0e-11*(1.0+abs(a))
        if nd == 4:
          var pi: PerfInfo
          # Its shared-staple contractions cancel U*U.adj and assume unitarity.
          let s = c.symanzik1loopAction(unit, pi)
          let off = lo.physVol.float*(6*c.plaq+12*c.rect+16*c.pgm)
          let au = c.gaugeAction1(unit, loop=work)
          check abs(au-(s.space+s.time-off)) < 1.0e-11*(1.0+abs(au))

    test "ambient action gradients, accumulation, and coefficient updates":
      for c in coeffs:
        (-1.0*c).gaugeActionDeriv(g, d, loop=work)
        threads:
          for mu in 0..<nd:
            e[mu] := 0
        c.gaugeDeriv2(g, e, work=work)
        check relativeDiff(d, e) < 2.0e-12
        let
          num = (c.gaugeAction1(gp, loop=work)-c.gaugeAction1(gm, loop=work))/(2*eps)
          ana = pair(d, h)
        check abs(num-ana) < 2.0e-8*(1.0+abs(ana))
        threads:
          for mu in 0..<nd:
            e[mu] := k[mu]
            tmp[mu] := k[mu]+d[mu]
        (-1.0*c).gaugeActionDeriv(g, e, accumulate=true, loop=work)
        check relativeDiff(tmp, e) < 2.0e-12
      for v in [0.0, 0.19, 0.0, -0.07]:
        let c = GaugeActionCoeffs(plaq: 0.8, rect: v, pgm: 2*v)
        let a = c.gaugeAction1(g, loop=work)
        check abs(a-c.gaugeAction2(g)) < 1.0e-11*(1.0+abs(a))

    test "ambient Hessians and symmetry":
      for c in coeffs:
        threads:
          for mu in 0..<nd:
            hv[mu] := 0
            kv[mu] := 0
        c.gaugeDerivDeriv2(g, h, hv, work=work)
        c.gaugeDerivDeriv2(g, k, kv, work=work)
        (-1.0*c).gaugeActionDeriv(gp, hp, loop=work)
        (-1.0*c).gaugeActionDeriv(gm, hm, loop=work)
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
        c.gaugeForce(g, d, loop=work)
        c.gaugeForce2(g, e, work=work)
        check relativeDiff(d, e) < 2.0e-12
        c.gaugeForce3(g, e, work=work)
        check relativeDiff(d, e) < 2.0e-12
        (-1.0*c).gaugeActionDeriv(g, d, loop=work)
        threads:
          for mu in 0..<nd:
            hv[mu] := 0
        c.gaugeDerivDeriv2(g, h, hv, work=work)
        c.gaugeForce(gp, hp, loop=work)
        c.gaugeForce(gm, hm, loop=work)
        threads:
          for mu in 0..<nd:
            tmp[mu] := (0.5/eps)*(hp[mu]-hm[mu])
            for x in lo:
              e[mu][x].projectTAH(-h[mu][x]*d[mu][x].adj-g[mu][x]*hv[mu][x].adj)
        check relativeDiff(tmp, e) < 2.0e-8

    test "subset gradient and Hessian masks, add, base, and sum":
      for c in @[GaugeActionCoeffs(plaq: 0.8)] & @coeffs:
        (-1.0*c).gaugeActionDeriv(g, d, loop=work)
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
            c.gaugeDeriv2Subset(g, hp[dir], newSubsetWork(g[0], parity, dir), work=work)
            check relativeDiff(hp, hm) < 2.0e-12
            let sw = newSubsetWork(g[0], parity, dir)
            threads:
              hp[dir] := k[dir]
              hm[dir] := k[dir]
              threadBarrier()
              for x in sub: hm[dir][x] := d[dir][x]
            c.gaugeDeriv2Subset(g, hp[dir], sw, false, tag=(1,1), work=work)
            check relativeDiff(hp, hm) < 2.0e-12
            # A Hessian at the derivative's tag reads the kept link boundaries
            # and equals one that exchanges them; a later one keeps the tag.
            let shifted = c.rect == 0 and c.pgm == 0
            check sw.tag == (if shifted: (1,1) else: default(LoopTag))
            for t in [(1,1), (1,1), (2,1)]:
              threads:
                for mu in 0..<nd:
                  hp[mu] := k[mu]
                  hm[mu] := k[mu]
              c.gaugeDerivDeriv2Subset(g, h[dir], hp, sw, t, work=work, accumulate=true)
              c.gaugeDerivDeriv2Subset(g, h[dir], hm, newSubsetWork(g[0], parity, dir), work=work, accumulate=true)
              check relativeDiff(hp, hm) == 0.0
              check sw.tag == (if shifted: t else: default(LoopTag))
            c.gaugeDerivDeriv2(g, tmp, hv, work=work)
            c.gaugeDerivDeriv2Subset(g, h[dir], e, newSubsetWork(g[0], parity, dir), work=work)
            check relativeDiff(hv, e) < 2.0e-12
            threads:
              for mu in 0..<nd:
                hp[mu] := k[mu]
                hm[mu] := k[mu]+hv[mu]
            c.gaugeDerivDeriv2Subset(g, h[dir], hp, newSubsetWork(g[0], parity, dir), work=work, accumulate=true)
            check relativeDiff(hp, hm) < 2.0e-12
            threads:
              for mu in 0..<nd:
                hp[mu] := 0
              threadBarrier()
              for x in sub:
                hp[dir][x] := h[dir][x]
                hm[dir][x] := h[dir][x]+hv[dir][x]
            c.gaugeDerivDeriv2SubsetAddBase(g, h[dir], k, hp, newSubsetWork(g[0], parity, dir), work=work)
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
            c.gaugeDerivDeriv2SubsetSum(g, @[h[dir],k[dir],h[dir]], w, hp, newSubsetWork(g[0], parity, dir), work=work)
            check relativeDiff(hp, hv) < 2.0e-12

    test "fundamental and adjoint families remain explicit":
      let c = GaugeActionCoeffs(rect: 0.1, pgm: 0.2, adjplaq: 0.3)
      expect(ValueError): discard c.gaugeAction1(g, loop=work)
      expect(ValueError): c.gaugeActionDeriv(g, d, loop=work)
      expect(ValueError): c.gaugeDerivDeriv2(g, h, d, work=work)
      expect(ValueError): c.gaugeDeriv2Subset(g, d[0], newSubsetWork(g[0], 0, 0), work=work)
      expect(ValueError): c.gaugeDerivDeriv2Subset(g, h[0], d, newSubsetWork(g[0], 0, 0), work=work)
      expect(ValueError): discard c.actionA(g)
      expect(ValueError): c.gaugeADeriv(g, d)

    test "public action and force dispatch by coefficient family":
      for c in [GaugeActionCoeffs(plaq: 0.8), coeffs[0], coeffs[1], coeffs[2],
                GaugeActionCoeffs(plaq: 0.8, adjplaq: 0.13)]:
        let a = c.action(g, work=work)
        # Both calls add the same rank partials; the MPI sum may associate differently.
        check c.action(g, work=nil) ~ a
        if c.adjplaq != 0:
          check abs(a-c.actionA(g)) < 1.0e-12*(1.0+abs(a))
          c.forceA(g, d)
        else:
          check abs(a-c.gaugeAction2(g)) < 1.0e-11*(1.0+abs(a))
          c.gaugeForce(g, d, loop=work)
        c.force(g, e, work=work)
        check relativeDiff(d, e) < 2.0e-12
      let mixed = GaugeActionCoeffs(plaq: 0.8, rect: 0.1, adjplaq: 0.2)
      expect(ValueError): discard mixed.action(g, work=work)
      expect(ValueError): mixed.force(g, e, work=work)

    test "independent workspaces refresh inputs and retain their own storage":
      let other = newLoopWork(g[0])
      proc buffers(w: auto): seq[pointer] =
        for row in w.halos:
          for slot in row.slots:
            for h in slot.fields:
              doAssert h.field == nil
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
        c.gaugeDerivDeriv2Subset(g, h[0], hp, newSubsetWork(g[0], 0, 0), work=work)
        let raw = getRawMemAllocated()
        c.gaugeDerivDeriv2Subset(gp, k[1], hp, newSubsetWork(gp[0], 1, 1), work=work)
        check getRawMemAllocated() == raw
        # The subset gradient should allocate only the existing action scratch.
        c.gaugeDeriv2Subset(g, hp[0], newSubsetWork(g[0], 0, 0), work=work)
        let before = getRawMemAllocated()
        (-1.0*c).gaugeActionDeriv(g, e, loop=work)
        let full = getRawMemAllocated()-before
        let after = getRawMemAllocated()
        c.gaugeDeriv2Subset(g, hp[1], newSubsetWork(g[0], 1, 1), work=work)
        check getRawMemAllocated()-after == full

    test "kept staples follow coefficient changes at the same tag":
      let plaq = GaugeActionCoeffs(plaq: 0.8)
      let rect = GaugeActionCoeffs(plaq: 0.8, rect: -0.3)
      let cache = newStapleWork(g[0])
      let tag: LoopTag = (17, 0)
      # Tagged plaquette actions also prepare the backward staples rectangles need.
      discard plaq.gaugeAction1(g, work=cache, tag=tag)
      var builds = cache.builds
      let want = rect.gaugeAction1(g)
      check abs(rect.gaugeAction1(g, work=cache, tag=tag)-want) < 1e-12*(1+abs(want))
      check cache.builds == builds
      builds = cache.builds
      discard plaq.gaugeAction1(g, work=cache, tag=tag)
      rect.gaugeActionDeriv(g, d, work=cache, tag=tag)
      rect.gaugeActionDeriv(g, e)
      check cache.builds == builds
      check relativeDiff(d, e) < 2e-14
      # An untagged plaquette action refreshes only the forward staples.
      discard plaq.gaugeAction1(g, work=cache)
      builds = cache.builds
      check abs(rect.gaugeAction1(g, work=cache, tag=tag)-want) < 1e-12*(1+abs(want))
      check cache.builds == builds+1

    test "tagged calls reuse kept staples at unchanged links":
      let c = GaugeActionCoeffs(plaq: 0.8, rect: -0.3)
      let cache = newStapleWork(g[0])
      var tag: LoopTag = (7, 1)
      for changed in [g, gp, gm, g]:
        inc tag.ver
        let native = c.gaugeAction1(changed)
        c.gaugeActionDeriv(changed, e)
        let a = c.gaugeAction1(changed, work=cache, tag=tag)
        # Streaming changes the order of direction-pair sums. Cache reuse
        # must still reproduce the cached value exactly.
        check abs(a-native) < 1e-12*(1+abs(native))
        let builds = cache.builds
        c.gaugeActionDeriv(changed, d, work=cache, tag=tag)  # reuses the action's build
        check cache.builds == builds
        check relativeDiff(d, e) < 2e-14
        check c.gaugeAction1(changed, work=cache, tag=tag) == a
        check cache.builds == builds
        threads:
          for mu in 0..<nd: tmp[mu] := 2.0*d[mu]
        c.gaugeActionDeriv(changed, d, accumulate=true, work=cache, tag=tag)
        check relativeDiff(d, tmp) < 2e-15
        # A tag of other links builds; an untagged call builds and clears the tag.
        check cache.reuses == 3
        cache.reuses = 0
      c.gaugeActionDeriv(g, d, work=cache)
      check cache.reuses == 0
      c.gaugeActionDeriv(g, e, work=cache, tag=tag)
      check cache.reuses == 0
      # Without an owner, tags cannot reuse data from earlier links.
      discard c.gaugeAction1(g, tag=tag)
      c.gaugeActionDeriv(gp, d, tag=tag)
      c.gaugeActionDeriv(gp, e)
      check relativeDiff(d, e) < 2e-14
      # Every tagged plaquette action prepares both sets for immediate force reuse.
      # Tagged builds fill both sets together, unretained derivatives sum the
      # boundary terms per batch: compare within rounding.
      let w = GaugeActionCoeffs(plaq: 0.8)
      let a = w.gaugeAction1(g)
      let fresh = newStapleWork(g[0])
      inc tag.ver
      check abs(w.gaugeAction1(g, work=fresh, tag=tag)-a) < 1e-14*(1+abs(a))
      var builds = fresh.builds
      w.gaugeActionDeriv(g, d, work=fresh, tag=tag)
      check fresh.builds == builds
      w.gaugeActionDeriv(g, e)
      check relativeDiff(d, e) < 2e-15
      inc tag.ver
      check abs(w.gaugeAction1(g, work=fresh, tag=tag)-a) < 1e-14*(1+abs(a))
      builds = fresh.builds
      w.gaugeActionDeriv(g, d, work=fresh, tag=tag)
      check fresh.builds == builds
      check relativeDiff(d, e) < 2e-15
      check w.gaugeAction1(g, work=fresh) == a  # untagged: stf only, no reuse
      builds = fresh.builds
      w.gaugeActionDeriv(g, d, work=fresh, tag=tag)
      check fresh.builds == builds+1
      check relativeDiff(d, e) < 2e-15

    test "staple and loop tags remain independent across interleaved calls":
      for c in [GaugeActionCoeffs(plaq: 0.8),
                GaugeActionCoeffs(plaq: 0.8, rect: -0.11, pgm: 0.07)]:
        let cache = newStapleWork(g[0])
        let lw = newLoopWork(g[0])
        for pass, src in [g, gp, gm, g]:
          let tag: LoopTag = (101, pass+1)
          let alt: LoopTag = (102, pass+1)
          let next = if pass mod 2 == 0: gp else: gm
          let a = c.gaugeAction1(src)
          # The loop call must not mark staples from the previous links as current.
          discard c.loopAction(src, work=lw, tag=tag)
          check abs(c.gaugeAction1(src, work=cache, loop=lw, tag=tag)-a) < 1e-11*(1+abs(a))
          c.gaugeActionDeriv(src, e)
          c.gaugeActionDeriv(src, d, work=cache, loop=lw, tag=tag)
          check relativeDiff(d, e) < 2e-12
          let builds = cache.builds
          c.loopDeriv(next, d, work=lw, tag=alt)
          check abs(c.gaugeAction1(src, work=cache, loop=lw, tag=tag)-a) < 1e-11*(1+abs(a))
          check cache.builds == builds
          # Force and kick use the loop cache too, while retained staples change independently.
          c.force(next, d, work=lw, tag=alt)
          c.gaugeActionDeriv(next, e)
          c.gaugeActionDeriv(next, d, work=cache, loop=lw, tag=alt)
          check relativeDiff(d, e) < 2e-12
          discard c.forceKick(src, h, tmp, 0.1, work=lw, tag=tag)
          c.gaugeActionDeriv(src, e)
          c.gaugeActionDeriv(src, d, work=cache, loop=lw, tag=tag)
          check relativeDiff(d, e) < 2e-12
          c.gaugeActionDeriv(src, d, work=cache, loop=lw)
          discard c.loopAction(src, work=lw)
          c.gaugeActionDeriv(src, d, work=cache, loop=lw, tag=tag)
          check relativeDiff(d, e) < 2e-12

    test "tagged loop calls reuse the kept halo at unchanged links":
      let c = GaugeActionCoeffs(plaq: 0.8, rect: -0.3, pgm: 0.1)
      let cache = newLoopWork(g[0])
      var tag: LoopTag = (9, 1)
      # Every tagged action fills the derivative halo; its force reuses it,
      # including the first pair through a fresh workspace.
      for changed in [g, gp, gm, g]:
        inc tag.ver
        let a = c.loopAction(changed)
        c.loopDeriv(changed, e)
        let updates = cache.updates
        check c.loopAction(changed, work=cache, tag=tag) == a
        check cache.updates == updates+nd
        c.loopDeriv(changed, d, work=cache, tag=tag)
        check cache.updates == updates+nd
        check relativeDiff(d, e) == 0
        c.loopDeriv(changed, d, work=cache, tag=tag)  # same links and halo: no exchange
        check cache.updates == updates+nd
        check relativeDiff(d, e) == 0
      # Tags belong to a halo row and slot: an untagged action exchanges its own
      # row and the tagged derivative still reuses; an untagged derivative clears
      # the derivative row's tag.
      var reuses = cache.reuses
      discard c.loopAction(g, work=cache)
      c.loopDeriv(g, d, work=cache, tag=tag)
      check cache.reuses == reuses+1
      c.loopDeriv(g, e)
      check relativeDiff(d, e) == 0
      c.loopDeriv(g, d, work=cache)
      reuses = cache.reuses
      c.loopDeriv(g, d, work=cache, tag=tag)
      check cache.reuses == reuses
      check relativeDiff(d, e) == 0

    test "halo reuse distinguishes plans and always refreshes seeds":
      let cache = newLoopWork(g[0])
      let tag: LoopTag = (19, 0)
      # These plans have the same rectangle extents but different link coverage.
      for c in [GaugeActionCoeffs(rect: -0.3),
                GaugeActionCoeffs(plaq: 0.8, rect: -0.3),
                GaugeActionCoeffs(rect: -0.3)]:
        var updates = cache.updates
        c.loopDeriv(g, d, work=cache, tag=tag)
        check cache.updates == updates+nd
        c.loopDeriv(g, e)
        check relativeDiff(d, e) == 0
        updates = cache.updates
        c.loopDeriv(g, d, work=cache, tag=tag)
        check cache.updates == updates
        for seed in [h, k, h]:
          c.loopDeriv(g, [seed], d, work=cache, tag=tag)
          check cache.updates == updates+nd
          updates = cache.updates
          c.loopDeriv(g, [seed], e)
          check relativeDiff(d, e) == 0

    test "retained staple storage survives collection between calls":
      let c = GaugeActionCoeffs(plaq: 0.8)
      let cache = newStapleWork(g[0])
      let threshold = getRawMemGcThreshold()
      setRawMemGcThreshold(0)
      defer: setRawMemGcThreshold(threshold)
      for changed in [g, gp, gm]:
        let a = c.gaugeAction1(changed)
        GC_fullCollect()
        check c.gaugeAction1(changed, work=cache) ~ a
        c.gaugeActionDeriv(changed, d, work=cache)
        GC_fullCollect()
        c.gaugeActionDeriv(changed, e)
        check relativeDiff(d, e) < 2e-15

    test "Wilson calls share staple storage and match native results":
      let c = GaugeActionCoeffs(plaq: 0.8)
      let cache = newStapleWork(g[0])
      let lw = newLoopWork(g[0])
      let saved = g.newOneOf
      let seed = h.newOneOf
      let blockBytes = nd*(nd-1)*lo.nSitesOuter*sizeof(typeof(g[0][0]))+g[0].s.align
      threads:
        for mu in 0..<nd:
          saved[mu] := g[mu]
          seed[mu] := h[mu]
      check cache.bytes == 0
      check c.gaugeAction1(g, work=cache) ~ c.gaugeAction1(g)
      check cache.bytes == blockBytes
      let before = getRawMemAllocated()
      c.gaugeActionDeriv(g, d, work=cache)
      check getRawMemAllocated()-before == blockBytes
      check cache.bytes == 2*blockBytes
      c.gaugeActionDeriv(g, e)
      check relativeDiff(d, e) < 2e-15
      threads:
        for mu in 0..<nd:
          hp[mu] := k[mu]
          hm[mu] := k[mu]
      c.gaugeDerivDeriv2(g, h, hp, work=lw)
      c.gaugeDerivDeriv2(g, h, hm)
      check relativeDiff(hp, hm) == 0
      let owned = 4*nd*g[0].s.mem.size
      check lw.bytes == owned
      check cache.bytes == 2*blockBytes
      for changed in [g, gp, gm]:
        let a = c.gaugeAction1(changed)
        c.gaugeActionDeriv(changed, e)
        threads:
          for mu in 0..<nd:
            hp[mu] := k[mu]
            hm[mu] := k[mu]
        c.gaugeDerivDeriv2(changed, k, hm)
        let raw = getRawMemAllocated()
        check c.gaugeAction1(changed, work=cache) ~ a
        c.gaugeActionDeriv(changed, d, work=cache)
        c.gaugeDerivDeriv2(changed, k, hp, work=lw)
        check c.gaugeAction1(changed, work=cache) ~ a
        check getRawMemAllocated() == raw
        check lw.bytes == owned
        check cache.bytes == 2*blockBytes
        check relativeDiff(d, e) < 2e-15
        check relativeDiff(hp, hm) == 0
        threads:
          for mu in 0..<nd:
            d[mu] := k[mu]
            tmp[mu] := k[mu]+e[mu]
        c.gaugeActionDeriv(changed, d, accumulate=true, work=cache)
        check relativeDiff(d, tmp) < 2e-15
      check relativeDiff(g, saved) == 0
      check relativeDiff(h, seed) == 0
      expect(ValueError):
        discard GaugeActionCoeffs(adjplaq: 1.0).gaugeAction1(g, work=cache)
      check c.gaugeAction1(g, work=cache) ~ c.gaugeAction1(g)
      let wrong = newStapleWork(newLayout(lo.physGeom).newGauge(nc)[0])
      expect(ValueError): discard c.gaugeAction1(g, work=wrong)
      expect(ValueError): c.gaugeActionDeriv(g, d, work=wrong)
      check wrong.bytes == 0
      let wrongLoop = newLoopWork(newLayout(lo.physGeom).newGauge(nc)[0])
      let builds = cache.builds
      expect(ValueError): discard coeffs[2].gaugeAction1(g, work=cache, loop=wrongLoop)
      expect(ValueError): coeffs[2].gaugeActionDeriv(g, d, work=cache, loop=wrongLoop)
      expect(ValueError): coeffs[2].gaugeForce(g, d, work=cache, loop=wrongLoop)
      check cache.builds == builds
      check wrongLoop.bytes == 0

    test "rectangle forces reuse the rect staple block":
      let c = GaugeActionCoeffs(plaq: 0.8, rect: -0.3)
      let cache = newStapleWork(g[0])
      let blockBytes = nd*(nd-1)*lo.nSitesOuter*sizeof(typeof(g[0][0]))+g[0].s.align
      c.gaugeActionDeriv(g, d, work=cache)
      check cache.bytes == 3*blockBytes
      c.gaugeActionDeriv(g, e)
      check relativeDiff(d, e) < 2e-14
      for changed in [gp, gm, g]:
        c.gaugeActionDeriv(changed, e)
        let raw = getRawMemAllocated()
        c.gaugeActionDeriv(changed, d, work=cache)
        check getRawMemAllocated() == raw
        check cache.bytes == 3*blockBytes
        check relativeDiff(d, e) < 2e-14

    test "workspace bytes follow owned halo product and scratch growth":
      let cache = newLoopWork(g[0])
      let c = GaugeActionCoeffs(pgm: 0.17)
      var raw = getRawMemAllocated()
      c.loopDeriv(g, [h], d, work=cache)
      check cache.bytes == getRawMemAllocated()-raw
      var owned = cache.bytes
      raw = getRawMemAllocated()
      discard c.loopAction(g, work=cache)
      check cache.bytes-owned == getRawMemAllocated()-raw
      # The higher-order product bank is retained for subsequent smaller plans.
      owned = cache.bytes
      raw = getRawMemAllocated()
      c.gaugeDeriv2Subset(g, d[0], newSubsetWork(g[0], 0, 0), work=cache)
      check cache.bytes-owned == getRawMemAllocated()-raw
      owned = cache.bytes
      raw = getRawMemAllocated()
      c.loopDeriv(gp, [k], d, work=cache)
      c.gaugeDeriv2Subset(gm, d[1], newSubsetWork(gm[0], 1, 1), work=cache)
      check cache.bytes == owned
      check getRawMemAllocated() == raw
      let st = newStapleWork(g[0])
      check c.gaugeAction1(g, work=st, loop=cache) ~ c.gaugeAction2(g)
      c.gaugeActionDeriv(g, d, work=st, loop=cache)
      check st.bytes == 0
      check st.builds == 0

    test "public staple reuse rejects overlapping or mismatched storage":
      let cs = startCornerShifts(g)
      let st = makeFwdStaples(g, cs)
      expect(ValueError): discard makeStaples(g, cs, st, st)
      var links = newSeq[typeof(g[0])](nd)
      for mu in 0..<nd: links[mu] = st[0,1]
      expect(ValueError): discard makeFwdStaples(links, cs, st)
      let other = newLayout(lo.physGeom).newGauge(nc)
      expect(ValueError): discard makeFwdStaples(other, cs, st)
      let wrong = newFieldArray(lo, typeof(g[0]), nd*(nd-1))
      expect(ValueError): discard makeFwdStaples(g, cs, wrong)


proc releaseLoops(lo: auto, nc: static int) =
  # refc scans live frames conservatively.  In runLoops, a stale word left by
  # earlier tests pointed into a freed cell that `transient` reused for a
  # field, which then survived GC_fullCollect: check from a small frame.
  let
    g = lo.newGauge(nc)
    h = lo.newGauge(nc)
    e = lo.newGauge(nc)
    work = newLoopWork(g[0])
  var r = newRNGField(RngMilc6, lo, 738142'u64)
  threads:
    g.random r
    h.gaussian r
  suite "Loop workspace release Nc=" & $nc & " lattice=" & $lo.physGeom:
    test "staple and loop owners release independently":
      let c = coeffs[2]
      let st = newStapleWork(g[0])
      discard c.loopAction(g, work=work)
      transientStaples(c, g, e, work)
      transientLoop(c, g, e, st)
      GC_fullCollect()
      GC_fullCollect()
      let before = getRawMemUsed()
      let sb = st.bytes
      let lb = work.bytes
      check sb > 0 and lb > 0
      transientStaples(c, g, e, work)
      transientLoop(c, g, e, st)
      GC_fullCollect()
      GC_fullCollect()
      check getRawMemUsed() <= before
      check st.bytes == sb
      check work.bytes == lb

    test "loop workspaces release temporary input fields":
      # A collection inside a call may also free unrelated garbage.
      for c in coeffs:
        discard c.gaugeAction1(g, loop=work)
        c.gaugeActionDeriv(g, e, loop=work)
        c.gaugeDerivDeriv2(g, h, e, work=work)
        c.gaugeDeriv2Subset(g, e[0], newSubsetWork(g[0], 0, 0), work=work)
        c.gaugeDerivDeriv2Subset(g, h[0], e, newSubsetWork(g[0], 0, 0), work=work)
        GC_fullCollect()
        GC_fullCollect()
        let before = getRawMemUsed()
        fresh(c, g, h, e, work)
        GC_fullCollect()
        GC_fullCollect()
        check getRawMemUsed() <= before

    #[ # fails with nvcc
    test "temporary loop workspaces do not accumulate storage":
      # refc may keep a dead workspace through conservative stack roots.
      # Warm that call path before checking repeated allocations for growth.
      for c in coeffs:
        var before = 0
        for rep in 0..<4:
          transient(c, g, h, e)
          GC_fullCollect()
          GC_fullCollect()
          let used = getRawMemUsed()
          if rep == 0: before = used
          else: check used <= before
    ]#

qexInit()
letParam:
  lat = latticeFromLocalLattice(@[4,4,4,4], nRanks)
let lo = newLayout(lat)
runLoops(lo, 3)
releaseLoops(lo, 3)
runLoops(lo, 1)
releaseLoops(lo, 1)
qexFinalize()
