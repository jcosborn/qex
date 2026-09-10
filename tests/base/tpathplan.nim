import testutils
import qex

qexInit()

proc fold(k:seq[int], a:bool):seq[int] =
  # flattened path of node key k with adjoint flag a
  if not a: return k
  result.newseq k.len
  for i in 0..<k.len:
    result[i] = -k[k.len-1-i]

proc shifted(f:auto, sh:seq[int]):auto =
  # f(x - sh) by unit shifts, as the plan defines sh
  var r = newOneOf f
  threads:
    r := f
  for mu,n in sh.pairs:
    if n != 0:
      var s = newShifter(f, mu, if n>0: -1 else: 1)
      for k in 1..abs(n):
        threads:
          r := s ^* r
  return r

suite "Test path plan":
  echo "rank ", myRank, "/", nRanks
  threads: echo "thread ",threadNum," / ",numThreads
  var (lo,g,_) = setupLattice([8,8,8,8])
  let
    plq = @[1,2,-1,-2]
    rct = @[1,1,2,-1,-1,-2]
    rot = [@[-2,1,2,-1], @[-2,1,1,2,-1,-1]]  # plq and rct started from the -2 link
    paths = @[plq, @[2,1,-2,-1], rot[0], rct, rot[1], @[1,2,-1,-2,1,-2,-1,2]]
  const nc = g[0][0].ncols
  let fac = 1.0/float(lo.physVol*nc)
  # staples from the origin, a = U1(x) U2(x+1) U1(x+2)^†, b = U1(x) U1(x+1) U2(x+1+1) U1(x+1+2)^† U1(x+2)^†
  var
    s1 = newShifter(g[0], 0, 1)  # y(x+1)
    s2 = newShifter(g[0], 1, 1)  # y(x+2)
    u2a = newOneOf g[0]  # U2(x+1)
    u1b = newOneOf g[0]  # U1(x+2)
    u1a = newOneOf g[0]  # U1(x+1)
    u2b = newOneOf g[0]  # U2(x+1+1)
    u1c = newOneOf g[0]  # U1(x+1+2)
    tmp = newOneOf g[0]
    a = newOneOf g[0]
    b = newOneOf g[0]
  threads:
    u2a := s1 ^* g[1]
    u1b := s2 ^* g[0]
    u1a := s1 ^* g[0]
    u2b := s1 ^* u2a
    u1c := s1 ^* u1b
    tmp := g[0] * u2a
    a := tmp * u1b.adj
    tmp := g[0] * u1a
    b := tmp * u2b
    tmp := b * u1c.adj
    b := tmp * u1b.adj

  test "deterministic":
    let t = paths.optimalPairs
    check(t.plan == t.plan)
    check(t.plan(false) == paths.optimalPairs.plan(false))
    check(t.plan != t.plan(false))

  test "self consistent":
    let pl = paths.optimalPairs.plan
    var keys = newseq[seq[int]]()
    for s in pl.steps:
      check(s.key == fold(s.l, s.la) & fold(s.r, s.ra))
      check(s.key.len > 1 and s.key notin keys)
      for k in [s.l, s.r]:
        check(k.len == 1 or k in keys)
        if k.len == 1: check(k[0] > 0)
      keys.add s.key
    check(pl.outs.len == paths.len)
    for i,o in pl.outs.pairs:
      check(fold(o.key, o.adj) == paths[i])
      check(o.key.len == 1 or o.key in keys)
      if o.key.len == 1: check(o.key[0] > 0)
    check(paths.optimalPairs.plan(false).steps == pl.steps)
    for o in paths.optimalPairs.plan(false).outs:
      check(o.sh.len == 0)

  test "products":
    let
      ps = @[plq, rct, rot[0], rot[1]]
      pt = ps.optimalPairs
      r1 = g.gaugeProd pt
      r0 = g.gaugeProd(pt, false)
      pl = pt.plan
      wl = g.wilsonLines ps
      p = g.plaq
    var
      tb2 = newTransporter(g[1], g[0], 1, -1)  # U2(x-2)^† y(x-2)
      tr: typeof(trace(g[0]))
    for i in 0..<ps.len:
      let st = if i mod 2 == 0: a else: b
      if i < 2:
        threads:
          tmp := st * g[1].adj
      else:
        threads:
          tmp := tb2 ^* st
      check(r1[i] ~ tmp)
      check(r1[i] ~ shifted(r0[i], pl.outs[i].sh))
      check(pl.outs[i].sh.len == (if i < 2: 0 else: 2))
      threads:
        tr = tmp.trace
      check(tr.re*fac ~ wl[i].re)
      check(tr.im*fac ~ wl[i].im)
      if i == 0:
        check(tr.re*fac/6.0 ~ p[0])

  test "adjoint output":
    let
      pt = optimalPairs([plq, @[2,1,-2,-1]])
      pl = pt.plan
      r = g.gaugeProd pt
    check(pl.outs[0].key == pl.outs[1].key)
    check(pl.outs[0].adj != pl.outs[1].adj)
    threads:
      tmp := r[0].adj
    check(r[1] ~ tmp)
    threads:
      tmp := a * g[1].adj
    check(r[0] ~ tmp)

qexFinalize()
