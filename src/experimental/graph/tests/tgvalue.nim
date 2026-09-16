#RUNCMD env OMP_NUM_THREADS=1 $RUNJOB
## Ordinary graph value lifetimes. QEX options configure the fixture; all cases run.
import base/globals
setVLENmax(4)

import math, unittest
import qex except epsilon
import base/alignedMem
import helpers
import ../[core, scalar, gauge, multi, functional]
import ../gauge/[types, basic_ops, field_ops, transport, stencil]

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))
qexInit()
letParam:
  expectRanks = nRanks
check nRanks == expectRanks
letParam:
  lat = latticeFromLocalLattice(@[4,4], nRanks)
let lo = lat.newLayout
var rng = lo.newRNGField(Philox4x64, 8091137u64)
let g = lo.newGauge
let u = lo.newGauge
let m = lo.newGauge
threads:
  g.random rng
  u.random rng
  m.randomTAH rng
  for f in m:
    f *= 0.1

proc distance(a, b: auto): float =
  let d = newOneOf(a)
  var value = 0.0
  threads:
    var n = 0.0
    for mu in 0..<a.len:
      d[mu] := a[mu] - b[mu]
      n += d[mu].norm2
    threadSingle: value = n
  value

proc same(x, y: Ggauge) =
  let a = x.eval.gval
  let b = y.eval.gval
  check distance(a, b) < 1e-18

proc same(x, y: Gscalar) =
  let a = x.eval.sval
  let b = y.eval.sval
  check abs(a-b) < 1e-10 * (1.0 + abs(b))

suite "graph value lifetime":
  setup:
    let rt = initGraphRuntime()
    let keep = initGraphRuntime()
    let x = rt.toGvalue(g)
    let y = rt.toGvalue(u)
    let b = rt.toGvalue(m)
    let kx = keep.toGvalue(g)
    let ky = keep.toGvalue(u)
    let kb = keep.toGvalue(m)
    defer:
      rt.resetGradCache
      rt.resetApplyCache
      rt.resetLdjCache
      keep.resetGradCache
      keep.resetApplyCache
      keep.resetLdjCache

  test "construction defers fields and generated identity values":
    let raw = getRawMemAllocated()
    let unit = x.unitGaugeLike
    let f = x*x + unit
    let moved = hop(x, shift(linkField(x, 0), 1, 1), 0, 1)
    let subset = gaugeActionDeriv(actWilson(rt.toGvalue(0.7)), x, 1, 0)
    let packed = multiValues("construction bundle", f, moved, subset)
    check getRawMemAllocated() == raw
    check not unit.hasStorage
    check not f.hasStorage
    check not moved.hasStorage
    check not subset.hasStorage
    check not packed.hasStorage
    check x.hasStorage
    check x.gval[0].s.data != g[0].s.data
    check f.nodeRepr.len > 0
    let score = norm2(f)
    same(score, norm2(kx*kx + kx.unitGaugeLike))
    check f.hasStorage
    check unit.hasStorage
    check x.hasStorage
    check score.hasStorage
    check score.valueReady
    let count = score.runCount
    let fieldCount = f.runCount
    discard score.eval
    check score.runCount == count
    check f.runCount == fieldCount

  test "repeated and alternating roots follow leaf updates":
    let a = x*y
    let f = norm2(a + x)
    let h = redot(a*a, y)
    let ka = kx*ky
    let kf = norm2(ka + kx)
    let kh = redot(ka*ka, ky)
    for i in 0..2:
      same(f, kf)
      same(h, kh)
      let count = a.runCount
      same(f, kf)
      check a.hasStorage
      check a.runCount == count
      if i == 0:
        x.update(u)
        kx.update(u)
      elif i == 1:
        y.update(g)
        ky.update(g)

  test "explicit writes to computed values remain authoritative":
    let a = x*x
    discard a.eval
    a.update(u)
    same(norm2(a + y), norm2(ky + ky))
    check a.hasStorage
    check a.valueOverride
    a.mutateGauge data:
      threads:
        for f in data:
          f *= 0.5
    same(norm2(a), norm2(0.5 * ky))
    check a.hasStorage
    x.update(u)
    kx.update(u)
    same(norm2(a), norm2(kx*kx))
    check not a.valueOverride
    check a.hasStorage

  test "partial mutations restore computed and identity values first":
    let a = x*x
    discard norm2(a).eval
    a.releaseStorage
    check not a.hasStorage
    a.mutateGauge data:
      threads:
        for f in data:
          f *= 0.5
    same(a, 0.5*(kx*kx))
    let unit = x.unitGaugeLike
    check not unit.hasStorage
    unit.mutateGauge data:
      threads:
        for f in data:
          f *= 2.0
    same(unit, 2.0*kx.unitGaugeLike)
    let f = linkField(x, 0)*linkField(y, 1)
    discard norm2(f).eval
    f.releaseStorage
    check not f.hasStorage
    f.mutateField data:
      threads: data *= 2.0
    let want = rt.toGvalue((2.0*linkField(kx, 0)*linkField(ky, 1)).eval.fval)
    check norm2(f-want).eval.sval < 1e-18

  test "lazy branches retain evaluated values and skip inactive forwards":
    let selector = rt.toGvalue(1)
    let a = x*x
    let z = x+y
    let score = norm2(cond(selector, a, z))
    same(score, norm2(kx*kx))
    check z.runCount == 0
    check not z.hasStorage
    check a.hasStorage
    selector.update 0
    same(score, norm2(kx+ky))
    check z.runCount == 1
    check z.hasStorage
    let count = a.runCount
    selector.update 1
    same(score, norm2(kx*kx))
    check a.runCount == count

  test "subset and injected complements are zero after reallocation":
    let a = maskSubset(1, 0, x*y)
    let i = injectLink(linkField(x, 1), 0, x)
    let ka = maskSubset(1, 0, kx*ky)
    let ki = injectLink(linkField(kx, 1), 0, kx)
    for _ in 0..1:
      same(a, ka)
      same(i, ki)
      a.releaseStorage
      i.releaseStorage
      check not a.hasStorage
      check not i.hasStorage
      x.update(u)
      kx.update(u)

  test "an expired injected override clears the untouched directions":
    let input = rt.toGvalue(g[0])
    let injected = injectLink(input, 0, x)
    discard injected.eval
    injected.update(u)
    check injected.valueOverride
    check norm2(injected.gval[1]) > 0.0
    input.update(u[0])
    discard injected.eval
    check not injected.valueOverride
    let delta = u[0].newOneOf
    threads: delta := injected.gval[0] - u[0]
    check norm2(delta) < 1e-18
    for mu in 1..<injected.gval.len:
      check norm2(injected.gval[mu]) == 0.0

  test "an accumulating forward clears a failed destination before retry":
    proc accumulate(v: Gvalue) =
      let z = Gfield(v)
      let input = Gfield(v.inputs[0])
      threads: z.fval += input.fval
      if Gint(v.inputs[1]).ival != 0:
        raiseError("requested accumulating forward failure")
    let input = rt.toGvalue(g[0])
    let flag = rt.toGvalue(1)
    let z = graphNode(input.fieldNodeLike, [Gvalue(input), Gvalue(flag)],
      Gfunc(forward: accumulate, bufferMode: bmZero, name: "lifetime accumulate"),
      "lifetime accumulate")
    expect(GraphError):
      discard z.eval
    check z.hasStorage
    check not z.valueReady
    flag.update(0)
    discard z.eval
    let delta = g[0].newOneOf
    threads: delta := z.fval - g[0]
    check norm2(delta) < 1e-18
    input.update(u[0])
    discard z.eval
    threads: delta := z.fval - u[0]
    check norm2(delta) < 1e-18

  test "shift and transport aliases retain data after owner release":
    let f = shift(linkField(x, 0), 1, 1)
    let h = hop(x, f, 0, 1)
    let kf = shift(linkField(kx, 0), 1, 1)
    let kh = hop(kx, kf, 0, 1)
    discard h.eval
    let snapshot = h.fval
    let reference = kh.eval.fval
    let want = reference.newOneOf
    threads: want := reference
    h.releaseStorage
    check not h.hasStorage
    let delta = snapshot.newOneOf
    threads: delta := snapshot - want
    check norm2(delta) < 1e-18
    discard h.eval
    x.update(u)
    kx.update(u)
    let got = rt.toGvalue(h.eval.fval)
    let refv = rt.toGvalue(kh.eval.fval)
    check norm2(got-refv).eval.sval < 1e-18
    threads: delta := snapshot - want
    check norm2(delta) < 1e-18

  test "halo gather scatter and products rebuild released work":
    let f = gather(linkField(x, 0), @[1, 0])
    let kf = gather(linkField(kx, 0), @[1, 0])
    let back = scatter(f, @[1, 0])
    let kback = scatter(kf, @[1, 0])
    let moved = gp(f, linkField(y, 1), @[1, -1])
    let kmoved = gp(kf, linkField(ky, 1), @[1, -1])
    let score = redot(moved + back, linkField(b, 0))
    let want = redot(kmoved + kback, linkField(kb, 0))
    for _ in 0..1:
      same(score, want)
      f.releaseStorage
      back.releaseStorage
      moved.releaseWork
      check not f.hasStorage
      check not back.hasStorage
      x.update(u)
      kx.update(u)

  test "field updates and mutations restore an explicitly released destination":
    let f = linkField(x, 0) * linkField(y, 1)
    discard norm2(f).eval
    f.releaseStorage
    check not f.hasStorage
    f.update(g[0])
    let want = rt.toGvalue(g[0])
    check norm2(f-want).eval.sval < 1e-18
    check f.hasStorage
    f.mutateField data:
      threads: data *= 0.5
    check norm2(f - 0.5*want).eval.sval < 1e-18
    check f.valueOverride

  test "packed values own wrappers and preserve explicitly released payloads":
    let a = x*y
    let args = multiValues("lifetime bundle", a, b)
    discard args.eval
    check a.hasStorage
    let packed = Ggauge(args.storedSlot(0))
    check packed.nodeKey != a.nodeKey
    let saved = packed.gval
    let want = (kx*ky).eval.gval
    check distance(saved, want) < 1e-18
    a.releaseStorage
    check not a.hasStorage
    check args.hasStorage
    args.releaseStorage
    check not args.hasStorage
    check distance(saved, want) < 1e-18
    discard args.eval
    check distance(Ggauge(args.storedSlot(0)).gval, want) < 1e-18
    let alias = slotVar(args)
    same(Ggauge(alias[0]), kx*ky)

  test "nested identity applications retain live inputs":
    let p = Ggauge(x.newOneOf)
    let identity = lambda(p, p)
    let fn = lambda(p, p*p + p.unitGaugeLike)
    let a = x*y
    let copied = Ggauge(apply(identity, a))
    let made = Ggauge(apply(fn, a))
    let score = norm2(made) + norm2(copied) + norm2(a)
    let ka = kx*ky
    let want = norm2(ka*ka + ka.unitGaugeLike) + 2.0*norm2(ka)
    for _ in 0..2:
      same(score, want)
      check a.hasStorage
      check copied.hasStorage
      check made.hasStorage
      let count = a.runCount
      discard score.eval
      check a.runCount == count
      x.update(u)
      kx.update(u)

  test "new derivatives can be built after earlier values are released":
    let a = x*y + x*x
    let ka = kx*ky + kx*kx
    let score = norm2(a)
    let want = norm2(ka)
    same(score, want)
    a.releaseStorage
    check not a.hasStorage
    let d = grad(score, x)
    let kd = grad(want, kx)
    same(d, kd)
    d.releaseStorage
    let dd = grad(redot(d, b), x)
    let kdd = grad(redot(kd, kb), kx)
    same(dd, kdd)
    x.update(u)
    kx.update(u)
    same(dd, kdd)

  test "stout packed cache payloads survive release and repeated pullbacks":
    let c = actWilson(rt.toGvalue(0.7))
    let kc = actWilson(keep.toGvalue(0.7))
    let alpha = rt.toGvalue(0.02)
    let kalpha = keep.toGvalue(0.02)
    let st = stoutUpdateLogDetJ(x, c, alpha, 1, 0)
    let ks = stoutUpdateLogDetJ(kx, kc, kalpha, 1, 0)
    let score = redot(st.Wnew, b) + st.lj
    let want = redot(ks.Wnew, kb) + ks.lj
    let dx = grad(score, x)
    let kd = grad(want, kx)
    let da = grad(score, alpha)
    let kad = grad(want, kalpha)
    for _ in 0..1:
      same(score, want)
      same(dx, kd)
      same(da, kad)
      same(st.Wnew, ks.Wnew)
      st.Wnew.releaseStorage
      dx.releaseStorage
      alpha.update 0.015
      kalpha.update 0.015
      x.update(u)
      kx.update(u)

  test "restoring unchanged storage preserves cached consumers and epochs":
    let a = x*x
    let score = norm2(a)
    same(score, norm2(kx*kx))
    let ep = a.epoch
    let ac = a.runCount
    let sc = score.runCount
    a.releaseStorage
    check a.epoch == ep
    same(score, norm2(kx*kx))
    check a.epoch == ep
    check a.runCount == ac+1
    check score.runCount == sc

  test "expired gauge overrides invalidate ordinary cached consumers":
    let source = lo.newGauge
    let overrideValue = lo.newGauge
    threads:
      for f in source: f := 2.0
      for f in overrideValue: f := 3.0
    let input = rt.toGvalue(source)
    let a = input*input
    discard a.eval
    a.update(overrideValue)
    let score = norm2(a)
    let oldValue = score.eval.sval
    let ep = a.epoch
    let sc = score.runCount
    let saved = a.gval
    a.releaseStorage
    check a.epoch == ep
    check distance(saved, overrideValue) < 1e-18
    let fresh = norm2(keep.toGvalue(source)*keep.toGvalue(source))
    same(score, fresh)
    check score.sval > oldValue
    check a.epoch > ep
    check score.runCount == sc+1
    check not a.valueOverride
    check distance(saved, overrideValue) < 1e-18

  test "generated identities restore with stable epochs":
    let unit = x.unitGaugeLike
    let field = x.unitFieldLike
    discard unit.eval
    discard field.eval
    let ep = unit.epoch
    let fp = field.epoch
    unit.releaseStorage
    field.releaseStorage
    same(unit, kx.unitGaugeLike)
    check norm2(field - rt.toGvalue(kx.unitFieldLike.eval.fval)).eval.sval < 1e-18
    check unit.epoch == ep
    check field.epoch == fp
    check unit.restoreValue != nil
    check field.restoreValue != nil

  test "subset producer work restores after storage and work release":
    let c = actWilson(rt.toGvalue(0.7))
    let kc = actWilson(keep.toGvalue(0.7))
    let a = gaugeActionDeriv(c, x, 1, 0)
    let ka = gaugeActionDeriv(kc, kx, 1, 0)
    let d = Ggauge(gradSeeded(a, x, b+b))
    let kd = Ggauge(gradSeeded(ka, kx, kb+kb))
    for _ in 0..1:
      same(a, ka)
      same(d, kd)
      a.releaseStorage
      d.releaseStorage
      check not a.hasStorage
      check not d.hasStorage
      a.releaseWork
      d.releaseWork
      x.update(u)
      kx.update(u)

  test "cond apply and slot results carry plain values from cached producers":
    let a = rt.toGvalue(0.04)
    let ka = keep.toGvalue(0.04)
    let sel = rt.toGvalue(1)
    let v = Ggauge(x.newOneOf)
    let cached = stoutUpdate(x, y, a, 1, 1)
    let kc = stoutUpdate(kx, ky, ka, 1, 1)
    let call = Ggauge(apply(lambda(v, stoutUpdate(v, y, a, 1, 1)), x))
    let choice = cond(sel, cached, x)
    let slot = slotVar(cached)
    check cached.bufferBytes > x.bufferBytes
    check call.bufferBytes == x.bufferBytes
    check choice.bufferBytes == x.bufferBytes
    check slot.bufferBytes == x.bufferBytes
    for flag in [1, 0, 1]:
      sel.update(flag)
      same(call, kc)
      same(slot, kc)
      same(choice, (if flag == 1: kc else: kx))
      call.releaseStorage
      choice.releaseStorage
      slot.releaseStorage
      cached.releaseStorage
      x.update(u)
      kx.update(u)


  test "nested mixed bundles release only their owned numerical wrappers":
    let p = rt.toGvalue(2.0)
    let square = lambda(p, p*p)
    let cube = lambda(p, p*p*p)
    let a = x*y
    let inner = multiValues("inner mixed bundle", square, a)
    proc forward(v: Gvalue) = discard
    let structure = newMultiStructureNode([Gvalue(a)], [Gvalue(a)],
      Gfunc(forward: forward, name: "lifetime structure"), "lifetime structure")
    let outer = multiValues("outer mixed bundle", inner, structure, x)
    discard outer.eval
    let owned = Gmulti(outer.storedSlot(0))
    let sourceField = a.gval[0]
    let innerField = Ggauge(inner.storedSlot(1)).gval[0]
    check owned.nodeKey != inner.nodeKey
    check owned.storedSlot(1).nodeKey != inner.storedSlot(1).nodeKey
    check outer.storedSlot(1).nodeKey == structure.nodeKey
    apply(owned.storedSlot(0), 3.0) :~ 9.0
    outer.releaseStorage
    outer.ensureStorage
    check a.hasStorage
    check inner.hasStorage
    check x.hasStorage
    check a.gval[0].s.data == sourceField.s.data
    check Ggauge(inner.storedSlot(1)).gval[0].s.data == innerField.s.data
    check not owned.hasStorage
    let replacement = multiValues("replacement inner", cube, y)
    let cloned = graphNode(Gmulti(outer.newOneOf),
      [Gvalue(replacement), Gvalue(structure), Gvalue(x)], outer.gfunc,
      "cloned outer mixed bundle")
    discard cloned.eval
    let nested = Gmulti(cloned.storedSlot(0))
    check nested.nodeKey != replacement.nodeKey
    check nested.storedSlot(0).nodeKey == cube.nodeKey
    apply(nested.storedSlot(0), 3.0) :~ 27.0
    same(Ggauge(nested.storedSlot(1)), ky)
    cloned.releaseStorage
    check replacement.hasStorage
    check y.hasStorage
    check a.hasStorage

qexFinalize()
