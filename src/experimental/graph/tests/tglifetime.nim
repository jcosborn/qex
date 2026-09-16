#RUNCMD env OMP_NUM_THREADS=1 $RUNJOB
## Graph buffer lifetimes on a small valid SIMD lattice.
import base/globals
setVLENmax(4)

import math, unittest
import qex except epsilon
import base/alignedMem
import helpers
import ../[core, scalar, gauge, multi, functional, plan]
import ../gauge/[types, field_ops, transport, stencil, matrix]

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

letParam:
  measure = false
if measure:
  letParam:
    mode = "planned"
  let planned = case mode
    of "direct": false
    of "planned": true
    else: raiseValueError("mode must be direct or planned")
  let rt = initGraphRuntime()
  let x = rt.toGvalue(g)
  let b = rt.toGvalue(m)
  let c = actWilson(rt.toGvalue(0.7))
  GC_fullCollect()
  let initial = getRawMemUsed()
  let f = gaugeActionDeriv(c, x)
  let first = grad(norm2(f), x)
  let second = grad(redot(first, b), x)
  let score = norm2(second)
  let p = if planned: plan(score) else: nil
  let built = getRawMemUsed() - initial
  var value: float
  if p == nil:
    value = score.eval.sval
  else:
    discard p.eval
    value = Gscalar(p[0]).sval
  let peak = getRawMemMaxUsed() - initial
  GC_fullCollect()
  let retained = getRawMemUsed() - initial
  echo "graph tower memory: mode=", mode, ", constructed=", built,
       ", peak=", peak, ", retained=", retained, ", value=", value
  if p != nil: p.clear
  rt.resetGradCache
  rt.resetApplyCache
  rt.resetLdjCache
else:
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
      let p = plan(f, moved, subset, packed)
      check getRawMemAllocated() == raw
      check not unit.hasStorage
      check not f.hasStorage
      check not moved.hasStorage
      check not subset.hasStorage
      check not packed.hasStorage
      check x.hasStorage
      check x.gval[0].s.data != g[0].s.data
      p.clear
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

    test "halo moves and fused plaquettes rebuild explicitly released work":
      let f = linkField(x, 0)
      let kf = linkField(kx, 0)
      let moved = gp(f, linkField(y, 1), @[1,-1])
      let kmoved = gp(kf, linkField(ky, 1), @[1,-1])
      let score = redot(moved, linkField(b, 0))
      let want = redot(kmoved, linkField(kb, 0))
      let a = x*x
      let p = stencil.plaqSum(a)
      let kp = stencil.plaqSum(kx*kx)
      for _ in 0..1:
        same(score, want)
        check moved.hasStorage
        same(p, kp)
        check a.hasStorage
        let count = p.runCount
        discard p.eval
        check p.runCount == count
        same(grad(p, x), grad(kp, kx))
        moved.releaseStorage
        p.releaseWork
        a.releaseStorage
        check not moved.hasStorage
        check not a.hasStorage
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

    test "functional clones and nested identity apply retain live inputs":
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

    test "matrix8 and scalar-field values restore through derivatives":
      let raw = lo.RealMatrix(8)
      threads:
        raw := 2.0
        for e in raw:
          raw[e][0,1] := 0.1
      let a = rt.toGvalue(raw)
      let ka = keep.toGvalue(raw)
      let prod = a*a
      let score = sumLogDet(prod) + norm2(inverse(prod))
      let want = sumLogDet(ka*ka) + norm2(inverse(ka*ka))
      let d = grad(score, a)
      let kd = grad(want, ka)
      same(score, want)
      prod.releaseStorage
      check not prod.hasStorage
      let dd = rt.toGvalue(kd.eval.fval)
      check norm2(d-dd).eval.sval < 1e-18
      d.releaseStorage
      same(score, want)
      check norm2(d-dd).eval.sval < 1e-18

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

qexFinalize()
