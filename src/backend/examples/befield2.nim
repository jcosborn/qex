## Gauge kernels in the style of field code against the same kernels on
## arrays of reals (gauge/gaugeGpu), on ColorMatrix fields in device memory
## in the layout of QEX fields:
##   axpy  p = p + t f
##   tah   p = p - t TAH(u w^+)
##   prod  x = u v w^+
##   exp   u = exp(t p) u
## Field style: QEX site types and operators on GpuField views, applied to
## the site values in place (direct) or to copies loaded first (load), with
## one site of the four directions per thread, or one link per thread on a
## field of all links (links).  The array style loads 18 reals per link and
## calls the procs of gaugeGpu, one link per thread; exp: the polynomial and
## squaring of gaugeGpu before it used QEX's expAH, against expAH.
##   -lat: lattice (default 16,16,16,16)  -reps: repetitions (default 20)
import qex, gauge, physics/qcdTypes
import backend/[accel, cgfield]
import gauge/gaugeGpu
import base/metaUtils
import times
getOptimPragmas()

qexInit()
let lat = intSeqParam("lat", @[16,16,16,16])
let reps = intParam("reps", 20)
let lo = lat.newLayout
const V = typeof(lo).V
let n = lo.nSites
let nl = 4*18*n
var rng = lo.newRNGField(RngMilc6, 987654321)
var u = lo.newGauge
var v = lo.newGauge
var w = lo.newGauge
var f = lo.newGauge
var p = lo.newGauge
u.random rng
v.random rng
w.random rng
threads:
  u.projectSU
  v.projectSU
  w.projectSU
  for mu in 0..3:
    for e in f[mu]: f[mu][e].projectTAH(v[mu][e])
    for e in p[mu]: p[mu][e].projectTAH(w[mu][e])

type
  Arr = ptr UncheckedArray[float]
  GF = gpuType(typeof(u[0]))
  GV = array[4, GF]

proc dev(x: seq): Arr =
  ## the 4 fields of x in one device array, [4][n div V][18][V]
  result = cast[Arr](gpuMalloc(nl*sizeof(float)))
  for mu in 0..3: gpuMemCpyToGpu(addr result[18*mu*n], addr x[mu][0], 18*n*sizeof(float))

proc host(x: seq; d: Arr) =
  for mu in 0..3: gpuMemCpyToCpu(addr x[mu][0], addr d[18*mu*n], 18*n*sizeof(float))

proc views(d: Arr): GV =
  ## field views of the 4 fields of d
  for mu in 0..3: result[mu] = gpuFieldView(lo, typeof(u[0][0]), addr d[18*mu*n])

proc allLinks(d: Arr): GF =
  ## the links of d as one field of 4 times the sites
  result = gpuFieldView(lo, typeof(u[0][0]), d)
  result.d.n = 4*result.d.n
  result.d.ne = 0

proc mexp(r: var array[18, float]; a: array[18, float]) {.gpuInline.} =
  ## r = exp(a) for anti-Hermitian a: exp(a/2^k) by its Taylor polynomial of
  ## degree 12 with |a/2^k|^2 <= 1/16, then squared k times (the link update
  ## of gaugeGpu before it used expAH)
  var n2 = 0.0
  forStatic i, 0, 17: n2 += a[i]*a[i]
  var s = 1.0
  var k = 0
  while n2 > 1.0/16.0:
    n2 *= 0.25
    s *= 0.5
    inc k
  var m {.noInit.}: array[18, float]
  forStatic i, 0, 17: m[i] = s*a[i]
  var m2 {.noInit.}, m3 {.noInit.}, t {.noInit.}: array[18, float]
  mmul(m2, m, m)
  mmul(m3, m2, m)
  # r = 1 + m + m2/2! + ... + m12/12!, Horner in m3 as expPoly12
  template poly(d, c2, c1, c0: untyped) =
    forStatic i, 0, 17: d[i] = c2*m2[i] + c1*m[i]
    forStatic i, 0, 2: d[8*i] += c0
  const f = [1.0, 1.0, 1.0/2, 1.0/6, 1.0/24, 1.0/120, 1.0/720, 1.0/5040, 1.0/40320,
             1.0/362880, 1.0/3628800, 1.0/39916800, 1.0/479001600]
  var c {.noInit.}: array[18, float]
  poly(c, f[11], f[10], f[9])
  forStatic i, 0, 17: c[i] += f[12]*m3[i]
  mmul(t, c, m3)
  poly(c, f[8], f[7], f[6])
  forStatic i, 0, 17: c[i] += t[i]
  mmul(t, c, m3)
  poly(c, f[5], f[4], f[3])
  forStatic i, 0, 17: c[i] += t[i]
  mmul(t, c, m3)
  poly(r, f[2], 1.0, 1.0)
  forStatic i, 0, 17: r[i] += t[i]
  for j in 0..<k:
    mmul(t, r, r)
    forStatic i, 0, 17: r[i] = t[i]

proc copy(d, s: Arr) =
  onGpu(nl):
    for i in gpuRange(nl): d[i] = s[i]

let du = dev(u)
let dv = dev(v)
let dw = dev(w)
let df = dev(f)
let dp = dev(p)
let dq = dev(p)  # the result of the array kernels
let dx = dev(p)
let dy = dev(p)
let uv = views(du)
let vv = views(dv)
let wv = views(dw)
let fv = views(df)
let pv = views(dp)
let xv = views(dx)
let ul = allLinks(du)
let vl = allLinks(dv)
let wl = allLinks(dw)
let xl = allLinks(dx)
let fl = allLinks(df)
let dz = dev(u)  # the links updated by exp
let zl = allLinks(dz)

# array style, one link per thread as in gaugeGpu
proc axpyA(p, f: Arr; t: float) =
  onGpu(nl):
    for i in gpuRange(nl): p[i] += t*f[i]

proc tahA(p, u, w: Arr; t: float) =
  onGpu(4*n):
    for i in gpuRange(4*n):
      let mu = i div n
      let o = 18*mu*n + lo18(V, i - mu*n)
      var a {.noInit.}, b {.noInit.}, c {.noInit.}, s {.noInit.}: array[18, float]
      mload(a, u, o, V)
      mload(b, w, o, V)
      mmulNA(c, a, b)
      mtah(s, c)
      forStatic e, 0, 17: p[o + e*V] -= t*s[e]

proc prodA(x, u, v, w: Arr) =
  onGpu(4*n):
    for i in gpuRange(4*n):
      let mu = i div n
      let o = 18*mu*n + lo18(V, i - mu*n)
      var a {.noInit.}, b {.noInit.}, c {.noInit.}, d {.noInit.}: array[18, float]
      mload(a, u, o, V)
      mload(b, v, o, V)
      mmul(c, a, b)
      mload(a, w, o, V)
      mmulNA(d, c, a)
      forStatic e, 0, 17: x[o + e*V] = d[e]

# field style, one site of each direction per thread
proc axpyF(p, f: GV; t: float) =
  onGpu(lo):
    for mu in 0..3: p[mu] += t*f[mu]

proc tahF(p, u, w: GV; t: float) =
  onGpu(lo):
    for mu in 0..3:
      for s in p[mu]:
        var m {.noInit.}: evalType(u[mu][s])
        m.projectTAH(u[mu][s] * w[mu][s].adj)
        p[mu][s] -= t*m

proc tahL(p, u, w: GV; t: float) =
  onGpu(lo):
    for mu in 0..3:
      for s in p[mu]:
        var a {.noInit.}, b {.noInit.}, m {.noInit.}: evalType(u[mu][s])
        a := u[mu][s]
        b := w[mu][s]
        m.projectTAH(a * b.adj)
        p[mu][s] -= t*m

proc prodF(x, u, v, w: GV) =
  onGpu(lo):
    for mu in 0..3: x[mu] := u[mu] * v[mu] * w[mu].adj

proc prodL(x, u, v, w: GV) =
  onGpu(lo):
    for mu in 0..3:
      for s in x[mu]:
        var a {.noInit.}, b {.noInit.}, c {.noInit.}: evalType(u[mu][s])
        a := u[mu][s]
        b := v[mu][s]
        c := w[mu][s]
        x[mu][s] := a * b * c.adj

proc tahK(p, u, w: GF; t: float) =
  onGpu(4*n):
    for s in p:
      var m {.noInit.}: evalType(u[s])
      m.projectTAH(u[s] * w[s].adj)
      p[s] -= t*m

proc prodK(x, u, v, w: GF) =
  onGpu(4*n):
    x := u * v * w.adj

proc expA(u, p: Arr; t: float) =
  onGpu(4*n):
    for i in gpuRange(4*n):
      let mu = i div n
      let o = 18*mu*n + lo18(V, i - mu*n)
      var a {.noInit.}, e {.noInit.}, x {.noInit.}, y {.noInit.}: array[18, float]
      forStatic j, 0, 17: a[j] = t*p[o + j*V]
      mexp(e, a)
      mload(x, u, o, V)
      mmul(y, e, x)
      forStatic j, 0, 17: u[o + j*V] = y[j]

proc expK(u, p: GF; t: float) =
  onGpu(4*n):
    for s in u:
      u[s] := expAH(t*p[s]) * u[s]

proc rel(a, b: Arr): float =
  ## |a-b|^2/|b|^2 of two device arrays like dp
  var x = lo.newGauge
  var y = lo.newGauge
  host(x, a)
  host(y, b)
  var d, s = 0.0
  threads:
    var dt, st = 0.0
    for mu in 0..3:
      dt += norm2(x[mu] - y[mu])
      st += y[mu].norm2
    threadMaster:
      d = dt
      s = st
  d/s

let t = 0.25
copy(dq, dp); axpyA(dq, df, t)
copy(dx, dp); axpyF(xv, fv, t)
echo "axpy field/array |d|^2/|a|^2: ", rel(dx, dq)
copy(dq, dp); tahA(dq, du, dw, t)
copy(dx, dp); tahF(xv, uv, wv, t)
echo "tah field direct/array: ", rel(dx, dq)
copy(dx, dp); tahL(xv, uv, wv, t)
echo "tah field load/array: ", rel(dx, dq)
prodA(dq, du, dv, dw)
prodF(xv, uv, vv, wv)
echo "prod field direct/array: ", rel(dx, dq)
prodL(xv, uv, vv, wv)
echo "prod field load/array: ", rel(dx, dq)
prodK(xl, ul, vl, wl)
echo "prod field links/array: ", rel(dx, dq)
copy(dq, dp); tahA(dq, du, dw, t)
copy(dx, dp); tahK(xl, ul, wl, t)
echo "tah field links/array: ", rel(dx, dq)
copy(dq, du); expA(dq, df, t)
copy(dz, du); expK(zl, fl, t)
echo "exp field links (expAH)/array: ", rel(dz, dq)

template bench(name: string; fields: int; body: untyped) =
  ## the best of reps runs; bytes: fields read or written of 4 directions
  body
  var best = 0.0
  for r in 0..<reps:
    let t0 = epochTime()
    body
    let dt = epochTime() - t0
    if r == 0 or dt < best: best = dt
  let bytes = float(fields*nl*sizeof(float))
  echo name, ": ", 1e6*best, " us  ", 1e-9*bytes/best, " GB/s"

echo "lattice ", lat, " V ", V
bench("axpy array", 3): axpyA(dq, df, 1e-9)
bench("axpy field", 3): axpyF(xv, fv, 1e-9)
bench("tah array", 4): tahA(dq, du, dw, 1e-9)
bench("tah field direct", 4): tahF(xv, uv, wv, 1e-9)
bench("tah field load", 4): tahL(xv, uv, wv, 1e-9)
bench("prod array", 4): prodA(dq, du, dv, dw)
bench("prod field direct", 4): prodF(xv, uv, vv, wv)
bench("prod field load", 4): prodL(xv, uv, vv, wv)
bench("tah field links", 4): tahK(xl, ul, wl, 1e-9)
bench("prod field links", 4): prodK(xl, ul, vl, wl)
bench("exp array", 3): expA(dq, df, 1e-3)
bench("exp field links", 3): expK(zl, fl, 1e-3)
when not backendIsGpu:
  var x = lo.newGauge
  bench("axpy threads", 3):
    threads:
      for mu in 0..3: p[mu] += 1e-9*f[mu]
  bench("tah threads", 4):
    threads:
      for mu in 0..3:
        for e in p[mu]:
          var m {.noInit.}: evalType(u[mu][e])
          m.projectTAH(u[mu][e] * w[mu][e].adj)
          p[mu][e] -= 1e-9*m
  bench("prod threads", 4):
    threads:
      for mu in 0..3:
        for e in x[mu]:
          x[mu][e] := u[mu][e] * v[mu][e] * w[mu][e].adj
  bench("exp threads", 3):
    threads:
      for mu in 0..3:
        for e in x[mu]:
          x[mu][e] := expAH(1e-3*f[mu][e]) * x[mu][e]
qexFinalize()
