## Gauge links on the GPU for HMC.
##
## Links and momenta of all directions are one device array each,
## [4][n div V][18][V] in the qex SIMD layout of the gauge field, so uploads
## are plain copies: real e = 6a+2b of U_ab (e+1 its imaginary part) at site
## k of direction mu is at 18*mu*n + (k div V)*18V + e*V + k mod V.
##
## Gauge action plaq + adjplaq, as actionA and forceA with nc = 3:
##   S = c_p sum_P (1 - Re tr P/3) + c_a sum_P (1 - |tr P|^2/9)
##   F_mu(x) = TAH sum_P conj(c_p/3 + 2 c_a/9 tr P) P
## over the plaquettes P = U_mu(x) U_nu(x+mu) U_mu(x+nu)^+ U_nu(x)^+ and
## U_mu(x) U_nu(x+mu-nu)^+ U_mu(x-nu)^+ U_nu(x-nu), nu != mu.
import qex
import gauge
import comms/[halo, halogpu]
import backend/accel
import base/metaUtils

type
  GpuGauge*[V: static int] = object
    lo*: Layout[V]
    n*: int  # local sites
    u*: ptr UncheckedArray[float]  # [4][n div V][18][V]
    ex*: array[4, GpuHaloEx[float]]  # halo of each direction, with the corners x+mu-nu
    nb*: ptr UncheckedArray[int32]  # [20][n]: x+mu, x-mu, x+mu-nu; j >= n is receive position j-n
    ord*: ptr UncheckedArray[int32]  # sites, the nin with all 20 neighbors local first
    nin*: int
    fresh*: bool  # the halo holds the current links

template lo18*(V, k: untyped): untyped =
  ## real 0 of site k in a link field, real e is e*V further
  (k div V)*(18*V) + k mod V

template nbF*(mu: untyped): untyped = mu
template nbB*(mu: untyped): untyped = 4 + mu
template nbD*(mu, nu: untyped): untyped = 8 + 3*mu + (if nu < mu: nu else: nu-1)  # x+mu-nu

proc newGpuGauge*[V: static int](lo: Layout[V]): GpuGauge[V] =
  tic("newGpuGauge")
  let nd = lo.nDim
  let no = lo.nSitesOuter
  let n = V*no
  let c = getDefaultComm()
  result.lo = lo
  result.n = n
  var w = newSeq[int](nd)
  for d in 0..<nd: w[d] = 1
  let hl = lo.makeHaloLayout(w, w)
  var offs = newSeq[seq[int32]](20)
  for mu in 0..<nd:
    offs[nbF(mu)] = newSeq[int32](nd)
    offs[nbF(mu)][mu] = 1
    offs[nbB(mu)] = newSeq[int32](nd)
    offs[nbB(mu)][mu] = -1
    for nu in 0..<nd:
      if nu != mu:
        var o = newSeq[int32](nd)
        o[mu] = 1
        o[nu] = -1
        offs[nbD(mu,nu)] = o
  let hm = hl.makeHaloMap(c, offs)
  for mu in 0..<nd:
    result.ex[mu] = newGpuHaloEx[float](hm.gather, 18, n, V, c)
  let src = hl.haloSource(hm.gather)
  var nb = newSeq[int32](20*n)
  for o in 0..<no:
    for mu in 0..<nd:
      let ef = int hl.neighborFwd[mu][o]
      for d in [nbF(mu), nbB(mu)]:
        let e = if d == nbF(mu): ef else: int hl.neighborBck[mu][o]
        for l in 0..<V:
          nb[d*n + V*o + l] = if e < no: int32(V*e + l) else: src[V*(e-no) + l]
      for nu in 0..<nd:
        if nu != mu:
          let e = int hl.neighborBck[nu][ef]
          for l in 0..<V:
            nb[nbD(mu,nu)*n + V*o + l] = if e < no: int32(V*e + l) else: src[V*(e-no) + l]
  result.nb = nb.toDevice
  var od = newSeq[int32](n)
  var m = 0
  for pass in 0..1:
    for k in 0..<n:
      var loc = true
      for d in 0..<20: loc = loc and nb[d*n + k] < int32(n)
      if loc == (pass == 0):
        od[m] = int32 k
        inc m
    if pass == 0: result.nin = m
  result.ord = od.toDevice
  result.u = cast[ptr UncheckedArray[float]](gpuMalloc(4*18*n*sizeof(float)))
  toc("done")

proc free*[V: static int](g: var GpuGauge[V]) =
  gpuFree(g.u)
  gpuFree(g.nb)
  gpuFree(g.ord)
  for mu in 0..3: g.ex[mu].free

proc upload*[V: static int](g: var GpuGauge[V]; d: ptr UncheckedArray[float]; f: openArray[Field]) =
  ## the links f to d, [4][n div V][18][V] as g.u
  for mu in 0..3:
    gpuMemCpyToGpu(addr d[18*mu*g.n], addr f[mu][0], 18*g.n*sizeof(float))
  if d == g.u: g.fresh = false

proc download*[V: static int](g: GpuGauge[V]; f: openArray[Field]; d: ptr UncheckedArray[float]) =
  for mu in 0..3:
    gpuMemCpyToCpu(addr f[mu][0], addr d[18*mu*g.n], 18*g.n*sizeof(float))

proc newLinks*[V: static int](g: GpuGauge[V]): ptr UncheckedArray[float] =
  ## a device array like g.u, for momenta or saved links
  cast[ptr UncheckedArray[float]](gpuMalloc(4*18*g.n*sizeof(float)))

proc copy*[V: static int](g: GpuGauge[V]; d, s: ptr UncheckedArray[float]) =
  ## d = s for arrays like g.u
  gpuFor(i, 4*18*g.n): d[i] = s[i]

proc update*[V: static int](g: var GpuGauge[V]) =
  ## Exchanges the halo of g.u.  The barrier keeps a peer from rewriting the
  ## halo while a kernel of this rank still reads the previous one.
  if g.fresh: return
  tic("gauge halo")
  getDefaultComm().barrier
  for mu in 0..3: g.ex[mu].pack(cast[ptr UncheckedArray[float]](addr g.u[18*mu*g.n]))
  for mu in 0..3: g.ex[mu].start
  for mu in 0..3: g.ex[mu].wait
  g.fresh = true
  toc("done")

# 3x3 complex matrices as 18 reals, U_ab at 6a+2b (re) and 6a+2b+1 (im)

template mload*(m, u, o, st: untyped) =
  forStatic i, 0, 17: m[i] = u[o + i*st]

template mmul*(r, a, b: untyped) =
  ## r = a b
  forStatic i, 0, 2:
    forStatic j, 0, 2:
      var re = a[6*i]*b[2*j] - a[6*i+1]*b[2*j+1]
      var im = a[6*i]*b[2*j+1] + a[6*i+1]*b[2*j]
      forStatic k, 1, 2:
        re += a[6*i+2*k]*b[6*k+2*j] - a[6*i+2*k+1]*b[6*k+2*j+1]
        im += a[6*i+2*k]*b[6*k+2*j+1] + a[6*i+2*k+1]*b[6*k+2*j]
      r[6*i+2*j] = re
      r[6*i+2*j+1] = im

template mmulNA*(r, a, b: untyped) =
  ## r = a b^+
  forStatic i, 0, 2:
    forStatic j, 0, 2:
      var re = a[6*i]*b[6*j] + a[6*i+1]*b[6*j+1]
      var im = a[6*i+1]*b[6*j] - a[6*i]*b[6*j+1]
      forStatic k, 1, 2:
        re += a[6*i+2*k]*b[6*j+2*k] + a[6*i+2*k+1]*b[6*j+2*k+1]
        im += a[6*i+2*k+1]*b[6*j+2*k] - a[6*i+2*k]*b[6*j+2*k+1]
      r[6*i+2*j] = re
      r[6*i+2*j+1] = im

template mtah*(r, a: untyped) =
  ## r = (a - a^+)/2 - tr(a - a^+)/6
  let t = (a[1] + a[9] + a[17])*(1.0/3.0)
  forStatic i, 0, 2:
    forStatic j, 0, 2:
      when i == j:
        r[8*i] = 0.0
        r[8*i+1] = a[8*i+1] - t
      else:
        r[6*i+2*j] = 0.5*(a[6*i+2*j] - a[6*j+2*i])
        r[6*i+2*j+1] = 0.5*(a[6*i+2*j+1] + a[6*j+2*i+1])

template mexp*(r, a: untyped) =
  ## r = exp(a) for anti-Hermitian a, as expAH: exp(a/2^k) by its Taylor
  ## polynomial of degree 12 with |a/2^k|^2 <= 1/16, then squared k times
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

template link*(m: untyped; nu: static int; j: untyped) =
  ## m = U_nu at site j, local or from the halo, in kernels with the locals
  ## n, u, nb, ro, rs, h0..h3 of a GpuGauge
  let jj = j
  if jj < n:
    mload(m, u, 18*nu*n + lo18(V, jj), V)
  else:
    let p = jj - n
    let o = int ro[p]
    let st = int rs[p]
    when nu == 0: mload(m, h0, o, st)
    elif nu == 1: mload(m, h1, o, st)
    elif nu == 2: mload(m, h2, o, st)
    else: mload(m, h3, o, st)

proc actionA*[V: static int](g: var GpuGauge[V]; c: GaugeActionCoeffs): float =
  ## the gauge action plaq + adjplaq of g.u, as actionA
  g.update
  let n = g.n
  let u = g.u
  let nb = g.nb
  let ro = g.ex[0].rofs
  let rs = g.ex[0].rstr
  let h0 = g.ex[0].rbuf
  let h1 = g.ex[1].rbuf
  let h2 = g.ex[2].rbuf
  let h3 = g.ex[3].rbuf
  var pl = gpuSum(k, n, 2):
    var a, b = 0.0
    forStatic mu, 1, 3:
      forStatic nu, 0, mu-1:
        var x {.noInit.}, y {.noInit.}, t {.noInit.}, p {.noInit.}: array[18, float]
        mload(x, u, 18*mu*n + lo18(V, k), V)
        link(y, nu, int nb[nbF(mu)*n + k])
        mmul(t, x, y)
        link(y, mu, int nb[nbF(nu)*n + k])
        mmulNA(p, t, y)
        mload(y, u, 18*nu*n + lo18(V, k), V)
        # tr P = tr(p y^+) = sum_ab p_ab conj(y_ab)
        var tr, ti = 0.0
        forStatic e, 0, 8:
          tr += p[2*e]*y[2*e] + p[2*e+1]*y[2*e+1]
          ti += p[2*e+1]*y[2*e] - p[2*e]*y[2*e+1]
        a += tr
        b += tr*tr + ti*ti
    [a, b]
  pl[0] /= 3.0
  pl[1] /= 9.0
  getDefaultComm().allReduce(addr pl[0], 2)
  let a0 = 6.0*float(g.lo.physVol)
  c.plaq*(a0 - pl[0]) + c.adjplaq*(a0 - pl[1])

proc plaq*[V: static int](g: var GpuGauge[V]): float =
  ## average Re tr P/3
  let c = GaugeActionCoeffs(plaq: 1.0)
  1.0 - g.actionA(c)/(6.0*float(g.lo.physVol))

proc forceA*[V: static int](g: var GpuGauge[V]; c: GaugeActionCoeffs; p: ptr UncheckedArray[float]; t: float) =
  ## p -= t F for the gauge force F of plaq + adjplaq, as forceA.  With a
  ## stale halo, the kernel over the sites with local neighbors also stores
  ## the send slots of the 4 halos, as update, and the kernel over the
  ## other sites follows the exchange.
  tic("gauge force")
  let n = g.n
  let u = g.u
  let nb = g.nb
  let od = g.ord
  let nin = g.nin
  let xch = not g.fresh
  let ns = if xch: g.ex[0].nsend else: 0  # send slots, the same in the 4 halos
  let si = g.ex[0].sidx
  let st = g.ex[0].sstr
  let sd0 = g.ex[0].sdst
  let sd1 = g.ex[1].sdst
  let sd2 = g.ex[2].sdst
  let sd3 = g.ex[3].sdst
  let ro = g.ex[0].rofs
  let rs = g.ex[0].rstr
  let h0 = g.ex[0].rbuf
  let h1 = g.ex[1].rbuf
  let h2 = g.ex[2].rbuf
  let h3 = g.ex[3].rbuf
  let cp = c.plaq/3.0
  let ca = 2.0*c.adjplaq/9.0
  template body(i, m, j0: untyped) =
    ## the link of direction i div m at site od[j0 + i mod m]
    let mu = i div m
    let k = int od[j0 + i - mu*m]
    var x {.noInit.}, y {.noInit.}, s {.noInit.}, q {.noInit.}, f {.noInit.}: array[18, float]
    forStatic e, 0, 17: f[e] = 0.0
    template plaqs(mu: static int) =
      mload(x, u, 18*mu*n + lo18(V, k), V)
      forStatic nu, 0, 3:
        when nu != mu:
          # U_mu(x) U_nu(x+mu) U_mu(x+nu)^+ U_nu(x)^+
          link(y, nu, int nb[nbF(mu)*n + k])
          mmul(s, x, y)
          link(y, mu, int nb[nbF(nu)*n + k])
          mmulNA(q, s, y)
          mload(y, u, 18*nu*n + lo18(V, k), V)
          mmulNA(s, q, y)
          acc(s)
          # U_mu(x) U_nu(x+mu-nu)^+ U_mu(x-nu)^+ U_nu(x-nu)
          link(y, nu, int nb[nbD(mu,nu)*n + k])
          mmulNA(s, x, y)
          link(y, mu, int nb[nbB(nu)*n + k])
          mmulNA(q, s, y)
          link(y, nu, int nb[nbB(nu)*n + k])
          mmul(s, q, y)
          acc(s)
    template acc(s: untyped) =
      ## f += conj(cp + ca tr s) s
      let cr = cp + ca*(s[0] + s[8] + s[16])
      let ci = -ca*(s[1] + s[9] + s[17])
      forStatic e, 0, 8:
        f[2*e] += cr*s[2*e] - ci*s[2*e+1]
        f[2*e+1] += cr*s[2*e+1] + ci*s[2*e]
    case mu
    of 0: plaqs(0)
    of 1: plaqs(1)
    of 2: plaqs(2)
    else: plaqs(3)
    mtah(s, f)
    let o = 18*mu*n + lo18(V, k)
    forStatic e, 0, 17: p[o + e*V] -= t*s[e]
  if xch: getDefaultComm().barrier  # as in update
  let mi = max(4*nin, 1)
  gpuFor(i, mi, 16):  # SIMD32 spills, 1.5-1.8x slower
    var q = i
    while q < 4*18*ns:  # the send slots, spread over the threads
      let mu = q div (18*ns)
      let r = q - mu*18*ns
      let f = cast[ptr UncheckedArray[float]](addr u[18*mu*n])
      case mu
      of 0: packAt(si, sd0, st, 18, ns, V, r, f)
      of 1: packAt(si, sd1, st, 18, ns, V, r, f)
      of 2: packAt(si, sd2, st, 18, ns, V, r, f)
      else: packAt(si, sd3, st, 18, ns, V, r, f)
      q += mi
    if i < 4*nin: body(i, nin, 0)
  toc("interior")
  if xch:
    for mu in 0..3: g.ex[mu].start
    for mu in 0..3: g.ex[mu].wait
    g.fresh = true
  toc("halo")
  if nin < n:
    gpuFor(i, 4*(n-nin), 16): body(i, n-nin, nin)
  toc("boundary")

proc expUpdate*[V: static int](g: var GpuGauge[V]; p: ptr UncheckedArray[float]; t: float) =
  ## g.u = exp(t p) g.u link by link, as axexpmuly
  tic("gauge exp")
  let n = g.n
  let u = g.u
  gpuFor(i, 4*n):
    let mu = i div n
    let k = i - mu*n
    let o = 18*mu*n + lo18(V, k)
    var a {.noInit.}, e {.noInit.}, x {.noInit.}, y {.noInit.}: array[18, float]
    forStatic j, 0, 17: a[j] = t*p[o + j*V]
    mexp(e, a)
    mload(x, u, o, V)
    mmul(y, e, x)
    forStatic j, 0, 17: u[o + j*V] = y[j]
  g.fresh = false
  toc("done")

proc norm2*[V: static int](g: GpuGauge[V]; p: ptr UncheckedArray[float]): float =
  ## global sum of |p|^2 over the links, p like g.u
  result = gpuSum(i, 4*18*g.n, 1, [p[i]*p[i]])[0]
  getDefaultComm().allReduce(result)
