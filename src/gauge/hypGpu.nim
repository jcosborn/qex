## HYP smearing of the links of a GpuGauge on the GPU, as gauge/hypsmear2.
##
## With a_1 = alpha1/2, a_2 = alpha2/4, a_3 = alpha3/6 and P the projection
## X (X^+ X)^-1/2 of projectU:
##   V1_mu,nu = P[(1-alpha1) U_mu + a_1 S_nu(U_nu, U_mu)]
##   V2_mu,nu = P[(1-alpha2) U_mu + a_2 sum_{a != mu,nu} S_a(V1_a,b, V1_mu,b)],  b != mu,nu,a
##   V3_mu = P[(1-alpha3) U_mu + a_3 sum_{nu != mu} S_nu(V2_nu,mu, V2_mu,nu)]
##   S_nu(X, Y)(x) = X(x) Y(x+nu) X(x+mu)^+ + X(x-nu)^+ Y(x-nu) X(x-nu+mu)
## Levels 1 and 2 run on the box sites: the local sites, then the sites of the
## width 1 shell around them in the dimensions split across ranks, whose thin
## links come in one exchange; level 3 runs on the local sites.  As in
## hypsmear2, a box site computes a level only where all its staples lie in
## the box, which covers every value the local level 3 reads.  The force runs
## the chain rule back through the levels as force of hypsmear2 and ends with
## one reverse exchange adding the force on the shell sites to their owners.
##
## Box site k < n is local site k of the fields, n + p the shell site at
## receive position p of the exchanges.  Box fields are [f][nb div V][18][V]
## like GpuGauge.u, the pairs mu != nu at f = ix(mu, nu).  A thread takes one
## pair of one site, in sub-groups of 16 (SIMD32 is 3-12% slower), which the
## kernels take by tiles of the local sites (ord), so that the links a
## sub-group reads at its neighbors, of either parity, are still in the L2
## cache when the neighbors read them.  The force computes the staple sums
## of level 2 and the derivative of its P in one kernel, those of level 1,
## with the staples of two directions, in two kernels, which is faster.
import qex
import gauge/[hypsmear, gaugeGpu]
import comms/[gather, halogpu]
import backend/accel
import base/metaUtils
import std/[algorithm, math]
getOptimPragmas()

type
  HypGpu*[V: static int] = object
    n*, nr*, nb*: int  # local sites, shell sites, box sites n + nr rounded up to V
    nbr*: ptr UncheckedArray[int32]  # [8][nb]: box index of x+mu, then of x-mu, -1 outside the box
    ord*: ptr UncheckedArray[int32]  # the sub-groups of 16 box sites in the order of the kernels, local ones first
    ex*: array[4, GpuHaloEx[float]]  # thin links of the shell sites
    rv*: array[4, GpuHaloEx[float]]  # force on the shell sites, to their owners
    v1x*, v1*, v2x*, v2*: ptr UncheckedArray[float]  # [12][nb div V][18][V]: levels 1 and 2, before and after P
    v3x*: ptr UncheckedArray[float]  # [4][n div V][18][V]: level 3 before P
    cf*, c3*: ptr UncheckedArray[float]  # [4][nb div V][18][V]: force on the thin links, derivative of P at level 3
    c1*, c2*: ptr UncheckedArray[float]  # [12][nb div V][18][V]: derivatives of P at levels 1 and 2

  M3 = array[18, float]  # 3x3 complex matrix, U_ab at 6a+2b (re) and 6a+2b+1 (im)

template ix(mu, nu: untyped): untyped = 3*mu + (if nu < mu: nu else: nu-1)

template bo(V, f, nb, k: untyped): untyped =
  ## real 0 of box site k of field f, real e is e*V further
  18*f*nb + lo18(V, k)

proc newHypGpu*[V: static int](lo: Layout[V]): HypGpu[V] =
  tic("newHypGpu")
  let nd = lo.nDim
  let n = lo.nSites
  let c = getDefaultComm()
  var b0, be = newSeq[int](nd)  # box origin and extent in local coordinates
  var nbox = 1
  for d in 0..<nd:
    let sp = lo.rankGeom[d] > 1
    b0[d] = if sp: -1 else: 0
    be[d] = lo.localGeom[d] + (if sp: 2 else: 0)
    nbox *= be[d]
  var bi = newSeq[int32](nbox)  # box index of each box site, lexicographic in be
  var rl = newSeq[RecvList]()
  var x, y = newSeq[int](nd)
  var nr = 0
  for b in 0..<nbox:
    var t = b
    var loc = true
    for d in 0..<nd:
      x[d] = b0[d] + t mod be[d]
      t = t div be[d]
      loc = loc and x[d] >= 0 and x[d] < lo.localGeom[d]
      y[d] = (lo.coordmin[d] + x[d] + lo.physGeom[d]) mod lo.physGeom[d]
    let ri = lo.rankIndex(y)
    if loc:
      bi[b] = int32 ri.index
    else:
      rl.add RecvList(didx: int32 nr, srank: int32 ri.rank, sidx: int32 ri.index)
      bi[b] = int32(-1 - nr)
      inc nr
  toc("box")
  let gm = c.makeGatherMap(rl)
  if gm.lidx.len > 0: qexError("newHypGpu: shell sites on this rank")
  var pos = newSeq[int32](nr)
  for p, d in gm.rdest: pos[d] = int32 p
  for b in 0..<nbox:
    if bi[b] < 0: bi[b] = int32(n) + pos[-1 - bi[b]]
  let nb = (n + nr + V - 1) div V * V
  result.n = n
  result.nr = nr
  result.nb = nb
  var nbr = newSeq[int32](2*nd*nb)
  for k in 0..<nbr.len: nbr[k] = -1
  for b in 0..<nbox:
    var t = b
    for d in 0..<nd:
      x[d] = t mod be[d]
      t = t div be[d]
    for mu in 0..<nd:
      for fb in 0..1:
        for d in 0..<nd: y[d] = x[d]
        y[mu] += 1 - 2*fb
        if lo.rankGeom[mu] == 1: y[mu] = (y[mu] + be[mu]) mod be[mu]
        if y[mu] >= 0 and y[mu] < be[mu]:
          var l = 0
          for d in countdown(nd-1, 0): l = l*be[d] + y[d]
          nbr[(fb*nd + mu)*nb + bi[b]] = bi[l]
  result.nbr = nbr.toDevice
  # the local sub-groups by tiles of 4 outer sites in each dimension but the
  # last, whole in it, lexicographic in a tile, then the shell ones: the
  # kernels read the links at the neighbors of a site again while they are
  # in the L2 cache, also those of the other parity, which qlayout keeps in
  # the other half
  let og = lo.outerGeom
  var ks = newSeq[(int, int, int)]((n + nr + 15) div 16)
  for g in 0..<ks.len:
    var t, w = 0
    if 16*g < n:
      for d in countdown(nd-1, 0):
        let e = if d == nd-1: og[d] else: min(4, og[d])
        let x = (lo.coords[d][16*g] - lo.coordmin[d]) mod og[d]
        t = t*((og[d] + e - 1) div e) + x div e
        w = w*e + x mod e
    else:
      t = int.high
    ks[g] = (t, w, g)
  ks.sort
  var ord = newSeq[int32](ks.len)
  for b, s in ks: ord[b] = int32 s[2]
  result.ord = ord.toDevice
  toc("neighbors")
  for mu in 0..<nd: result.ex[mu] = newGpuHaloEx[float](gm, 18, n, V, c)
  var gr = gm.reverse  # sends the shell sites, n + p at position p, back
  for p in 0..<nr: gr.sidx[p] = int32(n + p)
  for mu in 0..<nd: result.rv[mu] = newGpuHaloEx[float](gr, 18, nb, V, c)
  toc("exchanges")
  proc zeros(m: int): ptr UncheckedArray[float] =
    result = cast[ptr UncheckedArray[float]](gpuMalloc(m*sizeof(float)))
    let r = result
    gpuFor(i, m): r[i] = 0.0
  let m = 18*nb
  result.v1x = zeros(12*m)
  result.v1 = zeros(12*m)
  result.v2x = zeros(12*m)
  result.v2 = zeros(12*m)
  result.v3x = zeros(4*18*n)
  result.cf = zeros(4*m)
  result.c3 = zeros(4*m)
  result.c1 = zeros(12*m)
  result.c2 = zeros(12*m)
  toc("fields")

proc free*[V: static int](h: var HypGpu[V]) =
  gpuFree(h.nbr)
  gpuFree(h.ord)
  for p in [h.v1x, h.v1, h.v2x, h.v2, h.v3x, h.cf, h.c3, h.c1, h.c2]: gpuFree(p)
  for mu in 0..3:
    h.ex[mu].free
    h.rv[mu].free

# 3x3 complex matrices in kernels, loops unrolled so that the matrices stay
# in registers; the results may not alias the arguments

proc mul3(r: var M3; a, b: M3) {.alwaysInline.} =
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

proc mulNA3(r: var M3; a, b: M3) {.alwaysInline.} =
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

proc mulAN3(r: var M3; a, b: M3) {.alwaysInline.} =
  ## r = a^+ b
  forStatic i, 0, 2:
    forStatic j, 0, 2:
      var re = a[2*i]*b[2*j] + a[2*i+1]*b[2*j+1]
      var im = a[2*i]*b[2*j+1] - a[2*i+1]*b[2*j]
      forStatic k, 1, 2:
        re += a[6*k+2*i]*b[6*k+2*j] + a[6*k+2*i+1]*b[6*k+2*j+1]
        im += a[6*k+2*i]*b[6*k+2*j+1] - a[6*k+2*i+1]*b[6*k+2*j]
      r[6*i+2*j] = re
      r[6*i+2*j+1] = im

proc addc3(r: var M3; cr, ci: float; a: M3) {.alwaysInline.} =
  ## r += (cr + i ci) a
  forStatic e, 0, 8:
    r[2*e] += cr*a[2*e] - ci*a[2*e+1]
    r[2*e+1] += cr*a[2*e+1] + ci*a[2*e]

template adjEl(r, x: untyped; i, j, a, b, c, d, e, f, g, h: static int) =
  ## r_ij = x_ab x_cd - x_ef x_gh
  r[6*i+2*j] = x[6*a+2*b]*x[6*c+2*d] - x[6*a+2*b+1]*x[6*c+2*d+1] - x[6*e+2*f]*x[6*g+2*h] + x[6*e+2*f+1]*x[6*g+2*h+1]
  r[6*i+2*j+1] = x[6*a+2*b]*x[6*c+2*d+1] + x[6*a+2*b+1]*x[6*c+2*d] - x[6*e+2*f]*x[6*g+2*h+1] - x[6*e+2*f+1]*x[6*g+2*h]

proc adj3(r: var M3; x: M3) {.alwaysInline.} =
  ## r = the adjugate of x, as adjugate
  adjEl(r, x, 0, 0, 1, 1, 2, 2, 1, 2, 2, 1)
  adjEl(r, x, 0, 1, 2, 1, 0, 2, 2, 2, 0, 1)
  adjEl(r, x, 0, 2, 0, 1, 1, 2, 0, 2, 1, 1)
  adjEl(r, x, 1, 0, 1, 2, 2, 0, 1, 0, 2, 2)
  adjEl(r, x, 1, 1, 2, 2, 0, 0, 2, 0, 0, 2)
  adjEl(r, x, 1, 2, 0, 2, 1, 0, 0, 0, 1, 2)
  adjEl(r, x, 2, 0, 1, 0, 2, 1, 1, 1, 2, 0)
  adjEl(r, x, 2, 1, 2, 0, 0, 1, 2, 1, 0, 0)
  adjEl(r, x, 2, 2, 0, 0, 1, 1, 0, 1, 1, 0)

proc det3(dr, di: var float; x, a: M3) {.alwaysInline.} =
  ## dr + i di = det x = x_20 a_02 + x_21 a_12 + x_22 a_22 for a = adj x, as determinant
  dr = x[12]*a[4] - x[13]*a[5] + x[14]*a[10] - x[15]*a[11] + x[16]*a[16] - x[17]*a[17]
  di = x[12]*a[5] + x[13]*a[4] + x[14]*a[11] + x[15]*a[10] + x[16]*a[17] + x[17]*a[16]

proc inv3(r: var M3; x: M3) {.alwaysInline.} =
  ## r = x^-1 = adj x/det x, as inverse
  var dr, di: float
  adj3(r, x)
  det3(dr, di, x, r)
  let dn = 1.0/(dr*dr + di*di)
  let ir = dr*dn
  let ii = -di*dn
  forStatic e, 0, 8:
    let a = r[2*e]
    let b = r[2*e+1]
    r[2*e] = ir*a - ii*b
    r[2*e+1] = ir*b + ii*a

proc rsqrt3(z: var M3; x: M3) {.alwaysInline.} =
  ## z = (x^+ x + 1e-20)^-1/2, as projectUrsqrt: with t = x^+ x + 1e-20, its
  ## eigenvalues l_k from tr t, tr t^2 and det t as eigs3, and
  ## z = c0 + c1 t + c2 t^2 as rsqrtPHM3f
  var t {.noInit.}, t2 {.noInit.}: M3
  var det, deti: float
  mulAN3(t, x, x)
  forStatic i, 0, 2: t[8*i] += 1e-20
  mul3(t2, t, t)
  adj3(z, t)
  det3(det, deti, t, z)
  let tr = t[0] + t[8] + t[16]
  let p2 = t2[0] + t2[8] + t2[16]
  let tr3 = (1.0/3.0)*tr
  let p23 = (1.0/3.0)*p2
  let tr32 = tr3*tr3
  let q = abs(0.5*(p23-tr32))
  let r = 0.25*tr3*(5*tr32-p2) - 0.5*det
  let sq = sqrt(q)
  let sq3 = q*sq
  let isq3 = 1.0/max(sq3, 1.0/3e38)  # eigs3 clamps 1/sq3 to 3e38
  let rsq3 = min(1.0, max(-1.0, r*isq3))
  let th = (1.0/3.0)*arccos(rsq3)
  let st = sin(th)
  let ct = cos(th)
  let sqc = sq*ct
  let sqs = 1.73205080756887729352*sq*st  # sqrt(3)
  let ll = tr3 + sqc
  let l0 = tr3 - 2*sqc
  let l1 = ll + sqs
  let l2 = ll - sqs
  let sl0 = sqrt(abs(l0))
  let sl1 = sqrt(abs(l1))
  let sl2 = sqrt(abs(l2))
  let u = sl0 + sl1 + sl2
  let w = sl0 * sl1 * sl2
  let di = 1.0/(w*(sl0+sl1)*(sl0+sl2)*(sl1+sl2))
  let c0 = (w*u*u+l0*sl0*(l1+l2)+l1*sl1*(l0+l2)+l2*sl2*(l0+l1))*di
  let c1 = -(tr*u+w)*di
  let c2 = u*di
  forStatic e, 0, 17: z[e] = c1*t[e] + c2*t2[e]
  forStatic i, 0, 2: z[8*i] += c0

proc projU3(r: var M3; x: M3) {.alwaysInline.} =
  ## r = x (x^+ x)^-1/2, as projectU
  var z {.noInit.}: M3
  rsqrt3(z, x)
  mul3(r, x, z)

proc usdev(r: var float; x: M3) {.alwaysInline.} =
  ## r = |x^+ x - 1|^2 + |det x - 1|^2, as checkSU of a matrix
  var t {.noInit.}, a {.noInit.}: M3
  var dr, di: float
  mulAN3(t, x, x)
  forStatic i, 0, 2: t[8*i] -= 1.0
  r = 0.0
  forStatic e, 0, 17: r += t[e]*t[e]
  adj3(a, x)
  det3(dr, di, x, a)
  r += (dr - 1.0)*(dr - 1.0) + di*di

proc sylsolve3(x: var M3; a, c: M3) {.alwaysInline.} =
  ## x with a x + x a = c, as sylsolve: for d = adj a, t = tr a, s = tr d, r = det a,
  ##   x = c0 c - c4 (a c + c a) + c2 a c a + c1 d c d - c2 (d c + c d)
  ##   c2 = 1/(2 (s t - r)),  c0 = c2 (s + t^2),  c1 = c2 t/r,  c4 = c2 t
  ## d is computed again for its terms, which keeps five matrices live
  var d {.noInit.}, t {.noInit.}, w {.noInit.}: M3
  adj3(d, a)
  let tr = a[0] + a[8] + a[16]
  let ti = a[1] + a[9] + a[17]
  let sr = d[0] + d[8] + d[16]
  let si = d[1] + d[9] + d[17]
  let rr = a[0]*d[0] - a[1]*d[1] + a[2]*d[6] - a[3]*d[7] + a[4]*d[12] - a[5]*d[13]
  let ri = a[0]*d[1] + a[1]*d[0] + a[2]*d[7] + a[3]*d[6] + a[4]*d[13] + a[5]*d[12]
  let er = 2.0*(sr*tr - si*ti - rr)
  let ei = 2.0*(sr*ti + si*tr - ri)
  let en = 1.0/(er*er + ei*ei)
  let c2r = er*en
  let c2i = -ei*en
  let qr = sr + tr*tr - ti*ti
  let qi = si + 2.0*tr*ti
  let rn = 1.0/(rr*rr + ri*ri)
  let tor = (tr*rr + ti*ri)*rn
  let toi = (ti*rr - tr*ri)*rn
  let c4r = c2r*tr - c2i*ti
  let c4i = c2r*ti + c2i*tr
  forStatic e, 0, 17: x[e] = 0.0
  addc3(x, c2r*qr - c2i*qi, c2r*qi + c2i*qr, c)  # c0 c
  mul3(t, a, c)
  addc3(x, -c4r, -c4i, t)
  mul3(w, t, a)
  addc3(x, c2r, c2i, w)
  mul3(t, c, a)
  addc3(x, -c4r, -c4i, t)
  adj3(d, a)
  mul3(t, d, c)
  addc3(x, -c2r, -c2i, t)
  mul3(w, t, d)
  addc3(x, c2r*tor - c2i*toi, c2r*toi + c2i*tor, w)  # c1 d c d
  mul3(t, c, d)
  addc3(x, -c2r, -c2i, t)

template projUderiv3(r: untyped; lu, lx, lc, park, unpark: untyped) =
  ## r = the derivative of projectU at x for the chain c, u = projectU(x),
  ## as projectUderiv: z = (x^+ x)^-1/2, y = z^-1, s y + y s = u^+ c z,
  ## r = c z - x (s + s^+).  lu, lx, lc load u, x, c into a matrix, and
  ## r waits in memory through park and unpark while sylsolve3 runs, so
  ## that the matrices stay in registers.
  block:
    var z {.noInit.}, y {.noInit.}, t1 {.noInit.}, t2 {.noInit.}: M3
    lx(t1)
    rsqrt3(z, t1)
    inv3(y, z)
    lc(t1)
    mul3(r, t1, z)
    lu(t1)
    mulAN3(t2, t1, r)
    park(r)
    sylsolve3(t1, y, t2)
    forStatic i, 0, 2:
      forStatic j, 0, 2:
        t2[6*i+2*j] = t1[6*i+2*j] + t1[6*j+2*i]
        t2[6*i+2*j+1] = t1[6*i+2*j+1] - t1[6*j+2*i+1]
    lx(y)
    mul3(t1, y, t2)
    unpark(r)
    forStatic e, 0, 17: r[e] -= t1[e]

template planeSite(i: untyped; np: static int; ord, q, k: untyped) =
  ## plane q of np and site k of thread i: the planes of the 16 consecutive
  ## sites of sub-group ord[b] in consecutive sub-groups of 16, which then
  ## run close in time and share the loads of their neighbors in the caches,
  ## each sub-group still loading 16 consecutive sites
  let b = i div (16*np)
  let r = i - b*(16*np)
  let q = r div 16
  let k = 16*int(ord[b]) + r - q*16

template shellLinks(h, g: untyped) =
  ## the locals of thin, for the links of g at box sites
  let n {.inject.} = h.n
  let u {.inject.} = g.u
  let ro {.inject.} = h.ex[0].rofs
  let rs {.inject.} = h.ex[0].rstr
  let h0 {.inject.} = h.ex[0].rbuf
  let h1 {.inject.} = h.ex[1].rbuf
  let h2 {.inject.} = h.ex[2].rbuf
  let h3 {.inject.} = h.ex[3].rbuf

template thin(m: untyped; d, j: int) =
  ## m = the thin link of direction d at box site j, in kernels with the locals of shellLinks
  let jj = j
  if jj < n:
    mload(m, u, 18*d*n + lo18(V, jj), V)
  else:
    let p = jj - n
    let hb = if d == 0: h0 elif d == 1: h1 elif d == 2: h2 else: h3
    mload(m, hb, int ro[p], int rs[p])

template staple(s: untyped; al: float; lx, ly: untyped; i, fnu, fmu, bnu, fmubnu: int) =
  ## s += al S(X, Y) with lx(m, j), ly(m, j) loading X(j), Y(j) into m:
  ## X(i) Y(fnu) X(fmu)^+ + X(bnu)^+ Y(bnu) X(fmubnu)
  block:
    var a {.noInit.}, b {.noInit.}, t {.noInit.}: M3
    lx(a, i)
    ly(b, fnu)
    mul3(t, a, b)
    lx(a, fmu)
    mulNA3(b, t, a)
    forStatic e, 0, 17: s[e] += al*b[e]
    lx(a, bnu)
    ly(b, bnu)
    mulAN3(t, a, b)
    lx(a, fmubnu)
    mul3(b, t, a)
    forStatic e, 0, 17: s[e] += al*b[e]

template symderiv(hit, s: untyped; lx, ly, cx, cy: untyped; i, fnu, fmu, bnu, fmubnu, nl: int) =
  ## s += the derivative of the staples S(X, Y) with the chains cx, cy, as
  ## symderiv3 of hypsmear2, for lx, ly, cx, cy loading into m at j as in
  ## staple: the terms with a chain site j < nl, all for nl = int.high.
  ## hit = true if any staple is there.  A matrix is loaded again rather
  ## than kept, as registers are short.
  block:
    var a {.noInit.}, b {.noInit.}, d {.noInit.}, t {.noInit.}: M3
    if fnu >= 0 and fmu >= 0 and (i < nl or fnu < nl or fmu < nl):
      hit = true
      if i < nl:  # cx(i) Y(fnu) X(fmu)^+
        cx(a, i)
        ly(b, fnu)
        mul3(t, a, b)
        lx(a, fmu)
        mulNA3(d, t, a)
        forStatic e, 0, 17: s[e] += d[e]
      if fnu < nl:  # X(i) cy(fnu) X(fmu)^+
        lx(a, i)
        cy(b, fnu)
        mul3(t, a, b)
        lx(a, fmu)
        mulNA3(d, t, a)
        forStatic e, 0, 17: s[e] += d[e]
      if fmu < nl:  # X(i) Y(fnu) cx(fmu)^+
        lx(a, i)
        ly(b, fnu)
        mul3(t, a, b)
        cx(a, fmu)
        mulNA3(d, t, a)
        forStatic e, 0, 17: s[e] += d[e]
    if bnu >= 0 and fmubnu >= 0 and (bnu < nl or fmubnu < nl):
      hit = true
      if bnu < nl:  # (cx(bnu)^+ Y(bnu) + X(bnu)^+ cy(bnu)) X(fmubnu)
        cx(a, bnu)
        ly(b, bnu)
        mulAN3(t, a, b)
        lx(a, bnu)
        cy(b, bnu)
        mulAN3(d, a, b)
        forStatic e, 0, 17: t[e] += d[e]
        lx(a, fmubnu)
        mul3(d, t, a)
        forStatic e, 0, 17: s[e] += d[e]
      if fmubnu < nl:  # X(bnu)^+ Y(bnu) cx(fmubnu)
        lx(a, bnu)
        ly(b, bnu)
        mulAN3(t, a, b)
        cx(a, fmubnu)
        mul3(d, t, a)
        forStatic e, 0, 17: s[e] += d[e]

proc smear*[V: static int](h: HypGpu[V]; c: HypCoefs; g: GpuGauge[V]; fl: ptr UncheckedArray[float]) =
  ## fl = the smeared links of g.u, like g.u
  tic("hyp smear")
  getDefaultComm().barrier  # peers may still read the previous shell links
  for mu in 0..3: h.ex[mu].pack(cast[ptr UncheckedArray[float]](addr g.u[18*mu*g.n]))
  for mu in 0..3: h.ex[mu].start
  for mu in 0..3: h.ex[mu].wait
  toc("exchange")
  shellLinks(h, g)
  let nb = h.nb
  let nt = n + h.nr
  let ntp = 16*((nt + 15) div 16)  # nt padded to whole blocks of planeSite
  let nbr = h.nbr
  let ord = h.ord
  let v1x = h.v1x
  let v1 = h.v1
  let v2x = h.v2x
  let v2 = h.v2
  let v3x = h.v3x
  let a1 = c.alpha1/2
  let a2 = c.alpha2/4
  let a3 = c.alpha3/6
  let m1 = 1 - c.alpha1
  let m2 = 1 - c.alpha2
  let m3 = 1 - c.alpha3
  gpuFor(i, 12*ntp, 16):
    planeSite(i, 12, ord, q, k)
    if k < nt:
      let mu = q div 3  # q = ix(mu, nu)
      let nu = q - 3*mu + int(q - 3*mu >= mu)
      let fmu = int nbr[mu*nb + k]
      let fnu = int nbr[nu*nb + k]
      let bnu = int nbr[(4+nu)*nb + k]
      let fmubnu = if fmu >= 0 and bnu >= 0: int nbr[(4+nu)*nb + fmu] else: -1
      if fnu >= 0 and fmubnu >= 0:
        var s {.noInit.}, r {.noInit.}: M3
        thin(r, mu, k)
        forStatic e, 0, 17: s[e] = m1*r[e]
        template lx(m, j: untyped) = thin(m, nu, j)
        template ly(m, j: untyped) = thin(m, mu, j)
        staple(s, a1, lx, ly, k, fnu, fmu, bnu, fmubnu)
        let o = bo(V, q, nb, k)
        forStatic e, 0, 17: v1x[o + e*V] = s[e]
        projU3(r, s)
        forStatic e, 0, 17: v1[o + e*V] = r[e]
  toc("level 1")
  gpuFor(i, 12*ntp, 16):
    planeSite(i, 12, ord, q, k)
    if k < nt:
      let mu = q div 3  # q = ix(mu, nu)
      let nu = q - 3*mu + int(q - 3*mu >= mu)
      let fmu = int nbr[mu*nb + k]
      if fmu >= 0:
        var s {.noInit.}, r {.noInit.}: M3
        thin(r, mu, k)
        forStatic e, 0, 17: s[e] = m2*r[e]
        for a in 0..3:
          if a != mu and a != nu:
            let b = 6 - mu - nu - a
            let fa = int nbr[a*nb + k]
            let ba = int nbr[(4+a)*nb + k]
            let fmuba = if ba >= 0: int nbr[(4+a)*nb + fmu] else: -1
            if fa >= 0 and fmuba >= 0:
              template lx(m, j: untyped) = mload(m, v1, bo(V, ix(a, b), nb, j), V)
              template ly(m, j: untyped) = mload(m, v1, bo(V, ix(mu, b), nb, j), V)
              staple(s, a2, lx, ly, k, fa, fmu, ba, fmuba)
        let o = bo(V, q, nb, k)
        forStatic e, 0, 17: v2x[o + e*V] = s[e]
        projU3(r, s)
        forStatic e, 0, 17: v2[o + e*V] = r[e]
  toc("level 2")
  gpuFor(i, 4*16*((n + 15) div 16), 16):  # the local sub-groups come first in ord
    planeSite(i, 4, ord, mu, k)
    if k < n:
      let fmu = int nbr[mu*nb + k]
      var w {.noInit.}, s {.noInit.}, r {.noInit.}: M3
      thin(w, mu, k)
      forStatic e, 0, 17: s[e] = m3*w[e]
      for nu in 0..3:
        if nu != mu:
          let fnu = int nbr[nu*nb + k]
          let bnu = int nbr[(4+nu)*nb + k]
          let fmubnu = int nbr[(4+nu)*nb + fmu]
          template lx(m, j: untyped) = mload(m, v2, bo(V, ix(nu, mu), nb, j), V)
          template ly(m, j: untyped) = mload(m, v2, bo(V, ix(mu, nu), nb, j), V)
          staple(s, a3, lx, ly, k, fnu, fmu, bnu, fmubnu)
      let o = 18*mu*n + lo18(V, k)
      forStatic e, 0, 17: v3x[o + e*V] = s[e]
      projU3(r, s)
      forStatic e, 0, 17: fl[o + e*V] = r[e]
  toc("level 3")

proc force*[V: static int](h: HypGpu[V]; c: HypCoefs; g: GpuGauge[V]; fl, f, sg: ptr UncheckedArray[float];
                           ne: int; p: ptr UncheckedArray[float]) =
  ## p += TAH(F U^+) link by link, F the force on the thin links U = g.u of the
  ## chain e sg f on the smeared links fl of smear, e = 1 on the even and -1
  ## on the odd sites: for f the one link forces of outerM and sg of
  ## stagSigns, smearedOneLinkForce of staghmc_sh with the pullback of
  ## hypsmear2.  Uses the levels of the last smear, of g.u.
  tic("hyp force")
  shellLinks(h, g)
  let nb = h.nb
  let nt = n + h.nr
  let ntp = 16*((nt + 15) div 16)  # nt padded to whole blocks of planeSite
  let nbr = h.nbr
  let ord = h.ord
  let v1x = h.v1x
  let v1 = h.v1
  let v2x = h.v2x
  let v2 = h.v2
  let v3x = h.v3x
  let cf = h.cf
  let c3 = h.c3
  let c2 = h.c2
  let c1 = h.c1
  let a1 = c.alpha1/2
  let a2 = c.alpha2/4
  let a3 = c.alpha3/6
  let m1 = 1 - c.alpha1
  let m2 = 1 - c.alpha2
  let m3 = 1 - c.alpha3
  # a level l stores the derivative r of its P in cl, which its readers
  # scale: the chain of the level below is a_l r, the force on the thin
  # links through the same pair (1-alpha_l) r
  gpuFor(i, 4*ntp, 16):  # level 3: c3 = r
    planeSite(i, 4, ord, mu, k)
    if k < nt:
      let o = bo(V, mu, nb, k)
      if k < n:
        let ol = 18*mu*n + lo18(V, k)
        let sc = if k < ne: sg[mu*n + k] else: -sg[mu*n + k]
        var r {.noInit.}: M3
        template lu(m: untyped) = mload(m, fl, ol, V)
        template lx(m: untyped) = mload(m, v3x, ol, V)
        template lc(m: untyped) =
          mload(m, f, ol, V)
          forStatic e, 0, 17: m[e] *= sc
        template park(m: untyped) =
          forStatic e, 0, 17: c3[o + e*V] = m[e]
        template unpark(m: untyped) = mload(m, c3, o, V)
        projUderiv3(r, lu, lx, lc, park, unpark)
        forStatic e, 0, 17: c3[o + e*V] = r[e]
      else:
        forStatic e, 0, 17: c3[o + e*V] = 0.0
  toc("level 3")
  gpuFor(i, 12*ntp, 16):  # level 2: c2 = r of the chain a3 s, s the staple sums, where a staple is there, else 0
    planeSite(i, 12, ord, q, k)
    if k < nt:
      let mu = q div 3  # q = ix(mu, nu)
      let nu = q - 3*mu + int(q - 3*mu >= mu)
      let fmu = int nbr[mu*nb + k]
      let fnu = int nbr[nu*nb + k]
      let bnu = int nbr[(4+nu)*nb + k]
      let fmubnu = if fmu >= 0: int nbr[(4+nu)*nb + fmu] else: -1
      var s {.noInit.}: M3
      forStatic e, 0, 17: s[e] = 0.0
      var hit = false
      block:
        template lx(m, j: untyped) = mload(m, v2, bo(V, ix(nu, mu), nb, j), V)
        template ly(m, j: untyped) = mload(m, v2, bo(V, q, nb, j), V)
        template cx(m, j: untyped) = mload(m, c3, bo(V, nu, nb, j), V)
        template cy(m, j: untyped) = mload(m, c3, bo(V, mu, nb, j), V)
        symderiv(hit, s, lx, ly, cx, cy, k, fnu, fmu, bnu, fmubnu, n)
      let o = bo(V, q, nb, k)
      if hit:
        var r {.noInit.}: M3
        template lu(m: untyped) = mload(m, v2, o, V)
        template lx(m: untyped) = mload(m, v2x, o, V)
        template lc(m: untyped) =
          forStatic e, 0, 17: m[e] = a3*s[e]
        template park(m: untyped) =
          forStatic e, 0, 17: c2[o + e*V] = m[e]
        template unpark(m: untyped) = mload(m, c2, o, V)
        projUderiv3(r, lu, lx, lc, park, unpark)
        forStatic e, 0, 17: c2[o + e*V] = r[e]
      else:
        forStatic e, 0, 17: c2[o + e*V] = 0.0
  toc("level 2")
  # the chain a2 s of level 1, s its staple sums, goes to c1, then a second
  # kernel replaces it with r
  gpuFor(i, 12*ntp, 16):  # level 1: c1 = a2 s
    planeSite(i, 12, ord, q, k)
    if k < nt:
      let mu = q div 3  # q = ix(mu, nu)
      let nu = q - 3*mu + int(q - 3*mu >= mu)
      let fmu = int nbr[mu*nb + k]
      let fnu = int nbr[nu*nb + k]
      let bnu = int nbr[(4+nu)*nb + k]
      var s {.noInit.}: M3
      forStatic e, 0, 17: s[e] = 0.0
      var hit = false
      if fmu >= 0 and fnu >= 0 and bnu >= 0 and nbr[(4+nu)*nb + fmu] >= 0:
        for a in 0..3:
          if a != mu and a != nu:
            let b = 6 - mu - nu - a
            let fa = int nbr[a*nb + k]
            let ba = int nbr[(4+a)*nb + k]
            let fmuba = int nbr[(4+a)*nb + fmu]
            template lx(m, j: untyped) = mload(m, v1, bo(V, ix(a, nu), nb, j), V)
            template ly(m, j: untyped) = mload(m, v1, bo(V, q, nb, j), V)
            template cx(m, j: untyped) = mload(m, c2, bo(V, ix(a, b), nb, j), V)
            template cy(m, j: untyped) = mload(m, c2, bo(V, ix(mu, b), nb, j), V)
            symderiv(hit, s, lx, ly, cx, cy, k, fa, fmu, ba, fmuba, int.high)
      let o = bo(V, q, nb, k)
      forStatic e, 0, 17: c1[o + e*V] = a2*s[e]
  toc("level 1 sums")
  gpuFor(i, 12*ntp, 16):  # level 1: c1 = r where a staple is there, elsewhere s = 0
    planeSite(i, 12, ord, q, k)
    if k < nt:
      let mu = q div 3
      let nu = q - 3*mu + int(q - 3*mu >= mu)
      let fmu = int nbr[mu*nb + k]
      let fnu = int nbr[nu*nb + k]
      let bnu = int nbr[(4+nu)*nb + k]
      var hit = false
      if fmu >= 0 and fnu >= 0 and bnu >= 0 and nbr[(4+nu)*nb + fmu] >= 0:
        for a in 0..3:
          if a != mu and a != nu:
            hit = hit or nbr[a*nb + k] >= 0 or nbr[(4+a)*nb + k] >= 0 and nbr[(4+a)*nb + fmu] >= 0
      if hit:
        let o = bo(V, q, nb, k)
        var r {.noInit.}: M3
        template lu(m: untyped) = mload(m, v1, o, V)
        template lx(m: untyped) = mload(m, v1x, o, V)
        template lc(m: untyped) = mload(m, c1, o, V)
        template park(m: untyped) =
          forStatic e, 0, 17: c1[o + e*V] = m[e]
        template unpark(m: untyped) = mload(m, c1, o, V)
        projUderiv3(r, lu, lx, lc, park, unpark)
        forStatic e, 0, 17: c1[o + e*V] = r[e]
  toc("level 1")
  gpuFor(i, 4*ntp, 16):  # thin links: cf = a1 s + m3 c3 + sum_nu m2 c2 + m1 c1
    planeSite(i, 4, ord, mu, k)
    if k < nt:
      let fmu = int nbr[mu*nb + k]
      var s {.noInit.}, w {.noInit.}: M3
      forStatic e, 0, 17: s[e] = 0.0
      for nu in 0..3:
        if nu != mu:
          let fnu = int nbr[nu*nb + k]
          let bnu = int nbr[(4+nu)*nb + k]
          let fmubnu = if fmu >= 0: int nbr[(4+nu)*nb + fmu] else: -1
          template lx(m, j: untyped) = thin(m, nu, j)
          template ly(m, j: untyped) = thin(m, mu, j)
          template cx(m, j: untyped) = mload(m, c1, bo(V, ix(nu, mu), nb, j), V)
          template cy(m, j: untyped) = mload(m, c1, bo(V, ix(mu, nu), nb, j), V)
          var hit = false
          symderiv(hit, s, lx, ly, cx, cy, k, fnu, fmu, bnu, fmubnu, int.high)
      let o = bo(V, mu, nb, k)
      mload(w, c3, o, V)
      forStatic e, 0, 17: s[e] = a1*s[e] + m3*w[e]
      for nu in 0..3:
        if nu != mu:
          let oq = bo(V, ix(mu, nu), nb, k)
          mload(w, c2, oq, V)
          forStatic e, 0, 17: s[e] += m2*w[e]
          mload(w, c1, oq, V)
          forStatic e, 0, 17: s[e] += m1*w[e]
      forStatic e, 0, 17: cf[o + e*V] = s[e]
  toc("thin links")
  getDefaultComm().barrier  # peers may still read the previous force
  for mu in 0..3: h.rv[mu].pack(cast[ptr UncheckedArray[float]](addr cf[18*mu*nb]))
  for mu in 0..3: h.rv[mu].start
  for mu in 0..3: h.rv[mu].wait
  toc("exchange")
  let sl = h.ex[0].sslot
  let nsl = h.ex[0].nslot
  let qo = h.rv[0].rofs
  let qs = h.rv[0].rstr
  let r0 = h.rv[0].rbuf
  let r1 = h.rv[1].rbuf
  let r2 = h.rv[2].rbuf
  let r3 = h.rv[3].rbuf
  gpuFor(i, 4*n, 16):  # the force of the shell copies, then p += TAH(F U^+)
    let mu = i div n
    let k = i - mu*n
    var s {.noInit.}, w {.noInit.}, r {.noInit.}: M3
    mload(s, cf, bo(V, mu, nb, k), V)
    for q in 0..<nsl:
      let j = int sl[q*n + k]
      if j >= 0:
        let rb = if mu == 0: r0 elif mu == 1: r1 elif mu == 2: r2 else: r3
        mload(w, rb, int qo[j], int qs[j])
        forStatic e, 0, 17: s[e] += w[e]
    let ol = 18*mu*n + lo18(V, k)
    mload(w, u, ol, V)
    mulNA3(r, s, w)
    mtah(w, r)
    forStatic e, 0, 17: p[ol + e*V] += w[e]
  toc("combine")

proc reunit*[V: static int](g: var GpuGauge[V]): array[2, tuple[avg, max: float]] =
  ## g.u = W exp(-i arg(det W)/3) for W = U (U^+U)^-1/2, as projectSU of the
  ## host links, and checkSU of the links before and after: with d of
  ## usdev per link, avg = sqrt(sum d/(20 links)), max = sqrt(max d/20)
  const nm = 128  # links per partial maximum
  let n = g.n
  let u = g.u
  let nl = 4*n
  let np = (nl + nm - 1) div nm
  let d = cast[ptr UncheckedArray[float]](gpuMalloc(2*nl*sizeof(float)))
  let dp = cast[ptr UncheckedArray[float]](gpuMalloc(2*np*sizeof(float)))
  gpuFor(i, nl):
    let mu = i div n
    let k = i - mu*n
    let o = 18*mu*n + lo18(V, k)
    var x {.noInit.}, w {.noInit.}, a {.noInit.}: M3
    var dr, di: float
    forStatic e, 0, 17: x[e] = u[o + e*V]
    usdev(d[i], x)
    projU3(w, x)
    adj3(a, w)
    det3(dr, di, w, a)
    let p = -(1.0/3.0)*arctan2(di, dr)
    let cr = cos(p)
    let ci = sin(p)
    forStatic e, 0, 8:
      let re = cr*w[2*e] - ci*w[2*e+1]
      w[2*e+1] = cr*w[2*e+1] + ci*w[2*e]
      w[2*e] = re
    forStatic e, 0, 17: u[o + e*V] = w[e]
    usdev(d[nl + i], w)
  gpuFor(t, 2*np):
    let h = t div np
    let c = t - h*np
    var m = 0.0
    for j in c*nm ..< min(nl, (c+1)*nm): m = max(m, d[h*nl + j])
    dp[t] = m
  var hp = newSeq[float](2*np)
  gpuMemCpyToCpu(addr hp[0], dp, 2*np*sizeof(float))
  var s = gpuSum(i, nl, 2, [d[i], d[nl + i]])
  getDefaultComm().allReduce(addr s[0], 2)
  let c = 20.0
  let v = 4.0*float(g.lo.physVol)
  for h in 0..1:
    var m = 0.0
    for t in 0..<np: m = max(m, hp[h*np + t])
    rankMax(m)
    result[h] = (avg: sqrt(s[h]/(c*v)), max: sqrt(m/c))
  gpuFree(d)
  gpuFree(dp)
  g.fresh = false
