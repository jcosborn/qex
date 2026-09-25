## Staggered even-odd CG with the links, vectors and halo exchanges on the GPU.
##
## Device fields keep the qex SIMD layout of the gauge field, V lanes per
## outer site and even outer sites first: real c of site k at
## (k div V)*6V + c*V + k mod V, so neighboring threads read neighboring
## reals and the even part is the first 6*ne reals.  Links are blocked the
## same way with nl reals per link: 18; rows 0 and 1 only, 12, for links
## s W with W in SU(3) and s = +-1, rebuilding row 2 = s conj(row 0 x row 1);
## or rows 0 and 1 and d = det U, 14, for unitary links U, as HYP smeared
## ones, rebuilding row 2 = d conj(row 0 x row 1).
## Neighbor index j < n is a local site, possibly in another lane; j >= n is
## remote site j-n in the receive buffer of the exchange for that parity.
##
##   A x_e = 4 m^2 x_e - D_eo D_oe x_e,  D x = sum_mu U_mu(x) x(x+mu) - U_mu(x-mu)^+ x(x-mu)
##
## as stagD2ee.  The kernel producing a field also stores its boundary sites
## in the halo buffers of the neighbors, so a hop is: wait for the neighbors,
## then one kernel over the sites of the output parity.
##
## CG uses the Chronopoulos-Gear recurrence, one global sum per iteration:
##   p = r + beta p,  s = w + beta s,  x += alpha p,  r -= alpha s,  w = A r
##   gamma' = r.r,  delta = w.r,  beta = gamma'/gamma,
##   alpha = gamma'/(delta - beta gamma'/alpha)
## The second hop of A r sums both dot products and sends the boundary of w.
## Each rank then updates its halo copies of s and r with the same formulas,
## so the first hop needs no exchange, and the w exchange runs during the
## global sum.  With split, the second hop first runs on the even sites
## with local neighbors only, while the host waits for the exchange of t.
import std/[algorithm, macros]
import qex
import physics/qcdTypes
import solvers/solverBase
import backend/accel
import comms/[halo, halogpu]
import gauge/gaugeGpu
import base/metaUtils
import times

const nRed = 512  # slots per dot product of the second stage of the CG sums
const nBatch* {.intdefine.} = 3  # systems per hop in solveM of several systems; 3 beat 4, 6 and 8 on PVC
var hopSplit* = -1  ## CG second hop split around the exchange: 1 always, 0 never, -1 with off-node neighbors

type
  StagGpu*[V: static int; T] = object
    lo*: Layout[V]
    n*, ne*: int  # local sites, even sites
    nl*: int  # reals per link: 18, or 12 or 14 with row 2 rebuilt in the kernels
    lf*, lb*: ptr UncheckedArray[T]  # [4][n div V][nl][V]: U_mu(x), U_mu(x-mu)^+ or nil
    lh*: array[2, ptr UncheckedArray[T]]  # without lb, [nl][s.ex[1-q].nrecv]: U_mu(x-mu)^+ at the receive position of a remote x-mu, x of parity q
    nbr*: ptr UncheckedArray[int32]  # [8][n]: x+mu, then x-mu, with the link signs for nl = 12
    ex*: array[2, GpuHaloEx[T]]  # exchange of the sites of each parity
    red*: ptr UncheckedArray[float]  # [2][2*ne] dot products of the sites, two buffers, then [2*nRed] their slot sums
    hred*: ptr UncheckedArray[float]  # [2*nRed] pinned host copy of a buffer
    vec*: array[7, ptr UncheckedArray[T]]  # work vectors, 6*n reals each
    rh*, sh*: ptr UncheckedArray[T]  # CG halo copies of r and s, as s.ex[0].rbuf
    ord*: ptr UncheckedArray[int32]  # even sites, the ones with remote neighbors first, then the nin with local ones only
    nin*: int
    split*: bool  # second hop of the CG in two kernels around the exchange of t
    exB*: array[2, GpuHaloEx[T]]  # with batch, the exchanges of nBatch systems, system j in components 6j..6j+5
    redB*: ptr UncheckedArray[float]  # [2][2*nBatch][ne] dot products of the systems and sites, two buffers, then [2*nBatch*nRed] their slot sums
    hredB*: ptr UncheckedArray[float]  # [2*nBatch*nRed] pinned host copy of a buffer
    vecB*: array[7, ptr UncheckedArray[T]]  # work vectors of nBatch systems, 6*nBatch*n reals each
    rhB*, shB*: ptr UncheckedArray[T]  # CG halo copies of r and s of nBatch systems, as s.exB[0].rbuf

proc innerGeom*(lat, rg: seq[int]; v: int): seq[int] =
  ## Lanes in the dimensions not split across ranks first, then the largest
  ## local extent, at most 2 per dimension as qex lays out; the extent per
  ## lane stays even, so lanes share parity, and at least 4 in a dimension
  ## split across ranks, as qex halos go wrong with 2.
  result = newSeq[int](lat.len)
  for d in 0..<lat.len: result[d] = 1
  var k = v
  while k > 1:
    var best = -1
    for d in 0..<lat.len:
      let e = lat[d] div rg[d] div result[d]
      if result[d] == 1 and e mod (if rg[d] > 1: 8 else: 4) == 0 and (best < 0 or (rg[d] == 1, e) > (rg[best] == 1, lat[best] div rg[best] div result[best])):
        best = d
    if best < 0: qexError("no inner geometry for ", v, " lanes")
    result[best] *= 2
    k = k div 2

template vo(V, k: untyped): untyped =
  ## real 0 of site k in a vector field, real c is c*V further
  (k div V)*(6*V) + k mod V

template uo(V, k, nl: untyped): untyped =
  ## real 0 of site k in a link field, real e is e*V further
  (k div V)*(nl*V) + k mod V

template mdet(dr, di, m: untyped) =
  ## dr + i di = det m = sum_b m_2b (row 0 x row 1)_b, m as for load
  var dr, di = 0.0
  forStatic b, 0, 2:
    const b1 = (b+1) mod 3
    const b2 = (b+2) mod 3
    let xr = m[2*b1]*m[6+2*b2] - m[2*b1+1]*m[7+2*b2] - m[2*b2]*m[6+2*b1] + m[2*b2+1]*m[7+2*b1]
    let xi = m[2*b2]*m[7+2*b1] + m[2*b2+1]*m[6+2*b1] - m[2*b1]*m[7+2*b2] - m[2*b1+1]*m[6+2*b2]
    dr += m[12+2*b]*xr + m[13+2*b]*xi
    di += m[13+2*b]*xr - m[12+2*b]*xi

proc newStagGpu*[V: static int; E](g: openArray[Field[V,E]]; T: typedesc; reals = 12; fwd = -1;
                                   batch = false): StagGpu[V,T] =
  ## g are the phased links, as for newStag.  With reals = 12, if every link
  ## is s W with W in SU(3) and s = +-1, the device keeps rows 0 and 1 and
  ## the kernels rebuild row 2 = s conj(row 0 x row 1), with s in the sign of
  ## the neighbor index: j for s = 1, -1-j for s = -1.  With reals = 14 (or
  ## 12 when the links are not so), if every link is unitary, the device
  ## keeps rows 0 and 1 and d = det U, and the kernels rebuild row 2 =
  ## d conj(row 0 x row 1).  Otherwise all 18 reals.  With fwd = 1 the
  ## device keeps U_mu(x) only and a hop reads U_mu(x-mu) at x-mu, half the
  ## memory; fwd = -1 does so when both would take over 128 MB, as they
  ## would no longer stay in the L2 cache of a PVC tile between the hops.
  ## With batch, also the buffers of solveM of several systems.
  tic("newStagGpu")
  let lo = g[0].l
  let nd = lo.nDim
  let no = lo.nSitesOuter
  let n = V*no
  result.lo = lo
  result.n = n
  result.ne = V*lo.nEvenOuter
  let c = getDefaultComm()
  var w = newSeq[int](nd)
  for d in 0..<nd: w[d] = 1
  let hl = lo.makeHaloLayout(w, w)
  toc("layout")

  var off = newSeq[int32](nd)
  var offs = newSeq[seq[int32]](0)
  for mu in 0..<nd:
    off[mu] = 1
    offs.add off
    off[mu] = -1
    offs.add off
    off[mu] = 0
  var nb = newSeq[int32](2*nd*n)
  for p in 0..1:  # parity of the sites sent
    let hm = hl.makeHaloMap(c, offs, p)
    result.ex[p] = newGpuHaloEx[T](hm.gather, 6, n, V, c)
    if batch: result.exB[p] = newGpuHaloEx[T](hm.gather, 6*nBatch, n, V, c)
    let src = hl.haloSource(hm.gather)
    let q = 1 - p
    let o0 = if q == 0: 0 else: lo.nEvenOuter
    let o1 = if q == 0: lo.nEvenOuter else: no
    for o in o0..<o1:
      for mu in 0..<nd:
        for fb in 0..1:
          let e = int(if fb == 0: hl.neighborFwd[mu][o] else: hl.neighborBck[mu][o])
          for l in 0..<V:
            nb[(fb*nd+mu)*n + V*o + l] = if e < no: int32(V*e + l) else: src[V*(e-no) + l]
  var ord, bnd: seq[int32]
  for k in 0..<result.ne:
    var rem = false
    for d in 0..<2*nd: rem = rem or nb[d*n + k] >= n
    if rem: bnd.add int32(k) else: ord.add int32(k)
  result.nin = ord.len
  result.ord = (bnd & ord).toDevice  # remote neighbors first: the remote stores of a hop drain while the interior computes
  result.split = hopSplit > 0 or hopSplit < 0 and result.nin < result.ne and
                 false in result.ex[1].rpeer
  toc("neighbors")

  # lk[18*((fb*nd + mu)*n + k) + 6a+2b] (re), +1 (im): U_mu(x) for fb = 0,
  # U_mu(x-mu)^+ for fb = 1, their signs in sg and determinants in dt
  var lk = newSeq[float](2*nd*18*n)
  var sg = newSeq[int8](2*nd*n)
  var dt = newSeq[float](2*2*nd*n)
  var bad = [0.0, 0.0]  # links not s W, not unitary
  type U = eval(index(E, type(asSimd(0))))
  var u {.noInit.}: U
  template put(fb, mu, k: int) =
    let i = (fb*nd + mu)*n + k
    var m {.noInit.}: array[18, float]
    for a in 0..2:
      for b in 0..2:
        if fb == 0:
          m[6*a + 2*b] = float(u[a,b].re)
          m[6*a + 2*b + 1] = float(u[a,b].im)
        else:
          m[6*a + 2*b] = float(u[b,a].re)
          m[6*a + 2*b + 1] = -float(u[b,a].im)
    mdet(dr, di, m)
    var dp, dm, du, nn = 0.0
    for b in 0..2:  # xr + i xi = conj(row 0 x row 1)_b
      let b1 = (b+1) mod 3
      let b2 = (b+2) mod 3
      let xr = m[2*b1]*m[6+2*b2] - m[2*b1+1]*m[7+2*b2] - m[2*b2]*m[6+2*b1] + m[2*b2+1]*m[7+2*b1]
      let xi = m[2*b2]*m[7+2*b1] + m[2*b2+1]*m[6+2*b1] - m[2*b1]*m[7+2*b2] - m[2*b1+1]*m[6+2*b2]
      let pr = m[12+2*b] - xr
      let pi = m[13+2*b] - xi
      let mr = m[12+2*b] + xr
      let mi = m[13+2*b] + xi
      let ur = m[12+2*b] - (dr*xr - di*xi)
      let ui = m[13+2*b] - (dr*xi + di*xr)
      dp += pr*pr + pi*pi
      dm += mr*mr + mi*mi
      du += ur*ur + ui*ui
      nn += m[12+2*b]*m[12+2*b] + m[13+2*b]*m[13+2*b]
    for e in 0..17: lk[18*i + e] = m[e]
    dt[2*i] = dr
    dt[2*i+1] = di
    sg[i] = if dp <= dm: 1 else: -1
    if min(dp, dm) > 1e-24*nn: bad[0] += 1  # row 2 off by more than 1e-12
    if du > 1e-24*nn: bad[1] += 1
  for mu in 0..<nd:
    off[mu] = -1
    let hm = hl.makeHaloMap(c, @[off])
    off[mu] = 0
    let h = makeHalo(hl, g[mu])
    h.update(hm, c)
    for o in 0..<no:
      let e = hl.neighborBck[mu][o]
      for l in 0..<V:
        u := g[mu]{V*o+l}
        put(0, mu, V*o+l)
        u := h[e][asSimd(l)]
        put(1, mu, V*o+l)
  c.allReduce(addr bad[0], 2)
  let nl = if reals == 12 and bad[0] == 0: 12 elif reals <= 14 and bad[1] == 0: 14 else: 18
  result.nl = nl
  template lv(i, e: int): float =
    ## real e of link i with nl reals
    if nl == 14 and e >= 12: dt[2*i + e - 12] else: lk[18*i + e]
  let fw = fwd > 0 or fwd < 0 and 2*nd*nl*n*sizeof(T) > 128 shl 20
  var lf = newSeq[T](nd*nl*n)
  var lb = newSeq[T](if fw: 0 else: nd*nl*n)
  for mu in 0..<nd:
    for k in 0..<n:
      let o = nl*mu*n + uo(V, k, nl)
      for e in 0..<nl:
        lf[o + e*V] = T(lv(mu*n + k, e))
        if not fw: lb[o + e*V] = T(lv((nd + mu)*n + k, e))
  result.lf = lf.toDevice
  result.lb = lb.toDevice
  if fw:
    for q in 0..1:
      let nr = result.ex[1-q].nrecv
      var lh = newSeq[T](nl*nr)
      for k in (if q == 0: 0 else: result.ne) ..< (if q == 0: result.ne else: n):
        for mu in 0..<nd:
          let j = int nb[(nd+mu)*n + k]
          if j >= n:
            for e in 0..<nl: lh[e*nr + j-n] = T(lv((nd + mu)*n + k, e))
      result.lh[q] = lh.toDevice
  if nl == 12:
    for i in 0..<nb.len:
      if sg[i] < 0: nb[i] = -1 - nb[i]
  result.nbr = nb.toDevice
  toc("links")

  result.red = newSeq[float](4*result.ne + 2*nRed).toDevice
  result.hred = cast[ptr UncheckedArray[float]](gpuMallocHost(2*nRed*sizeof(float)))
  for i in 0..<result.vec.len:
    result.vec[i] = cast[ptr UncheckedArray[T]](gpuMalloc(6*n*sizeof(T)))
  let nh = max(1, 6*result.ex[0].nrecv)
  result.rh = cast[ptr UncheckedArray[T]](gpuMalloc(nh*sizeof(T)))
  result.sh = cast[ptr UncheckedArray[T]](gpuMalloc(nh*sizeof(T)))
  if batch:
    result.redB = newSeq[float](4*nBatch*result.ne + 2*nBatch*nRed).toDevice
    result.hredB = cast[ptr UncheckedArray[float]](gpuMallocHost(2*nBatch*nRed*sizeof(float)))
    for i in 0..<result.vecB.len:
      result.vecB[i] = cast[ptr UncheckedArray[T]](gpuMalloc(6*nBatch*n*sizeof(T)))
    result.rhB = cast[ptr UncheckedArray[T]](gpuMalloc(nBatch*nh*sizeof(T)))
    result.shB = cast[ptr UncheckedArray[T]](gpuMalloc(nBatch*nh*sizeof(T)))
  toc("vectors")

proc free*[V: static int; T](s: var StagGpu[V,T]) =
  for p in [s.lf, s.lb, s.lh[0], s.lh[1], s.rh, s.sh]:
    if p != nil: gpuFree(p)
  gpuFree(s.nbr)
  gpuFree(s.ord)
  gpuFree(s.red)
  gpuFreeHost(s.hred)
  for p in s.vec: gpuFree(p)
  for p in 0..1: s.ex[p].free
  if s.redB != nil:
    for p in [s.redB, s.rhB, s.shB]: gpuFree(p)
    gpuFreeHost(s.hredB)
    for p in s.vecB: gpuFree(p)
    for p in 0..1: s.exB[p].free

template load(m, u, o, st, e: untyped; nl: static int) =
  ## m[6a+2b] (re), m[6a+2b+1] (im) = U_ab from u[o + (6a+2b)*st], ...; for
  ## nl = 12 rows 0 and 1 only, row 2 = e conj(row 0 x row 1); for nl = 14
  ## rows 0 and 1 and d = det U at reals 12 and 13, row 2 = d conj(row 0 x row 1)
  when nl == 18:
    forStatic i, 0, 17: m[i] = u[o + i*st]
  else:
    forStatic i, 0, 11: m[i] = u[o + i*st]
    when nl == 14:
      let dr = u[o + 12*st]
      let di = u[o + 13*st]
    forStatic b, 0, 2:
      const b1 = (b+1) mod 3
      const b2 = (b+2) mod 3
      let xr = m[2*b1]*m[6+2*b2] - m[2*b1+1]*m[7+2*b2] - m[2*b2]*m[6+2*b1] + m[2*b2+1]*m[7+2*b1]
      let xi = m[2*b2]*m[7+2*b1] + m[2*b2+1]*m[6+2*b1] - m[2*b1]*m[7+2*b2] - m[2*b1+1]*m[6+2*b2]
      when nl == 12:
        m[12+2*b] = e*xr
        m[13+2*b] = e*xi
      else:
        m[12+2*b] = dr*xr - di*xi
        m[13+2*b] = dr*xi + di*xr

template hop(acc, u, o, v, st: untyped; sgn: static int; e: untyped; nl: static int) =
  ## acc += sgn U v, U as for load
  var m {.noInit.}: array[18, type(acc[0])]
  load(m, u, o, st, e, nl)
  forStatic r, 0, 2:
    var wr = m[6*r]*v[0] - m[6*r+1]*v[1]
    var wi = m[6*r]*v[1] + m[6*r+1]*v[0]
    forStatic b, 1, 2:
      wr += m[6*r+2*b]*v[2*b] - m[6*r+2*b+1]*v[2*b+1]
      wi += m[6*r+2*b]*v[2*b+1] + m[6*r+2*b+1]*v[2*b]
    when sgn > 0:
      acc[2*r] += wr
      acc[2*r+1] += wi
    else:
      acc[2*r] -= wr
      acc[2*r+1] -= wi

template hopA(acc, u, o, v, st: untyped; e: untyped; nl: static int) =
  ## acc -= U^+ v, U as for load
  var m {.noInit.}: array[18, type(acc[0])]
  load(m, u, o, st, e, nl)
  forStatic a, 0, 2:
    var wr = m[2*a]*v[0] + m[2*a+1]*v[1]
    var wi = m[2*a]*v[1] - m[2*a+1]*v[0]
    forStatic b, 1, 2:
      wr += m[6*b+2*a]*v[2*b] + m[6*b+2*a+1]*v[2*b+1]
      wi += m[6*b+2*a]*v[2*b+1] - m[6*b+2*a+1]*v[2*b]
    acc[2*a] -= wr
    acc[2*a+1] -= wi

template dots(i, d, y, yo, st, rs, ne: untyped) =
  ## rs[i] = y.y and rs[ne + i] = d.y of even site i, y at y[yo + c*st]
  var yy = 0.0
  var dy = 0.0
  forStatic c, 0, 5:
    let yc = float(y[yo + c*st])
    yy += yc*yc
    dy += float(d[c])*yc
  rs[i] = yy
  rs[ne + i] = dy

proc dslash[V: static int; T](s: StagGpu[V,T]; q: int; d, x, y, rb: ptr UncheckedArray[T];
                              a, b: T; rs, rz: ptr UncheckedArray[float]; dot, send: static bool;
                              j0 = 0; m = -1; nowait: static bool = false) =
  ## d = a y + b D x on the sites of parity q, the remote sites of x from rb,
  ## laid out as the receive buffer of s.ex[1-q]; for q = 0 on the sites
  ## s.ord[j0 ..< j0+m], all by default.  With send, stores d in the send
  ## slots of s.ex[q]; with dot, stores y.y and d.y of each site in rs.
  ## With nowait, returns before the kernel completes; gpuWaitAsync waits
  ## for it.
  let n = s.n
  let ne = s.ne
  let i0 = if q == 0: 0 else: s.ne
  let nk = if m >= 0: m elif q == 0: s.ne else: n - s.ne
  let od = s.ord
  let lf = s.lf
  let lb = s.lb
  let lh = s.lh[q]
  let nr = s.ex[1-q].nrecv
  let nb = s.nbr
  let ro = s.ex[1-q].rofs
  let rst = s.ex[1-q].rstr
  let sl = s.ex[q].sslot
  let sd = s.ex[q].sdst
  let st = s.ex[q].sstr
  let nsl = s.ex[q].nslot
  template fetch(v, j: untyped) =
    if j < n:
      let o = vo(V, j)
      forStatic c, 0, 5: v[c] = x[o + c*V]
    else:
      recvSite(ro, rst, rb, j-n, v)
  template body(i: untyped; nl: static int; fw: static bool) =
    let k = if q == 0: int od[j0 + i] else: i0 + i
    let ko = uo(V, k, nl)
    var acc {.noInit.}: array[6,T]
    forStatic c, 0, 5: acc[c] = T(0)
    forStatic mu, 0, 3:
      forStatic fb, 0, 1:
        let jj = int nb[(4*fb+mu)*n + k]
        let j = if nl == 12 and jj < 0: -1-jj else: jj
        let e = if jj < 0: T(-1) else: T(1)
        var v {.noInit.}: array[6,T]
        fetch(v, j)
        when fb == 0: hop(acc, lf, nl*mu*n + ko, v, V, 1, e, nl)
        elif fw:
          if j < n: hopA(acc, lf, nl*mu*n + uo(V, j, nl), v, V, e, nl)
          else: hop(acc, lh, j-n, v, nr, -1, e, nl)
        else: hop(acc, lb, nl*mu*n + ko, v, V, -1, e, nl)
    let yo = vo(V, k)
    if a == T(0):
      forStatic c, 0, 5: acc[c] = b*acc[c]
    else:
      forStatic c, 0, 5: acc[c] = a*y[yo + c*V] + b*acc[c]
    forStatic c, 0, 5: d[yo + c*V] = acc[c]
    when send:
      sendSite(sl, sd, st, nsl, n, k, 0, acc)
    when dot:
      dots(j0 + i, acc, y, yo, V, rs, ne)
  template kern(nl: static int; fw: static bool) =
    when nowait:
      gpuForAsync(i, nk): body(i, nl, fw)
    else:
      gpuFor(i, nk): body(i, nl, fw)
  if s.nl == 12:
    if lb == nil: kern(12, true) else: kern(12, false)
  elif s.nl == 14:
    if lb == nil: kern(14, true) else: kern(14, false)
  else:
    if lb == nil: kern(18, true) else: kern(18, false)

proc applyD2ee*[V: static int; T](s: StagGpu[V,T]; r, x, t: ptr UncheckedArray[T]; m2: float) =
  ## r_e = 4 m2 x_e - D_eo D_oe x_e, using t_o
  s.ex[0].pack(x)
  s.ex[0].start
  s.ex[0].wait
  s.dslash(1, t, x, x, s.ex[0].rbuf, T(0), T(1), nil, nil, dot = false, send = true)
  s.ex[1].start
  s.ex[1].wait
  s.dslash(0, r, t, x, s.ex[1].rbuf, T(4*m2), T(-1), nil, nil, dot = false, send = false)

proc applyD2eeCG[V: static int; T](s: StagGpu[V,T]; w, r, t: ptr UncheckedArray[T]; m2: float;
                                   rs, rz: ptr UncheckedArray[float]): array[2,float] =
  ## w_e = 4 m2 r_e - D_eo D_oe r_e with the remote sites of r from s.rh,
  ## returning the global r_e.r_e and w_e.r_e of the sites in rs, summed in
  ## a fixed order: slot k of nRed sums sites k, k+nRed, ..., the host the
  ## slots, so a solve repeats bit for bit.  Starts the
  ## exchange of the boundary of w in s.ex[0]; the caller waits for it.
  ## The first hop runs after the kernels already submitted, and the start
  ## of the exchange of t waits for them.
  tic("A r")
  s.dslash(1, t, r, r, s.rh, T(0), T(1), nil, nil, dot = false, send = true, nowait = true)
  toc("dslash oe")
  s.ex[1].start
  toc("start oe")
  if s.split:
    s.dslash(0, w, t, r, nil, T(4*m2), T(-1), rs, rz, dot = true, send = true, s.ne - s.nin, s.nin, nowait = true)
    toc("dslash eo local")
    s.ex[1].wait
    toc("wait oe")
    s.dslash(0, w, t, r, s.ex[1].rbuf, T(4*m2), T(-1), rs, rz, dot = true, send = true, 0, s.ne - s.nin)
    gpuWaitAsync()
  else:
    s.ex[1].wait(sync = false)
    toc("wait oe")
    s.dslash(0, w, t, r, s.ex[1].rbuf, T(4*m2), T(-1), rs, rz, dot = true, send = true, nowait = true)
  toc("dslash eo")
  s.ex[0].start
  toc("start eo")
  gpuWaitAsync()
  let ne = s.ne
  let sb = cast[ptr UncheckedArray[float]](addr s.red[4*ne])
  gpuFor(k, 2*nRed):
    let c = k div nRed
    var a = 0.0
    var j = k - c*nRed
    while j < ne:
      a += rs[c*ne + j]
      j += nRed
    sb[k] = a
  let h = s.hred
  gpuMemCpyToCpu(h, sb, 2*nRed*sizeof(float))
  for k in 0..<nRed:
    result[0] += h[k]
    result[1] += h[nRed+k]
  toc("dots")
  getDefaultComm().allReduce(addr result[0], 2)
  toc("sum")

proc redot[V: static int; T](s: StagGpu[V,T]; x, y: ptr UncheckedArray[T]): float =
  ## local sum of x_e.y_e
  gpuSum(i, 6*s.ne, 1, [float(x[i])*float(y[i])])[0]

proc update[V: static int; T](s: StagGpu[V,T]; x, r, p, sv, w: ptr UncheckedArray[T]; a, b: T) =
  ## p = r + b p, s = w + b s, x += a p, r -= a s on the even sites, and the
  ## same for the halo copies s.rh, s.sh from the w boundary received in
  ## s.ex[0]; a thread updates one real.  Returns before the kernel
  ## completes, the first hop of applyD2eeCG runs after it.
  let n6 = 6*s.ne
  let wh = s.ex[0].rbuf
  let rh = s.rh
  let sh = s.sh
  gpuForAsync(t, n6 + 6*s.ex[0].nrecv):
    if t < n6:
      let pk = r[t] + b*p[t]
      let sk = w[t] + b*sv[t]
      p[t] = pk
      sv[t] = sk
      x[t] += a*pk
      r[t] = r[t] - a*sk
    else:
      let k = t - n6
      let sk = wh[k] + b*sh[k]
      sh[k] = sk
      rh[k] = rh[k] - a*sk

proc upload*[V: static int; T](s: StagGpu[V,T]; d: ptr UncheckedArray[T]; f: Field) =
  ## even sites of f, on the layout of the links, to d
  when f.V != V or numberType(f[0]) isnot T:
    {.error: "upload: the field and StagGpu differ in lanes or precision".}
  gpuMemCpyToGpu(d, addr f[0], 6*s.ne*sizeof(T))

proc download*[V: static int; T](s: StagGpu[V,T]; f: Field; d: ptr UncheckedArray[T]) =
  ## d to the even sites of f, on the layout of the links
  when f.V != V or numberType(f[0]) isnot T:
    {.error: "download: the field and StagGpu differ in lanes or precision".}
  gpuMemCpyToCpu(addr f[0], d, 6*s.ne*sizeof(T))

proc convert[V: static int; T,U](s: StagGpu[V,T]; d: ptr UncheckedArray[T]; x: ptr UncheckedArray[U]; ns = 1) =
  ## d_e = x_e in the precision of d, for ns systems one after the other
  gpuFor(i, 6*ns*s.ne): d[i] = T(x[i])

proc addTo[V: static int; T,U](s: StagGpu[V,T]; d: ptr UncheckedArray[T]; x: ptr UncheckedArray[U]) =
  ## d_e += x_e
  gpuFor(i, 6*s.ne): d[i] += T(x[i])

proc resid[V: static int; T](s: StagGpu[V,T]; r, b, a: ptr UncheckedArray[T]) =
  ## r_e = b_e - a_e
  gpuFor(i, 6*s.ne): r[i] = b[i] - a[i]

proc cg[V: static int; T](s: StagGpu[V,T]; x, b: ptr UncheckedArray[T]; m, r2stop: float;
                          maxits, verb: int): tuple[its: int, r2: float] =
  ## Solves A x_e = b_e from x_e = 0 until the global |res|^2 <= r2stop,
  ## with s.vec[0..4] as work vectors.
  let c = getDefaultComm()
  let rr = s.vec[0]
  let p = s.vec[1]
  let sv = s.vec[2]
  let w = s.vec[3]
  let t = s.vec[4]
  gpuFor(i, 6*s.ne):
    x[i] = T(0)
    rr[i] = b[i]
    p[i] = T(0)
    sv[i] = T(0)
  var rs = s.red
  var rz = cast[ptr UncheckedArray[float]](addr s.red[2*s.ne])
  s.ex[0].pack(rr)
  s.ex[0].start
  s.ex[0].wait
  let hb = s.ex[0].rbuf
  let rh = s.rh
  let sh = s.sh
  gpuFor(k, 6*s.ex[0].nrecv):
    rh[k] = hb[k]
    sh[k] = T(0)
  var gd = s.applyD2eeCG(w, rr, t, m*m, rs, rz)
  swap(rs, rz)
  var r2 = gd[0]
  var alpha = gd[0]/gd[1]
  var beta = 0.0
  var itn = 0
  while itn < maxits and r2 > r2stop:
    tic("cg loop")
    s.ex[0].wait(sync = false)
    toc("wait w")
    s.update(x, rr, p, sv, w, T(alpha), T(beta))
    toc("update")
    gd = s.applyD2eeCG(w, rr, t, m*m, rs, rz)
    swap(rs, rz)
    toc("Ar")
    beta = gd[0]/r2
    alpha = gd[0]/(gd[1] - beta*gd[0]/alpha)
    r2 = gd[0]
    inc itn
    if verb > 1:
      echo "GPU CG iteration: ", itn, "  r2: ", r2
    if verb > 2:
      var rr2 = s.redot(rr, rr)
      var wr = s.redot(w, rr)
      c.allReduce(rr2)
      c.allReduce(wr)
      echo "  r.r: ", gd[0], " ", rr2, "  w.r: ", gd[1], " ", wr
  s.ex[0].wait
  (itn, r2)

proc solveEE*[V: static int; T](s: StagGpu[V,T]; r, x: Field; m: float; sp: var SolverParams) =
  ## Solves A r_e = x_e by CG from r_e = 0, as solveEE, stopping at
  ## |res|^2 <= sp.r2req |x_e|^2.
  tic("solveEE")
  let c = getDefaultComm()
  let b = s.vec[5]
  let xx = s.vec[6]
  s.upload(b, x)
  var b2 = s.redot(b, b)
  c.allReduce(b2)
  toc("setup")
  let t0 = epochTime()
  let (itn, r2) = s.cg(xx, b, m, sp.r2req*b2, sp.maxits, sp.verbosity)
  let secs = epochTime() - t0
  toc("cg")
  s.download(r, xx)
  sp.iterations = itn
  sp.seconds = secs
  sp.flops = float((2*8*72 + 60)*s.ne*itn)  # per even site: 2 hops of 8 links, 72 flops each, 60 for the vectors
  if sp.verbosity > 0:
    let gf = 1e-9*sp.flops*float(c.size)/secs
    echo "GPU CG iterations: ", itn, "  r2/b2: ", r2/b2, "  secs: ", secs, "  Gflops: ", gf
  toc("end")

proc cg[V: static int](s: StagGpu[V,float64]; ss: StagGpu[V,float32]; x, b: ptr UncheckedArray[float64];
                       m, r2stop, r2in: float; maxits, verb: int): tuple[its, nres: int, r2: float] =
  ## Solves A x_e = b_e from x_e = 0 until the global |res|^2 <= r2stop by
  ## single precision CGs in ss, each until |res|^2 drops by r2in, adding
  ## their solutions to x and restarting from the residual in double.
  let c = getDefaultComm()
  let rd = s.vec[0]
  let ad = s.vec[1]
  let t = s.vec[2]
  let rs = ss.vec[5]
  let e = ss.vec[6]
  gpuFor(i, 6*s.ne):
    x[i] = 0.0
    rd[i] = b[i]
  var r2 = s.redot(rd, rd)
  c.allReduce(r2)
  while result.its < maxits and r2 > r2stop:
    ss.convert(rs, rd)
    let (k, _) = ss.cg(e, rs, m, max(r2in*r2, 0.5*r2stop), maxits - result.its, verb)
    result.its += k
    s.addTo(x, e)
    s.applyD2ee(ad, x, t, m*m)
    s.resid(rd, b, ad)
    r2 = s.redot(rd, rd)
    c.allReduce(r2)
    inc result.nres
    if verb > 1:
      echo "GPU mixed CG restart: ", result.nres, "  iterations: ", result.its, "  r2: ", r2
  result.r2 = r2

proc solveEE*[V: static int](s: StagGpu[V,float64]; ss: StagGpu[V,float32]; r, x: Field; m: float;
                             sp: var SolverParams; r2in = 1e-6) =
  ## As solveEE, with the CG iterations in single precision: each restart
  ## solves A e = res in ss until |res|^2 drops by r2in, then r += e and
  ## res = x - A r in double.
  tic("solveEE mixed")
  let c = getDefaultComm()
  let b = s.vec[5]
  let xd = s.vec[6]
  s.upload(b, x)
  var b2 = s.redot(b, b)
  c.allReduce(b2)
  toc("setup")
  let t0 = epochTime()
  let (itn, nres, r2) = cg(s, ss, xd, b, m, sp.r2req*b2, r2in, sp.maxits, sp.verbosity)
  let secs = epochTime() - t0
  toc("cg")
  s.download(r, xd)
  sp.iterations = itn
  sp.seconds = secs
  sp.flops = float((2*8*72 + 60)*s.ne*(itn+nres))
  if sp.verbosity > 0:
    let gf = 1e-9*sp.flops*float(c.size)/secs
    echo "GPU mixed CG iterations: ", itn, "  restarts: ", nres, "  r2/b2: ", r2/b2, "  secs: ", secs, "  Gflops: ", gf
  toc("end")

# HMC with the links of a GpuGauge: s keeps 18 reals per link and both sets

proc stagSigns*[V: static int; E](g: openArray[Field[V,E]]): ptr UncheckedArray[float] =
  ## [8][n] on the device: the sign of U_mu(x), then of U_mu(x-mu), for g a
  ## unit gauge field after setBC and stagPhase
  let lo = g[0].l
  let nd = lo.nDim
  let no = lo.nSitesOuter
  let n = V*no
  let c = getDefaultComm()
  var w = newSeq[int](nd)
  for d in 0..<nd: w[d] = 1
  let hl = lo.makeHaloLayout(w, w)
  var off = newSeq[int32](nd)
  var sg = newSeq[float](2*nd*n)
  type U = eval(index(E, type(asSimd(0))))
  var u {.noInit.}: U
  for mu in 0..<nd:
    off[mu] = -1
    let hm = hl.makeHaloMap(c, @[off])
    off[mu] = 0
    let h = makeHalo(hl, g[mu])
    h.update(hm, c)
    for o in 0..<no:
      let e = hl.neighborBck[mu][o]
      for l in 0..<V:
        u := g[mu]{V*o+l}
        sg[mu*n + V*o+l] = float(u[0,0].re)
        u := h[e][asSimd(l)]
        sg[(nd+mu)*n + V*o+l] = float(u[0,0].re)
  sg.toDevice

proc setLinks*[V: static int; T](s: StagGpu[V,T]; g: var GpuGauge[V]; sg: ptr UncheckedArray[float]) =
  ## The links of s from those of g times the signs sg of stagSigns.  For
  ## links of 12 reals the signs in s.nbr are those of a unit gauge field
  ## with the phases, as newStagGpu gives, so the links of g must be SU(3).
  let n = g.n
  let u = g.u
  let nb = g.nb
  let ro = g.ex[0].rofs
  let rs = g.ex[0].rstr
  let h0 = g.ex[0].rbuf
  let h1 = g.ex[1].rbuf
  let h2 = g.ex[2].rbuf
  let h3 = g.ex[3].rbuf
  let nl = s.nl
  let ne = s.ne
  let fn = s.nbr
  let lf = s.lf
  let lb = s.lb
  let lh0 = s.lh[0]
  let lh1 = s.lh[1]
  let nr0 = s.ex[1].nrecv
  let nr1 = s.ex[0].nrecv
  template put(d, o, st, m, sc: untyped) =
    ## sc m to d[o + e*st], e < nl, for sc = +-1
    if nl == 14:
      forStatic e, 0, 11: d[o + e*st] = T(sc*m[e])
      mdet(dr, di, m)
      d[o + 12*st] = T(sc*dr)
      d[o + 13*st] = T(sc*di)
    else:
      for e in 0..<nl: d[o + e*st] = T(sc*m[e])
  forLinks(g, mu, k):
    let ol = nl*mu*n + uo(V, k, nl)
    var x {.noInit.}, y {.noInit.}: array[18, float]
    mload(x, u, 18*mu*n + lo18(V, k), V)
    let sf = sg[mu*n + k]
    put(lf, ol, V, x, sf)
    case mu
    of 0: link(x, 0, int nb[nbB(0)*n + k])
    of 1: link(x, 1, int nb[nbB(1)*n + k])
    of 2: link(x, 2, int nb[nbB(2)*n + k])
    else: link(x, 3, int nb[nbB(3)*n + k])
    let sb = sg[(4+mu)*n + k]
    forStatic a, 0, 2:  # y = sb U_mu(x-mu)^+
      forStatic b, 0, 2:
        y[6*a+2*b] = sb*x[6*b+2*a]
        y[6*a+2*b+1] = -sb*x[6*b+2*a+1]
    if lb != nil:
      put(lb, ol, V, y, 1.0)
    else:
      let jj = int fn[(4+mu)*n + k]
      let j = if jj < 0 and nl == 12: -1-jj else: jj
      if j >= n:
        if k < ne: put(lh0, j-n, nr0, y, 1.0)
        else: put(lh1, j-n, nr1, y, 1.0)

proc norm2*[V: static int; T](s: StagGpu[V,T]; x: ptr UncheckedArray[T]): float =
  ## global |x|^2 over all sites
  result = gpuSum(i, 6*s.n, 1, [float(x[i])*float(x[i])])[0]
  getDefaultComm().allReduce(result)

proc applyM*[V: static int; T](s: StagGpu[V,T]; d, x: ptr UncheckedArray[T]; m: float) =
  ## d_e = m x_e + D_eo x_o/2, the even sites of stag.D(d, x, m)
  getDefaultComm().barrier
  s.ex[1].pack(x)
  s.ex[1].start
  s.ex[1].wait
  s.dslash(0, d, x, x, s.ex[1].rbuf, T(m), T(0.5), nil, nil, dot = false, send = false)

proc applyMfull*[V: static int; T](s: StagGpu[V,T]; d, x: ptr UncheckedArray[T]; m: float) =
  ## d = (m + D/2) x on all sites, as stag.D
  let c = getDefaultComm()
  c.barrier
  s.ex[1].pack(x)
  s.ex[1].start
  s.ex[1].wait
  s.dslash(0, d, x, x, s.ex[1].rbuf, T(m), T(0.5), nil, nil, dot = false, send = false)
  c.barrier
  s.ex[0].pack(x)
  s.ex[0].start
  s.ex[0].wait
  s.dslash(1, d, x, x, s.ex[0].rbuf, T(m), T(0.5), nil, nil, dot = false, send = false)

proc addSolve(sp: var SolverParams; its: int; r2: float) =
  ## the statistics of one solve ending at |res|^2/|b|^2 = r2, as the CPU solvers keep them
  inc sp.calls
  sp.iterations += its
  sp.iterationsMax = max(sp.iterationsMax, its)
  sp.r2.push r2

proc solveM*[V: static int](s: StagGpu[V,float]; x, b: ptr UncheckedArray[float]; m: float;
                            sp: var SolverParams; ss: ptr StagGpu[V,float32] = nil; r2in = 1e-6;
                            full = false) =
  ## Solves (m + D/2) x = b for b_o = 0, as stag.solve, until the residual
  ## of the even system drops by sp.r2req, in mixed precision with ss:
  ##   x_e = A^-1 4m b_e,  x_o = -D_oe x_e/(2m)
  ## With full, b_o may be nonzero, as solveReconL of the CPU solve:
  ##   A x_e = 4m b_e - 2 D_eo b_o,  x_o = b_o/m - D_oe x_e/(2m)
  let c = getDefaultComm()
  let b4 = s.vec[5]
  if full:
    c.barrier
    s.ex[1].pack(b)
    s.ex[1].start
    s.ex[1].wait
    s.dslash(0, b4, b, b, s.ex[1].rbuf, 4.0*m, -2.0, nil, nil, dot = false, send = false)
  else:
    gpuFor(i, 6*s.ne): b4[i] = 4.0*m*b[i]
  var b2 = s.redot(b4, b4)
  c.allReduce(b2)
  if ss == nil:
    let (itn, r2) = s.cg(x, b4, m, sp.r2req*b2, sp.maxits, sp.verbosity)
    sp.addSolve(itn, r2/b2)
  else:
    let (itn, _, r2) = cg(s, ss[], x, b4, m, sp.r2req*b2, r2in, sp.maxits, sp.verbosity)
    sp.addSolve(itn, r2/b2)
  c.barrier
  s.ex[0].pack(x)
  s.ex[0].start
  s.ex[0].wait
  if full: s.dslash(1, x, x, b, s.ex[0].rbuf, 1.0/m, -0.5/m, nil, nil, dot = false, send = false)
  else: s.dslash(1, x, x, x, s.ex[0].rbuf, 0.0, -0.5/m, nil, nil, dot = false, send = false)

proc outerM*[V: static int](s: StagGpu[V,float]; f, x: ptr UncheckedArray[float]; t: float) =
  ## f_mu(y) += t x(y) x(y+mu)^+ over all sites, f like GpuGauge.u: the one
  ## link force before the phases and the smearing pullback, as fforce of
  ## staghmc_sh accumulates it over the mass terms
  getDefaultComm().barrier
  for q in 0..1:
    s.ex[q].pack(x)
    s.ex[q].start
  for q in 0..1: s.ex[q].wait
  let n = s.n
  let ne = s.ne
  let nl = s.nl
  let nb = s.nbr
  let ro0 = s.ex[0].rofs
  let rs0 = s.ex[0].rstr
  let rb0 = s.ex[0].rbuf
  let ro1 = s.ex[1].rofs
  let rs1 = s.ex[1].rstr
  let rb1 = s.ex[1].rbuf
  gpuFor(i, 4*n):
    let mu = i div n
    let k = i - mu*n
    let jj = int nb[mu*n + k]
    let j = if nl == 12 and jj < 0: -1-jj else: jj
    var v {.noInit.}: array[6, float]
    if j < n:
      let o = vo(V, j)
      forStatic c, 0, 5: v[c] = x[o + c*V]
    elif k < ne: recvSite(ro1, rs1, rb1, j-n, v)
    else: recvSite(ro0, rs0, rb0, j-n, v)
    let yo = vo(V, k)
    let o = 18*mu*n + lo18(V, k)
    forStatic a, 0, 2:
      forStatic b, 0, 2:
        f[o + (6*a+2*b)*V] += t*(x[yo + 2*a*V]*v[2*b] + x[yo + (2*a+1)*V]*v[2*b+1])
        f[o + (6*a+2*b+1)*V] += t*(x[yo + (2*a+1)*V]*v[2*b] - x[yo + 2*a*V]*v[2*b+1])

proc forceM*[V: static int](s: StagGpu[V,float]; g: GpuGauge[V]; sg: ptr UncheckedArray[float];
                            p, x: ptr UncheckedArray[float]; t: float) =
  ## p_mu(y) -= t e(y) TAH(x(y) (U_mu(y) x(y+mu))^+) with the links of g
  ## times the signs sg of stagSigns, e = 1 on even and -1 on odd sites,
  ## as oneLinkForce of staghmc
  getDefaultComm().barrier
  for q in 0..1:
    s.ex[q].pack(x)
    s.ex[q].start
  for q in 0..1: s.ex[q].wait
  let n = s.n
  let ne = s.ne
  let u = g.u
  let nl = s.nl
  let nb = s.nbr
  let ro0 = s.ex[0].rofs
  let rs0 = s.ex[0].rstr
  let rb0 = s.ex[0].rbuf
  let ro1 = s.ex[1].rofs
  let rs1 = s.ex[1].rstr
  let rb1 = s.ex[1].rbuf
  gpuFor(i, 4*n):
    let mu = i div n
    let k = i - mu*n
    let jj = int nb[mu*n + k]
    let j = if nl == 12 and jj < 0: -1-jj else: jj
    var v {.noInit.}: array[6, float]
    if j < n:
      let o = vo(V, j)
      forStatic c, 0, 5: v[c] = x[o + c*V]
    elif k < ne: recvSite(ro1, rs1, rb1, j-n, v)
    else: recvSite(ro0, rs0, rb0, j-n, v)
    let o = 18*mu*n + lo18(V, k)
    var m {.noInit.}: array[18, float]
    mload(m, u, o, V)
    let sf = sg[mu*n + k]
    var w {.noInit.}: array[6, float]  # U_mu(y) x(y+mu)
    forStatic a, 0, 2:
      var wr = 0.0
      var wi = 0.0
      forStatic b, 0, 2:
        wr += m[6*a+2*b]*v[2*b] - m[6*a+2*b+1]*v[2*b+1]
        wi += m[6*a+2*b]*v[2*b+1] + m[6*a+2*b+1]*v[2*b]
      w[2*a] = sf*wr
      w[2*a+1] = sf*wi
    let yo = vo(V, k)
    var y {.noInit.}: array[6, float]
    forStatic c, 0, 5: y[c] = x[yo + c*V]
    var f {.noInit.}, g {.noInit.}: array[18, float]  # f_ab = y_a conj(w_b)
    forStatic a, 0, 2:
      forStatic b, 0, 2:
        f[6*a+2*b] = y[2*a]*w[2*b] + y[2*a+1]*w[2*b+1]
        f[6*a+2*b+1] = y[2*a+1]*w[2*b] - y[2*a]*w[2*b+1]
    mtah(g, f)
    let sc = if k < ne: t else: -t
    forStatic e, 0, 17: p[o + e*V] -= sc*g[e]

# Several systems at once, nBatch per hop.  The vectors of the right-hand
# sides stay in their own arrays, 6n reals each as above, but one hop
# kernel loads each link once for the systems still iterating (bits of a
# mask) and one exchange of 6 nBatch reals per site carries them, system j
# in components 6j..6j+5.  With n_1 >= n_2 >= ... the iterations of the
# systems sorted by mass, groups of nBatch consecutive systems load the
# links n_1 + n_(nBatch+1) + ... times.  Needs newStagGpu(batch = true).

proc sysB[T](a: ptr UncheckedArray[T]; n6: int): array[nBatch, ptr UncheckedArray[T]] =
  ## nBatch vectors of n6 reals in one allocation
  for j in 0..<nBatch: result[j] = cast[ptr UncheckedArray[T]](addr a[j*n6])

# The kernels capture an array of the systems as nBatch variables a0, a1,
# ..., which both backends pass by value, the OpenMP one only as scalars.

macro unpackB(a: typed): untyped =
  ## let a0 = a[0], a1 = a[1], ...
  result = newStmtList()
  for j in 0..<nBatch:
    result.add newLetStmt(ident($a & $j), newTree(nnkBracketExpr, a, newLit(j)))

macro pickB(a, j: untyped): untyped =
  ## a_j of unpackB(a)
  result = nnkCaseStmt.newTree(j)
  for k in 0..<nBatch-1:
    result.add nnkOfBranch.newTree(newLit(k), ident($a & $k))
  result.add nnkElse.newTree(ident($a & $(nBatch-1)))

template hopm(acc, m, v: untyped; j, sgn: static int) =
  ## acc_j += sgn m v, m of 18 reals, v the 6 reals of system j
  forStatic r, 0, 2:
    var wr = m[6*r]*v[0] - m[6*r+1]*v[1]
    var wi = m[6*r]*v[1] + m[6*r+1]*v[0]
    forStatic b, 1, 2:
      wr += m[6*r+2*b]*v[2*b] - m[6*r+2*b+1]*v[2*b+1]
      wi += m[6*r+2*b]*v[2*b+1] + m[6*r+2*b+1]*v[2*b]
    when sgn > 0:
      acc[6*j+2*r] += wr
      acc[6*j+2*r+1] += wi
    else:
      acc[6*j+2*r] -= wr
      acc[6*j+2*r+1] -= wi

template hopAm(acc, m, v: untyped; j: static int) =
  ## acc_j -= m^+ v
  forStatic a, 0, 2:
    var wr = m[2*a]*v[0] + m[2*a+1]*v[1]
    var wi = m[2*a]*v[1] - m[2*a+1]*v[0]
    forStatic b, 1, 2:
      wr += m[6*b+2*a]*v[2*b] + m[6*b+2*a+1]*v[2*b+1]
      wi += m[6*b+2*a]*v[2*b+1] - m[6*b+2*a+1]*v[2*b]
    acc[6*j+2*a] -= wr
    acc[6*j+2*a+1] -= wi

proc dslashB[V: static int; T](s: StagGpu[V,T]; q: int; d, x, y: array[nBatch, ptr UncheckedArray[T]]; rb: ptr UncheckedArray[T];
                               a, b: array[nBatch, T]; act: int; rs, rz: ptr UncheckedArray[float]; dot, send: static bool;
                               j0 = 0; m = -1; nowait: static bool = false) =
  ## dslash for the systems j with bit j of act set: d_j = a_j y_j + b_j D x_j,
  ## the remote sites of x_j from components 6j.. of rb, laid out as the
  ## receive buffer of s.exB[1-q]; send and dot as in dslash with s.exB[q]
  ## and y_j.y_j, d_j.y_j of even site p at rs[2j ne + p], rs[(2j+1) ne + p]
  let n = s.n
  let ne = s.ne
  let i0 = if q == 0: 0 else: s.ne
  let nk = if m >= 0: m elif q == 0: s.ne else: n - s.ne
  let od = s.ord
  let lf = s.lf
  let lb = s.lb
  let lh = s.lh[q]
  let nr = s.ex[1-q].nrecv
  let nb = s.nbr
  let ro = s.exB[1-q].rofs
  let rst = s.exB[1-q].rstr
  let sl = s.exB[q].sslot
  let sd = s.exB[q].sdst
  let st = s.exB[q].sstr
  let nsl = s.exB[q].nslot
  unpackB(x)
  unpackB(y)
  unpackB(d)
  unpackB(a)
  unpackB(b)
  template fetch(v, j: untyped; xj: untyped; c0: static int) =
    if j < n:
      let o = vo(V, j)
      forStatic c, 0, 5: v[c] = xj[o + c*V]
    else:
      let o = int ro[j-n]
      let sk = int rst[j-n]
      forStatic c, 0, 5: v[c] = rb[o + (c0+c)*sk]
  template body(i: untyped; nl: static int; fw: static bool) =
    let k = if q == 0: int od[j0 + i] else: i0 + i
    let ko = uo(V, k, nl)
    var acc {.noInit.}: array[6*nBatch,T]
    forStatic c, 0, 6*nBatch-1: acc[c] = T(0)
    forStatic mu, 0, 3:
      forStatic fb, 0, 1:
        let jj = int nb[(4*fb+mu)*n + k]
        let j = if nl == 12 and jj < 0: -1-jj else: jj
        let e = if jj < 0: T(-1) else: T(1)
        var mm {.noInit.}: array[18,T]
        var v {.noInit.}: array[6,T]
        when fb == 0: load(mm, lf, nl*mu*n + ko, V, e, nl)
        elif fw:
          if j < n: load(mm, lf, nl*mu*n + uo(V, j, nl), V, e, nl)
          else: load(mm, lh, j-n, nr, e, nl)
        else: load(mm, lb, nl*mu*n + ko, V, e, nl)
        forStatic js, 0, nBatch-1:
          if (act and (1 shl js)) != 0:
            fetch(v, j, pickB(x, js), 6*js)
            when fb == 0: hopm(acc, mm, v, js, 1)
            elif fw:
              if j < n: hopAm(acc, mm, v, js)
              else: hopm(acc, mm, v, js, -1)
            else: hopm(acc, mm, v, js, -1)
    let yo = vo(V, k)
    forStatic js, 0, nBatch-1:
      if (act and (1 shl js)) != 0:
        let aj = pickB(a, js)
        let bj = pickB(b, js)
        let yj = pickB(y, js)
        let dj = pickB(d, js)
        if aj == T(0):
          forStatic c, 0, 5: acc[6*js+c] = bj*acc[6*js+c]
        else:
          forStatic c, 0, 5: acc[6*js+c] = aj*yj[yo + c*V] + bj*acc[6*js+c]
        forStatic c, 0, 5: dj[yo + c*V] = acc[6*js+c]
        when send:
          var w {.noInit.}: array[6,T]
          forStatic c, 0, 5: w[c] = acc[6*js+c]
          sendSite(sl, sd, st, nsl, n, k, 6*js, w)
        when dot:
          var yy = 0.0
          var dy = 0.0
          forStatic c, 0, 5:
            let yc = float(yj[yo + c*V])
            yy += yc*yc
            dy += float(acc[6*js+c])*yc
          let p = j0 + i
          rs[2*js*ne + p] = yy
          rs[(2*js+1)*ne + p] = dy
  template kern(nl: static int; fw: static bool) =
    when nowait:
      gpuForAsync(i, nk, 16): body(i, nl, fw)
    else:
      gpuFor(i, nk, 16): body(i, nl, fw)
  if s.nl == 12:
    if lb == nil: kern(12, true) else: kern(12, false)
  elif s.nl == 14:
    if lb == nil: kern(14, true) else: kern(14, false)
  else:
    if lb == nil: kern(18, true) else: kern(18, false)

proc packB[T](ex: GpuHaloEx[T]; f: array[nBatch, ptr UncheckedArray[T]]; act: int) =
  ## pack for nBatch vectors of 6 reals per site into components 6j.. of an
  ## exchange of 6 nBatch reals per site, the systems j in act
  let v = ex.v
  let ns = ex.nsend
  let si = ex.sidx
  let sd = ex.sdst
  let st = ex.sstr
  unpackB(f)
  gpuForAsync(t, 6*nBatch*ns):
    let c = t div ns
    let k = t - c*ns
    let js = c div 6
    if (act and (1 shl js)) != 0:
      let j = int si[k]
      sd[k][c*int st[k]] = pickB(f, js)[((j div v)*6 + c - 6*js)*v + j mod v]

proc applyD2eeB[V: static int; T](s: StagGpu[V,T]; r, x, t: array[nBatch, ptr UncheckedArray[T]]; m2: array[nBatch, float]; act: int) =
  ## r_j = 4 m2_j x_j - D_eo D_oe x_j on the even sites, using t_o, for the systems in act
  var a, z, o, mo: array[nBatch, T]
  for j in 0..<nBatch:
    a[j] = T(4*m2[j])
    o[j] = T(1)
    mo[j] = T(-1)
  s.exB[0].packB(x, act)
  s.exB[0].start
  s.exB[0].wait
  s.dslashB(1, t, x, x, s.exB[0].rbuf, z, o, act, nil, nil, dot = false, send = true)
  s.exB[1].start
  s.exB[1].wait
  s.dslashB(0, r, t, x, s.exB[1].rbuf, a, mo, act, nil, nil, dot = false, send = false)

proc applyD2eeCGB[V: static int; T](s: StagGpu[V,T]; w, r, t: array[nBatch, ptr UncheckedArray[T]]; m2: array[nBatch, float]; act: int;
                                    rs, rz: ptr UncheckedArray[float]): array[2*nBatch, float] =
  ## applyD2eeCG for the systems in act, returning r_j.r_j and w_j.r_j at 2j and 2j+1
  var a, z, o, mo: array[nBatch, T]
  for j in 0..<nBatch:
    a[j] = T(4*m2[j])
    o[j] = T(1)
    mo[j] = T(-1)
  s.dslashB(1, t, r, r, s.rhB, z, o, act, nil, nil, dot = false, send = true, nowait = true)
  s.exB[1].start
  if s.split:
    s.dslashB(0, w, t, r, nil, a, mo, act, rs, rz, dot = true, send = true, s.ne - s.nin, s.nin, nowait = true)
    s.exB[1].wait
    s.dslashB(0, w, t, r, s.exB[1].rbuf, a, mo, act, rs, rz, dot = true, send = true, 0, s.ne - s.nin)
    gpuWaitAsync()
  else:
    s.exB[1].wait(sync = false)
    s.dslashB(0, w, t, r, s.exB[1].rbuf, a, mo, act, rs, rz, dot = true, send = true, nowait = true)
  s.exB[0].start
  gpuWaitAsync()
  let ne = s.ne  # the site sums in a fixed order, as applyD2eeCG
  let sb = cast[ptr UncheckedArray[float]](addr s.redB[4*nBatch*ne])
  gpuFor(k, 2*nBatch*nRed):
    let c = k div nRed
    var sm = 0.0
    var j = k - c*nRed
    while j < ne:
      sm += rs[c*ne + j]
      j += nRed
    sb[k] = sm
  let h = s.hredB
  gpuMemCpyToCpu(h, sb, 2*nBatch*nRed*sizeof(float))
  for c in 0..<2*nBatch:
    if (act and (1 shl (c div 2))) != 0:
      for k in 0..<nRed: result[c] += h[c*nRed + k]
  getDefaultComm().allReduce(addr result[0], 2*nBatch)

proc updateB[V: static int; T](s: StagGpu[V,T]; x: array[nBatch, ptr UncheckedArray[T]]; r, p, sv, w: ptr UncheckedArray[T]; a, b: array[nBatch, T]; act: int) =
  ## update for the systems j in act with a_j and b_j; r, p, sv, w hold the
  ## systems one after the other, 6 ne reals each
  let n6 = 6*s.ne
  let nr = s.exB[0].nrecv
  let wh = s.exB[0].rbuf
  let rh = s.rhB
  let sh = s.shB
  let ro = s.exB[0].rofs
  let rst = s.exB[0].rstr
  unpackB(x)
  unpackB(a)
  unpackB(b)
  gpuForAsync(t, nBatch*n6 + 6*nBatch*nr):
    if t < nBatch*n6:
      let j = t div n6
      if (act and (1 shl j)) != 0:
        let aj = pickB(a, j)
        let bj = pickB(b, j)
        let pk = r[t] + bj*p[t]
        let sk = w[t] + bj*sv[t]
        p[t] = pk
        sv[t] = sk
        pickB(x, j)[t - j*n6] += aj*pk
        r[t] = r[t] - aj*sk
    else:
      let u = t - nBatch*n6
      let c = u div nr
      let j = c div 6
      if (act and (1 shl j)) != 0:
        let k = int ro[u - c*nr] + c*int rst[u - c*nr]
        let sk = wh[k] + pickB(b, j)*sh[k]
        sh[k] = sk
        rh[k] = rh[k] - pickB(a, j)*sk

proc cgB[V: static int; T](s: StagGpu[V,T]; x, b: array[nBatch, ptr UncheckedArray[T]]; m, r2stop: array[nBatch, float]; act0: int;
                           maxits, verb: int): tuple[its: array[nBatch, int], r2: array[nBatch, float]] =
  ## cg for the systems in act0 at once, each until its |res|^2 <= r2stop_j;
  ## a converged system leaves the kernels.  s.vecB[0..3] hold the systems
  ## one after the other, 6 ne reals each, s.vecB[4] with 6 n.
  let n6 = 6*s.ne
  let ra = s.vecB[0]
  let pa = s.vecB[1]
  let sa = s.vecB[2]
  let wa = s.vecB[3]
  let rr = sysB(ra, n6)
  let w = sysB(wa, n6)
  let t = sysB(s.vecB[4], 6*s.n)
  unpackB(x)
  unpackB(b)
  gpuFor(i, nBatch*n6):  # x = 0, r = b, p = s = 0
    let j = i div n6
    if (act0 and (1 shl j)) != 0:
      let o = i - j*n6
      pickB(x, j)[o] = T(0)
      ra[i] = pickB(b, j)[o]
      pa[i] = T(0)
      sa[i] = T(0)
  var rs = s.redB
  var rz = cast[ptr UncheckedArray[float]](addr s.redB[2*nBatch*s.ne])
  var act = act0
  s.exB[0].packB(rr, act)
  s.exB[0].start
  s.exB[0].wait
  let hb = s.exB[0].rbuf
  let rh = s.rhB
  let sh = s.shB
  gpuFor(k, 6*nBatch*s.exB[0].nrecv):
    rh[k] = hb[k]
    sh[k] = T(0)
  var m2: array[nBatch, float]
  for j in 0..<nBatch: m2[j] = m[j]*m[j]
  var gd = s.applyD2eeCGB(w, rr, t, m2, act, rs, rz)
  swap(rs, rz)
  var alpha, beta, r2: array[nBatch, float]
  for j in 0..<nBatch:
    if (act and (1 shl j)) != 0:
      r2[j] = gd[2*j]
      alpha[j] = gd[2*j]/gd[2*j+1]
      if r2[j] <= r2stop[j]: act = act and not (1 shl j)
  var itn = 0
  while itn < maxits and act != 0:
    s.exB[0].wait(sync = false)
    var a, bb: array[nBatch, T]
    for j in 0..<nBatch:
      a[j] = T(alpha[j])
      bb[j] = T(beta[j])
    s.updateB(x, ra, pa, sa, wa, a, bb, act)
    gd = s.applyD2eeCGB(w, rr, t, m2, act, rs, rz)
    swap(rs, rz)
    inc itn
    for j in 0..<nBatch:
      if (act and (1 shl j)) != 0:
        beta[j] = gd[2*j]/r2[j]
        alpha[j] = gd[2*j]/(gd[2*j+1] - beta[j]*gd[2*j]/alpha[j])
        r2[j] = gd[2*j]
        result.its[j] = itn
        if r2[j] <= r2stop[j]: act = act and not (1 shl j)
    if verb > 1:
      echo "GPU CGB iteration: ", itn, "  r2: ", r2
  s.exB[0].wait
  result.r2 = r2

proc cgB[V: static int](s: StagGpu[V,float64]; ss: StagGpu[V,float32]; x, b: array[nBatch, ptr UncheckedArray[float64]];
                        m, r2stop: array[nBatch, float]; act0: int; r2in: float; maxits, verb: int):
                        tuple[its: array[nBatch, int], r2: array[nBatch, float]] =
  ## the mixed precision cg for the systems in act0: single precision cgB in
  ## ss of the systems not converged, each until its |res|^2 drops by r2in,
  ## restarted in double
  let c = getDefaultComm()
  let n6 = 6*s.ne
  let rd = sysB(s.vecB[0], n6)
  let ad = sysB(s.vecB[1], n6)
  let t = sysB(s.vecB[2], 6*s.n)
  let rs = sysB(ss.vecB[5], n6)
  let e = sysB(ss.vecB[6], n6)
  var r2: array[nBatch, float]
  var act = 0
  for j in 0..<nBatch:
    if (act0 and (1 shl j)) != 0:
      let xj = x[j]
      let rj = rd[j]
      let bj = b[j]
      gpuFor(i, n6):
        xj[i] = 0.0
        rj[i] = bj[i]
      r2[j] = s.redot(rd[j], rd[j])
  c.allReduce(addr r2[0], nBatch)
  for j in 0..<nBatch:
    if (act0 and (1 shl j)) != 0 and r2[j] > r2stop[j]: act = act or (1 shl j)
  var its = 0
  while its < maxits and act != 0:
    ss.convert(rs[0], rd[0], nBatch)
    var r2s: array[nBatch, float]
    for j in 0..<nBatch: r2s[j] = max(r2in*r2[j], 0.5*r2stop[j])
    let (k, _) = ss.cgB(e, rs, m, r2s, act, maxits - its, verb)
    its += max(k)
    for j in 0..<nBatch:
      if (act and (1 shl j)) != 0:
        result.its[j] += k[j]
        s.addTo(x[j], e[j])
    var m2: array[nBatch, float]
    for j in 0..<nBatch: m2[j] = m[j]*m[j]
    s.applyD2eeB(ad, x, t, m2, act)
    var d: array[nBatch, float]
    for j in 0..<nBatch:
      if (act and (1 shl j)) != 0:
        s.resid(rd[j], b[j], ad[j])
        d[j] = s.redot(rd[j], rd[j])
    c.allReduce(addr d[0], nBatch)
    for j in 0..<nBatch:
      if (act and (1 shl j)) != 0:
        r2[j] = d[j]
        if r2[j] <= r2stop[j]: act = act and not (1 shl j)
    if verb > 1:
      echo "GPU mixed CGB restart, iterations: ", result.its, "  r2: ", r2
  result.r2 = r2

proc solveB[V: static int](s: StagGpu[V,float]; x, b: array[nBatch, ptr UncheckedArray[float]]; m: array[nBatch, float];
                           sp: var array[nBatch, SolverParams]; act0: int; ss: ptr StagGpu[V,float32]; r2in: float;
                           full: bool) =
  ## solveM of the systems in act0 at once, each until its residual drops
  ## by sp_j.r2req; a hop loads each link once for the systems still
  ## iterating
  let c = getDefaultComm()
  let n6 = 6*s.ne
  let b4 = sysB(s.vecB[5], n6)
  var b2: array[nBatch, float]
  var r2stop: array[nBatch, float]
  var maxits = 0
  var a4, t2, am, bm: array[nBatch, float]
  for j in 0..<nBatch:
    a4[j] = 4.0*m[j]
    t2[j] = -2.0
    am[j] = if full: 1.0/m[j] else: 0.0
    bm[j] = -0.5/m[j]
  if full:
    c.barrier
    s.exB[1].packB(b, act0)
    s.exB[1].start
    s.exB[1].wait
    s.dslashB(0, b4, b, b, s.exB[1].rbuf, a4, t2, act0, nil, nil, dot = false, send = false)
  for j in 0..<nBatch:
    if (act0 and (1 shl j)) != 0:
      let bj = b[j]
      let b4j = b4[j]
      let mj = a4[j]
      if not full: gpuFor(i, n6): b4j[i] = mj*bj[i]
      b2[j] = s.redot(b4j, b4j)
      maxits = max(maxits, sp[j].maxits)
  c.allReduce(addr b2[0], nBatch)
  for j in 0..<nBatch: r2stop[j] = sp[j].r2req*b2[j]
  let (its, r2) = if ss == nil: s.cgB(x, b4, m, r2stop, act0, maxits, sp[0].verbosity)
                  else: cgB(s, ss[], x, b4, m, r2stop, act0, r2in, maxits, sp[0].verbosity)
  for j in 0..<nBatch:
    if (act0 and (1 shl j)) != 0: sp[j].addSolve(its[j], r2[j]/b2[j])
  c.barrier
  s.exB[0].packB(x, act0)
  s.exB[0].start
  s.exB[0].wait
  s.dslashB(1, x, x, if full: b else: x, s.exB[0].rbuf, am, bm, act0, nil, nil, dot = false, send = false)

proc solveM*[V: static int](s: StagGpu[V,float]; x, b: openArray[ptr UncheckedArray[float]]; m: openArray[float];
                            sp: var openArray[SolverParams]; ss: ptr StagGpu[V,float32] = nil; r2in = 1e-6;
                            full = false) =
  ## solveM of the systems j, x_j = M(m_j)^-1 b_j, b_j,o = 0 unless full,
  ## lightest first in groups of nBatch: a hop loads each link once for the
  ## systems of a group still iterating, and each system stops at its own
  ## sp_j.r2req.  A group of one takes the solveM of one system.
  var o = newSeq[(float, int)](x.len)  # (|m_j|, j), lightest first
  for j in 0..<x.len: o[j] = (abs(m[j]), j)
  o.sort
  var k = 0
  while k < o.len:
    let nc = min(nBatch, o.len - k)
    if nc == 1:
      let i = o[k][1]
      s.solveM(x[i], b[i], m[i], sp[i], ss, r2in, full)
    else:
      var xg, bg: array[nBatch, ptr UncheckedArray[float]]
      var mg: array[nBatch, float]
      var pg: array[nBatch, SolverParams]
      for j in 0..<nBatch:  # unused slots repeat the last system, outside the mask
        let i = o[k + min(j, nc-1)][1]
        xg[j] = x[i]
        bg[j] = b[i]
        mg[j] = m[i]
        pg[j] = sp[i]
      s.solveB(xg, bg, mg, pg, (1 shl nc) - 1, ss, r2in, full)
      for j in 0..<nc: sp[o[k+j][1]] = pg[j]
    k += nc
