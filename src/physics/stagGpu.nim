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

const nRed = sumHost  # atomic slots, pinned host reals per dot product of the CG sums
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
    red*: ptr UncheckedArray[float]  # [2][2*max(ne, nRed)] dot products of the sites or slots, two buffers, then the workspace of sumFixed
    hred*: ptr UncheckedArray[float]  # [2*nRed] pinned host copy of a buffer
    vec*: array[7, ptr UncheckedArray[T]]  # work vectors, 6*n reals each
    rh*, sh*: ptr UncheckedArray[T]  # CG halo copies of r and s, as s.ex[0].rbuf
    ord*: ptr UncheckedArray[int32]  # even sites, the ones with remote neighbors first, then the nin with local ones only
    nin*: int
    split*: bool  # second hop of the CG in two kernels around the exchange of t
    fixed*: bool  # the CG sums its dot products in a fixed order (sumFixed), repeating bit for bit, else atomically
    exB*: array[2, GpuHaloEx[T]]  # with batch, the exchanges of nBatch systems, system j in components 6j..6j+5
    redB*: ptr UncheckedArray[float]  # [2][2*nBatch][max(ne, nRed)] dot products of the systems and sites or slots, two buffers, then the workspace of sumFixed
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

  result.red = newSeq[float](4*max(result.ne, nRed) + 2*(result.ne div 15 + 16)).toDevice
  result.hred = cast[ptr UncheckedArray[float]](gpuMallocHost(2*nRed*sizeof(float)))
  for i in 0..<result.vec.len:
    result.vec[i] = cast[ptr UncheckedArray[T]](gpuMalloc(6*n*sizeof(T)))
  let nh = max(1, 6*result.ex[0].nrecv)
  result.rh = cast[ptr UncheckedArray[T]](gpuMalloc(nh*sizeof(T)))
  result.sh = cast[ptr UncheckedArray[T]](gpuMalloc(nh*sizeof(T)))
  if batch:
    result.redB = newSeq[float](4*nBatch*max(result.ne, nRed) + 2*nBatch*(result.ne div 15 + 16)).toDevice
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
    gpuFree(s.redB)
    for p in [s.rhB, s.shB]: gpuFree(p)
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

template dots(i, d, y, yo, st, rs, rz, ne, fx: untyped) =
  ## y.y and d.y of even site i, y at y[yo + c*st]: with fx at rs[i] and
  ## rs[ne + i], else added to the slots i mod nRed of rs, zeroing those of
  ## rz, the buffer of the next sum
  var yy = 0.0
  var dy = 0.0
  forStatic c, 0, 5:
    let yc = float(y[yo + c*st])
    yy += yc*yc
    dy += float(d[c])*yc
  if fx:
    rs[i] = yy
    rs[ne + i] = dy
  else:
    if i < nRed:
      rz[i] = 0.0
      rz[nRed + i] = 0.0
    gpuAtomicAdd(rs, i mod nRed, yy)
    gpuAtomicAdd(rs, nRed + i mod nRed, dy)

proc dslash[V: static int; T](s: StagGpu[V,T]; q: int; d, x, y, rb: ptr UncheckedArray[T];
                              a, b: T; rs, rz: ptr UncheckedArray[float]; dot, send: static bool;
                              j0 = 0; m = -1; nowait: static bool = false) =
  ## d = a y + b D x on the sites of parity q, the remote sites of x from rb,
  ## laid out as the receive buffer of s.ex[1-q]; for q = 0 on the sites
  ## s.ord[j0 ..< j0+m], all by default.  With send, stores d in the send
  ## slots of s.ex[q]; with dot, y.y and d.y of each site in rs as dots.
  ## With nowait, returns before the kernel completes; gpuWaitAsync waits
  ## for it.
  let n = s.n
  let ne = s.ne
  let fx = s.fixed
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
      dots(j0 + i, acc, y, yo, V, rs, rz, ne, fx)
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
  ## returning the global r_e.r_e and w_e.r_e from rs, with s.fixed summed in
  ## a fixed order by sumFixed, so a solve repeats bit for bit.  Starts the
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
  if s.fixed:
    sumFixed(result, rs, cast[ptr UncheckedArray[float]](addr s.red[4*max(s.ne, nRed)]), s.hred, 2, s.ne)
  else:
    let h = s.hred
    gpuMemCpyToCpu(h, rs, 2*nRed*sizeof(float))
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
  let bo = 2*max(s.ne, nRed)
  let red = s.red
  gpuFor(k, 2*nRed):  # the atomic slots of both buffers
    red[k] = 0.0
    red[bo + k] = 0.0
  var rs = s.red
  var rz = cast[ptr UncheckedArray[float]](addr s.red[bo])
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

# Several systems at once, up to nBatch per hop.  The vectors of the
# right-hand sides stay in their own arrays, 6n reals each as above, but
# one hop kernel loads each link once for the C systems still iterating
# and one exchange carries them: slot j of the nBatch slots in the
# components 6j..6j+5 of an exchange of 6 nBatch reals per site.  The
# kernels come in widths C = 1..nBatch, so a slot without a system costs
# nothing.  The slots take the systems sorted by mass, lightest first, a
# slot the next one when its system stops, so the hops load the links as
# often as the busiest slot iterates: max(n_1, n_3 + n_4) times for four
# systems of n_1 >= n_2 >= n_3 >= n_4 iterations in three slots.  Needs
# newStagGpu(batch = true).

proc sysB[T](a: ptr UncheckedArray[T]; n6: int): array[nBatch, ptr UncheckedArray[T]] =
  ## nBatch vectors of n6 reals in one allocation
  for j in 0..<nBatch: result[j] = cast[ptr UncheckedArray[T]](addr a[j*n6])

# The kernels capture an array of the systems as the variables a0, a1,
# ..., which both backends pass by value, the OpenMP one only as scalars.

macro unpackB(a: typed; n: static int): untyped =
  ## let a0 = a[0], a1 = a[1], ..., a(n-1) = a[n-1]
  result = newStmtList()
  for j in 0..<n:
    result.add newLetStmt(ident($a & $j), newTree(nnkBracketExpr, a, newLit(j)))

macro pickS(a: untyped; j: static int): untyped =
  ## aj of unpackB(a, n) for a static j
  ident($a & $j)

macro pickB(a, j: untyped; n: static int): untyped =
  ## aj of unpackB(a, n)
  result = nnkCaseStmt.newTree(j)
  for k in 0..<n-1:
    result.add nnkOfBranch.newTree(newLit(k), ident($a & $k))
  result.add nnkElse.newTree(ident($a & $(n-1)))

template forActive(act: int; body: untyped) =
  ## body for the systems in act with C, their number, and so: array[C, int]
  ## their slots, C static
  var cnt = 0
  for j in 0..<nBatch:
    if (act and (1 shl j)) != 0: inc cnt
  forStatic cw, 1, nBatch:
    if cnt == cw:
      const C {.inject.} = cw
      var so {.inject, noInit.}: array[C, int]
      var q = 0
      for j in 0..<nBatch:
        if (act and (1 shl j)) != 0:
          so[q] = j
          inc q
      body

template sel(C: static int; a, so: untyped): untyped =
  ## a[so[q]] for q in 0..<C
  var r {.noInit.}: array[C, typeof(a[0])]
  for q in 0..<C: r[q] = a[so[q]]
  r

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

proc dslashB[V: static int; C: static int; T](s: StagGpu[V,T]; q: int; d, x, y: array[C, ptr UncheckedArray[T]];
                                              rb: ptr UncheckedArray[T]; a, b: array[C, T]; so: array[C, int];
                                              rs, rz: ptr UncheckedArray[float]; dot, send: static bool;
                                              j0 = 0; m = -1; nowait: static bool = false) =
  ## dslash for C systems, d_j = a_j y_j + b_j D x_j: system j of slot so_j
  ## reads the remote sites of x_j from the components 6 so_j.. of rb, laid
  ## out as the receive buffer of s.exB[1-q], sends as dslash with s.exB[q]
  ## to its components there, and stores y_j.y_j, d_j.y_j of even site p at
  ## rs[2 so_j ne + p] and rs[(2 so_j + 1) ne + p] with s.fixed, else adds
  ## them to the slots 2 so_j nRed + p mod nRed and those after nRed
  let n = s.n
  let ne = s.ne
  let fx = s.fixed
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
  unpackB(x, C)
  unpackB(y, C)
  unpackB(d, C)
  unpackB(a, C)
  unpackB(b, C)
  unpackB(so, C)
  template fetch(v, j: untyped; xj: untyped; c0: untyped) =
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
    var acc {.noInit.}: array[6*C,T]
    forStatic c, 0, 6*C-1: acc[c] = T(0)
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
        forStatic js, 0, C-1:
          fetch(v, j, pickS(x, js), 6*pickS(so, js))
          when fb == 0: hopm(acc, mm, v, js, 1)
          elif fw:
            if j < n: hopAm(acc, mm, v, js)
            else: hopm(acc, mm, v, js, -1)
          else: hopm(acc, mm, v, js, -1)
    let yo = vo(V, k)
    forStatic js, 0, C-1:
      let aj = pickS(a, js)
      let bj = pickS(b, js)
      let yj = pickS(y, js)
      let dj = pickS(d, js)
      let sj = pickS(so, js)
      if aj == T(0):
        forStatic c, 0, 5: acc[6*js+c] = bj*acc[6*js+c]
      else:
        forStatic c, 0, 5: acc[6*js+c] = aj*yj[yo + c*V] + bj*acc[6*js+c]
      forStatic c, 0, 5: dj[yo + c*V] = acc[6*js+c]
      when send:
        var w {.noInit.}: array[6,T]
        forStatic c, 0, 5: w[c] = acc[6*js+c]
        sendSite(sl, sd, st, nsl, n, k, 6*sj, w)
      when dot:
        var yy = 0.0
        var dy = 0.0
        forStatic c, 0, 5:
          let yc = float(yj[yo + c*V])
          yy += yc*yc
          dy += float(acc[6*js+c])*yc
        let p = j0 + i
        if fx:
          rs[2*sj*ne + p] = yy
          rs[(2*sj+1)*ne + p] = dy
        else:
          gpuAtomicAdd(rs, 2*sj*nRed + p mod nRed, yy)
          gpuAtomicAdd(rs, (2*sj+1)*nRed + p mod nRed, dy)
    when dot:
      let p = j0 + i
      if not fx and p < nRed:  # buffer of the next sum
        forStatic c, 0, 2*nBatch-1: rz[c*nRed + p] = 0.0
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

proc packB[C: static int; T](ex: GpuHaloEx[T]; f: array[C, ptr UncheckedArray[T]]; so: array[C, int]) =
  ## pack for C vectors of 6 reals per site, f_j into the components
  ## 6 so_j.. of an exchange of 6 nBatch reals per site
  let v = ex.v
  let ns = ex.nsend
  let si = ex.sidx
  let sd = ex.sdst
  let st = ex.sstr
  unpackB(f, C)
  unpackB(so, C)
  gpuForAsync(t, 6*C*ns):
    let c = t div ns
    let k = t - c*ns
    let js = c div 6
    let cc = c - 6*js
    let j = int si[k]
    sd[k][(6*pickB(so, js, C) + cc)*int st[k]] = pickB(f, js, C)[((j div v)*6 + cc)*v + j mod v]

proc applyD2eeB[V: static int; C: static int; T](s: StagGpu[V,T]; r, x, t: array[C, ptr UncheckedArray[T]];
                                                 m2: array[C, float]; so: array[C, int]) =
  ## r_j = 4 m2_j x_j - D_eo D_oe x_j on the even sites, using t_o
  var a, z, o, mo: array[C, T]
  for j in 0..<C:
    a[j] = T(4*m2[j])
    o[j] = T(1)
    mo[j] = T(-1)
  s.exB[0].packB(x, so)
  s.exB[0].start
  s.exB[0].wait
  s.dslashB(1, t, x, x, s.exB[0].rbuf, z, o, so, nil, nil, dot = false, send = true)
  s.exB[1].start
  s.exB[1].wait
  s.dslashB(0, r, t, x, s.exB[1].rbuf, a, mo, so, nil, nil, dot = false, send = false)

proc applyD2eeCGB[V: static int; C: static int; T](s: StagGpu[V,T]; w, r, t: array[C, ptr UncheckedArray[T]];
                                                   m2: array[C, float]; so: array[C, int];
                                                   rs, rz: ptr UncheckedArray[float]): array[2*nBatch, float] =
  ## applyD2eeCG for C systems, returning r_j.r_j and w_j.r_j at 2 so_j and 2 so_j + 1
  var a, z, o, mo: array[C, T]
  for j in 0..<C:
    a[j] = T(4*m2[j])
    o[j] = T(1)
    mo[j] = T(-1)
  s.dslashB(1, t, r, r, s.rhB, z, o, so, nil, nil, dot = false, send = true, nowait = true)
  s.exB[1].start
  if s.split:
    s.dslashB(0, w, t, r, nil, a, mo, so, rs, rz, dot = true, send = true, s.ne - s.nin, s.nin, nowait = true)
    s.exB[1].wait
    s.dslashB(0, w, t, r, s.exB[1].rbuf, a, mo, so, rs, rz, dot = true, send = true, 0, s.ne - s.nin)
    gpuWaitAsync()
  else:
    s.exB[1].wait(sync = false)
    s.dslashB(0, w, t, r, s.exB[1].rbuf, a, mo, so, rs, rz, dot = true, send = true, nowait = true)
  s.exB[0].start
  gpuWaitAsync()
  if s.fixed:  # the site sums in a fixed order, as applyD2eeCG
    var r: array[2*nBatch, float]
    sumFixed(r, rs, cast[ptr UncheckedArray[float]](addr s.redB[4*nBatch*max(s.ne, nRed)]), s.hredB, 2*nBatch, s.ne)
    for j in so:
      for c in 2*j..2*j+1: result[c] = r[c]
  else:
    let h = s.hredB
    gpuMemCpyToCpu(h, rs, 2*nBatch*nRed*sizeof(float))
    for j in so:
      for c in 2*j..2*j+1:
        for k in 0..<nRed: result[c] += h[c*nRed + k]
  getDefaultComm().allReduce(addr result[0], 2*nBatch)

proc updateB[V: static int; C: static int; T](s: StagGpu[V,T]; x: array[C, ptr UncheckedArray[T]]; so: array[C, int];
                                              r, p, sv, w: ptr UncheckedArray[T]; a, b: array[C, T]) =
  ## update for C systems with a_j and b_j; r, p, sv, w hold the slots one
  ## after the other, 6 ne reals each, system j in slot so_j
  let n6 = 6*s.ne
  let nr = s.exB[0].nrecv
  let wh = s.exB[0].rbuf
  let rh = s.rhB
  let sh = s.shB
  let ro = s.exB[0].rofs
  let rst = s.exB[0].rstr
  unpackB(x, C)
  unpackB(so, C)
  unpackB(a, C)
  unpackB(b, C)
  gpuForAsync(t, C*n6 + 6*C*nr):
    if t < C*n6:
      let js = t div n6
      let o = t - js*n6
      let i = pickB(so, js, C)*n6 + o
      let aj = pickB(a, js, C)
      let bj = pickB(b, js, C)
      let pk = r[i] + bj*p[i]
      let sk = w[i] + bj*sv[i]
      p[i] = pk
      sv[i] = sk
      pickB(x, js, C)[o] += aj*pk
      r[i] = r[i] - aj*sk
    else:
      let u = t - C*n6
      let c = u div nr
      let js = c div 6
      let cc = 6*pickB(so, js, C) + c - 6*js
      let k = int ro[u - c*nr] + cc*int rst[u - c*nr]
      let sk = wh[k] + pickB(b, js, C)*sh[k]
      sh[k] = sk
      rh[k] = rh[k] - pickB(a, js, C)*sk

type CgB[T] = object
  ## the systems in the nBatch slots of the batched CG: slot j iterates on
  ## x_j from startB until |res_j|^2 <= stop_j or its_j = mits_j
  x: array[nBatch, ptr UncheckedArray[T]]
  m2, stop, alpha, beta, r2: array[nBatch, float]
  its, mits: array[nBatch, int]
  act, fresh: int  # the slots iterating, those of them started since the last hop
  rs, rz: ptr UncheckedArray[float]  # the dot products of the next hop and of the one after

proc initB[V: static int; T](s: StagGpu[V,T]): CgB[T] =
  ## no systems, the atomic slots of both dot product buffers zero
  let bo = 2*nBatch*max(s.ne, nRed)
  let red = s.redB
  gpuFor(k, 2*nBatch*nRed):
    red[k] = 0.0
    red[bo + k] = 0.0
  result.rs = s.redB
  result.rz = cast[ptr UncheckedArray[float]](addr s.redB[bo])

proc startB[V: static int; T](s: StagGpu[V,T]; c: var CgB[T]; j: int; x, b: ptr UncheckedArray[T];
                              m, stop: float; mits: int) =
  ## slot j solves A x_e = b_e from x_e = 0: r = b, p = s = 0, the halo of r
  ## exchanged by the next stepB
  let n6 = 6*s.ne
  let o = j*n6
  let ra = s.vecB[0]
  let pa = s.vecB[1]
  let sa = s.vecB[2]
  gpuFor(i, n6):
    x[i] = T(0)
    ra[o + i] = b[i]
    pa[o + i] = T(0)
    sa[o + i] = T(0)
  c.x[j] = x
  c.m2[j] = m*m
  c.stop[j] = stop
  c.its[j] = 0
  c.mits[j] = mits
  c.act = c.act or (1 shl j)
  c.fresh = c.fresh or (1 shl j)

proc stepB[V: static int; T](s: StagGpu[V,T]; c: var CgB[T]; verb: int): int =
  ## Iterates the slots in c.act, the kernels as wide as the slots
  ## iterating, until some stop, and returns those; the next update of the
  ## others is queued.  s.vecB[0..3] hold r, p, s, w of the slots one after
  ## the other, 6 ne reals each, s.vecB[4] t with 6 n.
  let n6 = 6*s.ne
  let ra = s.vecB[0]
  let pa = s.vecB[1]
  let sa = s.vecB[2]
  let wa = s.vecB[3]
  let rr = sysB(ra, n6)
  let w = sysB(wa, n6)
  let t = sysB(s.vecB[4], 6*s.n)
  if c.fresh != 0:  # the halo copies r = b and s = 0 of the slots started
    forActive(c.fresh):
      s.exB[0].packB(sel(C, rr, so), so)
    s.exB[0].start
    s.exB[0].wait
    let nr = s.exB[0].nrecv
    let ro = s.exB[0].rofs
    let rst = s.exB[0].rstr
    let hb = s.exB[0].rbuf
    let rh = s.rhB
    let sh = s.shB
    forActive(c.fresh):
      unpackB(so, C)
      gpuFor(u, 6*C*nr):
        let cc = u div nr
        let js = cc div 6
        let k = int ro[u - cc*nr] + (6*pickB(so, js, C) + cc - 6*js)*int rst[u - cc*nr]
        rh[k] = hb[k]
        sh[k] = T(0)
  while true:
    var gd: array[2*nBatch, float]
    forActive(c.act):
      gd = s.applyD2eeCGB(sel(C, w, so), sel(C, rr, so), sel(C, t, so), sel(C, c.m2, so), so, c.rs, c.rz)
    swap(c.rs, c.rz)
    var d = 0
    for j in 0..<nBatch:
      if (c.act and (1 shl j)) != 0:
        if (c.fresh and (1 shl j)) != 0:
          c.beta[j] = 0.0
          c.alpha[j] = gd[2*j]/gd[2*j+1]
        else:
          c.beta[j] = gd[2*j]/c.r2[j]
          c.alpha[j] = gd[2*j]/(gd[2*j+1] - c.beta[j]*gd[2*j]/c.alpha[j])
          inc c.its[j]
        c.r2[j] = gd[2*j]
        if c.r2[j] <= c.stop[j] or c.its[j] >= c.mits[j]: d = d or (1 shl j)
    c.fresh = 0
    c.act = c.act and not d
    if verb > 1:
      echo "GPU CGB iterations: ", c.its, "  r2: ", c.r2
    s.exB[0].wait(sync = c.act == 0)
    if c.act != 0:
      forActive(c.act):
        var a, bb: array[C, T]
        for q in 0..<C:
          a[q] = T(c.alpha[so[q]])
          bb[q] = T(c.beta[so[q]])
        s.updateB(sel(C, c.x, so), so, ra, pa, sa, wa, a, bb)
    if d != 0: return d

proc rhsB[V: static int](s: StagGpu[V,float]; b: array[nBatch, ptr UncheckedArray[float]]; m: array[nBatch, float];
                         act: int; full: bool): array[nBatch, float] =
  ## b4_j = 4 m_j b_j,e (- 2 D_eo b_j,o with full) in slot j of s.vecB[5]
  ## for the slots j in act, returning the global |b4_j|^2
  let n6 = 6*s.ne
  let b4 = sysB(s.vecB[5], n6)
  var a4, t2: array[nBatch, float]
  for j in 0..<nBatch:
    a4[j] = 4.0*m[j]
    t2[j] = -2.0
  if full:
    getDefaultComm().barrier
    forActive(act):
      s.exB[1].packB(sel(C, b, so), so)
      s.exB[1].start
      s.exB[1].wait
      s.dslashB(0, sel(C, b4, so), sel(C, b, so), sel(C, b, so), s.exB[1].rbuf, sel(C, a4, so), sel(C, t2, so), so,
                nil, nil, dot = false, send = false)
  for j in 0..<nBatch:
    if (act and (1 shl j)) != 0:
      let bj = b[j]
      let dj = b4[j]
      let mj = a4[j]
      if not full: gpuFor(i, n6): dj[i] = mj*bj[i]
      result[j] = s.redot(dj, dj)
  getDefaultComm().allReduce(addr result[0], nBatch)

proc solveM*[V: static int](s: StagGpu[V,float]; x, b: openArray[ptr UncheckedArray[float]]; m: openArray[float];
                            sp: var openArray[SolverParams]; ss: ptr StagGpu[V,float32] = nil; r2in = 1e-6;
                            full = false) =
  ## solveM of the systems j, x_j = M(m_j)^-1 b_j, b_j,o = 0 unless full,
  ## each stopping at its own sp_j.r2req.  The slots of the batched CG take
  ## the systems lightest first, a slot the next one when its system stops.
  ## With ss, a slot whose single precision CG stops restarts its system in
  ## double at once; up to nBatch systems restart together once all their
  ## single precision CGs stopped.  One system takes the solveM of one
  ## system.
  let ns = x.len
  if ns == 1:
    s.solveM(x[0], b[0], m[0], sp[0], ss, r2in, full)
    return
  let c = getDefaultComm()
  let n6 = 6*s.ne
  let verb = sp[0].verbosity
  let b4 = sysB(s.vecB[5], n6)
  var o = newSeq[(float, int)](ns)  # (|m_j|, j), lightest first
  for j in 0..<ns: o[j] = (abs(m[j]), j)
  o.sort
  var q = 0  # the systems taken
  var sl: array[nBatch, int]  # the system in each slot
  var b2: array[nBatch, float]  # its |b4|^2
  var fr = (1 shl nBatch) - 1  # the free slots
  template take: int =
    ## the free slots take the next systems, b4 in s.vecB[5]; the slots taken
    var a = 0
    var bs: array[nBatch, ptr UncheckedArray[float]]
    var ms: array[nBatch, float]
    for j in 0..<nBatch:
      if (fr and (1 shl j)) != 0 and q < ns:
        let i = o[q][1]
        inc q
        sl[j] = i
        bs[j] = b[i]
        ms[j] = m[i]
        a = a or (1 shl j)
    if a != 0:
      let r = s.rhsB(bs, ms, a, full)
      for j in 0..<nBatch:
        if (a and (1 shl j)) != 0: b2[j] = r[j]
    fr = fr and not a
    a
  if ss == nil:
    var st = s.initB
    while q < ns or st.act != 0:
      let a = take
      for j in 0..<nBatch:
        if (a and (1 shl j)) != 0:
          let i = sl[j]
          s.startB(st, j, x[i], b4[j], m[i], sp[i].r2req*b2[j], sp[i].maxits)
      let d = s.stepB(st, verb)
      for j in 0..<nBatch:
        if (d and (1 shl j)) != 0: sp[sl[j]].addSolve(st.its[j], st.r2[j]/b2[j])
      fr = fr or d
  else:
    let rd = sysB(s.vecB[0], n6)
    let ad = sysB(s.vecB[1], n6)
    let t = sysB(s.vecB[2], 6*s.n)
    let rs = sysB(ss.vecB[5], n6)
    let e = sysB(ss.vecB[6], n6)
    var st = ss[].initB
    var its: array[nBatch, int]
    var r2: array[nBatch, float]
    var d = 0  # the slots whose single precision CG stopped
    while true:
      var go = 0  # the slots with a residual rd and its r2
      if d != 0:  # x += e, rd = b4 - A x in double
        var xs: array[nBatch, ptr UncheckedArray[float]]
        var m2: array[nBatch, float]
        for j in 0..<nBatch:
          if (d and (1 shl j)) != 0:
            its[j] += st.its[j]
            xs[j] = x[sl[j]]
            m2[j] = m[sl[j]]*m[sl[j]]
            s.addTo(xs[j], e[j])
        forActive(d):
          s.applyD2eeB(sel(C, ad, so), sel(C, xs, so), sel(C, t, so), sel(C, m2, so), so)
        var dd: array[nBatch, float]
        for j in 0..<nBatch:
          if (d and (1 shl j)) != 0:
            s.resid(rd[j], b4[j], ad[j])
            dd[j] = s.redot(rd[j], rd[j])
        c.allReduce(addr dd[0], nBatch)
        for j in 0..<nBatch:
          if (d and (1 shl j)) != 0: r2[j] = dd[j]
        go = d
      while true:  # the systems that stop leave their slots to the next ones
        for j in 0..<nBatch:
          if (go and (1 shl j)) != 0:
            let i = sl[j]
            if r2[j] <= sp[i].r2req*b2[j] or its[j] >= sp[i].maxits:
              sp[i].addSolve(its[j], r2[j]/b2[j])
              go = go and not (1 shl j)
              fr = fr or (1 shl j)
        let a = take
        if a == 0: break
        for j in 0..<nBatch:
          if (a and (1 shl j)) != 0:  # x_e = 0, rd = b4
            let xi = x[sl[j]]
            let rj = rd[j]
            let bj = b4[j]
            gpuFor(k, n6):
              xi[k] = 0.0
              rj[k] = bj[k]
            its[j] = 0
            r2[j] = b2[j]
        go = go or a
      for j in 0..<nBatch:  # A e = rd in single precision, until its |res|^2 drops by r2in
        if (go and (1 shl j)) != 0:
          let i = sl[j]
          let r2s = sp[i].r2req*b2[j]
          ss[].convert(rs[j], rd[j])
          ss[].startB(st, j, e[j], rs[j], m[i], max(r2in*r2[j], 0.5*r2s), sp[i].maxits - its[j])
      if st.act == 0: break
      d = ss[].stepB(st, verb)
      while ns <= nBatch and st.act != 0:  # started together: restart together, the hops shared
        d = d or ss[].stepB(st, verb)
  for k in countup(0, ns-1, nBatch):  # x_o = b_o/m - D_oe x_e/(2m), without b_o unless full
    let nc = min(nBatch, ns - k)
    var xg, bg: array[nBatch, ptr UncheckedArray[float]]
    var am, bm: array[nBatch, float]
    for j in 0..<nc:
      xg[j] = x[k+j]
      bg[j] = b[k+j]
      am[j] = if full: 1.0/m[k+j] else: 0.0
      bm[j] = -0.5/m[k+j]
    c.barrier
    forActive((1 shl nc) - 1):
      let xs = sel(C, xg, so)
      s.exB[0].packB(xs, so)
      s.exB[0].start
      s.exB[0].wait
      s.dslashB(1, xs, xs, if full: sel(C, bg, so) else: xs, s.exB[0].rbuf, sel(C, am, so), sel(C, bm, so), so,
                nil, nil, dot = false, send = false)

proc outerB[V: static int; C: static int](s: StagGpu[V,float]; f: ptr UncheckedArray[float];
                                         x: array[C, ptr UncheckedArray[float]]; t: array[C, float];
                                         so: array[C, int]) =
  ## f_mu(y) += sum_j t_j x_j(y) x_j(y+mu)^+ over all sites, the remote
  ## neighbors of x_j through slot so_j of both exchanges
  getDefaultComm().barrier
  for q in 0..1:
    s.exB[q].packB(x, so)
    s.exB[q].start
  for q in 0..1: s.exB[q].wait
  let n = s.n
  let ne = s.ne
  let nl = s.nl
  let nb = s.nbr
  let ro0 = s.exB[0].rofs
  let rs0 = s.exB[0].rstr
  let rb0 = s.exB[0].rbuf
  let ro1 = s.exB[1].rofs
  let rs1 = s.exB[1].rstr
  let rb1 = s.exB[1].rbuf
  unpackB(x, C)
  unpackB(t, C)
  unpackB(so, C)
  gpuFor(i, 4*n):
    let mu = i div n
    let k = i - mu*n
    let jj = int nb[mu*n + k]
    let j = if nl == 12 and jj < 0: -1-jj else: jj
    let yo = vo(V, k)
    var acc {.noInit.}: array[18, float]
    forStatic e, 0, 17: acc[e] = 0.0
    forStatic js, 0, C-1:
      let xj = pickS(x, js)
      let tj = pickS(t, js)
      var v {.noInit.}: array[6, float]
      if j < n:
        let oj = vo(V, j)
        forStatic c, 0, 5: v[c] = xj[oj + c*V]
      else:
        let p = j - n
        let c0 = 6*pickS(so, js)
        if k < ne:  # an odd neighbor
          let op = int ro1[p]
          let sk = int rs1[p]
          forStatic c, 0, 5: v[c] = rb1[op + (c0+c)*sk]
        else:
          let op = int ro0[p]
          let sk = int rs0[p]
          forStatic c, 0, 5: v[c] = rb0[op + (c0+c)*sk]
      forStatic a, 0, 2:
        forStatic b, 0, 2:
          acc[6*a+2*b] += tj*(xj[yo + 2*a*V]*v[2*b] + xj[yo + (2*a+1)*V]*v[2*b+1])
          acc[6*a+2*b+1] += tj*(xj[yo + (2*a+1)*V]*v[2*b] - xj[yo + 2*a*V]*v[2*b+1])
    let o = 18*mu*n + lo18(V, k)
    forStatic e, 0, 17: f[o + e*V] += acc[e]

proc outerM*[V: static int](s: StagGpu[V,float]; f: ptr UncheckedArray[float]; x: openArray[ptr UncheckedArray[float]];
                            t: openArray[float]) =
  ## outerM of the vectors x_j with t_j, nBatch in a kernel:
  ## f_mu(y) += sum_j t_j x_j(y) x_j(y+mu)^+.  Needs newStagGpu(batch = true).
  var k = 0
  while k < x.len:
    let nc = min(nBatch, x.len - k)
    var xg: array[nBatch, ptr UncheckedArray[float]]
    var tg: array[nBatch, float]
    for j in 0..<nc:
      xg[j] = x[k+j]
      tg[j] = t[k+j]
    forActive((1 shl nc) - 1):
      s.outerB(f, sel(C, xg, so), sel(C, tg, so), so)
    k += nc

proc hopOE*[V: static int](s: StagGpu[V,float]; y, x: openArray[ptr UncheckedArray[float]]; c: float) =
  ## y_j,o = c D_oe x_j,e for the vectors j, nBatch at a time, with the D of
  ## dslash, twice that of applyMfull; y_j,e unchanged.  Needs
  ## newStagGpu(batch = true).
  var k = 0
  while k < x.len:
    let nc = min(nBatch, x.len - k)
    var xg, yg: array[nBatch, ptr UncheckedArray[float]]
    var z, cg: array[nBatch, float]
    for j in 0..<nc:
      xg[j] = x[k+j]
      yg[j] = y[k+j]
      cg[j] = c
    getDefaultComm().barrier
    forActive((1 shl nc) - 1):
      let xs = sel(C, xg, so)
      let ys = sel(C, yg, so)
      s.exB[0].packB(xs, so)
      s.exB[0].start
      s.exB[0].wait
      s.dslashB(1, ys, xs, ys, s.exB[0].rbuf, sel(C, z, so), sel(C, cg, so), so, nil, nil, dot = false, send = false)
    k += nc
