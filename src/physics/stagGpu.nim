## Staggered even-odd CG with the links, vectors and halo exchanges on the GPU.
##
## Device fields keep the qex SIMD layout of the gauge field, V lanes per
## outer site and even outer sites first: real c of site k at
## (k div V)*6V + c*V + k mod V, so neighboring threads read neighboring
## reals and the even part is the first 6*ne reals.  Links are blocked the
## same way with nl reals per link: 18, or rows 0 and 1 only, 12, for links
## s W with W in SU(3) and s = +-1, rebuilding row 2 = s conj(row 0 x row 1).
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
import qex
import physics/qcdTypes
import solvers/solverBase
import backend/accel
import comms/[halo, halogpu]
import gauge/gaugeGpu
import base/metaUtils
import times

const nRed = 512  # atomic slots per dot product
var hopSplit* = -1  ## CG second hop split around the exchange: 1 always, 0 never, -1 with off-node neighbors

type
  StagGpu*[V: static int; T] = object
    lo*: Layout[V]
    n*, ne*: int  # local sites, even sites
    nl*: int  # reals per link, 12 with row 2 rebuilt in the kernels, else 18
    lf*, lb*: ptr UncheckedArray[T]  # [4][n div V][nl][V]: U_mu(x), U_mu(x-mu)^+ or nil
    lh*: array[2, ptr UncheckedArray[T]]  # without lb, [nl][s.ex[1-q].nrecv]: U_mu(x-mu)^+ at the receive position of a remote x-mu, x of parity q
    nbr*: ptr UncheckedArray[int32]  # [8][n]: x+mu, then x-mu, with the link signs for nl = 12
    ex*: array[2, GpuHaloEx[T]]  # exchange of the sites of each parity
    red*: ptr UncheckedArray[float]  # [2][2*nRed]: partial dot products, two buffers
    hred*: ptr UncheckedArray[float]  # [2*nRed] pinned host copy of a buffer
    vec*: array[7, ptr UncheckedArray[T]]  # work vectors, 6*n reals each
    rh*, sh*: ptr UncheckedArray[T]  # CG halo copies of r and s, as s.ex[0].rbuf
    ord*: ptr UncheckedArray[int32]  # even sites, the ones with remote neighbors first, then the nin with local ones only
    nin*: int
    split*: bool  # second hop of the CG in two kernels around the exchange of t

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

proc newStagGpu*[V: static int; E](g: openArray[Field[V,E]]; T: typedesc; recon = true; fwd = -1): StagGpu[V,T] =
  ## g are the phased links, as for newStag.  With recon, if every link is
  ## s W with W in SU(3) and s = +-1, the device keeps rows 0 and 1 and the
  ## kernels rebuild row 2 = s conj(row 0 x row 1), with s in the sign of
  ## the neighbor index: j for s = 1, -1-j for s = -1.  With fwd = 1 the
  ## device keeps U_mu(x) only and a hop reads U_mu(x-mu) at x-mu, half the
  ## memory; fwd = -1 does so when both would take over 128 MB, as they
  ## would no longer stay in the L2 cache of a PVC tile between the hops.
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
  # U_mu(x-mu)^+ for fb = 1, and their signs in sg
  var lk = newSeq[float](2*nd*18*n)
  var sg = newSeq[int8](2*nd*n)
  var bad = 0
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
    var dp, dm, nn = 0.0
    for b in 0..2:  # xr + i xi = conj(row 0 x row 1)_b
      let b1 = (b+1) mod 3
      let b2 = (b+2) mod 3
      let xr = m[2*b1]*m[6+2*b2] - m[2*b1+1]*m[7+2*b2] - m[2*b2]*m[6+2*b1] + m[2*b2+1]*m[7+2*b1]
      let xi = m[2*b2]*m[7+2*b1] + m[2*b2+1]*m[6+2*b1] - m[2*b1]*m[7+2*b2] - m[2*b1+1]*m[6+2*b2]
      let pr = m[12+2*b] - xr
      let pi = m[13+2*b] - xi
      let mr = m[12+2*b] + xr
      let mi = m[13+2*b] + xi
      dp += pr*pr + pi*pi
      dm += mr*mr + mi*mi
      nn += m[12+2*b]*m[12+2*b] + m[13+2*b]*m[13+2*b]
    for e in 0..17: lk[18*i + e] = m[e]
    sg[i] = if dp <= dm: 1 else: -1
    if min(dp, dm) > 1e-24*nn: inc bad  # row 2 off by more than 1e-12
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
  var nbad = float(bad)
  c.allReduce(nbad)
  let nl = if recon and nbad == 0: 12 else: 18
  result.nl = nl
  let fw = fwd > 0 or fwd < 0 and 2*nd*nl*n*sizeof(T) > 128 shl 20
  var lf = newSeq[T](nd*nl*n)
  var lb = newSeq[T](if fw: 0 else: nd*nl*n)
  for mu in 0..<nd:
    for k in 0..<n:
      let o = nl*mu*n + uo(V, k, nl)
      for e in 0..<nl:
        lf[o + e*V] = T(lk[18*(mu*n + k) + e])
        if not fw: lb[o + e*V] = T(lk[18*((nd + mu)*n + k) + e])
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
            for e in 0..<nl: lh[e*nr + j-n] = T(lk[18*((nd + mu)*n + k) + e])
      result.lh[q] = lh.toDevice
  if nl == 12:
    for i in 0..<nb.len:
      if sg[i] < 0: nb[i] = -1 - nb[i]
  result.nbr = nb.toDevice
  toc("links")

  result.red = newSeq[float](4*nRed).toDevice
  result.hred = cast[ptr UncheckedArray[float]](gpuMallocHost(2*nRed*sizeof(float)))
  for i in 0..<result.vec.len:
    result.vec[i] = cast[ptr UncheckedArray[T]](gpuMalloc(6*n*sizeof(T)))
  let nh = max(1, 6*result.ex[0].nrecv)
  result.rh = cast[ptr UncheckedArray[T]](gpuMalloc(nh*sizeof(T)))
  result.sh = cast[ptr UncheckedArray[T]](gpuMalloc(nh*sizeof(T)))
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

template load(m, u, o, st, e: untyped; nl: static int) =
  ## m[6a+2b] (re), m[6a+2b+1] (im) = U_ab from u[o + (6a+2b)*st], ...; for
  ## nl = 12 rows 0 and 1 only, row 2 = e conj(row 0 x row 1)
  forStatic i, 0, nl-1: m[i] = u[o + i*st]
  when nl == 12:
    forStatic b, 0, 2:
      const b1 = (b+1) mod 3
      const b2 = (b+2) mod 3
      m[12+2*b] = e*(m[2*b1]*m[6+2*b2] - m[2*b1+1]*m[7+2*b2] - m[2*b2]*m[6+2*b1] + m[2*b2+1]*m[7+2*b1])
      m[13+2*b] = e*(m[2*b2]*m[7+2*b1] + m[2*b2+1]*m[6+2*b1] - m[2*b1]*m[7+2*b2] - m[2*b1+1]*m[6+2*b2])

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

template dots(i, d, y, yo, st, rs: untyped) =
  ## rs slots get y.y and d.y, y of the site at y[yo + c*st]
  var yy = 0.0
  var dy = 0.0
  forStatic c, 0, 5:
    let yc = float(y[yo + c*st])
    yy += yc*yc
    dy += float(d[c])*yc
  gpuAtomicAdd(rs, i mod nRed, yy)
  gpuAtomicAdd(rs, nRed + i mod nRed, dy)

proc dslash[V: static int; T](s: StagGpu[V,T]; q: int; d, x, y, rb: ptr UncheckedArray[T];
                              a, b: T; rs, rz: ptr UncheckedArray[float]; dot, send: static bool;
                              j0 = 0; m = -1; nowait: static bool = false) =
  ## d = a y + b D x on the sites of parity q, the remote sites of x from rb,
  ## laid out as the receive buffer of s.ex[1-q]; for q = 0 on the sites
  ## s.ord[j0 ..< j0+m], all by default.  With send, stores d in the send
  ## slots of s.ex[q]; with dot, sums y.y and d.y in rs and zeros rz.  With
  ## nowait, returns before the kernel completes; gpuWaitAsync waits for it.
  let n = s.n
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
      let p = j0 + i
      if p < nRed:  # buffer of the next sum
        rz[p] = 0.0
        rz[nRed + p] = 0.0
      dots(p, acc, y, yo, V, rs)
  template kern(nl: static int; fw: static bool) =
    when nowait:
      gpuForAsync(i, nk): body(i, nl, fw)
    else:
      gpuFor(i, nk): body(i, nl, fw)
  if s.nl == 12:
    if lb == nil: kern(12, true) else: kern(12, false)
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
  ## returning the global r_e.r_e and w_e.r_e summed in rs.  Starts the
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
    s.ex[1].wait
    toc("wait oe")
    s.dslash(0, w, t, r, s.ex[1].rbuf, T(4*m2), T(-1), rs, rz, dot = true, send = true)
  toc("dslash eo")
  s.ex[0].start
  toc("start eo")
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

proc convert[V: static int; T,U](s: StagGpu[V,T]; d: ptr UncheckedArray[T]; x: ptr UncheckedArray[U]) =
  ## d_e = x_e in the precision of d
  gpuFor(i, 6*s.ne): d[i] = T(x[i])

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
  let red = s.red
  gpuFor(k, 4*nRed): red[k] = 0.0
  var rs = s.red
  var rz = cast[ptr UncheckedArray[float]](addr s.red[2*nRed])
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
    s.ex[0].wait
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
  forLinks(g, mu, k):
    let ol = nl*mu*n + uo(V, k, nl)
    var x {.noInit.}, y {.noInit.}: array[18, float]
    mload(x, u, 18*mu*n + lo18(V, k), V)
    let sf = sg[mu*n + k]
    for e in 0..<nl: lf[ol + e*V] = T(sf*x[e])
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
      for e in 0..<nl: lb[ol + e*V] = T(y[e])
    else:
      let jj = int fn[(4+mu)*n + k]
      let j = if jj < 0: -1-jj else: jj
      if j >= n:
        if k < ne:
          for e in 0..<nl: lh0[e*nr0 + j-n] = T(y[e])
        else:
          for e in 0..<nl: lh1[e*nr1 + j-n] = T(y[e])

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

proc solveM*[V: static int](s: StagGpu[V,float]; x, b: ptr UncheckedArray[float]; m: float;
                            sp: var SolverParams; ss: ptr StagGpu[V,float32] = nil; r2in = 1e-6) =
  ## Solves (m + D/2) x = b for b_o = 0, as stag.solve, until the residual
  ## of the even system drops by sp.r2req, in mixed precision with ss:
  ##   x_e = A^-1 4m b_e,  x_o = -D_oe x_e/(2m)
  let c = getDefaultComm()
  let b4 = s.vec[5]
  gpuFor(i, 6*s.ne): b4[i] = 4.0*m*b[i]
  var b2 = s.redot(b4, b4)
  c.allReduce(b2)
  if ss == nil:
    let (itn, _) = s.cg(x, b4, m, sp.r2req*b2, sp.maxits, sp.verbosity)
    sp.iterations += itn
  else:
    let (itn, _, _) = cg(s, ss[], x, b4, m, sp.r2req*b2, r2in, sp.maxits, sp.verbosity)
    sp.iterations += itn
  c.barrier
  s.ex[0].pack(x)
  s.ex[0].start
  s.ex[0].wait
  s.dslash(1, x, x, x, s.ex[0].rbuf, 0.0, -0.5/m, nil, nil, dot = false, send = false)

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
