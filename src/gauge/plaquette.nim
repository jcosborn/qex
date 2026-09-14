## P = sum_{x,mu>nu} Re tr U_mu(x) U_nu(x+mu) U_mu(x+nu)^dag U_nu(x)^dag.
## stapleSum(g, ds) = D^(ds.len) grad P(g), with no color/volume normalization.
## Each gauge bundle has one layout. Workspace order/action are fixed at construction.

import base, layout, field, maths, comms/halo

type PlaqWork*[L,F,T] = ref object
  ## Every evaluation rebinds h to its current inputs. Each workspace owns its
  ## halos; layout/map metadata are cached. Row zero is unused at order three.
  hl*: HaloLayout[L]
  hm*: HaloMap[L]
  h*: seq[seq[Halo[L,F,T]]]
  order*: int
  action*: bool

proc newPlaqWork*[F](f: F, order = 0, action = false): auto =
  if order < 0 or order > 3 or (action and order != 0):
    raise newException(ValueError, "plaquette workspace requires order 0..3, or order 0 for the action")
  type
    L = type(f.l)
    T = eval(F.type.index(int))
  let
    lo = f.l
    nd = lo.nDim
  var fw, bw = newSeq[int32](nd)
  var offsets: seq[seq[int32]]
  for mu in 0..<nd:
    fw[mu] = 1
    if not action: bw[mu] = 1
    var off = newSeq[int32](nd)
    off[mu] = 1
    offsets.add off
    if not action:
      off = newSeq[int32](nd)
      off[mu] = -1
      offsets.add off
      for nu in 0..<nd:
        if nu == mu: continue
        off = newSeq[int32](nd)
        off[mu] = 1
        off[nu] = -1
        offsets.add off
  let hl = haloLayout(lo, fw, bw)
  result = PlaqWork[L,F,T](hl: hl, hm: haloMap(hl, getDefaultComm(), offsets),
                          order: order, action: action)
  result.h.newSeq(order+1)
  for k in 0..order:
    if order == 3 and k == 0: continue
    result.h[k].newSeq(nd)
    for mu in 0..<nd:
      result.h[k][mu] = makeHalo(hl, f)

proc newOneOf*[L,F,T](w: PlaqWork[L,F,T]): PlaqWork[L,F,T] =
  let k = if w.order == 3: 1 else: 0
  newPlaqWork(w.h[k][0].field, w.order, w.action)

proc update[L,F,T](w: PlaqWork[L,F,T], g: openArray[F], ds: openArray[seq[F]]) =
  if ds.len != w.order:
    raise newException(ValueError, "plaquette seed count differs from workspace order")
  if g.len != w.hl.lo.nDim or g[0].l != w.hl.lo:
    raise newException(ValueError, "plaquette workspace requires its original gauge shape")
  for d in ds:
    if d[0].l != w.hl.lo:
      raise newException(ValueError, "plaquette seeds require the workspace layout")
  let comm = getDefaultComm()
  for k in 0..w.order:
    for mu in 0..<w.h[k].len:
      w.h[k][mu].field = if k == 0: g[mu] else: ds[k-1][mu]
      w.h[k][mu].update(w.hm, comm)

proc plaqSum*[L,F,T](w: PlaqWork[L,F,T], g: openArray[F]): float =
  if w.order != 0:
    raise newException(ValueError, "plaquette sum requires a workspace of order zero")
  let ds: seq[seq[F]] = @[]
  w.update(g, ds)
  let
    h = w.h[0]
    hl = w.hl
    nd = h.len
  var sums = newSeq[float](getMaxThreads())
  threads:
    var s = 0.0
    for i in h[0].field:
      for mu in 1..<nd:
        let im = hl.neighborFwd[mu][i]
        for nu in 0..<mu:
          let inn = hl.neighborFwd[nu][i]
          let a = h[mu][i] * h[nu][im]
          let b = h[nu][i] * h[mu][inn]
          s += simdSum(redot(a, b))
    sums[threadNum] = s
  for s in sums: result += s
  rankSum(result)

proc plaqSum*[F](g: openArray[F]): float =
  let w = newPlaqWork(g[0], action=true)
  w.plaqSum(g)

proc productJet(h: auto, mu, nu, i, j, k: int, aa, ca: static bool, order: static int): auto =
  mixin load1, adj
  type M = type(load1(h[0][0][0]))
  template factor(a, pos: untyped) =
    forStatic d, (if order == 3: 1 else: 0), order:
      when pos == 0:
        when aa: a[d] := h[d][nu][i].adj
        else: a[d] := h[d][nu][i]
      elif pos == 1:
        a[d] := h[d][mu][j]
      else:
        when ca: a[d] := h[d][nu][k].adj
        else: a[d] := h[d][nu][k]
  var s: M
  productJet(s, 3, order, factor)
  s

proc stapleSumImpl[L,F,T](w: PlaqWork[L,F,T], f: array|seq, order: static int) =
  let
    h = w.h
    hl = w.hl
    nd = hl.lo.nDim
    lo = hl.lo
  threads:
    for mu in 0..<nd:
      for i in lo:
        var s: type(load1(f[0][0]))
        s := 0
        let im = int(hl.neighborFwd[mu][i])
        for nu in 0..<nd:
          if nu == mu: continue
          let
            ip = int(hl.neighborFwd[nu][i])
            ib = int(hl.neighborBck[nu][i])
            imb = int(hl.neighborFwd[mu][ib])
          s += productJet(h, mu, nu, i, ip, im, false, true, order)
          s += productJet(h, mu, nu, ib, ib, imb, true, false, order)
        f[mu][i] := s

proc stapleSum*[L,F,T](w: PlaqWork[L,F,T], g: openArray[F], ds: openArray[seq[F]], f: array|seq) =
  ## Output storage must be disjoint from seeds and, below order three, g.
  if w.action:
    raise newException(ValueError, "staple sum requires a staple workspace")
  if f[0].l != w.hl.lo:
    raise newException(ValueError, "staple output requires the workspace layout")
  template disjoint(src: untyped) =
    for dst in f:
      for x in src:
        if dst.s.data == x.s.data:
          raise newException(ValueError, "staple output must not alias an input")
  if w.order < 3: disjoint(g)
  for d in ds: disjoint(d)
  w.update(g, ds)
  case w.order
  of 0: stapleSumImpl(w, f, 0)
  of 1: stapleSumImpl(w, f, 1)
  of 2: stapleSumImpl(w, f, 2)
  of 3: stapleSumImpl(w, f, 3)
  else: discard  # constructor validates the order

proc stapleSum*[F](g: openArray[F], ds: openArray[seq[F]], f: array|seq) =
  let nd = g.len
  if ds.len > 3:
    threads:
      for mu in 0..<nd:
        f[mu] := 0
    return
  let w = newPlaqWork(g[0], ds.len)
  w.stapleSum(g, ds, f)

proc stapleSum*[F](g: openArray[F], f: array|seq) =
  let ds: seq[seq[F]] = @[]
  stapleSum(g, ds, f)
