import base
#import stdUtils
#import profile
import maths
import layout
import layout/shifts
import gaugeUtils
import staples
import comms/halo
import std/[tables, bitops]
import base/alignedMem

func firstSetBit*(x: SomeInteger): int {.inline.} =  # bitops doesn't handle x==0 correctly
  result = 0
  if x != 0:
    result = bitops.firstSetBit(x)

type
  GaugeActionCoeffs* = object
    plaq*: float
    rect*: float
    pgm*: float
    adjplaq*: float

template requireFundamental(c: GaugeActionCoeffs, name: string) =
  if c.adjplaq != 0:
    raise newException(ValueError, name & " requires fundamental coefficients")

template requireAdjoint(c: GaugeActionCoeffs, name: string) =
  if c.rect != 0 or c.pgm != 0:
    raise newException(ValueError, name & " supports plaquette and adjoint-plaquette coefficients only")

proc `*`*(x: float, y: GaugeActionCoeffs): GaugeActionCoeffs =
  for r, v in fields(result, y):
    r = x * v

const
  C1Symanzik = -1.0/12.0  # tree-level
  C1Iwasaki = -0.331
  C1DBW2 = -1.4088

# for DBW2
#
# Table III in https://journals.aps.org/prd/pdf/10.1103/PhysRevD.54.1050
# beta11 7.986(13)  beta12 -0.9169(41)
# gives c1 = 1/(8+7.986/-0.9169) = -1.4088813767670574
# with standard error 0.08227016226543633 assuming no correlation
#
# https://journals.aps.org/prd/pdf/10.1103/PhysRevD.75.114501
# uses c1 = -1.4069
#
# https://journals.aps.org/prd/pdf/10.1103/PhysRevD.90.074502
# uses c1 = -1.4088
#
# The numbers are likely from truncations in the middle of the computation,
# r = beta12/beta11 = -0.9169/7.986 = -0.11481342349110946
# c1 = 1/(8+1/r) = 1/(8+1/-0.1148) = -1.4068627450980384
# c1 = 1/(8+1/r) = 1/(8+1/-0.114813) = -1.408817610680277

proc gaugeActRect*(beta:float, c1:float = C1Symanzik):GaugeActionCoeffs =
  result.plaq = (1.0-8.0*c1)*beta
  result.rect = c1*beta

proc Symanzik*(beta:float, c1:float = C1Symanzik):auto = gaugeActRect(beta, c1)
proc Iwasaki*(beta:float, c1:float = C1Iwasaki):auto = gaugeActRect(beta, c1)
proc DBW2*(beta:float, c1:float = C1DBW2):auto = gaugeActRect(beta, c1)

type
  LoopKind = enum lkPlaq, lkRect, lkPgm
  GaugeLoop = object
    path: array[5,int]
    len: int
    kind: LoopKind
  GaugeLoops = seq[seq[GaugeLoop]]
  GaugeLoopPlan = ref object
    paths: GaugeLoops
    fw, bw: seq[int32]
    offsets: seq[seq[seq[int32]]]  # base sites of the links of each direction
    products: array[6,LoopProducts]
    factors: seq[tuple[mu:int, offset:seq[int32]]]
    pairs: PathPlan
    outputs: seq[tuple[mu:int, kind:LoopKind]]
  LoopLoad = object
    dest, mu, seed: int
    walk: seq[tuple[d, n: int]]  # nonzero components of the factor's offset
  LoopMul = object
    dest, left, right: int
    la, ra, add: bool
  LoopOut = object
    source, mu: int
    kind: LoopKind
    adj, first: bool  # first: the first output of mu in sched
  LoopProducts = ref object
    loads: seq[LoopLoad]
    steps: seq[LoopMul]
    outs: seq[LoopOut]
    size: int  # values: loads and products
    sched: seq[int]  # k >= 0: steps[k]; k < 0: outs[-k-1], each after its source
    slot: seq[int]  # scratch slot of each product, -1 for a load
    slots: int

proc gaugeLoops(nd, mask, act: int): GaugeLoopPlan =
  # Each open path runs from x to x+mu. Differentiating a closed loop
  # contributes every occurrence, including repeated links after wrapping.
  proc makeLoops(): GaugeLoops =
    result.newSeq(nd)
    template addLoop(p: untyped, family: LoopKind) =
      block:
        let path = p
        for k in 0..(if act == 1: 0 else: path.len-1):
          var loop = GaugeLoop(len: path.len-1, kind: family)
          if path[k] > 0:
            for j in 0..<loop.len:
              loop.path[j] = -path[(k+path.len-1-j) mod path.len]
          else:
            for j in 0..<loop.len:
              loop.path[j] = path[(k+1+j) mod path.len]
          result[abs(path[k])-1].add loop
    for mu in 2..nd:
      for nu in 1..<mu:
        if (mask and 4) != 0:
          addLoop([mu,nu,-mu,-nu], lkPlaq)
        if (mask and 1) != 0:
          addLoop([mu,nu,nu,-mu,-nu,-nu], lkRect)
          addLoop([mu,mu,nu,-mu,-mu,-nu], lkRect)
        if (mask and 2) != 0:
          for sg in 1..<nu:
            # gaugeAction2's ts1, ts2, ts3, ts7, respectively.
            addLoop([mu,nu,sg,-mu,-nu,-sg], lkPgm)
            addLoop([mu,sg,nu,-mu,-sg,-nu], lkPgm)
            addLoop([nu,mu,sg,-nu,-mu,-sg], lkPgm)
            addLoop([mu,-nu,sg,-mu,nu,-sg], lkPgm)
  # Paths and offsets depend only on topology. Construct this metadata and its
  # derived product plans outside threads; evaluations share the cached plans.
  var cache {.global.}: Table[(int,int,int),GaugeLoopPlan]
  let key = (nd, mask, act)
  if key notin cache:
    let plan = GaugeLoopPlan(paths: makeLoops(), fw: newSeq[int32](nd), bw: newSeq[int32](nd),
                             offsets: newSeq[seq[seq[int32]]](nd))
    for row in plan.paths:
      for path in row:
        var x = newSeq[int32](nd)
        for k in 0..<path.len:
          let step = path.path[k]
          let mu = abs(step)-1
          if step < 0: dec x[mu]
          var off = newSeq[int32](nd)
          for d in 0..<nd:
            off[d] = x[d]
            plan.fw[d] = max(plan.fw[d], x[d])
            plan.bw[d] = max(plan.bw[d], -x[d])
          if off notin plan.offsets[mu]: plan.offsets[mu].add off
          if step > 0: inc x[mu]
    cache[key] = plan
  cache[key]

type
  LoopKeep* = enum
    ## What gaugeAction1 and gaugeActionDeriv leave in a LoopWork between calls.
    keepNone = "none"  # staple blocks allocated and freed per call
    keepStaples = "staples"  # stf, stu and ru kept: 2 or 3 arrays of nd*(nd-1) fields
  LoopTag* = tuple[id, ver: int]
    ## Identity and version of a link value, given by a caller that knows them
    ## (graph node id and epoch, a driver's configuration counter). id 0 is no tag.
  LoopHalos[V:static[int],T] = object
    layout*: HaloLayout[Layout[V]]
    maps*: seq[HaloMap[Layout[V]]]  # per direction: the sites the links of the last update reach
    fields*: seq[seq[Halo[Layout[V],Field[V,T],T]]]
    keys: seq[tuple[tag: LoopTag, plan: GaugeLoopPlan]]  # each slot's last update
  LoopWork*[V:static[int],T] = ref object
    ## Mutable loop buffers for one layout. Serialize calls using this object.
    ## Geometry initialization follows the halo API's threading requirements.
    lo*: Layout[V]
    halos*: seq[LoopHalos[V,T]]
    keep: LoopKeep
    stf, stu: FieldArray[V,T]  # staples of the last action/force call, passed back to makeStaples
    haveBackward: bool  # stu is current at stag; stf is always current there
    stag: LoopTag  # links of the kept staples
    spair: bool  # a staple force followed an action at the same tag: actions build stu too
    ltag: LoopTag  # links of the last loop call, independent of staple contents
    lpair: bool  # a loop derivative followed an action: actions fill the wider halo
    builds*, reuses*, updates*: int  # staple builds, reused builds, halo updates
    ru: FieldArray[V,T]  # rect parts of gaugeActionDeriv, zeroed on every use
    scratch: seq[Field[V,T]]
    trans: array[4,seq[Transporter[Field[V,T],Field[V,T],T]]]
    plan: LoopProducts
    tile: int  # sites per evaluation step of the loop kernels
    products: seq[alignedMem[T]]  # per thread: product slots, then nd sums, tile sites each
    ptrs: seq[alignedMem[ptr T]]  # per thread: each plan value at the tile's sites
    vals: seq[alignedMem[float]]  # per thread: action terms of the tile's sites

proc newLoopWork*[V:static[int],T](f: Field[V,T], keep = keepNone): LoopWork[V,T] =
  ## keep: see LoopKeep. Halos, scratch, Hessian transporters and loop products
  ## are made on first use and persist in every mode.
  LoopWork[V,T](lo: f.l, keep: keep)

proc bytes*(w: LoopWork): int =
  ## Owned raw-memory allocations; excludes borrowed inputs, shift buffers and metadata.
  template addArray(a: untyped) =
    if a.mem != nil: result += a.mem.size
  addArray(w.stf)
  addArray(w.stu)
  addArray(w.ru)
  for f in w.scratch:
    if f != nil and f.s.mem != nil: result += f.s.mem.size
  for row in w.trans:
    for t in row:
      if t.field != nil and t.field.s.mem != nil: result += t.field.s.mem.size
  for row in w.halos:
    for slot in row.fields:
      for h in slot:
        if h != nil and h.halo.mem != nil: result += h.halo.mem.size
  for p in w.products:
    if p.mem != nil: result += p.mem.size
  for p in w.ptrs:
    if p.mem != nil: result += p.mem.size
  for p in w.vals:
    if p.mem != nil: result += p.mem.size

proc useLoopWork(g: auto, work: LoopWork): auto =
  result = if work.isNil: newLoopWork(g[0]) else: work
  if result.lo != g[0].l:
    raise newException(ValueError, "loop workspace requires its original layout")

proc loopScratch(w: LoopWork): auto =
  if w.scratch.len == 0:
    w.scratch.newSeq(w.lo.nDim)
    for mu in 0..<w.scratch.len: w.scratch[mu].new(w.lo)
  w.scratch

proc hessianTrans(w: LoopWork, g, h: auto): auto =
  if g[0].l != w.lo or h[0].l != w.lo:
    raise newException(ValueError, "Hessian workspace requires its original layout")
  if w.trans[0].len == 0:
    w.trans[0] = newTransporters(g, g[0], 1)
    w.trans[1] = newTransporters(g, g[0], -1)
    w.trans[2] = newTransporters(h, g[0], 1)
    w.trans[3] = newTransporters(h, g[0], -1)
  else:
    w.trans[0].setLinks(g)
    w.trans[1].setLinks(g)
    w.trans[2].setLinks(h)
    w.trans[3].setLinks(h)
  (w.trans[0], w.trans[1], w.trans[2], w.trans[3])

proc clearInputs(w: LoopWork) =
  for row in w.trans.mitems: row.clearLinks
  for row in w.halos:
    for slot in row.fields:
      for h in slot:
        if h != nil: h.field = nil

proc gaugeLoopHalo[V:static[int],T](w: LoopWork[V,T], g: openArray[Field[V,T]], plan: GaugeLoopPlan, slot: int, tag = default(LoopTag)): auto =
  if g[0].l != w.lo:
    raise newException(ValueError, "loop inputs require the workspace layout")
  let
    lo = g[0].l
    nd = lo.nDim
    comm = getDefaultComm()
    hl = haloLayout(lo, plan.fw, plan.bw)
  var row = 0
  while row < w.halos.len and w.halos[row].layout != hl: inc row
  if row == w.halos.len: w.halos.add LoopHalos[V,T](layout: hl, maps: newSeq[HaloMap[Layout[V]]](nd))
  # Plans with the same extents share a row; each direction's map holds only
  # the sites that direction's links reach.
  for mu in 0..<nd: w.halos[row].maps[mu] = haloMap(hl, comm, plan.offsets[mu])
  if w.halos[row].fields.len <= slot:
    w.halos[row].fields.setLen(max(2,slot+1))
    w.halos[row].keys.setLen(w.halos[row].fields.len)
  if w.halos[row].fields[slot].len == 0:
    w.halos[row].fields[slot].newSeq(nd)
    for mu in 0..<nd:
      w.halos[row].fields[slot][mu] = makeHalo(hl, g[mu])
  let gh = w.halos[row].fields[slot]
  let same = tag.id != 0 and w.halos[row].keys[slot] == (tag, plan)
  for mu in 0..<nd:
    gh[mu].field = g[mu]
    if not same: gh[mu].update(w.halos[row].maps[mu], comm)
  if same: inc w.reuses
  else: w.updates += nd
  w.halos[row].keys[slot] = (tag, plan)
  gh

proc gaugeLoopProd(gh, hh: auto, path: GaugeLoop, i: int, deriv: static bool): auto =
  mixin adj, load1
  type M = type(load1(gh[0][0]))
  const order = if deriv: 1 else: 0
  var p: M
  var j = i
  template factor(a, k: untyped) =
    let step = path.path[k]
    let mu = abs(step)-1
    if step < 0: j = gh[mu].layout.neighborBck[mu][j]
    forStatic d, 0, order:
      when d == 0:
        if step > 0: a[d] := gh[mu][j]
        else: a[d] := gh[mu][j].adj
      else:
        if step > 0: a[d] := hh[mu][j]
        else: a[d] := hh[mu][j].adj
    if step > 0 and k < path.path.high: j = gh[mu].layout.neighborFwd[mu][j]
  productJet(p, path.path.len, order, factor)
  p

proc pgmAction[T](c: GaugeActionCoeffs, g: openArray[T], work: auto): float =
  let nd = g[0].l.nDim
  if nd < 3 or c.pgm == 0: return
  let w = useLoopWork(g, work)
  defer: w.clearInputs
  let
    plan = gaugeLoops(nd, 2, 1)
    gh = w.gaugeLoopHalo(g, plan, 0)
  var sums = newSeq[float](getMaxThreads())
  threads:
    var s = 0.0
    for mu in 0..<nd:
      for i in gh[0].field:
        for path in plan.paths[mu]:
          let p = gaugeLoopProd(gh, gh, path, i, false)
          s += simdSum(redot(gh[mu][i], p))
    sums[threadNum] = s
  for s in sums: result += s
  rankSum(result)
  result *= -c.pgm / float(g[0][0].nrows)

proc gaugeLoopDeriv[T](c: GaugeActionCoeffs, g, h: openArray[T], f: array|seq, deriv: static bool, work: auto) =
  # d(prod A_k) = sum_k A_0...dA_k...A_4, with no seed storage in
  # the ordinary gradient. All target links are local, so writes are disjoint.
  let nd = g[0].l.nDim
  let mask = (if c.rect != 0: 1 else: 0) + (if c.pgm != 0 and nd >= 3: 2 else: 0)
  if mask == 0: return
  let w = useLoopWork(g, work)
  defer: w.clearInputs
  let
    plan = gaugeLoops(nd, mask, 0)
    gh = w.gaugeLoopHalo(g, plan, 0)
    nc = float(g[0][0].nrows)
  when deriv:
    let hh = w.gaugeLoopHalo(h, plan, 1)
  threads:
    for mu in 0..<nd:
      for i in gh[0].field:
        var s: type(load1(gh[0][0]))
        s := 0
        for path in plan.paths[mu]:
          let a = (if path.kind == lkRect: c.rect else: c.pgm) / nc
          when deriv:
            s += a*gaugeLoopProd(gh, hh, path, i, true)
          else:
            s += a*gaugeLoopProd(gh, gh, path, i, false)
        f[mu][i] += s

proc loopMask(c: GaugeActionCoeffs, nd: int): int =
  (if c.plaq != 0 and nd > 1: 4 else: 0) or
  (if c.rect != 0 and nd > 1: 1 else: 0) or
  (if c.pgm != 0 and nd > 2: 2 else: 0)

proc loopProducts(plan: GaugeLoopPlan, order: int): LoopProducts =
  # A symbol identifies U_mu(x+offset); its sign denotes adjoint. All symbols
  # already refer to the same output site, so the pair planner needs no shifts.
  if plan.products[order] != nil: return plan.products[order]
  let p = LoopProducts()
  if plan.pairs.outs.len == 0:
    var
      symbols: Table[(int,seq[int32]),int]
      words: seq[seq[int]]
    for mu, row in plan.paths:
      for path in row:
        var x = newSeq[int32](plan.fw.len)
        var word: seq[int]
        for k in 0..<path.len:
          let step = path.path[k]
          let dir = abs(step)-1
          if step < 0: dec x[dir]
          var off = newSeq[int32](x.len)
          for d in 0..<x.len: off[d] = x[d]
          let symbol = (dir,off)
          if symbol notin symbols:
            plan.factors.add (dir,off)
            symbols[symbol] = plan.factors.len
          let id = symbols[symbol]
          word.add(if step > 0: id else: -id)
          if step > 0: inc x[dir]
        words.add word
        plan.outputs.add (mu,path.kind)
    plan.pairs = words.optimalPairs.plan(shifts=false)
  let pairs = plan.pairs
  var nodes: Table[seq[int],PathStep]
  for step in pairs.steps: nodes[step.key] = step
  var jets: Table[(seq[int],int),int]
  proc jet(word: seq[int], mask: int): int =
    let key = (word,mask)
    if key in jets: return jets[key]
    if word.len == 1:
      let f = plan.factors[word[0]-1]
      result = p.size
      inc p.size
      p.loads.add LoopLoad(dest:result,mu:f.mu,seed:firstSetBit(mask))
      for d, n in f.offset:
        if n != 0: p.loads[^1].walk.add (d, int(n))
    else:
      # (AB)_S = sum_{T subset S} A_T B_{S\T}. Only reachable coefficients
      # are built, including at the top degree where no primal loads remain.
      let step = nodes[word]
      var terms: seq[LoopMul]
      var sub = mask
      while true:
        if countSetBits(sub) <= step.l.len and countSetBits(mask xor sub) <= step.r.len:
          let left = jet(step.l,sub)
          let right = jet(step.r,mask xor sub)
          terms.add LoopMul(left:left,right:right,la:step.la,ra:step.ra,add:terms.len>0)
        if sub == 0: break
        sub = (sub-1) and mask
      result = p.size
      inc p.size
      for term in terms.mitems:
        term.dest = result
        p.steps.add term
    jets[key] = result
  for i, outp in pairs.outs:
    p.outs.add LoopOut(source:jet(outp.key,(1 shl order)-1),
      mu:plan.outputs[i].mu,kind:plan.outputs[i].kind,adj:outp.adj)
  # Each output follows the last step of its source and the previous output,
  # so sums keep the order of outs. Every step allocates its new destination
  # before releasing operands at their last use: a destination never aliases
  # an operand, and slots hold only live products.
  var done = newSeq[int](p.size)
  for v in done.mitems: v = -1
  for k, op in p.steps: done[op.dest] = k
  var next = 0
  # Active derivative plans cover every direction; empty/above-degree cases
  # return before execution. Each direction starts with its first output.
  var seen = newSeq[bool](plan.paths.len)
  for k, op in p.outs:
    while next <= done[op.source]:
      p.sched.add next
      inc next
    p.outs[k].first = not seen[op.mu]
    seen[op.mu] = true
    p.sched.add -k-1
  var last = newSeq[int](p.size)
  for pos, k in p.sched:
    if k >= 0:
      last[p.steps[k].left] = pos
      last[p.steps[k].right] = pos
    else: last[p.outs[-k-1].source] = pos
  p.slot = newSeq[int](p.size)
  for v in p.slot.mitems: v = -1
  var free: seq[int]
  template release(v: int, pos: int) =
    if p.slot[v] >= 0 and last[v] == pos and p.slot[v] notin free: free.add p.slot[v]
  for pos, k in p.sched:
    if k >= 0:
      let op = p.steps[k]
      if p.slot[op.dest] < 0:
        if free.len > 0: p.slot[op.dest] = free.pop
        else:
          p.slot[op.dest] = p.slots
          inc p.slots
      release(op.left, pos)
      release(op.right, pos)
    else: release(p.outs[-k-1].source, pos)
  plan.products[order] = p
  p

proc prepareProducts[V:static[int],T](w: LoopWork[V,T], plan: LoopProducts, action = false) =
  ## Per thread: a tile of 1 to 32 sites whose slots and sums take about 32 KB
  ## and whose value pointers number at most 8192 (Wilson force: 1 SU(3) site
  ## at 16 lanes, 6 at 4 lanes; 32 U(1) sites). The buffers are allocated for
  ## these bounds (and 32 sites of action terms), so a plan grows them only when
  ## a single site of it needs more, and a smaller plan after a larger one
  ## allocates nothing.
  w.plan = plan
  let m = plan.slots+w.lo.nDim
  let n = plan.size
  let o = plan.outs.len
  w.tile = max(1, min(32, min(32768 div (m*sizeof(T)), 8192 div n)))
  let nt = getMaxThreads()
  if w.products.len < nt:
    w.products.setLen(nt)
    w.ptrs.setLen(nt)
    w.vals.setLen(nt)
  for tid in 0..<nt:
    if w.products[tid].len < m*w.tile:
      w.products[tid].newU(max(m*w.tile, 32768 div sizeof(T)))
    if w.ptrs[tid].len < n*w.tile:
      w.ptrs[tid].newU(max(n*w.tile, 8192))
    if action and w.vals[tid].len < o*w.tile:
      w.vals[tid].newU(32*o)

template tileSites(n, tile: int; i0, nb, body: untyped) =
  ## The thread's share of outer sites 0..<n, as a field loop assigns them, in
  ## tiles of nb <= tile sites from i0.
  let ti1 = ((threadNum+1)*n) div numThreads
  var i0 = (threadNum*n) div numThreads
  while i0 < ti1:
    let nb = min(tile, ti1-i0)
    body
    i0 += nb

template evalProducts(plan: LoopProducts, h, hl: auto, i0, nb, tile: int, ptrs: auto, outp, src, b, body: untyped) =
  ## Products of plan at the outer sites i0+b, b < nb, each step for the whole
  ## tile; body runs for each output outp and site b, in the order of
  ## plan.outs, with src its value. Value v of site b is ptrs[v*tile+b]: loads
  ## point into the halos (layout hl), products into the thread's slots (bindSlots).
  mixin adj
  block:
    for op in plan.loads:
      for b in 0..<nb:
        var j = i0+b
        for st in op.walk: j = hl.step(st.d, j, st.n)
        ptrs[op.dest*tile+b] = h[op.seed][op.mu].indexPtr(j)
    for k in plan.sched:
      if k >= 0:
        let op = plan.steps[k]
        let r = op.dest*tile
        let x = op.left*tile
        let y = op.right*tile
        template run(xa, ya, ad: static bool) =
          for e in 0..<nb:
            let rp = ptrs[r+e]
            let xp = ptrs[x+e]
            let yp = ptrs[y+e]
            template mulStep(xv, yv: untyped) =
              when ad:
                var t {.noinit.}: type(rp[])
                mul(t, xv, yv)
                rp[] += t
              else: mul(rp[], xv, yv)
            when xa and ya: mulStep(xp[].adj, yp[].adj)
            elif xa: mulStep(xp[].adj, yp[])
            elif ya: mulStep(xp[], yp[].adj)
            else: mulStep(xp[], yp[])
        case 4*ord(op.la)+2*ord(op.ra)+ord(op.add)
        of 0: run(false, false, false)
        of 1: run(false, false, true)
        of 2: run(false, true, false)
        of 3: run(false, true, true)
        of 4: run(true, false, false)
        of 5: run(true, false, true)
        of 6: run(true, true, false)
        else: run(true, true, true)
      else:
        let outp = plan.outs[-k-1]
        let o = outp.source*tile
        for b in 0..<nb:
          template src: untyped = ptrs[o+b][]
          body

template bindSlots(plan: LoopProducts, sc, ptrs: auto, tile: int) =
  ## Point every product of plan at its scratch slots.
  for v, s in plan.slot:
    if s >= 0:
      for b in 0..<tile: ptrs[v*tile+b] = addr sc[s*tile+b]

template addOut(r: untyped, a: float, outp, src: untyped) =
  ## r = a src, or r += a src after the first output of the direction.
  if outp.first:
    if outp.adj: r := a*src.adj
    else: r := a*src
  else:
    if outp.adj: r += a*src.adj
    else: r += a*src

proc loopAction*[V:static[int],T](c: GaugeActionCoeffs, g: openArray[Field[V,T]]; work: LoopWork[V,T] = nil, tag = default(LoopTag)): float =
  ## Fused S = -sum(c_loop Re tr U_loop)/Nc for fundamental loop terms.
  ## tag: identity and version of g. A tagged call at the links of the kept halo
  ## skips the exchange; once a loopDeriv has followed an action at the same tag,
  ## tagged actions fill the derivative's wider halo so that the derivative reuses it.
  requireFundamental(c, "loopAction")
  let nd = g[0].l.nDim
  if g.len != nd: raise newException(ValueError, "loop action requires a complete gauge bundle")
  let mask = loopMask(c,nd)
  if mask == 0: return
  let w = useLoopWork(g,work)
  defer: w.clearInputs
  let plan = gaugeLoops(nd,mask,1)
  let products = loopProducts(plan,0)
  w.prepareProducts(products, action=true)
  let hplan = if tag.id != 0 and w.lpair: gaugeLoops(nd,mask,0) else: plan
  let h = [w.gaugeLoopHalo(g,hplan,0,tag)]
  w.ltag = tag
  var sums = newSeq[array[3,float]](getMaxThreads())
  let u = cast[ptr cArray[Field[V,T]]](unsafeAddr(g[0]))
  let hl = h[0][0].layout
  let tile = w.tile
  threads:
    let ptrs = w.ptrs[threadNum].data
    let vals = w.vals[threadNum].data
    products.bindSlots(w.products[threadNum].data, ptrs, tile)
    var s: array[3,float]
    tileSites(u[0].l.nSitesOuter, tile, i0, nb):
      var o = 0
      products.evalProducts(h,hl,i0,nb,tile,ptrs,outp,src,b):
        vals[o*tile+b] =
          if outp.adj: simdSum(redot(u[outp.mu][i0+b],src.adj))
          else: simdSum(redot(u[outp.mu][i0+b],src))
        if b == nb-1: inc o
      # Sum in site order, as one site at a time.
      for b in 0..<nb:
        for o, outp in products.outs: s[ord(outp.kind)] += vals[o*tile+b]
    sums[threadNum] = s
  var total: array[3,float]
  for s in sums:
    for k in 0..<3: total[k] += s[k]
  rankSum(total)
  -(c.plaq*total[0]+c.rect*total[1]+c.pgm*total[2])/float(g[0][0].nrows)

proc loopDeriv*[V:static[int],T,H; K:static int](c: GaugeActionCoeffs, g:openArray[Field[V,T]], ds:array[K,H], f:array|seq; work:LoopWork[V,T] = nil, tag = default(LoopTag)) =
  ## f = D^K grad S(g)[ds]. The seed count is static; orders above degree vanish.
  ## Output must be disjoint from seeds and, below the polynomial degree, g.
  ## tag: as in loopAction; seeds are always exchanged.
  requireFundamental(c, "loopDeriv")
  let lo = g[0].l
  let nd = lo.nDim
  if g.len != nd or f.len != nd or f[0].l != lo:
    raise newException(ValueError, "loop derivative requires matching gauge bundles")
  for d in ds:
    if d.len != nd or d[0].l != lo:
      raise newException(ValueError, "loop derivative seeds require matching gauge bundles")
  var mask = loopMask(c,nd)
  when K > 3: mask = mask and 3
  let degree = if (mask and 3) != 0: 5 else: 3
  if mask == 0 or K > degree:
    threads:
      for mu in 0..<nd: f[mu] := 0
    return
  template disjoint(src: untyped) =
    for a in f:
      for b in src:
        if a.s.data == b.s.data:
          raise newException(ValueError, "loop derivative output must not alias an input")
  if K < degree: disjoint(g)
  for d in ds: disjoint(d)
  let w = useLoopWork(g,work)
  defer: w.clearInputs
  let plan = gaugeLoops(nd,mask,0)
  let products = loopProducts(plan,K)
  w.prepareProducts(products)
  type HaloRow = type(w.gaugeLoopHalo(g,plan,0))
  var h: array[K+1,HaloRow]
  if tag.id != 0 and tag == w.ltag: w.lpair = true
  if K < degree: h[0] = w.gaugeLoopHalo(g,plan,0,tag)
  for k in 0..<K: h[k+1] = w.gaugeLoopHalo(ds[k],plan,k+1)
  w.ltag = tag
  let nc = float(g[0][0].nrows)
  let a = [-c.plaq/nc, -c.rect/nc, -c.pgm/nc]
  let hl = h[K][0].layout
  let tile = w.tile
  threads:
    let sc = w.products[threadNum].data
    let ptrs = w.ptrs[threadNum].data
    products.bindSlots(sc, ptrs, tile)
    let d = products.slots*tile  # sum of direction mu at site b: sc[d+mu*tile+b]
    tileSites(f[0].l.nSitesOuter, tile, i0, nb):
      products.evalProducts(h,hl,i0,nb,tile,ptrs,outp,src,b):
        addOut(sc[d+outp.mu*tile+b], a[ord(outp.kind)], outp, src)
      for b in 0..<nb:
        for mu in 0..<nd: f[mu][i0+b] := sc[d+mu*tile+b]

proc loopDeriv*[V:static[int],T](c: GaugeActionCoeffs, g:openArray[Field[V,T]], f:array|seq; work:LoopWork[V,T] = nil, tag = default(LoopTag)) =
  let ds = default(array[0,seq[Field[V,T]]])
  c.loopDeriv(g,ds,f,work,tag)

# plaq: 6 types
# rect: 12 types
# pgm: 32=4*2*4=4*3*2+4*2 types, via pgmAction
# shift corners: u[mu],nu mu != nu (12)
# make staples: s[mu][nu] mu != nu (12)
# plaq traces:
#  plaq: U[mu]^+ * sum_{nu!=mu} s[mu][nu] (6)
#  rect: shift(s[mu][nu], nu) (12 shifts)
proc gaugeAction1*[V:static[int],T](c: GaugeActionCoeffs, uu: openArray[Field[V,T]]; work: LoopWork[V,T] = nil, tag = default(LoopTag)): auto =
  ## tag: identity and version of uu. With keepStaples, a call whose tag equals
  ## the tag of the kept staples skips the link exchange and the staple build;
  ## a rect term still exchanges stu. Unretained rectangle calls reuse a bounded
  ## scratch batch across planes and release it before returning.
  requireFundamental(c, "gaugeAction1")
  if work != nil: discard useLoopWork(uu, work)
  defer:
    if work != nil: work.clearInputs
  tic("gaugeAction1")
  let lo = uu[0].l
  let nd = lo.nDim
  let nc = float(uu[0][0].ncols)
  let keep = work != nil and work.keep == keepStaples
  if not keep and c.rect != 0:
    result = stapleBatches(c.plaq, c.rect, uu, default(seq[Field[V,T]]), false)
    if work != nil: inc work.builds
    if c.pgm != 0: result += pgmAction(c, uu, work)
    toc("gaugeAction end")
    return
  let reusing = keep and tag.id != 0 and tag == work.stag and (c.rect == 0 or work.haveBackward)
  var
    stf, stu: FieldArray[V,T]
    ss: seq[seq[ShiftB[T]]]
  if reusing:
    stf = work.stf
    stu = work.stu
    if c.rect != 0:
      ss = newLinkShifts(lo, T, -1)
      threads:
        for mu in 1..<nd:
          for nu in 0..<mu:
            ss[mu][nu].startSB(stu[mu,nu][ix])
            ss[nu][mu].startSB(stu[nu,mu][ix])
    inc work.reuses
    toc("gaugeAction reuse")
  else:
    let cs = startCornerShifts(uu)
    toc("gaugeAction startCornerShifts")
    # Once a force has followed an action at the same links, tagged actions keep
    # both staple sets so that the next such force reuses them. Only a rect
    # term reads the backward set here; the force exchanges it itself.
    let both = c.rect != 0 or (keep and tag.id != 0 and work.spair)
    if not both:
      stf = if work == nil: makeFwdStaples(uu, cs)
        else: makeFwdStaples(uu, cs, work.stf)
    else:
      (stf, stu, ss) = if work == nil: makeStaples(uu, cs, post = c.rect != 0)
        else: makeStaples(uu, cs, work.stf, work.stu, c.rect != 0)
    if work != nil:
      inc work.builds
      work.haveBackward = both
      work.stag = tag
      if keep:
        work.stf = stf
        if both: work.stu = stu
    toc("gaugeAction makeStaples")
  var sums = newSeq[array[2,float]](getMaxThreads())
  addStapleActions(uu, stf, stu, ss, staplePlanes(nd), c.rect != 0, sums)
  var a = [0.0, 0.0]
  for part in sums:
    a[0] += part[0]
    a[1] += part[1]
  rankSum(a)
  result = -(c.plaq*a[0] + c.rect*a[1])/nc
  if c.pgm != 0:
    result += pgmAction(c, uu, work)
  if not keep:
    release(stf)
    release(stu)
  toc("gaugeAction end")

proc gaugeAction1*[T](uu: openarray[T]): auto =
  let gc = GaugeActionCoeffs(plaq:1.0)
  return gc.gaugeAction1(uu)

proc gaugeActionDeriv*[V:static[int],T](c: GaugeActionCoeffs, uu: openArray[Field[V,T]], f: array|seq, accumulate=false; work: LoopWork[V,T] = nil, tag = default(LoopTag)) =
  ## if accumulate, the derivatives will add to f.
  ## Otherwise the derivative overwrites f.
  ## tag: as in gaugeAction1; the kernel reads shifted links, so a reused build
  ## still exchanges the links and stu, and stf and ru for a rect term.
  requireFundamental(c, "gaugeActionDeriv")
  if work != nil: discard useLoopWork(uu, work)
  defer:
    if work != nil: work.clearInputs
  mixin load1, adj
  tic("gaugeActionDeriv")
  let u = cast[ptr cArray[Field[V,T]]](unsafeAddr(uu[0]))
  let lo = u[0].l
  let nd = lo.nDim
  #let np = (nd*(nd-1)) div 2
  let nc = u[0][0].ncols
  let cp = c.plaq / float(nc)
  let cr = c.rect / float(nc)
  let keep = work != nil and work.keep == keepStaples
  if not keep:
    discard stapleBatches(c.plaq, c.rect, uu, f, true, accumulate)
    if work != nil: inc work.builds
    if c.pgm != 0:
      gaugeLoopDeriv(GaugeActionCoeffs(pgm: c.pgm), uu, uu, f, false, work)
    toc("end")
    return
  let reusing = keep and tag.id != 0 and tag == work.stag and work.haveBackward
  if keep and tag.id != 0 and tag == work.stag: work.spair = true
  var ru:FieldArray[V,T]  # the rect parts of 3
  var sb:seq[seq[ShiftB[T]]]  # backward ru
  var sf:seq[seq[ShiftB[T]]]  # forward stf
  if cr!=0:
    # Temporary rectangle work returned above; this path retains its scratch.
    if work.ru.shape.len == 0: work.ru = newFieldArray2(lo,type(u[0]),[nd,nd],mu!=nu)
    ru = work.ru
    sb = newLinkShifts(lo, T, -1)
    sf = newLinkShifts(lo, T, 1)
  toc("init")
  let cs = startCornerShifts(uu)
  var stf, stu: FieldArray[V,T]
  var ss: seq[seq[ShiftB[T]]]
  if reusing:
    stf = work.stf
    stu = work.stu
    ss = newLinkShifts(lo, T, -1)  # started in the threads region
    inc work.reuses
    toc("reuse")
  else:
    (stf,stu,ss) = makeStaples(uu, cs, work.stf, work.stu)
    inc work.builds
    work.haveBackward = true
    work.stag = tag
    work.stf = stf
    work.stu = stu
    toc("makeStaples")
  threads:
    if reusing:
      for mu in 1..<nd:
        for nu in 0..<mu:
          ss[mu][nu].startSB(stu[mu,nu][ix])
          ss[nu][mu].startSB(stu[nu,mu][ix])
          boundaryWaitSB(cs[mu][nu]): discard
          boundaryWaitSB(cs[nu][mu]): discard
      boundarySyncSB()
    if not accumulate:
      for mu in 0..<nd: f[mu] := 0
  addStaplePlanes(uu, f, cp, cr, stf, stu, ru, cs, ss, sf, sb, staplePlanes(nd))
  if c.pgm != 0:
    gaugeLoopDeriv(GaugeActionCoeffs(pgm: c.pgm), uu, uu, f, false, work)
  toc("end")

proc gaugeForce*[V:static[int],T](c: GaugeActionCoeffs, uu: openArray[Field[V,T]], f: array|seq; work: LoopWork[V,T] = nil, tag = default(LoopTag)) =
  tic("gaugeForce")
  gaugeActionDeriv(c, uu, f, work=work, tag=tag)
  toc("gaugeActionDeriv")
  contractProjectTAH(uu, f)
  toc("gaugeForce end")

proc gaugeForce*[T](uu: openArray[T]): auto =
  let lo = uu[0].l
  var f = newOneOf @uu
  let gc = GaugeActionCoeffs(plaq:1.0)
  gc.gaugeForce(uu,f)
  return f

proc gaugeForce*(f,g: array|seq) =
  var c = GaugeActionCoeffs(plaq:1.0)
  gaugeForce(c,g,f)

proc gaugeAction2*(c: GaugeActionCoeffs, g: array|seq): auto =
  requireFundamental(c, "gaugeAction2")
  mixin redot
  tic("gaugeAction2")
  const nc = g[0][0].nrows
  let lo = g[0].l
  let nd = lo.nDim
  let t = newTransporters(g, g[0], 1)
  let t2 = newTransporters(g, g[0], 1)
  let td = newTransporters(g, g[0], -1)
  var pl = 0.0
  var rt = 0.0
  var pg = 0.0
  toc("gaugeAction2 setup")
  threads:
    tic()
    toc("gaugeAction2 zero")
    #var ip = 0
    for mu in 1..<nd:
      for nu in 0..<mu:
        tic()
        var tpl = redot(t[mu]^*g[nu], t[nu]^*g[mu])
        if threadNum==0:
          pl += tpl
        #echo mu, " ", nu, " ", trace(m)/nc
        toc("gaugeAction2 pl")
        if c.rect != 0:
          var tr1 = redot(t[mu]^*t[nu]^*g[nu], t2[nu]^*t[nu]^*g[mu])
          var tr2 = redot(t2[mu]^*t[mu]^*g[nu], t[nu]^*t[mu]^*g[mu])
          if threadNum==0:
            rt += tr1 + tr2
          toc("gaugeAction2 rt")
        if c.pgm != 0:
          for sg in 0..<nu:
            var ts1 = redot(t[mu]^*t[nu]^*g[sg], t[sg]^*t[nu]^*g[mu])
            var ts2 = redot(t[mu]^*t[sg]^*g[nu], t[nu]^*t[sg]^*g[mu])
            var ts3 = redot(t[nu]^*t[mu]^*g[sg], t[sg]^*t[mu]^*g[nu])
            #var ts4 = redot(t[nu]^*t[sg]^*g[mu], t[mu]^*t[sg]^*g[nu])
            #var ts5 = redot(t[sg]^*t[mu]^*g[nu], t[nu]^*t[mu]^*g[sg])
            #var ts6 = redot(t[sg]^*t[nu]^*g[mu], t[mu]^*t[nu]^*g[sg])
            #var ts7 = redot(td[sg]^*t[mu]^*td[nu]^*g[sg], td[nu]^*g[mu])
            #var ts8 = redot(td[sg]^*t[nu]^*td[mu]^*g[sg], td[mu]^*g[nu])
            var ts7 = redot(t[mu]^*td[nu]^*g[sg], t[sg]^*td[nu]^*g[mu])
            #var ts8 = redot(t[mu]^*td[sg]^*g[nu], t[nu]^*td[sg]^*g[mu])
            if threadNum==0:
              #pg += ts1 + ts2 + ts3 + ts4 + ts5 + ts6 + ts7 + ts8
              pg += ts1 + ts2 + ts3 + ts7
          toc("gaugeAction2 pg")
    toc("gaugeAction2 work")
  toc("gaugeAction2 threads")
  #echo "plaq: ", pl, "  rect: ", rt, "  pgm: ", pg
  #result = (pl,rt,pg)
  result = (-1.0/nc.float) * (c.plaq*pl + c.rect*rt + c.pgm*pg)
template gaugeAction2*(g: array|seq, c: GaugeActionCoeffs): untyped =
  gaugeAction2(c, g)
proc gaugeAction2*(g: array|seq): auto =
  var c = GaugeActionCoeffs(plaq:1.0)
  gaugeAction2(c, g)

proc gaugeDeriv2*[V:static[int],T](c: GaugeActionCoeffs, g: openArray[Field[V,T]], f: array|seq; work: LoopWork[V,T] = nil) =
  requireFundamental(c, "gaugeDeriv2")
  mixin adj
  tic("gaugeDeriv2")
  let lo = g[0].l
  let nd = lo.nDim
  const nc = g[0][0].nrows
  let cp = - c.plaq / float(nc)
  let cr = - c.rect / float(nc)
  let t = newTransporters(g, g[0], 1)
  let t2 = newTransporters(g, g[0], 1)
  let tg = newTransporters(g, g[0], 1)
  let td = newTransporters(g, g[0], -1)
  let td2 = newTransporters(g, g[0], -1)
  toc("gaugeForce2 setup")
  let u = cast[ptr cArray[Field[V,T]]](unsafeAddr(g[0]))
  threads:
    for mu in 0..<nd:
      #let mu = (mux + 1) mod nd
      #f[mu] := 0
      for nu in 0..<nd:
        if nu==mu: continue
        discard t[nu] ^* u[mu]
        shiftExpr(t[mu].sb, f[mu][ir] += cp * t[nu].field[ir]*adj(it), u[nu][ix])
        f[mu] += cp * td[nu] ^* t[mu] ^* u[nu]
        if cr != 0:
          discard t2[nu] ^* t[nu] ^* u[mu]
          discard tg[nu] ^* u[nu]
          shiftExpr(t[mu].sb, f[mu][ir] += cr * t2[nu].field[ir]*adj(it), tg[nu].field[ix])
          f[mu] += cr * td2[nu] ^* td[nu] ^* tg[mu] ^* t[nu] ^* u[nu]
          f[mu] += cr * td2[mu] ^* td[nu] ^* tg[mu] ^* t[mu] ^* u[nu]
          discard td[nu] ^* tg[mu] ^* t[mu] ^* u[nu]
          shiftExpr(t2[mu].sb, f[mu][ir] += cr * td[nu].field[ir]*adj(it), u[mu][ix])
          discard td[mu] ^* t[nu] ^* tg[mu] ^* u[mu]
          shiftExpr(t2[mu].sb, f[mu][ir] += cr * td[mu].field[ir]*adj(it), u[nu][ix])
          # Back-to-back shifts through one buffer need every thread done with
          # the previous receive before the next one is posted (see gaugeForce2).
          threadBarrier()
          shiftExpr(t2[mu].sb, f[mu][ir] += cr * t[nu].field[ir]*adj(it), t[mu].field[ix])
  if c.pgm != 0:
    gaugeLoopDeriv(GaugeActionCoeffs(pgm: -c.pgm), g, g, f, false, work)
  toc("end")

proc gaugeDerivSubset[V:static[int],T](c: GaugeActionCoeffs, g: openArray[Field[V,T]], fd: Field[V,T], parity, dir: int, clear: static bool; work: LoopWork[V,T] = nil) =
  # Full derivative into workspace scratch, then select the requested links.
  let
    w = useLoopWork(g, work)
    d = w.loopScratch
    sub = g[0].l.getSubset(if parity == 0: "even" else: "odd")
  gaugeActionDeriv(-1.0*c, g, d, work=w)
  threads:
    for x in fd:
      if x >= sub.lowOuter and x < sub.highOuter:
        fd[x] := d[dir][x]
      else:
        when clear: fd[x] := 0

type
  SubsetShifts*[V:static[int],T] = ref object
    ## Shift buffers that each masked kernel call completes, shared by the
    ## works of all (parity, dir) and made on first use: dn[nd*p+nu] moves by
    ## -nu onto parity p, up[nd*p+nu] by +nu, bk[nd*dir+nu] by -dir.
    lo: Layout[V]
    dn, up, bk: seq[ShiftB[T]]
  SubsetWork*[V:static[int],T] = ref object
    ## Shift buffers of the masked derivative and Hessian at one (parity, dir),
    ## entries nu != dir. Derivative: sd[nu] moves u[nu] by +dir, sf[nu] moves
    ## u[dir] by +nu onto the parity P, sb[nu] moves the lower staples by -nu
    ## onto P. Hessian, made on its first call: so[nu] moves the seed by +nu off
    ## P, sq[nu] moves its products by -dir, sr[nu] by -nu off P. sd and sf keep
    ## the boundaries of their last exchange, which came from the links of tag;
    ## a masked Hessian at the same tag reads them back. sb, so, sq and sr
    ## share their buffers with sh.
    lo: Layout[V]
    parity, dir*: int  # fixed at construction
    sd, sf, sb, so, sq, sr: seq[ShiftB[T]]
    sh: SubsetShifts[V,T]
    tag*: LoopTag

proc newSubsetShifts*[V:static[int],T](f: Field[V,T]): SubsetShifts[V,T] =
  let nd = f.l.nDim
  SubsetShifts[V,T](lo: f.l, dn: newSeq[ShiftB[T]](2*nd), up: newSeq[ShiftB[T]](2*nd), bk: newSeq[ShiftB[T]](nd*nd))

proc entry[T](s: var seq[ShiftB[T]]; i: int; f: Field; dir, len: int; sub: string): ShiftB[T] =
  ## The ShiftB at s[i], made on first use; a copy shares its buffers.
  if s[i].si.isNil: s[i].initShiftB(f, dir, len, sub)
  s[i]

proc newSubsetWork*[V:static[int],T](f: Field[V,T], parity, dir: int; sh: SubsetShifts[V,T] = nil): SubsetWork[V,T] =
  ## Geometry is fixed; buffers are made on first use. Shared calls are serialized.
  ## parity is 0 (even) or 1 (odd), and 0 <= dir < nDim.
  if parity notin 0..1 or dir notin 0..<f.l.nDim:
    raise newException(ValueError, "subset work requires parity 0 or 1 and a lattice direction")
  if sh != nil and sh.lo != f.l:
    raise newException(ValueError, "subset shifts require their original layout")
  SubsetWork[V,T](lo: f.l, parity: parity, dir: dir,
    sh: if sh.isNil: newSubsetShifts(f) else: sh)

proc useSubsetWork(f: Field, sw: SubsetWork, parity, dir: int): auto =
  result = if sw.isNil: newSubsetWork(f, parity, dir) else: sw
  if result.lo != f.l or result.parity != parity or result.dir != dir:
    raise newException(ValueError, "subset workspace requires its original layout, parity and direction")

proc initShifts[V:static[int],T](sw: SubsetWork[V,T], f: Field[V,T], hessian: static bool) =
  let nd = sw.lo.nDim
  if sw.sd.len == 0:
    let ps = if sw.parity == 0: "even" else: "odd"
    sw.sd.newSeq(nd)
    sw.sf.newSeq(nd)
    sw.sb.newSeq(nd)
    for nu in 0..<nd:
      if nu != sw.dir:
        sw.sd[nu].initShiftB(f, sw.dir, 1)
        sw.sf[nu].initShiftB(f, nu, 1, ps)
        sw.sb[nu] = sw.sh.dn.entry(nd*sw.parity+nu, f, nu, -1, ps)
  when hessian:
    if sw.so.len != 0: return
    let po = 1 - sw.parity
    let ps = if po == 0: "even" else: "odd"
    sw.so.newSeq(nd)
    sw.sq.newSeq(nd)
    sw.sr.newSeq(nd)
    for nu in 0..<nd:
      if nu != sw.dir:
        sw.so[nu] = sw.sh.up.entry(nd*po+nu, f, nu, 1, ps)
        sw.sq[nu] = sw.sh.bk.entry(nd*sw.dir+nu, f, sw.dir, -1, "all")
        sw.sr[nu] = sw.sh.dn.entry(nd*po+nu, f, nu, -1, ps)

proc gaugeDeriv2SubsetWork*[V:static[int],T](c: GaugeActionCoeffs, g: openArray[Field[V,T]], fd: Field[V,T], sw: SubsetWork[V,T], clear: static bool; tag = default(LoopTag); work: LoopWork[V,T] = nil) =
  ## fd|P = D(g)[dir]|P for the direction field fd, shifting through sw.
  ## tag: identity and version of g, which sw keeps for a Hessian at g.
  ## clear=false preserves fd|~P; clear=true sets it to zero.
  ## rect/pgm coefficients take the full-derivative route and leave sw untagged.
  requireFundamental(c, "gaugeDeriv2SubsetWork")
  if g[0].l != sw.lo or fd.l != sw.lo:
    raise newException(ValueError, "subset fields require the workspace layout")
  let parity = sw.parity
  let dir = sw.dir
  sw.tag = default(LoopTag)
  if c.rect != 0 or c.pgm != 0:
    gaugeDerivSubset(c, g, fd, parity, dir, clear, work)
    return
  mixin adj, load1
  tic("gaugeDeriv2Subset")
  sw.initShifts(g[0], false)
  let lo = g[0].l
  let nd = lo.nDim
  const nc = g[0][0].nrows
  let cp = -c.plaq / float(nc)
  let sub = lo.getSubset(if parity == 0: "even" else: "odd")
  when clear:
    let other = lo.getSubset(if parity == 0: "odd" else: "even")
  let u = cast[ptr cArray[Field[V,T]]](unsafeAddr(g[0]))
  type M = type(load1(u[0][0]))
  # One pass over P after the link exchanges. With a_nu(y) = u_nu(y+dir) and
  # the lower staple l_nu(y) = u_nu(y)^dag u_dir(y) a_nu(y),
  #   fd(x) = cp sum_nu [u_nu(x) u_dir(x+nu) a_nu(x)^dag + l_nu(x-nu)],  x in P.
  template lower(nu, y0: untyped): untyped =
    block:
      let y = y0  # bound before getSB injects its own ix
      var a {.noinit.}: M
      getSB(sw.sd[nu], y, assign(a, it), u[nu][ix])
      u[nu][y].adj * (u[dir][y] * a)
  threads:
    when clear:
      for x in other:
        fd[x] := 0
    for nu in 0..<nd:
      if nu != dir:
        startSB(sw.sd[nu], u[nu][ix])
        startSB(sw.sf[nu], u[dir][ix])
    for nu in 0..<nd:
      if nu != dir:
        waitRecvSB(sw.sd[nu])
        waitRecvSB(sw.sf[nu])
    boundarySyncSB()
    # The lower staples at the send sites read the a_nu boundary.
    for nu in 0..<nd:
      if nu != dir:
        startSB(sw.sb[nu], lower(nu, ix))
    toc("gaugeDeriv2Subset exchange")
    for x in sub:
      var s {.noinit.}: M
      s := 0
      for nu in 0..<nd:
        if nu != dir:
          var a, b {.noinit.}: M
          getSB(sw.sd[nu], x, assign(a, it), u[nu][ix])
          getSB(sw.sf[nu], x, assign(b, it), u[dir][ix])
          s += (u[nu][x] * b) * a.adj
          localSB(sw.sb[nu], x, s += it, lower(nu, ix))
      fd[x] := cp * s
    toc("gaugeDeriv2Subset local")
    # The boundary loops split P like the local loop.
    for nu in 0..<nd:
      if nu != dir:
        boundarySB(sw.sb[nu], fd[ir] += cp*it)
        doneSB(sw.sd[nu])
        doneSB(sw.sf[nu])
  sw.tag = tag
  toc("gaugeDeriv2Subset end")

proc gaugeDeriv2Subset*[V:static[int],T](c: GaugeActionCoeffs, g: openArray[Field[V,T]], f: array|seq, parity, dir: int; work: LoopWork[V,T] = nil) =
  ## f[dir]|P = D(g)|P, f[dir]|~P = 0; other directions are unchanged.
  if c.rect != 0 or c.pgm != 0:
    gaugeDerivSubset(c, g, f[dir], parity, dir, true, work)
    return
  c.gaugeDeriv2SubsetWork(g, f[dir], newSubsetWork(g[0], parity, dir), true, work=work)

proc gaugeDerivDeriv2*[V:static[int],T](c: GaugeActionCoeffs, g: openArray[Field[V,T]], h, f: array|seq; work: LoopWork[V,T] = nil) =
  ## f += H_g(h), under the ambient real Frobenius pairing.
  requireFundamental(c, "gaugeDerivDeriv2")
  mixin adj
  tic("gaugeDeriv2")
  let lo = g[0].l
  let nd = lo.nDim
  const nc = g[0][0].nrows
  let cp = - c.plaq / float(nc)
  var t, td: typeof(newTransporters(g, g[0], 1))
  var th, thd: typeof(newTransporters(h, g[0], 1))
  defer:
    if work != nil: work.clearInputs
  when typeof(g[0]) is typeof(h[0]):
    if work != nil:
      (t, td, th, thd) = work.hessianTrans(g, h)
  if t.len == 0:
    t = newTransporters(g, g[0], 1)
    td = newTransporters(g, g[0], -1)
    th = newTransporters(h, g[0], 1)
    thd = newTransporters(h, g[0], -1)
  toc("gaugeForce2 setup")
  let u = cast[ptr cArray[Field[V,T]]](unsafeAddr(g[0]))
  threads:
    for mu in 0..<nd:
      for nu in 0..<nd:
        if nu==mu: continue
        discard t[nu] ^* u[mu]
        shiftExpr(t[mu].sb, f[mu][ir] += cp * t[nu].field[ir]*adj(it), h[nu][ix])
        discard t[nu] ^* h[mu]
        shiftExpr(t[mu].sb, f[mu][ir] += cp * t[nu].field[ir]*adj(it), u[nu][ix])
        discard th[nu] ^* u[mu]
        shiftExpr(t[mu].sb, f[mu][ir] += cp * th[nu].field[ir]*adj(it), u[nu][ix])
        f[mu] += cp * td[nu] ^* t[mu] ^* h[nu]
        f[mu] += cp * td[nu] ^* th[mu] ^* u[nu]
        f[mu] += cp * thd[nu] ^* t[mu] ^* u[nu]
  if c.rect != 0 or c.pgm != 0:
    gaugeLoopDeriv(GaugeActionCoeffs(rect: -c.rect, pgm: -c.pgm), g, h, f, true, work)
  toc("end")

type SubsetHessianMode = enum hmOverwrite, hmAccumulate, hmAddBase, hmSum

proc gaugeDerivDeriv2SubsetImpl[V:static[int],T](c: GaugeActionCoeffs, g: openArray[Field[V,T]], hs, hdir: auto, f: array|seq; mode: static SubsetHessianMode; sw: SubsetWork[V,T]; tag: LoopTag; work: LoopWork[V,T] = nil) =
  # S=(parity,dir). hs supplies the summands for hmSum or the base for hmAddBase.
  # hmSum forms hdir|P=sum(hs)|P, then overwrites f with its Hessian.
  # hmAddBase: f[S]+=H(hdir|P); f[~S]=base[~S]+H(hdir|P).
  # The link shifts go through sw; at sw's tag they reuse its boundaries.
  requireFundamental(c, "gaugeDerivDeriv2Subset")
  let parity = sw.parity
  let dir = sw.dir
  if c.rect != 0 or c.pgm != 0:
    sw.tag = default(LoopTag)
    # Apply the full Hessian to a masked seed kept in workspace scratch.
    let
      w = useLoopWork(g, work)
      hm = w.loopScratch
      sub = g[0].l.getSubset(if parity == 0: "even" else: "odd")
    when mode == hmAddBase:
      let other = g[0].l.getSubset(if parity == 0: "odd" else: "even")
    threads:
      for mu in 0..<hm.len:
        hm[mu] := 0
      threadBarrier()
      for x in sub:
        when mode == hmSum:
          hdir[x] := hs[0][x]
          for k in 1..<hs.len:
            hdir[x] += hs[k][x]
        hm[dir][x] := hdir[x]
      threadBarrier()
      when mode == hmAddBase:
        for mu in 0..<f.len:
          if mu != dir:
            f[mu] := hs[mu]
        for x in other:
          f[dir][x] := hs[dir][x]
      elif mode in {hmOverwrite, hmSum}:
        for mu in 0..<f.len:
          f[mu] := 0
    c.gaugeDerivDeriv2(g, hm, f, work=w)
    return
  mixin adj, load1
  tic("gaugeDerivDeriv2Subset")
  let lo = g[0].l
  let nd = lo.nDim
  const nc = g[0][0].nrows
  let cp = -c.plaq / float(nc)
  let sub = lo.getSubset(if parity == 0: "even" else: "odd")
  let other = lo.getSubset(if parity == 0: "odd" else: "even")
  sw.initShifts(g[0], true)
  let kept = tag.id != 0 and tag == sw.tag
  sw.tag = default(LoopTag)
  let u = cast[ptr cArray[Field[V,T]]](unsafeAddr(g[0]))
  type M = type(load1(u[0][0]))
  # Each (dir,nu) plaquette holds one dir link on P, so with H = h on P and
  # u_dir off P the Hessian is the (dir,nu) plaquette derivative at u_dir -> H.
  # One pass over all sites after the seed exchange; with a_nu(y) = u_nu(y+dir),
  #   f_nu(y) += cp [H(y) a_nu(y) H(y+nu)^dag + q_nu(y-dir)],
  #   f_dir(y) += cp sum_nu [u_nu(y) H(y+nu) a_nu(y)^dag + r_nu(y-nu)],  y off P,
  #   q_nu(x) = H(x)^dag u_nu(x) H(x+nu),  r_nu(x) = u_nu(x)^dag h(x) a_nu(x).
  # Off-rank q and r come from their send sites; the rest is read in place.
  template inP(x: untyped): bool = x >= sub.lowOuter and x < sub.highOuter
  template anbr(nu, x0: untyped): untyped =
    block:
      let x = x0  # bound before getSB injects its own ix
      var t {.noinit.}: M
      getSB(sw.sd[nu], x, assign(t, it), u[nu][ix])
      t
  template hnbr(nu, x0: untyped): untyped =
    block:
      let x = x0
      var t {.noinit.}: M
      if inP(x): getSB(sw.sf[nu], x, assign(t, it), u[dir][ix])
      else: getSB(sw.so[nu], x, assign(t, it), hdir[ix])
      t
  template q(nu, x0: untyped): untyped =
    block:
      let x = x0
      let n = hnbr(nu, x)
      var t {.noinit.}: M
      if inP(x): t := hdir[x].adj * (u[nu][x] * n)
      else: t := u[dir][x].adj * (u[nu][x] * n)
      t
  template r(nu, x0: untyped): untyped =
    block:
      let x = x0
      u[nu][x].adj * (hdir[x] * anbr(nu, x))
  template put(mu, y, t: untyped) =
    when mode == hmAddBase: f[mu][y] := hs[mu][y] + cp*t
    elif mode == hmAccumulate: f[mu][y] += cp*t
    else: f[mu][y] := cp*t
  threads:
    when mode == hmSum:
      for x in sub:
        hdir[x] := hs[0][x]
        for k in 1..<hs.len:
          hdir[x] += hs[k][x]
      threadBarrier()
    for nu in 0..<nd:
      if nu != dir:
        if not kept:
          startSB(sw.sd[nu], u[nu][ix])
          startSB(sw.sf[nu], u[dir][ix])
        startSB(sw.so[nu], hdir[ix])
    for nu in 0..<nd:
      if nu != dir:
        if not kept:
          waitRecvSB(sw.sd[nu])
          waitRecvSB(sw.sf[nu])
        waitRecvSB(sw.so[nu])
    boundarySyncSB()
    for nu in 0..<nd:
      if nu != dir:
        startSB(sw.sq[nu], q(nu, ix))
        startSB(sw.sr[nu], r(nu, ix))
    toc("gaugeDerivDeriv2Subset exchange")
    # The sites of both parities differ in work: each iteration takes one of
    # each at the same position in its half, adjacent sites sharing neighbors.
    for y in sub:
      for nu in 0..<nd:
        if nu != dir:
          var a, n, t {.noinit.}: M
          getSB(sw.sd[nu], y, assign(a, it), u[nu][ix])
          getSB(sw.sf[nu], y, assign(n, it), u[dir][ix])
          t := (hdir[y] * a) * n.adj
          localSB(sw.sq[nu], y, t += it, q(nu, ix))
          put(nu, y, t)
      when mode in {hmOverwrite, hmSum}: f[dir][y] := 0
      let z = y - sub.lowOuter + other.lowOuter
      var e {.noinit.}: M
      e := 0
      for nu in 0..<nd:
        if nu != dir:
          var a, n, k, t {.noinit.}: M
          getSB(sw.sd[nu], z, assign(a, it), u[nu][ix])
          getSB(sw.so[nu], z, assign(n, it), hdir[ix])
          k := a * n.adj
          t := u[dir][z] * k
          e += u[nu][z] * k.adj
          localSB(sw.sr[nu], z, e += it, r(nu, ix))
          localSB(sw.sq[nu], z, t += it, q(nu, ix))
          put(nu, z, t)
      put(dir, z, e)
    toc("gaugeDerivDeriv2Subset local")
    # The boundary loops of sq split all sites unlike the parity loops.
    threadBarrier()
    for nu in 0..<nd:
      if nu != dir:
        boundarySB(sw.sq[nu], f[nu][ir] += cp*it)
        boundarySB(sw.sr[nu], f[dir][ir] += cp*it)
        doneSB(sw.so[nu])
        if not kept:
          doneSB(sw.sd[nu])
          doneSB(sw.sf[nu])
  sw.tag = tag
  toc("gaugeDerivDeriv2Subset end")

# The masked Hessians below take the derivative's shift buffers sw (a new set
# when nil) and the tag of g: at sw's tag they reuse its link boundaries.

proc gaugeDerivDeriv2Subset*[V:static[int],T](c: GaugeActionCoeffs, g: openArray[Field[V,T]], h, f: array|seq, parity, dir: int; sw: SubsetWork[V,T] = nil; tag = default(LoopTag); work: LoopWork[V,T] = nil) =
  ## f = H_g(h[dir]|P), including all affected links.
  c.gaugeDerivDeriv2SubsetImpl(g, h, h[dir], f, hmOverwrite, useSubsetWork(g[0], sw, parity, dir), tag, work=work)

proc gaugeDerivDeriv2Subset*[V:static[int],T](c: GaugeActionCoeffs, g: openArray[Field[V,T]], hdir: Field[V,T], f: array|seq, parity, dir: int; sw: SubsetWork[V,T] = nil; tag = default(LoopTag); work: LoopWork[V,T] = nil) =
  ## f = H_g(hdir|P in direction dir), including all affected links.
  c.gaugeDerivDeriv2SubsetImpl(g, hdir, hdir, f, hmOverwrite, useSubsetWork(g[0], sw, parity, dir), tag, work=work)

proc gaugeDerivDeriv2SubsetAdd*[V:static[int],T](c: GaugeActionCoeffs, g: openArray[Field[V,T]], hdir: auto, f: array|seq, parity, dir: int; sw: SubsetWork[V,T] = nil; tag = default(LoopTag); work: LoopWork[V,T] = nil) =
  ## f += H_g(hdir|P).
  c.gaugeDerivDeriv2SubsetImpl(g, hdir, hdir, f, hmAccumulate, useSubsetWork(g[0], sw, parity, dir), tag, work=work)

proc gaugeDerivDeriv2SubsetAddBase*[V:static[int],T](c: GaugeActionCoeffs, g: openArray[Field[V,T]], hdir: auto, base, f: array|seq, parity, dir: int; sw: SubsetWork[V,T] = nil; tag = default(LoopTag); work: LoopWork[V,T] = nil) =
  ## S=(P,dir): f[S] += H_g(hdir|P)[S].
  ## f[~S] = base[~S] + H_g(hdir|P)[~S].
  c.gaugeDerivDeriv2SubsetImpl(g, base, hdir, f, hmAddBase, useSubsetWork(g[0], sw, parity, dir), tag, work=work)

proc gaugeDerivDeriv2SubsetSum*[V:static[int],T](c: GaugeActionCoeffs, g: openArray[Field[V,T]], h: array|seq, w: auto, f: array|seq, parity, dir: int; sw: SubsetWork[V,T] = nil; tag = default(LoopTag); work: LoopWork[V,T] = nil) =
  ## w|P = sum(h)|P; f = H_g(w|P).
  c.gaugeDerivDeriv2SubsetImpl(g, h, w, f, hmSum, useSubsetWork(g[0], sw, parity, dir), tag, work=work)

proc gaugeForce2*[V:static[int],T](c: GaugeActionCoeffs, g: openArray[Field[V,T]], f: array|seq; work: LoopWork[V,T] = nil) =
  requireFundamental(c, "gaugeForce2")
  if c.pgm != 0:
    c.gaugeForce(g, f, work=work)
    return
  mixin adj,projectTAH
  tic("gaugeForce2")
  let lo = g[0].l
  let nd = lo.nDim
  const nc = g[0][0].nrows
  let cp = c.plaq / float(nc)
  let cr = c.rect / float(nc)
  let t = newTransporters(g, g[0], 1)
  let t2 = newTransporters(g, g[0], 1)
  let tg = newTransporters(g, g[0], 1)
  let td = newTransporters(g, g[0], -1)
  let td2 = newTransporters(g, g[0], -1)
  toc("gaugeForce2 setup")
  let u = cast[ptr cArray[Field[V,T]]](unsafeAddr(g[0]))
  threads:
    for mu in 0..<nd:
      #let mu = (mux + 1) mod nd
      f[mu] := 0
      for nu in 0..<nd:
        if nu==mu: continue
        discard t[nu] ^* u[mu]
        shiftExpr(t[mu].sb, f[mu][ir] += cp * t[nu].field[ir]*adj(it), u[nu][ix])
        f[mu] += cp * td[nu] ^* t[mu] ^* u[nu]
        if cr != 0:
          discard t2[nu] ^* t[nu] ^* u[mu]
          discard tg[nu] ^* u[nu]
          shiftExpr(t[mu].sb, f[mu][ir] += cr * t2[nu].field[ir]*adj(it), tg[nu].field[ix])
          f[mu] += cr * td2[nu] ^* td[nu] ^* tg[mu] ^* t[nu] ^* u[nu]
          f[mu] += cr * td2[mu] ^* td[nu] ^* tg[mu] ^* t[mu] ^* u[nu]
          discard td[nu] ^* tg[mu] ^* t[mu] ^* u[nu]
          shiftExpr(t2[mu].sb, f[mu][ir] += cr * td[nu].field[ir]*adj(it), u[mu][ix])
          discard td[mu] ^* t[nu] ^* tg[mu] ^* u[mu]
          shiftExpr(t2[mu].sb, f[mu][ir] += cr * td[mu].field[ir]*adj(it), u[nu][ix])
          threadBarrier()
          shiftExpr(t2[mu].sb, f[mu][ir] += cr * t[nu].field[ir]*adj(it), t[mu].field[ix])
    for mu in 0..<nd: contractProjectTAH(f[mu], u[mu], f[mu])
  toc("end")
proc gaugeForce2*(f,g: array|seq) =
  var c = GaugeActionCoeffs(plaq:1.0)
  gaugeForce2(c,g,f)

proc gaugeAction3*[V:static[int],T](c: GaugeActionCoeffs, g: openArray[Field[V,T]]; work: LoopWork[V,T] = nil): auto =
  requireFundamental(c, "gaugeAction3")
  tic("gaugeAction3")
  const nc = g[0][0].nrows
  let lo = g[0].l
  let gs = @g  # wilsonLines captures its links in closures
  var pl = 0.0
  var rt = 0.0
  for mu in 1..<lo.nDim:
    for nu in 0..<mu:
      var ls = newseq[seq[int]]()
      block:
        let
          mu = mu+1
          nu = nu+1
        if c.plaq!=0:
          ls.add @[mu, nu, -mu, -nu]
        if c.rect!=0:
          ls.add [@[mu, nu, nu, -mu, -nu, -nu], @[mu, mu, nu, -mu, -mu, -nu]]
      let ws = gs.wilsonLines ls
      var i = 0
      if c.plaq!=0:
        pl += ws[0].re
        i = 1
      if c.rect!=0:
        rt += ws[i].re + ws[i+1].re
  result = (-lo.physVol.float) * (c.plaq*pl + c.rect*rt)
  if c.pgm != 0:
    result += pgmAction(c, g, work)
  toc("end")
proc gaugeAction3*(g: array|seq): auto =
  var c = GaugeActionCoeffs(plaq:1.0)
  gaugeAction3(c, g)

proc plaqRectPath_fun(c:GaugeActionCoeffs, mu,nu:int):auto =
  let
    mu = mu+1
    nu = nu+1
  var ls = newseq[seq[int]]()
  if c.plaq!=0:
    ls.add [ @[nu, mu, -nu], @[-nu, mu, nu]
           , @[mu, nu, -mu], @[-mu, nu, mu]
           ]
  if c.rect!=0:
    ls.add [ @[nu, nu, mu, -nu, -nu], @[-nu, -nu, mu, nu, nu]
           , @[nu, mu, mu, -nu, -mu], @[-mu, nu, mu, mu, -nu]
           , @[-nu, mu, mu, nu, -mu], @[-mu, -nu, mu, mu, nu]
           , @[mu, mu, nu, -mu, -mu], @[-mu, -mu, nu, mu, mu]
           , @[mu, nu, nu, -mu, -nu], @[-nu, mu, nu, nu, -mu]
           , @[-mu, nu, nu, mu, -nu], @[-nu, -mu, nu, nu, mu]
           ]
  ls.optimalPairs

proc plaqRectPath(c:GaugeActionCoeffs, mu,nu:int):auto =
  var j = 0
  if c.plaq!=0:
    inc j
  if c.rect!=0:
    j += 2
  memoize(j,mu,nu):
    c.plaqRectPath_fun(mu,nu)

proc gaugeForce3*[V:static[int],T](c: GaugeActionCoeffs, g: openArray[Field[V,T]], f: auto; work: LoopWork[V,T] = nil) =
  requireFundamental(c, "gaugeForce3")
  if c.pgm != 0:
    c.gaugeForce(g, f, work=work)
    return
  tic("gaugeForce3")
  const nc = g[0][0].nrows
  let nd = g[0].l.nDim
  let cp = c.plaq / nc.float
  let cr = c.rect / nc.float
  let gs = @g  # gaugeProd captures its links in closures
  threads:
    for mu in 0..<nd:
      f[mu] := 0
  for mu in 1..<nd:
    for nu in 0..<mu:
      let ptree = c.plaqRectPath(mu,nu)
      let ws = gs.gaugeProd ptree
      threads:
        for ir in f[mu]:
          var pmu,rmu,pnu,rnu: type(load1(f[0][0]))
          var i = 0
          if c.plaq!=0:
            for j in 0..<2:
              pmu += ws[j][ir]
              pnu += ws[2+j][ir]
            i = 4
            f[mu][ir] += cp * pmu
            f[nu][ir] += cp * pnu
          if c.rect!=0:
            for j in 0..<6:
              rmu += ws[i+j][ir]
              rnu += ws[6+i+j][ir]
            f[mu][ir] += cr * rmu
            f[nu][ir] += cr * rnu
  threads:
    for mu in 0..<nd: contractProjectTAH(f[mu], gs[mu], f[mu])
  toc("end")
proc gaugeForce3*(f,g: array|seq) =
  var c = GaugeActionCoeffs(plaq:1.0)
  gaugeForce3(c,g,f)

proc actionA*(c: GaugeActionCoeffs, g: auto): auto =
  ## Specialized gauge action for plaq + adjplaq
  requireAdjoint(c, "actionA")
  mixin mul, load1, createShiftBufs, re
  tic("actionA")
  let lo = g[0].l
  let nd = lo.nDim
  let nc = g[0][0].ncols
  var sf = newSeq[type(createShiftBufs(g[0],1,"all"))](nd)
  for i in 0..<nd-1:
    sf[i] = createShiftBufs(g[0], 1, "all")
  sf[nd-1].newSeq(nd)
  for i in 0..<nd-1: sf[nd-1][i] = sf[i][i]
  var pl = [0.0, 0.0]
  toc("plaq setup")
  threads:
    tic()
    var plt = [0.0, 0.0]
    var umunu,unumu: type(load1(g[0][0]))
    for mu in 0..<nd:
      for nu in 0..<nd:
        if mu != nu:
          startSB(sf[mu][nu], g[mu][ix])
    toc("plaq start shifts")
    for ir in g[0]:
      for mu in 1..<nd:
        for nu in 0..<mu:
          if isLocal(sf[mu][nu],ir) and isLocal(sf[nu][mu],ir):
            localSB(sf[mu][nu], ir, mul(unumu,g[nu][ir],it), g[mu][ix])
            localSB(sf[nu][mu], ir, mul(umunu,g[mu][ir],it), g[nu][ix])
            let dt = dot(umunu,unumu)
            plt[0] += simdSum(dt.re)
            plt[1] += simdSum(dt.norm2)
    toc("plaq local")
    var needBoundary = false
    for mu in 0..<nd:
      for nu in 0..<nd:
        if mu != nu:
          boundaryWaitSB(sf[mu][nu]): needBoundary = true
    toc("plaq wait")
    if needBoundary:
      boundarySyncSB()
      for ir in g[0]:
        for mu in 1..<nd:
          for nu in 0..<mu:
            if not isLocal(sf[mu][nu],ir) or not isLocal(sf[nu][mu],ir):
              if isLocal(sf[mu][nu], ir):
                localSB(sf[mu][nu], ir, mul(unumu,g[nu][ir],it), g[mu][ix])
              else:
                boundaryGetSB(sf[mu][nu], ir):
                  mul(unumu, g[nu][ir], it)
              if isLocal(sf[nu][mu], ir):
                localSB(sf[nu][mu], ir, mul(umunu,g[mu][ir],it), g[nu][ix])
              else:
                boundaryGetSB(sf[nu][mu], ir):
                  mul(umunu, g[mu][ir], it)
              let dt = dot(umunu,unumu)
              plt[0] += simdSum(dt.re)
              plt[1] += simdSum(dt.norm2)
    toc("plaq boundary")
    threadSum(plt)
    if threadNum == 0:
      pl[0] = plt[0] / float(nc)
      pl[1] = plt[1] / float(nc*nc)
      rankSum(pl)
    toc("plaq sum")
  let a0 = 0.5 * float(nd*(nd-1)*lo.physVol)
  result = c.plaq*(a0-pl[0]) + c.adjplaq*(a0-pl[1])
  toc("plaq end", flops=lo.nSites.float*float(2*8*nc*nc*nc-1))

proc gaugeADeriv*(c: GaugeActionCoeffs, g,f: auto, accumulate=false) =
  requireAdjoint(c, "gaugeADeriv")
  ## Specialized gauge force for plaq + adjplaq
  ## if accumulate, the derivatives will add to f.
  ## if not, f is set to 0 first.
  mixin load1, adj
  tic("gaugeADeriv")
  let lo = g[0].l
  let nd = lo.nDim
  let nc = g[0][0].ncols
  let cp = c.plaq / float(nc)
  let ca = 2.0 * c.adjplaq / float(nc*nc)
  var cs = startCornerShifts(g)
  toc("gaugeADeriv startCornerShifts")
  var (stf,stu,ss) = makeStaples(g, cs)
  toc("gaugeADeriv makeStaples")
  if not accumulate:
    for i in 0..<nd:
      f[i] := 0
  threads:
    tic()
    for ir in g[0]:
      for mu in 1..<nd:
        for nu in 0..<mu:
          let tmn = dot(stf[mu,nu][ir], g[mu][ir])
          f[mu][ir] += (cp+ca*tmn) * stf[mu,nu][ir]
          let tnm = dot(stf[nu,mu][ir], g[nu][ir])
          f[nu][ir] += (cp+ca*tnm) * stf[nu,mu][ir]
          if isLocal(ss[mu][nu],ir):
            var bmu: type(load1(g[0][0]))
            localSB(ss[mu][nu], ir, assign(bmu,it), stu[mu,nu][ix])
            let tmu = dot(bmu, g[mu][ir])
            f[mu][ir] += (cp+ca*tmu) * bmu
          if isLocal(ss[nu][mu],ir):
            var bnu: type(load1(g[0][0]))
            localSB(ss[nu][mu], ir, assign(bnu,it), stu[nu,mu][ix])
            let tnu = dot(bnu, g[nu][ir])
            f[nu][ir] += (cp+ca*tnu) * bnu
    toc("gaugeADeriv local")
    for mu in 1..<nd:
      for nu in 0..<mu:
        var needBoundary = false
        boundaryWaitSB(ss[mu][nu]): needBoundary = true
        boundaryWaitSB(ss[nu][mu]): needBoundary = true
        if needBoundary:
          boundarySyncSB()
          for ir in lo:
            if not isLocal(ss[mu][nu],ir):
              var bmu: type(load1(g[0][0]))
              getSB(ss[mu][nu], ir, assign(bmu,it), stu[mu,nu][ix])
              let tmu = dot(bmu, g[mu][ir])
              f[mu][ir] += (cp+ca*tmu) * bmu
            if not isLocal(ss[nu][mu],ir):
              var bnu: type(load1(g[0][0]))
              getSB(ss[nu][mu], ir, assign(bnu,it), stu[nu,mu][ix])
              let tnu = dot(bnu, g[nu][ir])
              f[nu][ir] += (cp+ca*tnu) * bnu
    #toc("gaugeADeriv boundary")
  toc("gaugeADeriv threads")

proc forceA*(c: GaugeActionCoeffs, g,f: auto) =
  tic("forceA")
  gaugeADeriv(c, g, f)
  toc("gaugeADeriv")
  contractProjectTAH(g, f)
  #contractProjectTAH(f, g)
  toc("forceA end")

proc action*[V:static[int],T](c: GaugeActionCoeffs, g: openArray[Field[V,T]]; work: LoopWork[V,T] = nil, tag = default(LoopTag)): float =
  ## Dispatch by coefficient family; mixed adjoint/improved terms are rejected.
  ## Fundamental terms run the loop kernels: one link exchange per call and no
  ## staple fields; gaugeAction1 and gaugeForce remain the staple path.
  if c.adjplaq != 0: c.actionA(@g)  # actionA captures its links in threads
  else: c.loopAction(g, work=work, tag=tag)

proc force*[V:static[int],T](c: GaugeActionCoeffs, g: openArray[Field[V,T]], f: auto; work: LoopWork[V,T] = nil, tag = default(LoopTag)) =
  ## Projected force for action(c,g). Keep f disjoint from g.
  if c.adjplaq != 0: c.forceA(@g, f)  # forceA captures its links in threads
  else:
    (-1.0*c).loopDeriv(g, f, work=work, tag=tag)  # f = -grad S, as gaugeActionDeriv
    contractProjectTAH(g, f)

proc forceKick*[V:static[int],T](c: GaugeActionCoeffs, g, p, q: openArray[Field[V,T]]; a: float; work: LoopWork[V,T] = nil, tag = default(LoopTag)): Norm2Stats =
  ## q = a F + p with F the projected force of `force`, in loopDeriv's site pass:
  ##   d = -grad S at the site, F = projectTAH(U d†),
  ## so F is never stored. Returns the Norm2Stats of F. Fundamental coefficients;
  ## q may alias p, not g. g, p and q may be temporaries. tag: as in loopAction.
  requireFundamental(c, "forceKick")
  let lo = g[0].l
  let nd = lo.nDim
  if g.len != nd or p.len != nd or q.len != nd or p[0].l != lo or q[0].l != lo:
    raise newException(ValueError, "force kick requires matching gauge bundles")
  for x in q:
    for y in g:
      if x.s.data == y.s.data:
        raise newException(ValueError, "force kick output must not alias the links")
  # Copies keep the arrays referenced while the kernel allocates and through
  # the threads regions; refc does not see raw pointers into them.
  let
    u = @g
    pp = @p
    qq = @q
    mask = loopMask(c,nd)
  if mask == 0:
    threads:
      for mu in 0..<nd: qq[mu] := pp[mu]
    return
  let w = useLoopWork(g,work)
  defer: w.clearInputs
  let plan = gaugeLoops(nd,mask,0)
  let products = loopProducts(plan,0)
  w.prepareProducts(products)
  if tag.id != 0 and tag == w.ltag: w.lpair = true
  let h = [w.gaugeLoopHalo(g,plan,0,tag)]
  w.ltag = tag
  let
    nc = float(g[0][0].nrows)
    cs = [c.plaq/nc, c.rect/nc, c.pgm/nc]  # the sums are -grad S
  var r: Norm2Stats
  let hl = h[0][0].layout
  let tile = w.tile
  threads:
    let sc = w.products[threadNum].data
    let ptrs = w.ptrs[threadNum].data
    products.bindSlots(sc, ptrs, tile)
    let d = products.slots*tile  # sum of direction mu at site b: sc[d+mu*tile+b]
    var s: Norm2Acc
    tileSites(qq[0].l.nSitesOuter, tile, i0, nb):
      products.evalProducts(h,hl,i0,nb,tile,ptrs,outp,src,b):
        addOut(sc[d+outp.mu*tile+b], cs[ord(outp.kind)], outp, src)
      for b in 0..<nb:
        let i = i0+b
        for mu in 0..<nd:
          var f {.noinit.}: evalType(u[mu][i])
          let m = u[mu][i] * sc[d+mu*tile+b].adj
          f.projectTAH m
          s.add f.norm2
          qq[mu][i] := a*f + pp[mu][i]
    let t = s.total
    threadSingle: r = t
  r

when isMainModule:
  import qex
  import physics/qcdTypes
  #import matrixFunctions
  qexInit()
  var defaultGaugeFile = "l88.scidac"
  #let defaultLat = @[2,2,2,2]
  let defaultLat = @[8,8,8,8]
  #let defaultLat = @[8,8,8]
  #let defaultLat = @[8,8]
  defaultSetup()
  #for mu in 0..<g.len: g[mu] := 1
  g.random

  proc test(g:auto) =
    tic("test")
    echo "Test C_plaq = 1"
    var pl = plaq(g)
    echo "plaq:"
    echo pl
    echo pl.sum
    var f = newOneOf g
    var f2 = newOneOf g
    var f3 = newOneOf g
    var fa = newOneOf g
    var gc = GaugeActionCoeffs(plaq:1.0)
    #var ga = gaugeAction(g)
    toc("plaq")
    var ga = gaugeAction.gaugeAction1(g)
    toc("ga1")
    var ga2 = gaugeAction.gaugeAction2(gc,g)
    toc("ga2")
    var ga3 = gaugeAction.gaugeAction3(gc,g)
    toc("ga3")
    var gaA = gaugeAction.actionA(gc,g)
    toc("aA")
    gaugeAction.gaugeForce(gc,g,f)
    toc("gf")
    gaugeAction.gaugeForce2(gc,g,f2)
    toc("gf2")
    gaugeAction.gaugeForce3(gc,g,f3)
    toc("gf3")
    gaugeAction.forceA(gc,g,fa)
    toc("fA")
    echo "ga: ", ga, "\t", ga2, "\t", ga3, "\t", gaA
    for i in 0..<f.len:
      echo "f[", i, "]: ", f[i].norm2, "\t", f2[i].norm2, "\t", f3[i].norm2, "\t", fa[i].norm2
    toc("end")

  proc testR(g:auto) =
    tic("testR")
    echo "Test Rectangle"
    var gc = gaugeAction.DBW2(0.7796)
    echo gc
    var ga = gaugeAction.gaugeAction1(gc,g)
    toc("ga1")
    var ga2 = gaugeAction.gaugeAction2(gc,g)
    toc("ga2")
    var ga3 = gaugeAction.gaugeAction3(gc,g)
    toc("ga3")
    echo "ga: ", ga, "\t", ga2, "\t", ga3
    var f = newOneOf g
    var f2 = newOneOf g
    var f3 = newOneOf g
    var fA = newOneOf g
    toc("init f")
    gaugeAction.gaugeForce(gc,g,f)
    toc("gf")
    gaugeAction.gaugeForce2(gc,g,f2)
    toc("gf2")
    gaugeAction.gaugeForce3(gc,g,f3)
    toc("gf3")
    gaugeAction.forceA(gc,g,fA)
    toc("gfA")
    echo "gf: \t",  f[0].norm2, "\t",  f[1].norm2, "\t",  f[2].norm2, "\t",  f[3].norm2
    echo "gf2:\t", f2[0].norm2, "\t", f2[1].norm2, "\t", f2[2].norm2, "\t", f2[3].norm2
    echo "gf3:\t", f3[0].norm2, "\t", f3[1].norm2, "\t", f3[2].norm2, "\t", f3[3].norm2
    toc("end")

  test(g)
  testR(g)

  proc updateX(g,p,eps:auto) =
    axexpahmuly(g, eps, p, g)

  proc updateP(c:GaugeActionCoeffs, g,p,eps:auto) =
    var f = newOneOf g
    gaugeAction.gaugeForce(c, g, f)
    for mu in 0..<f.len:
      #echo "f[", mu, "]: ", f[mu].norm2
      p[mu] += (-eps)*f[mu]

  var g0 = g[0].l.newGauge
  for mu in 0..<g.len: g0[mu] := g[mu]
  proc test2(steps:int, lambda=0.1931833):auto {.discardable.} =
    const t = 0.02
    let eps = t/steps.float
    echo "eps: ",eps
    var p = newSeq[type(g[0])](g.len)
    for mu in 0..<p.len:
      g[mu] := g0[mu]
      p[mu].new(g[0].l)
      for e in p[mu]:
        when p[mu][e].nrows==1:
          #p[mu][e] := asImag(1)
          p[mu][e] := newComplex(0,1)
        else:
          p[mu][e] := 0
          let t = (2*(e mod 2)-1).float
          #let t = 1.0
          p[mu][e][0,1] := t
          p[mu][e][1,0] := -t
    var gc = gaugeAction.DBW2(0.7796)
    let ga = gaugeAction.gaugeAction1(gc,g)
    var p2 = 0.0
    for mu in 0..<p.len: p2 += p[mu].norm2
    let s0 = ga + 0.5*p2
    echo "ACT: ", ga, "\t", 0.5*p2, "\t", s0

    for n in 1..steps:
      #echo "pdiff: ", (p[0]-p[1]).norm2
      #echo "gdiff: ", (g[0]-g[1]).norm2
      #echo "ga: ", gaugeAction.gaugeAction1(gc,g)
      updateX(g,p,lambda*eps)
      #echo "pdiff: ", (p[0]-p[1]).norm2
      #echo "gdiff: ", (g[0]-g[1]).norm2
      #echo "ga: ", gaugeAction.gaugeAction1(gc,g)
      gc.updateP(g,p,0.5*eps)
      #echo "pdiff: ", (p[0]-p[1]).norm2
      #echo "gdiff: ", (g[0]-g[1]).norm2
      #echo "ga: ", gaugeAction.gaugeAction1(gc,g)
      updateX(g,p,(1.0-2.0*lambda)*eps)
      #echo "pdiff: ", (p[0]-p[1]).norm2
      #echo "gdiff: ", (g[0]-g[1]).norm2
      #echo "ga: ", gaugeAction.gaugeAction1(gc,g)
      gc.updateP(g,p,0.5*eps)
      #echo "pdiff: ", (p[0]-p[1]).norm2
      #echo "gdiff: ", (g[0]-g[1]).norm2
      #echo "ga: ", gaugeAction.gaugeAction1(gc,g)
      updateX(g,p,lambda*eps)

    let ga2 = gaugeAction.gaugeAction1(gc,g)
    p2 = 0.0
    for mu in 0..<p.len: p2 += p[mu].norm2
    let s2 = ga2 + 0.5*p2
    echo "ACT2: ", ga2, "\t", 0.5*p2, "\t", s2
    echo "dH: ", s2 - s0
    let sr = (s2-s0)/(eps*eps)
    echo "error rate: ", sr
    return (s2-s0)/s0

  proc testE4(lambda:float, steps=4):auto =
    # e_n = a*t^2/n^2 + b*t^4/n^4 + ...
    let
      e1 = test2(steps, lambda)
      e4 = test2(4*steps, lambda)
    return e1-16.0*e4

  var lambda = 0.23748

  when false:
    # Search for the lambda that cancels higher order terms.
    let tol = 1e-14
    var
      xlo = 0.15
      xhi = 0.25
      elo = testE4(xlo)
      ehi = testE4(xhi)
      x,e:float
    while elo>0 xor ehi>0:
      x = (ehi*xlo-elo*xhi)/(ehi-elo)
      e = testE4(x)
      echo "lambda: ",x," err_4/s: ",e
      if abs(e)<tol:
        break
      if e>0 xor ehi>0:
        xlo = x
        elo = e
      else:
        xhi = x
        ehi = e
    lambda = x

  test2(200, lambda)
  test2(20, lambda)
  test2(2, lambda)

  let dev = testE4(lambda,10)
  echo "Relative deviation from dt^2 scaling: ",dev
  if abs(dev)>1e-13:
    qexError "Large deviation."

  #echoTimers()
  echoProf()
  qexFinalize()
