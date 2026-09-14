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
    offsets: seq[seq[int32]]
    products: array[6,LoopProducts]
    factors: seq[tuple[mu:int, offset:seq[int32]]]
    pairs: PathPlan
    outputs: seq[tuple[mu:int, kind:LoopKind]]
  LoopLoad = object
    dest, mu, seed: int
    offset: seq[int32]
  LoopMul = object
    dest, left, right: int
    la, ra, add: bool
  LoopOut = object
    source, mu: int
    kind: LoopKind
    adj: bool
  LoopProducts = ref object
    loads: seq[LoopLoad]
    steps: seq[LoopMul]
    outs: seq[LoopOut]
    size: int

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
    let plan = GaugeLoopPlan(paths: makeLoops(), fw: newSeq[int32](nd), bw: newSeq[int32](nd))
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
          if off notin plan.offsets: plan.offsets.add off
          if step > 0: inc x[mu]
    cache[key] = plan
  cache[key]

type
  LoopHalos[V:static[int],T] = object
    layout*: HaloLayout[Layout[V]]
    fields*: seq[seq[Halo[Layout[V],Field[V,T],T]]]
  LoopWork*[V:static[int],T] = ref object
    ## Mutable loop buffers for one layout. Serialize calls using this object.
    ## Geometry initialization follows the halo API's threading requirements.
    lo*: Layout[V]
    halos*: seq[LoopHalos[V,T]]
    scratch: seq[Field[V,T]]
    trans: array[4,seq[Transporter[Field[V,T],Field[V,T],T]]]
    plan: LoopProducts
    products: seq[alignedMem[T]]

proc newLoopWork*[V:static[int],T](f: Field[V,T]): LoopWork[V,T] =
  LoopWork[V,T](lo: f.l)

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

proc gaugeLoopHalo[V:static[int],T](w: LoopWork[V,T], g: openArray[Field[V,T]], plan: GaugeLoopPlan, slot: int): auto =
  if g[0].l != w.lo:
    raise newException(ValueError, "loop inputs require the workspace layout")
  let
    lo = g[0].l
    nd = lo.nDim
    comm = getDefaultComm()
    hl = haloLayout(lo, plan.fw, plan.bw)
    hm = haloMap(hl, comm, plan.offsets)
  var row = 0
  while row < w.halos.len and w.halos[row].layout != hl: inc row
  if row == w.halos.len: w.halos.add LoopHalos[V,T](layout: hl)
  if w.halos[row].fields.len <= slot:
    w.halos[row].fields.setLen(max(2,slot+1))
  if w.halos[row].fields[slot].len == 0:
    w.halos[row].fields[slot].newSeq(nd)
    for mu in 0..<nd:
      w.halos[row].fields[slot][mu] = makeHalo(hl, g[mu])
  let gh = w.halos[row].fields[slot]
  for mu in 0..<nd:
    gh[mu].field = g[mu]
    gh[mu].update(hm, comm)
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

proc loopCoeff(c: GaugeActionCoeffs, kind: LoopKind): float =
  case kind
  of lkPlaq: c.plaq
  of lkRect: c.rect
  of lkPgm: c.pgm

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
      p.loads.add LoopLoad(dest:result,mu:f.mu,seed:firstSetBit(mask),offset:f.offset)
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
  plan.products[order] = p
  p

proc prepareProducts(w: LoopWork, plan: LoopProducts) =
  w.plan = plan
  let nt = getMaxThreads()
  if w.products.len < nt: w.products.setLen(nt)
  for tid in 0..<nt:
    if w.products[tid].len < plan.size+w.lo.nDim:
      w.products[tid].newU(plan.size+w.lo.nDim)

proc evalProducts(plan: LoopProducts, h: auto, i: int, values: auto) =
  for op in plan.loads:
    let hl = h[op.seed][op.mu].layout
    var j = i
    for d, n in op.offset:
      if n > 0:
        for _ in 0..<int(n): j = hl.neighborFwd[d][j]
      elif n < 0:
        for _ in 0..<int(-n): j = hl.neighborBck[d][j]
    values[op.dest] := h[op.seed][op.mu][j]
  for op in plan.steps:
    template multiply(a,b: untyped) =
      if op.add: values[op.dest] += a*b
      else: values[op.dest] := a*b
    if op.la:
      if op.ra: multiply(values[op.left].adj,values[op.right].adj)
      else: multiply(values[op.left].adj,values[op.right])
    else:
      if op.ra: multiply(values[op.left],values[op.right].adj)
      else: multiply(values[op.left],values[op.right])

proc loopAction*[G:array|seq](c: GaugeActionCoeffs, g: G; work: typeof(newLoopWork(g[0])) = nil): float =
  ## Fused S = -sum(c_loop Re tr U_loop)/Nc for fundamental loop terms.
  requireFundamental(c, "loopAction")
  let nd = g[0].l.nDim
  if g.len != nd: raise newException(ValueError, "loop action requires a complete gauge bundle")
  let mask = loopMask(c,nd)
  if mask == 0: return
  let w = useLoopWork(g,work)
  defer: w.clearInputs
  let plan = gaugeLoops(nd,mask,1)
  let products = loopProducts(plan,0)
  w.prepareProducts(products)
  let h = [w.gaugeLoopHalo(g,plan,0)]
  var sums = newSeq[array[3,float]](getMaxThreads())
  threads:
    let values = w.products[threadNum]
    var s: array[3,float]
    for i in g[0]:
      products.evalProducts(h,i,values)
      for outp in products.outs:
        if outp.adj:
          s[ord(outp.kind)] += simdSum(redot(g[outp.mu][i],values[outp.source].adj))
        else:
          s[ord(outp.kind)] += simdSum(redot(g[outp.mu][i],values[outp.source]))
    sums[threadNum] = s
  var total: array[3,float]
  for s in sums:
    for k in 0..<3: total[k] += s[k]
  rankSum(total)
  -(c.plaq*total[0]+c.rect*total[1]+c.pgm*total[2])/float(g[0][0].nrows)

proc loopDeriv*[G:array|seq,H; K:static int](c: GaugeActionCoeffs, g:G, ds:array[K,H], f:array|seq; work:typeof(newLoopWork(g[0])) = nil) =
  ## f = D^K grad S(g)[ds]. The seed count is static; orders above degree vanish.
  ## Output must be disjoint from seeds and, below the polynomial degree, g.
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
  if K < degree: h[0] = w.gaugeLoopHalo(g,plan,0)
  for k in 0..<K: h[k+1] = w.gaugeLoopHalo(ds[k],plan,k+1)
  let nc = float(g[0][0].nrows)
  threads:
    let values = w.products[threadNum]
    for i in f[0]:
      products.evalProducts(h,i,values)
      for mu in 0..<nd: values[products.size+mu] := 0
      for outp in products.outs:
        let a = -loopCoeff(c,outp.kind)/nc
        if outp.adj: values[products.size+outp.mu] += a*values[outp.source].adj
        else: values[products.size+outp.mu] += a*values[outp.source]
      for mu in 0..<nd: f[mu][i] := values[products.size+mu]

proc loopDeriv*[G:array|seq](c: GaugeActionCoeffs, g:G, f:array|seq; work:typeof(newLoopWork(g[0])) = nil) =
  let ds = default(array[0,G])
  c.loopDeriv(g,ds,f,work)

# plaq: 6 types
# rect: 12 types
# pgm: 32=4*2*4=4*3*2+4*2 types, via pgmAction
# shift corners: u[mu],nu mu != nu (12)
# make staples: s[mu][nu] mu != nu (12)
# plaq traces:
#  plaq: U[mu]^+ * sum_{nu!=mu} s[mu][nu] (6)
#  rect: shift(s[mu][nu], nu) (12 shifts)
proc gaugeAction1*[T](c: GaugeActionCoeffs, uu: openarray[T]; work: typeof(newLoopWork(uu[0])) = nil): auto =
  requireFundamental(c, "gaugeAction1")
  mixin mul, redot, load1
  tic("gaugeAction1")
  let u = cast[ptr cArray[T]](unsafeAddr(uu[0]))
  let lo = u[0].l
  let nd = lo.nDim
  #let np = (nd*(nd-1)) div 2
  let nc = u[0][0].ncols
  var cs = startCornerShifts(uu)
  toc("gaugeAction startCornerShifts")
  var
    stf, stu: FieldArray[type(u[0]).V, type(u[0]).T]
    ss: seq[seq[ShiftB[type(uu[0][0])]]]
  if c.rect == 0:
    stf = makeFwdStaples(uu, cs)
  else:
    (stf, stu, ss) = makeStaples(uu, cs)
  toc("gaugeAction makeStaples")
  #var ss = startStapleShifts(st)
  #toc("gaugeAction startStapleShifts")
  let maxThreads = getMaxThreads()
  var nth = 0
  var act = newSeq[float](2*maxThreads)
  toc("gaugeAction setup")
  threads:
    tic()
    var plaq = 0.0
    var rect = 0.0
    for ir in u[0]:
      for mu in 1..<nd:
        for nu in 0..<mu:
          # plaq
          let p1 = redot(u[mu][ir], stf[mu,nu][ir])
          plaq += simdSum(p1)
          if c.rect!=0:
            if isLocal(ss[mu][nu],ir):
              var bmu: type(load1(u[0][0]))
              localSB(ss[mu][nu], ir, assign(bmu,it), stu[mu,nu][ix])
              # rect
              let r = redot(bmu, stf[mu,nu][ir])
              rect += simdSum(r)
            if isLocal(ss[nu][mu],ir):
              var bnu: type(load1(u[0][0]))
              localSB(ss[nu][mu], ir, assign(bnu,it), stu[nu,mu][ix])
              # rect
              let r = redot(bnu, stf[nu,mu][ir])
              rect += simdSum(r)
    toc("gaugeAction local")
    if c.rect != 0:
      for mu in 1..<nd:
        for nu in 0..<mu:
          var needBoundary = false
          boundaryWaitSB(ss[mu][nu]): needBoundary = true
          boundaryWaitSB(ss[nu][mu]): needBoundary = true
          if c.rect != 0 and needBoundary:
            boundarySyncSB()
            for ir in lo:
              if not isLocal(ss[mu][nu],ir):
                var bmu: type(load1(u[0][0]))
                getSB(ss[mu][nu], ir, assign(bmu,it), stu[mu,nu][ix])
                # rect
                let r = redot(bmu, stf[mu,nu][ir])
                rect += simdSum(r)
              if not isLocal(ss[nu][mu],ir):
                var bnu: type(load1(u[0][0]))
                getSB(ss[nu][mu], ir, assign(bnu,it), stu[nu,mu][ix])
                # rect
                let r = redot(bnu, stf[nu,mu][ir])
                rect += simdSum(r)
    act[threadNum*2]   = plaq
    act[threadNum*2+1] = rect
    if threadNum==0: nth = numThreads
    # toc("gaugeAction boundary")
  toc("gaugeAction threads")
  var a = [0.0, 0.0]
  for i in 0..<nth:
    a[0] += act[i*2]
    a[1] += act[i*2+1]
  rankSum(a)
  result = (-1.0/nc.float) * (c.plaq*a[0] + c.rect*a[1])
  if c.pgm != 0:
    result += pgmAction(c, uu, work)
  toc("gaugeAction end")

proc gaugeAction1*[T](uu: openarray[T]): auto =
  let gc = GaugeActionCoeffs(plaq:1.0)
  return gc.gaugeAction1(uu)

proc gaugeActionDeriv*[T](c: GaugeActionCoeffs, uu: openArray[T], f: array|seq, accumulate=false; work: typeof(newLoopWork(uu[0])) = nil) =
  ## if accumulate, the derivatives will add to f.
  ## if not, f is set to 0 first.
  requireFundamental(c, "gaugeActionDeriv")
  mixin load1, adj
  tic("gaugeActionDeriv")
  let u = cast[ptr cArray[T]](unsafeAddr(uu[0]))
  let lo = u[0].l
  let nd = lo.nDim
  #let np = (nd*(nd-1)) div 2
  let nc = u[0][0].ncols
  let cp = c.plaq / float(nc)
  let cr = c.rect / float(nc)
  var cs = startCornerShifts(uu)
  var ru:FieldArray[type(u[0]).V,type(u[0]).T]  # the rect parts of 3
  var sb:seq[seq[ShiftB[type(u[0][0])]]]  # backward ru
  var sf:seq[seq[ShiftB[type(u[0][0])]]]  # forward stf
  if cr!=0:
    ru = newFieldArray2(lo,type(u[0]),[nd,nd],mu!=nu)
    sb.newseq(nd)
    for mu in 0..<nd:
      sb[mu].newseq(nd)
      for nu in 0..<nd:
        if mu==nu: continue
        sb[mu][nu].initShiftB(ru[mu,nu], nu, -1, "all")
    sf.newseq(nd)
    for mu in 0..<nd:
      sf[mu].newseq(nd)
      for nu in 0..<nd:
        if mu==nu: continue
        sf[mu][nu].initShiftB(u[mu], nu, 1, "all")
  toc("init")
  var (stf,stu,ss) = makeStaples(uu, cs)
  toc("makeStaples")
  threads:
    tic("gaugeActionDerivThreads")
    if cr!=0:
      for mu in 1..<nd:
        for nu in 0..<mu:
          sf[mu][nu].startSB(stf[mu,nu][ix])
          sf[nu][mu].startSB(stf[nu,mu][ix])
    for mu in 0..<nd:
      if not accumulate:
        f[mu] := 0
      if cr!=0:
        for nu in 0..<nd:
          if mu!=nu:
            ru[mu,nu] := 0
    for ir in u[0]:
      for mu in 1..<nd:
        for nu in 0..<mu:
          # plaq
          f[mu][ir] += cp * stf[mu,nu][ir]
          f[nu][ir] += cp * stf[nu,mu][ir]
          if isLocal(ss[mu][nu],ir):
            var bmu: type(load1(u[0][0]))
            localSB(ss[mu][nu], ir, assign(bmu,it), stu[mu,nu][ix])
            f[mu][ir] += cp * bmu
            if cr!=0:
              var umu,unu,bmunu: type(load1(u[0][0]))
              getSB(cs[nu][mu], ir, assign(unu,it), u[nu][ix])
              getSB(cs[mu][nu], ir, assign(umu,it), u[mu][ix])
              bmunu := bmu * unu
              f[nu][ir] += cr * bmunu * umu.adj
              ru[nu,mu][ir] += bmu.adj * u[nu][ir] * umu
              ru[mu,nu][ir] += u[nu][ir].adj * bmunu
          if isLocal(ss[nu][mu],ir):
            var bnu: type(load1(u[0][0]))
            localSB(ss[nu][mu], ir, assign(bnu,it), stu[nu,mu][ix])
            f[nu][ir] += cp * bnu
            if cr!=0:
              var unu,umu,bnumu: type(load1(u[0][0]))
              getSB(cs[mu][nu], ir, assign(umu,it), u[mu][ix])
              getSB(cs[nu][mu], ir, assign(unu,it), u[nu][ix])
              bnumu := bnu * umu
              f[mu][ir] += cr * bnumu * unu.adj
              ru[mu,nu][ir] += bnu.adj * u[mu][ir] * unu
              ru[nu,mu][ir] += u[mu][ir].adj * bnumu
          if cr!=0:
            if isLocal(sf[mu][nu],ir):
              var smu,unu,smunu: type(load1(u[0][0]))
              localSB(sf[mu][nu], ir, assign(smu,it), stf[mu,nu][ix])
              getSB(cs[nu][mu], ir, assign(unu,it), u[nu][ix])
              smunu := smu * unu.adj
              f[mu][ir] += cr * u[nu][ir] * smunu
              f[nu][ir] += cr * u[mu][ir] * smunu.adj
              ru[nu,mu][ir] += u[mu][ir].adj * u[nu][ir] * smu
            if isLocal(sf[nu][mu],ir):
              var snu,umu,snumu: type(load1(u[0][0]))
              localSB(sf[nu][mu], ir, assign(snu,it), stf[nu,mu][ix])
              getSB(cs[mu][nu], ir, assign(umu,it), u[mu][ix])
              snumu := snu * umu.adj
              f[nu][ir] += cr * u[mu][ir] * snumu
              f[mu][ir] += cr * u[nu][ir] * snumu.adj
              ru[mu,nu][ir] += u[nu][ir].adj * u[mu][ir] * snu
    toc("local")
    for mu in 1..<nd:
      for nu in 0..<mu:
        var needBoundary = false
        boundaryWaitSB(ss[mu][nu]): needBoundary = true
        boundaryWaitSB(ss[nu][mu]): needBoundary = true
        if needBoundary:
          boundarySyncSB()
          for ir in lo:
            if not isLocal(ss[mu][nu],ir):
              var bmu: type(load1(u[0][0]))
              getSB(ss[mu][nu], ir, assign(bmu,it), stu[mu,nu][ix])
              f[mu][ir] += cp * bmu
              if cr!=0:
                var umu,unu,bmunu: type(load1(u[0][0]))
                getSB(cs[nu][mu], ir, assign(unu,it), u[nu][ix])
                getSB(cs[mu][nu], ir, assign(umu,it), u[mu][ix])
                bmunu := bmu * unu
                f[nu][ir] += cr * bmunu * umu.adj
                ru[nu,mu][ir] += bmu.adj * u[nu][ir] * umu
                ru[mu,nu][ir] += u[nu][ir].adj * bmunu
            if not isLocal(ss[nu][mu],ir):
              var bnu: type(load1(u[0][0]))
              getSB(ss[nu][mu], ir, assign(bnu,it), stu[nu,mu][ix])
              f[nu][ir] += cp * bnu
              if cr!=0:
                var unu,umu,bnumu: type(load1(u[0][0]))
                getSB(cs[mu][nu], ir, assign(umu,it), u[mu][ix])
                getSB(cs[nu][mu], ir, assign(unu,it), u[nu][ix])
                bnumu := bnu * umu
                f[mu][ir] += cr * bnumu * unu.adj
                ru[mu,nu][ir] += bnu.adj * u[mu][ir] * unu
                ru[nu,mu][ir] += u[mu][ir].adj * bnumu
    if cr!=0:
      for mu in 1..<nd:
        for nu in 0..<mu:
          var needBoundary = false
          boundaryWaitSB(sf[mu][nu]): needBoundary = true
          boundaryWaitSB(sf[nu][mu]): needBoundary = true
          if needBoundary:
            boundarySyncSB()
            for ir in lo:
              if not isLocal(sf[mu][nu],ir):
                var smu,unu,smunu: type(load1(u[0][0]))
                getSB(sf[mu][nu], ir, assign(smu,it), stf[mu,nu][ix])
                getSB(cs[nu][mu], ir, assign(unu,it), u[nu][ix])
                smunu := smu * unu.adj
                f[mu][ir] += cr * u[nu][ir] * smunu
                f[nu][ir] += cr * u[mu][ir] * smunu.adj
                ru[nu,mu][ir] += u[mu][ir].adj * u[nu][ir] * smu
              if not isLocal(sf[nu][mu],ir):
                var snu,umu,snumu: type(load1(u[0][0]))
                getSB(sf[nu][mu], ir, assign(snu,it), stf[nu,mu][ix])
                getSB(cs[mu][nu], ir, assign(umu,it), u[mu][ix])
                snumu := snu * umu.adj
                f[nu][ir] += cr * u[mu][ir] * snumu
                f[mu][ir] += cr * u[nu][ir] * snumu.adj
                ru[mu,nu][ir] += u[nu][ir].adj * u[mu][ir] * snu
          threadBarrier()
          sb[mu][nu].startSB(ru[mu,nu][ix])
          sb[nu][mu].startSB(ru[nu,mu][ix])
      toc("staple boundary")
      for ir in u[0]:
        for mu in 1..<nd:
          for nu in 0..<mu:
            if isLocal(sb[mu][nu],ir):
              var b: type(load1(u[0][0]))
              localSB(sb[mu][nu], ir, assign(b,it), ru[mu,nu][ix])
              f[mu][ir] += cr * b
            if isLocal(sb[nu][mu],ir):
              var b: type(load1(u[0][0]))
              localSB(sb[nu][mu], ir, assign(b,it), ru[nu,mu][ix])
              f[nu][ir] += cr * b
      toc("back rect local")
      for mu in 1..<nd:
        for nu in 0..<mu:
          var needBoundary = false
          boundaryWaitSB(sb[mu][nu]): needBoundary = true
          boundaryWaitSB(sb[nu][mu]): needBoundary = true
          if needBoundary:
            boundarySyncSB()
            for ir in lo:
              if not isLocal(sb[mu][nu],ir):
                var b: type(load1(u[0][0]))
                getSB(sb[mu][nu], ir, assign(b,it), ru[mu,nu][ix])
                f[mu][ir] += cr * b
              if not isLocal(sb[nu][mu],ir):
                var b: type(load1(u[0][0]))
                getSB(sb[nu][mu], ir, assign(b,it), ru[nu,mu][ix])
                f[nu][ir] += cr * b
  if c.pgm != 0:
    gaugeLoopDeriv(GaugeActionCoeffs(pgm: c.pgm), uu, uu, f, false, work)
  toc("end")

proc gaugeForce*[T](c: GaugeActionCoeffs, uu: openArray[T], f: array|seq; work: typeof(newLoopWork(uu[0])) = nil) =
  tic("gaugeForce")
  gaugeActionDeriv(c, uu, f, work=work)
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

proc gaugeDeriv2*[G:array|seq](c: GaugeActionCoeffs, g: G, f: array|seq; work: typeof(newLoopWork(g[0])) = nil) =
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
  threads:
    for mu in 0..<nd:
      #let mu = (mux + 1) mod nd
      #f[mu] := 0
      for nu in 0..<nd:
        if nu==mu: continue
        discard t[nu] ^* g[mu]
        shiftExpr(t[mu].sb, f[mu][ir] += cp * t[nu].field[ir]*adj(it), g[nu][ix])
        f[mu] += cp * td[nu] ^* t[mu] ^* g[nu]
        if cr != 0:
          discard t2[nu] ^* t[nu] ^* g[mu]
          discard tg[nu] ^* g[nu]
          shiftExpr(t[mu].sb, f[mu][ir] += cr * t2[nu].field[ir]*adj(it), tg[nu].field[ix])
          f[mu] += cr * td2[nu] ^* td[nu] ^* tg[mu] ^* t[nu] ^* g[nu]
          f[mu] += cr * td2[mu] ^* td[nu] ^* tg[mu] ^* t[mu] ^* g[nu]
          discard td[nu] ^* tg[mu] ^* t[mu] ^* g[nu]
          shiftExpr(t2[mu].sb, f[mu][ir] += cr * td[nu].field[ir]*adj(it), g[mu][ix])
          discard td[mu] ^* t[nu] ^* tg[mu] ^* g[mu]
          shiftExpr(t2[mu].sb, f[mu][ir] += cr * td[mu].field[ir]*adj(it), g[nu][ix])
          shiftExpr(t2[mu].sb, f[mu][ir] += cr * t[nu].field[ir]*adj(it), t[mu].field[ix])
  if c.pgm != 0:
    gaugeLoopDeriv(GaugeActionCoeffs(pgm: -c.pgm), g, g, f, false, work)
  toc("end")

proc gaugeDerivSubset[G:array|seq](c: GaugeActionCoeffs, g: G, f: array|seq, parity, dir: int, clear: static bool; work: typeof(newLoopWork(g[0])) = nil) =
  # Full derivative into workspace scratch, then select the requested links.
  let
    w = useLoopWork(g, work)
    d = w.loopScratch
    sub = g[0].l.getSubset(if parity == 0: "even" else: "odd")
  gaugeActionDeriv(-1.0*c, g, d, work=w)
  threads:
    for x in f[dir]:
      if x >= sub.lowOuter and x < sub.highOuter:
        f[dir][x] := d[dir][x]
      else:
        when clear: f[dir][x] := 0

proc gaugeDeriv2SubsetWork*[G:array|seq](c: GaugeActionCoeffs, g: G, f: array|seq, sd, sf, sb: auto, parity, dir: int, clear: static bool; work: typeof(newLoopWork(g[0])) = nil) =
  ## f[dir]|P = D(g)|P; other directions are unchanged.
  ## clear=false preserves f[dir]|~P; clear=true sets it to zero.
  ## rect/pgm coefficients take the full-derivative route and ignore sd, sf, sb.
  requireFundamental(c, "gaugeDeriv2SubsetWork")
  if c.rect != 0 or c.pgm != 0:
    gaugeDerivSubset(c, g, f, parity, dir, clear, work)
    return
  mixin adj
  tic("gaugeDeriv2Subset")
  let lo = g[0].l
  let nd = lo.nDim
  const nc = g[0][0].nrows
  let cp = -c.plaq / float(nc)
  when clear:
    let other = lo.getSubset(if parity == 0: "odd" else: "even")
  let firstNu = if dir == 0: 1 else: 0
  toc("gaugeDeriv2Subset setup")
  threads:
    when clear:
      for x in other:
        f[dir][x] := 0
    for nu in 0..<nd:
      if nu == dir: continue

      discard sd ^*! g[nu]
      threadBarrier()
      if nu == firstNu:
        shiftExpr(sf[nu], f[dir][ir] := cp * (g[nu][ir] * it) * adj(sd.field[ir]), g[dir][ix])
      else:
        shiftExpr(sf[nu], f[dir][ir] += cp * (g[nu][ir] * it) * adj(sd.field[ir]), g[dir][ix])
      # := and += may use different shift partitions.
      threadBarrier()
      toc("gaugeDeriv2Subset forward")
      shiftExpr(sb[nu], f[dir][ir] += cp*it, g[nu][ix].adj * (g[dir][ix] * sd.field[ix]))
      threadBarrier()
      toc("gaugeDeriv2Subset backward")
  toc("gaugeDeriv2Subset end")

proc gaugeDeriv2Subset*[G:array|seq](c: GaugeActionCoeffs, g: G, f: array|seq, parity, dir: int; work: typeof(newLoopWork(g[0])) = nil) =
  ## f[dir]|P = D(g)|P; other directions are unchanged.
  if c.rect != 0 or c.pgm != 0:
    gaugeDerivSubset(c, g, f, parity, dir, true, work)
    return
  let ps = if parity == 0: "even" else: "odd"
  let sd = newShifter(g[0], dir, 1)
  let sf = createShiftBufs(g[0], 1, ps)
  let sb = createShiftBufs(g[0], -1, ps)
  c.gaugeDeriv2SubsetWork(g, f, sd, sf, sb, parity, dir, true, work=work)

proc gaugeDerivDeriv2*[G:array|seq](c: GaugeActionCoeffs, g: G, h, f: array|seq; work: typeof(newLoopWork(g[0])) = nil) =
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
  threads:
    for mu in 0..<nd:
      for nu in 0..<nd:
        if nu==mu: continue
        discard t[nu] ^* g[mu]
        shiftExpr(t[mu].sb, f[mu][ir] += cp * t[nu].field[ir]*adj(it), h[nu][ix])
        discard t[nu] ^* h[mu]
        shiftExpr(t[mu].sb, f[mu][ir] += cp * t[nu].field[ir]*adj(it), g[nu][ix])
        discard th[nu] ^* g[mu]
        shiftExpr(t[mu].sb, f[mu][ir] += cp * th[nu].field[ir]*adj(it), g[nu][ix])
        f[mu] += cp * td[nu] ^* t[mu] ^* h[nu]
        f[mu] += cp * td[nu] ^* th[mu] ^* g[nu]
        f[mu] += cp * thd[nu] ^* t[mu] ^* g[nu]
  if c.rect != 0 or c.pgm != 0:
    gaugeLoopDeriv(GaugeActionCoeffs(rect: -c.rect, pgm: -c.pgm), g, h, f, true, work)
  toc("end")

proc gaugeDerivDeriv2SubsetImpl[G:array|seq](c: GaugeActionCoeffs, g: G, hs, hdir: auto, f: array|seq, parity, dir: int, sum, add, addBase: static bool; work: typeof(newLoopWork(g[0])) = nil) =
  # S=(parity,dir); sum: hdir|P=sum(hs)|P; add: f+=H(hdir|P).
  # addBase: f[S]+=H(hdir|P); f[~S]=base[~S]+H(hdir|P).
  requireFundamental(c, "gaugeDerivDeriv2Subset")
  if c.rect != 0 or c.pgm != 0:
    # Apply the full Hessian to a masked seed kept in workspace scratch.
    let
      w = useLoopWork(g, work)
      hm = w.loopScratch
      sub = g[0].l.getSubset(if parity == 0: "even" else: "odd")
    when addBase:
      let other = g[0].l.getSubset(if parity == 0: "odd" else: "even")
    threads:
      for mu in 0..<g.len:
        hm[mu] := 0
      threadBarrier()
      for x in sub:
        when sum:
          hdir[x] := hs[0][x]
          for k in 1..<hs.len:
            hdir[x] += hs[k][x]
        hm[dir][x] := hdir[x]
      threadBarrier()
      when addBase:
        for mu in 0..<g.len:
          if mu != dir:
            f[mu] := hs[mu]
        for x in other:
          f[dir][x] := hs[dir][x]
      elif not add:
        for mu in 0..<g.len:
          f[mu] := 0
    c.gaugeDerivDeriv2(g, hm, f, work=w)
    return
  mixin adj
  tic("gaugeDerivDeriv2Subset")
  let lo = g[0].l
  let nd = lo.nDim
  const nc = g[0][0].nrows
  let cp = -c.plaq / float(nc)
  let ps = if parity == 0: "even" else: "odd"
  let po = if parity == 0: "odd" else: "even"
  when sum:
    let sub = lo.getSubset(ps)
  when addBase:
    let other = lo.getSubset(po)
    let firstNu = if dir == 0: 1 else: 0
  let sd = newShifter(g[0], dir, 1)
  let sfSame = createShiftBufs(g[0], 1, ps)
  let sfOther = createShiftBufs(g[0], 1, po)
  let sb = createShiftBufs(g[0], -1, po)
  let sq = g[0].newOneOf
  let bqSame = createShiftB(g[0], dir, -1, ps)
  let bqOther = createShiftB(g[0], dir, -1, po)
  toc("gaugeDerivDeriv2Subset setup")
  threads:
    when sum:
      for x in sub:
        hdir[x] := hs[0][x]
        for k in 1..<hs.len:
          hdir[x] += hs[k][x]
      threadBarrier()
      toc("gaugeDerivDeriv2Subset sum")
    when addBase:
      if nd == 1:
        for x in other:
          f[dir][x] := hs[dir][x]
    elif not add:
      f[dir] := 0

    for nu in 0..<nd:
      if nu == dir: continue

      discard sd ^*! g[nu]
      threadBarrier()
      template put(i, y: untyped) =
        when addBase:
          f[nu][i] := hs[nu][i] + y
        elif add:
          f[nu][i] += y
        else:
          f[nu][i] := y
      template setOtherAdd(i, t: untyped) =
        put(i, cp * (g[dir][i] * sd.field[i]) * adj(t))
        sq[i] := g[nu][i] * t
        f[dir][i] += cp * sq[i] * adj(sd.field[i])
      when addBase:
        template setOtherBase(i, t: untyped) =
          put(i, cp * (g[dir][i] * sd.field[i]) * adj(t))
          sq[i] := g[nu][i] * t
          f[dir][i] := hs[dir][i] + cp * sq[i] * adj(sd.field[i])
      template setSame(i, t: untyped) =
        sd.field[i] := hdir[i] * sd.field[i]
        put(i, cp * sd.field[i] * adj(t))
        sq[i] := g[nu][i] * t
      when addBase:
        if nu == firstNu:
          shiftExpr(sfOther[nu], setOtherBase(ir, it), hdir[ix])
        else:
          shiftExpr(sfOther[nu], setOtherAdd(ir, it), hdir[ix])
      else:
        shiftExpr(sfOther[nu], setOtherAdd(ir, it), hdir[ix])
      shiftExpr(sfSame[nu], setSame(ir, it), g[dir][ix])
      threadBarrier()
      toc("gaugeDerivDeriv2Subset forward")

      shiftExpr(bqSame, f[nu][ir] += cp*it, g[dir][ix].adj * sq[ix])
      shiftExpr(bqOther, f[nu][ir] += cp*it, hdir[ix].adj * sq[ix])
      shiftExpr(sb[nu], f[dir][ir] += cp*it, g[nu][ix].adj * sd.field[ix])
      threadBarrier()
      toc("gaugeDerivDeriv2Subset backward")
  toc("gaugeDerivDeriv2Subset end")

proc gaugeDerivDeriv2Subset*[G:array|seq](c: GaugeActionCoeffs, g: G, h, f: array|seq, parity, dir: int; work: typeof(newLoopWork(g[0])) = nil) =
  ## f = H_g(h[dir]|P), including all affected links.
  c.gaugeDerivDeriv2SubsetImpl(g, h, h[dir], f, parity, dir, false, false, false, work=work)

proc gaugeDerivDeriv2SubsetAdd*[G:array|seq](c: GaugeActionCoeffs, g: G, hdir: auto, f: array|seq, parity, dir: int; work: typeof(newLoopWork(g[0])) = nil) =
  ## f += H_g(hdir|P).
  c.gaugeDerivDeriv2SubsetImpl(g, hdir, hdir, f, parity, dir, false, true, false, work=work)

proc gaugeDerivDeriv2SubsetAddBase*[G:array|seq](c: GaugeActionCoeffs, g: G, hdir: auto, base, f: array|seq, parity, dir: int; work: typeof(newLoopWork(g[0])) = nil) =
  ## S=(P,dir): f[S] += H_g(hdir|P)[S].
  ## f[~S] = base[~S] + H_g(hdir|P)[~S].
  c.gaugeDerivDeriv2SubsetImpl(g, base, hdir, f, parity, dir, false, false, true, work=work)

proc gaugeDerivDeriv2SubsetSum*[G:array|seq](c: GaugeActionCoeffs, g: G, h: array|seq, w: auto, f: array|seq, parity, dir: int; work: typeof(newLoopWork(g[0])) = nil) =
  ## w|P = sum(h)|P; f = H_g(w|P).
  c.gaugeDerivDeriv2SubsetImpl(g, h, w, f, parity, dir, true, false, false, work=work)

proc gaugeForce2*[G:array|seq](c: GaugeActionCoeffs, g: G, f: array|seq; work: typeof(newLoopWork(g[0])) = nil) =
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
  threads:
    for mu in 0..<nd:
      #let mu = (mux + 1) mod nd
      f[mu] := 0
      for nu in 0..<nd:
        if nu==mu: continue
        discard t[nu] ^* g[mu]
        shiftExpr(t[mu].sb, f[mu][ir] += cp * t[nu].field[ir]*adj(it), g[nu][ix])
        f[mu] += cp * td[nu] ^* t[mu] ^* g[nu]
        if cr != 0:
          discard t2[nu] ^* t[nu] ^* g[mu]
          discard tg[nu] ^* g[nu]
          shiftExpr(t[mu].sb, f[mu][ir] += cr * t2[nu].field[ir]*adj(it), tg[nu].field[ix])
          f[mu] += cr * td2[nu] ^* td[nu] ^* tg[mu] ^* t[nu] ^* g[nu]
          f[mu] += cr * td2[mu] ^* td[nu] ^* tg[mu] ^* t[mu] ^* g[nu]
          discard td[nu] ^* tg[mu] ^* t[mu] ^* g[nu]
          shiftExpr(t2[mu].sb, f[mu][ir] += cr * td[nu].field[ir]*adj(it), g[mu][ix])
          discard td[mu] ^* t[nu] ^* tg[mu] ^* g[mu]
          shiftExpr(t2[mu].sb, f[mu][ir] += cr * td[mu].field[ir]*adj(it), g[nu][ix])
          threadBarrier()
          shiftExpr(t2[mu].sb, f[mu][ir] += cr * t[nu].field[ir]*adj(it), t[mu].field[ix])
    for mu in 0..<nd:
      for e in f[mu]:
        let s = g[mu][e] * f[mu][e].adj
        f[mu][e].projectTAH s
  toc("end")
proc gaugeForce2*(f,g: array|seq) =
  var c = GaugeActionCoeffs(plaq:1.0)
  gaugeForce2(c,g,f)

proc gaugeAction3*[G:array|seq](c: GaugeActionCoeffs, g: G; work: typeof(newLoopWork(g[0])) = nil): auto =
  requireFundamental(c, "gaugeAction3")
  tic("gaugeAction3")
  const nc = g[0][0].nrows
  let lo = g[0].l
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
      let ws = g.wilsonLines ls
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

proc gaugeForce3*(c: GaugeActionCoeffs, g,f: auto; work: typeof(newLoopWork(g[0])) = nil) =
  requireFundamental(c, "gaugeForce3")
  if c.pgm != 0:
    c.gaugeForce(g, f, work=work)
    return
  tic("gaugeForce3")
  const nc = g[0][0].nrows
  let nd = g[0].l.nDim
  let cp = c.plaq / nc.float
  let cr = c.rect / nc.float
  threads:
    for mu in 0..<nd:
      f[mu] := 0
  for mu in 1..<nd:
    for nu in 0..<mu:
      let ptree = c.plaqRectPath(mu,nu)
      let ws = g.gaugeProd ptree
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
    for mu in 0..<nd:
      for e in f[mu]:
        let s = g[mu][e] * f[mu][e].adj
        f[mu][e].projectTAH s
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

proc action*(c: GaugeActionCoeffs, g: auto; work: typeof(newLoopWork(g[0])) = nil): float =
  ## Dispatch by coefficient family; mixed adjoint/improved terms are rejected.
  if c.adjplaq != 0: c.actionA(g)
  else: c.gaugeAction1(g, work=work)

proc force*(c: GaugeActionCoeffs, g, f: auto; work: typeof(newLoopWork(g[0])) = nil) =
  ## Projected force for action(c,g). Keep f disjoint from g.
  if c.adjplaq != 0: c.forceA(g, f)
  else: c.gaugeForce(g, f, work=work)

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
    mixin exp
    for mu in 0..<g.len:
      for e in g[mu]:
        let t = exp(eps*p[mu][e])*g[mu][e]
        g[mu][e] := t
      #echo "g[", mu, "]: ", g[mu].norm2

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
