## Integer path factorization and gauge-loop schedules. Build outside threads;
## consumers read the cached plans and own their numerical storage separately.
import std/[hashes, tables, bitops]

type
  Coord = object
    x: seq[int]

proc `+`(x,y:Coord):Coord =
  let
    nx = x.x.len
    ny = y.x.len
    n = max(nx,ny)
    m = min(nx,ny)
  result.x.newseq n
  for i in 0..<m:
    result.x[i] = x.x[i] + y.x[i]
  if nx>ny:
    for i in m..<n:
      result.x[i] = x.x[i]
  else:
    for i in m..<n:
      result.x[i] = y.x[i]

proc `-`(x,y:Coord):Coord =
  let
    nx = x.x.len
    ny = y.x.len
    n = max(nx,ny)
    m = min(nx,ny)
  result.x.newseq n
  for i in 0..<m:
    result.x[i] = x.x[i] - y.x[i]
  if nx>ny:
    for i in m..<n:
      result.x[i] = x.x[i]
  else:
    for i in m..<n:
      result.x[i] = -y.x[i]

proc `-`(x:Coord):Coord =
  let
    n = x.x.len
  result.x.newseq n
  for i in 0..<n:
    result.x[i] = -x.x[i]

proc `+=`(x:var Coord,y:Coord) =
  let
    nx = x.x.len
    ny = y.x.len
  if nx<ny:
    x.x.setLen ny
  for i in 0..<ny:
    x.x[i] += y.x[i]

proc `[]`(x:Coord, i:SomeInteger):int =
  if i>=x.x.len:
    result = 0
  else:
    result = x.x[i]

proc `[]`(x:var Coord, i:SomeInteger):var int =
  if i>=x.x.len:
    x.x.setLen(i+1)
  return x.x[i]

proc `[]=`(x:var Coord, i:SomeInteger, d:int) =
  if i>=x.x.len:
    x.x.setLen(i+1)
  x.x[i] = d

type
  OrdPathKind = enum
    opInt, opPair, opList, opAdj
  OrdPath = ref object
    case k: OrdPathKind
    of opInt: d:int
    of opPair: l,r:OrdPath
    of opList: s:seq[OrdPath]
    of opAdj: p:OrdPath
  OrdPathTree* = ref object
    paths: seq[OrdPath]
    segments: seq[OrdPath]
    counts: seq[int]

proc `$`(path:OrdPath):string =
  return
    case path.k:
    of opInt: $path.d
    of opPair: "(" & $path.l & " " & $path.r & ")"
    of opList: $path.s
    of opAdj: "Adj" & $path.p

proc newOrdPath(path:openarray[int]):OrdPath =
  let n = path.len
  var op = OrdPath(k:opList, s:newseq[OrdPath](n))
  for i in 0..<n:
    let d = path[i]
    if d==0:
      raise newException(ValueError, "path directions must be nonzero: " & $path)
    op.s[i] = OrdPath(k:opInt, d:d)
  return op

proc adjoint(p:OrdPath):OrdPath =
  case p.k:
  of opInt:
    return OrdPath(k:opInt, d: -p.d)
  of opPair:
    return OrdPath(k:opPair, l:p.r.adjoint, r:p.l.adjoint)
  of opList:
    var
      l = p.s.len
      l1 = l-1
      n = newseq[OrdPath](l)
    for i in 0..<l:
      n[i] = p.s[l1-i].adjoint
    return OrdPath(k:opList, s:n)
  of opAdj:
    return p.p

proc adjointOp(p:OrdPath):OrdPath =
  ## Like adjoint(), but create new OrdPath(k:opAdj) if it's not already an opInt or opAdj.
  case p.k:
  of opInt:
    return p.adjoint
  of opAdj:
    return p.p
  else:
    return OrdPath(k:opAdj, p:p)

proc flatten(path:OrdPath):seq[int] =
  case path.k:
  of opInt:
    return @[path.d]
  of opPair:
    var r = path.l.flatten
    r.add path.r.flatten
    return r
  of opList:
    var r = newseq[int]()
    for l in path.s:
      r.add l.flatten
    return r
  of opAdj:
    return path.p.adjoint.flatten

proc deltaX(path:OrdPath):Coord =
  # return the Coord difference between the end of the path to the beginning.
  case path.k:
  of opInt:
    let d = abs(path.d)-1
    if path.d>0:
      result[d] = 1
    else:
      result[d] = -1
  of opPair:
    result = path.l.deltaX + path.r.deltaX
  of opList:
    for l in path.s:
      result += l.deltaX
  of opAdj:
    result = -path.p.deltaX

proc position(path:OrdPath):Coord =
  # return the Coord relative to the left most starting point.
  case path.k:
  of opInt:
    if path.d<0:
      result[-1-path.d] = -1
  of opPair:
    result = path.l.position
  of opList:
    if path.s.len>0:
      result = path.s[0].position
  of opAdj:
    result = path.p.position - path.p.deltaX

proc `$`*(path:OrdPathTree):string =
  result = "OrdPathTree:\n  paths:"
  for p in path.paths:
    result.add "\n    " & $p.flatten & " \t " & $p & " \t shift " & $p.position.x
  result.add "\n  segments:"
  for i in 0..<path.segments.len:
    let
      l = path.segments[i].l
      r = path.segments[i].r
    result.add "\n    " & $path.counts[i] & "  " & $l & " \t " & $r & " \t shift " & $(l.position - l.deltaX - r.position).x

proc hash(path:OrdPath):Hash =
  return
    case path.k:
    of opInt: !$(opInt.hash !& path.d.hash)
    of opPair: !$(opPair.hash !& path.l.hash !& path.r.hash)
    of opList: !$(opList.hash !& path.s.hash)
    of opAdj: !$(opAdj.hash !& path.p.hash)

proc `==`(x,y:OrdPath):bool {.noSideEffect.} =  # noSideEffect needed for Nim devel
  # strictly follows hash(path:OrdPath)
  return
    if x.k != y.k:
      false
    else:
      case x.k:
      of opInt: x.d == y.d
      of opPair: x.l == y.l and x.r == y.r
      of opList: x.s == y.s
      of opAdj: x.p == y.p

proc lexLess(a, b: seq[int]): bool =
  for i in 0..<min(a.len, b.len):
    if a[i] != b[i]: return a[i] < b[i]
  a.len < b.len

proc mostSharedPair(paths:openarray[OrdPath]):(OrdPath,int) =
  ## Receives paths and search each element of OrdPath(k:opList).
  ## Return an OrdPath(k:opPair) occured most frequently among all the paths,
  ## and its count.
  ## If there are no pairs, return OrdPath(k:opList) with len 0, and 0.
  var pc = initTable[OrdPath,int]()
  var p,pa: OrdPath
  var fp,fpa: seq[int]
  var c,ca: int
  while true:
    c = 0
    ca = 0
    for ps in paths:
      if not (ps.k==opList and ps.s.len>1):
        continue
      var i = 0
      while i<ps.s.len-1:
        i.inc
        let t = OrdPath(k:opPair, l:ps.s[i-1], r:ps.s[i])
        let ft = t.flatten
        if c==0:
          p = t
          pa = OrdPath(k:opPair, l:t.r.adjointOp, r:t.l.adjointOp)
          if p in pc or pa in pc:
            continue
          fp = p.flatten()
          fpa.newseq(fp.len)
          for i in 0..<fp.len:
            fpa[i] = -fp[fp.len-i-1]
          c = 1
          i.inc
        elif ft == fp:
          c.inc
          i.inc
        elif ft == fpa:
          ca.inc
          i.inc
    if c==0:
      # no new pairs
      break
    else:
      let ct = c+ca
      if c>=ca:
        pc[p] = ct
      else:
        pc[pa] = ct
  c = 0
  var fp0: seq[int]
  for k,v in pc.pairs:
    # Ties resolve to the smallest flattened path, so the plan is deterministic.
    let fk = k.flatten
    if v>c or (v==c and fk.lexLess(fp0)):
      c = v
      p = k
      fp0 = fk
  if c==0:
    return (OrdPath(k:opList, s: @[]), 0)
  else:
    return (p, c)

proc groupByPair(paths:openarray[OrdPath], pair:OrdPath):seq[OrdPath] =
  ## Receives paths and group each element of OrdPath(k:opList) by the pair OrdPath(k:opPair).
  let
    fp = pair.flatten
    pairAdj = OrdPath(k:opAdj, p:pair)
  var afp = newseq[int](fp.len)
  for i in 0..<fp.len:
    afp[i] = -fp[fp.len-i-1]
  var newpaths = newseq[OrdPath]()
  for ps in paths:
    if ps.k!=opList:
      newpaths.add ps
      continue
    var
      p = newseq[OrdPath]()
      i = 0
      n = ps.s.len
    while i<n-1:
      var t = ps.s[i].flatten
      t.add ps.s[i+1].flatten
      if t == fp:
        p.add pair
        # echo "groupByPair: add ",p[^1]
        i.inc 2
      elif t==afp:
        p.add pairAdj
        # echo "groupByPair: add ",p[^1]
        i.inc 2
      else:
        p.add ps.s[i]
        i.inc
    if i==n-1:
      p.add ps.s[i]
    newpaths.add(
      if p.len==1:
        p[0]
      else:
        OrdPath(k:opList, s:p)
    )
  return newpaths

proc optimalPairs(paths:openarray[OrdPath]):OrdPathTree =
  ## Runtime O(N^3) with N=paths.len
  var
    pgroup = @paths
    segments = newseq[OrdPath]()
    counts = newseq[int]()
    (p,c) = pgroup.mostSharedPair
  while p.k==opPair:
    # echo "optimalPairs: add ",p
    segments.add p
    counts.add c
    pgroup = pgroup.groupByPair p
    (p,c) = pgroup.mostSharedPair
  return OrdPathTree(paths:pgroup, segments:segments, counts:counts)

proc optimalPairs*(paths:openarray[seq[int]]):OrdPathTree =
  ## Runtime O(N^3) with N=paths.len
  let n = paths.len
  var ps = newseq[OrdPath](n)
  for i,p in paths.pairs:
    ps[i] = p.newOrdPath
  ps.optimalPairs

type
  PathStep* = object
    ## nodes[key](x) = L(x) * R(x - sh),  L = adjIf(nodes[l], la),  R = adjIf(nodes[r], ra)
    key*, l*, r*: seq[int]
    la*, ra*: bool
    sh*: seq[int]
  PathOut* = object
    ## adjIf(nodes[key], adj)(x - sh)
    key*: seq[int]
    adj*: bool
    sh*: seq[int]
  PathPlan* = object
    ## Evaluation plan over nodes keyed by flattened paths.
    ## A key of length one is the raw link mu+1 (never negative, the sign is in la/ra/adj);
    ## longer keys are step results.  Steps are in evaluation order, so every operand
    ## is a link or an earlier step.  sh lists per-direction offsets, trailing zeros dropped.
    steps*: seq[PathStep]
    outs*: seq[PathOut]

proc operand(s:OrdPath):(seq[int],bool) =
  # Node key and adjoint flag of a step operand or an output path.
  case s.k:
  of opInt:
    return (@[abs s.d], s.d<0)
  of opPair:
    return (s.flatten, false)
  of opAdj:
    let (k,a) = operand s.p
    return (k, not a)
  of opList:
    if s.s.len==1:
      return operand s.s[0]
    raise newException(ValueError, "expected a path operand: " & $s)

proc trim(c:Coord):seq[int] =
  result = c.x
  while result.len>0 and result[^1]==0:
    result.setLen(result.len-1)

proc plan*(t:OrdPathTree, origin=true, shifts=true):PathPlan =
  ## Steps in the order gaugeProd evaluates them, then one output per input path.
  ## origin=true shifts each output back to the path's starting site.
  ## shifts=false treats integers as independent matrix symbols at one site.
  result.steps.newseq t.segments.len
  for i,s in t.segments.pairs:
    let
      (l,la) = operand s.l
      (r,ra) = operand s.r
    result.steps[i] = PathStep(key:s.flatten, l:l, r:r, la:la, ra:ra)
    if shifts:
      result.steps[i].sh = trim(s.l.position - s.l.deltaX - s.r.position)
  result.outs.newseq t.paths.len
  for i,p in t.paths.pairs:
    let (k,a) = operand p
    result.outs[i] = PathOut(key:k, adj:a)
    if origin and shifts:
      result.outs[i].sh = trim(-p.position)

func firstSetBit(x: SomeInteger): int {.inline.} =  # bitops doesn't handle x==0 correctly
  result = 0
  if x != 0:
    result = bitops.firstSetBit(x)

type
  LoopKind* = enum lkPlaq, lkRect, lkPgm
  GaugeLoop* = object
    path*: array[5,int]
    len*: int
    kind*: LoopKind
  GaugeLoops* = seq[seq[GaugeLoop]]
  GaugeLoopPlan* = ref object
    paths*: GaugeLoops
    fw*, bw*: seq[int32]
    offsets*: seq[seq[seq[int32]]]  # base sites of the links of each direction
    products: array[6,LoopProducts]
    factors: seq[tuple[mu:int, offset:seq[int32]]]
    pairs: PathPlan
    outputs: seq[tuple[mu:int, kind:LoopKind]]
  LoopLoad* = object
    dest*, mu*, seed*: int
    walk*: seq[tuple[d, n: int]]  # nonzero components of the factor's offset
  LoopMul* = object
    dest*, left*, right*: int
    la*, ra*, add*: bool
  LoopOut* = object
    source*, mu*: int
    kind*: LoopKind
    adj*, first*: bool  # first: the first output of mu in sched
  LoopProducts* = ref object
    loads*: seq[LoopLoad]
    steps*: seq[LoopMul]
    outs*: seq[LoopOut]
    size*: int  # values: loads and products
    sched*: seq[int]  # k >= 0: steps[k]; k < 0: outs[-k-1], each after its source
    slot*: seq[int]  # scratch slot of each product, -1 for a load
    slots*: int

proc gaugeLoops*(nd: int, kinds: set[LoopKind], action = false): GaugeLoopPlan =
  # Each open path runs from x to x+mu. Differentiating a closed loop
  # contributes every occurrence, including repeated links after wrapping.
  proc makeLoops(): GaugeLoops =
    result.newSeq(nd)
    template addLoop(p: untyped, family: LoopKind) =
      block:
        let path = p
        for k in 0..(if action: 0 else: path.len-1):
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
        if lkPlaq in kinds:
          addLoop([mu,nu,-mu,-nu], lkPlaq)
        if lkRect in kinds:
          addLoop([mu,nu,nu,-mu,-nu,-nu], lkRect)
          addLoop([mu,mu,nu,-mu,-mu,-nu], lkRect)
        if lkPgm in kinds:
          for sg in 1..<nu:
            # gaugeAction2's ts1, ts2, ts3, ts7, respectively.
            addLoop([mu,nu,sg,-mu,-nu,-sg], lkPgm)
            addLoop([mu,sg,nu,-mu,-sg,-nu], lkPgm)
            addLoop([nu,mu,sg,-nu,-mu,-sg], lkPgm)
            addLoop([mu,-nu,sg,-mu,nu,-sg], lkPgm)
  # Paths and offsets depend only on topology. Construct this metadata and its
  # derived product plans outside threads; evaluations share the cached plans.
  var cache {.global.}: Table[(int,set[LoopKind],bool),GaugeLoopPlan]
  let key = (nd, kinds, action)
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

proc loopProducts*(plan: GaugeLoopPlan, order: int): LoopProducts =
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

