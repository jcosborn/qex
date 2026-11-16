import std/unittest
import gauge/loopPlan

const families = [{lkPlaq}, {lkRect}, {lkPgm}, {lkPlaq,lkRect},
  {lkPlaq,lkPgm}, {lkRect,lkPgm}, {lkPlaq,lkRect,lkPgm}]

suite "Loop plans":
  test "schedules initialize values and reuse only retired slots":
    for nd in 2..4:
      for kinds in families:
        for action in [false,true]:
          for order in 0..5:
            if action and order != 0: continue
            var active = kinds
            if nd == 2: active.excl lkPgm
            if order > 3: active.excl lkPlaq
            if active == {}: continue
            checkpoint "dimensions=" & $nd & " kinds=" & $active & " action=" & $action & " order=" & $order
            let p = gaugeLoops(nd, active, action).loopProducts(order)
            var ready = newSeq[bool](p.size)
            var last = newSeq[int](p.size)
            var owner = newSeq[int](p.slots)
            for x in owner.mitems: x = -1
            for op in p.loads:
              check not ready[op.dest]
              check op.seed in 0..order and op.mu in 0..<nd
              ready[op.dest] = true
              check p.slot[op.dest] == -1
            for pos, k in p.sched:
              if k >= 0:
                last[p.steps[k].left] = pos
                last[p.steps[k].right] = pos
              else: last[p.outs[-k-1].source] = pos
            var output = 0
            var seen = newSeq[bool](nd)
            for pos, k in p.sched:
              if k >= 0:
                let op = p.steps[k]
                for v in [op.left,op.right]:
                  check ready[v]
                  if p.slot[v] >= 0: check owner[p.slot[v]] == v
                let slot = p.slot[op.dest]
                check slot in 0..<p.slots
                if op.add:
                  check ready[op.dest] and owner[slot] == op.dest
                else:
                  check not ready[op.dest] and owner[slot] == -1
                  ready[op.dest] = true
                  owner[slot] = op.dest
                for v in [op.left,op.right]:
                  if p.slot[v] >= 0 and last[v] == pos: owner[p.slot[v]] = -1
              else:
                check -k-1 == output
                inc output
                let op = p.outs[-k-1]
                check ready[op.source] and op.kind in active
                check op.first == not seen[op.mu]
                seen[op.mu] = true
                if p.slot[op.source] >= 0:
                  check owner[p.slot[op.source]] == op.source
                  if last[op.source] == pos: owner[p.slot[op.source]] = -1
            check output == p.outs.len
            for x in owner: check x == -1

  test "shared products reduce work for active loop families":
    for nd in 2..4:
      let planes = nd*(nd-1) div 2
      let triples = nd*(nd-1)*(nd-2) div 6
      for kinds in families:
        let np = if lkPlaq in kinds: planes else: 0
        let nr = if lkRect in kinds: 2*planes else: 0
        let ng = if lkPgm in kinds: 4*triples else: 0
        if np+nr+ng == 0: continue
        let a = gaugeLoops(nd, kinds, true).loopProducts(0)
        let d = gaugeLoops(nd, kinds).loopProducts(0)
        let h = gaugeLoops(nd, kinds).loopProducts(1)
        check a.outs.len == np+nr+ng
        check a.steps.len <= 2*np+4*(nr+ng)
        check d.outs.len == 4*np+6*(nr+ng)
        check d.steps.len < 2*4*np+4*6*(nr+ng)
        check h.steps.len < 5*4*np+11*6*(nr+ng)
        check d.slots < d.size-d.loads.len

  test "top jets load only their seeds":
    for kinds in [{lkPlaq}, {lkRect}, {lkPgm}, {lkRect,lkPgm}]:
      let order = if kinds == {lkPlaq}: 3 else: 5
      let p = gaugeLoops(4,kinds).loopProducts(order)
      for op in p.loads: check op.seed > 0
      check p == gaugeLoops(4,kinds).loopProducts(order)

  test "path operands reject zero directions and empty paths":
    expect ValueError: discard optimalPairs([@[1,0,2]])
    expect ValueError: discard optimalPairs([newSeq[int]()]).plan
