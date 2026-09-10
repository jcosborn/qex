## Halo-based lattice moves of a single field: no interior copy, boundary
## exchange only, with the halo layout narrowed to the offset's directions.
##
##   gather(f, sh)(x)  = f(x - sh)                        adjoint: scatter
##   scatter(f, sh)(y) = f(y + sh)                        adjoint: gather
##   gp(a, b, sh, fa, fb)(x) = A(x) B(x - sh),  A = adjIf(a, fa), B = adjIf(b, fb)
##     a-bar = adjIf(gp(u, b, sh, false, not fb), fa)     (= adjIf(u B^dag, fa))
##     b-bar = scatter(adjIf(adjIf(a, not fa) u, fb), sh) (= scatter(adjIf(A^dag u, fb), sh))
## under the pairing dL = Re tr(G^dag dF). The three ops close on each other
## and the site algebra, so a lineProducts step is one gp node.
##
## A node owns its halo buffers and index table and rebinds the halo's field
## to its input on every evaluation (clone rule, DESIGN section 6); layouts
## and maps come from the halo cache. Halo updates run outside `threads:`.

import ../core
import ../support/op
import layout, gauge, physics/qcdTypes
import comms/[commsTypes, halo]
import types, basic_ops

type
  Lo = Layout[VLEN]
  FieldHalo = Halo[Lo, DLatticeColorMatrixV, DColorMatrixV]
  Gmove = ref object of Gfield
    ## Shared plumbing of gather, scatter, and gp: the layout/map for one
    ## offset and the extended index of x - sh for every local site.
    sh: seq[int]
    hl: HaloLayout[Lo]
    hm: HaloMap[Lo]
    idx: seq[int32]
    h: FieldHalo          # over the gathered input (gather, gp) or the output (scatter)
    fa, fb: bool          # gp flags

proc moveNodeLike(x: Gfield, sh: seq[int], onOutput: bool): Gmove =
  ## Halo widths only along the offset's directions: x - sh needs bck width
  ## sh_d for sh_d > 0 and fwd width -sh_d for sh_d < 0.
  let lo = x.fval.l
  let nd = lo.nDim
  var fwd, bck = newSeq[int32](nd)
  var off = newSeq[int32](nd)
  for d in 0 ..< min(nd, sh.len):
    off[d] = int32(-sh[d])
    if sh[d] > 0: bck[d] = int32(sh[d]) else: fwd[d] = int32(-sh[d])
  let hl = haloLayout(lo, fwd, bck)
  result = Gmove(runtime: x.runtime, sh: sh, hl: hl,
                 hm: haloMap(hl, getDefaultComm(), @[off]))
  result.fval = x.fval.newOneOf
  result.fval.zeroFieldStorage
  result.h = makeHalo(hl, (if onOutput: result.fval else: x.fval))
  result.idx = newSeq[int32](hl.nOut)
  for x in 0 ..< hl.nOut:
    var i = int32 x
    for d in 0 ..< min(nd, sh.len):
      for k in 1 .. abs(sh[d]):
        i = (if sh[d] > 0: hl.neighborBck[d][i] else: hl.neighborFwd[d][i])
    result.idx[x] = i
  result.assignStableNodeId

proc gather*(f: Gfield, sh: seq[int]): Gfield
proc scatter*(f: Gfield, sh: seq[int]): Gfield
proc gp*(a, b: Gfield, sh: seq[int], fa = false, fb = false): Gfield

proc gatherf(v: Gvalue) =
  let z = Gmove(v)
  let f = Gfield(v.inputs[0])
  z.h.field = f.fval
  z.h.update(z.hm, getDefaultComm())
  threads:
    for x in z.fval:
      z.fval[x] := z.h[z.idx[x]]

proc gatherb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  scatter(requireUpstream(zb, "gather backward", Gfield), Gmove(z).sh)

let gatherg = Gfunc(forward: gatherf, backward: gatherb, name: "gather")

proc scatterf(v: Gvalue) =
  # Every local x has one target x - sh: local targets are written directly,
  # shell targets go to their owners through the reverse exchange.
  let z = Gmove(v)
  let f = Gfield(v.inputs[0])
  threads:
    # The reverse exchange accumulates, so the output and the shell start
    # from zero on every evaluation.
    z.fval := 0.0
    tfor c, 0 ..< z.h.halo.len:
      z.h.halo[c] := 0.0
    threadBarrier()
    for x in z.fval:
      let i = z.idx[x]
      if i < z.hl.nOut:
        z.fval[i] := f.fval[x]
      else:
        z.h.halo[i - z.hl.nOut] := f.fval[x]
  z.h.updateRev(z.hm, getDefaultComm())

proc scatterb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  gather(requireUpstream(zb, "scatter backward", Gfield), Gmove(z).sh)

let scatterg = Gfunc(forward: scatterf, backward: scatterb, name: "scatter")

proc gpf(v: Gvalue) =
  let z = Gmove(v)
  let a = Gfield(v.inputs[0])
  let b = Gfield(v.inputs[1])
  z.h.field = b.fval
  z.h.update(z.hm, getDefaultComm())
  let fa = z.fa
  let fb = z.fb
  threads:
    for x in z.fval:
      var t {.noinit.}: evalType(z.fval[x])
      t := z.h[z.idx[x]]
      if fa:
        if fb: z.fval[x] := a.fval[x].adj * t.adj
        else: z.fval[x] := a.fval[x].adj * t
      else:
        if fb: z.fval[x] := a.fval[x] * t.adj
        else: z.fval[x] := a.fval[x] * t

proc gpb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let n = Gmove(z)
  let a = Gfield(z.inputs[0])
  let b = Gfield(z.inputs[1])
  let u = requireUpstream(zb, "gp backward", Gfield)
  if i == 0:
    let t = gp(u, b, n.sh, false, not n.fb)
    return (if n.fa: t.adj else: t)
  var t = (if n.fa: a else: a.adj) * u
  if n.fb: t = t.adj
  scatter(t, n.sh)

let gpg = Gfunc(forward: gpf, backward: gpb, name: "gp")

method newOneOf(x: Gmove): Gvalue =
  var r = moveNodeLike(x, x.sh, x.gfunc == scatterg)
  r.fa = x.fa
  r.fb = x.fb
  r

proc gather*(f: Gfield, sh: seq[int]): Gfield =
  ## f(x - sh) through a halo of f.
  graphNode(moveNodeLike(f, sh, false), @[Gvalue(f)], gatherg, "gather")

proc scatter*(f: Gfield, sh: seq[int]): Gfield =
  ## f(y + sh), the adjoint of gather, through the reverse exchange.
  graphNode(moveNodeLike(f, sh, true), @[Gvalue(f)], scatterg, "scatter")

proc gp*(a, b: Gfield, sh: seq[int], fa = false, fb = false): Gfield =
  ## adjIf(a, fa)(x) * adjIf(b, fb)(x - sh) in one node: the product of a
  ## Wilson-line step, with b gathered rather than copied.
  a.requireSameFieldShape(b, "gp")
  var n = moveNodeLike(b, sh, false)
  n.fa = fa
  n.fb = fb
  graphNode(n, @[Gvalue(a), Gvalue(b)], gpg, "gp")
