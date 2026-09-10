## Covariant transport (hop) and Wilson lines.
##
## `hop` is the closed transport primitive: one covariant parallel-transport
## step of a general matrix field,
##
##   hop(g, f, mu, +1)(x) = U_mu(x) * f(x+mu)
##   hop(g, f, mu, -1)(x) = U_mu(x-mu)^dag * f(x-mu)
##
## Its backward closes on the basic set: the field cotangent is the reverse
## hop (transport is unitary, so the adjoint is the return trip) and the
## gauge cotangent is a single-direction product injected into slot mu.
##
## `wilsonLine`/`transport` chain hops, so path-ordered link products and all
## their derivatives (to any order) come from the basic-op closure.

import ../core
import ../support/op
import layout, physics/qcdTypes
import types, basic_ops, field_ops

type Ghop = ref object of Gfield
  ## Owns the fused shift-multiply comm buffers (cloned by newOneOf).
  tr: Transporter[DLatticeColorMatrixV, DLatticeColorMatrixV, DColorMatrixV]
  mu, sgn: int

proc hopNodeLike(x: Gfield, mu, sgn: int): Ghop =
  # Node storage aliases the transporter's receive buffer, so the fused
  # shift-multiply writes the result directly into this node's value;
  # newOneOf clones them as a pair. The initial link binding is shape only:
  # hopf rebinds it from the gauge input on every evaluation.
  result = Ghop(
    runtime: x.runtime,
    tr: newTransporter(x.fval, x.fval, mu, sgn),
    mu: mu,
    sgn: sgn)
  result.fval = result.tr.field
  result.fval.zeroFieldStorage
  result.assignStableNodeId

method newOneOf(x: Ghop): Gvalue =
  hopNodeLike(x, x.mu, x.sgn)

proc hop*(g: Ggauge, x: Gfield, mu: int, sgn: int): Gfield

proc hopb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let h = Ghop(z)
  let g = Ggauge(z.inputs[0])
  let x = Gfield(z.inputs[1])
  let u = requireUpstream(zb, "hop backward", Gfield)
  if i == 0:
    if h.sgn > 0:
      # z(x) = U(x) f(x+mu): Ubar(x) = u(x) f(x+mu)^dag
      return injectLink(u * shift(x, h.mu, 1).adj, h.mu, g)
    # z(x) = U(x-mu)^dag f(x-mu): Ubar(y) = f(y) u(y+mu)^dag
    return injectLink(x * shift(u, h.mu, 1).adj, h.mu, g)
  # Transport is unitary: the field cotangent is the reverse hop.
  hop(g, u, h.mu, -h.sgn)

proc hopf(v: Gvalue) =
  let g = Ggauge(v.inputs[0])
  let x = Gfield(v.inputs[1])
  let z = Ghop(v)
  z.tr.setLink g.gval[z.mu]
  threads:
    discard z.tr ^* x.fval

let hopg = Gfunc(forward: hopf, backward: hopb, name: "hop")

proc hop*(g: Ggauge, x: Gfield, mu: int, sgn: int): Gfield =
  ## One covariant transport step of x by the mu links of g (see module doc).
  if sgn != 1 and sgn != -1:
    raiseValueError("hop sign must be +1 or -1")
  g.requireLinkShape(mu, x.fval, "hop")
  if mu >= x.fval.l.nDim:
    raiseValueError("hop direction out of range")
  graphNode(hopNodeLike(x, mu, sgn), @[Gvalue(g), Gvalue(x)], hopg, "hop")

proc transport*(g: Ggauge, x: Gfield, path: openArray[int]): Gfield =
  ## Path-ordered transport of x along `path`. Entries are +-(mu+1); the
  ## first entry is the first link of the line as seen from the base site:
  ##   transport(g, x, [p1, p2, ..., pk]) = hop_p1(hop_p2(... hop_pk(x)))
  for p in path:
    if p == 0 or p < -g.gval.len or p > g.gval.len:
      raiseValueError(
        "transport path entry must be in [-" & $g.gval.len & ",-1] or [1," &
        $g.gval.len & "]; got " & $p)
  result = x
  for k in countdown(path.high, 0):
    let p = path[k]
    result = hop(g, result, abs(p) - 1, (if p > 0: 1 else: -1))

proc wilsonLine*(g: Ggauge, path: openArray[int]): Gfield =
  ## Path-ordered product of links along `path` starting at each site, e.g.
  ## wilsonLine(g, [mu+1, nu+1, -(mu+1), -(nu+1)])(x) is the mu-nu plaquette
  ##   U_mu(x) U_nu(x+mu) U_mu(x+nu)^dag U_nu(x)^dag.
  transport(g, g.unitFieldLike, path)

proc plaqPath*(mu, nu: int): array[4, int] =
  if mu < 0 or nu < 0 or mu == nu:
    raiseValueError("plaquette directions must be distinct and nonnegative")
  [mu + 1, nu + 1, -(mu + 1), -(nu + 1)]
