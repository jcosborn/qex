## The 16-cell honeycomb as a QEX layout (task **L**).
##
## The honeycomb is a hypercubic lattice of **cells** (`N_s^3 x N_t`, an ordinary
## QEX `Layout`) with a two-site basis:
##
##   sublattice A (sub = 0) at cell `y`,
##   sublattice B (sub = 1) at `y + (1/2,1/2,1/2,1/2)`.
##
## See doc/FORMULATION.md section 1.2.  Nothing in QEX understands non-hypercubic
## geometry, so everything here is built on `Layout` plus **single-axis** shifts.
##
## The other half of this module is `HcShift16`, the binary-tree 16-way shift:
## the diagonal neighbours need all 16 copies `f(y +- delta)`, `delta in {0,1}^4`,
## and this builds them with **15 single-axis shifts** instead of 15 general
## gathers (PLAN.md section 1.1).  Each individual shift is a plain QEX `Shifter`,
## hence SIMD-safe and MPI-safe; `makeShiftSubQ` with a multi-axis displacement is
## *not* (it mis-permutes SIMD lanes) and must not be used.

import base, layout, field, maths
import physics/qcdTypes
import hcgeom

export hcgeom

type
  HcLayout*[V: static[int]] = ref object
    lo*: Layout[V]        ## the CELL layout, N_s^3 x N_t
    ns*, nt*: int
    geom*: seq[int]       ## the full cell geometry

proc newHcLayoutX*(geom: openArray[int], V: static[int]): HcLayout[V] =
  ## Cell layout from an explicit 4d geometry.  `V` is the SIMD length; use
  ## `V = 1` for a non-vectorised layout (handy when isolating lane bugs).
  doAssert geom.len == nDim, "honeycomb cell geometry must be 4 dimensional"
  result = HcLayout[V](ns: geom[0], nt: geom[nDim-1], geom: @geom)
  result.lo = newLayout(geom, V)

template newHcLayout*(geom: openArray[int]): untyped =
  newHcLayoutX(geom, VLEN)
template newHcLayout*(geom: openArray[int]; V: static[int]): untyped =
  newHcLayoutX(geom, V)
template newHcLayout*(ns, nt: int): untyped =
  newHcLayoutX([ns, ns, ns, nt], VLEN)
template newHcLayout*(ns, nt: int; V: static[int]): untyped =
  newHcLayoutX([ns, ns, ns, nt], V)

template nCells*(hl: HcLayout): int = hl.lo.physVol
  ## total number of cells (global)
template nSites*(hl: HcLayout): int = 2*hl.lo.physVol
  ## total number of honeycomb sites (global) = 2 per cell
template nLinks*(hl: HcLayout): int = nDirs*hl.lo.physVol
  ## total number of links (global) = 24 per cell = 12 per site

# ---------------------------------------------------------------------------
# the 16-way binary-tree shift
# ---------------------------------------------------------------------------

func topBit*(delta: int): int {.inline.} =
  ## index of the highest set bit of `delta` (delta > 0)
  var d = delta
  result = -1
  while d != 0:
    inc result
    d = d shr 1

type
  HcShift16*[F, S] = object
    ## `f[delta]` holds `src(y + sign*delta)` for all 16 `delta in {0,1}^4`.
    ##
    ## `f[0]` **aliases** `src`; `f[1..15]` are the private buffers of the 15
    ## shifters.  Two type parameters, not one: `F` is the field type and `S` the
    ## `Shifter[F,T]` type, which Nim cannot derive from `F` alone.
    f*: array[nDiag, F]
    sh*: array[nDiag, S]  ## sh[0] is unused
    sign*: int            ## +1 (forward, y+delta) or -1 (backward, y-delta)

proc newHcShift16*[F](src: F; sign: int = 1): auto =
  ## Allocate the 15 shifters of the binary tree.  Call **outside** `threads:`.
  ##
  ##   level 0:  f[0]                       = src               (delta = 0000)
  ##   level mu: f[delta or 2^mu](y) = f[delta](y + sign*e_mu)   for delta < 2^mu
  ##
  ## `sign = +1` gives `f[delta](y) = src(y+delta)` (the B -> A diagonal
  ## neighbours, `uD[delta]: B(y) -> A(y+delta)`);
  ## `sign = -1` gives `f[delta](y) = src(y-delta)` (the A -> B ones,
  ## `A(y) -> B(y-dbar)`).
  doAssert sign == 1 or sign == -1, "HcShift16 sign must be +-1"
  type S = type(newShifter(src, 0, 1))
  var r: HcShift16[F, S]
  r.sign = sign
  r.f[0] = src
  for d in 1..<nDiag:
    let mu = topBit(d)
    r.sh[d] = newShifter(src, mu, sign)
    r.f[d] = r.sh[d].field
  r

proc run*(s: var HcShift16) =
  ## Refresh `f[1..15]` from `f[0]`.  Call inside a `threads:` block.
  ## The usual pattern is to bind the source once and call this whenever the
  ## contents of the source field change.
  for d in 1..<nDiag:
    let mu = topBit(d)
    let p = d and not (1 shl mu)
    discard s.sh[d] ^* s.f[p]

proc setSrc*[F, S](s: var HcShift16[F, S], src: F) =
  ## Point the tree at a different source field.  This is a plain reference
  ## assignment, so call it **outside** a `threads:` block.
  s.f[0] = src

proc run*[F, S](s: var HcShift16[F, S], src: F) =
  ## Rebind the source and refresh.  Call **outside** a `threads:` block (it is
  ## `setSrc` + `run`; wrap the `run` in `threads:` yourself for the threaded
  ## version).
  s.setSrc src
  s.run

template `[]`*(s: HcShift16, delta: int): untyped = s.f[delta]
