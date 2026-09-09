## Gauge fields on the 16-cell honeycomb.
##
## Sites Z^4 + (Z+1/2)^4 are a hypercubic lattice of cells (a QEX `Layout`)
## with sublattices A (cell y) and B (y + 1/2).  24 link fields per cell:
##   uA[mu]     A(y) -> A(y+e_mu)         mu = 0..3
##   uB[mu]     B(y) -> B(y+e_mu)
##   uD[delta]  B(y) -> A(y+delta)        delta in {0,1}^4, direction delta - 1/2
## `links` holds the same 24 fields flat (uA, uB, uD order) for QEX's
## openArray routines.  Sites per cell: 2; links per cell: 24.
##
## HcShift16 builds all 16 copies f(y +- delta) with 15 single-axis shifts
## (level mu: f[delta or 2^mu](y) = f[delta](y + sign e_mu)).  Multi-axis
## `makeShiftSubQ` displacements mis-permute SIMD lanes; do not use them.
##
## Configuration files: one LIME file with metadata <hcGauge> (geom, beta,
## traj, info) and one record of the 24 fields in `links` order, cell
## geometry as the lattice size, so plain QEX readers order the sites.

import std/[math, strutils, os]
import base, layout, field, maths, rng, io
import physics/qcdTypes
import gauge
import hcgeom

export hcgeom

# ---------------------------------------------------------------------------
# 16-way binary-tree shift
# ---------------------------------------------------------------------------

func topBit*(delta: int): int {.inline.} =
  ## index of the highest set bit of delta > 0
  var d = delta
  result = -1
  while d != 0:
    inc result
    d = d shr 1

type
  HcShift16*[F, S] = object
    ## f[delta] = src(y + sign*delta); f[0] aliases src, f[1..15] are the
    ## shifter buffers.  S is the Shifter type (not derivable from F).
    f*: array[nDiag, F]
    sh*: array[nDiag, S]
    sign*: int

proc newHcShift16*[F](src: F; sign: int = 1): auto =
  ## allocates 15 shifters; call outside `threads:`
  type S = type(newShifter(src, 0, 1))
  var r: HcShift16[F, S]
  r.sign = sign
  r.f[0] = src
  for d in 1..<nDiag:
    r.sh[d] = newShifter(src, topBit(d), sign)
    r.f[d] = r.sh[d].field
  r

proc run*(s: var HcShift16) =
  ## refresh f[1..15] from f[0]; inside `threads:`
  for d in 1..<nDiag:
    let mu = topBit(d)
    discard s.sh[d] ^* s.f[d and not (1 shl mu)]

proc setSrc*[F, S](s: var HcShift16[F, S], src: F) =
  ## ref assignment; outside `threads:`
  s.f[0] = src

template `[]`*(s: HcShift16, delta: int): untyped = s.f[delta]

# ---------------------------------------------------------------------------
# gauge field
# ---------------------------------------------------------------------------

type
  HcGauge*[F] = object
    uA*, uB*: array[nDim, F]
    uD*: array[nDiag, F]
    links*: seq[F]       ## uA[0..3], uB[0..3], uD[0..15]

template lo*(g: HcGauge): untyped = g.uA[0].l

proc setLinks[F](g: var HcGauge[F]) =
  g.links = newSeq[F](nDirs)
  for mu in 0..<nDim:
    g.links[mu] = g.uA[mu]
    g.links[nDim+mu] = g.uB[mu]
  for d in 0..<nDiag:
    g.links[2*nDim+d] = g.uD[d]

proc newHcGauge*(lo: Layout, nc: static[int] = getDefaultNc()): auto =
  ## 24 unit link fields on the cell layout `lo`; allocates
  type F = type(lo.ColorMatrix(nc))
  var g: HcGauge[F]
  for mu in 0..<nDim:
    g.uA[mu] = lo.ColorMatrix(nc)
    g.uB[mu] = lo.ColorMatrix(nc)
  for d in 0..<nDiag:
    g.uD[d] = lo.ColorMatrix(nc)
  g.setLinks
  for u in g.links: u := 1
  g

proc newOneOf*[F](g: HcGauge[F]): HcGauge[F] =
  ## same shape, contents undefined; allocates
  for mu in 0..<nDim:
    result.uA[mu] = g.uA[mu].newOneOf
    result.uB[mu] = g.uB[mu].newOneOf
  for d in 0..<nDiag:
    result.uD[d] = g.uD[d].newOneOf
  result.setLinks

# The following are allocation free and run inside a `threads:` block.

proc unit*(g: HcGauge) =
  for u in g.links: u := 1

proc random*(g: HcGauge, r: var RNGField) =
  random(g.links, r)

proc warm*(g: HcGauge, s: float, r: var RNGField) =
  warm(g.links, s, r)

proc reunit*(g: HcGauge) =
  when g.uA[0][0].nrows == 1:
    projectU(g.links)
  else:
    projectSU(g.links)

proc checkSU*(g: HcGauge): tuple[avg, max: float] =
  checkSU(g.links)

proc randomTAH*(g: HcGauge, r: RNGField) =
  ## gaussian traceless anti-Hermitian matrices on all 24 fields
  for u in g.links: randomTAH(u, r)

proc link*[F](g: HcGauge[F], k: LinkKind, idx: int): F =
  ## the field of hcgeom.LinkRef (kind, idx)
  case k
  of lkA: g.uA[idx]
  of lkB: g.uB[idx]
  of lkD: g.uD[idx]

proc `:=`*(a: HcGauge, b: HcGauge) =
  for i in 0..<nDirs:
    a.links[i] := b.links[i]

proc redot*(a, b: HcGauge): float =
  ## sum over all 24 fields and sites of Re tr(a^dag b); outside `threads:`
  var res = 0.0
  threads:
    var s = 0.0
    for i in 0..<nDirs:
      s += redot(a.links[i], b.links[i])
    threadMaster: res = s
  res

proc norm2*(g: HcGauge): float =
  ## sum over all 24 fields and sites of |u|^2; inside `threads:`
  for i in 0..<nDirs:
    result += norm2(g.links[i])

proc norm2diff*(a, b: HcGauge): float =
  ## sum_l |a_l - b_l|^2; inside `threads:`
  for i in 0..<nDirs:
    result += norm2diff(a.links[i], b.links[i])

proc expMul*[F](r, g, p: HcGauge[F], t: float) =
  ## r := exp(t p) g on all 24 fields (the HMC link update); r may alias g.
  ## Outside `threads:`.
  threads:
    for i in 0..<nDirs:
      let u = r.links[i]
      let u0 = g.links[i]
      let pu = p.links[i]
      for e in u:
        u[e] := exp(t*pu[e])*u0[e]

template toF*(x: untyped): float =
  ## a single-site lane accessor as a float (QEX `f{i}` returns proxies)
  block:
    var v: float
    v := x
    v

template setC*(dest: untyped; a, b: float) =
  ## dest := a + i b on a lane accessor
  dest.re := a
  dest.im := b

# ---------------------------------------------------------------------------
# triangle loops
# ---------------------------------------------------------------------------

proc triangleTrace*[F](g: HcGauge[F]): tuple[re, im: float] =
  ## sum_x sum_{i=1..32} Tr P_i(x) over all sites, each triangle once at its
  ## apex.  B apex at y: uD[d](y) uA[mu](y+d) uD[d'](y)^dag, d' = d or 2^mu.
  ## A apex, re-based at z = y-db (db = d xor 15, db' = db-2^mu):
  ## uB[mu](z) uD[db'](z+e_mu) uD[db](z)^dag.  Allocates shifters.
  type SH = type(newHcShift16(g.uA[0], 1))
  var shA: array[nDim, SH]
  for mu in 0..<nDim:
    shA[mu] = newHcShift16(g.uA[mu], 1)
  type SS = type(newShifter(g.uD[0], 0, 1))
  var sD: array[nDim, SS]
  for mu in 0..<nDim:
    sD[mu] = newShifter(g.uD[0], mu, 1)
  var m = g.uA[0].newOneOf
  var tr: type(trace(m))
  threads:
    m := 0
    for mu in 0..<nDim:
      shA[mu].run
    threadBarrier()
    for t in apexTris:
      m += (g.uD[t.delta] * shA[t.mu].f[t.delta]) * g.uD[t.deltaP].adj
    for t in apexTris:
      let
        db = t.delta xor 15
        dbp = t.deltaP xor 15
        s = sD[t.mu] ^* g.uD[dbp]
      m += (g.uB[t.mu] * s) * g.uD[db].adj
    threadBarrier()
    tr = trace(m)
  (tr.re, tr.im)

proc triangleSum*(g: HcGauge): float =
  ## (1/(32 N_sites)) sum_x sum_i Re Tr P_i(x)/N; 1 for the unit gauge.
  ## The imaginary part is not identically zero (fixed apex orientation).
  const nc = g.uA[0][0].nrows
  triangleTrace(g).re/float(nTriPerSite*2*g.lo.physVol*nc)

# ---------------------------------------------------------------------------
# gauge transformations
# ---------------------------------------------------------------------------

proc gaugeTransform*[F](g: HcGauge[F], vA, vB: F) =
  ## u -> V(start) u V(end)^dag on all 24 fields; vA on A sites, vB on B.
  ## Allocates shifters.
  var shA = newHcShift16(vA, 1)
  type SS = type(newShifter(vB, 0, 1))
  var sB: array[nDim, SS]
  for mu in 0..<nDim:
    sB[mu] = newShifter(vB, mu, 1)
  var t = vA.newOneOf
  threads:
    shA.run
    for mu in 0..<nDim:
      discard sB[mu] ^* vB
    threadBarrier()
    for mu in 0..<nDim:
      t := vA * g.uA[mu]
      g.uA[mu] := t * shA.f[1 shl mu].adj
      t := vB * g.uB[mu]
      g.uB[mu] := t * sB[mu].field.adj
    for d in 0..<nDiag:
      t := vB * g.uD[d]
      g.uD[d] := t * shA.f[d].adj

# ---------------------------------------------------------------------------
# configuration files
# ---------------------------------------------------------------------------

type
  HcMeta* = object
    version*: int
    geom*: seq[int]
    beta*: float
    traj*: int
    info*: string

const recordMd = "<?xml version=\"1.0\"?><note>16-cell honeycomb gauge: " &
  "uA[0..3], uB[0..3], uD[0..15]</note>"

proc metaXml(m: HcMeta): string =
  var gs = newSeq[string](m.geom.len)
  for i in 0..<m.geom.len: gs[i] = $m.geom[i]
  "<?xml version=\"1.0\"?>\n<hcGauge>\n" &
    "  <version>" & $m.version & "</version>\n" &
    "  <geom>" & gs.join(" ") & "</geom>\n" &
    "  <beta>" & $m.beta & "</beta>\n" &
    "  <traj>" & $m.traj & "</traj>\n" &
    "  <info>" & m.info & "</info>\n" &
    "</hcGauge>\n"

proc xmlTag(s, tag: string): string =
  let a = s.find("<" & tag & ">")
  if a < 0: return ""
  let b = s.find("</" & tag & ">", a)
  if b < 0: return ""
  s[a + tag.len + 2 ..< b].strip

proc parseMeta*(xml: string): HcMeta =
  result.version = -1
  result.traj = -1
  let v = xmlTag(xml, "version")
  if v.len > 0: result.version = parseInt(v)
  for t in xmlTag(xml, "geom").splitWhitespace:
    result.geom.add parseInt(t)
  let b = xmlTag(xml, "beta")
  if b.len > 0: result.beta = parseFloat(b)
  let tr = xmlTag(xml, "traj")
  if tr.len > 0: result.traj = parseInt(tr)
  result.info = xmlTag(xml, "info")

proc save*(g: HcGauge, fn: string, beta = 0.0, traj = -1, info = ""): int =
  ## returns 0 on success; outside `threads:`
  let m = HcMeta(version: 1, geom: g.lo.physGeom, beta: beta, traj: traj,
                 info: info)
  saveGauge(g.links, fn, "", metaXml(m), recordMd)

proc load*(g: HcGauge, fn: string): tuple[status: int, meta: HcMeta] =
  ## `g` must have the file's cell layout; returns (0, metadata) on success
  var rd = g.lo.newReader(fn)
  if rd.status != 0: return (rd.status, HcMeta())
  let fmd = rd.fileMetadata
  rd.read(g.links)
  if rd.status != 0: return (rd.status, HcMeta())
  rd.close()
  if rd.status != 0: return (rd.status, HcMeta())
  let m = parseMeta(fmd)
  if m.version != 1:
    qexWarn "unexpected hcGauge metadata version in ", fn, ": ", fmd
  (0, m)

proc fileGeom*(fn: string): seq[int] =
  ## cell geometry stored in a configuration file; empty if unreadable
  if fileExists(fn): getFileLattice(fn) else: @[]

template withLayout*(simdlen: int, geom: openArray[int], lo, body: untyped) =
  ## run `body` with `lo` a cell layout of SIMD length simdlen
  ## (0 = build default; 1 or 2 for geometries the default rejects)
  case simdlen
  of 0:
    let lo {.inject.} = newLayout(geom)
    body
  of 1:
    let lo {.inject.} = newLayout(geom, 1)
    body
  of 2:
    let lo {.inject.} = newLayout(geom, 2)
    body
  else:
    qexError "unsupported simdlen ", simdlen, " (0 = default, 1, 2)"
