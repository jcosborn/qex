## Configuration I/O for the 16-cell honeycomb gauge field (task **M**).
##
## Format (documented here, pinned by the round-trip test in tests/thmc.nim):
##
##   One SciDAC/LIME file per configuration, written with QEX's standard
##   `Writer` (the same machinery as `saveGauge`), containing
##
##   * file metadata: the XML string below, carrying the cell geometry, beta,
##     the update/trajectory number and a free-form info string:
##
##       <?xml version="1.0"?>
##       <hcGauge>
##         <version>1</version>
##         <geom>ns ns ns nt</geom>
##         <beta>...</beta>
##         <traj>...</traj>
##         <info>...</info>
##       </hcGauge>
##
##   * ONE binary record with datacount = 24: the link fields in the canonical
##     order  uA[0..3], uB[0..3], uD[0..15]  (i.e. `hcgauge.allLinks` order),
##     each a colour matrix per CELL site, written in the field's native
##     precision (double), so a save/load round trip is bit exact.
##
## The lattice dimensions stored in the file are the CELL geometry, so
## `getFileLattice` returns the cell layout size and a plain QEX reader/writer
## does all the site ordering (lexicographic by cell coordinate, rank
## independent).  RNG state is NOT stored (runs record their seed and update
## number; that is enough to document a stream, not to splice it).

import std/[strutils, os]
import base, layout, field, io
import physics/qcdTypes
import gauge
import hcgeom, hclayout, hcgauge

export hcgauge

type
  HcMeta* = object
    version*: int
    geom*: seq[int]
    beta*: float
    traj*: int
    info*: string

const hcRecordMd = "<?xml version=\"1.0\"?><note>16-cell honeycomb gauge: " &
  "uA[0..3], uB[0..3], uD[0..15]</note>"

proc metaXml*(m: HcMeta): string =
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
  ## the text between <tag> and </tag>, or "" if absent
  let a = s.find("<" & tag & ">")
  if a < 0: return ""
  let b = s.find("</" & tag & ">", a)
  if b < 0: return ""
  s[a + tag.len + 2 ..< b].strip

proc parseHcMeta*(xml: string): HcMeta =
  result.version = -1
  let v = xmlTag(xml, "version")
  if v.len > 0:
    try: result.version = parseInt(v)
    except ValueError: discard
  for t in xmlTag(xml, "geom").splitWhitespace:
    try: result.geom.add parseInt(t)
    except ValueError: discard
  let b = xmlTag(xml, "beta")
  if b.len > 0:
    try: result.beta = parseFloat(b)
    except ValueError: discard
  let tr = xmlTag(xml, "traj")
  result.traj = -1
  if tr.len > 0:
    try: result.traj = parseInt(tr)
    except ValueError: discard
  result.info = xmlTag(xml, "info")

proc saveHcGauge*[V: static[int], F](g: HcGauge[V, F], fn: string,
                                     beta = 0.0, traj = -1,
                                     info = ""): int =
  ## Save `g` to `fn` in the format above.  Returns 0 on success.
  ## Call outside `threads:`.
  let m = HcMeta(version: 1, geom: g.hl.geom, beta: beta, traj: traj,
                 info: info)
  let links = allLinks(g)
  saveGauge(links, fn, "", metaXml(m), hcRecordMd)

proc loadHcGauge*[V: static[int], F](g: HcGauge[V, F], fn: string):
    tuple[status: int, meta: HcMeta] =
  ## Load `fn` into `g` (which must already have the matching cell layout;
  ## the reader errors out on a geometry mismatch).  Returns (0, metadata)
  ## on success.  Call outside `threads:`.
  var links = allLinks(g)
  var rd = g.lo.newReader(fn)
  if rd.status != 0: return (rd.status, HcMeta())
  let fmd = rd.fileMetadata
  rd.read(links)
  if rd.status != 0: return (rd.status, HcMeta())
  rd.close()
  if rd.status != 0: return (rd.status, HcMeta())
  var m = parseHcMeta(fmd)
  if m.version != 1:
    qexWarn "hcio: unexpected hcGauge metadata version in ", fn, ": ", fmd
  (0, m)

proc hcFileGeom*(fn: string): seq[int] =
  ## cell geometry stored in a configuration file (empty if unreadable)
  if fileExists(fn): getFileLattice(fn) else: @[]

when isMainModule:
  import std/strformat
  import rng
  import hcaction
  qexInit()
  echo "hcio smoke test: save + load round trip"
  let hl = newHcLayout([4, 4, 4, 6])
  var r = hl.lo.newRNGField(MRG32k3a, 777'u64)
  var g = newHcGauge(hl)
  threads:
    g.random r
  let ts0 = g.triangleSum
  let fn = getTempDir() / "hcio_smoke.lime"
  doAssert 0 == g.saveHcGauge(fn, beta = 6.25, traj = 42, info = "smoke")
  var g2 = newHcGauge(hl)
  let (st, meta) = g2.loadHcGauge(fn)
  doAssert st == 0
  echo &"meta: version {meta.version} geom {meta.geom} beta {meta.beta} traj {meta.traj} info '{meta.info}'"
  let ts1 = g2.triangleSum
  var d2 = 0.0
  block:
    var scr = newOneOf(g)
    threads:
      for mu in 0..<nDim:
        scr.uA[mu] := g.uA[mu] - g2.uA[mu]
        scr.uB[mu] := g.uB[mu] - g2.uB[mu]
      for d in 0..<nDiag:
        scr.uD[d] := g.uD[d] - g2.uD[d]
    d2 = redot(scr, scr)
  echo &"triangleSum {ts0} -> {ts1}  (diff {abs(ts1-ts0):.3e})"
  echo &"sum|U-U'|^2 = {d2:.3e} (must be exactly 0)"
  doAssert d2 == 0.0 and ts1 == ts0
  removeFile fn
  echo "hcio round trip OK"
  qexFinalize()
