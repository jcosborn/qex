## Plain-text TSV with `#key=value` header metadata.
##
##   # lattice=L4
##   # at=0.095238095238095233
##   # columns=t g e
##   0	0.1	0.001
##
## Numbers are written with %.17g, so a write/read round trip is bit identical.
## `columns` is the one reserved key: it carries the column names and is returned
## as `names`, not in `meta`.  A `#` line with no `=` is a comment and is dropped.
## Blank lines are skipped, so a file may carry gnuplot dataset separators.

import std/[os, strformat, strutils, tables]

const colKey = "columns"

proc writeTsv*(path: string, header: openArray[(string, string)],
               colNames: openArray[string], cols: openArray[seq[float]]) =
  ## `cols` is column major: cols[j][i] is row i of column j.  Every column must
  ## have the length of cols[0].
  let dir = path.parentDir
  if dir.len > 0: createDir dir
  var f = open(path, fmWrite)
  defer: f.close()
  for (k, v) in header: f.writeLine("# " & k & "=" & v)
  if colNames.len > 0: f.writeLine("# " & colKey & "=" & colNames.join(" "))
  let nrow = if cols.len == 0: 0 else: cols[0].len
  for i in 0..<nrow:
    var line = ""
    for j in 0..<cols.len:
      if j > 0: line.add '\t'
      line.add &"{cols[j][i]:.17g}"
    f.writeLine line

proc readTsv*(path: string): tuple[meta: Table[string, string], names: seq[string],
                                   cols: seq[seq[float]]] =
  ## Inverse of `writeTsv`.  The column count is set by the first data line.
  result.meta = initTable[string, string]()
  for raw in path.lines:
    let s = raw.strip
    if s.len == 0: continue
    if s[0] == '#':
      let h = s[1..^1].strip
      let i = h.find '='
      if i < 0: continue                      # free-form comment
      let
        k = h[0..<i].strip
        v = h[i+1..^1].strip
      if k == colKey: result.names = v.splitWhitespace
      else: result.meta[k] = v
    else:
      let w = s.splitWhitespace
      if result.cols.len == 0: result.cols.setLen w.len
      for j in 0..<w.len: result.cols[j].add parseFloat(w[j])
