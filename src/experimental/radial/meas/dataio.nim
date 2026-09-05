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

import std/[os, strformat, strutils, tables, tempfiles]

const colKey = "columns"

proc readTsvMeta*(path: string): Table[string, string] =
  ## Read only the leading metadata header.  This is used by restart/analysis
  ## compatibility checks without parsing a potentially large data body.
  result = initTable[string, string]()
  for raw in path.lines:
    let s = raw.strip
    if s.len == 0: continue
    if s[0] != '#': break
    let h = s[1..^1].strip
    let i = h.find '='
    if i < 0: continue
    let
      k = h[0..<i].strip
      v = h[i+1..^1].strip
    if k != colKey: result[k] = v

proc requireTsvMeta*(path: string, expected: openArray[(string, string)]) =
  ## Extra descriptive metadata is allowed; every identifying field must match.
  let meta = readTsvMeta(path)
  for (k, v) in expected:
    if not meta.hasKey(k) or meta[k] != v:
      let got = if meta.hasKey(k): meta[k] else: "<missing>"
      raise newException(ValueError, &"{path}: incompatible {k}: expected {v}, got {got}; " &
        "remeasure with -skipDone:false or use a fresh output directory")

proc outputsDone*(paths: openArray[string], expected: openArray[(string, string)]): bool =
  ## An interrupted group is incomplete. Check all existing members for conflicts.
  result = true
  for path in paths:
    if fileExists(path): requireTsvMeta(path, expected)
    else: result = false

proc writeTsv*[T: float|string = float](path: string, header: openArray[(string, string)],
               colNames: openArray[string], cols: openArray[seq[T]]) =
  ## `cols` is column major: cols[j][i] is row i of column j.  Every column must
  ## have the length of cols[0]. String cells are preformatted fields and must
  ## contain no tabs or line breaks. Publish only a closed, complete file.
  let nrow = if cols.len == 0: 0 else: cols[0].len
  if colNames.len > 0 and colNames.len != cols.len:
    raise newException(ValueError, path & ": column name count mismatch")
  for col in cols:
    if col.len != nrow:
      raise newException(ValueError, path & ": column length mismatch")
    when T is string:
      for cell in col:
        if '\t' in cell or '\n' in cell or '\r' in cell:
          raise newException(ValueError, path & ": tab or line break in TSV field")
  let dir = path.parentDir
  if dir.len > 0: createDir dir
  let (f, tmp) = createTempFile(path.extractFilename & ".", ".tmp",
                               if dir.len > 0: dir else: ".")
  try:
    block:
      defer: f.close()
      for (k, v) in header: f.writeLine("# " & k & "=" & v)
      if colNames.len > 0: f.writeLine("# " & colKey & "=" & colNames.join(" "))
      for i in 0..<nrow:
        var line = ""
        for j in 0..<cols.len:
          if j > 0: line.add '\t'
          when T is string: line.add cols[j][i]
          else: line.add &"{cols[j][i]:.17g}"
        f.writeLine line
      f.flushFile()
    moveFile(tmp, path)
  finally:
    if fileExists(tmp): removeFile tmp

proc readTsv*(path: string): tuple[meta: Table[string, string], names: seq[string],
                                   cols: seq[seq[float]]] =
  ## Read numeric TSV fields. The column count is set by the first data line.
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
