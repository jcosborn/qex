import testutils
import qex, io/parallelIo, base/hyper
import os, sequtils

# Parallel write and read back under the GC stress and borrow check of
# tpiogc.nims.  The IO buffers and the write map's send lists reach the gathers
# only as raw pointers.  If one is not kept referenced, the check aborts;
# without it, a build whose compiler keeps no stack copy reads back other data.

# [0,0].re holds the site's global lex index
proc fill(f: auto) =
  let lo = f.l
  var x = newSeq[int32](lo.nDim)
  for i in 0..<lo.nSites:
    lo.coord(x, i)
    f[i] := 0
    f[i][0,0].re := float lexIndex(x, lo.physGeom)

qexInit()

suite "Parallel IO":
  echo "rank ", myRank, "/", nRanks
  let lo = newLayout(latticeFromLocalLattice(@[4,4,4,4], nRanks), 1)
  let fn = "tpiogc.bin"
  var f = lo.ColorMatrix1()
  var g = lo.ColorMatrix1()
  fill f
  for ioranks in deduplicate(@[@[0], @[nRanks-1], deduplicate(@[0, nRanks-1])]):
    test "write and read back, io ranks " & $ioranks:
      let wm = lo.setupWrite(lo.physGeom, newSeq[int](lo.nDim), ioranks)
      var pw = openCreate(fn)
      pw.write(f, wm)
      pw.close
      g := 0
      var pr = openRead(fn)
      pr.read(g, wm)
      pr.close
      check(diffNorm2(f, g) == 0.0)
  if myRank == 0: removeFile(fn)

qexFinalize()
