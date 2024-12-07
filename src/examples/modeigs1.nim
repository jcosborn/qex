import qex
import io / [parallelIo, timesliceIo, modfile]
import xmltree, strutils, times

template TAG(x: untyped): untyped = xmltree.`<>`(x)
template TXT(x: untyped): untyped = xmltree.newText(x)

proc makeTestEigFile*(fn: string, f: Field, n: int, tsio: TimesliceIo) =
  let lattSize = f.l.physGeom.join(" ")
  let rankGeom = f.l.rankGeom.join(" ")
  let localGeom = f.l.localGeom.join(" ")
  let numVecs = $n
  #let runDate = "05 Jul 13 16:45:15 EDT"
  let runDate = $now()
  let totalVolume = $f.l.physVol
  var ud =
    TAG MODMetaData(
      TAG id(TXT "eigenVecsTimeSlice"),
      TAG lattSize(TXT lattSize),
      TAG decay_dir(TXT "3"),
      TAG num_vecs(TXT numVecs),
      TAG ProgramInfo(
        TAG code_version(
          TAG basePrecision(TXT "32")
        ),
        TAG run_date(TXT runDate),
        TAG Setgeom(
          TAG latt_size(TXT lattSize),
          TAG logical_size(TXT rankGeom),
          TAG subgrid_size(TXT localGeom),
          TAG total_volume(TXT totalVolume),
          TAG subgrid_volume(TXT totalVolume)
        )
      )
    )
  var weights = TAG Weights()
  let nt = f.l[^1]
  var ws = newSeq[float](nt)
  for i in 0..<n:
    for j in 0..<nt:
      ws[j] = float(i+j+1)
    let w = ws.join(" ")
    weights.add TAG elem(TXT w)
  ud.add weights
  #echo ud

  var mw = newModFileWriter(fn, $ud)
  var r = newRngField(f.l, RngMilc6)
  for i in 0..<n:
    f.gaussian(r)
    for t in 0..<nt:
      #echo "begin write: ", mw.w.pos
      mw.beginWrite(packKey(@[t,i]))
      tsio.write(mw.w, f, t)
      mw.endWrite()
      #echo "end write: ", mw.w.pos, "  cksum: ", mw.w.crc32
  mw.close()

