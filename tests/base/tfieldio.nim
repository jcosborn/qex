import testutils
import qex, os

# TODO array of field?

qexInit()

threads: echo "thread ",threadNum," / ",numThreads
const
  fn = "tmpfield.lime"
  filemd = "test filemd"
  recordmd = "test recordmd"
var
  lat = latticeFromLocalLattice([8,8,8,8], nRanks)
  (l,g,_) = setupLattice(lat)
  rng = l.newRNGField(RngMilc6, 987654321)
  fr = l.Real
  fc = l.Complex
  fv = l.ColorVector
  fm = l.ColorMatrix
  fh = l.HalfFermion
  fd = l.DiracFermion
  c = getDefaultComm()

proc writer(): auto =
  result = l.newWriter(fn, filemd)
  check(result.status==0)

proc finish(w: var Writer) =
  w.close
  check(w.status==0)
  #c.barrier

proc save(w: var Writer, f: var Field) =
  f.gaussian rng
  w.write(f, recordmd)
  check(w.status==0)

proc save(f: var Field) =
  var w = writer()
  w.save f
  w.finish

proc reader(): auto =
  result = l.newReader(fn)
  check(result.status==0)
  check(result.fileMetadata==filemd)

proc finish(r: var Reader) =
  r.close
  check(r.status==0)
  #c.barrier

proc load(r: var Reader, f: Field) =
  var f2 = f.newOneOf
  r.read(f2)
  check(r.status==0)
  check(r.recordMetadata==recordmd)
  f2 -= f
  let n2 = f2.norm2
  check(n2==0)

proc load(f: Field) =
  var r = reader()
  r.load f
  r.finish

suite "Test field IO":

  test "save real":
    save fr
  test "load real":
    load fr

  test "save complex":
    save fc
  test "load complex":
    load fc

  test "save color vector":
    save fv
  test "load color vector":
    load fv

  test "save color matrix":
    save fm
  test "load color matrix":
    load fm

  test "save half fermion":
    save fh
  test "load half fermion":
    load fh

  test "save dirac fermion":
    save fd
  test "load dirac fermion":
    load fd

  test "save all":
    var w = writer()
    w.save fr
    w.save fc
    w.save fv
    w.save fm
    w.save fh
    w.save fd
    w.finish

  test "load all":
    var r = reader()
    r.load fr
    r.load fc
    r.load fv
    r.load fm
    r.load fh
    r.load fd
    r.finish

removeFile fn
qexFinalize()

