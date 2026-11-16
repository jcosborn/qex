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
  (l,_,_) = setupLattice(lat)
  rng = l.newRNGField(RngMilc6, 987654321)
  fr = l.Real
  fc = l.Complex
  fv = l.ColorVector
  fm = l.ColorMatrix
  fh = l.HalfFermion
  fd = l.DiracFermion
  #c = getDefaultComm()

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
    subtest save(fr)
  test "load real":
    subtest load(fr)

  test "save complex":
    subtest save(fc)
  test "load complex":
    subtest load(fc)

  test "save color vector":
    subtest save(fv)
  test "load color vector":
    subtest load(fv)

  test "save color matrix":
    subtest save(fm)
  test "load color matrix":
    subtest load(fm)

  test "save half fermion":
    subtest save(fh)
  test "load half fermion":
    subtest load(fh)

  test "save dirac fermion":
    subtest save(fd)
  test "load dirac fermion":
    subtest load(fd)

  test "save all":
    var w: typeof(writer())
    subtest:
      w = writer()
    subtest w.save(fr)
    subtest w.save(fc)
    subtest w.save(fv)
    subtest w.save(fm)
    subtest w.save(fh)
    subtest w.save(fd)
    subtest w.finish

  test "load all":
    var r: typeof(reader())
    subtest:
      r = reader()
    subtest r.load(fr)
    subtest r.load(fc)
    subtest r.load(fv)
    subtest r.load(fm)
    subtest r.load(fh)
    subtest r.load(fd)
    subtest r.finish

  when defined(QioReader) or not defined(QioliteReader):
    # QIO may leave peers inside a broadcast after a header error.
    if nRanks == 1:
      test "reject empty file on open":
        writeFile(fn, "")
        let r = l.newReader(fn)
        check(r.status != 0)

removeFile fn
qexFinalize()

