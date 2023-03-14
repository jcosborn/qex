import qex
import testutils

qexInit()

suite "Test field IO":
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

  proc save(f: var auto) =
    f.gaussian rng
    var w = l.newWriter(fn, filemd)
    check(w.status==0)
    w.write(f, recordmd)
    check(w.status==0)

  proc load(f: auto) =
    var f2 = f.newOneOf
    var r = l.newReader(fn)
    check(r.status==0)
    r.read(f2)
    f2 -= f
    let n2 = f2.norm2
    check(n2==0)

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

qexFinalize()
