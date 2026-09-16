import qex, base/alignedMem
import testutils

qexInit()
letParam:
  lat = latticeFromLocalLattice(@[4,4,4,4], nRanks)
let lo = lat.newLayout
let x = lo.ColorMatrix(3)
let u = lo.ColorMatrix(3)
let v = lo.ColorMatrix(3)
let y = lo.ColorMatrix(3)
let z = lo.ColorMatrix(3)
var rng = lo.newRNGField(Philox4x64, 940260914'u64)
threads:
  x.gaussian rng
  u.gaussian rng
  v.gaussian rng

suite "Field storage":
  test "shape descriptors preserve data owned by existing aliases":
    var f = x
    let old = f
    let raw = getRawMemAllocated()
    f = f.newShape
    let typed = lo.newShape(DColorMatrixV)
    check getRawMemAllocated() == raw
    check f.l == lo and typed.l == lo
    check f.s.data == nil and typed.s.data == nil
    check f.s.bytes == x.s.bytes and typed.s.bytes == x.s.bytes
    check old.s.data == x.s.data
    check diffNorm2(old, x) == 0.0

  test "a shape from a field array describes one field":
    let fa = newFieldArray(lo, type(x), 3)
    threads:
      fa[1] := x
    let f = fa[1]
    let raw = getRawMemAllocated()
    let sh = f.newShape
    check getRawMemAllocated() == raw
    check sh.s.len == lo.nSitesOuter
    check sh.s.bytes == x.s.bytes
    check sh.s.data == nil
    check f.s.data == fa[1].s.data
    check diffNorm2(f, x) == 0.0

  test "shifters reuse supplied outputs after input and destination changes":
    for mu in 0..<lo.nDim:
      for sgn in [-1,1]:
        let refsh = newShifter(x, mu, sgn)
        let raw = getRawMemAllocated()
        var sh = newShifter(x, mu, sgn, dest=y)
        check getRawMemAllocated() == raw
        check sh.field.s.data == y.s.data
        for pass in 0..2:
          let src = if pass == 1: u else: x
          sh.field = if pass == 1: z else: y
          threads:
            discard refsh ^* src
            discard sh ^* src
          check diffNorm2(sh.field, refsh.field) == 0.0
        check sh.field.s.data == y.s.data

  test "transporters reuse outputs with current links and both signs":
    for mu in 0..<lo.nDim:
      for sgn in [-1,1]:
        var reftr = newTransporter(u, x, mu, sgn)
        let raw = getRawMemAllocated()
        var tr = newTransporter(u, x, mu, sgn, dest=y)
        check getRawMemAllocated() == raw
        check tr.field.s.data == y.s.data
        for pass in 0..2:
          let link = if pass == 1: v else: u
          let src = if pass == 2: u else: x
          reftr.setLink(link)
          tr.setLink(link)
          tr.field = if pass == 1: z else: y
          threads:
            discard reftr ^* src
            discard tr ^* src
          check diffNorm2(tr.field, reftr.field) == 0.0
        check tr.field.s.data == y.s.data

  test "supplied outputs require the operand layout":
    var lat = newSeq[int](lo.nDim)
    for i, n in lo.physGeom: lat[i] = int(n)
    lat[0] *= 2
    let wrong = lat.newLayout.ColorMatrix(3)
    expect(ValueError):
      discard newShifter(x, 0, 1, dest=wrong)
    expect(ValueError):
      discard newTransporter(u, x, 0, 1, dest=wrong)

qexFinalize()
