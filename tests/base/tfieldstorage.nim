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

  test "array ownership survives aliases and preserves sparse shapes on release":
    let fa = newFieldArray2(lo, type(x), [3,3], mu != nu)
    let alias = fa
    let view = fa[0,1]
    let copied = fa.newOneOf
    check fa.mem != copied.mem
    for mu in 0..<3:
      for nu in 0..<3:
        if mu == nu:
          check fa[mu,nu] == nil and copied[mu,nu] == nil
        else:
          threads:
            fa[mu,nu] := x
            copied[mu,nu] := u
    let raw = getRawMemUsed()
    let bytes = fa.mem.size
    fa.release
    check getRawMemUsed() == raw-bytes
    check alias.mem.size == 0 and alias.mem.p == nil
    check view.s.data == nil and view.s.len == lo.nSitesOuter
    alias.release
    check getRawMemUsed() == raw-bytes
    let fresh = fa.newOneOf
    for mu in 0..<3:
      for nu in 0..<3:
        if mu != nu:
          check diffNorm2(copied[mu,nu], u) == 0
          threads: fresh[mu,nu] := x
          check diffNorm2(fresh[mu,nu], x) == 0
    copied.release
    fresh.release

  test "empty field arrays release their backing allocation":
    let fa = newFieldArray2(lo, type(x), [2,2], false)
    let copied = fa.newOneOf
    check copied.shape == fa.shape and copied.mem == nil
    copied.release
    let raw = getRawMemUsed()
    let bytes = fa.mem.size
    fa.release
    fa.release
    check getRawMemUsed() == raw-bytes

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
