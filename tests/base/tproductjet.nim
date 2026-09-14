import qex
import testutils

type M = MatrixArray[2,2,float64]

proc checkProduct(n, order: static int) =
  var a: array[n,array[order+1,M]]
  for k in 0..<n:
    for d in 0..order:
      for i in 0..<2:
        for j in 0..<2:
          a[k][d][i,j] = float((k+1)*(3*i-j+2)+(d+1)*(i+2*j-1))/13.0
  var got, want: M
  template factor(x, k: untyped) =
    forStatic d, (if order == n: 1 else: 0), order:
      x[d] := a[k][d]
  productJet(got, n, order, factor)

  # Independent expansion: assign each distinct seed to a different factor.
  var pos: array[order,int]
  proc expand(d: int) =
    if d == order:
      var p: M
      for k in 0..<n:
        var slot = 0
        for j in 0..<order:
          if pos[j] == k: slot = j+1
        if k == 0: p := a[k][slot]
        else:
          var t: M
          t := p*a[k][slot]
          p := t
      want += p
    else:
      for k in 0..<n:
        var used = false
        for j in 0..<d:
          if pos[j] == k: used = true
        if not used:
          pos[d] = k
          expand(d+1)
  expand(0)
  check norm2(got-want) < 1e-24*(1.0+norm2(want))

# Count actual element multiplications to guard the stencil cost contract.
type Counted = object
  value: float
var products, loads: int
template eval(t: typedesc[Counted]): typedesc = Counted
proc `:=`(r: var Counted, x: Counted) = r = x
proc `:=`(r: var Counted, x: SomeNumber) = r.value = float(x)
proc `+=`(r: var Counted, x: Counted) = r.value += x.value
proc `*`(x, y: Counted): Counted =
  inc products
  Counted(value: x.value*y.value)

proc checkCount(order: static int) =
  products = 0
  loads = 0
  var r: Counted
  template factor(a, k: untyped) =
    inc loads
    forStatic d, (if order == 3: 1 else: 0), order: a[d] := 1
  productJet(r, 3, order, factor)
  check products == [2,5,9,9,0][order]
  check r.value == [1.0,3.0,6.0,6.0,0.0][order]
  check loads == (if order > 3: 0 else: 3)

qexInit()
suite "Affine matrix product jets":
  test "noncommuting products and mixed derivatives":
    forStatic order, 0, 1: checkProduct(1, order)
    forStatic order, 0, 4: checkProduct(3, order)
    forStatic order, 0, 6: checkProduct(5, order)
  test "three-factor jets keep their multiplication counts":
    forStatic order, 0, 4: checkCount(order)
qexFinalize()
