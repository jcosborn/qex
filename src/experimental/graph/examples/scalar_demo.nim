import math
import std/assertions

import ../[core, scalar]

let grt = initGraphRuntime()
grt.graphDebug = true

let x = toGvalue(grt, 0.0)
let y = toGvalue(grt, 0.0)
let w = x - 2.0
let v = w + y
let z = v * (-v) / w
let dzdy = z.grad y

func f(a, b: float): float = (a + b - 2.0) * (2.0 - a - b) / (a - 2.0)
func dfdb(a, b: float): float = -2.0 * (a + b - 2.0) / (a - 2.0)

let a = 1.1
let b = 3.7
let c = 1.3

x.update a
y.update b
echo z.treeRepr
echo dzdy.treeRepr
z.eval
dzdy.eval
echo "z = ", z
echo z.treeRepr
echo "dzdy = ", dzdy
echo dzdy.treeRepr

dumpGradientList(grt)

doAssert almostEqual(z.getfloat, f(a, b))
doAssert almostEqual(dzdy.getfloat, dfdb(a, b))

y.update c
z.eval
dzdy.eval
echo "z = ", z
echo z.treeRepr
echo "dzdy = ", dzdy
echo dzdy.treeRepr
doAssert almostEqual(z.getfloat, f(a, c))
doAssert almostEqual(dzdy.getfloat, dfdb(a, c))
