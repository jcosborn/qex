import testutils
import base
import maths/complexNumbers
import maths/matrixConcept
import maths/matrixFunctions
import maths/types
import algorithms/numdiff
import math

# expDeriv(m, c) is the adjoint (reverse-mode VJP) derivative of exp: D exp(m†)[c]
# It satisfies
#   redot(c, D exp(m)[e]) == redot(expDeriv(m, c), e),
# with redot the real inner product Re tr(A† B).

proc fillAH(m: var auto, s: float) =
  ## anti-Hermitian (m† = -m): the stout case m = αF
  const nc = m.nrows
  m := 0
  for i in 0..<nc:
    m[i,i].im := 0.2*cos(s + 0.5*i.float)
    for j in (i+1)..<nc:
      let re = 0.2*sin(s + 1.7*i.float - 0.9*j.float)
      let im = 0.2*cos(s + 0.5*i.float + 1.3*j.float)
      m[i,j].re :=  re;  m[i,j].im := im
      m[j,i].re := -re;  m[j,i].im := im

proc fillGen(m: var auto, s: float) =
  const nc = m.nrows
  for i in 0..<nc:
    for j in 0..<nc:
      m[i,j].re := 0.3*sin(s + 1.1*i.float - 0.7*j.float)
      m[i,j].im := 0.3*cos(s + 0.6*i.float + 1.2*j.float)

proc runExpDerivTests(nc: static int) =
  type M = MatrixArray[nc, nc, ComplexType[float]]
  var mAH, mGen, c, e: M
  fillAH(mAH, 1.0)
  fillGen(mGen, 4.0)
  fillGen(c, 2.0)
  fillGen(e, 3.0)
  proc fwd(m: M): float =          # redot(c, D exp(m)[e]) via ndiff
    var err: float
    ndiff(result, err, proc(t: float): float = redot(c, exp(m + t*e)), 0.0, 0.1)
  suite "expDeriv Nc=" & $nc:
    test "anti-Hermitian m":
      withCT(1e-9): check fwd(mAH) ~ redot(expDeriv(mAH, c), e)
    test "general m":
      withCT(1e-9): check fwd(mGen) ~ redot(expDeriv(mGen, c), e)

when isMainModule:
  runExpDerivTests(1)
  runExpDerivTests(3)
