import qex, nn
import std/[math, unittest]
export unittest

proc fields*[T: SomeFloat](lo: Layout[VLEN]; channels: int; value: T = T(0)): seq[RealField[T]] =
  result = newSeq[RealField[T]](channels)
  for c in 0..<channels: result[c] = realField(lo,T)
  let fs = result
  threads:
    for f in fs: f := value

proc sample*[F](f: F; site: int): auto =
  var r: numberType(F)
  r := f{site}
  r

proc fill*[F; T: SomeFloat](fs: seq[F]; value: T) =
  threads:
    for f in fs: f := value

proc dot*[F](a, b: seq[F]): float64 =
  ## Global sum of a_c(s) b_c(s), accumulated in double precision.
  var acc: DLatticeRealV.T
  acc := 0
  for c in 0..<a.len:
    for e in a[c]:
      var x, y: DLatticeRealV.T
      x := a[c][e]
      y := b[c][e]
      acc += x*y
  result = simdSum(acc)
  a[0].l.comm.rankSum(result)

proc close*[T: SomeFloat](x, y: T): bool =
  let tol = when T is float32: T(2e-6) else: T(2e-12)
  abs(x-y) <= tol*max(T(1),abs(y))

proc sameFields*[F](a, b: seq[F]): bool =
  ## The RMS difference is within the precision's tolerance of the RMS of b (at least one).
  if a.len != b.len: return false
  type T = numberType(F)
  let tol = when T is float32: 2e-6 else: 2e-12
  let n = float(a[0].l.physVol)
  for c in 0..<a.len:
    let d = sqrt((a[c]-b[c]).norm2/n)
    if d > tol*max(1.0, sqrt(b[c].norm2/n)): return false
  true
