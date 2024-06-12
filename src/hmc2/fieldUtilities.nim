import qex
import base
import gauge

proc zeroOdd*(phi: auto) =
  threads: phi.odd := 0

proc zero*(u: auto) =
  threads:
    for mu in 0..<u.len:
      for s in u[mu]: u[mu][s] := 0

proc pdagp*(p: auto): float = 
  var p2: float
  threads:
    var p2t = 0.0
    for mu in 0..<p.len: p2t += p[mu].norm2
    threadBarrier()
    threadMaster: p2 = p2t
  result = p2

proc reunit*(u: auto) =
  threads: u.projectSU

proc updateU*(u: auto; p: auto; c: float) =
  threads:
    for mu in 0..<p.len:
      for s in p[mu]:
        u[mu][s] := exp(c*p[mu][s])*u[mu][s]

proc setU*(u: auto; v: auto) =
  threads:
    for mu in 0..<v.len:
      for s in v[mu]: u[mu][s] := v[mu][s]

proc subtractForce*(p: auto; f: auto) =
  threads:
    for mu in 0..<f.len:
      for s in f[mu]: p[mu][s] -= f[mu][s]

if isMainModule:
  qexinit()

  var
    lat = intSeqParam("lat", @[4, 4, 4, 4])
    lo = lat.newLayout(@[1, 1, 1, 1])

    u1 = lo.newGauge
    u2 = lo.newGauge

  qexfinalize()