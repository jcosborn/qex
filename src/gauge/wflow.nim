import base
import gaugeUtils, gaugeAction

#[

d/dt Vt = Z(Vt) Vt

Runge-Kutta:

W0 <- Vt
W1 <- exp(1/4 Z0) W0
W2 <- exp(8/9 Z1 - 17/36 Z0) W1
V(t+eps) <- exp(3/4 Z2 - 8/9 Z1 + 17/36 Z0) W2

where

Zi = eps Z(Wi)

]#

proc gaugeFlow*[T](g:openArray[T], tmax = 1.0, eps = 0.01) =
  ## Wilson flow.
  ## The input gauge field will be modified.
  var
    p = g[0].l.newGauge  # mom
    f = g[0].l.newGauge  # force
  threads:
    var t = 0.0
    for i in 0..<p.len: p[i] := 0  # FIX needed?
    while t < tmax:
      threadBarrier()
      f.gaugeForce2 g
      threadBarrier()
      for mu in 0..<f.len:
        echo t, " f[", mu, "]: ", f[mu].norm2
        p[mu] += (1.0/4.0)*eps*f[mu]
        for e in g[mu]:
          let t = exp(p[mu][e])*g[mu][e]
          g[mu][e] := t
      threadBarrier()
      f.gaugeForce2 g
      threadBarrier()
      for mu in 0..<f.len:
        echo t, " f[", mu, "]: ", f[mu].norm2
        p[mu] := (8.0/9.0)*eps*f[mu] + (-17.0/9.0)*p[mu]
        for e in g[mu]:
          let t = exp(p[mu][e])*g[mu][e]
          g[mu][e] := t
      threadBarrier()
      f.gaugeForce2 g
      threadBarrier()
      for mu in 0..<f.len:
        echo t, " f[", mu, "]: ", f[mu].norm2
        p[mu] := (3.0/4.0)*eps*f[mu] - p[mu]
        for e in g[mu]:
          let t = exp(p[mu][e])*g[mu][e]
          g[mu][e] := t
      echo t, " act: ", g.gaugeAction2

when isMainModule:
  import qex, gauge, physics/qcdTypes

  proc printPlaq(g: any) =
    let
      p = g.plaq
      sp = 2.0*(p[0]+p[1]+p[2])
      tp = 2.0*(p[3]+p[4]+p[5])
    echo "plaq ",p
    echo "plaq ss: ",sp," st: ",tp," tot: ",p.sum

  qexInit()

  var
    lat = [8,8,8,8]
    lo = lat.newLayout
    g = lo.newGauge
    r = newRNGField(RngMilc6, lo)
  threads: g.random r
  g.printPlaq

  g.gaugeFlow

  qexFinalize()
