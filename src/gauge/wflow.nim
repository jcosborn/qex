import base, maths
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

template gaugeFlow*(g: array|seq, steps: int, eps: float, measure: untyped): untyped =
  ## Wilson flow.
  ## The input gauge field will be modified.
  ## `wflowT` is injected for `measure`.
  proc flowProc {.gensym.} =
    tic("flowProc")
    const nc = g[0][0].nrows.float
    var
      p = g[0].l.newGauge  # mom
      f = g[0].l.newGauge  # force
      n = 1
    while true:
      let t = n * eps
      let epsnc = eps * nc  # compensate force normalization
      f.gaugeForce g
      threads:
        for mu in 0..<f.len:
          for e in g[mu]:
            var v {.noinit.}:type(load1(f[0][0]))
            v := (-1.0/4.0)*epsnc*f[mu][e]
            let t = exp(v)*g[mu][e]
            p[mu][e] := v
            g[mu][e] := t
      f.gaugeForce g
      threads:
        for mu in 0..<f.len:
          for e in g[mu]:
            var v {.noinit.}:type(load1(f[0][0]))
            v := (-8.0/9.0)*epsnc*f[mu][e] + (-17.0/9.0)*p[mu][e]
            let t = exp(v)*g[mu][e]
            p[mu][e] := v
            g[mu][e] := t
      f.gaugeForce g
      threads:
        for mu in 0..<f.len:
          for e in g[mu]:
            var v {.noinit.}:type(load1(f[0][0]))
            v := (-3.0/4.0)*epsnc*f[mu][e] - p[mu][e]
            let t = exp(v)*g[mu][e]
            g[mu][e] := t
      let wflowT {.inject,used.} = t
      measure
      inc n
      if steps>0 and n>steps:
        break
    toc("end")
  flowProc()
template gaugeFlow*(g: array|seq, eps: float, measure: untyped): untyped =
  gaugeFlow(g, 0, eps):
    measure

when isMainModule:
  import qex, gauge, physics/qcdTypes
  import os, sequtils

  proc printPlaq(g: auto) =
    let
      p = g.plaq
      sp = 2.0*(p[0]+p[1]+p[2])
      tp = 2.0*(p[3]+p[4]+p[5])
    echo "plaq ",p
    echo "plaq ss: ",sp," st: ",tp," tot: ",p.sum

  qexInit()

  let
    fn = if paramCount() > 0: paramStr 1 else: ""
    lat = if fn.len == 0: @[8,8,8,8] else: fn.getFileLattice
    lo = lat.newLayout
  var g = lo.newGauge
  if fn.len == 0:
    g.random
  elif 0 != g.loadGauge fn:
    echo "ERROR: couldn't load gauge file: ",fn
    qexFinalize()
    quit(-1)
  g.printPlaq

  g.gaugeFlow(6, 0.01):
    echo "WFLOW ",wflowT
    g.printPlaq

  when g[0][0].nrows == 1 or g[0][0].nrows == 3:
    if fn.len == 0:
      when g[0][0].nrows == 1:
        when defined(FUELCompat):
          const
            p0 = @[
              0.01554978485185948,
              0.01864915539989097,
              0.02185632226074754,
              0.01631465138128996,
              0.01599186401283072,
              0.02143528647233249]
        else:
          const
            p0 = @[
              0.02071819961326493,
              0.01982444094659434,
              0.01998727952235725,
              0.01837579824880835,
              0.01528029895697336,
              0.01782410235480616]
      when g[0][0].nrows == 3:
        when defined(FUELCompat):
          const
            p0 = @[
              0.01984173412117658,
              0.01775826434617258,
              0.01814446406048498,
              0.01960829801483532,
              0.01913142046361973,
              0.01953858380137068]
        else:
          const
            p0 = @[
              0.01960725848281519,
              0.01982378149813489,
              0.01938877647467847,
              0.0185899778070918,
              0.0180821938831715,
              0.01876842496122964]
      let p = g.plaq
      let d = zip(p,p0).foldl(a + abs(b[0]-b[1]), 0.0) / sum(p0)
      if d > 2e-14:
        echo "Test failed, relative diff: ", d
        echo "Expected:\t",p0
        echo "Actual:\t",p
        qexExit -1

  qexFinalize()
