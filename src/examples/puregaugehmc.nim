import qex
import gauge, physics/qcdTypes
import mdevolve

qexinit()

let
  lat = @[8,8,8,8]
  #lat = @[8,8,8]
  #lat = @[32,32]
  #lat = @[1024,1024]
  lo = lat.newLayout
  gc = GaugeActionCoeffs(plaq:6)
var r = lo.newRNGField RngMilc6
var R:RngMilc6  # global RNG
R.seed(7^7, 0)

var g = lo.newgauge
g.random r

echo g.plaq
echo g.gaugeAction2 gc

var
  p = lo.newgauge
  f = lo.newgauge
  g0 = lo.newgauge

let
  tau = 1.0
  steps = 20
  trajs = 10

proc mdt(t: float) =
  threads:
    for mu in 0..<g.len:
      for s in g[mu]:
        g[mu][s] := exp(t*p[mu][s])*g[mu][s]
proc mdv(t: float) =
  f.gaugeforce2(g, gc)
  threads:
    for mu in 0..<f.len:
      p[mu] -= t*f[mu]

# For FGYin11
proc fgv(t: float) =
  f.gaugeforce2(g, gc)
  threads:
    for mu in 0..<g.len:
      for s in g[mu]:
        g[mu][s] := exp((-t)*f[mu][s])*g[mu][s]
var gg = lo.newgauge
proc fgsave =
  threads:
    for mu in 0..<g.len:
      gg[mu] := g[mu]
proc fgload =
  threads:
    for mu in 0..<g.len:
      g[mu] := gg[mu]

let
  # H = mkLeapfrog(steps = steps, V = mdv, T = mdt)
  # H = mkSW92(steps = steps, V = mdv, T = mdt)
  # H = mkOmelyan2MN(steps = steps, V = mdv, T = mdt)
  # H = mkOmelyan4MN4FP(steps = steps, V = mdv, T = mdt)
  # H = mkOmelyan4MN5FV(steps = steps, V = mdv, T = mdt)
  H = mkFGYin11(steps = steps, V = mdv, T = mdt, Vfg = fgv, save = fgsave(), load = fgload())

for n in 1..trajs:
  var p2 = 0.0
  threads:
    p.randomTAH r
    var p2t = 0.0
    for i in 0..<p.len:
      p2t += p[i].norm2
      g0[i] := g[i]
    threadMaster: p2 = p2t
  let
    ga0 = g0.gaugeAction2 gc
    t0 = 0.5*p2
    h0 = ga0 + t0
  echo "Begin H: ",h0,"  Sg: ",ga0,"  T: ",t0

  H.evolve tau

  threads:
    var p2t = 0.0
    for i in 0..<p.len:
      p2t += p[i].norm2
    threadMaster: p2 = p2t
  let
    ga1 = g.gaugeAction2 gc
    t1 = 0.5*p2
    h1 = ga1 + t1
  echo "End H: ",h1,"  Sg: ",ga1,"  T: ",t1

  #when true:
  when false:
    block:
      var g1 = lo.newgauge
      var p1 = lo.newgauge
      threads:
        for i in 0..<g1.len:
          g1[i] := g[i]
          p1[i] := p[i]
          p[i] := -1*p[i]
      H.evolve tau
      threads:
        var p2t = 0.0
        for i in 0..<p.len:
          p2t += p[i].norm2
        threadMaster: p2 = p2t
      let
        ga1 = g.gaugeAction2 gc
        t1 = 0.5*p2
        h1 = ga1 + t1
      echo "Reversed H: ",h1,"  Sg: ",ga1,"  T: ",t1
      echo "Reversibility: dH: ",h1-h0,"  dSg: ",ga1-ga0,"  dT: ",t1-t0
      #echo p[0][0]
      for i in 0..<g1.len:
        g[i] := g1[i]
        p[i] := p1[i]

  let
    dH = h1 - h0
    acc = exp(-dH)
    accr = R.uniform
  if accr <= acc:  # accept
    echo "ACCEPT:  dH: ",dH,"  exp(-dH): ",acc,"  r: ",accr
  else:  # reject
    echo "REJECT:  dH: ",dH,"  exp(-dH): ",acc,"  r: ",accr
    threads:
      for i in 0..<g.len:
        g[i] := g0[i]

  echo g.plaq

echoTimers()
qexfinalize()
