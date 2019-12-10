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
  beta = floatParam("beta", 6.0)
  adjFac = floatParam("adjFac", -0.25)
  gc = GaugeActionCoeffs(plaq: beta, adjplaq: beta*adjFac)
var r = lo.newRNGField(RngMilc6, 987654321)
var R:RngMilc6  # global RNG
R.seed(987654321, 987654321)

var g = lo.newgauge
#g.random r
g.unit

echo g.plaq
echo g.gaugeAction2 gc
echo gc.actionA g

var
  p = lo.newgauge
  f = lo.newgauge
  g0 = lo.newgauge

let
  tau = 2.0
  steps = 4
  trajs = 10

proc mdt(t: float) =
  threads:
    for mu in 0..<g.len:
      for s in g[mu]:
        g[mu][s] := exp(t*p[mu][s])*g[mu][s]
proc mdv(t: float) =
  #f.gaugeforce2(g, gc)
  gc.forceA(g, f)
  threads:
    for mu in 0..<f.len:
      p[mu] -= t*f[mu]

# For force gradient update
const useFG = true
const useApproxFG2 = false
proc fgv(t: float) =
  #f.gaugeforce2(g, gc)
  gc.forceA(g, f)
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
proc updatefga(ts,gs:openarray[float]) =
  let
    t = ts[0]
    g = gs[0]
  if g != 0:
    if t != 0:
      # Approximate the force gradient update with a Taylor expansion.
      let (tf,tg) = approximateFGcoeff(t,g)
      fgsave()
      fgv tg
      mdv tf
      fgload()
    else:
      quit("Force gradient without the force update.")
  elif t != 0:
    mdv t
  else:
    quit("No updates required.")
proc updatefga2(ts,gs:openarray[float]) =
  let
    t = ts[0]
    g = gs[0]
  if g != 0:
    if t != 0:
      # Approximate the force gradient update with two Taylor expansions.
      let (tf,tg) = approximateFGcoeff2(t,g)
      fgsave()
      fgv tg[0]
      mdv tf[0]
      fgload()
      fgv tg[1]
      mdv tf[1]
      fgload()
    else:
      quit("Force gradient without the force update.")
  elif t != 0:
    mdv t
  else:
    quit("No updates required.")

let
  # Omelyan's triple star integrators, see Omelyan et. al. (2003)
  H =
    when useFG:
      let
        (VG,T) =
          if useApproxFG2: newIntegratorPair(updatefga2, mdt)
          else: newIntegratorPair(updatefga, mdt)
        V = VG[0]
      # mkOmelyan4MN4F2GVG(steps = steps, V = V, T = T)
      # mkOmelyan4MN4F2GV(steps = steps, V = V, T = T)
      # mkOmelyan4MN5F1GV(steps = steps, V = V, T = T)
      # mkOmelyan4MN5F1GP(steps = steps, V = V, T = T)
      # mkOmelyan4MN5F2GV(steps = steps, V = V, T = T)
      mkOmelyan4MN5F2GP(steps = steps, V = V, T = T)
      # mkOmelyan6MN5F3GP(steps = steps, V = V, T = T)
    else:
      let (V,T) = newIntegratorPair(mdv, mdt)
      # mkOmelyan2MN(steps = steps, V = V, T = T)
      # mkOmelyan4MN5FP(steps = steps, V = V, T = T)
      # mkOmelyan4MN5FV(steps = steps, V = V, T = T)
      mkOmelyan6MN7FV(steps = steps, V = V, T = T)

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
    #ga0 = g0.gaugeAction2 gc
    ga0 = gc.actionA g0
    t0 = 0.5*p2
    h0 = ga0 + t0
  echo "Begin H: ",h0,"  Sg: ",ga0,"  T: ",t0

  H.evolve tau
  H.finish

  threads:
    var p2t = 0.0
    for i in 0..<p.len:
      p2t += p[i].norm2
    threadMaster: p2 = p2t
  let
    #ga1 = g.gaugeAction2 gc
    ga1 = gc.actionA g
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
      H.finish
      threads:
        var p2t = 0.0
        for i in 0..<p.len:
          p2t += p[i].norm2
        threadMaster: p2 = p2t
      let
        ga1 = gc.actionA g
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
