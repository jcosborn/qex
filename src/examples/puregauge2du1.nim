import qex
import gauge, physics/qcdTypes
import mdevolve
import os, strutils

template checkgc =
  let fl = instantiationInfo(fullPaths=true)
  echo "CHECK GC: ",fl.filename,":",fl.line
  GC_fullCollect()
  echo GC_getStatistics()
  #dumpNumberOfInstances()  # remember to -d:nimTypeNames

proc setup(defaultLat:openarray[int],
    fseed:uint64 = uint64(11^11), gseed:uint64 = uint64(7^7)):auto =
  var lat:seq[int]
  let pc = paramCount()
  if pc > 0 and paramStr(1).isDigit:
    lat = @[]
    for i in 1..pc:
      if not paramStr(i).isDigit: break
      lat.add paramStr(i).parseInt
  else:
    lat = @defaultLat
  let lo = lat.newLayout
  var
    g = lo.newGauge
    r = newRNGField(RngMilc6, lo, fseed)
    R:RngMilc6  # global RNG
  R.seed(gseed, 0)
  return (g, r, R)

proc topo2DU1(g:array or seq):float =
  tic()
  const nc = g[0][0].nrows
  let
    lo = g[0].l
    nd = lo.nDim
    t = newTransporters(g, g[0], 1)
  var p = 0.0
  toc("topo2DU1 setup")
  threads:
    tic()
    var tp:type(atan2(g[0][0][0,0].im, g[0][0][0,0].re))
    for mu in 1..<nd:
      for nu in 0..<mu:
        let tpl = (t[mu]^*g[nu]) * (t[nu]^*g[mu]).adj
        for i in tpl:
          tp += atan2(tpl[i][0,0].im, tpl[i][0,0].re)
    var v = tp.simdSum
    v.threadRankSum
    threadSingle: p += v
    toc("topo2DU1 work")
  toc("topo2DU1 threads")
  p/TAU

proc maxTreeFix(f:seq, val:float) =
  ## Set link to `val`, for those links on the maximal tree.
  let
    nd = f.len
    lo = f[0].l
    lat = lo.physGeom
  threads:
    var co = newseq[cint](nd)
    for i in 0..<nd:
      for j in lo.sites:
        lo.coord(co,(lo.myRank,j))
        var fix = co[i] < lat[i]-1
        for k in i+1 .. nd-1: fix = fix and co[k] == 0
        if fix: f[i]{j} := val

qexinit()
threads: echo "thread ",threadNum," / ",numThreads

let
  beta = floatParam("beta", 5.0)
  trajs = intParam("trajs", 512)
  tau = floatParam("tau", 2.0)
  steps = intParam("steps", 5)
  fseed = intParam("fseed", 11^11).uint64
  gseed = intParam("gseed", 7^7).uint64
  gfix = intParam("gfix", 0).bool
  gfixunit = intParam("gfixunit", 1).bool

var (g,r,R) = setup([64,64],fseed,gseed)

let
  lo = g[0].l
  lat = lo.physGeom
  nd = lat.len
  gc = GaugeActionCoeffs(plaq:beta)

echo "latsize = ",lat
echo "volume = ",lo.physVol
echo "beta = ",beta
echo "trajs = ",trajs
echo "tau = ",tau
echo "steps = ",steps
echo "fseed = ",fseed
echo "gseed = ",gseed
echo "gfix = ",gfix.int
echo "gfixunit = ",gfixunit.int

g.random r
if gfix and gfixunit: g.maxTreeFix 1.0

#echo g.plaq
#echo g.plaq2
echo "Initial plaq: ",g.plaq3
#echo g.gaugeAction2 gc

var
  p = lo.newgauge
  f = lo.newgauge
  g0 = lo.newgauge

proc mdt(t:float) =
  threads:
    for i in 0..<g.len:
      for e in g[i]:
        let etpg = exp(t*p[i][e])*g[i][e]
        g[i][e] := etpg
proc mdv(t:float) =
  f.gaugeforce2(g, gc)
  if gfix: f.maxTreeFix 0.0
  threads:
    for i in 0..<f.len:
      for e in f[i]:
        let tf = (-t)*f[i][e]
        p[i][e] += tf

# For force gradient update
const useFG = true
const useApproxFG2 = false
when useFG:
  proc fgv(t:float) =
    f.gaugeforce2(g, gc)
    if gfix: f.maxTreeFix 0.0
    threads:
      for i in 0..<g.len:
        for e in g[i]:
          let etfg = exp((-t)*f[i][e])*g[i][e]
          g[i][e] := etfg
  var gg = lo.newgauge
  proc fgsave =
    threads:
      for i in 0..<g.len:
        gg[i] := g[i]
  proc fgload =
    threads:
      for i in 0..<g.len:
        g[i] := gg[i]
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
      mkOmelyan2MN(steps = steps, V = V, T = T)
      # mkOmelyan4MN5FP(steps = steps, V = V, T = T)
      # mkOmelyan4MN5FV(steps = steps, V = V, T = T)
      # mkOmelyan6MN7FV(steps = steps, V = V, T = T)

for n in 1..trajs:
  echo "Begin traj: ",n
  var p2 = 0.0
  threads:
    p.randomTAH r
  if gfix: p.maxTreeFix 0.0
  threads:
    var p2t = 0.0
    for i in 0..<p.len:
      p2t += p[i].norm2
      g0[i] := g[i]
    threadMaster: p2 = p2t
  let
    ga0 = g0.gaugeAction2 gc
    t0 = 0.5*p2
    h0 = ga0 + t0
  #echo "Begin H: ",h0,"  Sg: ",ga0,"  T: ",t0

  H.evolve tau
  H.finish

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
      H.finish
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

  #echo g.plaq
  #echo g.plaq2
  let pl = g.plaq3
  echo "plaq: ",pl.re," ",pl.im
  echo "topo: ",g.topo2DU1

echoTimers()
checkgc()
qexfinalize()
