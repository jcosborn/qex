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
  steps = intParam("steps", 10)
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

# For FGYin11
const useFGYin11 = false
when useFGYin11:
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

let
  # H = mkLeapfrog(steps = steps, V = mdv, T = mdt)
  # H = mkSW92(steps = steps, V = mdv, T = mdt)
  # H = mkOmelyan2MN(steps = steps, V = mdv, T = mdt)
  # H = mkOmelyan4MN4FP(steps = steps, V = mdv, T = mdt)
  # H = mkOmelyan4MN5FV(steps = steps, V = mdv, T = mdt)
  H =
    when useFGYin11:
      mkFGYin11(steps = steps, V = mdv, T = mdt, Vfg = fgv, save = fgsave(), load = fgload())
    else:
      mkOmelyan4MN5FV(steps = steps, V = mdv, T = mdt)

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

  #echo g.plaq
  #echo g.plaq2
  let pl = g.plaq3
  echo "plaq: ",pl.re," ",pl.im
  echo "topo: ",g.topo2DU1

echoTimers()
checkgc()
qexfinalize()
