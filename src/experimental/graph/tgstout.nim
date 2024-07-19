import qex, algorithms/numdiff, gauge/stoutsmear
import core, scalar, gauge

qexInit()

letParam:
  lat = @[12,12,12,24]
  dt = 0.1
  eps = 0.004
  nstep = 3
  beta = 5.4
  seed:uint = 1234567891

echoParams()
echo "rank ", myRank, "/", nRanks
threads: echo "thread ", threadNum, "/", numThreads

let
  lo = lat.newLayout
  vol = lo.physVol
  gc = GaugeActionCoeffs(plaq: beta)
  g = lo.newgauge
  u = lo.newgauge

var
  r = lo.newRNGField(RngMilc6, seed)
  ss = lo.newStoutSmear(dt)

for i in 0..3:
  ss.smear(g, g)

g.random r
g.echoPlaq
for i in 1..nstep:
  if i==1:
    ss.smear(g, u)
  else:
    ss.smear(u, u)
u.echoPlaq
let sgs = gc.gaugeAction1 u
echo "smear S: ",sgs

proc act(t: float): float =
  var ss = lo.newStoutSmear(t)
  for i in 1..nstep:
    if i==1:
      ss.smear(g, u)
    else:
      ss.smear(u, u)
  gc.gaugeAction1 u

var ndt, err: float
ndiff(ndt, err, act, dt, eps, ordMax=3)
echo "numdiff smear dS/dt: ",ndt," +/- ",err

proc stout(g, t: Gvalue, n: int): Gvalue =
  var g = g
  for i in 1..n:
    g = axexpmuly(t, gaugeForce(actWilson(-3.0), g), g)
  g

let
  gdt = toGvalue dt
  gg = toGvalue g
  gs = gg.stout(gdt, nstep)
  s = gc.gaugeAction gs
  ddt = s.grad gdt

# echo ddt.treeRepr

gs.eval.getgauge.echoPlaq
let sgg = s.eval.getfloat
echo "graph S: ",sgg

let gddt = ddt.eval.getfloat
echo "graph dS/dt: ",gddt
# echo ddt.treeRepr

proc gact(t: float): float =
  gdt.update t
  s.eval.getfloat
var gndt, gerr: float
ndiff(gndt, gerr, gact, dt, eps, ordMax=4)
echo "numdiff graph dS/dt: ",gndt," +/- ",gerr

let
  rds = abs((sgs-sgg)/(sgs+sgg))
  rgdt = abs((ndt-gddt)/(ndt+gddt))
  rndt = abs((ndt-gndt)/(ndt+gndt))
echo "rel dS: ",rds
echo "rel graph dS/dt: ",rgdt
echo "rel ndiff dS/dt: ",rndt

doassert rds < 1e-11
doassert rgdt < 1e-11
doassert rndt < 1e-11

qexFinalize()
