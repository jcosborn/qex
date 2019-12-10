#[
  Port of my_conn4d_meas.lua
]#

import qex, gauge/hypsmear, physics/stagSolve
import times, strformat

qexInit()

# Accept command line arguments and set up parameters
let
  # Required parameters
  inlat = strParam("inlat")   # Input gauge file name
  outfn = strParam("outfn", "scprop.out") # Save results to outfn
  mass = floatParam("mass", 0.1)   # The quark mass

  # optional parameters
  cg_prec = floatParam("cg_prec", 1e-9) # Max residual with default
  cg_max = intParam("cg_max", 10_000)  # Max number of iterations
  num_source = intParam("num_source", 16)  # Number of stochastic sources
  sq_min_distance0 = intParam("sq_min_distance", 0)
  seed = intParam("seed", int(1000*epochTime())).uint64
  write_group_size = intParam("write_group_size", 128)

var sp = initSolverParams()
sp.r2req = cg_prec*cg_prec
sp.maxits = cg_max
sp.sloppySolve = intParam("sloppy", 2).SloppyType # 0: None, 1: Single, 2: Half

# Load lattice and determine size and parameters
template getLat(fn:string): seq[int] =
  if inlat.len == 0:
    echo "WARNING: no input lattice file specicified."
    @[4,4,4,8]
  else:
    let lat = inlat.getFileLattice
    if lat.len == 0:
      echo "ERROR: getFileLattice failed on '", inlat, "'."
      qexAbort()
    lat
let
  lat = inlat.getLat
  nt = lat[^1]
  metadata_prefix = "l" & $lat[0] & ".t" & $nt & ".m" & $mass

var
  lo = lat.newLayout
  g = lo.newGauge
  r = newRNGField(RngMilc6, lo, seed)
  R: RngMilc6  # global RNG
  nc = g[0][0].nrows
R.seed(seed, 987654321)
if inlat.len == 0:
  threads:
    g.unit
    #g.random r
else:
  if 0 != g.loadGauge(inlat):
    echo "ERROR: loadGauge failed on '", inlat, "'."
    qexAbort()

echo "latsize = ",lo.physGeom
echo "volume = ",lo.physVol
echo "mass = ",mass
echo "seed = ",seed
echo "cg_prec = ",cg_prec
echo "num_source = ",num_source

var sq_min_distance = sq_min_distance0
if sq_min_distance == 0:
  let v = lo.physVol.float / num_source.float
  sq_min_distance = floor(0.9*sqrt(v)).int
echo "sq_min_distance = ", sq_min_distance

var spatv = 1
for i in 0..<lat.len-1: spatv *= lat[i]

proc printPlaq(g: any) =
  let
    p = g.plaq
    sp = 2.0*(p[0]+p[1]+p[2])
    tp = 2.0*(p[3]+p[4]+p[5])
  echo "plaq ",p
  echo "plaq ss: ",sp," st: ",tp," tot: ",p.sum

threads:
  let d = g.checkSU
  echo "unitary deviation avg: ",d.avg," max: ",d.max
g.printPlaq

threads: g.projectSU

threads:
  let d = g.checkSU
  echo "new unitary deviation avg: ",d.avg," max: ",d.max
g.printPlaq

var
  info: PerfInfo
  coef = HypCoefs(alpha1:0.4, alpha2:0.5, alpha3:0.5)
echo "smear = ",coef
var sg = lo.newGauge
coef.smear(g, sg, info)
#for i in 0..<g.len:
#  sg[i] := g[i]

threads:
  sg.setBC
  sg.stagPhase
var s = sg.newStag

var pts = newSeq[seq[int]]()
proc randomPoint(): seq[int] =
  let nd = lo.nDim
  result.newSeq(nd)
  while true:
    for i in 0..<nd:
      let li = lo[i]
      result[i] = floor(li * R.uniform()).int
      doAssert(result[i]>=0 and result[i]<li)
    #echo result
    var far = true
    for l in 0..<pts.len:
      var d = 0
      for i in 0..<nd:
        let dx1 = abs(pts[l][i] - result[i])
        let dx2 = lo[i] - dx1
        let dx = min(dx1, dx2)
        #echo " ", dx
        d += dx*dx
      #echo d
      if d < sq_min_distance:
        far = false
        break
    # Force the point to lay on a unit corner.
    for j in 0..<result.len:
      result[j] = result[j] - (result[j] mod 2)
    if far:
      pts.add result
      break

proc pointSource(r: Field; c: openArray[int]; ic: int) =
  let (ptRank,ptIndex) = r.l.rankIndex(c)
  threads:
    r := 0
    threadBarrier()
    #echo "point: ", r.norm2
    if myRank==ptRank and threadNum==0:
      r{ptIndex}[ic] := 1
    threadBarrier()
    #echo "point: ", r.norm2

proc translate(r: Field, x: Field2, pt: openArray[int]) =
  r := x
  for mu in 0..<pt.len:
    let n = pt[mu]
    #echo &"translate: {mu} {n} ", r.norm2
    # Works around a bug.
    # var s = newShifter(r, mu, n)
    # r := s ^* r
    var s = newShifter(r, mu, 1)
    for i in 1..n:
      r := s ^* r
    #echo &"translate:     ", r.norm2

var
  eta = lo.ColorVector        # Random source
  phi = lo.ColorVector        # Propagator
  n2phi = lo.Real             # Local norm2 of prop
  n2phit = lo.Real            # Local norm2 of prop after translating
  locmes = lo.Real            # Local meson field
  scal = 1.0/num_source.float
locmes := 0

tic()
for i in 0..<num_source:
  # Generate a random point to place the point source.
  var pt = randomPoint()
  #var pt = @[0,0,0,0]
  echo &"Starting inversion work on point source {i}, living at {pt}."

  for c in 0..<nc:
    echo "Color ", c
    tic()
    # Generate a point source.
    pointSource(eta, pt, c)

    # Invert
    phi := 0
    s.solve(phi, eta, mass, sp)
    echo phi.norm2

    threads:
      for i in phi:
        n2phi[i] := scal * phi[i].norm2

    translate(n2phit, n2phi, pt)

    locmes += n2phit

    var indiv_time = getElapsedTime2()
    echo "one inversion and shift time: ", indiv_time

# rephase and get correlator
var est = newSeq[float](nt)
for i in locmes.sites:
  var t: float
  t := locmes{i}
  let par = (lo.coords[0][i]+lo.coords[1][i]+lo.coords[2][i]) mod 2
  if par!=0:
    t := -t
    locmes{i} := t
  est[lo.coords[3][i]] += t
est.ranksum
for t in 0..<nt:
  echo t, "  ", est[t]

var w = lo.newWriter(outfn, metadata_prefix)
w.write(locmes, metadata_prefix, "D")
w.close

#[

# Test loading local meson
echo "Loading local meson"
var rdr = lo.newReader outfn
rdr.read n2phi
rdr.close

echo "File metadata: ", rdr.fileMetadata
echo "Record metadata: ", rdr.recordMetadata

# per-slice estimate at 0 momentum
var est2 = newseq[float](nt)
for i in n2phi.sites:
  var t: float
  t := n2phi{i}
  est2[lo.coords[3][i]] += t
est2.ranksum
for t in 0..<nt:
  echo t, "  ", est2[t]

]#


qexFinalize()
