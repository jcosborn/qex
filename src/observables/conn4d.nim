#[
  Port of my_conn4d_meas.lua
]#

import qex, gauge/hypsmear, physics/stagSolve
import times, strformat, os

qexInit()

# Accept command line arguments and set up parameters
letParam:
  # Required parameters
  inlat = ""   # Input gauge file name
  outfn = "scprop.out" # Save results to outfn
  mass = 0.1   # The quark mass

  # optional parameters
  lat =
    if existsFile(inlat):
      getFileLattice inlat
    else:
      qexWarn "Nonexistent gauge file: ", inlat
      @[4,4,4,8]
  cg_prec = 1e-9 # Max residual with default
  cg_max = 100_000  # Max number of iterations
  num_source = 16  # Number of stochastic sources
  sq_min_distance =  # squared minimum distance
    floor(0.25*sqrt(float(lat[0]*lat[1]*lat[2]*lat[3])/num_source.float)).int
  seed:uint64 = int(1000*epochTime())
  ## write_group_size = 128
  showTimers:bool = 0 # print out the timers in the end

echoParams()

echo "rank ", myRank, "/", nRanks
threads: echo "thread ", threadNum, "/", numThreads

var sp = initSolverParams()
sp.r2req = cg_prec*cg_prec
sp.maxits = cg_max
sp.sloppySolve = intParam("sloppy", 2).SloppyType # 0: None, 1: Single, 2: Half

# Load lattice and determine size and parameters
let
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
    qexError "loadGauge failed on '", inlat, "'."

echo "latsize = ",lo.physGeom
echo "volume = ",lo.physVol

qexLog "Finished loading conifguration."

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

qexLog "Finished re-unitarization."

proc smear(sg,g:any) =
  var
    info: PerfInfo
    coef = HypCoefs(alpha1:0.4, alpha2:0.5, alpha3:0.5)
  echo "smear = ",coef
  coef.smear(g, sg, info)

var sg = lo.newGauge
sg.smear g

#echo GC_getStatistics()
GC_fullCollect()
#echo GC_getStatistics()

qexLog "Finished smearing."

threads:
  sg.setBC
  sg.stagPhase
var s = sg.newStag

qexLog "Finished stagPhase."

var pts = newSeq[seq[int]]()
proc randomPoint(): seq[int] =
  let nd = lo.nDim
  result.newSeq(nd)
  const maxiter = 1_000_000_000
  for iter in 0..<maxiter:
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
      return
  qexError "max iteration reached without finding a random point.\n",
    "Perhaps sq_min_distance=", sq_min_distance, " is too large."

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

for i in 0..<num_source:
  # Generate a random point to place the point source.
  var pt = randomPoint()
  #var pt = @[0,0,0,0]
  qexLog &"Starting inversion work on point source {i}, living at {pt}."

  for c in 0..<nc:
    echo "Color ", c
    tic()
    # Generate a point source.
    pointSource(eta, pt, c)

    # Invert
    phi := 0
    s.solve(phi, eta, mass, sp)
    echo "norm2: ", phi.norm2

    threads:
      for i in phi:
        n2phi[i] := scal * phi[i].norm2

    translate(n2phit, n2phi, pt)

    locmes += n2phit

    var indiv_time = getElapsedTime()
    echo "one inversion and shift time: ", indiv_time

qexLog "Finished all sources."

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

if showTimers: echoTimers()
qexFinalize()
