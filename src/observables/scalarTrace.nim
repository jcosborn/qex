#[
  Port of scalar propagator code from FUEL.

  Only 4D SU(3) for now.

  Original comments:

  JXY 2017-09-20
  Scale trace by 1/Nc
  output PBP spatial sum normalized by SpatialVolume

  JXY 2016-12-07
  Combine even-odd and corner to one file
  Incorporating James's single file writer.

  ESW 2016-11-14
  Script file that generates noise sources, dilutes in time and corners, then saves packed traces.
  Optionally, this script saves the original noise sources and propagators.
]#

import qex, gauge/hypsmear, physics/stagSolve
import times, os

qexInit()

# Accept command line arguments and set up parameters
letParam:
  # Required parameters
  inlat = ""   # Input gauge file name
  outfn = "output" # Save results to outfn
  mass = 0.1   # The quark mass

  # Intermediate definitions
  lat =
    if fileExists(inlat):
      getFileLattice inlat
    else:
      qexWarn "Nonexistent gauge file: ", inlat
      @[4,4,4,8]

  # optional parameters
  sloppy:SloppyType = 2 # 0: None, 1: Single, 2: Half
  cg_prec = 1e-9 # Max residual with default
  cg_max = 100_000  # Max number of iterations
  num_stoch = 1  # Number of stochastic sources
  improved_trace:bool = 1
  save_props:bool = 0
  source_type = "Z4" # Z4, Z2, U1, Gauss
  dilute_type:parseDilution = "EO" # EO, CORNER
  seed:uint64 = int(1000*epochTime())
  ## write_group_size = 128
  showTimers:bool = 0 # print out the timers in the end

echoParams()

echo "rank ", myRank, "/", nRanks
threads: echo "thread ", threadNum, "/", numThreads

var sp = initSolverParams()
sp.r2req = cg_prec*cg_prec
sp.maxits = cg_max
sp.sloppySolve = sloppy

# Load lattice and determine size and parameters
let
  nt = lat[^1]
  metadata_prefix = "l" & $lat[0] & ".t" & $nt & ".m" & $mass & ".cfg" & inlat

proc trace_file(i:int):string =
  let imp = if improved_trace: "1" else: "0"
  outfn & ".trace" & $i & "." & source_type & ".imp" & imp & "." & $dilute_type
proc trace_meta(i:int):string =
  let imp = if improved_trace: "1" else: "0"
  metadata_prefix & ".type" & source_type & ".src" & $i & ".imp" & imp & "." & $dilute_type

var
  lo = lat.newLayout
  g = lo.newGauge
  r = newRNGField(RngMilc6, lo, seed)
if inlat.len == 0:
  threads: g.random r
else:
  if 0 != g.loadGauge(inlat):
    echo "ERROR: loadGauge failed on '", inlat, "'."
    qexAbort()

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

var
  eta = lo.ColorVector        # Noise source
  phi = lo.ColorVector        # Propagator
  trce = lo.Complex     # Trace
  tmps = lo.ColorVector

for i in 0..<num_stoch:
  qexLog "Starting work on noise source ",i," ."
  phi := 0
  trce := 0

  echo "Generating a ",source_type," noise source."
  threads:
    case source_type
    of "Z4": eta.z4 r
    of "Z2": eta.z2 r
    of "U1": eta.u1 r
    of "Gauss":
      eta.gaussian r
      threadBarrier()
      eta *= 1.0/sqrt(2.0)
    else:
      qexError "Invalid noise type ",source_type,"."
    threadBarrier()
    echo "noise norm2: ",eta.norm2

  if save_props:              # XXX save eta
    qexWarn "Save not implemented"

  for t in 0..<nt:
    for dl in dilution(dilute_type):
      qexLog "Source ",i," dilution pattern ",dl," timeslice ",t," ."
      threads:
        tmps := 0
        threadBarrier()

        # XXX implement subset later
        for i in tmps.sites(dl):
          if lo.coords[^1][i] == t:
            tmps{i} := eta{i}
            #forO c, 0, tmps{0}.len-1:
            #  tmps{i}[c].re := eta{i}[c].re
            #  tmps{i}[c].im := eta{i}[c].im
        phi := 0
        threadBarrier()
        echo "src norm2: ",tmps.norm2
      s.solve(phi, tmps, mass, sp)
      threads:
        echo "dest norm2: ",phi.norm2
      qexLog "Finished one solve."

      # if save_props:        # XXX save phi
      #   echo "WARNING: save not implemented"

      threads:
        if improved_trace:
          echo "Computing the improved trace."
          for i in trce:
            trce[i] += mass * phi[i].dot phi[i]
        else:
          echo "Computing the unimproved trace."
          for i in trce:
            trce[i] += tmps[i].dot phi[i]
      qexLog "Finished one trace."
  qexLog "Finished source ",i," ."

  threads:
    let invnc = 1.0 / tmps[0].len.float
    trce *= invnc

  var est = newseq[float](nt)
  for i in trce.sites:
    var t:float
    t := trce{i}.re
    est[lo.coords[3][i]] += t
  est.ranksum
  for t in 0..<nt:
    echo "initsrc ",i," timeslice ",t," pbp ",est[t] / spatv.float
  qexLog "Finished computing source ",i," ."

  var trceWriter = lo.newWriter(trace_file(i), trace_meta(i))
  trceWriter.write(trce, trace_meta(i), "D")
  trceWriter.close
  qexLog "Finished saving trace."

# Test loading traces
for i in 0..<num_stoch:
  qexLog "Loading trace from source ", i

  var trceReader = lo.newReader trace_file(i)
  trceReader.read trce
  trceReader.close

  echo "File metadata for trace ", i, ": ", trceReader.fileMetadata
  echo "Trace metadata for trace ", i, ": ", trceReader.recordMetadata

  # per-slice estimate at 0 momentum
  var est = newseq[float](nt)
  for i in trce.sites:
    var t:float
    t := trce{i}.re
    est[lo.coords[3][i]] += t
  est.ranksum
  for t in 0..<nt:
    echo "loadsrc ",i," mom 0 0 0 timeslice ",t," pbp ",est[t] / spatv.float

  # per-slice estimate at +z momentum XXX

if showTimers: echoTimers()
qexFinalize()
