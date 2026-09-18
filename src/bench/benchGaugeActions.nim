import qex
import base/[alignedMem, profile]

qexInit()
letParam:
  lat = latticeFromLocalLattice(@[4,4,4,4], nRanks)
  reps = 512
  trials = 5
installHelpParam()
processHelpParam()
let
  lo = lat.newLayout
  cases = [
    ("Wilson", GaugeActionCoeffs(plaq: 5.4)),
    ("Symanzik", Symanzik(5.4)),
    ("Iwasaki", Iwasaki(5.4)),
    ("DBW2", DBW2(5.4)),
    ("adjoint", GaugeActionCoeffs(plaq: 5.4, adjplaq: 0.2))]
var
  rng = lo.newRNGField(Philox4x64, 129911u64)
  g = lo.newGauge
  f = lo.newGauge
  total = 0.0
let work = newLoopWork(g[0])
threads:
  g.random rng
setRawMemGcThreshold(high(int))
echo "family,operation,trial,seconds_per_call,allocated_bytes_per_call"
for (name,c) in cases:
  for force in [false,true]:
    template step =
      if force:
        c.force(g, f, work=work)
      else:
        total += c.action(g, work=work)
    for _ in 0..<4: step
    for trial in 0..<trials:
      GC_fullCollect()
      let before = getRawMemAllocated()
      tic("gauge action sample")
      for _ in 0..<reps: step
      let secs = getElapsedTime()/float(reps)
      toc("sample")
      let bytes = float(getRawMemAllocated()-before)/float(reps)
      echo name, ",", (if force: "force" else: "action"), ",", trial, ",", secs, ",", bytes
echo "action checksum: ", total
qexFinalize()
