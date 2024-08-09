# coupled harmonic oscillators example (currently only 1d degree of freedom)
# H = (1/2) sum_<i,j> (x_i - x_j)^2
import qex, strformat

qexInit()

let
  lat = intSeqParam("lat", @[16])
  lo = lat.newLayout
  seed = uint64 intParam("seed", 1)
  ntraj = intParam("ntraj", 2)
  nsteps = intParam("nsteps", 2)
  tau = floatParam("tau", 1)

var
  nAccept = 0
  x = lo.Real
  p = lo.Real
  xSave = lo.Real
  #pSave = lo.Real
  rng = lo.newRNGField(MRG32k3a, seed)
  globalRng: MRG32k3a  # global RNG
globalRng.seed(seed, 987654321)

echo "ntraj: ", ntraj
echo "nsteps: ", nsteps
echo "tau: ", tau

proc refreshMomentum(p: auto) =
  threads:
    p.gaussian rng

proc action(p,x: auto): float =
  var sp, sx: float
  threads:
    sp := 0.5*p.norm2
  var t = newShifters(x, 1)
  let nd = lo.nDim
  for mu in 0..<nd:
    sx += 0.5*norm2diff(x, t[mu] ^* x)
  result = sp + sx

proc updateX(x,p: auto, s: float) =
  threads:
    x += s * p

proc updateP(p,x: auto, s: float) =
  var tf = newShifters(x, 1)
  var tb = newShifters(x, -1)
  let nd = lo.nDim
  threads:
    for mu in 0..<nd:
      discard tf[mu] ^* x
      discard tb[mu] ^* x
      p += s * ( tf[mu].field + tb[mu].field - 2*x)

proc evolve(p,x: auto) =
  let eps = tau / nsteps
  for i in 0..<nsteps:
    updateX(x, p, 0.5*eps)
    updateP(p, x, eps)
    updateX(x, p, 0.5*eps)

proc recenter(x: auto) =
  let xs = x.sum / lo.physVol
  x -= xs

proc printObservables(x: auto) =
  let xs = x.sum / lo.physVol
  let x2s = x.norm2 / lo.physVol
  echo "  ave x: ", xs
  echo "  ave x2: ", x2s

threads:
  x := 0

var ds2 = 0.0
for traj in 1..ntraj:
  refreshMomentum(p)
  threads:
    xSave := x
    #pSave := p
  let s0 = action(p, x)
  evolve(p, x)
  let s1 = action(p, x)
  let ds = s1 - s0
  ds2 += ds*ds
  let pacc = exp(-ds)
  let r = globalRng.uniform
  if r <= pacc: # accept
    echo &"Accept: {ds} {pacc} {r}"
    inc nAccept
  else: # reject
    echo &"Reject: {ds} {pacc} {r}"
    threads:
      x := xSave

  recenter(x)
  printObservables(x)

echo "Acceptance ratio: ", nAccept/ntraj
echo "ds2: ", ds2/ntraj
qexFinalize()
