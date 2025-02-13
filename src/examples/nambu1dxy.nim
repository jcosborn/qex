import qex, strformat
import utils/resample

qexInit()

letParam:
  beta = 2.0
  xlen = 16
  seed:uint64 = 4321
  ntrajThermo = 4
  ntraj = 8
  nsteps = 16
  tau = 2.0
  revCheckFreq = max(1, ntraj div 64)

  hmc:bool = 0
  gamma = 1.0

  gpChoice = "none"
  gpFactor = 1.0
  gpParam  = @[0.0]

  grChoice = "identity"
  grFactor = 1.0
  grParam  = @[0.0]

  gxChoice = "dfieldaction"
  gxFactor = 0.5*beta
  gxParam  = @[0.0]

  jkBlockSize = max(1, ntraj div 64)
echoParams()

echo "rank ", myRank, "/", nRanks
threads: echo "thread ", threadNum, "/", numThreads

let
  lat = @[xlen]
  lo = lat.newLayout

var
  nAccept = 0
  x = lo.Real
  p = lo.Real
  r = lo.Real
  xSave = lo.Real
  rng = lo.newRNGField(MRG32k3a, seed)
  globalRng: MRG32k3a  # global RNG
globalRng.seed(seed, 987654321)

proc refreshMomentum(p: auto) =
  threads:
    p.gaussian rng

proc fieldAction(x: auto, beta: float): float =
  var sx: float
  var t = newShifter(x, 0, 1)
  threads:
    var c: type(x[0])
    c := 0.0
    discard t ^* x
    for i in x:
      c += cos(t.field[i]-x[i])
    var cs = simdSum(c)
    x.l.threadRankSum(cs)
    threadSingle:
      sx = cs
  (-beta)*sx

proc dFieldAction(f,x: auto, beta: float) =
  var tf = newShifter(x, 0, 1)
  var tb = newShifter(x, 0, -1)
  threads:
    discard tf ^* x
    discard tb ^* x
    for i in f:
      f[i] := beta * ( sin(x[i]-tf.field[i]) + sin(x[i]-tb.field[i]) )

proc dTopo(f,x: auto, beta: float) =
  # \beta [\cos(x_i - x_{i-1}) - \cos(x_i - x_{i+1})]
  var tf = newShifter(x, 0, 1)
  var tb = newShifter(x, 0, -1)
  threads:
    discard tf ^* x
    discard tb ^* x
    for i in f:
      f[i] := beta * ( cos(x[i]-tb.field[i]) - cos(x[i]-tf.field[i]) )

proc sinsum(f,x: auto, beta: float) =
  # \beta \sin(\sum_j x_j)
  threads:
    let s = x.sum
    for i in f:
      f[i] := beta * sin(s)

proc cossum(f,x: auto, beta: float) =
  # \beta \sin(\sum_j x_j)
  threads:
    let s = x.sum
    for i in f:
      f[i] := beta * cos(s)

proc laplacian(f, x: auto, beta: float) =
  # \beta (x_{i+1} - 2x_i + x_{i-1}).
  let sF = newShifter(x, 0, +1)
  let sB = newShifter(x, 0, -1)
  threads:
    discard sF ^* x
    discard sB ^* x
    for i in f:
      let lap = sF.field[i] - 2.0*x[i] + sB.field[i]
      f[i] = beta * lap

proc longRange(f, x: auto, beta: float, maxr: int) =
  # \beta \sum_{r=1}^{\text{maxr}}
  #     (\sin(x_i - x_{i+r}) + \sin(x_i - x_{i-r}))
  var sF = newSeq[Shifter[typeof(x), typeof(x[0])]](maxr)
  var sB = newSeq[Shifter[typeof(x), typeof(x[0])]](maxr)
  for r in 1..maxr:
    sF[r-1] = newShifter(x, 0, +r)
    sB[r-1] = newShifter(x, 0, -r)
  threads:
    for r in 0..<maxr:
      discard sF[r] ^* x
      discard sB[r] ^* x
    for i in f:
      var xi:typeof(x[i])
      xi := 0.0
      for r in 0..<maxr:
        xi += sin(x[i] - sF[r].field[i])
        xi += sin(x[i] - sB[r].field[i])
      f[i] = beta * xi / float(maxr)

proc kick(f, x: auto, beta: float, threshold, slope: float) =
  # beta * (0.5 + 0.5 \tanh(slope * (threshold - \cos(x_{i+1}-x_i))))
  let sF = newShifter(x, 0, +1)
  threads:
    discard sF ^* x
    for i in f:
      let d = slope * (threshold - cos(sF.field[i] - x[i]))
      f[i] = beta * (0.5 + 0.5*tanh(d))

proc blockAverage(f, x: auto, beta: float, maxr: int) =
  # "Block": average +/- r neighbors
  var sF = newSeq[Shifter[typeof(x), typeof(x[0])]](maxr)
  var sB = newSeq[Shifter[typeof(x), typeof(x[0])]](maxr)
  for r in 1..maxr:
    sF[r-1] = newShifter(x, 0, +r)
    sB[r-1] = newShifter(x, 0, -r)
  let nr = float(2*maxr+1)
  threads:
    for r in 0..<maxr:
      discard sF[r] ^* x
      discard sB[r] ^* x
    for i in f:
      var xi = x[i]
      for r in 1..maxr:
        xi += sF[r-1].field[i] + sB[r-1].field[i]
      f[i] = beta * (xi/nr - x[i])

proc binomial(n, k: int): int =
  if k > n-k: 
    return binomial(n, n-k)
  result = 1
  for i in 0..<k:
    result = result*(n-i) div (i+1)

proc binomial(f, x: auto, beta: float, radius: int) =
  ## Perform a binomial smoothing of width nPoints around each site,
  ## then the derivative is beta*(weightedAverage - x[i]).
  ##
  ## E.g. if radius=2, and the binomial kernel is
  ##   [1,4,6,4,1] => sum=16 => each site i is replaced by
  ##   (1*x[i-2] + 4*x[i-1] + 6*x[i] + 4*x[i+1] + 1*x[i+2]) / 16
  ## The derivative is beta*( that average - x[i] ).
  let nPoints = radius*2+1
  var coefs = newSeq[int](nPoints)
  for k in 0..<nPoints:
    coefs[k] = binomial(nPoints-1, k)
  let sumC = float(1 shl (nPoints-1))
  var sF = newSeq[Shifter[typeof(x), typeof(x[0])]](radius)
  var sB = newSeq[Shifter[typeof(x), typeof(x[0])]](radius)
  for r in 1..radius:
    sF[r-1] = newShifter(x, 0, +r)
    sB[r-1] = newShifter(x, 0, -r)
  threads:
    for r in 0..<radius:
      discard sF[r] ^* x
      discard sB[r] ^* x
    for i in f:
      var xi = coefs[radius]*x[i]
      for r in 1..radius:
        xi += coefs[radius+r]*sF[r-1].field[i]
        xi += coefs[radius-r]*sB[r-1].field[i]
      f[i] = beta * (xi/sumC - x[i])

proc applyG(gfield, z: auto, choice: string, factor: float, param: openarray[float]) =
  case choice:
  of "none":
    threads:
      for i in gfield:
        gfield[i] := 0.0
  of "identity":
    threads:
      for i in gfield:
        gfield[i] = factor * z[i]
  of "dfieldaction":
    gfield.dFieldAction(z, factor)
  of "dtopo":
    gfield.dTopo(z, factor)
  of "sinsum":
    gfield.sinsum(z, factor)
  of "cossum":
    gfield.cossum(z, factor)
  of "laplacian":
    gfield.laplacian(z, factor)
  of "longrange":
    gfield.longRange(z, factor, max(1, int(param[0])))
  of "kick":
    gfield.kick(z, factor, param[0], param[1])
  of "block":
    gfield.blockAverage(z, factor, max(1, int(param[0])))
  of "binomial":
    gfield.binomial(z, factor, max(1, int(param[0])))
  else:
    echo "ERROR: unknown choice: ", choice
    qexError("Invalid function selection")

proc momAction(p: auto): float =
  var sp: float
  threads:
    sp := 0.5*p.norm2
  sp

proc action(p,x: auto): float =
  var sp, sx: float
  sp = momAction(p)
  sx = fieldAction(x, beta)
  sp+sx

proc action(p,x,r: auto): float =
  var sp, sr, sx: float
  sp = momAction(p)
  sr = gamma*momAction(r)
  sx = fieldAction(x, beta)
  sp+sr+sx

proc updateX(x,dp: auto, s: float) =
  threads:
    x += s * dp

proc updateP(p,f: auto, s: float) =
  threads:
    p += s * f

proc evoHMC(p,x: auto) =
  let f = lo.Real
  let eps = tau / nsteps
  updateX(x, p, 0.5*eps)
  for i in 0..<nsteps:
    if i>0:
      updateX(x, p, eps)
    f.dFieldAction(x, beta)
    updateP(p, f, -eps)
  updateX(x, p, 0.5*eps)

proc poissonX(f, p, r, gp, gr: auto) =
  # dh/dp dg/dr - dh/dr dg/dp
  threads:
    f := p*gr - gamma*r*gp

proc poissonP(f, r, dhdx, gr, gx: auto) =
  # dh/dr dg/dx - dh/dx dg/dr
  threads:
    f := gamma*r*gx - dhdx*gr

proc poissonR(f, p, dhdx, gp, gx: auto) =
  # dh/dx dg/dp - dh/dp dg/dx
  threads:
    f := dhdx*gp - p*gx

proc evoNambu(p,x,r: auto) =
  let
    dhdx = lo.Real
    gx = lo.Real  # dg/dx
    gp = lo.Real  # dg/dp
    gr = lo.Real  # dg/dr
    fx = lo.Real  # dh/dp dg/dr - dh/dr dg/dp
    fp = lo.Real  # dh/dr dg/dx - dh/dx dg/dr
    fr = lo.Real  # dh/dx dg/dp - dh/dp dg/dx
  let eps = tau / nsteps

  gp.applyG(p, gpChoice, gpFactor, gpParam)
  gr.applyG(r, grChoice, grFactor, grParam)
  fx.poissonX(p, r, gp, gr)
  updateX(x, fx, 0.5*eps)

  for i in 0..<nsteps:
    if i>0:
      fx.poissonX(p, r, gp, gr)
      updateX(x, fx, eps)

    dhdx.dFieldAction(x, beta)
    gx.applyG(x, gxChoice, gxFactor, gxParam)

    fr.poissonR(p, dhdx, gp, gx)
    updateP(r, fr, 0.5*eps)
    gr.applyG(r, grChoice, grFactor, grParam)

    fp.poissonP(r, dhdx, gr, gx)
    updateP(p, fp, eps)
    gp.applyG(p, gpChoice, gpFactor, gpParam)

    fr.poissonR(p, dhdx, gp, gx)
    updateP(r, fr, 0.5*eps)
    gr.applyG(r, grChoice, grFactor, grParam)

  fx.poissonX(p, r, gp, gr)
  updateX(x, fx, 0.5*eps)

var Evals = newSeq[float](ntraj)
var Qvals = newSeq[float](ntraj)
var Avals = newSeq[float](ntraj)

proc printObservables(traj:int, acc:bool, ds,pacc,r:float, x: auto) =
  var t = newShifter(x, 0, 1)
  var rc,rs,rq: float
  threads:
    var c,s,q: type(x[0])
    c := 0.0
    s := 0.0
    q := 0.0
    discard t ^* x
    for i in x:
      let
        d = t.field[i]-x[i]
        cd = cos(d)
        sd = sin(d)
      c += cd
      s += sd
      q += atan2(sd,cd)
    var
      cs = simdSum(c)
      ss = simdSum(s)
      qs = simdSum(q)
    x.l.threadRankSum(cs)
    x.l.threadRankSum(ss)
    x.l.threadRankSum(qs)
    threadSingle:
      rc = cs/float(lo.physVol)
      rs = ss/TAU
      rq = round(qs/TAU)
  let ar = if acc: "Accept" else: "Reject"
  echo &"traj: {traj}  {ar}  dH: {ds}  exp(-dH): {pacc}  r: {r}  energy: {rc}  qr: {rs}  q: {rq}"
  if traj>0:
    Evals[traj-1] = rc
    Qvals[traj-1] = rq

threads:
  x := 0

var revCorrect = true
for traj in 1-ntrajThermo..ntraj:
  threads:
    xSave := x
  var s0, s1: float
  if hmc:
    refreshMomentum(p)
    s0 = action(p, x)
    evoHMC(p, x)
    s1 = action(p, x)
    if revCheckFreq>0 and traj>0 and traj mod revCheckFreq == 0:
      let
        x1 = lo.Real
        p1 = lo.Real
      threads:
        x1 := x
        p1 := p
        p := -p
      evoHMC(p, x)
      let ds = action(p, x)-s0
      var dx: float
      threads:
        x := x - xSave
        dx = x.norm2/float(xlen)
      echo "Reversibility: dx2: ", dx,"  dH: ",ds
      if ds>1e-7 or dx>1e-14:
        revCorrect = false
        echo "WARNING: broken reversibility."
      threads:
        x := x1
        p := p1
  else:
    refreshMomentum(p)
    refreshMomentum(r)
    let isqrtgamma = 1.0/sqrt(gamma)
    threads:
      r := isqrtgamma*r
    s0 = action(p, x, r)
    evoNambu(p, x, r)
    s1 = action(p, x, r)
    if revCheckFreq>0 and traj>0 and traj mod revCheckFreq == 0:
      let
        x1 = lo.Real
        p1 = lo.Real
        r1 = lo.Real
      threads:
        x1 := x
        p1 := p
        r1 := r
        p := -p
      evoNambu(p, x, r)
      let ds = action(p, x, r)-s0
      var dx: float
      threads:
        x := x - xSave
        dx = x.norm2/float(xlen)
      echo "Reversibility: dx2: ", dx,"  dH: ",ds
      if ds>1e-7 or dx>1e-14:
        revCorrect = false
        echo "WARNING: failed reversibility."
      threads:
        x := x1
        p := p1
        r := r1
  let ds = s1 - s0
  let pacc = exp(-ds)
  if traj>0:
    Avals[traj-1] = pacc
  let racc = globalRng.uniform
  let accept = racc <= pacc
  if accept:
    if traj>0:
      inc nAccept
  else:
    threads:
      x := xSave

  printObservables(traj, accept, ds, pacc, racc, x)
  threads:
    for i in x:
      let t = x[i]
      x[i] := atan2(sin(t), cos(t))

const EPS = 1e-15

# I0(x) = sum_{k=0 to infinity} ( (x/2)^(2k) ) / (k!)^2
proc besselI0(x: float): float =
  let halfX = x / 2.0
  var
    term = 1.0   # term for k = 0
    sumVal = 1.0 # sum starts with k = 0
    k = 1
  while abs(term) > EPS * abs(sumVal):
    # term *= ( (x/2)^2 ) / (k^2 )
    term *= (halfX * halfX) / (float(k) * float(k))
    sumVal += term
    inc k
  sumVal

# I1(x) = sum_{k=0 to infinity} ( (x/2)^(2k+1) ) / ( k! * (k+1)! )
proc besselI1(x: float): float =
  let halfX = x / 2.0
  var
    term = halfX # term for k = 0
    sumVal = halfX
    k = 1
  while abs(term) > EPS * abs(sumVal):
    # term *= ( (x/2)^2 ) / [ k * (k+1) ]
    term *= (halfX * halfX) / (float(k) * float(k + 1))
    sumVal += term
    inc k
  sumVal

proc mean(xs: Ensemble[seq[float]]): float =
  var m = 0.0
  for i in 0..<xs.len:
    let x = xs[i]
    m += (x - m)/float(i+1)
  m

proc naiveIntAutocorr(xs: Ensemble[seq[float]], maxlen: int): float =
  ## Very rough integrated autocorrelation: sum of correlation out to some cut.
  let N = xs.len
  let m = mean(xs)
  # c(0)
  var c0 = 0.0
  for i in 0..<xs.len:
    let d = xs[i]-m
    c0 += d*d
  c0 /= float(N)
  if c0 <= 1e-14: return 1.0

  var sumRho = 1.0
  for lag in 1..(min(maxlen, N-1)):
    var cLag = 0.0
    for i in 0..<(N-lag):
      cLag += (xs[i] - m)*(xs[i+lag] - m)
    cLag /= float(N-lag)
    let rho = cLag/c0
    if rho < 0.0:
      break
    sumRho += 2.0*rho
  sumRho

proc tunnelingRate(xs: Ensemble[seq[float]]): float =
  var r = 0.0
  for i in 1..<xs.len:
    let d = xs[i]-xs[i-1]
    r += (d*d-r)/float(i)
  sqrt(r)

if not revCorrect:
  echo "WARNING: Reversibility check failed."
echo "Acceptance ratio: ", nAccept/ntraj
echo "Jackknife Block Size: ", jkBlockSize

let expmdh = Avals.jackknife(jkBlockSize, mean)
for a in Avals.mitems:
  if a>1.0:
    a = 1.0
let pacc = Avals.jackknife(jkBlockSize, mean)
let Emean = Evals.jackknife(jkBlockSize, mean)
let Qmean = Qvals.jackknife(jkBlockSize, mean)
let Qac = Qvals.jackknife(jkBlockSize, naiveIntAutocorr, 200)
let dQrms = Qvals.jackknife(jkBlockSize, tunnelingRate)

var qmax = 0
for qv in Qvals:
  let q = int(round(qv))
  if q>qmax:
    qmax = q
  elif -q>qmax:
    qmax = -q
var qdist = newseq[JackknifeStat[float]]()
proc count(xs:Ensemble[seq[float]], x:int):float =
  for i in 0..<xs.len:
    if x == int(abs(round(xs[i]))):
      result += 1.0
  if x!=0:
    result *= 0.5
for i in 0..qmax:
  qdist.add Qvals.jackknife(jkBlockSize, count, i)

for a in Qvals.mitems:
  a = a*a
let Q2 = Qvals.jackknife(jkBlockSize, mean)
let Q2ac = Qvals.jackknife(jkBlockSize, naiveIntAutocorr, 200)

echo "exp(-dH) = ", expmdh.mean, " ± ", expmdh.stdev
echo "Pacc = ",pacc.mean, " ± ", pacc.stdev
echo "Emean = ", Emean.mean, " ± ", Emean.stdev, "  dE = ",Emean.mean-besselI1(beta)/besselI0(beta)
echo "Qmean = ", Qmean.mean, " ± ", Qmean.stdev
echo "Tau_Q = ", Qac.mean, " ± ", Qac.stdev
echo "dQrms = ", dQrms.mean, " ± ", dQrms.stdev
echo "Q2/L = ", Q2.mean/float(xlen), " ± ", Q2.stdev/float(xlen)
echo "Tau_Q2 = ", Q2ac.mean, " ± ", Q2ac.stdev
for i in 0..qmax:
  echo "P(Q=",i,") = ",qdist[i].mean/float(ntraj), " ± ", qdist[i].stdev/float(ntraj)

qexFinalize()
