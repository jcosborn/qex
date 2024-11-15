import qex
import gauge, physics/qcdTypes
import os, strutils, times

const
  twistDir = 0
  twistBoundaryIndex = 0
  NMax = 8

proc sumEnergy(fr,fi:auto, J, h:auto, g,b:auto, bb,sf,sb:auto) =
  fr := 0
  fi := 0
  for nu in 0..<g.l.nDim:
    discard sf[nu] ^* g
    discard sb[nu] ^* g
    var pf,cpf,spf,pb,cpb,spb:typeof(g[0])
    threadBarrier()
    for i in fr:
      pf := sf[nu].field[i]
      pb := sb[nu].field[i]
      if nu == twistDir:
        pf -= b[i]
        pb += bb.field[i]
      cpf = cos pf
      spf = sin pf
      cpb = cos pb
      spb = sin pb
      fr[i] += cpf + cpb
      fi[i] += spf + spb
  fr := J*fr + h
  fi := J*fi

proc pickClock[D](r:var RNG, N:float, lambda,phi:D):D =
  var p {.noinit,align(64).}:array[NMax,D]
  let n = int ceil N
  for i in 0..<n:
    p[i] = lambda*cos(2.0*PI*float(i)/N - phi)
  var pmax = p[0]
  for i in 1..<n:
    if p[i]>pmax:
      pmax = p[i]
  for i in 0..<n:
    p[i] -= pmax
  var psum = 0.0
  for i in 0..<n:
    psum += exp(p[i])
    p[i] = psum
  let u = psum * r.uniform
  for i in 0..<n:
    if u<p[i]:
      return 2.0*PI*float(i)/N

type PhaseDiff[F,E] = object
  cosd,sind:seq[float]
  f: seq[Shifter[F,E]]

proc phaseDiffB(del:var PhaseDiff,g:auto):auto =
  let
    # del cannot be captured by nim in threads
    f = del.f
    cosd = cast[ptr UnCheckedArray[float]](del.cosd[0].addr)
    sind = cast[ptr UnCheckedArray[float]](del.sind[0].addr)
  threads:
    var d,t,s:typeof(g{0}[][])
    discard f[twistDir] ^* g
    threadBarrier()
    for i in g.sites:
      if g.l.coords[twistDir][i] == twistBoundaryIndex:
        d = f[twistDir].field{i}[][] - g{i}[][]
        t += cos(d)
        s += sin(d)
    t.threadRankSum
    s.threadRankSum
    threadSingle:
      cosd[twistDir] = t
      sind[twistDir] = s

proc phaseDiff(del:var PhaseDiff,g,b:auto):auto =
  let
    # del cannot be captured by nim in threads
    f = del.f
    cosd = cast[ptr UnCheckedArray[float]](del.cosd[0].addr)
    sind = cast[ptr UnCheckedArray[float]](del.sind[0].addr)
  threads:
    for nu in 0..<g.l.nDim:
      var d,t,s:typeof(g[0])
      discard f[nu] ^* g
      threadBarrier()
      for i in g:
        d := f[nu].field[i] - g[i]
        if nu == twistDir:
          d -= b[i]
        t += cos(d)
        s += sin(d)
      var
        v = t.simdSum
        u = s.simdSum
      v.threadRankSum
      u.threadRankSum
      threadSingle:
        cosd[nu] = v
        sind[nu] = u

type HeatBath[F,E] = object
  fr,fi: F
  sf,sb: array[2,seq[Shifter[F,E]]]
  subs: array[2,Subset]
  del: PhaseDiff[F,E]

proc newHeatBath(lo:auto):auto =
  let
    nd = lo.nDim
    fr = lo.Real
    fi = lo.Real
  type
    F = typeof(fr)
    E = typeof(fr[0])
  var r = HeatBath[F,E](fr:fr, fi:fi)
  const p = ["even","odd"]
  for j in 0..1:
    r.sf[j] = newseq[Shifter[F,E]](nd)
    r.sb[j] = newseq[Shifter[F,E]](nd)
    for i in 0..<nd:
      r.sf[j][i] = newShifter(fr, i, 1, p[j])
      r.sb[j][i] = newShifter(fr, i, -1, p[j])
  r.subs[0].layoutSubset(lo,"e")
  r.subs[1].layoutSubset(lo,"o")
  r.del.cosd = newseq[float](nd)
  r.del.sind = newseq[float](nd)
  r.del.f = newseq[Shifter[F,E]](nd)
  for i in 0..<nd:
    r.del.f[i] = newShifter(fr, i, 1)
  r

proc evolve(H:HeatBath, g,b:auto, bb:auto, d:var float, gc:auto, r:auto, R:var RngMilc6,
    sample = true, twistSample = true) =
  tic("heatbath")
  let
    lo = g.l
    nd = lo.nDim
    (beta, N, J, h) = gc
    p = d
    z = newseq[float](nd)
  if H.subs.len != 2:
    qexError "HeatBath only works with even-odd subsets for now."
  if sample:
    tic("threads")
    threads:
      discard bb ^* b
      threadBarrier()
      for j in 0..<H.subs.len:
        let
          s = H.subs[j]
          so = H.subs[(j+1) mod 2]
        sumEnergy(H.fr[s], H.fi[s], J, h, g,b,bb, H.sf[j], H.sb[j])
        threadBarrier()
        for i in g[s].sites:
          let
            yr = H.fr{i}[][]
            yi = H.fi{i}[][]
            lambda = beta*hypot(yi, yr)
            phi = arctan2(yi, yr)
          g{i} := pickClock(r{i}, N, lambda, phi)
    toc("sample")
  if twistSample:
    tic()
    var del = H.del
    del.phaseDiffB g
    let
      yr = del.cosd[twistDir]
      yi = del.sind[twistDir]
      phi = arctan2(yi,yr)
    d = floormod(vonMises(R, beta*J*hypot(yi,yr))+phi+PI,2*PI) - PI
    let c = d
    threads:
      for i in b.sites:
        if b.l.coords[twistDir][i] == twistBoundaryIndex:
          b{i} := c
    toc("twist sample")
  toc("end")

proc magnet(g:auto):auto =
  tic("magnet")
  var mr,mi = 0.0
  threads:
    var t,s:typeof(g[0])
    for i in g:
      t += cos(g[i])
      s += sin(g[i])
    var
      v = t.simdSum
      u = s.simdSum
    v.threadRankSum
    u.threadRankSum
    threadSingle:
      mr = v
      mi = u
  toc("done")
  (mr, mi)

proc showMeasure[F,E](del:var PhaseDiff[F,E],g,b:F,label="") =
  let
    (mr,mi) = g.magnet
    v = 1.0/g.l.physVol.float
    s = (mr*mr+mi*mi)*v
    nd = g.l.nDim
  del.phaseDiff(g,b)
  echo label,"magnet: ",mr," ",mi," ",s
  var diff = ""
  for i in 0..<nd:
    diff &= "CosSinDel" & $i & ": " & $(del.cosd[i]*v) & " " & $(del.sind[i]*v) & "\t"
  diff.setlen(diff.len-1)
  echo label,diff

proc showTwist[T](d:T,label="") =
  echo label, "twist: ", d

proc hitFreq(num, freq:int):bool = freq>0 and 0==num mod freq

qexinit()
tic()

letParam:
  #lat = @[8,8,8,8]
  #lat = @[8,8,8]
  lat = @[32,32]
  #lat = @[1024,1024]
  beta = 1.0
  N = 4.0
  J = 1.0
  h = 0.0
  sweeps = 10
  sampleFreq = 1
  twistSampleFreq = 1
  twistAngle = 0.0
  measureFreq = 1
  seed:uint64 = int(1000*epochTime())

echoParams()
echo "rank ", myRank, "/", nRanks
threads: echo "thread ", threadNum, "/", numThreads

if N>NMax:
  qexError "N(" & $N & ") is larger than maximum allowed NMax(" & $NMax & ")."

let
  gc = (beta:beta, N:N, J:J, h:h)
  lo = lat.newLayout
  vol = lo.physVol

var
  r = lo.newRNGField(RngMilc6, seed)
  R:RngMilc6  # global RNG for the twisting angle
  g = lo.Real
  d = twistAngle
  b = lo.Real
  bb = newShifter(b, twistDir, -1)
  H = lo.newHeatBath
R.seed(seed,987654321)

threads:
  b := 0.0
  for i in g.sites:
    let u = floor(ceil(N)*uniform(r{i}))
    g{i} := 2.0*PI*u/N
    if lo.coords[twistDir][i] == twistBoundaryIndex:
      b{i} := d

d.showTwist("Initial: ")
H.del.showMeasure(g,b, "Initial: ")

toc("init")

for n in 1..sweeps:
  tic("sweep")
  echo "Begin sweep: ",n

  H.evolve(g,b,bb,d,gc,r,R,
    hitFreq(n,sampleFreq),
    hitFreq(n,twistSampleFreq))
  toc("evolve")

  d.showTwist
  if hitFreq(n,measureFreq): H.del.showMeasure(g,b)
  toc("measure")

toc("done")
echoTimers()
qexfinalize()
