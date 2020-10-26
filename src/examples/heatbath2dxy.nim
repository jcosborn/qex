import qex
import gauge, physics/qcdTypes
import os, strutils, times

proc sumEnergy(fr,fi:any, J, h:any, g:any, sf,sb:any) =
  fr := 0
  fi := 0
  for nu in 0..<g.l.nDim:
    discard sf[nu] ^* g
    discard sb[nu] ^* g
    threadBarrier()
    for i in fr:
      fr[i] += cos(sf[nu].field[i]) + cos(sb[nu].field[i])
    for i in fi:
      fi[i] += sin(sf[nu].field[i]) + sin(sb[nu].field[i])
  fr := J*fr + h
  fi := J*fi

type HeatBath[F,E] = object
  fr,fi: F
  sf,sb: array[2,seq[Shifter[F,E]]]
  subs: array[2,Subset]

proc newHeatBath(lo:any):auto =
  let
    fr = lo.Real
    fi = lo.Real
  type
    F = typeof(fr)
    E = typeof(fr[0])
  var r = HeatBath[F,E](fr:fr, fi:fi)
  const p = ["even","odd"]
  for j in 0..1:
    r.sf[j] = newseq[Shifter[F,E]](lo.nDim)
    r.sb[j] = newseq[Shifter[F,E]](lo.nDim)
    for i in 0..<lo.nDim:
      r.sf[j][i] = newShifter(fr, i, 1, p[j])
      r.sb[j][i] = newShifter(fr, i, -1, p[j])
  r.subs[0].layoutSubset(lo,"e")
  r.subs[1].layoutSubset(lo,"o")
  r

proc evolve(H:HeatBath, g:any, gc:any, r:any) =
  tic("heatbath")
  let
    lo = g.l
    nd = lo.nDim
    (beta, J, h) = gc
  if H.subs.len != 2:
    qexError "HeatBath only works with even-odd subsets for now."
  threads:
    tic("threads")
    # sample
    for j in 0..<H.subs.len:
      let
        s = H.subs[j]
        so = H.subs[(j+1) mod 2]
      sumEnergy(H.fr[s], H.fi[s], J, h, g, H.sf[j], H.sb[j])
      threadBarrier()
      for i in g[s].sites:
        let
          yr = H.fr{i}[][]
          yi = H.fi{i}[][]
          lambda = beta*hypot(yi, yr)
          phi = arctan2(yi, yr)
        g{i} := vonMises(r{i}, lambda)+phi
    threadBarrier()
    toc("sample")
    # over relaxation: flip
    for j in 0..<H.subs.len:
      let
        s = H.subs[j]
        so = H.subs[(j+1) mod 2]
      sumEnergy(H.fr[s], H.fi[s], J, h, g, H.sf[j], H.sb[j])
      threadBarrier()
      for i in g[s].sites:
        let
          yr = H.fr{i}[][]
          yi = H.fi{i}[][]
        g{i} := 2.0*arctan2(yi,yr)-g{i}[][]
    toc("flip")
  toc("end")

proc magnet(g:any):auto =
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
  let vol = g.l.physVol.float
  (mr/vol, mi/vol)

proc showMeasure(g:any,label="") =
  let
    (mr,mi) = g.magnet
    v = g.l.physVol.float
    s = (mr*mr+mi*mi)*v
  echo label,"magnet: ",mr," ",mi," ",s

qexinit()
tic()

letParam:
  #lat = @[8,8,8,8]
  #lat = @[8,8,8]
  lat = @[32,32]
  #lat = @[1024,1024]
  beta = 1.0
  J = 1.0
  h = 0.0
  sweeps = 10
  seed:uint64 = int(1000*epochTime())

echoParams()
echo "rank ", myRank, "/", nRanks
threads: echo "thread ", threadNum, "/", numThreads

let
  gc = (beta:beta, J:J, h:h)
  lo = lat.newLayout
  vol = lo.physVol

var r = lo.newRNGField(RngMilc6, seed)

var g = lo.Real

threads:
  for i in g.sites:
    let u = uniform r{i}
    g{i} := PI*(2.0*u-1.0)

g.showMeasure "Initial: "

var H = lo.newHeatBath

toc("init")

for n in 1..sweeps:
  tic("sweep")
  echo "Begin sweep: ",n

  H.evolve(g,gc,r)
  toc("evolve")

  g.showMeasure
  toc("measure")

toc("done")
echoTimers()
qexfinalize()
