import qex
import gauge, physics/qcdTypes
import os, strutils, times

proc sumStaples(f:auto, g:auto, mu:int, sg:auto, ftmp:auto, sftmp:auto) =
  f := 0
  for nu in 0..<g.len:
    if mu==nu: continue
    discard sg[nu] ^* g[mu]
    discard sg[mu] ^* g[nu]
    threadBarrier()
    ftmp := g[nu].adj * g[mu] * sg[mu].field
    threadBarrier()
    discard sftmp[nu] ^* ftmp
    f += g[nu] * sg[nu].field * sg[mu].field.adj
    threadBarrier()
    f += sftmp[nu].field

type HeatBath[F,E] = object
  f: F
  ftmp: F
  sf: seq[Shifter[F,E]]
  sftmp: array[2,seq[Shifter[F,E]]]
  subs: array[2,Subset]

proc newHeatBath(lo:auto):auto =
  type F = typeof(lo.ColorMatrix())
  var f = lo.ColorMatrix()
  type E = type(f[0])
  var r = HeatBath[F,E](f: f, ftmp: newOneOf(f))
  r.sf = newseq[Shifter[F,E]](lo.nDim)
  for i in 0..<lo.nDim:
    r.sf[i] = newShifter(r.f, i, 1)
  const p = ["even","odd"]
  for j in 0..1:
    r.sftmp[j] = newseq[Shifter[F,E]](lo.nDim)
    for i in 0..<lo.nDim:
      r.sftmp[j][i] = newShifter(r.ftmp, i, -1, p[j])
  r.subs[0].layoutSubset(lo,"e")
  r.subs[1].layoutSubset(lo,"o")
  r

proc evolve(h:HeatBath, g:array or seq, gc:GaugeActionCoeffs, r:auto) =
  tic("heatbath")
  const nc = g[0][0].nrows
  let
    lo = g[0].l
    nd = lo.nDim
    bn = gc.plaq/nc.float
  when nc!=1:
    qexError "HeatBath only works with nc = 1 for now."
  if gc.rect!=0 or gc.pgm!=0 or gc.adjplaq!=0:
    qexError "HeatBath only works with plaq action for now."
  if h.subs.len != 2:
    qexError "HeatBath only works with even-odd subsets for now."
  threads:
    # sample
    for j in 0..<h.subs.len:
      let
        s = h.subs[j]
        so = h.subs[(j+1) mod 2]
      for mu in 0..<nd:
        sumStaples(h.f[s], g, mu, h.sf, h.ftmp[so], h.sftmp[j])
        threadBarrier()
        for i in g[mu][s].sites:
          let  # TODO only works with nc=1
            y = h.f{i}[0,0]
            yr = y.re[][]
            yi = y.im[][]
            lambda = bn*hypot(yi, yr)
            phi = arctan2(yi, yr)
            x = vonMises(r{i}, lambda)+phi
          g[mu]{i}[0,0].re = cos(x)
          g[mu]{i}[0,0].im = sin(x)
    threadBarrier()
    # over relaxation: flip
    for j in 0..<h.subs.len:
      let
        s = h.subs[j]
        so = h.subs[(j+1) mod 2]
      for mu in 0..<nd:
        sumStaples(h.f[s], g, mu, h.sf, h.ftmp[so], h.sftmp[j])
        threadBarrier()
        for i in g[mu][s].sites:
          let  # TODO only works with nc=1
            y = h.f{i}[0,0]
            yr = y.re[][]
            yi = -y.im[][]
            yr2 = yr*yr
            yi2 = yi*yi
            yn2 = 1.0/(yr2+yi2)
            y2r = (yr2-yi2)*yn2
            y2i = (-2.0)*yi*yr*yn2
            z = g[mu]{i}[0,0]
            zr = z.re[][]
            zi = z.im[][]
            # x <- z.adj * y.adj / y
            xr = y2r*zr+y2i*zi
            xi = y2i*zr-y2r*zi
          g[mu]{i}[0,0].re = xr
          g[mu]{i}[0,0].im = xi
  toc("end")

proc topo2DU1(g:array or seq):float =
  tic()
  #const nc = g[0][0].nrows
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

qexinit()

letParam:
  #lat = @[8,8,8,8]
  #lat = @[8,8,8]
  lat = @[32,32]
  #lat = @[1024,1024]
  beta = 6.0
  sweeps = 10
  seed:uint64 = int(1000*epochTime())

echoParams()
echo "rank ", myRank, "/", nRanks
threads: echo "thread ", threadNum, "/", numThreads

let
  gc = GaugeActionCoeffs(plaq: beta)
  lo = lat.newLayout
  #vol = lo.physVol

var r = lo.newRNGField(RngMilc6, seed)

var g = lo.newgauge
g.random r
#g.unit

#echo g.plaq
#echo g.plaq2
echo "Initial plaq: ",g.plaq3
#echo g.gaugeAction2 gc

var H = lo.newHeatBath

for n in 1..sweeps:
  echo "Begin sweep: ",n

  H.evolve(g,gc,r)

  let ga = g.gaugeAction2 gc
  echo "Sg: ",ga

  #echo g.plaq
  #echo g.plaq2
  let pl = g.plaq3
  echo "plaq: ",pl.re," ",pl.im
  echo "topo: ",g.topo2DU1

echoTimers()
qexfinalize()
