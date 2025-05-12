import base
import layout
import gauge

proc plaqEvenOddS4*(g:auto):auto =
  let
    lo = g[0].l
    nd = lo.nDim
    nc = g[0][0].ncols
  var
    pl = lo.Real() 
    peo = newseq[array[2,float]](nd)
    t = newTransporters(g, g[0], 1) 
  threads:
    var peot = newseq[array[2,float]](nd)
    for mu in 1..<nd:
      for nu in 0..<mu:
        discard t[mu]^*g[nu]
        discard t[nu]^*g[mu]
        threadBarrier()
        for i in g[mu]:
          pl[i] := redot(t[mu].field[i], t[nu].field[i])
        threadBarrier()
        for site in pl.sites:
          let ps = pl{site}
          peot[mu][lo.coords[mu][site] mod 2] += ps
          peot[nu][lo.coords[nu][site] mod 2] += ps
        threadBarrier()
    peot.threadRankSum
    threadSingle:
      for dir in 0..<nd:
        peo[dir][0] += peot[dir][0]
        peo[dir][1] += peot[dir][1]
  let n = 1.0 / (lo.physVol.float*0.5*float((nd-1)*nc))
  for dir in 0..<nd:
    peo[dir][0] *= n
    peo[dir][1] *= n
  peo

when isMainModule:
  import qex, gauge, physics/qcdTypes
  import os, sequtils

  qexInit()

  let
    fn = if paramCount() > 0: paramStr 1 else: ""
    lat = if fn.len == 0: @[8,8,8,8] else: fn.getFileLattice
    lo = lat.newLayout
  var g = lo.newGauge
  if fn.len == 0:
    g.random
  elif 0 != g.loadGauge fn:
    echo "ERROR: couldn't load gauge file: ",fn
    qexFinalize()
    quit(-1)
  echo "plaq: ",g.plaq
  let plaqeo = g.plaqEvenOddS4
  echo "plaq_eo: ",plaqeo

  qexFinalize()
