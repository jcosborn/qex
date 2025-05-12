#[ ~~~~ Imports ~~~~ ]#

import qex # QEX
import gauge # For gauge measurements
import sequtils, layout, field # For some operators w/ gauge field
import physics/qcdTypes # For generic QCD types
import base/hyper

#[ ~~~~ Various measurement functions ~~~~ ]#

#[ For measuring plaquette ]#
proc mplaq*(g: auto) =
   # Calculate plaquette
   let
      pl = g.plaq
      nl = pl.len div 2
      ps = pl[0..<nl].sum * 2.0
      pt = pl[nl..^1].sum * 2.0

   #echo pl[0..<nl]
   #echo pl[nl..^1]

   # Print information about plaquette
   echo "MEASplaq ss: ",ps,"  st: ",pt,"  tot: ",0.5*(ps+pt)

#[ For measuring gauge S4 order parameter ]#
proc s4_gauge*(g: auto) = 
  #[ Calculates "S4" (pure gauge) order parameter from arXiv:1111.2317 
     Written by Xiaoyong Jin ]#
   
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

  for dir in 0..<nd:
   echo "MEASplaq ", dir, "-dir even/odd: ", peo[dir][0], " ", peo[dir][1]

#[ For measuring polyakov loop ]#
proc ploop*(g: auto) =
   # Calculate Polyakov loop
   let pg = g[0].l.physGeom
   var pl = newseq[typeof(g.wline @[1])](pg.len)
   for i in 0..<pg.len:
      pl[i] = g.wline repeat(i+1, pg[i])
   let
      pls = pl[0..^2].sum / float(pl.len-1)
      plt = pl[^1]

   # Print information about Polyakov loop
   echo "MEASploop spatial: ",pls.re," ",pls.im," temporal: ",plt.re," ",plt.im

#[ Testing ]#
when isMainModule:
   #[ Imports ]#
   import qex
   import physics/qcdTypes

   #[ Start QEX ]#

   # Start QEX
   qexInit()

   # Setup gauge field
   let defaultLat = @[8, 8, 8, 8]
   defaultSetup()
   g.random()

   #[ Gauge field w/o single-site symmetry ]#

   # Print result of s4 order parameter
   g.s4_gauge()
   g.mplaq()

   qexFinalize()