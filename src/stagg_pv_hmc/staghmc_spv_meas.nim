#[ ~~~~ Imports ~~~~ ]#

import qex # QEX
import gauge # For gauge measurements
import sequtils # For simple seq operations

#[ ~~~~ Various measurement functions ~~~~ ]#

#[ For measuring plaquette ]#
proc mplaq*(g: auto) =
   # Calculate plaquette
   let
      pl = g.plaq
      nl = pl.len div 2
      ps = pl[0..<nl].sum * 2.0
      pt = pl[nl..^1].sum * 2.0

   # Print information about plaquette
   echo "MEASplaq ss: ",ps,"  st: ",pt,"  tot: ",0.5*(ps+pt)

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