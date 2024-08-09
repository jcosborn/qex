import qex
import sequtils

proc plaquette*(g: auto) =
   let
      pl = g.plaq
      nl = pl.len div 2
      ps = pl[0..<nl].sum * 2.0
      pt = pl[nl..^1].sum * 2.0
   echo "MEASplaq ss: ",ps,"  st: ",pt,"  tot: ",0.5*(ps+pt)

proc polyakov*(g: auto) =
   let pg = g[0].l.physGeom
   var pl = newseq[typeof(g.wline @[1])](pg.len)
   for i in 0..<pg.len:
      pl[i] = g.wline repeat(i+1, pg[i])
   let
      pls = pl[0..^2].sum / float(pl.len-1)
      plt = pl[^1]
   echo "MEASploop spatial: ",pls.re," ",pls.im," temporal: ",plt.re," ",plt.im

proc reunit*(g: auto) =
    threads:
      let d = g.checkSU
      threadBarrier()
      echo "unitary deviation avg: ",d.avg," max: ",d.max
      g.projectSU
      threadBarrier()
      let dd = g.checkSU
      echo "new unitary deviation avg: ",dd.avg," max: ",dd.max

proc updateGauge*(u: auto, p: auto; dtau: float) =
  threads:
    for mu in 0..<u.len:
      for s in u[mu]: u[mu][s] := exp(dtau*p[mu][s])*u[mu][s]

proc setGauge*(u: auto; v: auto) =
   threads:
    for mu in 0..<u.len:
      for s in u[mu]: u[mu][s] := v[mu][s]