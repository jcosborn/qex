#[ Imports from the original gaugeAction.nim code ]#
import base
import maths
import layout
import layout/shifts
import gauge/gaugeUtils
import gauge/staples
import gauge/gaugeAction

#[

This code modifies the forces in <>/src/gauge/gaugeAction.nim written
by James Osborn and Xiaoyong Jin so that smearing can be applied to the gauge force

]#

proc gaugeForceCust*[T](c: GaugeActionCoeffs; uu: openArray[T]; f: array|seq) =
  mixin load1, adj
  tic("gaugeForce")
  let u = cast[ptr cArray[T]](unsafeAddr(uu[0]))
  let lo = u[0].l
  let nd = lo.nDim
  #let np = (nd*(nd-1)) div 2
  let nc = u[0][0].ncols
  let cp = c.plaq / float(nc)
  let cr = c.rect / float(nc)
  var cs = startCornerShifts(uu)
  var ru:FieldArray[type(u[0]).V,type(u[0]).T]  # the rect parts of 3
  var sb:seq[seq[ShiftB[type(u[0][0])]]]  # backward ru
  var sf:seq[seq[ShiftB[type(u[0][0])]]]  # forward stf
  if cr!=0:
    ru = newFieldArray2(lo,type(u[0]),[nd,nd],mu!=nu)
    sb.newseq(nd)
    for mu in 0..<nd:
      sb[mu].newseq(nd)
      for nu in 0..<nd:
        if mu==nu: continue
        sb[mu][nu].initShiftB(ru[mu,nu], nu, -1, "all")
    sf.newseq(nd)
    for mu in 0..<nd:
      sf[mu].newseq(nd)
      for nu in 0..<nd:
        if mu==nu: continue
        sf[mu][nu].initShiftB(u[mu], nu, 1, "all")
  toc("gaugeForce init")
  var (stf,stu,ss) = makeStaples(uu, cs)
  toc("gaugeForce makeStaples")
  threads:
    tic()
    if cr!=0:
      for mu in 1..<nd:
        for nu in 0..<mu:
          sf[mu][nu].startSB(stf[mu,nu][ix])
          sf[nu][mu].startSB(stf[nu,mu][ix])
    for mu in 0..<nd:
      f[mu] := 0
      if cr!=0:
        for nu in 0..<nd:
          if mu!=nu:
            ru[mu,nu] := 0
    for ir in u[0]:
      for mu in 1..<nd:
        for nu in 0..<mu:
          # plaq
          f[mu][ir] += cp * stf[mu,nu][ir]
          f[nu][ir] += cp * stf[nu,mu][ir]
          if isLocal(ss[mu][nu],ir):
            var bmu: type(load1(u[0][0]))
            localSB(ss[mu][nu], ir, assign(bmu,it), stu[mu,nu][ix])
            f[mu][ir] += cp * bmu
            if cr!=0:
              var umu,unu,bmunu: type(load1(u[0][0]))
              getSB(cs[nu][mu], ir, assign(unu,it), u[nu][ix])
              getSB(cs[mu][nu], ir, assign(umu,it), u[mu][ix])
              bmunu := bmu * unu
              f[nu][ir] += cr * bmunu * umu.adj
              ru[nu,mu][ir] += bmu.adj * u[nu][ir] * umu
              ru[mu,nu][ir] += u[nu][ir].adj * bmunu
          if isLocal(ss[nu][mu],ir):
            var bnu: type(load1(u[0][0]))
            localSB(ss[nu][mu], ir, assign(bnu,it), stu[nu,mu][ix])
            f[nu][ir] += cp * bnu
            if cr!=0:
              var unu,umu,bnumu: type(load1(u[0][0]))
              getSB(cs[mu][nu], ir, assign(umu,it), u[mu][ix])
              getSB(cs[nu][mu], ir, assign(unu,it), u[nu][ix])
              bnumu := bnu * umu
              f[mu][ir] += cr * bnumu * unu.adj
              ru[mu,nu][ir] += bnu.adj * u[mu][ir] * unu
              ru[nu,mu][ir] += u[mu][ir].adj * bnumu
          if cr!=0:
            if isLocal(sf[mu][nu],ir):
              var smu,unu,smunu: type(load1(u[0][0]))
              localSB(sf[mu][nu], ir, assign(smu,it), stf[mu,nu][ix])
              getSB(cs[nu][mu], ir, assign(unu,it), u[nu][ix])
              smunu := smu * unu.adj
              f[mu][ir] += cr * u[nu][ir] * smunu
              f[nu][ir] += cr * u[mu][ir] * smunu.adj
              ru[nu,mu][ir] += u[mu][ir].adj * u[nu][ir] * smu
            if isLocal(sf[nu][mu],ir):
              var snu,umu,snumu: type(load1(u[0][0]))
              localSB(sf[nu][mu], ir, assign(snu,it), stf[nu,mu][ix])
              getSB(cs[mu][nu], ir, assign(umu,it), u[mu][ix])
              snumu := snu * umu.adj
              f[nu][ir] += cr * u[mu][ir] * snumu
              f[mu][ir] += cr * u[nu][ir] * snumu.adj
              ru[mu,nu][ir] += u[nu][ir].adj * u[mu][ir] * snu
    toc("gaugeForce local")
    for mu in 1..<nd:
      for nu in 0..<mu:
        var needBoundary = false
        boundaryWaitSB(ss[mu][nu]): needBoundary = true
        boundaryWaitSB(ss[nu][mu]): needBoundary = true
        if needBoundary:
          boundarySyncSB()
          for ir in lo:
            if not isLocal(ss[mu][nu],ir):
              var bmu: type(load1(u[0][0]))
              getSB(ss[mu][nu], ir, assign(bmu,it), stu[mu,nu][ix])
              f[mu][ir] += cp * bmu
              if cr!=0:
                var umu,unu,bmunu: type(load1(u[0][0]))
                getSB(cs[nu][mu], ir, assign(unu,it), u[nu][ix])
                getSB(cs[mu][nu], ir, assign(umu,it), u[mu][ix])
                bmunu := bmu * unu
                f[nu][ir] += cr * bmunu * umu.adj
                ru[nu,mu][ir] += bmu.adj * u[nu][ir] * umu
                ru[mu,nu][ir] += u[nu][ir].adj * bmunu
            if not isLocal(ss[nu][mu],ir):
              var bnu: type(load1(u[0][0]))
              getSB(ss[nu][mu], ir, assign(bnu,it), stu[nu,mu][ix])
              f[nu][ir] += cp * bnu
              if cr!=0:
                var unu,umu,bnumu: type(load1(u[0][0]))
                getSB(cs[mu][nu], ir, assign(umu,it), u[mu][ix])
                getSB(cs[nu][mu], ir, assign(unu,it), u[nu][ix])
                bnumu := bnu * umu
                f[mu][ir] += cr * bnumu * unu.adj
                ru[mu,nu][ir] += bnu.adj * u[mu][ir] * unu
                ru[nu,mu][ir] += u[mu][ir].adj * bnumu
    if cr!=0:
      for mu in 1..<nd:
        for nu in 0..<mu:
          var needBoundary = false
          boundaryWaitSB(sf[mu][nu]): needBoundary = true
          boundaryWaitSB(sf[nu][mu]): needBoundary = true
          if needBoundary:
            boundarySyncSB()
            for ir in lo:
              if not isLocal(sf[mu][nu],ir):
                var smu,unu,smunu: type(load1(u[0][0]))
                getSB(sf[mu][nu], ir, assign(smu,it), stf[mu,nu][ix])
                getSB(cs[nu][mu], ir, assign(unu,it), u[nu][ix])
                smunu := smu * unu.adj
                f[mu][ir] += cr * u[nu][ir] * smunu
                f[nu][ir] += cr * u[mu][ir] * smunu.adj
                ru[nu,mu][ir] += u[mu][ir].adj * u[nu][ir] * smu
              if not isLocal(sf[nu][mu],ir):
                var snu,umu,snumu: type(load1(u[0][0]))
                getSB(sf[nu][mu], ir, assign(snu,it), stf[nu,mu][ix])
                getSB(cs[mu][nu], ir, assign(umu,it), u[mu][ix])
                snumu := snu * umu.adj
                f[nu][ir] += cr * u[mu][ir] * snumu
                f[mu][ir] += cr * u[nu][ir] * snumu.adj
                ru[mu,nu][ir] += u[nu][ir].adj * u[mu][ir] * snu
          threadBarrier()
          sb[mu][nu].startSB(ru[mu,nu][ix])
          sb[nu][mu].startSB(ru[nu,mu][ix])
      toc("gaugeForce staple boundary")
      for ir in u[0]:
        for mu in 1..<nd:
          for nu in 0..<mu:
            if isLocal(sb[mu][nu],ir):
              var b: type(load1(u[0][0]))
              localSB(sb[mu][nu], ir, assign(b,it), ru[mu,nu][ix])
              f[mu][ir] += cr * b
            if isLocal(sb[nu][mu],ir):
              var b: type(load1(u[0][0]))
              localSB(sb[nu][mu], ir, assign(b,it), ru[nu,mu][ix])
              f[nu][ir] += cr * b
      toc("gaugeForce back rect local")
      for mu in 1..<nd:
        for nu in 0..<mu:
          var needBoundary = false
          boundaryWaitSB(sb[mu][nu]): needBoundary = true
          boundaryWaitSB(sb[nu][mu]): needBoundary = true
          if needBoundary:
            boundarySyncSB()
            for ir in lo:
              if not isLocal(sb[mu][nu],ir):
                var b: type(load1(u[0][0]))
                getSB(sb[mu][nu], ir, assign(b,it), ru[mu,nu][ix])
                f[mu][ir] += cr * b
              if not isLocal(sb[nu][mu],ir):
                var b: type(load1(u[0][0]))
                getSB(sb[nu][mu], ir, assign(b,it), ru[nu,mu][ix])
                f[nu][ir] += cr * b
  toc("gaugeForce end")

proc forceACust*(c: GaugeActionCoeffs; g,f: auto) =
  ## Specialized gauge force for plaq + adjplaq
  mixin load1, adj
  tic("forceA")
  let lo = g[0].l
  let nd = lo.nDim
  let nc = g[0][0].ncols
  let cp = c.plaq / float(nc)
  let ca = 2.0 * c.adjplaq / float(nc*nc)
  var cs = startCornerShifts(g)
  toc("gaugeForce startCornerShifts")
  var (stf,stu,ss) = makeStaples(g, cs)
  toc("gaugeForce makeStaples")
  for i in 0..<nd:
    f[i] := 0
  threads:
    tic()
    for ir in g[0]:
      for mu in 1..<nd:
        for nu in 0..<mu:
          let tmn = dot(stf[mu,nu][ir], g[mu][ir])
          f[mu][ir] += (cp+ca*tmn) * stf[mu,nu][ir]
          let tnm = dot(stf[nu,mu][ir], g[nu][ir])
          f[nu][ir] += (cp+ca*tnm) * stf[nu,mu][ir]
          if isLocal(ss[mu][nu],ir):
            var bmu: type(load1(g[0][0]))
            localSB(ss[mu][nu], ir, assign(bmu,it), stu[mu,nu][ix])
            let tmu = dot(bmu, g[mu][ir])
            f[mu][ir] += (cp+ca*tmu) * bmu
          if isLocal(ss[nu][mu],ir):
            var bnu: type(load1(g[0][0]))
            localSB(ss[nu][mu], ir, assign(bnu,it), stu[nu,mu][ix])
            let tnu = dot(bnu, g[nu][ir])
            f[nu][ir] += (cp+ca*tnu) * bnu
    toc("gaugeForce local")
    for mu in 1..<nd:
      for nu in 0..<mu:
        var needBoundary = false
        boundaryWaitSB(ss[mu][nu]): needBoundary = true
        boundaryWaitSB(ss[nu][mu]): needBoundary = true
        if needBoundary:
          boundarySyncSB()
          for ir in lo:
            if not isLocal(ss[mu][nu],ir):
              var bmu: type(load1(g[0][0]))
              getSB(ss[mu][nu], ir, assign(bmu,it), stu[mu,nu][ix])
              let tmu = dot(bmu, g[mu][ir])
              f[mu][ir] += (cp+ca*tmu) * bmu
            if not isLocal(ss[nu][mu],ir):
              var bnu: type(load1(g[0][0]))
              getSB(ss[nu][mu], ir, assign(bnu,it), stu[nu,mu][ix])
              let tnu = dot(bnu, g[nu][ir])
              f[nu][ir] += (cp+ca*tnu) * bnu
  toc("gaugeForce end")

#[ Project to traceless/anti-Hermitian ]#
proc projTAH*(f, g: auto; adj = "no_adj") = 
   #[ Projects to traceless/anti-Hermitian component. Silly case statement
      takes care of appropriate adjoint ordering for gauge and matter
      field forces as they are currently coded up ]#

   # Start case
   case adj: # Appropriate for matter fields
      of "no_adj":
         # Start thread block
         threads:
            # Cycle through direction
            for mu in 0..<f.len:
               # Cycle through lattice sites
               for i in f[mu]:
                  # Create useful variable for product
                  var s {.noinit.}: typeof(f[0][0])

                  # Calculate product with gauge field link
                  s := f[mu][i] * g[mu][i].adj

                  # Project to traceless/anti-Hermitian component
                  projectTAH(f[mu][i], s)
      of "adj": # Appropriate for gauge fields
         # Start thread block
         threads:
            # Cycle through direction
            for mu in 0..<f.len:
               # Cycle through lattice sites
               for i in f[mu]:
                  # Create useful variable for product
                  var s {.noinit.}: typeof(f[0][0])

                  # Calculate product with gauge field link
                  s := g[mu][i] * f[mu][i].adj

                  # Project to traceless/anti-Hermitian component
                  projectTAH(f[mu][i], s)