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

#[ Custom loop structure for going through lattice ]#

# Define custom object to represent lattice sites
type
   # Lattice sites object
   LatSites* = object
      # Lower bound on sites loop
      low: int

      # Upper bound on sites loop
      high: int

      # Physical geometry of lattice
      lat: seq[int]

      # Lattice coordinate
      coord: seq[int]

#[ Not used in any code, but kept for potential later use
# LatSites object constructor
proc AllLatSites(n_sites: int, lat: seq[int]): LatSites =
   # Return LatSites object
   result = LatSites(low: 0, high: n_sites, lat: lat)

# Iterator to go through lattice sites
iterator items(range: LatSites): int =
   # Create variable for site
   var site = range.low

   # Go through sites
   while site < range.high:
      # Pass lattice site on
      yield site

      # Increment site
      inc site

# Iterator to go over lattice sites & coordinates
iterator pairs(range: LatSites): tuple[site: int, coord: seq[int]] =
   # Cycle through lattice sites
   for site in range:
      # Grab cartesian coordinate
      range.coord.lexCoord(site, range.lat)

      # Return site and its cartesian coordinate
      yield (site, range.coord)
]#

#[ For printing plaquette information out ]#
proc print_plaq[T](plqs: seq[T]) =
   # Define strings for output
   let plq_str = @["MEASEvenplaq: ", "MEASOddplaq: "]

   # Cycle through plaquettes
   for ind in 0..<plqs.len:
      # Define string
      var output = plq_str[ind]

      # Cycle through components
      for comp in plqs[ind]:
         # Add to string
         output = output & $comp & " "

      # Print result of plaquette
      echo output

#[ For measuring gauge S4 order parameter - modified from <>/src/gaugeUtils.nim]#
proc s4_gauge*[T](g: seq[T]) = 
   #[ Calculates plaquette and "S4" order parameter from arXiv:1111.2317 and 
      laid out more explicitly in https://github.com/daschaich/KS_nHYP_FA ]#
   
   #[ Set things up ]#

   # Get appropriate procs
   mixin adj, newTransporters

   # Define immutable variables to be used for calculation
   let
      # Lattice layout
      lo = g[0].l

      # Set nd
      nd = lo.nDim 

      # Set nc
      nc = g[0][0].ncols

      # Calculate normalization
      norm = 4.0 / float(lo.physVol*(nd-1)*nd*nc)

   # Define mutable variables to be used for calculation
   var
      # Plaquette
      pl = lo.ColorMatrix()

      # Trace on plaquette on all lattice sites
      tr = lo.Complex()

      # Plaquettes on even sites
      plaq_e = newseq[float](lo.nDim)

      # Plaquettes on odd sites
      plaq_o = newseq[float](lo.nDim)

      # Shifters
      t = newTransporters(g, g[0], 1)

      # Create lat sites object
      sites = LatSites(low: 0, high: lo.nSites, lat: lo.physGeom, 
                       coord: newseq[int](lo.physGeom.len))

   #[ Calculate S4 order parameter ]#

   # Start thread block
   threads:
      # Cycle through mu < nu
      for mu in 1..<nd:
         for nu in 0..<mu:
            # Temporarily store plaquette and shift by n_mu
            pl := (t[mu]^*g[nu]) * (t[nu]^*g[mu]).adj

            # Lay down thread barrier
            threadBarrier()

            # Cycle plaquette lattice sites
            for site in pl:
               # Calculate trace
               tr[site] := trace(pl[site])

            # Lay down another thread barrier
            threadBarrier()

            # Cycle through lattice sites again
            for site in tr.sites:
               # Get cartesian coordinate
               sites.coord.lexCoord(site, sites.lat)

               # Cycle through directions
               for dir in 0..<sites.lat.len:
                  # Check if even sites to be increments
                  if ((mu == dir) or (nu == dir)) and (sites.coord[dir] mod 2 == 0):
                     # Increment plaquette on even n_mu, n_nu sites
                     plaq_e[dir] += tr{site}.re
                  elif (mu == dir) or (nu == dir):
                     # Increment plaquette on odd n_mu, n_nu sites
                     plaq_o[dir] += tr{site}.re

            # Lay down last thread barrier
            threadBarrier()

   # Cycle through directions
   for dir in 0..<sites.lat.len:
      # Normalize plaquettes
      plaq_e[dir] *= norm; plaq_o[dir] *= norm;

   # Print output
   print_plaq(@[plaq_e, plaq_o])

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