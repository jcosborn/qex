import qex # QEX
import physics/qcdTypes # Types
import strutils, times # OS operations, string manip., timing
import xy_cluster_init_and_io # For IO
import base/hyper # "hyper"
import streams # For reading/writing
import arraymancer # For doing gradient flow with vectorized operations

#[ Initialize QEX ]#

# Initialize QEX
qexinit()

# Start timer
tic()

#[ Grab command line and XML inputs ]#

# Create command line inputs
let (start_config, end_config, config_space, save_freq,
     xml_file, rank_geom, def_fn, io_path) = read_cmd()

# Read XML file input
let (int_prms, flt_prms, seed_prms, str_prms) = read_xml(xml_file)

#[ Set simuation up ]#

# Set parameters
let
  # Lattice dimensions
  lat = @[int_prms["Nx"], int_prms["Ny"]]

  # Number of configurations
  nconfig = end_config - start_config

  # Set default lattice name
  lat_fn = io_path & def_fn & "."

# Maximum flow time
var max_flt = 0.0

# Check if maximum flow time set to zero
if flt_prms["max_flt"] == 0.0:
   # If so, set to (0.5 * L)^2 / 4.0
   max_flt = round(pow(0.5*lat[0], 2.0) / 4.0, 2)
else:
   # Otherwise, set to value specified
   max_flt = flt_prms["max_flt"]

# Print number of ranks (and current rank)
echo "rank ", myRank, "/", nRanks

# Print number of threads (and current thread)
threads: echo "thread ", threadNum, "/", numThreads

# Set lattice layout
let lo = lat.newLayout(rank_geom)

# Create lattice data structures
var g = lo.Real

# Create serial RNG generator
var R: RngMilc6

# Seeed serial RNG
R.seed(seed_prms["serial_seed"], 987654321)

# Check if starting off of old configuration
if (int_prms["start_config"] == 0) and (start_config == 0):
  # Check start
  if str_prms["start"] == "hot":
    # Cycle through lattice sites
    for i in g.sites:
       # Set lattice site to random angle
       g{i} := PI*(2.0*R.uniform - 1.0)
  elif str_prms["start"] == "cold":
    # Cycle through lattice sites
    for i in g.sites:
       # Set angle at this site to zero
       g{i} := 0.0
else:
  # Define filename
  let in_fn = lat_fn & intToStr(int_prms["start_config"]) & ".lat"

  # Create reader
  var reader = g.l.newReader(in_fn)

  # Read new spin field in
  reader.read(g)

  # Close reader
  reader.close()

#[ Define procs for Monte Carlo & measurements ]#

# Get (x,y) coordinate
proc get_cart_coord(site: int): seq =
   # Create variable for coordinate
   var coord = newseq[int](2)

   # Grab cartesian coordinate
   coord.lexCoord(site, lat)

   # Return (x,y) coordinate
   result = coord

# Get neighbor lex coordinate
proc get_nghbr_lex_coord(site, mu, dir: int): int =
   # Get (x,y) coordinate of "site"
   var coord = get_cart_coord(site)

   # Define (x,y) coord of neighbor
   var nghbr_coord = coord

   # Update "site" coord to neighbor coord
   nghbr_coord[mu] = coord[mu] + dir

   # Enforce boundary conditions
   if (nghbr_coord[mu] < 0):
      # Wrap around lattice
      nghbr_coord[mu] = lat[mu] - 1
   elif (nghbr_coord[mu] == lat[mu]):
      # Wrap around lattice
      nghbr_coord[mu] = 0

   # Return lex coordinate of neighbor
   result = lexIndex(nghbr_coord, lat)

# Calculate local energy
proc meas_loc_energ(g: any, site: int): float =
  # Initialize local energy
  var loc_energ = 0.0

  # Go through positive x and y directions
  for mu in 0..<g.l.nDim:
     # Get lexicographical coordinate
     let nghbr_site = get_nghbr_lex_coord(site, mu, 1)

     # Calculate energy
     loc_energ = loc_energ - cos(g{site}[][] - g{nghbr_site}[][]) + 1.0

  # Return result
  result = loc_energ

# Measure energy
proc meas_energ(g: any): float =
  # Initialize local energy
  var energ = 0.0
  
  # Cycle through sites
  for site in g.sites:
     # Calculate local energy
     energ = energ + g.meas_loc_energ(site)
  
  # Return energy
  result = energ / lat[0] / lat[1]

# Measure magnetization
proc meas_mag(g: any): seq =
  # Initialize magnetization
  var mag = @[0.0, 0.0]

  # Go through each site
  for site in g.sites:
     # Increment x magnetization
     mag[0] = mag[0] + cos(g{site}[][]) / lat[0] / lat[1]

     # Increment y magnetization
     mag[1] = mag[1] + sin(g{site}[][]) / lat[0] / lat[1]

  # Return magnetization
  result = mag

# For getting random sites
proc rand_site(): int =
   # Return random lattice lex index for lattice site
   result = int round((lat[0]*lat[1] - 0.5)*R.uniform() - 0.5)

#[ For growing cluster ]#

# Proc for growing cluster
proc grow_cluster(site: int, rangle: float, clstr: var seq[int]) =
   #[ Do calculations at this site and flip it ]#

   # Get dot product
   let st_dot = cos(g{site}[][] - rangle)

   #[ Flip vector ]#

   # Get x coordinate
   let x = cos(g{site}[][]) - 2.0 * st_dot * cos(rangle)
   
   # Get y coordinate
   let y = sin(g{site}[][]) - 2.0 * st_dot * sin(rangle)

   # Flip spin at site
   g{site} := arctan2(y, x)

   # Mark spin in cluster
   clstr.add(site)

   #[ Recursively go through neighbors ]#

   # Cycle through neighbors
   for mu in 0..<g.l.nDim:
      # Cycle through directions
      for dir in [-1, 1]:
         #[ Take care of indexing ]#

         # Get site index
         let nghbr_site = get_nghbr_lex_coord(site, mu, dir)

         #[ Do calculations if neighbor not in cluster ]#

         # Check if site already in cluster
         if nghbr_site notin clstr:
            #[ Calculate change in energy ]#

            # Get dot product
            let nb_dot = cos(g{nghbr_site}[][] - rangle)

            # Get change in energy
            let dH = -2.0 * flt_prms["J"] * st_dot * nb_dot

            #[ Check if lattice site to be added to cluster & do recursion ]#

            # Define probability
            var prob = R.uniform()

            # Check if member to be added
            if 1.0 - exp(min(0.0, dH)) >= prob:
               # Recursively call "grow cluster" for this neighbor
               grow_cluster(nghbr_site, rangle, clstr)

#[ For doing cluster update ]#
proc evolve_cluster(): float =
  #[ Small initial setup ]#

  # Initialize cluster updates
  var clstr = newSeq[int]()

  # Get random lattice site
  let st = rand_site()

  #[ Grow cluster ]#

  let rangle = 2.0*PI*R.uniform() - PI

  # Grow cluster
  grow_cluster(st, rangle, clstr)

  # Return size of cluster
  result = float(clstr.len)

#[ Do gradient flow ]#

# Complex GF proc
proc gflow() = 
   #[ Create appropriate data structures ]#

   # Initialize array of lattice variables
   var
      # Initialize complex spins, "momentum", and Z = X
      z = newTensor[type(complex(0.0, 0.0))]([lo.nSites])
      p = newTensor[type(complex(0.0, 0.0))]([lo.nSites])
      Z = newTensor[type(complex(0.0, 0.0))]([lo.nSites])

      # Forward shifter indices
      fi = newTensor[int]([lo.nDim, lo.nSites])

      # Backward shifter indices
      bi = newTensor[int]([lo.nDim, lo.nSites])

   # Set some coefficients
   let
      # Set GF epsilon
      eps = complex(flt_prms["gf_eps"], 0.0)

      # First RK coeff
      c1 = complex(0.25, 0.0)

      # Second RK coeff
      c2 = complex(8.0/9.0, 0.0)

      # Third RK coeff
      c3 = complex(17.0/9.0, 0.0)

      # Fourth RK coeff
      c4 = complex(0.75, 0.0)

   # Cycle through sites
   for site in 0..<lo.nSites:
      # Fill complex z
      z[site] = complex(cos(g{site}[][]), sin(g{site}[][]))

      # Cycle through lattice directions
      for dir in 0..<lo.nDim:
         # Get forward neighbor
         fi[dir, site] = get_nghbr_lex_coord(site, dir, 1)

         # Get backward neighbor
         bi[dir, site] = get_nghbr_lex_coord(site, dir, -1)

   #[ Functions for doing gradient flow ]#
   
   # Fast proc for shift
   proc shift(x: int): Complex[system.float64] =
      # z[x+/-mu]
      result = z[x]

   # Proc for calculating real part
   proc Re(z: Complex[system.float64]): Complex[system.float64] =
      # Real part
      complex(z.re, 0.0)

   # Proc for calculating imaginary part
   proc Im(z: Complex[system.float64]): Complex[system.float64] =
      # Imaginary part (with "i")
      complex(0.0, z.im)

   # Proc for calculating norm
   proc norm2[T](x: Tensor[T]): T =
      # Calculate norm
      result = reduce(z *. z.map(conjugate), `+`)

   # Calculate "X" in dot(z) = X(z) * z
   proc X() =
      # Cycle through directions in lattice
      for mu in 0..<lo.nDim:
         # Start case
         case mu:                     # Forward             # Backward
            of 0: Z[_] =        fi[mu, _].map(shift) + bi[mu, _].map(shift)
            else: Z[_] = Z[_] + fi[mu, _].map(shift) + bi[mu, _].map(shift)

      # Multiply by z^dagger and get i * Im(X) (*. = Hadamard product)
      Z = map(z.map(conjugate) *. Z[_], Im)

   # Regular norm
   proc norm[T](x, y: Tensor[T]): T =
      # Return norm
      result = reduce(x[_] *. y[_], `+`)

   # For doing calculations of gf quantities
   proc gf_calcs(flt: float; itn: int) =
      #[ Calculation of standard gradient flow measurements ]#
      proc do_calc() = 
         # Initialize energy density
         var ce = 0.0

         # Cycle through lattice directions
         for mu in 0..<lo.nDim:
            # Get shifted field
            Z[_] = fi[mu, _].map(shift)

            # Forward direction
            ce = ce - norm(z.map(conjugate), Z).re / lo.nSites + 1.0

         # Flow time multiplied by ce
         let tce = flt * ce

         # Calculate norm deviation
         let ndev = z.norm2().re / lo.nSites - 1.0

         # Print result of calculations
         echo "flt, E, tE, dev: ", round(flt, 2), " ", ce, " ", tce, " ", ndev

      #[ Do calculations ]#

      # Do calculation
      do_calc()

   proc integrate(flt: var float; itn: var int) = 
      #[ Initial calculations for previous iteration ]#

      # Calculations
      gf_calcs(flt, itn)

      #[ Do integration ]#

      # "Z0" and "W1"
      X(); p = c1 * Z; z = map(eps * p, exp) *. z[_];

      # "Z1" and "W2"
      X(); p = c2 * Z - c3 * p; z = map(eps * p, exp) *. z[_];

      # "Z2" and "W3"
      X(); p = c4 * Z  - p; z = map(eps * p, exp) *. z[_];

      # Increment flow time
      flt = flt + flt_prms["gf_eps"]

      # Increment iteration
      itn = itn + 1

   #[ Do gradient flow integration ]#

   # Define flow time
   var flow_time = 0.0

   # Define iteration
   var itn = 0

   # Cycle through flow times
   while flow_time <= max_flt:
      # Do gradient flow integration
      integrate(flow_time, itn)

#[ Perform Monte Carlo w/ cluster update ]#

# Make initial measurements
let
  # Get magnetization
  magi = g.meas_mag()

  # Get squared mag
  tmagi = magi[0]*magi[0] + magi[1]*magi[1]

  # Get energy
  energi = g.meas_energ()

# Print out initial values
echo "initial H, Mx, My, M^2: ", energi, " ", magi[0], " ", magi[1], " ", tmagi, "\n"

# Cycle through configurations
for config in 1..<nconfig + 1:
   #[ Do cluster update ]#

   # Print out configuration number
   echo "Configuration ", config

   # Average cluster size
   var avg_size = 0.0

   # Cycle through sweeps
   for sweep in 1..<config_space: 
      # Evolve cluster
      let size = evolve_cluster()

      # Calculate average cluster size
      avg_size = avg_size + size / config_space

   #[ Measurements and IO ]#

   # Check if measurements to be performed
   if (int_prms["meas_freq"] > 0):
      # Check if measurement to be made
      if (0 == config mod int_prms["meas_freq"]):
         #[ Energy information ]#

         # Get magnetization
         let mag = g.meas_mag()

         # Get the total magnetization
         let tmag = mag[0]*mag[0] + mag[1]*mag[1]

         # Initialize energy
         let energ = g.meas_energ()

         # Print value of energy density
         echo "H, Mx, My, M^2, <size>: ",energ," ",mag[0]," ",mag[1]," ",tmag," ",avg_size

   # Check if gradient flow to be performed
   if (int_prms["gf_freq"] > 0):
      # Check if gauge flow to be performed
      if (0 == config mod int_prms["gf_freq"]):
         #[ Do flow ]#

         # Evolve with gradient flow
         gflow()

   # Check if configuration to be saved
   if (save_freq > 0):
      # Check if config to be saved
      if (0 == config mod save_freq):
          # Filename
          let fn = lat_fn & intToStr(config) & ".lat"

          # Create new file
          var file = newFileStream(fn, fmWrite)

          # Make sure the configuration is going to be able to be saved
          if not file.isNil:
             # Save spin field
             file.write g

             # Tell user what you did
             echo "Wrote " & fn
          else:
             # Tell user that configuration did not save successfully
             quit("Was not able to write " & fn)

          # Flush file
          file.flush

   # Create space in output
   echo "\n"

#[ Finalize QEX ]#

# End timer
toc("done")

# Finalize QEX
qexfinalize()
