import qex # QEX
import physics/qcdTypes # Gauge & types
import strutils, times # OS operations, string manip., timing
import xy_cluster_init_and_io # For IO
import base/hyper # "hyper"
import streams # For reading/writing

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
    # Start thread block
    threads:
       # Cycle through lattice sites
       for i in g.sites:
          # Set lattice site to random angle
          g{i} := PI*(2.0*R.uniform - 1.0)
  elif str_prms["start"] == "cold":
    # Start thread block
    threads:
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

proc vector(angle: float): seq = 
   result = @[cos(angle), sin(angle)]

proc dot(v1, v2: auto): float =
   result = v1[0] * v2[0] + v1[1] * v2[1]

# Calculate local energy
proc meas_loc_energ(g: any, site: int): float =
  # Initialize local energy
  var loc_energ = 0.0

  # Get vector at site
  let site_vec = vector(g{site}[][])
  
  # Initialize coordinate
  var coord = newseq[int](2)

  # Get coordinate
  coord.lexCoord(site, lat)

  # Go through positive x and y directions
  for nu in 0..<g.l.nDim:
     # Get neighbor coordinate
     var nghbr_coord = [coord[0], coord[1]]

     # Update neighbor coordinate
     nghbr_coord[nu] = nghbr_coord[nu] + 1

     # Get lexicographical coordinate
     let lex_coord = lexIndex(nghbr_coord, lat)

     # Calculate energy
     loc_energ = loc_energ + dot(vector(g{lex_coord}[][]), site_vec)

  # Return result
  result = -flt_prms["J"] * loc_energ

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
   #[ Get vector coordinates and angles ]#

   # Initialize coordinates
   var coord = newseq[int](2)

   # Get coordinate of this site
   coord.lexCoord(site, lat)

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
   for nu in 0..<g.l.nDim:
      # Cycle through directions
      for dir in [-1, 1]:
         #[ Take care of indexing ]#

         # Define neighbor coordinate
         var nghbr_coord = [coord[0], coord[1]]

         # Increment neighbor coordinate
         nghbr_coord[nu] = nghbr_coord[nu] + dir

         # Enforce boundary condition
         if (nghbr_coord[nu] < 0):
            # Wrap around lattice
            nghbr_coord[nu] = lat[nu] - 1
         elif  (nghbr_coord[nu] == lat[nu]):
            # Wrap around lattice
            nghbr_coord[nu] = 0

         # Get site index
         let nghbr_site = lexIndex(nghbr_coord, lat)

         # Check if site already in cluster
         if nghbr_site notin clstr:
            #[ Calculate change in energy ]#

            # Get dot product
            let nb_dot = cos(g{nghbr_site}[][] - rangle)

            # Get change in energy
            let dH = 2.0 * flt_prms["J"] * st_dot * nb_dot

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

  # Get norm
  #let norm = vector(2.0*PI*R.uniform() - PI)

  #[ Grow cluster ]#

  let rangle = 2.0*PI*R.uniform() - PI

  # Grow cluster
  grow_cluster(st, rangle, clstr)

  # Return size of cluster
  result = float(clstr.len)

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
echo "initial H, Mx, My, M^2: ", energi, " ", magi[0], " ", magi[1], " ", tmagi

# Cycle through configurations
for config in 1..<nconfig:
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
         echo "H, Mx, My, M^2, <size>: ", energ, " ", mag[0], " ", mag[1], " ", tmag, " ", avg_size

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

#[ Finalize QEX ]#

# End timer
toc("done")

# Finalize QEX
qexfinalize()
