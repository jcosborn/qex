#[ Imports ]#

import qex # QEX
import physics/qcdTypes # Types
import streams # For reading/writing
import arraymancer # For doing gradient flow with vectorized operations
import xy_cluster_init_and_io # For IO
import parseopt # For parsing command line arguments
import tables # For organizing data
import streams, parsexml, strutils # For parsing XML
import system # For system-specific operations
import times # Timing for generating seeds
import os # For checking if files exist
import base/hyper # For translating between Cart. & lex. coords.

#[

   This is just a copy of the U(1) XY model gradient flow code that is already
   in xy_cluster.nim; however, this code can be used to just to run gradient flow.

]#

#[ Read command line inputs ]#
proc read_inputs*(): auto =
   #[ Initialize variables ]#

   # Set variables
   var
      # Starting configuration
      start_config = 0

      # Ending configuration
      end_config = 0

      # XML file name
      xml_file = "input.xml"

      # Structure for all file names
      def_fn = "checkpoint"

      # Path for all IO
      io_path = "./"

   #[ Parse inputs and return ]#

   # Print what you're doing
   echo "\n ~~~~ Command line options ~~~~\n"

   # Create parser for command line options
   var cm_opts = initOptParser()

   # Cycle through command line arguments
   while true:
      # Go to next option
      cm_opts.next()

      # Start case
      case cm_opts.kind
         of cmdEnd: break # Exit of options
         of cmdShortOption, cmdLongOption, cmdArgument:
            # Check if starting config.
            if cm_opts.key == "start_config":

               # Set ending config.
               start_config = parseInt(cm_opts.val)

               # Print ending config.
               echo "start config: " & cm_opts.val

            # Check if ending config.
            if cm_opts.key == "end_config":
               # Set ending config.
               end_config = parseInt(cm_opts.val)

               # Print ending config.
               echo "end config: " & cm_opts.val


            # Check if xml file
            if cm_opts.key == "xml":
               # Set xml file
               xml_file = cm_opts.val

               # Tell user where information is being read from
               echo "XML file: " & cm_opts.val

            # Check if user wants a different filename
            if cm_opts.key == "filename":
               # Set filename
               def_fn = cm_opts.val

            # Check if user wants a different path
            if cm_opts.key == "path":
               # Set path
               io_path = cm_opts.val 

   # Return parsed options as tuple
   result = (start_config, end_config, xml_file, def_fn, io_path)

#[ Run U(1) gradient flow ]#

# Initialize QEX
qexinit()

#[ Grab command line and XML inputs ]#

# Create command line inputs
let (start_config, end_config, xml_file, def_fn, io_path) = read_inputs()

# Read XML file input
let (int_prms, flt_prms, seed_prms, str_prms) = read_xml(xml_file)

#[ Set lattice up ]#

# Set parameters
let
  # Lattice dimensions
  lat = @[int_prms["Nx"], int_prms["Ny"]]

  # Set default lattice name
  lat_fn = io_path & def_fn & "."

# Set lattice layout
let lo = lat.newLayout(@[1, 1])

# Create lattice data structures
var g = lo.Real

#[ Set gradient flow up ]#

# Maximum flow time
var max_flt = 0.0

# Check if maximum flow time set to zero
if flt_prms["max_flt"] == 0.0:
   # If so, set to (0.5 * L)^2 / 4.0
   max_flt = round(pow(0.5*lat[0], 2.0) / 4.0, 2)
else:
   # Otherwise, set to value specified
   max_flt = flt_prms["max_flt"]

#[ Gradient flow procs ]#

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

#[ Do gradient flow ]#

# Cycle through configurations
for config in start_config..<end_config + 1:
   #[ Load configuration ]#
   
   # Define file name
   let fn = lat_fn & intToStr(config) & ".lat"

   # Open file
   var file = newFileStream(fn, fmRead)

   # Check if able to read
   if file == nil:
      # Quit program
      quit("Was not able to read " & fn & ". Exiting.")
   else:
      # Cycle through sites
      for site in 0..<lo.nSites:
         # Temporary variable
         var site_val: typeof(g{site}[][])

         # Read new spin field in
         discard file.readData(site_val.addr, site_val.sizeof)

         # Save to g
         g{site} := site_val


   #[ Do gradient flow ]#

   # Gradient flow measurement
   gflow()

# Finalize QEX
qexfinalize()