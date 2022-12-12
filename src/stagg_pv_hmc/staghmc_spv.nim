#[ ~~~~ Import necessary modules ~~~~ ]#

import qex # QEX
import gauge # Gauge field
import gauge / hypsmear # Import nHYP smearing
import physics / [qcdTypes, stagSolve] # Staggered fermions
import mdevolve # MD evolution
import times # Timing
import macros # Useful macros
import os # For operating system-specific tasks
import streams, parsexml, strutils # For parsing XML
import parseopt # For parsing command line arguments
import tables # For organizing data
import sequtils # For dealing with sequences
import math # Basic mathematical operations
import streams # For reading MILC6 global RNG file
import std/monotimes # For some timing applications
import gauge/wflow # For intermediate measurements

#[ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Authors: Curtis Taylor Peterson, James Osborn, Xiaoyong Jin

Description:

   Hybrid Monte Carlo with staggered fermions and staggered Pauli
   -Villars bosons with nHYP smearing.

Credits:

   Thank you to James Osborn Xiaoyong Jin for developing QEX and helping
   develop this program.

Running program:
   
   To run the program, make sure that you have the approprate XML file
   and run the binary file as follows.

   <parallel_exec> <location_of_bin_file> --start_config=<first_config> ...
   --end_config=<last_config> --save_freq=<frequency_of_data> ...
   --config_space=<spacing_between_numbered_trajs>
   --xml=<location_of_xml_file> --rank_geom=<layout_of_rank_geometry> ...
   --path=<input/output_path> --filename=<base_name_for_files>

   If "start_config" is 0, then a new gauge configuration will be created
   according to specifications (type of start) in XML file. Otherwise,
   program will attempt to look in "path" for files with the base
   filename "filename" (e.g., "checkpoint") for the starting trajectory.

   If "save_freq" is 0, then gauge/rng fields and global rng object will
   not be saved to disk.

Layout of program:

   The following lists the function of each section of the program.

   ~ Define function for integration algorithm: Defines a function that
     takes in a string to specify the HMC integrator
   
   ~ Initialize info for whole file: Defines various variables that 
     control behavior of program; namely, those that are specified 
     in the input.xml file and as command line options
   
   ~ Parallel information: Prints out number of ranks and number
     of threads

   ~ Functions for timing: Custom functions for silly timing 
     operations
   
   ~ Functions for gauge measurements: Gauge measurements
     + mplaq: Measures plaquette
     + ploop: Measures Polyakov loop

   ~ Functions for IO: Functions for readings/writing gauge fields,
     RNG fields and global RNG object to disk.
   
   ~ Initialize fields and RNG: Grabs command line options, uses
     command line options to grab information from XML file, 
     initializes gauge fields, fermion fields, RNG fields and
     global RNG (for accept/reject). Gauge fields, RNG fields
     and global RNG are either grabbed from disk or initialized
     according to command line options and XML file inputs.
   
   ~ Initialize gauge action and smearing: Sets up variables
     for manipulating smeared Dirac operator, gauge action,
     and nHYP smearing itself

   ~ Initialize CG: Sets up conjugate gradient solvers for
     calculating action, doing force computation, and 
     calculating the chiral condensate

   ~ Define functions for Wilson flow measurements: For simple
     Wilson flow measurements

   ~ Define functions setting boundary conditions: Sets boundary
     conditions according to user specification in XML file
     (e.g., pppa for periodic in space and antiperiodic in time,
      aaaa for fully antiperiodic, etc.)

   ~ Define function for pbp measurement: Define functions to 
     measure the chiral condensate "pbp"

   ~ Define functions for staggered operations: Various functions 
     for dealing with fermions/bosons. Most notably, if the fermion
     fields are massless, only the solver dealing with the
     even-even block of the Dirac operator is used; otherwise,
     the full solver is used to calculate the inverse of the
     (massive) Dirac operator.

   ~ Define functions for momentum and field generation: Generates
     random momenta and pseudofermion/boson fields.

   ~ Define functions for action calculation: Calculates action
     and Hamiltonian (there are some poor name choices here)

   ~ Define functions for force calculation: Various functions 
     for calculating the contribution of the fermions and 
     Pauli-Villars bosons to the force (and smears the force)

   ~ Define functions HMC integration: Various key functions for 
     doing HMC integration with Omelyan integrator that shares
     smeared computations

   ~ Define functions for checks of HMC: Function(s) for check
     of hybrid Monte Carlo
     + rev_check: Checks reversibility of integrator

   ~ Do HMC: Main code dealing with explicit hybrid Monte Carlo.
     Cycles through and performs trajectories w/ Metropolis step;
     saves gauge/RNG fields and global RNG object according to 
     specifications in XML file; makes measurements of plaquette,
     Polyakov loop and chiral condensate according to specifications
     in XML file; finally makes appropriate checks, again according
     to specifications in XML file

   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ]#

#[ ~~~~ Define function for integration algorithm ~~~~ ]#

type IntProc = proc(T,V:Integrator; steps:int):Integrator
converter toIntProc(s:string):IntProc =
  template mkProc1(s:untyped):IntProc =
    proc mkInt(T,V:Integrator; steps:int):Integrator {.gensym.} =
      `mk s`(T = T, V = V, steps = steps)
    mkInt
  template mkProc2(s:untyped):IntProc =
    proc mkInt(T,V:Integrator; steps:int):Integrator {.gensym.} =
      `mk s`(T = T, V = V, steps = steps, ss[1].parseFloat)
    mkInt
  template mkProc3(s:untyped):IntProc =
    proc mkInt(T,V:Integrator; steps:int):Integrator {.gensym.} =
      `mk s`(T = T, V = V, steps = steps, ss[1].parseFloat, ss[2].parseFloat)
    mkInt
  template mkProc4(s:untyped):IntProc =
    proc mkInt(T,V:Integrator; steps:int):Integrator {.gensym.} =
      `mk s`(T = T, V = V, steps = steps, ss[1].parseFloat, ss[2].parseFloat, ss[3].parseFloat)
    mkInt
  template mkProc5(s:untyped):IntProc =
    proc mkInt(T,V:Integrator; steps:int):Integrator {.gensym.} =
      `mk s`(T = T, V = V, steps = steps, ss[1].parseFloat, ss[2].parseFloat, ss[3].parseFloat, ss[4].parseFloat)
    mkInt
  let ss = s.split(',')
  # Omelyan's triple star integrators, see Omelyan et. al. (2003)
  case ss[0]:
  of "2MN":
    if ss.len == 1: return mkProc1(Omelyan2MN)
    else: return mkProc2(Omelyan2MN)
  of "4MN5FP":
    if ss.len == 1: return mkProc1(Omelyan4MN5FP)
    elif ss.len == 2: return mkProc2(Omelyan4MN5FP)
    elif ss.len == 3: return mkProc3(Omelyan4MN5FP)
    elif ss.len == 4: return mkProc4(Omelyan4MN5FP)
    elif ss.len == 5: return mkProc5(Omelyan4MN5FP)
    else: return mkProc2(Omelyan4MN5FP)
  of "4MN5FV":
    if ss.len == 1: return mkProc1(Omelyan4MN5FV)
    elif ss.len == 2: return mkProc2(Omelyan4MN5FV)
    elif ss.len == 3: return mkProc3(Omelyan4MN5FV)
    elif ss.len == 4: return mkProc4(Omelyan4MN5FV)
    elif ss.len == 5: return mkProc5(Omelyan4MN5FV)
    else: return mkProc2(Omelyan4MN5FV)
  of "6MN7FV": return mkProc1(Omelyan6MN7FV)
  of "4MN3F1GP":  # lambda = 0.2725431326761773  is  FUEL f3g a0=0.109
    if ss.len == 1: return mkProc1(Omelyan4MN3F1GP)
    else: return mkProc2(Omelyan4MN3F1GP)
  of "4MN4F2GVG": return mkProc1(Omelyan4MN4F2GVG)
  of "4MN4F2GV": return mkProc1(Omelyan4MN4F2GV)
  of "4MN5F1GV": return mkProc1(Omelyan4MN5F1GV)
  of "4MN5F1GP": return mkProc1(Omelyan4MN5F1GP)
  of "4MN5F2GV": return mkProc1(Omelyan4MN5F2GV)
  of "4MN5F2GP": return mkProc1(Omelyan4MN5F2GP)
  of "6MN5F3GP": return mkProc1(Omelyan6MN5F3GP)
  else:
    qexError "Cannot parse integrator: '", s, "'\n",
      """Available integrators (with default parameters):
      2MN,0.1931833275037836
      4MN5FP,0.2750081212332419,-0.1347950099106792,-0.08442961950707149,0.3549000571574260
      4MN5FV,0.2539785108410595,-0.03230286765269967,0.08398315262876693,0.6822365335719091
      6MN7FV
      4MN3F1GP,0.2470939580390842
      4MN4F2GVG
      4MN4F2GV
      4MN5F1GV
      4MN5F1GP
      4MN5F2GV
      4MN5F2GP
      6MN5F3GP"""


#[ ~~~~ Initialize QEX ~~~~ ]#

# Initialize QEX
qexinit()

# Print start date
echo "\nStart: ", now().utc

#[ ~~~~ Initialize info for whole file ~~~~ ]#

# For RNG file header information
const fileMd = "<?xml version=\"1.0\"?>\n<note>generated by QEX</note>\n"

# For RNG file header information
const recordMd = "<?xml version=\"1.0\"?>\n<note>RNG field</note>\n"

# Initialize starting trajectory
var start_config = 0

# Initialize ending trajectory
var end_config = 0

# Initialize spacing between names in saved configurations
var config_space = 1

# Initialize save frequency
var save_freq = 0

# Define integer parameters
var int_prms = {"Ns" : 0, "Nt" : 0, "num_pv" : 0, 
                "Nf" : 0, "a_maxits" : 0, "f_maxits" : 0,
                "pbp_maxits" : 0, 
                "g_steps" : 0, "f_steps" : 0, "pv_steps" : 0,
                "no_metropolis_until" : 0, "num_stoch_srcs" : 0,
                "plaq_freq" : 0, "ploop_freq" : 0, "pbp_freq" : 0, 
                "num_pbp_stoch_srcs" : 0, "approx_order" : 1,
                "rev_check_freq" : 0, "check_solvers" : 0,
                "f_munu_loop" : 0, "wflow_freq" : 0,
                "start_config" : 0}.toTable

# Define float parameters
var flt_prms = {"tau" : 0.0, "beta" : 0.0, "adj_fac" : 0.0, "mass" : 0.0,
                "alpha_1" : 0.0, "alpha_2" : 0.0, "alpha_3" : 0.0, 
                "mass_pv" : 0.0, "a_tol" : 0.0, "f_tol" : 0.0,
                "pbp_tol" : 0.0, "timer_waste_ratio" : 0.05,
                "timer_expand_ratio" : 0.05, "dt" : 0.1,
                "t_max" : 0.0}.toTable

# Define seed parameters
var seed_prms = {"parallel_seed" : intParam("seed", int(1000 * epochTime())).uint64,
                 "serial_seed" : intParam("seed", int(1000 * epochTime())).uint64,
                 "pbp_seed" : intParam("seed", int(1000 * epochTime())).uint64}.toTable

# Define string parameters
var str_prms = {"bc" : "pppa", "start" : "hot", 
                "gauge_int_alg" : "2MN", 
                "ferm_int_alg" : "2MN",
                "pv_int_alg" : "2MN"}.toTable

# Define options
var options = {"verbose_cg_stats" : false, "verbose_timer" : false,
               "show_timers" : false}.toTable

# Initialize attribute name
var attrName = ""

# Initialize attribute value
var attrVal = ""

# Initialize xml file
var xml_file = ""

# Initialize rank geometry
var rank_geom = @[1, 1, 1, 1]

# Initialize volume
var vol = 0

# Initialize default filename
var def_fn = "checkpoint"

# Initialize IO path
var io_path = "./"

# Initialize variable to set update number
var update_num = 0

#[ ~~~~ Parallel information ~~~~ ]#

echo "\n ~~~~ Parallel information  ~~~~\n"

# Print number of ranks
echo "# ranks: ", nRanks

# Print number of threads
threads: echo "# threads: ", numThreads

#[ ~~~~ Functions for extra timing ~~~~ ]#

proc ticc(): float = 
   # Return t0
   result = cpuTime()

proc tocc(message: string, t0: float) =
   # Print timing
   echo message, " ", cpuTime() - t0, "s"

#[ ~~~~ Functions for gauge measurements ~~~~ ]#

#[ For reunitarization ]#
proc reunit(g: auto) =
   # Start timer
   tic()

   # Create separator
   echo ""

   # Get initial time
   let t0 = ticc()

   # Start thread block and reunitarize
   threads:
      let d = g.checkSU
      threadBarrier()
      echo "unitary deviation avg: ",d.avg," max: ",d.max
      g.projectSU
      threadBarrier()
      let dd = g.checkSU
      echo "new unitary deviation avg: ",dd.avg," max: ",dd.max

   # Print timing information
   tocc("Reunitarization:", t0)

   # End timer
   toc("reunit")

#[ For measuring plaquette ]#
proc mplaq(g: auto) =
   # Start timer
   tic()

   # Create separator
   echo ""

   # Get initial time
   let t0 = ticc()

   # Calculate plaquette
   let
      pl = g.plaq
      nl = pl.len div 2
      ps = pl[0..<nl].sum * 2.0
      pt = pl[nl..^1].sum * 2.0

   # Print information about plaquette
   echo "MEASplaq ss: ",ps,"  st: ",pt,"  tot: ",0.5*(ps+pt)

   # Print timing information
   tocc("Plaquette measurement:", t0)

   # End timer
   toc("plaq")

#[ For measuring polyakov loop ]#
proc ploop(g: auto) =
   # Start timer
   tic()

   # Create separator
   echo ""

   # Get initial time
   let t0 = ticc()

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
   
   # Print timing information
   tocc("Polyakov loop measurement:", t0)

   # Create separator
   echo ""

   # End timer
   toc("ploop")

#[ ~~~~ Functions for IO ~~~~ ]#

#[ For reading RngMilc6 ]#
proc read_rng_milc6(filename: string): auto =
   # Get iniial time
   let t0 = ticc()

   # Tell user what you're doing
   echo "loading global rng file: " & filename

   # Create RngMilc6 object
   var rng: RngMilc6

   # Create new file stream
   var file = newFileStream(filename, fmRead)

   # Read file data and put in rng object
   discard file.readData(rng.addr, rng.sizeof)

   # Tell user what you did
   echo "wrote global rng file: " & filename

   # Print timing information
   tocc("Global RNG read:", t0)

   # Return rng
   result = rng

#[ For reading RngMilc6 ]#
proc write_rng_milc6(filename: string, rng: auto) =
   # Get initial time
   let t0 = ticc()

   # Tell user what you're doing
   echo "writing global rng file: " & filename

   # Create new file stream
   var file = newFileStream(filename, fmWrite)

   # Check if nil
   if not file.isNil:
      # Write RngMilc6 object to file
      file.write rng

   # Flush
   file.flush

   # Tell user what you did
   echo "wrote global rng file: " & filename

   # Print timing information
   tocc("Global RNG write:", t0)

#[ For reading in gauge field and RNG ]#
proc initialize_params_fields_and_rngs(): auto =
   #[ ~~~~ Read command line and XML information ~~~~ ]#

   # Get initial time
   let t0 = ticc()

   # Create XML parser
   var x: XmlParser

   # Print out separator
   echo "\n ~~~~ Command line and XML information  ~~~~\n"

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
            # Check if starting trajectory
            if cm_opts.key == "start_config":
               # Set ending trajectory
               start_config = parseInt(cm_opts.val)

               # Print ending trajectory
               echo "start config: " & cm_opts.val

            # Check if ending trajectory
            if cm_opts.key == "end_config":
               # Set ending trajectory
               end_config = parseInt(cm_opts.val)

               # Print ending trajectory
               echo "end config: " & cm_opts.val

            # Check if save frequency
            if cm_opts.key == "save_freq":
               # Set save frequency
               save_freq = parseInt(cm_opts.val)

               # Print save frequency
               echo "save frequency: " & cm_opts.val

            # Check if configuration spacing
            if cm_opts.key == "config_space":
               # Save configuration spacing
               config_space = parseInt(cm_opts.val)

               # Print configuration spacing
               echo "config. spacing: " & cm_opts.val & " trajs."

            # Check if xml file
            if cm_opts.key == "xml":
               # Set xml file
               xml_file = cm_opts.val

               # Tell user where information is being read from
               echo "XML file: " & cm_opts.val

            # Check if MPI layout for physical geometry
            if cm_opts.key == "rank_geom":
               # Split geometry string
               let gm_str_splt = cm_opts.val.split({'.'})

               # Cycle through entries in rank geometry
               for ind in 0..<gm_str_splt.len:
                  # Fill rank geometry string
                  rank_geom[ind] = parseInt(gm_str_splt[ind]) 

            # Check if user wants a different filename
            if cm_opts.key == "filename":
               # Set filename
               def_fn = cm_opts.val

            # Check if user wants a different path
            if cm_opts.key == "path":
               # Set path
               io_path = cm_opts.val   

   # Define xml filename
   var file_stream = newFileStream(xml_file, fmRead)

   # Check if file exists
   if file_stream == nil:
      # If file does not exist, exit
      quit("Cannot open " & xml_file)
   else:
      # Open file
      open(x, file_stream, xml_file)

      # Cycle through data
      while true:
         # Go to next option
         x.next()

         # Start case
         case x.kind
            # Do checks
            of xmlElementStart: # If element name
               # Define name
               attrName = x.elementName
            of xmlCharData: # If element attribute
               # Define value
               attrVal = x.charData

               # Check if attribute of parameter tables
               if int_prms.hasKey(attrName):
                  # Save parameter
                  int_prms[attrName] = parseInt(attrVal)
               elif flt_prms.hasKey(attrName):
                  # Save parameter
                  flt_prms[attrName] = parseFloat(attrVal)
               elif seed_prms.hasKey(attrName):
                  # Save seed
                  seed_prms[attrName] = parseInt(attrVal).uint64
               elif options.hasKey(attrName):
                  # Save option
                  options[attrName] = parseBool(attrVal)
               elif str_prms.hasKey(attrName):
                  # Save string information
                  str_prms[attrName] = attrVal

               # Print variable
               echo attrName & ": " & attrVal
            of xmlEof: break # If end of file, exit
            else: discard # Otherwise, do nothing

      # Close
      x.close()

   # Finally, user tell where gauge fields/rngs will be saved and what name they will go by
   echo "Gauge configs/RNG file names: " & io_path & def_fn & ".<traj #>.{lat, rng, global_rng}"

   #[ ~~~~ Initialize gauge field variables ~~~~ ]#

   # Print out separator
   echo "\n~~~~ Gauge information ~~~~\n"

   # Define various variables for gauge action
   let
      # Define lattice
      lat = intSeqParam("lat", @[int_prms["Ns"], int_prms["Ns"],
                                 int_prms["Ns"], int_prms["Nt"]])

      # Define new lattice layout
      lo = lat.newLayout(rank_geom)

   # Create new gauge
   var
      # Starting gauge field
      g0 = lo.newgauge

      # Gauge field
      g = lo.newgauge

      # Backup gauge field
      gg = lo.newgauge

      # Define gauge momenta
      p = lo.newgauge

      # Gauge force
      f = lo.newgauge

      # Smeared gauge field
      sg = lo.newGauge

      # Backup smeared gauge field
      sgg = lo.newGauge

   # Set physical volume
   vol = lo.physVol

   #[ ~~~~ Load global RNG and RNG field ~~~ ]#

   # Set name for global RNG file
   let global_rng_file = io_path & def_fn & "_" & intToStr(start_config) & ".global_rng"   

   # Set name for RNG field field
   let field_rng_file = io_path & def_fn & "_" & intToStr(start_config) & ".rng"

   # Define new RNG field for pbp
   var r_pbp = lo.newRNGField(RngMilc6, seed_prms["pbp_seed"])   

   # Define new RNG field for fermions
   var r = lo.newRNGField(RngMilc6, seed_prms["parallel_seed"])

   # Create global RNG for HMC
   var R: RngMilc6

   # Seed RNG
   R.seed(seed_prms["serial_seed"], 987654321)

   # Check if RNG is to be loaded
   if start_config != 0:
      # Tell user what you're doing
      echo "\n~~~~ RNG IO ~~~~\n"

      #[ Take care of global RNG ]#
      
      # Read and set global RNG
      R = read_rng_milc6(global_rng_file)

      #[ Take care of RNG field ]#
      # Tell user what you're doing      
      echo "loading rng field: " & field_rng_file

      # Create reader
      var reader = r.l.newReader(field_rng_file)

      # Read RNG file and store info in RNG field
      reader.read(r)

      # Close reader file
      reader.close()

      # Tell user what you did
      echo "read rng field: " & field_rng_file

   #[ ~~~~ Load gauge field ~~~~ ]#

   # Set gauge file name
   let gaugefile = io_path & def_fn & "_" & intToStr(start_config) & ".lat"

   # Check if gauge field is to be loaded
   if start_config != 0:
      # Tell user what you're doing
      echo "\n~~~~ Gauge field IO ~~~~\n"

      # Check if gauge field file exists
      if fileExists(gaugefile):
         # Start timer
         tic("Loading gauge file")

         # Check if gauge file can be read
         if 0 != g.loadGauge gaugefile:
            # If not, throw qex error
            qexError "failed to load gauge file: ", gaugefile

         # Set output to qexLog
         qexLog "loaded gauge from file: ", gaugefile," secs: ", getElapsedTime()

         # End timer
         toc("read")

         # Check if starting config
         if int_prms["start_config"] == start_config:
            # Tell user what you're doing
            echo "\nCurrent config. is same as start config. in XML file. Reunitarizing."

            # Reunitarize
            g.reunit

         # Measure plaquette
         g.mplaq

         # Measure polyakov loop
         g.ploop
   else:
      # Check kind of start that user wants
      if str_prms["start"] == "rand":
         # Do hot start
         g.random r
      elif str_prms["start"] == "unit":
         # Do cold start
         g.unit

   #[ ~~~~ Initialize fermion fields ~~~~ ]#

   # Print out separator
   echo "\n~~~~ Fermion information ~~~~\n"

   # Set up array of fermion and PV fermion masses
   var masses = newSeq[float](int_prms["Nf"] + int_prms["num_pv"])

   # Create fermion field for temporary operations
   var ftmp = lo.ColorVector()

   # Get type of fermion variable
   type f_type = typeof(ftmp)

   # Define fermion fields
   var
      # Define phi field
      phi = newSeq[f_type](masses.len)

      # Define psi field
      psi = newSeq[f_type](masses.len)

   #[ Set fermion mass ]#
   proc set_ferm_mass(ind: int) = 
      # Check if first in series of fermion masses
      if ind == 0:
         # Print separator
         echo "\nStaggered fermion masses:"

      # Set regular fermion mass
      masses[ind] = flt_prms["mass"]

      # Print mass
      echo "Staggered fermion mass: " & flt_prms["mass"].formatFloat(ffDecimal, 4)

   #[ Set boson mass ]#
   proc set_bosn_mass(ind: int) =
      # Check if first in series of PV boson masses
      if ind == int_prms["Nf"]:
         # Print separator
         echo "\nStaggered PV boson masses:"
 
      # Set boson mass
      masses[ind] = flt_prms["mass_pv"]

      # Print boson mass
      echo "Staggered PV boson mass: " & flt_prms["mass_pv"].formatFloat(ffDecimal, 4)

   # Tell user some information about all fields
   echo "All staggered fields in action are half-fields (defined on even sites)"
   echo int_prms["Nf"] * 4, " continuum fermions"
   echo int_prms["num_pv"] * 4, " continuum PV bosons"

   # Add and fermion fields
   for mass_ind in 0..<masses.len:
      # Check through pseudofermion cases
      if mass_ind < int_prms["Nf"]:
         # Set regular fermion mass
         set_ferm_mass(mass_ind)
      else:
         # Set PV boson mass
         set_bosn_mass(mass_ind)

      # Create phi field variable
      phi[mass_ind] = lo.ColorVector()

      # Create psi field variable
      psi[mass_ind] = lo.ColorVector()

   #[ ~~~~ Print timing information ~~~~ ]#

   # Print timing information
   tocc("Initial IO:", t0)

   #[ ~~~~ Return fields and RNG's ~~~~ ]#
   # Return fields
   result = (lo,                      # Layout
             g0, g, sg,               # Gauge fields 
             gg, sgg,                 # Backup gauge fields
             p, f,                    # Momentum and force
             r, r_pbp, R,             # Random number generators
             phi, psi, ftmp, masses)  # Fermion fields/masses

proc save_gaugefield_and_rng(base_filename: string, gauge: auto, 
                             rng_global: auto, rng_field: auto) =
   #[ Save global rng to file ]#

   # Tell user what you're doing
   echo "\n~~~~ RNG IO ~~~~\n"

   # Get initial time
   let t0 = ticc()

   # Call proc to save global rng
   write_rng_milc6(base_filename & &".global_rng", rng_global)

   #[ Save rng field to file ]#
   
   # Tell user what you're doing
   echo "writing rng field: " & base_filename & &".rng"

   # Create writer
   var writer = rng_field.l.newWriter(base_filename & &".rng", fileMd)

   # Write RNG field
   writer.write(rng_field, recordMd)

   # Close writer
   writer.close()

   # Tell user what you did
   echo "wrote rng field: " & base_filename & &".rng"

   #[ Save gauge field ]#

   # Tell user what you're doing
   echo "\n~~~~ Gauge field IO ~~~~\n"

   # Attempt to save gauge configuration
   if 0 != gauge.saveGauge(base_filename & &".lat"):
      # If not able to, throw qex error
      qexError "Failed to save gauge to file: ", base_filename & &".lat"

   # Create separator
   echo "wrote gauge file: " & base_filename & &".lat"

   # Print timing information
   tocc("Gauge/RNG field IO:", t0)

#[ ~~~~ Initialize fields and RNG ~~~~ ]#

# Get fields and RNG
var (lo,               # Layout
     g0, g, sg,        # Gauge fields
     gg, sgg,          # Backup gauge fields
     p, f,             # Momentum and force
     r, r_pbp, R,      # Random number generators
     phi, psi, ftmp,   # Fermion fields
     masses) = initialize_params_fields_and_rngs()

#[ ~~~~ Initialize gauge action and smearing ~~~~ ]#

# Initialize smeared Dirac operator and (unsmeared) gauge action
let
   # For operations involving (smeared) Dirac operator
   stag = newStag(sg)

   # For operations involving backup (smeared) Dirac operator
   stagg = newStag(sgg)

   # Define gauge action coefficients
   gc = GaugeActionCoeffs(plaq: flt_prms["beta"],
                          adjplaq: flt_prms["beta"] * flt_prms["adj_fac"])

# Initialize nHYP smearing parameters
var
   # Get info
   info: PerfInfo

   # Get nHYP coefficients
   coef = HypCoefs(alpha1: flt_prms["alpha_1"], 
                   alpha2: flt_prms["alpha_2"], 
                   alpha3: flt_prms["alpha_3"])

#[ ~~~ Initialize CG ~~~~ ]#

#[ Set up parameters of pbp solver ]#

# Set up force solver
var sppbp = initSolverParams()

# Set force solver tolerance
sppbp.r2req = flt_prms["pbp_tol"]

# Set max number of CG iterations
sppbp.maxits = int_prms["pbp_maxits"]

# Set verbosity of force solver
sppbp.verbosity = 1

#[ Set up parameters for action solver ]#

# Create array of action solves
var spa = newseq[typeof(sppbp)](masses.len)

# Create array of force solves
var spf = newseq[typeof(sppbp)](int_prms["Nf"])

# Cycle through masses
for mass_ind in 0..<masses.len:
   #[ Set up parameters for action solver ]#

   # Set up action solver
   spa[mass_ind] = initSolverParams()

   # Set action solver tolerance
   spa[mass_ind].r2req = flt_prms["a_tol"]

   # Set max number of CG iterations
   spa[mass_ind].maxits = int_prms["a_maxits"]

   # Set verbosite
   spa[mass_ind].verbosity = 1

   # Check mass index
   if mass_ind < int_prms["Nf"]:
      #[ Set up parameters for force solver ]#

      # Set up force solver
      spf[mass_ind] = initSolverParams()

      # Set force solver tolerance
      spf[mass_ind].r2req = flt_prms["f_tol"]

      # Set max number of CG iterations
      spf[mass_ind].maxits = int_prms["f_maxits"]

      # Set verbosity of force solver
      spf[mass_ind].verbosity = 1

#[ ~~~~ For information regarding testing ~~~~]#

# Set time waster ratio
DropWasteTimerRatio = flt_prms["timer_waste_ratio"]

# Set option for CG stats
VerboseGCStats = options["verbose_cg_stats"]

# Set option for verbose timers
VerboseTimer = options["verbose_timer"]

#[ ~~~~ Define functions for Wilson flow measurements ~~~~ ]#

#[ For calculating flowed measurements ]#
proc EQ(gauge: auto, loop: int):auto =
   # Calculate measurements
   let
      # Get loop
      f = gauge.fmunu loop

      # Calculate Yang-Mills density
      (es, et) = f.densityE

      # Calculate topological charge
      q = f.topoQ

   # Return Yang-Mills density and topology
   return (es, et, q)

#[ For flow and flow measurements ]#
proc wflow(gauge: auto) =
   # Start timer
   tic()

   # Get initial time
   let t0 = ticc()

   # Create smeared gauge field
   var flowed_gauge = lo.newgauge

   # Start thread block
   threads:
      # Cycle through gauge field
      for i in 0..<gauge.len:
         # Set starting flowed gauge to gauge field
         flowed_gauge[i] := gauge[i]

   # Initialize measurements
   let (es, et, q) = flowed_gauge.EQ int_prms["f_munu_loop"]

   # Initialize t^2 E(t)
   var t2E = 0.0

   # Initialize parameters for flow
   var ot2E, dt2E, tdt2E: float

   # Print initial measurements
   echo "WFLOW ", 0.0, " ", et, " ", es, " ", q

   # Start flow loop
   flowed_gauge.gaugeFlow(flt_prms["dt"]):
      # Calculate measurement
      let (es, et, q) = flowed_gauge.EQ int_prms["f_munu_loop"]

      # Set original t^2 E(t)
      ot2E = t2E

      # Calculate t^2 E(t)
      t2E = wflowT * wflowT * (es+et)

      # Calculate change in t^2 E(t)
      dt2E = (t2E - ot2E) / flt_prms["dt"]

      # Update change
      tdt2E = wflowT * dt2E

      # Print result of measurement
      echo "WFLOW ", wflowT, " ", et, " ", es, " ", q, " ", t2E, " ", tdt2E

      # Breaking condition
      if (flt_prms["t_max"] > 0) and (wflowT > flt_prms["t_max"]):
         # Exit flow loop
         break

   # End timer
   toc("Wflow")

   # Print timing to output file
   tocc("Wilson flow:", t0)

#[ ~~~~ Define functions setting boundary conditions ~~~~ ]#

#[ Set condition in specific bounary ]#
proc BC(g: openArray[Field], mu: int) =
   # Get boundary condition
   let bc = str_prms["bc"][mu]

   # Get g along direction
   let gt = g[mu]

   # Check if boundary is anti-periodic
   if ($bc == $"a"):
      # Cycle through coordinates
      tfor i, 0..<gt.l.nSites:
         # Check if boundary
         if gt.l.coords[mu][i] == gt.l.physGeom[mu] - 1:
            # Multiply by -1
            gt{i} *= -1.0

#[ Set boundary conditions ]#
proc setBC*(g: openArray[Field]) =
   # Cycle through BC
   for mu in 0..<g.len:
      # Set boundary condition
      BC(g, mu)

#[ ~~~~ Define function for pbp measurement ~~~~ ]#

#[ For measuring pbp ]#
proc pbp(s: auto) =
   # Require number of stoch. srcs. > 0
   if int_prms["num_pbp_stoch_srcs"] > 0:
      # Start timer
      tic()

      # Get initial time
      let t0 = ticc()

      # Create separator
      echo ""

      # Cycle through pbp masses
      for mass_ind in 0..<int_prms["Nf"]:
         # Tell user what you're doing
         echo "Working on fermion ", mass_ind

         # Cycle through sources
         for src in 0..<int_prms["num_pbp_stoch_srcs"]:
            # Start thread block
            threads:
               # Create U(1) color source
               phi[mass_ind].u1 r_pbp

            # Calculate inverse
            s.solve_fermion(psi[mass_ind], phi[mass_ind], masses[mass_ind], sppbp)

            # Check if massless
            if masses[mass_ind] == 0:
               # Apply D^{d} and fill odd entries appropriately
               s.apply_massless_Ddag(psi[mass_ind], psi[mass_ind], "action")

            # Create thread block
            threads:
               # Measure norm
               var pbp = psi[mass_ind].norm2()

               # Set task for master thread
               threadMaster:
                  echo "MEASpbp mass, pbp/mass: ", masses[mass_ind], ", ", pbp/vol.float
      
      # Print timing information
      tocc("Chiral condensate measurement:", t0)

      # Create separator
      echo ""

      # End timer
      toc("Measure pbp")
   else:
      # Tell user the pbp not being measured
      echo "pbp not measured. Set # srcs > 0"

#[ ~~~~ Define functions for staggered operations ~~~~ ]#

#[ Fermion solve ]#
proc solve_fermion(s: Staggered; x,b: Field; mass: float; sp0: var SolverParams) =
   #[ Start timers ]#

   # Start timer
   tic()

   # Get initial time
   let t0 = ticc()

   #[ Do solve ]#
   
   # Check if fermion is massless
   if mass != 0:
      #[ Do full solve ]#

      # Do solve
      s.solve(x, b, mass, sp0)

      # End timer
      toc("Massive force solve")
   else:
      #[ Initialize CG parameters ]#

      # Create temporary solver params
      var sp = sp0

      # Reset stats
      sp.resetStats()

      # Set verbosity
      dec sp.verbosity

      # Set previous solution variable
      sp.usePrevSoln = false

      #[ Do solve and reconstruct ]#

      # Do solve of just even part
      s.solveEE(x, b, 0, sp)

      # Start thread block
      threads:
         # Copy even 4X even sites
         x.even := 4*x

      #[ Take care of stats ]#

      # Set calls
      sp.calls = 1

      # Set elapsed time
      sp.seconds = cpuTime() - t0

      # Add more flops to count (not sure if correct)
      sp.flops += float(24 * x.l.nEven)

      # Check verbosity
      if sp0.verbosity > 0:
         # Print information
         echo "stagSolve: ", sp.getStats()

      # Add stats to sp0
      sp0.addStats(sp)

   # End timer
   toc("Solve")

#[ Apply massless D^{dagger} ]#
proc apply_massless_Ddag(s: Staggered; x,b: Field; option: string) =
   # Start timer
   tic()

   # Set initial time
   let t0 = ticc()

   # Start thread block
   threads:
      # Apply 2Doe to phi_even and update x
      stagD2(s.so, x, s.g, b, 0, 0)

      # Set thread barrier
      threadBarrier()

      # Check option
      if option == "action":
         # Offset factor of 2 in 2Doe (also needs negative sign)
         x.odd := -0.5*x

         # Set even sites to zero
         x.even := 0
      elif option == "pv_force":
         # Set even sites to phi
         x.even := b

   # Print time taken
   tocc("D app.:", t0)

   # End timer
   toc("D app.")

#[ ~~~~ Define functions for momentum and field generation ~~~~ ]#

#[ Generate momenta ]#
proc generate_momenta(p: any, g: any) =
   # Start timer
   tic()
 
   # Set initial time
   let t0 = ticc()

   # Start thread block
   threads:
      # Create random momentum
      p.randomTAH r

      # Cycle through all momenta
      for i in 0..<g.len:
         # Set initial gauge variable
         g0[i] := g[i]

   # Print time taken
   tocc("Generate momenta and save gauge field:", t0)

   # End timer
   toc("Heatbath")

#[ Generate pseudofermion fields ]#
proc generate_pseudoferms(s: Staggered) =
   # Start timer
   tic()

   # Set initial time
   let t0 = ticc()

   # Cycle through fermion fields
   for mass_ind in 0..<masses.len:
      # Start thread block
      threads:
         # Create random fermion field
         psi[mass_ind].gaussian r

      # Check if regular fermion field
      if mass_ind < int_prms["Nf"]:
         # Apply staggered Dirac operator (D^{d})
         s.D(phi[mass_ind], psi[mass_ind], -masses[mass_ind])
      else:
         # Apply inverted staggered Dirac operator
         s.solve(phi[mass_ind], psi[mass_ind], -masses[mass_ind], spa[mass_ind])

      # Create thread block
      threads:
         # Set odd sites to zero
         phi[mass_ind].odd := 0

   # Print timing information
   tocc("Generate fermion/boson fields:", t0)

   # End timer
   toc("Generate fermion/boson fields")

#[ ~~~~ Define functions for action calculation ~~~~ ]#

#[ For calculating gauge action ]#
proc gaction(g: any, f2: seq, p2: float): auto =
   # Set variables
   var
      # Set gauge action
      ga = gc.actionA g

      # Set fermion action
      fa = sum(f2)

      # Get momentum from gauge field
      t = 0.5*p2 - float(16*vol)

   # Set hamiltonian
   let h = ga + fa + t

   # Return variables
   result = (ga, fa, t, h)

#[ Calculate action ]#
proc calc_action(s: Staggered, g: any, p: any): auto =
   #[ Initial timing info ]#

   # Start timer
   tic()

   # Get initial time
   let t0 = ticc()

   #[ Calculate momenta ]#
   
   # Initialize p2 to be zero
   var p2 = 0.0

   # Create thread block
   threads:
      # Create temporary p2 variable
      var p2_c = 0.0

      # Cycle throug sites
      for i in 0..<p.len:
         # Increment p2
         p2_c += p[i].norm2

      # Set p2 to p2t in master thread
      threadMaster: p2 = p2_c

   #[ Take care of fermion part of the action ]#
   
   # Initialize fermion action
   var f2 = newSeq[float](masses.len)

   # Cycle through fermion fields
   for mass_ind in 0..<masses.len:
      # Check if regular fermion
      if mass_ind < int_prms["Nf"]:
         # Do fermion solve
         s.solve_fermion(psi[mass_ind], phi[mass_ind], masses[mass_ind], spa[mass_ind])

         # Check if massless
         if masses[mass_ind] == 0:
            # Apply D^{d} and fill odd entries appropriately
            s.apply_massless_Ddag(psi[mass_ind], psi[mass_ind], "action")
      else:
         # Create thread block
         threads:
            # Apply D to phi to get psi for PV fermion
            s.D(psi[mass_ind], phi[mass_ind], masses[mass_ind])

      # Create temporary psi^{d}*psi variable
      var psi2 = 0.0

      # Start thread block
      threads:
         # Increment f2 (reuses psi from HMC trajectory)
         psi2 = psi[mass_ind].norm2()

      # Increment f2
      f2[mass_ind] = 0.5 * psi2

   #[ Calculate gauge action and put everything together ]#

   # Calculate gauge action and return whole action
   let (ga, fa, t, h) = g.gaction(f2, p2)

   #[ End timing information ]#
   
   # Print timing info
   tocc("Calculate action:", t0)

   # End timer
   toc("Calculate action")

   # Return results for gauge action
   result = (ga, fa, f2, t, h)

#[ ~~~~ Define functions for force calculation ~~~~ ]#

#[ For smearing force ]#
proc smearRephase(gauge: any, smeared_gauge: any): auto  =
   # Start timer
   tic()

   # Get initial time
   let t0 = ticc()

   # Get smeared force
   let smearedForce = coef.smearGetForce(gauge, smeared_gauge, info)

   # End timer
   toc("smear")

   # Start thread block
   threads:
      # Set boundary conditions
      smeared_gauge.setBC

      # Create thread barrier 
      threadBarrier()

      # Set staggered phases
      smeared_gauge.stagPhase
   
   # Print timing information
   tocc("Smear gauge and prepare smeared force:", t0)

   # End timer
   toc("BC & Phase")

   # Return force proc
   result = smearedForce

#[ For smearing gauge field ]#
proc smearRephaseDiscardForce(gauge: auto, smeared_gauge: auto) =
   # Start timer
   tic()

   # Get initial time
   let t0 = ticc()

   # Smear gauge field
   coef.smear(gauge, smeared_gauge, info)

   # Update user
   qexGC "smear done"

   # End timer
   toc("smear w/o force")

   # Start thread block
   threads:
      # Set boundary conditions
      smeared_gauge.setBC

      # Set thread barrier
      threadBarrier()

      # Set staggered phases
      smeared_gauge.stagPhase

   # Get timing information
   tocc("Smear fields:", t0)

   # End timer
   toc("BC & Phase")

#[ Rescaling for different fields ]#
proc rescale(index: int, t: float): float =
   # Get value of s
   var s = -0.5 * t

   # Check type of fermion
   if (index < int_prms["Nf"]) and (masses[index] != 0):
      # Add factors appropriate to regular fermion
      s = s / masses[index]
   elif index < int_prms["Nf"]:
      # Add factors appropriate to massess fermion
      s = -0.5 * s
   else:
      # Multiply by appropriate factor
      s = 0.5 * s

   # Return rescale
   result = s

#[ Smeared one link force ]#
proc smeared_one_link_force(f: auto, smeared_force: proc, g: auto) =
   
   #[ Correcting phase ]#
   
   # Start timer
   tic("One link force")

   # Get initial time
   let t0 = ticc()

   # Start thread block
   threads:
      # Set boundary conditions
      f.setBC

      # Set thread barrier
      threadBarrier()

      # Staggered phases
      f.stagPhase()

      # Set thread barrier
      threadBarrier()

      # Cycle through directions
      for mu in 0..<f.len:
         # Cycle through odd lattice sites
         for i in f[mu].odd:
            # Resign odd sites
            f[mu][i] *= -1

   # Print timing information
   tocc("Rephase force/set BC's:", t0)

   # End timer
   toc("Phase")

   #[ Smear force ]#

   # Smear
   f.smeared_force f

   # Print timing information
   tocc("Smear fermion/boson force:", t0)

   # End timer
   toc("Smear")

   #[ Multiply by link and project to traceless/anti-Hermitian component ]#

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

   # Print timing information
   tocc("Multiply force by gauge link and project TA:", t0)

#[ For initializing Dslash ]#
proc init_Dslash(f: auto, p: auto, t: auto; n: int, scale: float) =
   # Get initial time
   let t0 = ticc()

   # Start thread block
   threads:
      # Cycle through directions
      for mu in 0..<f.len:
         # Cycle through sites
         for i in f[mu]:
            # Cycle through n
            forO a, 0, n-1:
               # Cycle through n
               forO b, 0, n-1:
                  # Outer product with D on the right
                  f[mu][i][a, b] := scale * p[i][a] * t[mu].field[i][b].adj

   # Print timing information
   tocc("First contribution to Dslash:", t0)

#[ For adding to existing calculation of Dslash ]#
proc append_Dslash(f: auto, p: auto, t: auto; n: int, scale: float) =
   # Get initial time
   let t0 = ticc()

   # Start thread block
   threads:
      # Cycle through directions
      for mu in 0..<f.len:
         # Cycle through sites
         for i in f[mu]:
            # Cycle through n
            forO a, 0, n-1:
               # Cycle through n
               forO b, 0, n-1:
                  # Outer product with D on the right
                  f[mu][i][a, b] += scale * p[i][a] * t[mu].field[i][b].adj

   # Print timing
   tocc("Appended onto Dslash", t0) 

#[ Calculate smeared force ]#
proc fforce(s: auto, f: auto, sf: proc, g: auto,
            ix: openarray[int], ts: openarray[float]) =

   #[ Do Dslash ]#

   # Start timer
   tic("Fermion force")

   # Get initial time
   let t0 = ticc()

   # Create shifter
   var t: array[4, Shifter[typeof(psi[0]), typeof(psi[0][0])]]

   # Cycle through directions
   for mu in 0..<f.len:
      # Set shifter
      t[mu] = newShifter(psi[0], mu, 1)

   # Create variable controlling f updating behavior
   var first = true

   # Cycle through indices
   for f_ind in ix:
      #[ Solves or application of massive D ]#

      # Start timer
      tic("Pseudofermion loop")

      # Temporary copy of f_ind for threads (fix issue with memory safety)
      var f_ind_copy = f_ind

      # Check if regular fermion (place where things can go wrong)
      if f_ind < int_prms["Nf"]:
         # Do solve
         s.solve_fermion(psi[f_ind], phi[f_ind], masses[f_ind], spf[f_ind])

         # End timer
         toc("Force solve")
      
         # Check if massless
         if masses[f_ind] == 0:
            # Apply D^{d} and fill odd entries appropriately
            s.apply_massless_Ddag(psi[f_ind], psi[f_ind], "ferm_force")
      else:
         # Apply D^{d} and fill odd entries appropriately
         s.apply_massless_Ddag(psi[f_ind], phi[f_ind], "pv_force")

      #[ Application of massless D and creation of outer product ]#

      # Get rescale
      let scale = rescale(f_ind, ts[f_ind])

      # Cycle through directions
      for mu in 0..<f.len:
         # Essentially apply staggered Dirac operator
         discard t[mu] ^* psi[f_ind]

      # Create variable for convenience
      const n = psi[0][0].len

      # Check if first computation
      if first: # Create force
         # Set first variable to false
         first = false

         # Get first contribution to Dslash
         init_Dslash(f, psi[f_ind], t, n, scale)
      else: # Update already created force
         # Calculate next contribution to Dslash
         append_Dslash(f, psi[f_ind], t, n, scale)

      # End timer
      toc("Outer")

   # Print timing information
   tocc("Dslash and outer product", t0)

   # End timer
   toc("Outer product")

   #[ Smear and rephase ]#

   # Calculate the smeared one link force
   f.smeared_one_link_force(sf, g)

   # Print timing information
   tocc("Full fermion/boson force calculation:", t0)

   # End timer
   toc("Smeared one link force")

#[ ~~~~ Define functions HMC integration ~~~~ ]#

#[ For saving backup gauge field ]#
proc fgsave =
   # Get initial time
   let t0 = ticc()

   # Start thread block
   threads:
      # Cycle through directions
      for mu in 0..<g.len:
         # Save gauge field
         gg[mu] := g[mu]

   # Print timing information
   tocc("Set backup gauge field", t0)

#[ For loading backup gauge field ]#
proc fgload = 
   # Get initial time
   let t0 = ticc()

   # Start thread block
   threads:
      # Cycle through directions
      for mu in 0..<g.len:
         # Load gauge field
         g[mu] := gg[mu]
 
   # Print timing information
   tocc("Loaded gauge field from backup:", t0)

#[ Gauge field update from momenta ]#
proc mdt(t: float) =
   # Start timer
   tic()

   # Get t0
   let t0 = ticc()

   # Start thread block
   threads:
      # Cycle through directions
      for mu in 0..<g.len:
         # Cycle through lattice sites
         for s in g[mu]:
            # Update gauge field
            g[mu][s] := exp(t*p[mu][s])*g[mu][s]
  
   # Print timing information
   tocc("Gauge field update:", t0)

   # End timer
   toc("mdt")

#[ Momentum update from gauge sector ]#
proc mdv(t: float) =
   # Start timer
   tic()

   # Get initial time
   let t0 = ticc()

   # Get gauge force
   gc.forceA(g, f)

   # Start thread block
   threads:
      # Cycle through directions     
      for mu in 0..<f.len:
         # Update momenta
         p[mu] -= t*f[mu]

   # Print timing information
   tocc("Gauge field momentum update from gauge sector:", t0)

   # End timer
   toc("mdv")

#[ Momentum update w/ fermions/bosons ]#
proc mdvf(ix: openarray[int], sf: proc, ts: openarray[float]) =
   # Start timer
   tic()

   # Get initial time
   let t0 = ticc()

   # Calculate force
   stag.fforce(f, sf, g, ix, ts)

   # Update user
   qexGC "mdvf fforce"

   # Start thread block
   threads:
      # Cycle through directions
      for mu in 0..<f.len:
         # Update momenta
         p[mu] -= f[mu]

   # Print timing information
   tocc("Gauge field momentum update from fermion/boson sector:", t0)

   # End timer
   toc("mdvf")

#[ Force gradient update from gauge sector ]#
proc fgv(t: float) = 
   # Start timer
   tic()

   # Get initial time
   let t0 = ticc()

   # Calculate gauge force
   gc.forceA(gg, f)

   # Update user
   qexGC "fgv forceA"

   # Start thread block
   threads:
      # Cycle through directions
      for mu in 0..<g.len:
         # Cycle through lattice sites
         for s in g[mu]:
            # Update gauge field from backup
            g[mu][s] := exp((-t)*f[mu][s])*g[mu][s]

   # Print timing information
   tocc("Gauge field force-gradient update from gauge sector:", t0)

   # End timer
   toc("fgv")

#[ Force gradient update from fermion sector ]#
proc fgvf(ix: openarray[int], sf: proc, ts: openarray[float]) =
   # Start timer
   tic()

   # Get initial time
   let t0 = ticc()

   # Calculate force
   stagg.fforce(f, sf, gg, ix, ts)

   # Update user
   qexGC "fgvf fforce"

   # Start thread block
   threads:
      # Cycle through directions
      for mu in 0..<f.len:
         # Cycle through lattice sites
         for s in g[mu]:
            # Update gauge field (sign may introduce error spot)
            g[mu][s] := exp((-1)*f[mu][s])*g[mu][s]

   # Print timing information
   tocc("Gauge field force-gradient update from fermion/boson sector:", t0)

   # End timer
   toc("fgvf")

#[ Integration with shared update ]#
proc mdvAllfga(ts, gs: openarray[float]) =
   #[ Initial setup ]#
   
   # Start timer
   tic("Combined update")

   # Get initial time
   let t0 = ticc()

   # First, define variables
   var
      # Shared force proc
      sforceShared: typeof(g.smearRephase sg) = nil

      # Option for updating gauge
      updateG = false

      # Option for updating backup gauge
      updateGG = false

      # Option for updating fermions
      updateF = newseq[int](0)

      # Option for updating fermions w/ backup gauge
      updateFG = newseq[int](0)

   # Define type for steps
   type
      # Backup gauge step type
      GGstep = tuple[t, g: float]

      # Fermion w/ backup gauge type
      FGstep = tuple[t, g: seq[float]]

   # Define arrays for approximation coefficients
   var
      # For backup gauge
      ggs: array[2, GGstep]

      # For fermion w/ backup gauge
      fgs: array[2, FGstep]

   # Set up approximation coefficient sequences
   for order in 0..1:
      # Fill fgs time
      fgs[order].t = newseq[float](gs.len - 1)

      # Fill fgs g time
      fgs[order].g = newseq[float](gs.len - 1)

   #[ Determine whether or not to update gauge ]#
   
   # Check if to update backup gauge
   if gs[0] != 0:
      # Set update for backup gauge to True
      updateGG = true
      
      # Check if non-backup gauge not updated
      if ts[0] == 0:
         # If so, throw error
         qexError: "Force gradient without the force update."

      # Check approximation order
      if int_prms["approx_order"] == 2:
         # Get coefficients
         let (tf, tg) = approximateFGcoeff2(ts[0], gs[0])

         # Cycle through orders
         for order in 0..1:
            # Fill backup gauge coefficients
            ggs[order] = (t: tf[order], g: tg[order])
      else:
         # Otherwise, get first-order coefficients only
         let (tf, tg) = approximateFGcoeff(ts[0], gs[0])

         # Fill backup gauge coefficients
         ggs[0] = (t: tf, g: tg)
   elif ts[0] != 0:
      # Otherwise, just update gauge
      updateG = true

   #[ Determine whether or not to update fermion ]#
   
   # Cycle through pseudofermion/boson fields
   for k in 0..<gs.len-1:   
      # Update for correct indexing
      let i = k + 1

      # Check if fermion w/ backup to be updated
      if gs[i] != 0:
         # Set update for backup gauge to True
         updateFG.add k

         # Check if non-backup gauge not updated
         if ts[i] == 0:
            # If so, throw error
            qexError: "Force gradient without the force update."

         # Check approximation order
         if int_prms["approx_order"] == 2:
            # Get coefficients
            let (tf, tg) = approximateFGcoeff2(ts[i], gs[i])

            # Cycle through orders
            for order in 0..1:
               # Fill fermion w/ backup gauge coefficients
               fgs[order].t[k] = tf[order]

               # Fill fermion w/ backup gauge coefficients
               fgs[order].g[k] = tg[order]
         else:
            # Otherwise, get first-order coefficients only
            let (tf, tg) = approximateFGcoeff(ts[i], gs[i])

            # Fill fermion w/ backup gauge coefficients
            fgs[0].t[k] = tf

            # Fill fermion w/ backup gauge coefficients
            fgs[0].g[k] = tg
      elif ts[i] != 0:
         # Otherwise, just update fermions
         updateF.add k

   #[ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
     Note from the QEX developers:
     
     ~ Note that the smeared force proc retains a reference to the input
     gauge field.  In order to reuse the force proc, we need to use the
     correct gauge field that would remain the same.
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ]#

   #[ Prepare for update ]#
   
   # Check if backup field to be used
   if updateGG or updateFG.len > 0:
      # Start timer
      tic()

      # Save backup gauge field
      fgsave()

      # End timer
      toc("FG save")

      # Check if fermions to be updated with backup gauge
      if updateFG.len > 0:
         # Create shared force
         sforceShared = gg.smearRephase sgg
         
         # End timer
         toc("FG copy smeared")

         # Check if fermions with regular gauge to be updated
         if updateF.len > 0:
            # Start thread block
            threads:
               # Cycle through directions
               for mu in 0..<sg.len:
                  # Fill gauge field with backup
                  sg[mu] := sgg[mu]
            
            # End timer
            toc("Done")
         # End timer
         toc("Done")       
      # End timer
      toc("Done")

   #[ Do molecular dynamics ]#

   # Check if gauge field to be updated
   if updateG:
      # Start timer
      tic("MD mdv")

      # Update gauge sector
      mdv ts[0]

      # End timer
      toc("Done")

   # Check if fermions to updated
   if updateF.len > 0:
      # Start timer
      tic("MD mdvf")
 
      # Check if to be updated with backup
      if updateFG.len == 0:
         # Create smeared sforce proc
         var sforcef = g.smearRephase sg

         # Update fermion sector
         mdvf(updateF, sforcef, ts[1..^1])

         # Set sforce proc to nil
         sforcef = nil

         # Update user
         qexGC "MD mdvf w/ smear done"
      else:
         # Do molecular dynamics w/ shared force proc
         mdvf(updateF, sforceShared, ts[1..^1])

         # Update user
         qexGC "MD mdvf done"

      # End timer
      toc("Done")

   #[ Do force gradient step ]#

   # Check if backup gauge to be updated
   if updateGG or updateFG.len > 0:
      # Start timer
      tic("updateFG")

      # Cycle through approximation orders
      for order in 0..<int_prms["approx_order"]:
         # Start timer
         tic()

         # Check if backup gauge to be updated
         if updateGG:
            # Start timer
            tic("fgv")

            # Do force gradient
            fgv ggs[order].g

            # End timer
            toc("Done")

         # Check if fermion w/ backup to be updated
         if updateFG.len > 0:
            # Start timer
            tic("fgvf")

            # Do force gradient
            fgvf(updateFG, sforceShared, fgs[order].g)

            # Check order
            if order + 1 == int_prms["approx_order"]: sforceShared = nil

            # Update user
            qexGC "fgvf done"

            # End timer
            toc("FG")
       
         # Check if backup gauge to be updated
         if updateGG:
            # Start timer
            tic("FG mdv")

            # Do molecular dynamics
            mdv ggs[order].t

            # End timer
            toc("Done")

         # Check again if fermion w/ backup gauge to be updated
         if updateFG.len > 0:
            # Start timer
            tic("FG mdvf")
     
            # Create force proc
            var sforceg = g.smearRephase sg

            # End timer
            toc("FG smear rephase temp")

            # Do molecular dynamics
            mdvf(updateFG, sforceg, fgs[order].t)

            # Set sforce proc to nil
            sforceg = nil

            # Update user
            qexGC "FG mdvf done"

            # End timer
            toc("FG MD")

         # Load backup gauge
         fgload()

         # End timer
         toc("Load")
      # End timer
      toc("Done")      
   # End timer
   toc("Done")

   #[ Print timing information ]#

   # Update update number
   update_num = update_num + 1

   # Print timing information
   tocc("Integrator update " & intToStr(update_num) & ": ", t0)

#[ Set up integrator ]#
let
   # Define gauge integration algorithm
   gauge_int_alg: IntProc = str_prms["gauge_int_alg"]

   # Define fermion integration algorithm
   ferm_int_alg: IntProc = str_prms["ferm_int_alg"]

   # Define Pauli-Villars integration algorithm
   pv_int_alg: IntProc = str_prms["pv_int_alg"]

   # Create integration pair
   (V, T) = newIntegratorPair(mdvAllfga, mdt)

   # Set integrator for gauge
   H = newParallelEvolution gauge_int_alg(T = T, V = V[0], steps = int_prms["g_steps"])
block:
   # Cycle through masses
   for mass_ind in 0..<masses.len:
      # Check if dealing with regular fermions
      if mass_ind < int_prms["Nf"]:
         # Set integrator for regular fermion
         H.add ferm_int_alg(T = T, V = V[mass_ind + 1], steps = int_prms["f_steps"])
      else:
         # Add integrator for pv fermion
         H.add pv_int_alg(T = T, V = V[mass_ind + 1], steps = int_prms["pv_steps"])

#[ ~~~~ Define functions for checks of HMC ~~~~ ]#

#[ For checking individual solvers ]#
proc checkStats(label: string, sp: var SolverParams) =
   # Print stats
   echo label, sp.getAveStats

   # Check r2 max
   if sp.r2.max > sp.r2req:
      # Throw qex error
      qexError "Max r2 larger than requested"

   # Reset stats
   sp.resetStats

#[ For checking solvers ]#
proc checkSolvers(traj: int) =
   # Check if pbp to be measured
   if (int_prms["pbp_freq"] > 0) and ((traj + 1) mod int_prms["pbp_freq"] == 0):
      # Print information about pbp solver
      checkStats("Solver [pbp]: ", sppbp)

   # Cycle through fermions
   for mass_ind in 0..<masses.len:
      # Print stats for action solver
      checkStats("Solver " & intToStr(mass_ind) & " [action]: ", spa[mass_ind])

      # Check mass index
      if mass_ind < int_prms["Nf"]:
         # Print stats for force solver
         checkStats("Solver " & intToStr(mass_ind) & " [force]: ", spf[mass_ind])

#[ Reversibility check ]#
proc rev_check(evol: auto; h0, ga0, T0, fa0: float; 
               h1, ga1, t1, fa1: float; f20, f21: seq) =
   #[ Define/fill temporary variables ]#
 
   # Get initial time
   let t_0 = ticc()

   # Tell user what you're doing
   echo "\n~~~~ Reversibility check ~~~~\n"

   # Create temporary variables for fields
   var
      # New gauge field
      g1 = lo.newgauge

      # New momentum
      p1 = lo.newgauge

   # Create thread block
   threads:
      # Cycle through gauge field components
      for i in 0..<g1.len:
         # Fill temporary gauge field
         g1[i] := g[i]

         # Fill temporary momentum
         p1[i] := p[i]

         # Fill reversed momentum
         p[i] := -1*p[i]

   #[ Evolve gauge field ]#

   # Evolve gauge field
   evol.evolve flt_prms["tau"]

   # Finish evolution
   evol.finish

   #[ Smear and calculate action ]#

   # Now resmear gauge field
   g.smearRephaseDiscardForce sg

   # Calculate action after evolution
   let (gar, far, f2r, tr, hr) = stag.calc_action(g, p)

   # Create new sequence for change in ferm. from rev
   var df2r = newseq[float](f2r.len)

   # Create new sequence for change in ferm. from init.
   var df20 = newseq[float](f2r.len)

   # Fill df2
   for ind in 0..<df2r.len:
      # Calculate difference from rev.
      df2r[ind] = f2r[ind] - f21[ind]

      # Calculate difference from init
      df20[ind] = f2r[ind] - f20[ind]

   #[ Show results ]#
   
   # Print information about reversed H
   echo "Reversed H: ", hr, " Sg: ", gar, ", Sf: ", far, ", far (indiv.): ", f2r, ", T: ", tr

   # Print change in Hamiltonian from before and after reversed traj.
   echo "dH: ",hr-h1, " dSg: ",gar-ga1, " dSf: ",far-fa1, ", dSf (indiv.) ",df2r, " dT ",tr-t1

   # Print changes from initial Hamiltonian
   echo "dH0: ",hr-h0," dSg0: ",gar-ga0," dSf0: ",far-fa0,", dSf0 (indiv.) ",df20," dT0 ",tr-T0

   #[ Restore state and print timing ]#

   # Start thread block
   threads:
      # Cycle through gauge field components
      for i in 0..<g1.len:
         # Fill temporary gauge field
         g[i] := g1[i]

         # Fill temporary momentum
         p[i] := p1[i]

   # Print timing information
   tocc("Reversibility check:", t_0)

   # Create space
   echo ""

#[ ~~~~ Do HMC ~~~~ ]#

# Create separator
echo "\n~~~~ Integrator information ~~~~\n"

# Print information about gauge field evolution
echo H

#[ Cycle through configurations ]#
for config in start_config..<end_config:
   
   # Get initial time for config
   let config_init_time = ticc()

   #[ Cycle through trajectories ]#
   for traj in 0..<config_space:
      
      #[ Start HMC ]#

      # Start timer
      tic()

      # Get initial time for trajectory
      let traj_init_time = ticc()

      # Define ratio string
      let rat_str = intToStr(traj + 1) & "/" & intToStr(config_space)

      # Tell user what trajectory that you're on
      echo "\n~~~~ Trajectory #: " & intToStr(config) & " " & rat_str & " ~~~~\n"

      # Initialize update number
      update_num = 0

      #[ Generate field momenta and save ]#

      # Generate random momenta and save initial gauge field
      p.generate_momenta(g)

      # Finalize timer for setting up initial momenta
      toc("HB momenta for gauge created")

      #[ Smear fields and save ]#

      # Now smear and rephase
      g.smearRephaseDiscardForce sg

      # Tell user that you smeared gauge field
      toc("Smeared gauge field")

      #[ Generate pseudofermion/boson fields ]#

      # Generate pseudofermion/boson fields
      stag.generate_pseudoferms()

      #[ Calculate initial action ]#

      # Calculate action
      let (ga0, fa0, f20, t0, h0) = stag.calc_action(g, p)

      # Print information about initial action out
      echo "Beginning H: ", h0, " Sg: ", ga0, ", Sf: ", fa0, ", Sf (indiv.): ", f20, ", T: ", t0

      #[ Do trajectory ]#

      # Evolve gauge field
      H.evolve flt_prms["tau"]

      # Finish Evolution
      H.finish

      # Stop timer
      toc("evolve")    

      #[ Smear gauge field again ]#
   
      # Now resmear gauge field
      g.smearRephaseDiscardForce sg

      #[ Calculate final action ]#

      # Calculate action
      let (ga1, fa1, f21, t1, h1) = stag.calc_action(g, p)

      # Create array for getting individual changes in fermion/boson action
      var df2 = newseq[float](f21.len)

      # Cycle through individual fermions/bosons
      for ind in 0..<df2.len:
         # Calculate difference
         df2[ind] = f21[ind] - f20[ind]

      # Print information about final action out
      echo "Ending H: ", h1, " Sg: ", ga1, ", Sf: ", fa1, ", Sf (indiv.): ", f21, ", T: ", t1

      # Print information about change
      echo "dg, d(t + g), df2, df2 (indiv.): ",ga1-ga0, ", ",t1-t0+ga1-ga0, ", ",fa1-fa0, ", ",df2

      #[ Checks ]#

      # Check if reversibility to be checked
      if (int_prms["rev_check_freq"] > 0) and (traj mod int_prms["rev_check_freq"] == 0):
         # Do reversibility check
         H.rev_check(h0, ga0, t0, fa0, h1, ga1, t1, fa1, f20, f21)

      #[ Metropolis step ]#
   
      # Set dH
      var dH = h1 - h0

      # Set acceptance probability
      var acc = exp(-dH)

      # For random number
      var accr = 0.0

      # Check if NaN
      if dH.classify == fcNaN:
         # Create separator
         echo "\n\n"

         # Print warning message
         echo "WARNING! WARNING! WARNING! dH is Nan!"

         # Create separator
         echo "\n\n"

         # Set dH to infinity
         dH = Inf

         # Set acceptance to zero
         acc = 0

      # Check to make sure above user-set Metropolis threshold
      if config >= int_prms["no_metropolis_until"]:
         # Draw random number 
         accr = R.uniform

         # Do Metropolis test
         if accr <= acc:
            # Tell user that new configuration has been accepted
            echo "ACCEPT: dH: ", dH,"  exp(-dH): ", acc,"  r: ", accr
         else:
            # Otherwise, tell user that trajectory rejected
            echo "REJECT:  dH: ", dH, " exp(-dH): ", acc, "  r: ",accr

            # Create thread block
            threads:
               # Cycle through gauge field
               for i in 0..<g.len:
                  # Set gauge field to original gauge field
                  g[i] := g0[i]
      else:
         # Tell user that new configuration has been accepted
         echo "ACCEPT (no metrop. test performed): dH: ", dH,"  exp(-dH): ",acc,"  r: ",accr

      #[ Save information ]#

      # Check if gauge config to be saved
      if (save_freq > 0) and ((config + 1) mod save_freq == 0) and (traj == config_space - 1):
         # Start timer
         tic("save")

         # Create filename
         let fn = io_path & def_fn & "_" & intToStr(config + 1)
   
         # Save global rng, rng field, and gauge configuration
         save_gaugefield_and_rng(fn, g, R, r) 

         # End timer
         toc("Finished saving configuration")

      #[ Do measurements ]#

      # Check if trajectory was accepted again
      if accr <= acc:
         # Start timer
         tic("Reunit and smear")

         # Reunitarize gauge field
         g.reunit

         # Smear gauge field
         g.smearRephaseDiscardForce sg   

         # End timer
         toc("Reunit and smear")

      # Check if plaquette to be measured
      if (int_prms["plaq_freq"] > 0) and ((traj + 1) mod int_prms["plaq_freq"] == 0):
         # Measure plaquette
         g.mplaq

      # Check if Polyakov loop to be measured
      if (int_prms["ploop_freq"] > 0) and ((traj + 1) mod int_prms["ploop_freq"] == 0):
         # Measure polyakov loop
         g.ploop

      # Check if pbp to be measured
      if (int_prms["pbp_freq"] > 0) and ((traj + 1) mod int_prms["pbp_freq"] == 0):
         # Measure pbp on g (potential error spot)
         stag.pbp

      # Check if Wilson flow measurements to be performed
      if (int_prms["wflow_freq"] > 0) and ((traj + 1) mod int_prms["wflow_freq"] == 0):
         # Do Wilson flow and Wilson flow measurements
         g.wflow()

      #[ Show final information about trajectory ]#

      # Check solvers
      if (int_prms["check_solvers"] > 0) and ((traj + 1) mod int_prms["check_solvers"] == 0):
         # Check solvers
         checkSolvers(traj)

      # Give user timing information for trajectory
      tocc("\nTotal time for config. #" & intToStr(config) & " " & rat_str & ":", traj_init_time)

   # Give user timing information for configuration
   tocc("\nTotal time for config. #" & intToStr(config) & ":", config_init_time)

#[ ~~~~ Finalize QEX ~~~~ ]#

# Check if user wants timers shown
if options["show_timers"]:
   # Print information about timers
   echoTimers(flt_prms["timer_expand_ratio"], 
              options["timer_echo_dropped"])

# Print end date
echo "\nEnd: ", now().utc, "\n" 

# Finalize QEX
qexfinalize()