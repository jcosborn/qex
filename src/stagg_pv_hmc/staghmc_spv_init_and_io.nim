#[ ~~~~ Imports ~~~~ ]#

import qex # Import qex
import gauge/gaugeUtils # For handling gauge fields
import layout # For initializing lattice layout
import rng # For handing RNG
import parseopt # For parsing command line arguments
import tables # For organizing data
import streams, parsexml, strutils # For parsing XML
import system # For system-specific operations
import times # Timing for generating seeds
import os # For checking if files exist
import staghmc_spv_rng # For certain RNG operations
import options # For optional IO behavior with fields

#[ ~~~~ Command line and XML inputs ~~~~ ]#

proc print_banner() =
   # Create banner
   let banner = """
                                .oooooo.      oooooooooooo ooooooo  ooooo                                
                               d8P'  `Y8b     `888'     `8  `8888    d8'                                 
                              888      888     888            Y888..8P                                   
                              888      888     888oooo8        `8888'                                    
                              888      888     888    "       .8PY888.                                   
                              `88b    d88b     888       o   d8'  `888b                                  
                               `Y8bood8P'Ybd' o888ooooood8 o888o  o88888o                                
                                                                                                         
                                                                                                         
                                                                                                         
 .oooooo..o     .                                            ooooo   ooooo ooo        ooooo   .oooooo.   
d8P'    `Y8   .o8                                            `888'   `888' `88.       .888'  d8P'  `Y8b  
Y88bo.      .o888oo  .oooo.    .oooooooo  .oooooooo           888     888   888b     d'888  888          
 `"Y8888o.    888   `P  )88b  888' `88b  888' `88b            888ooooo888   8 Y88. .P  888  888          
     `"Y88b   888    .oP"888  888   888  888   888            888     888   8  `888'   888  888          
oo     .d8P   888 . d8(  888  `88bod8P'  `88bod8P'  .o.       888     888   8    Y     888  `88b    ooo  
8""88888P'    "888" `Y888""8o `8oooooo.  `8oooooo.  Y8P      o888o   o888o o8o        o888o  `Y8bood8P'  
                              d"     YD  d"     YD                                                       
                              "Y88888P'  "Y88888P'                                                       
                
Quantum EXpressions (QEX) authors: James Osborn & Xiao-Yong Jin
QEX Staggered HMC authors: James Osborn, Xiao-Yong Jin, Curtis Taylor Peterson, Anna Hasenfratz
QEX GitHub: https://github.com/jcosborn/qex
QEX Staggered HMC GitHub: https://github.com/ctpeterson/qex_staghmc 
   """

   # Print banner out
   echo banner

#[ Read command line inputs ]#
proc read_cmd*(): auto =
   #[ Print banner ]#

   # Print banner at top of output
   print_banner()

   #[ Initialize variables ]#

   # Set variables
   var
      # Starting configuration
      start_config = 0

      # Ending configuration
      end_config = 0

      # Frequency to save configurations
      save_freq = 1

      # Spacing (in trajectories) between configurations
      config_space = 1
      
      # XML file name
      xml_file = "input.xml"

      # MPI rank geometry
      rank_geom = newseq[int](0)

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
                  rank_geom.add parseInt(gm_str_splt[ind])

            # Check if user wants a different filename
            if cm_opts.key == "filename":
               # Set filename
               def_fn = cm_opts.val

            # Check if user wants a different path
            if cm_opts.key == "path":
               # Set path
               io_path = cm_opts.val 

   # Return parsed options as tuple
   result = (start_config, end_config, config_space, save_freq, 
             xml_file, rank_geom, def_fn, io_path)

#[ Read XML inputs ]#
proc read_xml*(xml_file: string): auto = 
   #[ Initialize variables ]#

   # Set variables
   var
      # Integer parameters
      int_prms = {"Ns" : 0, "Nt" : 0, "num_pv": 0,
                  "Nf" : 0, "a_maxits" : 0, "f_maxits" : 0,
                  "g_steps" : 0, "sg_steps": 0, "f_steps" : 0, "pv_steps" : 0,
                  "no_metropolis_until" : 0, "start_config" : 0,
                  "plaq_freq" : 0, "ploop_freq" : 0,
                  "rev_check_freq" : 0, "check_solvers" : 0, "sg_opt": 0,
                  "num_Ns" : 3, "num_Nt" : 1}.toTable

      # Float parameters
      flt_prms = {"beta" : 1.0, "adj_fac" : -1.0/4.0, "mass" : 0.0, 
                  "c1" : -1.0/12.0, "tau" : 0.0, "sm_adj_fac" : -1.0/4.0,
                  "sm_c1": -1.0/12.0, "sm_beta": 1.0,
                  "alpha_1" : 0.0, "alpha_2" : 0.0, "alpha_3" : 0.0,
                  "mass_pv" : 0.0, "a_tol" : 0.0, "f_tol" : 0.0}.toTable

      # Seed parameters
      seed_prms = {"parallel_seed" : intParam("seed", int(1000 * epochTime())).uint64,
                   "serial_seed" : intParam("seed", int(1000 * epochTime())).uint64}.toTable

      # String parameters
      str_prms = {"bc" : "pppa", "start" : "unit",
                  "gauge_int_alg" : "2MN", "ferm_int_alg" : "2MN",
                  "pv_int_alg" : "2MN", "smeared_gauge_int_alg": "2MN",
                  "rng_type" : "RngMilc6",
                  "gauge_act": "adjoint",
                  "smeared_gauge_act": "Wilson",
                  "gauge_smearing": "nhyp", "matter_smearing": "nhyp",
                  "plaq_verbose": "true"}.toTable

      # Initialize XML attribute name
      attrName = ""

      # Initialize XML attribute value
      attrVal = ""

   #[ Read XML file and return ]#   

   # Print what you're doing
   echo "\n ~~~~ XML information ~~~~\n"

   # Create XML parser
   var x: XmlParser

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
               elif str_prms.hasKey(attrName):
                  # Save string information
                  str_prms[attrName] = attrVal

               # Print variable
               echo attrName & ": " & attrVal
            of xmlEof: break # If end of file, exit
            else: discard # Otherwise, do nothing

      # Close
      x.close()

   # Return result at tuple
   result = (int_prms, flt_prms, seed_prms, str_prms)

#[ ~~~~ Gauge/RNG field initialization and IO ~~~~ ]#

#[ For reunitarization ]#
proc reunit*(g: auto) =
   # Start thread block and reunitarize
   threads:
      let d = g.checkSU
      threadBarrier()
      echo "unitary deviation avg: ",d.avg," max: ",d.max
      g.projectSU
      threadBarrier()
      let dd = g.checkSU
      echo "new unitary deviation avg: ",dd.avg," max: ",dd.max

#[ Initialize lattice layout ]#
proc init_layout*(phys_geom, rank_geom: seq): auto =
   # Define lattice
   let lat = intSeqParam("lat", phys_geom)

   # Return layout
   result = lat.newLayout(rank_geom)

#[ For initializing fields ]#
proc init_fields*(lo: Layout; gauge_start, rng_type: string; seed: uint64): auto =
   #[ Initialize RNG field ]#

   # Say what you're doing
   echo "\n ~~~~ Initialize gauge/RNG fields ~~~~\n"

   # Initialize RNG field
   var rng_field: ParallelRNG

   # Initialize and seed RNG field
   rng_field.init_parallel_rng(lo, rng_type, seed)

   #[ Initialize gauge field ]#

   # Initialize gauge field
   var g = lo.newgauge

   # Give gauge field particular start
   case gauge_start:
      of "cold": g.unit() # "Unit" start
      of "hot":
         # Start thread block
         threads:
            # Fill gauge field with random numbers 
            g.random(rng_field) # "Random" start
      of "read": discard # Do nothing
      else: quit("Unrecognized option for starting config. Exiting.")

   #[ Return field as tuple ]#   

   # Return fields at tuple
   result = (rng_field, g)

#[ For reading fields ]#
proc read_fields*(lo: Layout; base_fn, rng_type: string; seed: uint64;
                  current_config, start_config = none(int)): auto =
   #[ Initiaize fields and variables ]#

   # Say what you're doing
   echo "\n ~~~~ gauge/RNG IO ~~~~\n"

   # Initialize fields
   var (rng_field, g) = lo.init_fields("read", rng_type, seed)

   # Define filenames
   let
      # Name for RNG field
      rng_fn = base_fn & ".rng"

      # Name for gauge field
      gauge_fn = base_fn & ".lat"

   #[ Read RNG field ]#

   # Check to make sure that file exists
   if fileExists(rng_fn):
      # Read file
      rng_field.read_rng(rng_fn)
   elif (current_config.isSome and start_config.isSome) and (current_config == start_config):
      # Tell user what you're doing
      echo rng_fn & " does not exist." & " Creating new RNG field for just this config."
   else:
      # Exit program
      quit(rng_fn & " does not exist. Exiting.")

   #[ Read gauge field ]#

   # Check to make sure that file exists
   if fileExists(gauge_fn):
      # Read file and throw error if failure
      if 0 != g.loadGauge(gauge_fn):
         # Exit program there is an issue with reading file
         quit("Error reading gauge file. Exiting.")
   else:
      # Exit program
      quit(gauge_fn & " does not exist. Exiting.")

   # Check if gauge configuration needs to be reunitarized
   if (current_config.isSome and start_config.isSome) and (current_config == start_config):
      # Tell user what you're doing
      echo "Reunitarizing just this gauge configuration."

      # Reunitarize gauge field
      g.reunit()

   #[ Return fields ]#

   # Return fields as tuple
   result = (rng_field, g)

#[ For writing fields ]#
proc write_fields*(base_fn: string; rng_field: ParallelRNG; gf: auto) =
   #[ Initialize variables ]#

   # Tell user what you're doing
   echo "\n ~~~~ gauge/RNG IO ~~~~\n"

   # Initialize variables
   let
      # Name for RNG field
      rng_fn = base_fn & ".rng"

      # Name for gauge field
      gauge_fn = base_fn & ".lat"

   #[ Write RNG field ]#
   
   # Write RNG field
   rng_field.write_rng(rng_fn)

   #[ Write gauge field ]#

   # Attempt to save gauge configuration
   if 0 != gf.saveGauge(gauge_fn):
      # If not able to write gauge field, exit
      quit("Unable to write gauge field. Exiting.")

#[~~~~ Initialize pseudofermion/boson fields ~~~~ ]#

#[ Initialize pseudofermion and boson fields ]#
proc init_matter_fields*(lo: Layout, int_prms: Table[string, int]): auto =
   #[ Initialize matter fields ]#

   # Initialize number of fields
   let n_fields = int_prms["Nf"] + int_prms["num_pv"]

   # Initialize fermion field for temporary computations
   var psi = lo.ColorVector()

   # Initialize fermion type
   type clr_vec_type = typeof(psi)

   # Define generic "phi" fields (held constant in mol. dynmcs.)
   var phi = newseq[clr_vec_type](n_fields)

   # Cycle through field types
   for fld_ind in 0..<n_fields:
      # Initialize color vector
      phi[fld_ind] = lo.ColorVector()

   #[ Return fields ]#

   # Return fields as tuple
   result = (psi, phi)

#[ ~~~~ Initialize CG ~~~~ ]#

#[ Initialize CG ]#
proc init_cg*(flt_prms: Table[string, float], int_prms: Table[string, int]): auto =
   #[ Set variables up ]#
   
   # Set number of each field type
   let n_fields = int_prms["Nf"] + int_prms["num_pv"]

   # Create sequences of solver parameters
   var
      # Create dummy solverParams for type specification
      spd = initSolverParams()

      # Action solver parameters
      spa = newseq[typeof(spd)](n_fields)

      # Force solver parameters
      spf = newseq[typeof(spd)](int_prms["Nf"])

   #[ Initialize solvers ]#
   
   # Cycle through fields
   for fld_ind in 0..<n_fields:
      #[ Initialize action solver parameters ]#
         
      # Initialize solver parameters
      spa[fld_ind] = initSolverParams()

      # Set tolerance
      spa[fld_ind].r2req = flt_prms["a_tol"]

      # Set maximum number of CG iterations
      spa[fld_ind].maxits = int_prms["a_maxits"]

      # Set verbosity
      spa[fld_ind].verbosity = 1

      #[ Initialize force solver parameters ]#

      # Check is pseudofermion field
      if fld_ind < int_prms["Nf"]:
         # Initialize solver parameters
         spf[fld_ind] = initSolverParams()

         # Set tolerance
         spf[fld_ind].r2req = flt_prms["f_tol"]

         # Set maximum number of CG iterations
         spf[fld_ind].maxits = int_prms["f_maxits"]

         # Set verbosity
         spf[fld_ind].verbosity = 1

   #[ Return solver parameters ]#

   # Return solver parameters as tuple
   result = (spa, spf)