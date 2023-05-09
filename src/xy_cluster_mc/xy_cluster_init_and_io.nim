#[ ~~~~ Imports ~~~~ ]#
import qex # Import qex
import parseopt # For parsing command line arguments
import tables # For organizing data
import streams, parsexml, strutils # For parsing XML
import system # For system-specific operations
import times # Timing for generating seeds
import os # For checking if files exist

#[ ~~~~ Command line and XML inputs ~~~~ ]#

#[ Read command line inputs ]#
proc read_cmd*(): auto =
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
      rank_geom = @[1, 1]

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
                  rank_geom[ind] = parseInt(gm_str_splt[ind])

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
      int_prms = {"Nx" : 0, "Ny" : 0, "num_pv": 0,
                  "no_metropolis_until" : 0, "start_config" : 0,
                  "meas_freq" : 0, "start_config" : 0}.toTable

      # Float parameters
      flt_prms = {"J" : 0.0}.toTable

      # Seed parameters
      seed_prms = {"serial_seed" : intParam("seed", int(1000 * epochTime())).uint64}.toTable

      # String parameters
      str_prms = {"start" : "hot"}.toTable

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