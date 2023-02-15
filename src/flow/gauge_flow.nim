#[ ~~~~ Import necessary modules ~~~~ ]#

import qex # QEX
import gauge # Gauge field
import times # Timing
import macros # Useful macros
import os # For operating system-specific tasks
import streams, parsexml, strutils # For parsing XML
import parseopt # For parsing command line arguments
import tables # For organizing data
import sequtils # For dealing with sequences
import math # Basic mathematical operations
import flow # For intermediate measurements

#[ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Authors: Curtis Taylor Peterson, James Osborn, Xiaoyong Jin

Description:

   Gauge flow measurements. Currently only uses Wilson flow. Other
   flows/measurements to be added in the future.

Credits:

   Thank you to James Osborn Xiaoyong Jin for developing QEX and helping
   develop this program.

   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ]#

#[ ~~~~ Initialize QEX ~~~~ ]#

# Initialize QEX
qexinit()

# Print start date
echo "\nStart: ", now().utc

#[ ~~~~ Initialize info for whole file ~~~~ ]#

# Initialize attribute name
var attrName = ""

# Initialize attribute value
var attrVal = ""

# Define maximum flow times
var max_flts = newseq[float](0)

# Define dt's
var dts = newseq[float](0)

# Define integer parameters
var int_prms = {"Ns" : 0, "Nt" : 0, "f_munu_loop" : 0}.toTable

# Define float parameters
var flt_prms = {"t_max" : 0.0, "beta_w" : 1.0,
                "beta_r" : 1.0, "c1" : 1.0, "beta_a" : 1.0,
                "adj_plaq" : -0.25}.toTable

# Define string parameters
var str_prms = {"flow_act" : "Wilson"}.toTable

# Initialize starting trajectory
var start_config = 0

# Initialize ending trajectory
var end_config = 0

# Initialize filename
var fn = "checkpoint"

# Initialize xml file
var xml_file = ""

# Initialize rank geometry
var rank_geom = @[1, 1, 1, 1]

#[ ~~~~ Parallel information ~~~~ ]#

# Print number of ranks
echo "# ranks: ", nRanks

# Print number of threads
threads: echo "# threads: ", numThreads

#[ ~~~~ For timing ~~~~ ]#

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

   # Start thread block and reunitarize
   threads:
      let d = g.checkSU
      threadBarrier()
      echo "unitary deviation avg: ",d.avg," max: ",d.max
      g.projectSU
      threadBarrier()
      let dd = g.checkSU
      echo "new unitary deviation avg: ",dd.avg," max: ",dd.max

   # End timer
   toc("reunit")

#[ For measuring plaquette ]#
proc meas_plaq(g: auto): auto =
   # Start timer
   tic()

   # Calculate plaquette
   let
      pl = g.plaq
      nl = pl.len div 2
      ps = pl[0..<nl].sum * 2.0
      pt = pl[nl..^1].sum * 2.0

   # End timer
   toc("plaq")

   # Return result
   result = (ps, pt)

#[ For measuring polyakov loop ]#
proc meas_ploop(g: auto): auto =
   # Start timer
   tic()

   # Calculate Polyakov loop
   let pg = g[0].l.physGeom
   var pl = newseq[typeof(g.wline @[1])](pg.len)
   for i in 0..<pg.len:
      pl[i] = g.wline repeat(i+1, pg[i])
   let
      pls = pl[0..^2].sum / float(pl.len-1)
      plt = pl[^1]

   # End timer
   toc("ploop")

   # Return Polyakov loop
   result = (pls, plt)

#[ ~~~~ Functions for initialization ~~~~ ]#

#[ For reading flow parameters ]#
proc initialize_gauge_field_and_params(): auto = 
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
            
            # Check if starting config.
            if cm_opts.key == "start_config":
               # Set ending config.
               start_config = parseInt(cm_opts.val)

               # Print ending config.
               echo "start config: " & cm_opts.val

            # Check if ending config.
            if cm_opts.key == "end_config":
               # Set ending trajectory
               end_config = parseInt(cm_opts.val)

               # Print ending config.
               echo "end config: " & cm_opts.val

            # Check if base file name
            if cm_opts.key == "filename":
               # Set filename
               fn = cm_opts.val

               # Print filename
               echo "flowing gauge file: " & cm_opts.val
           
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
               elif str_prms.hasKey(attrName):
                  # Save parameter
                  str_prms[attrName] = attrVal
               elif attrName[0..3] == "time":
                  # Add to list of times
                  max_flts.add(parseFloat(attrVal))
               elif attrName[0..1] == "dt":
                  # Add to list of increments
                  dts.add(parseFloat(attrVal))

               # Print variable
               echo attrName & ": " & attrVal

            of xmlEof: break # If end of file, exit
            else: discard # Otherwise, do nothing

      # Close
      x.close()

   #[ ~~~~ Initialize lattice and gauge field ~~~~ ]#
   
   # Define various variables for gauge action
   let
      # Define lattice
      lat = intSeqParam("lat", @[int_prms["Ns"], int_prms["Ns"],
                                 int_prms["Ns"], int_prms["Nt"]])

      # Define new lattice layout
      lo = lat.newLayout(rank_geom)

   # Create new gauge
   var g = lo.newgauge

   #[ Initialize gauge action coefficients ]#

   # Initialize gauge action coefficients
   var gc: GaugeActionCoeffs

   # Start case for setting coefficeints
   case str_prms["flow_act"]
      of "Wilson": gc = GaugeActionCoeffs(plaq: flt_prms["beta_w"])
      of "adj": gc = GaugeActionCoeffs(plaq: flt_prms["beta_a"], 
                                       adjplaq: flt_prms["beta_a"] * flt_prms["adj_plaq"])
      of "rect": gc = gaugeActRect(flt_prms["beta_r"], flt_prms["c1"])

   #[ ~~~~ Print summary of information ~~~~ ]#

   # Print header
   echo "\n ~~~~ SUMMARY OF FLOW INFORMATION  ~~~~\n"

   # Start case
   case str_prms["flow_act"]:
      of "Wilson":
         # Print output
         echo "Wilson flow w/ beta_w = ", flt_prms["beta_w"]
      of "adj":
         # Temporary definitions
         let beta = flt_prms["beta_a"]
         let adjplaq = flt_prms["adj_plaq"]

         # Print output
         echo "Adj. plaq. flow w/ beta_f = ", beta, " & beta_f/beta_a = ", adj_plaq
      of "rect":
         # Temporary definitions
         let beta = flt_prms["beta_r"]
         let c1 = flt_prms["c1"]

         # Print output
         echo "Rect. flow w/ beta_r = ", beta, " & c1 = ", c1

   # Cycle through dt's
   for dt_ind in 0..<max_flts.len:
      # Print information
      echo "dt = ", dts[dt_ind], " up to t/a^2 = ", max_flts[dt_ind]

   #[ ~~~~ Return gauge field ~~~~ ]#
   # Return result
   result = (g, gc)

#[ ~~~~ Functions for gauge field IO ~~~~ ]#

#[ For reading in gauge flows ]#
proc read_gauge_file(gaugefile: string; g: auto) =
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

      # Reunitarize
      g.reunit

#[ ~~~~ Functions for gauge flow ~~~~ ]#

#[ For calculating flowed measurements ]#
proc EQ(gauge: auto, loop: int): auto =
   # Calculate measurements
   let
      # Get loop
      f = gauge.fmunu loop

      # Calculate Yang-Mills density
      (es, et) = f.densityE

      # Calculate topological charge
      q = f.topoQ

      # Measure Plaquette
      (ss, st) = gauge.meas_plaq()

      # Measure Polyakov loop
      (pls, plt) = gauge.meas_ploop()

   # Return Yang-Mills density and topology
   return (es, et, ss, st, q, pls, plt)

#[ Helper proc. for formatting floats ]#
proc format_float(num: float, precis: int): string =
   # Return formatted string
   result = num.formatFloat(ffDecimal, precis) & " "

#[ Helper proc. for printout ]#
proc print_info(dtau, tau, es, et, ss, st, q: float; pls, plt: auto; old_t2E: float): float =
   #[ Take care of some preliminary computations ]#

   let
      #[ Plaquettes ]#

      # Multiply spatial plaquette by 3
      plaq_ss = 3.0 * ss

      # Multiply tempoary plaquette by 3
      plaq_st = 3.0 * st

      # Average spatial and temporal plaquettes
      plaq = (plaq_ss + plaq_st) / 2.0

      # Calculate "check"
      check = 12.0 * tau * tau * (3.0 - plaq)

      #[ Clover operator ]#

      # Calculate clover operator averaged over space/time directions
      clov = es + et

      # Calculate t^2 E(t)
      t2E = tau * tau * clov

      # Calculate t^2 E(t) for ss
      t2E_ss = tau * tau * es

      # Calculate t^2 E(t) for st
      t2E_st = tau * tau * et

      # Calculate derivative
      der_t2E = (t2E - old_t2E) / dtau

      #[ Polyakov loops ]#

      # Get temporal Polyakov loop normalized by 3
      poly_t = plt * 3.0 

      # Get spatial Polyakov loop normalized by 3
      poly_s = pls * 3.0

   #[ Take care of extra details ]#

   # Define default pricision
   let def_pres = 13

   #[ Print out appropriate information ]#

   # Define string
   var printout = "FLOW "

   # Add flow time
   printout = printout & format_float(tau, 2)

   # Add plaquette (normalized to 3)
   printout = printout & format_float(plaq, def_pres)

   # Add clover [E(t)]
   printout = printout & format_float(clov, def_pres)

   # Add t^2 E(t)
   printout = printout & format_float(t2E, def_pres)

   # Add d t^2 E(t) / dt
   printout = printout & format_float(der_t2E, def_pres)

   # Add "check"
   printout = printout & format_float(check, def_pres)

   # Add topology
   printout = printout & format_float(q, def_pres)

   # Add t^2 E(t) for ss
   printout = printout & format_float(t2E_ss, def_pres)

   # Add t^2 E(t) for st
   printout = printout & format_float(t2E_st, def_pres)

   # Add real part of Polyakov loop
   printout = printout & format_float(poly_t.re, def_pres)

   # Add imaginary part of Polyakov loop
   printout = printout & format_float(poly_t.im, def_pres)

   # Add imaginary part of Polyakov loop
   printout = printout & format_float(poly_s.re, def_pres)

   # Add imaginary part of Polyakov loop
   printout = printout & format_float(poly_s.im, def_pres)

   #[ Print information out and return t^2 E(t) ]#

   # Print information out
   echo printout

   # Return t^2 E(t)
   result = t2E

#[ For flow and flow measurements ]#
proc flow(gc: GaugeActionCoeffs; flowed_gauge: auto; dt_ind: int; t2E: float): float =
   #[ Set things up ]#

   # Get dt
   let dt = dts[dt_ind]

   # Get maximum flow time
   let max_flt = max_flts[dt_ind]

   # Set last maximum flow time
   var last_max_flt = 0.0

   # Check flow time 
   if dt_ind != 0:
      # Reset last maximum flow time
      last_max_flt = max_flts[dt_ind - 1]

   # Set temprary t2E
   var t2E_temp = t2E

   #[ Do gauge flow ]#
   
   # Start flow loop
   gc.gaugeFlow(str_prms["flow_act"], flowed_gauge, last_max_flt, dt):
      # Breaking condition
      if (wflowT > max_flt) or ((flt_prms["t_max"] > 0) and (wflowT > flt_prms["t_max"])):
         # Exit flow loop
         break

      # Calculate measurement
      let (es, et, ss, st, q, pls, plt) = flowed_gauge.EQ int_prms["f_munu_loop"]

      # Print result of measurement
      t2E_temp = print_info(dt, wflowT, es, et, ss, st, q, pls, plt, t2E_temp)

   # Return t2E
   result = t2E_temp

#[ ~~~~ Perform gauge flow ~~~~ ]#

# Read in appropriate information and initialize gauge field
var (g, gc) = initialize_gauge_field_and_params()

#[ Cycle through configurations ]#
for config in start_config..<end_config + 1:

   #[ Take care of initial information ]#

   # Get initial time for config
   let config_time = ticc()
   
   # Tell user what config. that you're on
   echo "\n~~~~ Flowing config. # " & intToStr(config) & " ~~~~\n"

   #[ Read in appropriate gauge field ]#

   # Name gauge file
   let filename = fn & "_" & intToStr(config) & ".lat"

   # Load gauge field
   filename.read_gauge_file(g)

   #[ Flow gauge field ]#

   # Initialize measurements
   let (es, et, ss, st, q, pls, plt) = g.EQ int_prms["f_munu_loop"]

   # Define default string
   var t2E = print_info(dts[0], 0.0, es, et, ss, st, q, pls, plt, 0.0)

   # Cycle through different incremenets
   for dt_ind in 0..<max_flts.len:
      # Do gauge flow
      t2E = gc.flow(g, dt_ind, t2E)

   #[ Finalize this gauge flow ]#

   # Print out timing for this flow
   tocc("Flow for config. # " & intToStr(config) & ":", config_time)

# Print end date
echo "\nEnd: ", now().utc, "\n" 

# Finalize QEX
qexfinalize()