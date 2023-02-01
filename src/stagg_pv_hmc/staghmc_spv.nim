#[ ~~~~ Import necessary modules ~~~~ ]#

import qex # QEX
import staghmc_spv_modules # For custom functions
import gauge # Gauge field
import gauge/hypsmear # Import nHYP smearing
import physics/[qcdTypes, stagSolve] # Staggered fermions
import mdevolve # MD evolution
import times # Timing
import macros # Useful macros
import tables # For organizing data
import sequtils # For dealing with sequences
import strutils # For manipulating strings
import math # Basic mathematical operations
import algorithms/integrator # For integrator options

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

   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ]#

#[ ~~~~ Initialize QEX ~~~~ ]#

# Initialize QEX
qexinit()

# Print start date
echo "\nStart: ", now().utc

echo "\n ~~~~ Parallel information  ~~~~\n"

# Print number of ranks
echo "# ranks: ", nRanks

# Print number of threads
threads: echo "# threads: ", numThreads

#[ ~~~~ Initial IO ~~~~ ]#

# Create command line inputs
let (start_config, end_config, config_space, save_freq,
     xml_file, rank_geom, def_fn, io_path) = read_cmd()

# Read XML file input
let (int_prms, flt_prms, seed_prms, str_prms) = read_xml(xml_file)

# Define base filename
var base_fn = io_path & def_fn & "_" & intToStr(start_config)

#[ ~~~~ Setup global RNG and CG ~~~~ ]#

# Create RNG object (custom, not native to QEX)
var R: GlobalRNG

# Set rng type
R = GlobalRNG(rng_type: str_prms["rng_type"])

# Check if starting new ensemble
if start_config == 0:
   # Seed RNG
   R.seed(seed_prms["serial_seed"])
else:
   # If not, read previous RNG file
   R.read_rng(base_fn)

# Initialize CG
var (spa, spf) = init_cg(flt_prms, int_prms)

#[ ~~~~ Set fields up ~~~~ ]#

# Initialize lattice layout
var lo = init_layout(@[int_prms["Ns"], int_prms["Ns"],
                       int_prms["Ns"], int_prms["Nt"]],
                     rank_geom)

# Initialize fields
var
   # RNG field
   r: ParallelRNG

   # Gauge field
   g: type(lo.newgauge)

   # Backup copy of gauge field
   g0 = lo.newgauge

   # Smeared gauge field
   sg = lo.newgauge

   # Gauge momentum
   p = lo.newgauge

   # Force f
   f = lo.newgauge

# Set up adjoint gauge action coefficients
let gc = GaugeActionCoeffs(plaq: flt_prms["beta"],
                           adjplaq: flt_prms["beta"] * flt_prms["adj_fac"])

# Check if starting w/ first configuration
if start_config == 0:
   # Set initialize gauge and RNG field
   (r, g) = lo.init_fields(str_prms["start"], str_prms["rng_type"],
                           seed_prms["parallel_seed"])
else:
   # Otherwise, pick up gauge configuration and RNG field fork disk
   (r, g) = lo.read_fields(base_fn, str_prms["rng_type"],
                           seed_prms["parallel_seed"])

# Initialize nHYP smearing parameters
var
   # Get info
   info: PerfInfo

   # Get nHYP coefficients
   coef = HypCoefs(alpha1: flt_prms["alpha_1"],
                   alpha2: flt_prms["alpha_2"],
                   alpha3: flt_prms["alpha_3"])

# Initialize fermion fields
var (psi, phi) = lo.init_matter_fields(int_prms)

# Set masses
let masses = @[flt_prms["mass"], flt_prms["mass_pv"]]

# Initialize stag for operations involving smeared stag op
let stag = newStag(sg)

# Calculate physical volume
let vol = lo.physVol

#[ ~~~~ Functions for extra timing ~~~~ ]#

proc ticc(): float = 
   # Return t0
   result = cpuTime()

proc tocc(message: string, t0: float) =
   # Print timing
   echo message, " ", cpuTime() - t0

#[ ~~~~ Define functions setting boundary conditions ~~~~ ]#

#[ Set boundary conditions ]#
proc setBC_cust(gf: auto; bc: string) =
   # Cycle through BC
   for mu in 0..<gf.len:
      # Check if boundary is anti-periodic
      if ($bc[mu] == $"a"):
         # Cycle through coordinates
         tfor i, 0..<gf[mu].l.nSites:
            # Check if boundary
            if gf[mu].l.coords[mu][i] == gf[mu].l.physGeom[mu] - 1:
               # Multiply by -1
               gf[mu]{i} *= -1.0

#[ ~~~~ Define functions for staggered operations ~~~~ ]#

#[ Fermion solve ]#
proc solve_fermion(s: Staggered; x, b: auto; 
                   mass: float; sp0: var SolverParams) =
   #[ Start timers ]#

   # Get initial time
   let t0 = ticc()

   #[ Do solve ]#
   
   # Check if fermion is massless
   if mass != 0:
      # Do solve
      s.solve(x, b, mass, sp0)
   else:
      #[ Initialize CG parameters ]#

      # Create temporary solver params
      var sp = sp0

      # Reset stats
      sp.resetStats()

      # Set verbosity
      sp.verbosity = 1

      # Set previous solution variable
      sp.usePrevSoln = false

      #[ Do solve and reconstruct ]#

      # Start thread block
      threads:
         # Set x to zero
         x := 0

      # Do solve of just even part
      s.solveEE(x, b, 0, sp)

      # Set elapsed time
      sp.seconds = cpuTime() - t0

      # Start thread block
      threads:
         # Copy even 4X even sites
         x.even := 4*x

      #[ Take care of stats ]#

      # Set calls
      sp.calls = 1

      # Add stats to sp0
      sp0.addStats(sp)

#[ Apply massless D^{dagger} ]#
proc apply_massless_Ddag(s: Staggered; x, b: auto; option: string) =
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

         # Prevent race condition
         threadBarrier()

         # Set even sites to zero
         x.even := 0
      elif option == "force":
         # Set even sites to phi
         x.even := b

   # Print time taken
   tocc("D app.:", t0)

#[ ~~~~ Define functions for momentum and field generation ~~~~ ]#

#[ Generate momenta ]#
proc generate_momenta(p, gf: auto) =
   # Set initial time
   let t0 = ticc()

   # Start case (workaround)
   case r.rng_type:
      of "RngMilc6":
         # Start thread block
         threads:
            # Sample
            p.randomTAH(r.milc)
      of "MRG32k3a":
        # Start thread block
        threads:
            # Sample
            p.randomTAH(r.mrg32k3a)

   # Start thread block
   threads:
      # Cycle through all momenta
      for i in 0..<g.len:
         # Set initial gauge variable
         g0[i] := gf[i]

   # Print time taken
   tocc("Generate momenta and save backup gauge field:", t0)

#[ Generate pseudofermion fields ]#
proc generate_pseudoferms(s: Staggered) =
   # Set initial time
   let t0 = ticc()

   # Cycle through fermion fields
   for fld_ind in 0..<phi.len:
      # Start case  (work around)
      case r.rng_type:
         of "RngMilc6":
            # Start thread block
            threads:
               # Sample
               psi.gaussian(r.milc)
         of "MRG32k3a":
            # Start thread block
            threads:
               # Sample
               psi.gaussian(r.mrg32k3a)

      # Check if regular fermion field
      if fld_ind < int_prms["Nf"]:
         # Apply staggered Dirac operator (D^{d})
         s.Ddag(phi[fld_ind], psi, masses[0])
      else:
         # Apply inverted staggered Dirac operator
         s.solve(phi[fld_ind], psi, masses[1], spa[fld_ind])

      # Create thread block
      threads:
         # Set odd sites to zero
         phi[fld_ind].odd := 0

   # Print timing information
   tocc("Generate fermion/boson fields:", t0)

#[ ~~~~ Define functions for action calculation ~~~~ ]#

#[ For calculating gauge action ]#
proc gaction(gf: auto; f2: seq; p2: float): auto =
   # Set variables
   var
      # Set gauge action
      ga = gc.actionA(gf)

      # Set fermion action
      fa = sum(f2)

      # Get momentum from gauge field
      t = 0.5*p2 - float(16*vol)

   # Set hamiltonian
   let h = ga + fa + t

   # Return variables
   result = (ga, fa, t, h)

#[ Calculate action ]#
proc calc_action(s: Staggered; gf, p: auto): auto =
   #[ Initial timing info ]#

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
   var f2 = newSeq[float](phi.len)

   # Cycle through fermion fields
   for fld_ind in 0..<phi.len:
      # Start thread block
      threads:
         # Set psi to zero
         psi := 0

      # Check if regular fermion
      if fld_ind < int_prms["Nf"]:
         # Do fermion solve
         s.solve_fermion(psi, phi[fld_ind], -masses[0], spa[fld_ind])

         # Check if massless
         if masses[0] == 0:
            # Apply D^{d} and fill odd entries appropriately
            s.apply_massless_Ddag(psi, psi, "action")
      else:
         # Create thread block
         threads:
            # Apply D to phi to get psi for PV fermion
            s.D(psi, phi[fld_ind], masses[1])

      # Start thread block
      threads:
         # Increment f2 (reuses psi from HMC trajectory)
         let psi2 = psi.norm2()

         # Prevent race condition
         threadBarrier()

         # Increment f2
         threadMaster: f2[fld_ind] = 0.5 * psi2

   #[ Calculate gauge action and put everything together ]#

   # Calculate gauge action and return whole action
   let (ga, fa, t, h) = gf.gaction(f2, p2)

   #[ End timing information ]#
   
   # Print timing info
   tocc("Calculate action:", t0)

   # Return results for gauge action
   result = (ga, fa, f2, t, h)

#[ ~~~~ Define functions for force calculation ~~~~ ]#

#[ For smearing force ]#
proc smearRephase(gf, sgf: auto): proc =
   # Get initial time
   let t0 = ticc()

   # Get smeared force
   let smearedForce = coef.smearGetForce(gf, sgf, info)

   # Start thread block
   threads:
      # Set boundary conditions
      sgf.setBC_cust(str_prms["bc"])

      # Create thread barrier 
      threadBarrier()

      # Set staggered phases
      sgf.stagPhase
   
   # Print timing information
   tocc("Smear gauge and prepare smeared force:", t0)

   # Return force proc
   result = smearedForce

#[ For smearing gauge field ]#
proc smearRephaseDiscardForce(gf, sgf: auto) =
   # Get initial time
   let t0 = ticc()

   # Smear gauge field
   coef.smear(gf, sgf, info)

   # Start thread block
   threads:
      # Set boundary conditions
      sgf.setBC_cust(str_prms["bc"])

      # Set thread barrier
      threadBarrier()

      # Set staggered phases
      sgf.stagPhase

   # Get timing information
   tocc("Smear fields:", t0)

#[ Rescaling for different fields ]#
proc rescale(index: int, t: float): float =
   # Get value of s
   var s = -0.5 * t

   # Check type of fermion
   if (index < int_prms["Nf"]) and (masses[0] != 0):
      # Add factors appropriate to regular fermion
      s = s / masses[0]
   elif index < int_prms["Nf"]:
      # Add factors appropriate to massess fermion
      s = -0.5 * s
   else:
      # Multiply by appropriate factor
      s = 0.5 * s

   # Return rescale
   result = s

#[ Smeared one link force ]#
proc smeared_one_link_force(f: auto; gf: auto; smeared_force: proc) =
   #[ Correcting phase ]#
   
   # Get initial time
   let t0 = ticc()

   # Start thread block
   threads:
      # Set boundary conditions
      f.setBC_cust(str_prms["bc"])

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
            # Re-sign odd sites
            f[mu][i] *= -1

   # Print timing information
   tocc("Rephase force/set BC's:", t0)

   #[ Smear force ]#

   # Smear
   f.smeared_force(f)

   # Print timing information
   tocc("Smear fermion/boson force:", t0)

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
            s := f[mu][i] * gf[mu][i].adj

            # Project to traceless/anti-Hermitian component
            projectTAH(f[mu][i], s)

   # Print timing information
   tocc("Multiply force by gauge link and project TA:", t0)

#[ Calculate smeared force ]#
proc fforce(s: Staggered; f: auto; gf: auto;
            smeared_force: proc;
            ix: openarray[int]; ts: openarray[float]) =

   #[ Do Dslash ]#

   # Get initial time
   let t0 = ticc()

   # Create shifter
   var t: array[4, Shifter[typeof(psi), typeof(psi[0])]]

   # Cycle through directions
   for mu in 0..<f.len:
      # Set shifter
      t[mu] = newShifter(psi, mu, 1)

   # Define convenient variables
   var
      # Create variable controlling f updating behavior
      frc_ind = 0

      # Initialize scale (for proper normalization on force)
      scale = 0.0

      # Define update sequence (to tell the updater what fields to include)
      update_seq = newseq[int](0)

      # Define number of fields (for cycling through all fields)
      n_fields = int_prms["Nf"] + int_prms["num_pv"]

   # Cycle through indices
   for f_ind in 0..<n_fields:
      # Check if fermions to be updated
      if (f_ind < int_prms["Nf"]) and (0 in ix):
         # Add fermion field
         update_seq.add f_ind
      elif (f_ind >= int_prms["Nf"]) and (1 in ix):
         # Add pv field
         update_seq.add f_ind

   #[ Calculate outer product for Dslash ]#

   # Cycle through indices
   for f_ind in update_seq:
      #[ Solves or application of massive D ]#

      # Start thread block
      threads:
         # Set psi to zero
         psi := 0

      # Check if regular fermion (place where things can go wrong)
      if f_ind < int_prms["Nf"]:
         # Do solve
         s.solve_fermion(psi, phi[f_ind], masses[0], spf[f_ind])
 
         # Check if massless
         if masses[0] == 0:
            # Apply D^{d} and fill odd entries appropriately
            s.apply_massless_Ddag(psi, psi, "force")

         # Set scale
         scale = rescale(f_ind, ts[0])
      else:
         # Apply D^{d} and fill odd entries appropriately
         s.apply_massless_Ddag(psi, phi[f_ind], "force")

         # Set scale
         scale = rescale(f_ind, ts[1])

      #[ Application of massless D and creation of outer product ]#

      # Cycle through directions
      for mu in 0..<f.len:
         # Essentially apply staggered Dirac operator
         discard t[mu] ^* psi

      # Create variable for convenience
      let n = psi[0].len

      # Calculate outer product for Dslash
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
                     # Start case
                     case frc_ind:
                        # Either initialize outer product or append to it
                        of 0: f[mu][i][a, b] := scale * psi[i][a] * t[mu].field[i][b].adj
                        else: f[mu][i][a, b] += scale * psi[i][a] * t[mu].field[i][b].adj

      # Append force index
      frc_ind = frc_ind + 1

   # Print timing information
   tocc("Dslash and outer product", t0)

   #[ Smear and rephase ]#

   # Smear force
   f.smeared_one_link_force(gf, smeared_force)

   # Print timing information
   tocc("Full fermion/boson force calculation:", t0)

#[ ~~~~ Define functions HMC integration ~~~~ ]#

#[ Gauge field update from momenta ]#
proc mdt(t: float) =
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

#[ Momentum update from gauge sector ]#
proc mdv(t: float) =
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

#[ Momentum update w/ fermions/bosons ]#
proc mdvf(ix: openarray[int], ts: openarray[float], smeared_force: proc) =
   # Get initial time
   let t0 = ticc()

   # Calculate force
   stag.fforce(f, g, smeared_force, ix, ts)

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

#[ Integration with shared update ]#
proc mdvAllfga(ts: openarray[float]) =
   #[ Initial setup ]#
   
   # Get initial time
   let t0 = ticc()

   # First, define variables
   var
      # Option for updating gauge
      updateG = false

      # Option for updating fermions
      updateF = newseq[int](0)

   #[ Determine what fields are to be updated ]#   

   # Check if gauge field to be updated
   if ts[0] != 0:
      # Otherwise, just update gauge
      updateG = true
   
   # Check if pseudofermion/boson fields to be updated
   for k in 0..<ts.len - 1:
      # Define convenient index
      let i = k + 1

      # Check if field to be updated
      if ts[i] != 0:
         # Otherwise, just update fermions
         updateF.add k

   #[ Do molecular dynamics ]#

   # Check if gauge field to be updated
   if updateG:
      # Update gauge sector
      mdv(ts[0])

   # Check if fermions to updated
   if updateF.len > 0:
      # Create smeared sforce proc
      var smeared_force = g.smearRephase(sg)

      # Update fermion sector
      mdvf(updateF, ts[1..^1], smeared_force)

      # Set sforce proc to nil
      smeared_force = nil

   #[ Print timing information ]#

   # Print timing information
   tocc("Integrator update: ", t0)

#[ Set up integrator ]#
let
   # Define gauge integration algorithm
   gauge_int_alg: IntegratorProc = str_prms["gauge_int_alg"]

   # Define fermion integration algorithm
   ferm_int_alg: IntegratorProc = str_prms["ferm_int_alg"]

   # Define Pauli-Villars integration algorithm
   pv_int_alg: IntegratorProc = str_prms["pv_int_alg"]

   # Create integration pair
   (V, T) = newIntegratorPair(mdvAllfga, mdt)

   # Set integrator for gauge
   H = newParallelEvolution gauge_int_alg(T = T, V = V[0], steps = int_prms["g_steps"])
block:
   # Add fermions to be updated
   H.add ferm_int_alg(T = T, V = V[1], steps = int_prms["f_steps"])

   # Check if PV bosons to be updated
   H.add pv_int_alg(T = T, V = V[2], steps = int_prms["pv_steps"])

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
   # Cycle through fermions
   for fld_ind in 0..<phi.len:
      # Print stats for action solver
      checkStats("Solver " & intToStr(fld_ind) & " [action]: ", spa[fld_ind])

      # Check mass index
      if fld_ind < int_prms["Nf"]:
         # Print stats for force solver
         checkStats("Solver " & intToStr(fld_ind) & " [force]: ", spf[fld_ind])

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

      # Get initial time for trajectory
      let traj_init_time = ticc()

      # Define ratio string
      let rat_str = intToStr(traj + 1) & "/" & intToStr(config_space)

      # Tell user what trajectory that you're on
      echo "\n~~~~ Trajectory #: " & intToStr(config) & " " & rat_str & " ~~~~\n"

      #[ Generate field momenta and save ]#

      # Generate random momenta and save initial gauge field
      generate_momenta(p, g)

      #[ Smear fields and save ]#

      # Now smear and rephase
      g.smearRephaseDiscardForce sg

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

      # Check to make sure above user-set Metropolis threshold
      if (config + 1) >= int_prms["no_metropolis_until"]:
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
         # Create filename
         let fn = io_path & def_fn & "_" & intToStr(config + 1)
   
         # Write gauge and RNG field
         fn.write_fields(r, g)

         # Write global RNG
         R.write_rng(fn)

      #[ Do measurements ]#

      # Check if trajectory was accepted again
      if accr <= acc:
         # Reunitarize gauge field
         g.reunit

         # Smear gauge field
         g.smearRephaseDiscardForce sg   

      # Check if plaquette to be measured
      if (int_prms["plaq_freq"] > 0) and ((traj + 1) mod int_prms["plaq_freq"] == 0):
         # Measure plaquette
         g.mplaq

      # Check if Polyakov loop to be measured
      if (int_prms["ploop_freq"] > 0) and ((traj + 1) mod int_prms["ploop_freq"] == 0):
         # Measure polyakov loop
         g.ploop

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

# Print end date
echo "\nEnd: ", now().utc, "\n" 

# Finalize QEX
qexfinalize()