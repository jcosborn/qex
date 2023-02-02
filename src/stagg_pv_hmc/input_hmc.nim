<?xml version="1.0"?>
<!-- Credit to Oliver Witzel for formatting of XML file -->
<qex>
  <hmc><!-- Hybrid Monte Carlo parameters -->
    <tau>0.7071067811865475</tau>
    <g_steps>60</g_steps> <!--Gauge field steps -->
    <f_steps>35</f_steps> <!--Fermion field steps -->
    <pv_steps>35</pv_steps> <!--PV boson field steps -->
    <gauge_int_alg>2MN</gauge_int_alg> <!--Integration algorithm for gauge field -->
    <ferm_int_alg>2MN</ferm_int_alg> <!--Integration algorithm for fermion fields -->
    <pv_int_alg>2MN</pv_int_alg> <!--Integration algorithm for Pauli-Villars fields -->
    <no_metropolis_until>1</no_metropolis_until> <!-- When to start Metropolis test -->
  </hmc>
  <config_opts><!-- Options for gauge configs -->
    <start_config>137</start_config> <!-- Starting config. is reunitarized -->
    <start>unit</start> <!-- Cold start = unit, hot start = rand -->
  </config_opts>
  <rng><!-- Random number seeds -->
    <rng_type>RngMilc6</rng_type> <!-- RNG type: MRG32k3a or RngMilc6 -->
    <parallel_seed>987654321</parallel_seed> <!-- Seed for RNG fields -->
    <serial_seed>987654321</serial_seed> <!-- Seed for global RNG (Metrop. step) -->
  </rng>
  <action><!-- Action parameters -->
    <geom><!-- Specify geometry, Ns^3 X Nt -->
      <Ns>32</Ns>
      <Nt>32</Nt>
      <bc>aaaa</bc> <!-- Boundary conditions; e.g., pppa, aaaa, etc. -->
    </geom>
    <gauge><!-- Gauge action parameters (adjoint gauge action) -->
      <beta>9.2</beta>
      <adj_fac>-0.25</adj_fac>
    </gauge>
    <ferm><!-- Action parameters for fermions -->
      <Nf>2</Nf>
      <mass>0.0</mass>
    </ferm>
    <pv><!-- Action parameters for Pauli-Villars -->
      <num_pv>16</num_pv> <!-- num_pv = 0, then no Pauli-Villars bosons -->
      <mass_pv>0.75</mass_pv>
    </pv>
  </action>
  <nhyp_smearing><!-- Parameters for nHYP smearing -->
    <alpha_1>0.4</alpha_1>
    <alpha_2>0.5</alpha_2>
    <alpha_3>0.5</alpha_3>
  </nhyp_smearing>
  <solver><!-- Solver parameters for action, force and pbp -->
    <a_tol>1e-20</a_tol> <!-- Action solver tolerance -->
    <a_maxits>10000</a_maxits> <!-- Action solver maximum itns. -->
    <f_tol>1e-12</f_tol> <!-- Force solver " " -->
    <f_maxits>10000</f_maxits> <!-- Force solver " " -->
    <check_solvers>1</check_solvers> <!-- Frequency of checking solvers -->
  </solver>
  <extra><!-- Extra optional parameters -->
    <basic_meas> <!-- For basic measurement of observables -->
      <plaq_freq>1</plaq_freq> <!-- Plaquette; if zero, no measurement -->
      <ploop_freq>1</ploop_freq> <!-- Polyakov loop; if zero, no measurement -->
    </basic_meas>
    <hmc_checks> <!-- For checks of HMC -->
      <rev_check_freq>0</rev_check_freq> <!-- Frequency of reversibility check for integrator -->
    </hmc_checks>
  </extra>
</qex>
