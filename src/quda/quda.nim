import enum_quda, quda_constants
##
##  @file  quda.h
##  @brief Main header file for the QUDA library
##
##  Note to QUDA developers: When adding new members to QudaGaugeParam
##  and QudaInvertParam, be sure to update lib/check_params.h as well
##  as the Fortran interface in lib/quda_fortran.F90.
##

type
  ConstInt* {.importc:"const int".} = cint
  double_complex* {.importc:"double _Complex".} = object
converter toDoubleComplex*(x: array[2,float]): double_complex =
  var r = cast[ptr array[2,float]](addr result)
  r[] = x

##
##  Parameters having to do with the gauge field or the
##  interpretation of the gauge field by various Dirac operators
##

type
  INNER_C_UNION_quda_2* {.importc: "QudaGaugeParam::no_name", header: "quda.h",
                         bycopy, union.} = object
    use_split_gauge_bkup* {.importc: "use_split_gauge_bkup".}: bool
    ## < Used by gauge split buffers (default=true keep split gauge after usage)
    pad* {.importc: "pad".}: cint
    ## < Forces 4-byte alignment

  QudaGaugeParam* {.importc: "QudaGaugeParam", header: "quda.h", bycopy.} = object
    struct_size* {.importc: "struct_size".}: csize_t
    ## < Size of this struct in bytes.  Used to ensure that the host application and QUDA see the same struct size
    location* {.importc: "location".}: QudaFieldLocation
    ## < The location of the gauge field
    X* {.importc: "X".}: array[4, cint]
    ## < The local space-time dimensions (without checkboarding)
    anisotropy* {.importc: "anisotropy".}: cdouble
    ## < Used for Wilson and Wilson-clover
    tadpole_coeff* {.importc: "tadpole_coeff".}: cdouble
    ## < Used for staggered only
    scale* {.importc: "scale".}: cdouble
    ## < Used by staggered long links
    `type`* {.importc: "type".}: QudaLinkType
    ## < The link type of the gauge field (e.g., Wilson, fat, long, etc.)
    gauge_order* {.importc: "gauge_order".}: QudaGaugeFieldOrder
    ## < The ordering on the input gauge field
    t_boundary* {.importc: "t_boundary".}: QudaTboundary
    ## < The temporal boundary condition that will be used for fermion fields
    cpu_prec* {.importc: "cpu_prec".}: QudaPrecision
    ## < The precision used by the caller
    cuda_prec* {.importc: "cuda_prec".}: QudaPrecision
    ## < The precision of the cuda gauge field
    reconstruct* {.importc: "reconstruct".}: QudaReconstructType
    ## < The reconstruction type of the cuda gauge field
    cuda_prec_sloppy* {.importc: "cuda_prec_sloppy".}: QudaPrecision
    ## < The precision of the sloppy gauge field
    reconstruct_sloppy* {.importc: "reconstruct_sloppy".}: QudaReconstructType
    ## < The recontruction type of the sloppy gauge field
    cuda_prec_refinement_sloppy* {.importc: "cuda_prec_refinement_sloppy".}: QudaPrecision
    ## < The precision of the sloppy gauge field for the refinement step in multishift
    reconstruct_refinement_sloppy* {.importc: "reconstruct_refinement_sloppy".}: QudaReconstructType
    ## < The recontruction type of the sloppy gauge field for the refinement step in multishift
    cuda_prec_precondition* {.importc: "cuda_prec_precondition".}: QudaPrecision
    ## < The precision of the preconditioner gauge field
    reconstruct_precondition* {.importc: "reconstruct_precondition".}: QudaReconstructType
    ## < The recontruction type of the preconditioner gauge field
    cuda_prec_eigensolver* {.importc: "cuda_prec_eigensolver".}: QudaPrecision
    ## < The precision of the eigensolver gauge field
    reconstruct_eigensolver* {.importc: "reconstruct_eigensolver".}: QudaReconstructType
    ## < The recontruction type of the eigensolver gauge field
    gauge_fix* {.importc: "gauge_fix".}: QudaGaugeFixed
    ## < Whether the input gauge field is in the axial gauge or not
    ga_pad* {.importc: "ga_pad".}: cint
    ## < The pad size that native GaugeFields will use (default=0)
    site_ga_pad* {.importc: "site_ga_pad".}: cint
    ## < Used by link fattening and the gauge and fermion forces
    staple_pad* {.importc: "staple_pad".}: cint
    ## < Used by link fattening
    llfat_ga_pad* {.importc: "llfat_ga_pad".}: cint
    ## < Used by link fattening
    mom_ga_pad* {.importc: "mom_ga_pad".}: cint
    ## < Used by the gauge and fermion forces
    ano_quda_3* {.importc: "ano_quda_3".}: INNER_C_UNION_quda_2
    staggered_phase_type* {.importc: "staggered_phase_type".}: QudaStaggeredPhase
    ## < Set the staggered phase type of the links
    staggered_phase_applied* {.importc: "staggered_phase_applied".}: cint
    ## < Whether the staggered phase has already been applied to the links
    i_mu* {.importc: "i_mu".}: cdouble
    ## < Imaginary chemical potential
    overlap* {.importc: "overlap".}: cint
    ## < Width of overlapping domains
    overwrite_gauge* {.importc: "overwrite_gauge".}: cint
    ## < When computing gauge, should we overwrite it or accumulate to it
    overwrite_mom* {.importc: "overwrite_mom".}: cint
    ## < When computing momentum, should we overwrite it or accumulate to it
    use_resident_gauge* {.importc: "use_resident_gauge".}: cint
    ## < Use the resident gauge field as input
    use_resident_mom* {.importc: "use_resident_mom".}: cint
    ## < Use the resident momentum field as input
    make_resident_gauge* {.importc: "make_resident_gauge".}: cint
    ## < Make the result gauge field resident
    make_resident_mom* {.importc: "make_resident_mom".}: cint
    ## < Make the result momentum field resident
    return_result_gauge* {.importc: "return_result_gauge".}: cint
    ## < Return the result gauge field
    return_result_mom* {.importc: "return_result_mom".}: cint
    ## < Return the result momentum field
    gauge_offset* {.importc: "gauge_offset".}: csize_t
    ## < Offset into MILC site struct to the gauge field (only if gauge_order=MILC_SITE_GAUGE_ORDER)
    mom_offset* {.importc: "mom_offset".}: csize_t
    ## < Offset into MILC site struct to the momentum field (only if gauge_order=MILC_SITE_GAUGE_ORDER)
    site_size* {.importc: "site_size".}: csize_t
    ## < Size of MILC site struct (only if gauge_order=MILC_SITE_GAUGE_ORDER)


##
##  Parameters relating to the solver and the choice of Dirac operator.
##

type
  QudaInvertParam* {.importc: "QudaInvertParam", header: "quda.h", bycopy.} = object
    ##  Size of this struct in bytes.  Used to ensure that the host application and QUDA see the same struct size
    struct_size* {.importc: "struct_size".}: csize_t
    input_location* {.importc: "input_location".}: QudaFieldLocation
    ## < The location of the input field
    output_location* {.importc: "output_location".}: QudaFieldLocation
    ## < The location of the output field
    dslash_type* {.importc: "dslash_type".}: QudaDslashType
    ## < The Dirac Dslash type that is being used
    inv_type* {.importc: "inv_type".}: QudaInverterType
    ## < Which linear solver to use
    mass* {.importc: "mass".}: cdouble
    ## < Used for staggered only
    kappa* {.importc: "kappa".}: cdouble
    ## < Used for Wilson and Wilson-clover
    m5* {.importc: "m5".}: cdouble
    ## < Domain wall height
    Ls* {.importc: "Ls".}: cint
    ## < Extent of the 5th dimension (for domain wall)
    b_5* {.importc: "b_5".}: array[QUDA_MAX_DWF_LS, double_complex]
    ## < Mobius coefficients - only real part used if regular Mobius
    c_5* {.importc: "c_5".}: array[QUDA_MAX_DWF_LS, double_complex]
    ## < Mobius coefficients - only real part used if regular Mobius
    ## <
    ##  The following specifies the EOFA parameters. Notation follows arXiv:1706.05843
    ##  eofa_shift: the "\beta" in the paper
    ##  eofa_pm: plus or minus for the EOFA operator
    ##  mq1, mq2, mq3 are the three masses corresponds to Hasenbusch mass spliting.
    ##  As far as I know mq1 is always the same as "mass" but it's here just for consistence.
    ##
    eofa_shift* {.importc: "eofa_shift".}: cdouble
    eofa_pm* {.importc: "eofa_pm".}: cint
    mq1* {.importc: "mq1".}: cdouble
    mq2* {.importc: "mq2".}: cdouble
    mq3* {.importc: "mq3".}: cdouble
    mu* {.importc: "mu".}: cdouble
    ## < Twisted mass parameter
    tm_rho* {.importc: "tm_rho".}: cdouble
    ## < Hasenbusch mass shift applied like twisted mass to diagonal (but not inverse)
    epsilon* {.importc: "epsilon".}: cdouble
    ## < Twisted mass parameter
    evmax* {.importc: "evmax".}: cdouble
    ##  maximum of the eigenvalues of the ndeg twisted mass operator needed for fermionic forces  *
    twist_flavor* {.importc: "twist_flavor".}: QudaTwistFlavorType
    ## < Twisted mass flavor
    laplace3D* {.importc: "laplace3D".}: cint
    ## < omit this direction from laplace operator: x,y,z,t -> 0,1,2,3 (-1 is full 4D)
    covdev_mu* {.importc: "covdev_mu".}: cint
    ## < Apply forward/backward covariant derivative in direction mu(mu<=3)/mu-4(mu>3)
    tol* {.importc: "tol".}: cdouble
    ## < Solver tolerance in the L2 residual norm
    tol_restart* {.importc: "tol_restart".}: cdouble
    ## < Solver tolerance in the L2 residual norm (used to restart InitCG)
    tol_hq* {.importc: "tol_hq".}: cdouble
    ## < Solver tolerance in the heavy quark residual norm
    compute_true_res* {.importc: "compute_true_res".}: cint
    ##  Whether to compute the true residual post solve
    true_res* {.importc: "true_res".}: array[QUDA_MAX_MULTI_SRC, cdouble]
    ## < Actual L2 residual norm achieved in the solver
    true_res_hq* {.importc: "true_res_hq".}: array[QUDA_MAX_MULTI_SRC, cdouble]
    ## < Actual heavy quark residual norm achieved in the solver
    maxiter* {.importc: "maxiter".}: cint
    ## < Maximum number of iterations in the linear solver
    reliable_delta* {.importc: "reliable_delta".}: cdouble
    ## < Reliable update tolerance
    reliable_delta_refinement* {.importc: "reliable_delta_refinement".}: cdouble
    ## < Reliable update tolerance used in post multi-shift solver refinement
    use_alternative_reliable* {.importc: "use_alternative_reliable".}: cint
    ## < Whether to use alternative reliable updates
    use_sloppy_partial_accumulator* {.importc: "use_sloppy_partial_accumulator".}: cint
    ## < Whether to keep the partial solution accumuator in sloppy precision
    ## < This parameter determines how often we accumulate into the
    ##        solution vector from the direction vectors in the solver.
    ##        E.g., running with solution_accumulator_pipeline = 4, means we
    ##        will update the solution vector every four iterations using the
    ##        direction vectors from the prior four iterations.  This
    ##        increases performance of mixed-precision solvers since it means
    ##        less high-precision vector round-trip memory travel, but
    ##        requires more low-precision memory allocation.
    solution_accumulator_pipeline* {.importc: "solution_accumulator_pipeline".}: cint
    ## < This parameter determines how many consecutive reliable update
    ##     residual increases we tolerate before terminating the solver,
    ##     i.e., how long do we want to keep trying to converge
    max_res_increase* {.importc: "max_res_increase".}: cint
    ## < This parameter determines how many total reliable update
    ##     residual increases we tolerate before terminating the solver,
    ##     i.e., how long do we want to keep trying to converge
    max_res_increase_total* {.importc: "max_res_increase_total".}: cint
    ## < This parameter determines how many consecutive heavy-quark
    ##     residual increases we tolerate before terminating the solver,
    ##     i.e., how long do we want to keep trying to converge
    max_hq_res_increase* {.importc: "max_hq_res_increase".}: cint
    ## < This parameter determines how many total heavy-quark residual
    ##     restarts we tolerate before terminating the solver, i.e., how long
    ##     do we want to keep trying to converge
    max_hq_res_restart_total* {.importc: "max_hq_res_restart_total".}: cint
    ## < After how many iterations shall the heavy quark residual be updated
    heavy_quark_check* {.importc: "heavy_quark_check".}: cint
    pipeline* {.importc: "pipeline".}: cint
    ## < Whether to use a pipelined solver with less global sums
    num_offset* {.importc: "num_offset".}: cint
    ## < Number of offsets in the multi-shift solver
    num_src* {.importc: "num_src".}: cint
    ## < Number of sources in the multiple source solver
    num_src_per_sub_partition* {.importc: "num_src_per_sub_partition".}: cint
    ## < Number of sources in the multiple source solver, but per sub-partition
    ## < The grid of sub-partition according to which the processor grid will be partitioned.
    ##     Should have:
    ##       split_grid[0] * split_grid[1] * split_grid[2] * split_grid[3] * num_src_per_sub_partition == num_src. *
    split_grid* {.importc: "split_grid".}: array[QUDA_MAX_DIM, cint]
    overlap* {.importc: "overlap".}: cint
    ## < Width of domain overlaps
    ##  Offsets for multi-shift solver
    offset* {.importc: "offset".}: array[QUDA_MAX_MULTI_SHIFT, cdouble]
    ##  Solver tolerance for each offset
    tol_offset* {.importc: "tol_offset".}: array[QUDA_MAX_MULTI_SHIFT, cdouble]
    ##  Solver tolerance for each shift when refinement is applied using the heavy-quark residual
    tol_hq_offset* {.importc: "tol_hq_offset".}: array[QUDA_MAX_MULTI_SHIFT, cdouble]
    ##  Actual L2 residual norm achieved in solver for each offset
    true_res_offset* {.importc: "true_res_offset".}: array[QUDA_MAX_MULTI_SHIFT,
        cdouble]
    ##  Iterated L2 residual norm achieved in multi shift solver for each offset
    iter_res_offset* {.importc: "iter_res_offset".}: array[QUDA_MAX_MULTI_SHIFT,
        cdouble]
    ##  Actual heavy quark residual norm achieved in solver for each offset
    true_res_hq_offset* {.importc: "true_res_hq_offset".}: array[
        QUDA_MAX_MULTI_SHIFT, cdouble]
    ##  Residuals in the partial faction expansion
    residue* {.importc: "residue".}: array[QUDA_MAX_MULTI_SHIFT, cdouble]
    ##  Whether we should evaluate the action after the linear solver
    compute_action* {.importc: "compute_action".}: cint
    ##  Computed value of the bilinear action (complex-valued)
    ## 	invert: \phi^\dagger A^{-1} \phi
    ## 	multishift: \phi^\dagger r(x) \phi = \phi^\dagger (sum_k residue[k] * (A + offset[k])^{-1} ) \phi
    action* {.importc: "action".}: array[2, cdouble]
    solution_type* {.importc: "solution_type".}: QudaSolutionType
    ## < Type of system to solve
    solve_type* {.importc: "solve_type".}: QudaSolveType
    ## < How to solve it
    matpc_type* {.importc: "matpc_type".}: QudaMatPCType
    ## < The preconditioned matrix type
    dagger* {.importc: "dagger".}: QudaDagType
    ## < Whether we are using the Hermitian conjugate system or not
    mass_normalization* {.importc: "mass_normalization".}: QudaMassNormalizationT
    ## < The mass normalization is being used by the caller
    solver_normalization* {.importc: "solver_normalization".}: QudaSolverNormalization
    ## < The normalization desired in the solver
    preserve_source* {.importc: "preserve_source".}: QudaPreserveSource
    ## < Preserve the source or not in the linear solver (deprecated)
    cpu_prec* {.importc: "cpu_prec".}: QudaPrecision
    ## < The precision used by the input fermion fields
    cuda_prec* {.importc: "cuda_prec".}: QudaPrecision
    ## < The precision used by the QUDA solver
    cuda_prec_sloppy* {.importc: "cuda_prec_sloppy".}: QudaPrecision
    ## < The precision used by the QUDA sloppy operator
    cuda_prec_refinement_sloppy* {.importc: "cuda_prec_refinement_sloppy".}: QudaPrecision
    ## < The precision of the sloppy gauge field for the refinement step in multishift
    cuda_prec_precondition* {.importc: "cuda_prec_precondition".}: QudaPrecision
    ## < The precision used by the QUDA preconditioner
    cuda_prec_eigensolver* {.importc: "cuda_prec_eigensolver".}: QudaPrecision
    ## < The precision used by the QUDA eigensolver
    dirac_order* {.importc: "dirac_order".}: QudaDiracFieldOrder
    ## < The order of the input and output fermion fields
    gamma_basis* {.importc: "gamma_basis".}: QudaGammaBasis
    ## < Gamma basis of the input and output host fields
    clover_location* {.importc: "clover_location".}: QudaFieldLocation
    ## < The location of the clover field
    clover_cpu_prec* {.importc: "clover_cpu_prec".}: QudaPrecision
    ## < The precision used for the input clover field
    clover_cuda_prec* {.importc: "clover_cuda_prec".}: QudaPrecision
    ## < The precision used for the clover field in the QUDA solver
    clover_cuda_prec_sloppy* {.importc: "clover_cuda_prec_sloppy".}: QudaPrecision
    ## < The precision used for the clover field in the QUDA sloppy operator
    clover_cuda_prec_refinement_sloppy* {.
        importc: "clover_cuda_prec_refinement_sloppy".}: QudaPrecision
    ## < The precision of the sloppy clover field for the refinement step in multishift
    clover_cuda_prec_precondition* {.importc: "clover_cuda_prec_precondition".}: QudaPrecision
    ## < The precision used for the clover field in the QUDA preconditioner
    clover_cuda_prec_eigensolver* {.importc: "clover_cuda_prec_eigensolver".}: QudaPrecision
    ## < The precision used for the clover field in the QUDA eigensolver
    clover_order* {.importc: "clover_order".}: QudaCloverFieldOrder
    ## < The order of the input clover field
    use_init_guess* {.importc: "use_init_guess".}: QudaUseInitGuess
    ## < Whether to use an initial guess in the solver or not
    clover_csw* {.importc: "clover_csw".}: cdouble
    ## < Csw coefficient of the clover term
    clover_coeff* {.importc: "clover_coeff".}: cdouble
    ## < Coefficient of the clover term
    clover_rho* {.importc: "clover_rho".}: cdouble
    ## < Real number added to the clover diagonal (not to inverse)
    compute_clover_trlog* {.importc: "compute_clover_trlog".}: cint
    ## < Whether to compute the trace log of the clover term
    trlogA* {.importc: "trlogA".}: array[2, cdouble]
    ## < The trace log of the clover term (even/odd computed separately)
    compute_clover* {.importc: "compute_clover".}: cint
    ## < Whether to compute the clover field
    compute_clover_inverse* {.importc: "compute_clover_inverse".}: cint
    ## < Whether to compute the clover inverse field
    return_clover* {.importc: "return_clover".}: cint
    ## < Whether to copy back the clover matrix field
    return_clover_inverse* {.importc: "return_clover_inverse".}: cint
    ## < Whether to copy back the inverted clover matrix field
    verbosity* {.importc: "verbosity".}: QudaVerbosity
    ## < The verbosity setting to use in the solver
    iter* {.importc: "iter".}: cint
    ## < The number of iterations performed by the solver
    gflops* {.importc: "gflops".}: cdouble
    ## < The Gflops rate of the solver
    secs* {.importc: "secs".}: cdouble
    ## < The time taken by the solver
    energy* {.importc: "energy".}: cdouble
    ## < The energy consumed by the solver
    power* {.importc: "power".}: cdouble
    ## < The mean power of the solver
    temp* {.importc: "temp".}: cdouble
    ## < The mean temperature of the device for the duration of the solve
    clock* {.importc: "clock".}: cdouble
    ## < The mean clock frequency of the device for the duration of the solve
    ##  Number of steps in s-step algorithms
    Nsteps* {.importc: "Nsteps".}: cint
    ##  Maximum size of Krylov space used by solver
    gcrNkrylov* {.importc: "gcrNkrylov".}: cint
    ##
    ##  The following parameters are related to the solver
    ##  preconditioner, if enabled.
    ##
    ##
    ##  The inner Krylov solver used in the preconditioner.  Set to
    ##  QUDA_INVALID_INVERTER to disable the preconditioner entirely.
    ##
    inv_type_precondition* {.importc: "inv_type_precondition".}: QudaInverterType
    ##  Preconditioner instance, e.g., multigrid
    preconditioner* {.importc: "preconditioner".}: pointer
    ##  Deflation instance
    deflation_op* {.importc: "deflation_op".}: pointer
    ##  defines deflation
    eig_param* {.importc: "eig_param".}: pointer
    ##  If true, deflate the initial guess
    deflate* {.importc: "deflate".}: QudaBoolean
    ##  Dirac Dslash used in preconditioner
    dslash_type_precondition* {.importc: "dslash_type_precondition".}: QudaDslashType
    ##  Verbosity of the inner Krylov solver
    verbosity_precondition* {.importc: "verbosity_precondition".}: QudaVerbosity
    ##  Tolerance in the inner solver
    tol_precondition* {.importc: "tol_precondition".}: cdouble
    ##  Maximum number of iterations allowed in the inner solver
    maxiter_precondition* {.importc: "maxiter_precondition".}: cint
    ##  Relaxation parameter used in GCR-DD (default = 1.0)
    omega* {.importc: "omega".}: cdouble
    ##  Basis for CA algorithms
    ca_basis* {.importc: "ca_basis".}: QudaCABasis
    ##  Minimum eigenvalue for Chebyshev CA basis
    ca_lambda_min* {.importc: "ca_lambda_min".}: cdouble
    ##  Maximum eigenvalue for Chebyshev CA basis
    ca_lambda_max* {.importc: "ca_lambda_max".}: cdouble
    ##  Basis for CA algorithms in a preconditioned solver
    ca_basis_precondition* {.importc: "ca_basis_precondition".}: QudaCABasis
    ##  Minimum eigenvalue for Chebyshev CA basis in a preconditioner solver
    ca_lambda_min_precondition* {.importc: "ca_lambda_min_precondition".}: cdouble
    ##  Maximum eigenvalue for Chebyshev CA basis in a preconditioner solver
    ca_lambda_max_precondition* {.importc: "ca_lambda_max_precondition".}: cdouble
    ##  Number of preconditioner cycles to perform per iteration
    precondition_cycle* {.importc: "precondition_cycle".}: cint
    ##  Whether to use additive or multiplicative Schwarz preconditioning
    schwarz_type* {.importc: "schwarz_type".}: QudaSchwarzType
    ##  The type of accelerator type to use for preconditioner
    accelerator_type_precondition* {.importc: "accelerator_type_precondition".}: QudaAcceleratorType
    ##
    ##  The following parameters are the ones used to perform the adaptive MADWF in MSPCG
    ##  See section 3.3 of [arXiv:2104.05615]
    ##
    ##  The diagonal constant to suppress the low modes when performing 5D transfer
    madwf_diagonal_suppressor* {.importc: "madwf_diagonal_suppressor".}: cdouble
    ##  The target MADWF Ls to be used in the accelerator
    madwf_ls* {.importc: "madwf_ls".}: cint
    ##  The minimum number of iterations after which to generate the null vectors for MADWF
    madwf_null_miniter* {.importc: "madwf_null_miniter".}: cint
    ##  The maximum tolerance after which to generate the null vectors for MADWF
    madwf_null_tol* {.importc: "madwf_null_tol".}: cdouble
    ##  The maximum number of iterations for the training iterations
    madwf_train_maxiter* {.importc: "madwf_train_maxiter".}: cint
    ##  Whether to load the MADWF parameters from the file system
    madwf_param_load* {.importc: "madwf_param_load".}: QudaBoolean
    ##  Whether to save the MADWF parameters to the file system
    madwf_param_save* {.importc: "madwf_param_save".}: QudaBoolean
    ##  Path to load from the file system
    madwf_param_infile* {.importc: "madwf_param_infile".}: array[256, char]
    ##  Path to save to the file system
    madwf_param_outfile* {.importc: "madwf_param_outfile".}: array[256, char]
    ##
    ##  Whether to use the L2 relative residual, Fermilab heavy-quark
    ##  residual, or both to determine convergence.  To require that both
    ##  stopping conditions are satisfied, use a bitwise OR as follows:
    ##
    ##  p.residual_type = (QudaResidualType) (QUDA_L2_RELATIVE_RESIDUAL
    ##                                      | QUDA_HEAVY_QUARK_RESIDUAL);
    ##
    residual_type* {.importc: "residual_type".}: QudaResidualType
    ## Parameters for deflated solvers
    ##  The precision of the Ritz vectors
    cuda_prec_ritz* {.importc: "cuda_prec_ritz".}: QudaPrecision
    ##  How many vectors to compute after one solve
    ##   for eigCG recommended values 8 or 16
    ##
    n_ev* {.importc: "n_ev".}: cint
    ##  EeigCG  : Search space dimension
    ##   gmresdr : Krylov subspace dimension
    ##
    max_search_dim* {.importc: "max_search_dim".}: cint
    ##  For systems with many RHS: current RHS index
    rhs_idx* {.importc: "rhs_idx".}: cint
    ##  Specifies deflation space volume: total number of eigenvectors is n_ev*deflation_grid
    deflation_grid* {.importc: "deflation_grid".}: cint
    ##  eigCG: selection criterion for the reduced eigenvector set
    eigenval_tol* {.importc: "eigenval_tol".}: cdouble
    ##  mixed precision eigCG tuning parameter:  minimum search vector space restarts
    eigcg_max_restarts* {.importc: "eigcg_max_restarts".}: cint
    ##  initCG tuning parameter:  maximum restarts
    max_restart_num* {.importc: "max_restart_num".}: cint
    ##  initCG tuning parameter:  tolerance for cg refinement corrections in the deflation stage
    inc_tol* {.importc: "inc_tol".}: cdouble
    ##  Whether to make the solution vector(s) after the solve
    make_resident_solution* {.importc: "make_resident_solution".}: cint
    ##  Whether to use the resident solution vector(s)
    use_resident_solution* {.importc: "use_resident_solution".}: cint
    ##  Whether to use the solution vector to augment the chronological basis
    chrono_make_resident* {.importc: "chrono_make_resident".}: cint
    ##  Whether the solution should replace the last entry in the chronology
    chrono_replace_last* {.importc: "chrono_replace_last".}: cint
    ##  Whether to use the resident chronological basis
    chrono_use_resident* {.importc: "chrono_use_resident".}: cint
    ##  The maximum length of the chronological history to store
    chrono_max_dim* {.importc: "chrono_max_dim".}: cint
    ##  The index to indicate which chrono history we are augmenting
    chrono_index* {.importc: "chrono_index".}: cint
    ##  Precision to store the chronological basis in
    chrono_precision* {.importc: "chrono_precision".}: QudaPrecision
    ##  Which external library to use in the linear solvers (Eigen)
    extlib_type* {.importc: "extlib_type".}: QudaExtLibType
    ##  Whether to use the platform native or generic BLAS / LAPACK
    native_blas_lapack* {.importc: "native_blas_lapack".}: QudaBoolean
    ##  Whether to use fused kernels for mobius
    use_mobius_fused_kernel* {.importc: "use_mobius_fused_kernel".}: QudaBoolean
    ##
    ##  Parameters for distance preconditioning algorithm proposed in arXiv:1006.4028,
    ##  which is useful to solve a precise heavy quark propagator.
    ##  alpha0 and t0 follow Eq.(9) in the article.
    ##
    ##  The alpha0 parameter for distance preconditioning, related to the pseudoscalar meson mass
    distance_pc_alpha0* {.importc: "distance_pc_alpha0".}: cdouble
    ##  The t0 parameter for distance preconditioning, the timeslice where the source is located
    distance_pc_t0* {.importc: "distance_pc_t0".}: cint


##  Parameter set for solving eigenvalue problems.

type
  QudaEigParam* {.importc: "QudaEigParam", header: "quda.h", bycopy.} = object
    ##  Size of this struct in bytes.  Used to ensure that the host application and QUDA see the same struct size
    struct_size* {.importc: "struct_size".}: csize_t
    ##  EIGENSOLVER PARAMS
    ## -------------------------------------------------
    ##  Used to store information pertinent to the operator *
    invert_param* {.importc: "invert_param".}: ptr QudaInvertParam
    ##  Type of eigensolver algorithm to employ *
    eig_type* {.importc: "eig_type".}: QudaEigType
    ##  Use Polynomial Acceleration *
    use_poly_acc* {.importc: "use_poly_acc".}: QudaBoolean
    ##  Degree of the Chebysev polynomial *
    poly_deg* {.importc: "poly_deg".}: cint
    ##  Range used in polynomial acceleration *
    a_min* {.importc: "a_min".}: cdouble
    a_max* {.importc: "a_max".}: cdouble
    ##  Whether to preserve the deflation space between solves.  If
    ##         true, the space will be stored in an instance of the
    ##         deflation_space struct, pointed to by preserve_deflation_space
    preserve_deflation* {.importc: "preserve_deflation".}: QudaBoolean
    ##  This is where we store the deflation space.  This will point
    ##         to an instance of deflation_space. When a deflated solver is enabled, the deflation space will be obtained from this.
    preserve_deflation_space* {.importc: "preserve_deflation_space".}: pointer
    ##  If we restore the deflation space, this boolean indicates
    ##         whether we are also preserving the evalues or recomputing
    ##         them.  For example if a different mass shift is being used
    ##         than the one used to generate the space, then this should be
    ##         false, but preserve_deflation would be true
    preserve_evals* {.importc: "preserve_evals".}: QudaBoolean
    ##  Whether to use the smeared gauge field for the Dirac operator
    ##         for whose eigenvalues are are computing.
    use_smeared_gauge* {.importc: "use_smeared_gauge".}: bool
    ##  What type of Dirac operator we are using *
    ##  If !(use_norm_op) && !(use_dagger) use M. *
    ##  If use_dagger, use Mdag *
    ##  If use_norm_op, use MdagM *
    ##  If use_norm_op && use_dagger use MMdag. *
    ##  If use_pc for any, then use the even-odd pc version *
    use_dagger* {.importc: "use_dagger".}: QudaBoolean
    use_norm_op* {.importc: "use_norm_op".}: QudaBoolean
    use_pc* {.importc: "use_pc".}: QudaBoolean
    ##  Use Eigen routines to eigensolve the upper Hessenberg via QR *
    use_eigen_qr* {.importc: "use_eigen_qr".}: QudaBoolean
    ##  Performs an MdagM solve, then constructs the left and right SVD. *
    compute_svd* {.importc: "compute_svd".}: QudaBoolean
    ##  Performs the \gamma_5 OP solve by Post multipling the eignvectors with
    ##         \gamma_5 before computing the eigenvalues
    compute_gamma5* {.importc: "compute_gamma5".}: QudaBoolean
    ##  If true, the solver will error out if the convergence criteria are not met *
    require_convergence* {.importc: "require_convergence".}: QudaBoolean
    ##  Which part of the spectrum to solve *
    spectrum* {.importc: "spectrum".}: QudaEigSpectrumType
    ##  Size of the eigenvector search space *
    n_ev* {.importc: "n_ev".}: cint
    ##  Total size of Krylov space *
    n_kr* {.importc: "n_kr".}: cint
    ##  Max number of locked eigenpairs (deduced at runtime) *
    nLockedMax* {.importc: "nLockedMax".}: cint
    ##  Number of requested converged eigenvectors *
    n_conv* {.importc: "n_conv".}: cint
    ##  Number of requested converged eigenvectors to use in deflation *
    n_ev_deflate* {.importc: "n_ev_deflate".}: cint
    ##  Tolerance on the least well known eigenvalue's residual *
    tol* {.importc: "tol".}: cdouble
    ##  Tolerance on the QR iteration *
    qr_tol* {.importc: "qr_tol".}: cdouble
    ##  For IRLM/IRAM, check every nth restart *
    check_interval* {.importc: "check_interval".}: cint
    ##  For IRLM/IRAM, quit after n restarts *
    max_restarts* {.importc: "max_restarts".}: cint
    ##  For the Ritz rotation, the maximal number of extra vectors the solver may allocate *
    batched_rotate* {.importc: "batched_rotate".}: cint
    ##  For block method solvers, the block size *
    block_size* {.importc: "block_size".}: cint
    ##  The batch size used when computing eigenvalues *
    compute_evals_batch_size* {.importc: "compute_evals_batch_size".}: cint
    ##  For block method solvers, quit after n attempts at block orthonormalisation *
    max_ortho_attempts* {.importc: "max_ortho_attempts".}: cint
    ##  For hybrid modifeld Gram-Schmidt orthonormalisations *
    ortho_block_size* {.importc: "ortho_block_size".}: cint
    ##  In the test function, cross check the device result against ARPACK *
    arpack_check* {.importc: "arpack_check".}: QudaBoolean
    ##  For Arpack cross check, name of the Arpack logfile *
    arpack_logfile* {.importc: "arpack_logfile".}: array[512, char]
    ##  Name of the QUDA logfile (residua, upper Hessenberg/tridiag matrix updates) *
    QUDA_logfile* {.importc: "QUDA_logfile".}: array[512, char]
    ##  The orthogonal direction in the 3D eigensolver *
    ortho_dim* {.importc: "ortho_dim".}: cint
    ##  The size of the orthogonal direction in the 3D eigensolver, local *
    ortho_dim_size_local* {.importc: "ortho_dim_size_local".}: cint
    ## -------------------------------------------------
    ##  EIG-CG PARAMS
    ## -------------------------------------------------
    nk* {.importc: "nk".}: cint
    np* {.importc: "np".}: cint
    ##  Whether to load eigenvectors
    import_vectors* {.importc: "import_vectors".}: QudaBoolean
    ##  The precision of the Ritz vectors
    cuda_prec_ritz* {.importc: "cuda_prec_ritz".}: QudaPrecision
    ##  The memory type used to keep the Ritz vectors
    mem_type_ritz* {.importc: "mem_type_ritz".}: QudaMemoryType
    ##  Location where deflation should be done
    location* {.importc: "location".}: QudaFieldLocation
    ##  Whether to run the verification checks once set up is complete
    run_verify* {.importc: "run_verify".}: QudaBoolean
    ##  Filename prefix where to load the null-space vectors
    vec_infile* {.importc: "vec_infile".}: array[256, char]
    ##  Filename prefix for where to save the null-space vectors
    vec_outfile* {.importc: "vec_outfile".}: array[256, char]
    ##  The precision with which to save the vectors
    save_prec* {.importc: "save_prec".}: QudaPrecision
    ##  Whether to inflate single-parity eigen-vector I/O to a full
    ##         field (e.g., enabling this is required for compatability with
    ##         MILC I/O)
    io_parity_inflate* {.importc: "io_parity_inflate".}: QudaBoolean
    ##  Whether to save eigenvectors in QIO singlefile or partfile format
    partfile* {.importc: "partfile".}: QudaBoolean
    ##  Which external library to use in the deflation operations (Eigen)
    extlib_type* {.importc: "extlib_type".}: QudaExtLibType
    ## -------------------------------------------------

  QudaMultigridParam* {.importc: "QudaMultigridParam", header: "quda.h", bycopy.} = object
    ##  Size of this struct in bytes.  Used to ensure that the host application and QUDA see the same struct size
    struct_size* {.importc: "struct_size".}: csize_t
    invert_param* {.importc: "invert_param".}: ptr QudaInvertParam
    eig_param* {.importc: "eig_param".}: array[QUDA_MAX_MG_LEVEL, ptr QudaEigParam]
    ##  Number of multigrid levels
    n_level* {.importc: "n_level".}: cint
    ##  Geometric block sizes to use on each level
    geo_block_size* {.importc: "geo_block_size".}: array[QUDA_MAX_MG_LEVEL,
        array[QUDA_MAX_DIM, cint]]
    ##  Spin block sizes to use on each level
    spin_block_size* {.importc: "spin_block_size".}: array[QUDA_MAX_MG_LEVEL, cint]
    ##  Number of null-space vectors to use on each level
    n_vec* {.importc: "n_vec".}: array[QUDA_MAX_MG_LEVEL, cint]
    ##  Precision to store the null-space vectors in (post block orthogonalization)
    precision_null* {.importc: "precision_null".}: array[QUDA_MAX_MG_LEVEL,
        QudaPrecision]
    ##  Number of times to repeat Gram-Schmidt in block orthogonalization
    n_block_ortho* {.importc: "n_block_ortho".}: array[QUDA_MAX_MG_LEVEL, cint]
    ##  Whether to do passes at block orthogonalize in fixed point for improved accuracy
    block_ortho_two_pass* {.importc: "block_ortho_two_pass".}: array[
        QUDA_MAX_MG_LEVEL, QudaBoolean]
    ##  Verbosity on each level of the multigrid
    verbosity* {.importc: "verbosity".}: array[QUDA_MAX_MG_LEVEL, QudaVerbosity]
    ##  Setup MMA usage on each level of the multigrid
    setup_use_mma* {.importc: "setup_use_mma".}: array[QUDA_MAX_MG_LEVEL,
        QudaBoolean]
    ##  Dslash MMA usage on each level of the multigrid
    dslash_use_mma* {.importc: "dslash_use_mma".}: array[QUDA_MAX_MG_LEVEL,
        QudaBoolean]
    ##  Transfer MMA usage on each level of the multigrid
    transfer_use_mma* {.importc: "transfer_use_mma".}: array[QUDA_MAX_MG_LEVEL,
        QudaBoolean]
    ##  Inverter to use in the setup phase
    setup_inv_type* {.importc: "setup_inv_type".}: array[QUDA_MAX_MG_LEVEL,
        QudaInverterType]
    ##  Solver batch size to use in the setup phase
    n_vec_batch* {.importc: "n_vec_batch".}: array[QUDA_MAX_MG_LEVEL, cint]
    ##  Number of setup iterations
    num_setup_iter* {.importc: "num_setup_iter".}: array[QUDA_MAX_MG_LEVEL, cint]
    ##  Tolerance to use in the setup phase
    setup_tol* {.importc: "setup_tol".}: array[QUDA_MAX_MG_LEVEL, cdouble]
    ##  Maximum number of iterations for each setup solver
    setup_maxiter* {.importc: "setup_maxiter".}: array[QUDA_MAX_MG_LEVEL, cint]
    ##  Maximum number of iterations for refreshing the null-space vectors
    setup_maxiter_refresh* {.importc: "setup_maxiter_refresh".}: array[
        QUDA_MAX_MG_LEVEL, cint]
    ##  Basis to use for CA solver setup
    setup_ca_basis* {.importc: "setup_ca_basis".}: array[QUDA_MAX_MG_LEVEL,
        QudaCABasis]
    ##  Basis size for CA solver setup
    setup_ca_basis_size* {.importc: "setup_ca_basis_size".}: array[
        QUDA_MAX_MG_LEVEL, cint]
    ##  Minimum eigenvalue for Chebyshev CA basis
    setup_ca_lambda_min* {.importc: "setup_ca_lambda_min".}: array[
        QUDA_MAX_MG_LEVEL, cdouble]
    ##  Maximum eigenvalue for Chebyshev CA basis
    setup_ca_lambda_max* {.importc: "setup_ca_lambda_max".}: array[
        QUDA_MAX_MG_LEVEL, cdouble]
    ##  Null-space type to use in the setup phase
    setup_type* {.importc: "setup_type".}: QudaSetupType
    ##  Pre orthonormalize vectors in the setup phase
    pre_orthonormalize* {.importc: "pre_orthonormalize".}: QudaBoolean
    ##  Post orthonormalize vectors in the setup phase
    post_orthonormalize* {.importc: "post_orthonormalize".}: QudaBoolean
    ##  The solver that wraps around the coarse grid correction and smoother
    coarse_solver* {.importc: "coarse_solver".}: array[QUDA_MAX_MG_LEVEL,
        QudaInverterType]
    ##  Tolerance for the solver that wraps around the coarse grid correction and smoother
    coarse_solver_tol* {.importc: "coarse_solver_tol".}: array[QUDA_MAX_MG_LEVEL,
        cdouble]
    ##  Maximum number of iterations for the solver that wraps around the coarse grid correction and smoother
    coarse_solver_maxiter* {.importc: "coarse_solver_maxiter".}: array[
        QUDA_MAX_MG_LEVEL, cint]
    ##  Basis to use for CA coarse solvers
    coarse_solver_ca_basis* {.importc: "coarse_solver_ca_basis".}: array[
        QUDA_MAX_MG_LEVEL, QudaCABasis]
    ##  Basis size for CA coarse solvers
    coarse_solver_ca_basis_size* {.importc: "coarse_solver_ca_basis_size".}: array[
        QUDA_MAX_MG_LEVEL, cint]
    ##  Minimum eigenvalue for Chebyshev CA basis
    coarse_solver_ca_lambda_min* {.importc: "coarse_solver_ca_lambda_min".}: array[
        QUDA_MAX_MG_LEVEL, cdouble]
    ##  Maximum eigenvalue for Chebyshev CA basis
    coarse_solver_ca_lambda_max* {.importc: "coarse_solver_ca_lambda_max".}: array[
        QUDA_MAX_MG_LEVEL, cdouble]
    ##  Smoother to use on each level
    smoother* {.importc: "smoother".}: array[QUDA_MAX_MG_LEVEL, QudaInverterType]
    ##  Tolerance to use for the smoother / solver on each level
    smoother_tol* {.importc: "smoother_tol".}: array[QUDA_MAX_MG_LEVEL, cdouble]
    ##  Number of pre-smoother applications on each level
    nu_pre* {.importc: "nu_pre".}: array[QUDA_MAX_MG_LEVEL, cint]
    ##  Number of post-smoother applications on each level
    nu_post* {.importc: "nu_post".}: array[QUDA_MAX_MG_LEVEL, cint]
    ##  Basis to use for CA smoother solvers
    smoother_solver_ca_basis* {.importc: "smoother_solver_ca_basis".}: array[
        QUDA_MAX_MG_LEVEL, QudaCABasis]
    ##  Minimum eigenvalue for Chebyshev CA smoother basis
    smoother_solver_ca_lambda_min* {.importc: "smoother_solver_ca_lambda_min".}: array[
        QUDA_MAX_MG_LEVEL, cdouble]
    ##  Maximum eigenvalue for Chebyshev CA smoother basis
    smoother_solver_ca_lambda_max* {.importc: "smoother_solver_ca_lambda_max".}: array[
        QUDA_MAX_MG_LEVEL, cdouble]
    ##  Over/under relaxation factor for the smoother at each level
    omega* {.importc: "omega".}: array[QUDA_MAX_MG_LEVEL, cdouble]
    ##  Precision to use for halo communication in the smoother
    smoother_halo_precision* {.importc: "smoother_halo_precision".}: array[
        QUDA_MAX_MG_LEVEL, QudaPrecision]
    ##  Whether to use additive or multiplicative Schwarz preconditioning in the smoother
    smoother_schwarz_type* {.importc: "smoother_schwarz_type".}: array[
        QUDA_MAX_MG_LEVEL, QudaSchwarzType]
    ##  Number of Schwarz cycles to apply
    smoother_schwarz_cycle* {.importc: "smoother_schwarz_cycle".}: array[
        QUDA_MAX_MG_LEVEL, cint]
    ##  The type of residual to send to the next coarse grid, and thus the
    ## 	type of solution to receive back from this coarse grid
    coarse_grid_solution_type* {.importc: "coarse_grid_solution_type".}: array[
        QUDA_MAX_MG_LEVEL, QudaSolutionType]
    ##  The type of smoother solve to do on each grid (e/o preconditioning or not)
    smoother_solve_type* {.importc: "smoother_solve_type".}: array[
        QUDA_MAX_MG_LEVEL, QudaSolveType]
    ##  The type of multigrid cycle to perform at each level
    cycle_type* {.importc: "cycle_type".}: array[QUDA_MAX_MG_LEVEL,
        QudaMultigridCycleType]
    ##  Whether to use global reductions or not for the smoother / solver at each level
    global_reduction* {.importc: "global_reduction".}: array[QUDA_MAX_MG_LEVEL,
        QudaBoolean]
    ##  Location where each level should be done
    location* {.importc: "location".}: array[QUDA_MAX_MG_LEVEL, QudaFieldLocation]
    ##  Location where the coarse-operator construction will be computedn
    setup_location* {.importc: "setup_location".}: array[QUDA_MAX_MG_LEVEL,
        QudaFieldLocation]
    ##  Whether to use eigenvectors for the nullspace or, if the coarsest instance deflate
    use_eig_solver* {.importc: "use_eig_solver".}: array[QUDA_MAX_MG_LEVEL,
        QudaBoolean]
    ##  Whether to compute the null vectors or reload them
    compute_null_vector* {.importc: "compute_null_vector".}: QudaComputeNullVector
    ##  Whether to generate on all levels or just on level 0
    generate_all_levels* {.importc: "generate_all_levels".}: QudaBoolean
    ##  Whether to run the verification checks once set up is complete
    run_verify* {.importc: "run_verify".}: QudaBoolean
    ##  Whether to run null Vs eigen vector overlap checks once set up is complete
    run_low_mode_check* {.importc: "run_low_mode_check".}: QudaBoolean
    ##  Whether to run null vector oblique checks once set up is complete
    run_oblique_proj_check* {.importc: "run_oblique_proj_check".}: QudaBoolean
    ##  Whether to load the null-space vectors to disk (requires QIO)
    vec_load* {.importc: "vec_load".}: array[QUDA_MAX_MG_LEVEL, QudaBoolean]
    ##  Filename prefix where to load the null-space vectors
    vec_infile* {.importc: "vec_infile".}: array[QUDA_MAX_MG_LEVEL, array[256, char]]
    ##  Whether to store the null-space vectors to disk (requires QIO)
    vec_store* {.importc: "vec_store".}: array[QUDA_MAX_MG_LEVEL, QudaBoolean]
    ##  Filename prefix for where to save the null-space vectors
    vec_outfile* {.importc: "vec_outfile".}: array[QUDA_MAX_MG_LEVEL,
        array[256, char]]
    ##  Whether to store the null-space vectors in singlefile or partfile format
    mg_vec_partfile* {.importc: "mg_vec_partfile".}: array[QUDA_MAX_MG_LEVEL,
        QudaBoolean]
    ##  Whether to use and initial guess during coarse grid deflation
    coarse_guess* {.importc: "coarse_guess".}: QudaBoolean
    ##  Whether to preserve the deflation space during MG update
    preserve_deflation* {.importc: "preserve_deflation".}: QudaBoolean
    ##  Multiplicative factor for the mu parameter
    mu_factor* {.importc: "mu_factor".}: array[QUDA_MAX_MG_LEVEL, cdouble]
    ##  Boolean for aggregation type, implies staggered or not
    transfer_type* {.importc: "transfer_type".}: array[QUDA_MAX_MG_LEVEL,
        QudaTransferType]
    ##  Whether or not to let MG coarsening drop improvements, for ex dropping long links in small aggregation dimensions
    allow_truncation* {.importc: "allow_truncation".}: QudaBoolean
    ##  Whether or not to use the dagger approximation for the KD preconditioned operator
    staggered_kd_dagger_approximation* {.importc: "staggered_kd_dagger_approximation".}: QudaBoolean
    ##  Whether to do a full (false) or thin (true) update in the context of updateMultigridQuda
    thin_update_only* {.importc: "thin_update_only".}: QudaBoolean

  QudaGaugeObservableParam* {.importc: "QudaGaugeObservableParam",
                             header: "quda.h", bycopy.} = object
    struct_size* {.importc: "struct_size".}: csize_t
    ## < Size of this struct in bytes.  Used to ensure that the host application and QUDA see the same struct
    su_project* {.importc: "su_project".}: QudaBoolean
    ## < Whether to project onto the manifold prior to measurement
    compute_plaquette* {.importc: "compute_plaquette".}: QudaBoolean
    ## < Whether to compute the plaquette
    plaquette* {.importc: "plaquette".}: array[3, cdouble]
    ## < Total, spatial and temporal field energies, respectively
    compute_rectangle* {.importc: "compute_rectangle".}: QudaBoolean
    ## < Whether to compute the rectangle
    rectangle* {.importc: "rectangle".}: array[3, cdouble]
    ## < Total, spatial and temporal rectangle, respectively
    compute_polyakov_loop* {.importc: "compute_polyakov_loop".}: QudaBoolean
    ## < Whether to compute the temporal Polyakov loop
    ploop* {.importc: "ploop".}: array[2, cdouble]
    ## < Real and imaginary part of temporal Polyakov loop
    compute_gauge_loop_trace* {.importc: "compute_gauge_loop_trace".}: QudaBoolean
    ## < Whether to compute gauge loop traces
    traces* {.importc: "traces".}: ptr double_complex
    ## < Individual complex traces of each loop
    input_path_buff* {.importc: "input_path_buff".}: ptr ptr cint
    ## < Array of paths
    path_length* {.importc: "path_length".}: ptr cint
    ## < Length of each path
    loop_coeff* {.importc: "loop_coeff".}: ptr cdouble
    ## < Multiplicative factor for each loop
    num_paths* {.importc: "num_paths".}: cint
    ## < Total number of paths
    max_length* {.importc: "max_length".}: cint
    ## < Maximum length of any path
    factor* {.importc: "factor".}: cdouble
    ## < Global multiplicative factor to apply to each loop trace
    compute_qcharge* {.importc: "compute_qcharge".}: QudaBoolean
    ## < Whether to compute the topological charge and field energy
    qcharge* {.importc: "qcharge".}: cdouble
    ## < Computed topological charge
    energy* {.importc: "energy".}: array[3, cdouble]
    ## < Total, spatial and temporal field energies, respectively
    compute_qcharge_density* {.importc: "compute_qcharge_density".}: QudaBoolean
    ## < Whether to compute the topological charge density
    qcharge_density* {.importc: "qcharge_density".}: pointer
    ## < Pointer to host array of length volume where the q-charge density will be copied
    remove_staggered_phase* {.importc: "remove_staggered_phase".}: QudaBoolean
    ## < Whether or not the resident gauge field has staggered phases applied and if they should
    ##                                  be removed; this was needed for the Polyakov loop calculation when called through MILC,
    ##                                  with the underlying issue documented https://github.com/lattice/quda/issues/1315

  QudaGaugeSmearParam* {.importc: "QudaGaugeSmearParam", header: "quda.h", bycopy.} = object
    struct_size* {.importc: "struct_size".}: csize_t
    ## < Size of this struct in bytes.  Used to ensure that the host application and QUDA see the same struct
    n_steps* {.importc: "n_steps".}: cuint
    ## < The total number of smearing steps to perform.
    epsilon* {.importc: "epsilon".}: cdouble
    ## < Serves as one of the coefficients in Over Improved Stout smearing, or as the step size in
    ##                              Wilson/Symanzik flow
    smear_anisotropy* {.importc: "smear_anisotropy".}: cdouble
    ##  Used in anisotropic Wilson/Symanzik flow and APE, STOUT, and OvrimpSTOUT *
    rk_order* {.importc: "rk_order".}: cuint
    ##  Order of the Runga-Kutta integrator: 3 or 4 *
    alpha* {.importc: "alpha".}: cdouble
    ## < The single coefficient used in APE smearing
    rho* {.importc: "rho".}: cdouble
    ## < Serves as one of the coefficients used in Over Improved Stout smearing, or as the single coefficient used in Stout
    alpha1* {.importc: "alpha1".}: cdouble
    ## < The coefficient used in HYP smearing step 3 (will not be used in 3D smearing)
    alpha2* {.importc: "alpha2".}: cdouble
    ## < The coefficient used in HYP smearing step 2
    alpha3* {.importc: "alpha3".}: cdouble
    ## < The coefficient used in HYP smearing step 1
    meas_interval* {.importc: "meas_interval".}: cuint
    ## < Perform the requested measurements on the gauge field at this interval
    smear_type* {.importc: "smear_type".}: QudaGaugeSmearType
    ## < The smearing type to perform
    adj_n_save* {.importc: "adj_n_save".}: cuint
    ## < How many intermediate gauge fields to save at each large nblock to perform adj flow
    hier_threshold* {.importc: "hier_threshold".}: cuint
    ## < Minimum *hierarchical* threshold for adj gradient flow
    restart* {.importc: "restart".}: QudaBoolean
    ## < Used to restart the smearing from existing gaugeSmeared
    t0* {.importc: "t0".}: cdouble
    ## < Starting flow time for Wilson flow
    dir_ignore* {.importc: "dir_ignore".}: cint
    ## < The direction to be ignored by the smearing algorithm
    ##                                         A negative value means 3D for APE/STOUT and 4D for OVRIMP_STOUT/HYP

  QudaBLASParam* {.importc: "QudaBLASParam", header: "quda.h", bycopy.} = object
    struct_size* {.importc: "struct_size".}: csize_t
    ## < Size of this struct in bytes.  Used to ensure that the host application and QUDA see the same struct
    blas_type* {.importc: "blas_type".}: QudaBLASType
    ## < Type of BLAS computation to perfrom
    ##  GEMM params
    trans_a* {.importc: "trans_a".}: QudaBLASOperation
    ## < operation op(A) that is non- or (conj.) transpose.
    trans_b* {.importc: "trans_b".}: QudaBLASOperation
    ## < operation op(B) that is non- or (conj.) transpose.
    m* {.importc: "m".}: cint
    ## < number of rows of matrix op(A) and C.
    n* {.importc: "n".}: cint
    ## < number of columns of matrix op(B) and C.
    k* {.importc: "k".}: cint
    ## < number of columns of op(A) and rows of op(B).
    lda* {.importc: "lda".}: cint
    ## < leading dimension of two-dimensional array used to store the matrix A.
    ldb* {.importc: "ldb".}: cint
    ## < leading dimension of two-dimensional array used to store matrix B.
    ldc* {.importc: "ldc".}: cint
    ## < leading dimension of two-dimensional array used to store matrix C.
    a_offset* {.importc: "a_offset".}: cint
    ## < position of the A array from which begin read/write.
    b_offset* {.importc: "b_offset".}: cint
    ## < position of the B array from which begin read/write.
    c_offset* {.importc: "c_offset".}: cint
    ## < position of the C array from which begin read/write.
    a_stride* {.importc: "a_stride".}: cint
    ## < stride of the A array in strided(batched) mode
    b_stride* {.importc: "b_stride".}: cint
    ## < stride of the B array in strided(batched) mode
    c_stride* {.importc: "c_stride".}: cint
    ## < stride of the C array in strided(batched) mode
    alpha* {.importc: "alpha".}: double_complex
    ## < scalar used for multiplication.
    beta* {.importc: "beta".}: double_complex
    ## < scalar used for multiplication. If beta==0, C does not have to be a valid input.
    ##  LU inversion params
    inv_mat_size* {.importc: "inv_mat_size".}: cint
    ## < The rank of the square matrix in the LU inversion
    ##  Common params
    batch_count* {.importc: "batch_count".}: cint
    ## < number of pointers contained in arrayA, arrayB and arrayC.
    data_type* {.importc: "data_type".}: QudaBLASDataType
    ## < Specifies if using S(C) or D(Z) BLAS type
    data_order* {.importc: "data_order".}: QudaBLASDataOrder
    ## < Specifies if using Row or Column major


##
##  Interface functions, found in interface_quda.cpp
##
##
##  Set parameters related to status reporting.
##
##  In typical usage, this function will be called once (or not at
##  all) just before the call to initQuda(), but it's valid to call
##  it any number of times at any point during execution.  Prior to
##  the first time it's called, the parameters take default values
##  as indicated below.
##
##  @param verbosity  Default verbosity, ranging from QUDA_SILENT to
##                    QUDA_DEBUG_VERBOSE.  Within a solver, this
##                    parameter is overridden by the "verbosity"
##                    member of QudaInvertParam.  The default value
##                    is QUDA_SUMMARIZE.
##
##  @param prefix     String to prepend to all messages from QUDA.  This
##                    defaults to the empty string (""), but you may
##                    wish to specify something like "QUDA: " to
##                    distinguish QUDA's output from that of your
##                    application.
##
##  @param outfile    File pointer (such as stdout, stderr, or a handle
##                    returned by fopen()) where messages should be
##                    printed.  The default is stdout.
##

proc setVerbosityQuda*(verbosity: QudaVerbosity; prefix: ptr char; outfile: ptr FILE) {.
    importc: "setVerbosityQuda", header: "quda.h".}
##
##  initCommsGridQuda() takes an optional "rank_from_coords" argument that
##  should be a pointer to a user-defined function with this prototype.
##
##  @param coords  Node coordinates
##  @param fdata   Any auxiliary data needed by the function
##  @return        MPI rank or QMP node ID cooresponding to the node coordinates
##
##  @see initCommsGridQuda
##

type
  QudaCommsMap* = proc (coords: ptr ConstInt; fdata: pointer): cint {.cdecl.}

##
##  @param mycomm User provided MPI communicator in place of MPI_COMM_WORLD
##

proc qudaSetCommHandle*(mycomm: pointer) {.importc: "qudaSetCommHandle",
                                        header: "quda.h".}
##
##  Declare the grid mapping ("logical topology" in QMP parlance)
##  used for communications in a multi-GPU grid.  This function
##  should be called prior to initQuda().  The only case in which
##  it's optional is when QMP is used for communication and the
##  logical topology has already been declared by the application.
##
##  @param nDim   Number of grid dimensions.  "4" is the only supported
##                value currently.
##
##  @param dims   Array of grid dimensions.  dims[0]*dims[1]*dims[2]*dims[3]
##                must equal the total number of MPI ranks or QMP nodes.
##
##  @param func   Pointer to a user-supplied function that maps coordinates
##                in the communication grid to MPI ranks (or QMP node IDs).
##                If the pointer is NULL, the default mapping depends on
##                whether QMP or MPI is being used for communication.  With
##                QMP, the existing logical topology is used if it's been
##                declared.  With MPI or as a fallback with QMP, the default
##                ordering is lexicographical with the fourth ("t") index
##                varying fastest.
##
##  @param fdata  Pointer to any data required by "func" (may be NULL)
##
##  @see QudaCommsMap
##

proc initCommsGridQuda*(nDim: cint; dims: ptr cint; `func`: QudaCommsMap; fdata: pointer) {.
    importc: "initCommsGridQuda", header: "quda.h".}
##
##  Initialize the library.  This is a low-level interface that is
##  called by initQuda.  Calling initQudaDevice requires that the
##  user also call initQudaMemory before using QUDA.
##
##  @param device CUDA device number to use.  In a multi-GPU build,
##                this parameter may either be set explicitly on a
##                per-process basis or set to -1 to enable a default
##                allocation of devices to processes.
##

proc initQudaDevice*(device: cint) {.importc: "initQudaDevice", header: "quda.h".}
##
##  Initialize the library persistant memory allocations (both host
##  and device).  This is a low-level interface that is called by
##  initQuda.  Calling initQudaMemory requires that the user has
##  previously called initQudaDevice.
##

proc initQudaMemory*() {.importc: "initQudaMemory", header: "quda.h".}
##
##  Initialize the library.  This function is actually a wrapper
##  around calls to initQudaDevice() and initQudaMemory().
##
##  @param device  CUDA device number to use.  In a multi-GPU build,
##                 this parameter may either be set explicitly on a
##                 per-process basis or set to -1 to enable a default
##                 allocation of devices to processes.
##

proc initQuda*(device: cint) {.importc: "initQuda", header: "quda.h".}
##
##  Finalize the library.
##

proc endQuda*() {.importc: "endQuda", header: "quda.h".}
##
##  @brief update the radius for halos.
##  @details This should only be needed for automated testing when
##  different partitioning is applied within a single run.
##

proc updateR*() {.importc: "updateR", header: "quda.h".}
##
##  A new QudaGaugeParam should always be initialized immediately
##  after it's defined (and prior to explicitly setting its members)
##  using this function.  Typical usage is as follows:
##
##    QudaGaugeParam gauge_param = newQudaGaugeParam();
##

proc newQudaGaugeParam*(): QudaGaugeParam {.importc: "newQudaGaugeParam",
    header: "quda.h".}
##
##  A new QudaInvertParam should always be initialized immediately
##  after it's defined (and prior to explicitly setting its members)
##  using this function.  Typical usage is as follows:
##
##    QudaInvertParam invert_param = newQudaInvertParam();
##

proc newQudaInvertParam*(): QudaInvertParam {.importc: "newQudaInvertParam",
    header: "quda.h".}
##
##  A new QudaMultigridParam should always be initialized immediately
##  after it's defined (and prior to explicitly setting its members)
##  using this function.  Typical usage is as follows:
##
##    QudaMultigridParam mg_param = newQudaMultigridParam();
##

proc newQudaMultigridParam*(): QudaMultigridParam {.
    importc: "newQudaMultigridParam", header: "quda.h".}
##
##  A new QudaEigParam should always be initialized immediately
##  after it's defined (and prior to explicitly setting its members)
##  using this function.  Typical usage is as follows:
##
##    QudaEigParam eig_param = newQudaEigParam();
##

proc newQudaEigParam*(): QudaEigParam {.importc: "newQudaEigParam",
                                     header: "quda.h".}
##
##  A new QudaGaugeObservableParam should always be initialized
##  immediately after it's defined (and prior to explicitly setting
##  its members) using this function.  Typical usage is as follows:
##
##    QudaGaugeObservalbeParam obs_param = newQudaGaugeObservableParam();
##

proc newQudaGaugeObservableParam*(): QudaGaugeObservableParam {.
    importc: "newQudaGaugeObservableParam", header: "quda.h".}
##
##  A new QudaGaugeSmearParam should always be initialized
##  immediately after it's defined (and prior to explicitly setting
##  its members) using this function.  Typical usage is as follows:
##
##    QudaGaugeSmearParam smear_param = newQudaGaugeSmearParam();
##

proc newQudaGaugeSmearParam*(): QudaGaugeSmearParam {.
    importc: "newQudaGaugeSmearParam", header: "quda.h".}
##
##  A new QudaBLASParam should always be initialized immediately
##  after it's defined (and prior to explicitly setting its members)
##  using this function.  Typical usage is as follows:
##
##    QudaBLASParam blas_param = newQudaBLASParam();
##

proc newQudaBLASParam*(): QudaBLASParam {.importc: "newQudaBLASParam",
                                       header: "quda.h".}
##
##  Print the members of QudaGaugeParam.
##  @param param The QudaGaugeParam whose elements we are to print.
##

proc printQudaGaugeParam*(param: ptr QudaGaugeParam) {.
    importc: "printQudaGaugeParam", header: "quda.h".}
##
##  Print the members of QudaInvertParam.
##  @param param The QudaInvertParam whose elements we are to print.
##

proc printQudaInvertParam*(param: ptr QudaInvertParam) {.
    importc: "printQudaInvertParam", header: "quda.h".}
##
##  Print the members of QudaMultigridParam.
##  @param param The QudaMultigridParam whose elements we are to print.
##

proc printQudaMultigridParam*(param: ptr QudaMultigridParam) {.
    importc: "printQudaMultigridParam", header: "quda.h".}
##
##  Print the members of QudaEigParam.
##  @param param The QudaEigParam whose elements we are to print.
##

proc printQudaEigParam*(param: ptr QudaEigParam) {.importc: "printQudaEigParam",
    header: "quda.h".}
##
##  Print the members of QudaGaugeObservableParam.
##  @param param The QudaGaugeObservableParam whose elements we are to print.
##

proc printQudaGaugeObservableParam*(param: ptr QudaGaugeObservableParam) {.
    importc: "printQudaGaugeObservableParam", header: "quda.h".}
##
##  Print the members of QudaBLASParam.
##  @param param The QudaBLASParam whose elements we are to print.
##

proc printQudaBLASParam*(param: ptr QudaBLASParam) {.importc: "printQudaBLASParam",
    header: "quda.h".}
##
##  Load the gauge field from the host.
##  @param h_gauge Base pointer to host gauge field (regardless of dimensionality)
##  @param param   Contains all metadata regarding host and device storage
##

proc loadGaugeQuda*(h_gauge: pointer; param: ptr QudaGaugeParam) {.
    importc: "loadGaugeQuda", header: "quda.h".}
##
##  Free QUDA's internal copy of the gauge field.
##

proc freeGaugeQuda*() {.importc: "freeGaugeQuda", header: "quda.h".}
##
##  Free a unique type (Wilson, HISQ fat, HISQ long, smeared) of internal gauge field.
##  @param link_type[in] Type of link type to free up
##

proc freeUniqueGaugeQuda*(link_type: QudaLinkType) {.
    importc: "freeUniqueGaugeQuda", header: "quda.h".}
##
##  Free QUDA's internal smeared gauge field.
##

proc freeGaugeSmearedQuda*() {.importc: "freeGaugeSmearedQuda", header: "quda.h".}
##
##  Free QUDA's internal two-link gauge field.
##

proc freeGaugeTwoLinkQuda*() {.importc: "freeGaugeTwoLinkQuda", header: "quda.h".}
##
##  Save the gauge field to the host.
##  @param h_gauge Base pointer to host gauge field (regardless of dimensionality)
##  @param param   Contains all metadata regarding host and device storage
##

proc saveGaugeQuda*(h_gauge: pointer; param: ptr QudaGaugeParam) {.
    importc: "saveGaugeQuda", header: "quda.h".}
##
##  Load the clover term and/or the clover inverse from the host.
##  Either h_clover or h_clovinv may be set to NULL.
##  @param h_clover    Base pointer to host clover field
##  @param h_cloverinv Base pointer to host clover inverse field
##  @param inv_param   Contains all metadata regarding host and device storage
##

proc loadCloverQuda*(h_clover: pointer; h_clovinv: pointer;
                    inv_param: ptr QudaInvertParam) {.importc: "loadCloverQuda",
    header: "quda.h".}
##
##  Free QUDA's internal copy of the clover term and/or clover inverse.
##

proc freeCloverQuda*() {.importc: "freeCloverQuda", header: "quda.h".}
##
##  Perform the solve, according to the parameters set in param.  It
##  is assumed that the gauge field has already been loaded via
##  loadGaugeQuda().
##  @param h_x    Solution spinor field
##  @param h_b    Source spinor field
##  @param param  Contains all metadata regarding host and device
##                storage and solver parameters
##

proc lanczosQuda*(k0: cint; m: cint; hp_Apsi: pointer; hp_r: pointer; hp_V: pointer;
                 hp_alpha: pointer; hp_beta: pointer; eig_param: ptr QudaEigParam) {.
    importc: "lanczosQuda", header: "quda.h".}
##
##  Perform the eigensolve. The problem matrix is defined by the invert param, the
##  mode of solution is specified by the eig param. It is assumed that the gauge
##  field has already been loaded via  loadGaugeQuda().
##  @param h_evecs  Array of pointers to application eigenvectors
##  @param h_evals  Host side eigenvalues
##  @param param Contains all metadata regarding the type of solve.
##

proc eigensolveQuda*(h_evecs: ptr pointer; h_evals: ptr double_complex;
                    param: ptr QudaEigParam) {.importc: "eigensolveQuda",
    header: "quda.h".}
##
##  Perform the solve, according to the parameters set in param.  It
##  is assumed that the gauge field has already been loaded via
##  loadGaugeQuda().
##  @param h_x    Solution spinor field
##  @param h_b    Source spinor field
##  @param param  Contains all metadata regarding host and device
##                storage and solver parameters
##

proc invertQuda*(h_x: pointer; h_b: pointer; param: ptr QudaInvertParam) {.
    importc: "invertQuda", header: "quda.h".}
##
##  @brief Perform the solve like @invertQuda but for multiple rhs by spliting the comm grid into
##  sub-partitions: each sub-partition invert one or more rhs'.
##  The QudaInvertParam object specifies how the solve should be performed on each sub-partition.
##  Unlike @invertQuda, the interface also takes the host side gauge as input. The gauge pointer and
##  gauge_param are used if for inv_param split_grid[0] * split_grid[1] * split_grid[2] * split_grid[3]
##  is larger than 1, in which case gauge field is not required to be loaded beforehand; otherwise
##  this interface would just work as @invertQuda, which requires gauge field to be loaded beforehand,
##  and the gauge field pointer and gauge_param are not used.
##  @param _hp_x       Array of solution spinor fields
##  @param _hp_b       Array of source spinor fields
##  @param param       Contains all metadata regarding host and device storage and solver parameters
##

proc invertMultiSrcQuda*(hp_x: ptr pointer; hp_b: ptr pointer;
                        param: ptr QudaInvertParam) {.
    importc: "invertMultiSrcQuda", header: "quda.h".}
##
##  Solve for multiple shifts (e.g., masses).
##  @param _hp_x    Array of solution spinor fields
##  @param _hp_b    Source spinor fields
##  @param param  Contains all metadata regarding host and device
##                storage and solver parameters
##

proc invertMultiShiftQuda*(hp_x: ptr pointer; hp_b: pointer;
                          param: ptr QudaInvertParam) {.
    importc: "invertMultiShiftQuda", header: "quda.h".}
##
##  Setup the multigrid solver, according to the parameters set in param.  It
##  is assumed that the gauge field has already been loaded via
##  loadGaugeQuda().
##  @param param  Contains all metadata regarding host and device
##                storage and solver parameters
##

proc newMultigridQuda*(param: ptr QudaMultigridParam): pointer {.
    importc: "newMultigridQuda", header: "quda.h".}
##
##  @brief Free resources allocated by the multigrid solver
##  @param mg_instance Pointer to instance of multigrid_solver
##  @param param Contains all metadata regarding host and device
##  storage and solver parameters
##

proc destroyMultigridQuda*(mg_instance: pointer) {.importc: "destroyMultigridQuda",
    header: "quda.h".}
##
##  @brief Updates the multigrid preconditioner for the new gauge / clover field
##  @param mg_instance Pointer to instance of multigrid_solver
##  @param param Contains all metadata regarding host and device
##  storage and solver parameters, of note contains a flag specifying whether
##  to do a full update or a thin update.
##

proc updateMultigridQuda*(mg_instance: pointer; param: ptr QudaMultigridParam) {.
    importc: "updateMultigridQuda", header: "quda.h".}
##
##  @brief Dump the null-space vectors to disk
##  @param[in] mg_instance Pointer to the instance of multigrid_solver
##  @param[in] param Contains all metadata regarding host and device
##  storage and solver parameters (QudaMultigridParam::vec_outfile
##  sets the output filename prefix).
##

proc dumpMultigridQuda*(mg_instance: pointer; param: ptr QudaMultigridParam) {.
    importc: "dumpMultigridQuda", header: "quda.h".}
##
##  Apply the Dslash operator (D_{eo} or D_{oe}).
##  @param[out] h_out  Result spinor field
##  @param[in] h_in   Input spinor field
##  @param[in] param  Contains all metadata regarding host and device
##                storage
##  @param[in] parity The destination parity of the field
##

proc dslashQuda*(h_out: pointer; h_in: pointer; inv_param: ptr QudaInvertParam;
                parity: QudaParity) {.importc: "dslashQuda", header: "quda.h".}
##
##  Apply the covariant derivative.
##  @param[out] h_out  Result spinor field
##  @param[in] h_in   Input spinor field
##  @param[in] dir    Direction of application
##  @param[in] param  Metadata for host and device storage
##

proc covDevQuda*(h_out: pointer; h_in: pointer; dir: cint; param: ptr QudaInvertParam) {.
    importc: "covDevQuda", header: "quda.h".}
##
##  Apply the covariant derivative.
##  @param[out] h_out  Result spinor field
##  @param[in] h_in   Input spinor field
##  @param[in] dir    Direction of application
##  @param[in] sym    Apply forward=2, backward=2 or symmetric=3 shift
##  @param[in] param  Metadata for host and device storage
##

proc shiftQuda*(h_out: pointer; h_in: pointer; dir: cint; sym: cint;
               param: ptr QudaInvertParam) {.importc: "shiftQuda",
    header: "quda.h".}
##
##  Apply the spin-taste operator.
##  @param[out] h_out  Result spinor field
##  @param[in] h_in   Input spinor field
##  @param[in] spin   Spin gamma structure
##  @param[in] taste  Taste gamma structure
##  @param[in] param  Metadata for host and device storage
##

proc spinTasteQuda*(h_out: pointer; h_in: pointer; spin: cint; taste: cint;
                   param: ptr QudaInvertParam) {.importc: "spinTasteQuda",
    header: "quda.h".}
##
##  @brief Perform the solve like @dslashQuda but for multiple rhs by spliting the comm grid into
##  sub-partitions: each sub-partition does one or more rhs'.
##  The QudaInvertParam object specifies how the solve should be performed on each sub-partition.
##  Unlike @invertQuda, the interface also takes the host side gauge as
##  input - gauge field is not required to be loaded beforehand.
##  @param _hp_x       Array of solution spinor fields
##  @param _hp_b       Array of source spinor fields
##  @param param       Contains all metadata regarding host and device storage and solver parameters
##  @param parity      Parity to apply dslash on
##

proc dslashMultiSrcQuda*(hp_x: ptr pointer; hp_b: ptr pointer;
                        param: ptr QudaInvertParam; parity: QudaParity) {.
    importc: "dslashMultiSrcQuda", header: "quda.h".}
##
##  Apply the clover operator or its inverse.
##  @param h_out  Result spinor field
##  @param h_in   Input spinor field
##  @param param  Contains all metadata regarding host and device
##                storage
##  @param parity The source and destination parity of the field
##  @param inverse Whether to apply the inverse of the clover term
##

proc cloverQuda*(h_out: pointer; h_in: pointer; inv_param: ptr QudaInvertParam;
                parity: QudaParity; inverse: cint) {.importc: "cloverQuda",
    header: "quda.h".}
##
##  Apply the full Dslash matrix, possibly even/odd preconditioned.
##  @param h_out  Result spinor field
##  @param h_in   Input spinor field
##  @param param  Contains all metadata regarding host and device
##                storage
##

proc MatQuda*(h_out: pointer; h_in: pointer; inv_param: ptr QudaInvertParam) {.
    importc: "MatQuda", header: "quda.h".}
##
##  Apply M^{\dag}M, possibly even/odd preconditioned.
##  @param h_out  Result spinor field
##  @param h_in   Input spinor field
##  @param param  Contains all metadata regarding host and device
##                storage
##

proc MatDagMatQuda*(h_out: pointer; h_in: pointer; inv_param: ptr QudaInvertParam) {.
    importc: "MatDagMatQuda", header: "quda.h".}
##
##  The following routines are temporary additions used by the HISQ
##  link-fattening code.
##

proc set_dim*(a1: ptr cint) {.importc: "set_dim", header: "quda.h".}
proc pack_ghost*(cpuLink: ptr pointer; cpuGhost: ptr pointer; nFace: cint;
                precision: QudaPrecision) {.importc: "pack_ghost",
    header: "quda.h".}
proc computeKSLinkQuda*(fatlink: pointer; longlink: pointer; ulink: pointer;
                       inlink: pointer; path_coeff: ptr cdouble;
                       param: ptr QudaGaugeParam) {.importc: "computeKSLinkQuda",
    header: "quda.h".}
##
##  Compute two-link field
##
##  @param[out] twolink computed two-link field
##  @param[in] inlink  the external field
##  @param[in] param  Contains all metadata regarding host and device
##                storage
##

proc computeTwoLinkQuda*(twolink: pointer; inlink: pointer; param: ptr QudaGaugeParam) {.
    importc: "computeTwoLinkQuda", header: "quda.h".}
##
##  Either downloads and sets the resident momentum field, or uploads
##  and returns the resident momentum field
##
##  @param[in,out] mom The external momentum field
##  @param[in] param The parameters of the external field
##

proc momResidentQuda*(mom: pointer; param: ptr QudaGaugeParam) {.
    importc: "momResidentQuda", header: "quda.h".}
##
##  Compute the gauge force and update the momentum field
##
##  @param[in,out] mom The momentum field to be updated
##  @param[in] sitelink The gauge field from which we compute the force
##  @param[in] input_path_buf[dim][num_paths][path_length]
##  @param[in] path_length One less that the number of links in a loop (e.g., 3 for a staple)
##  @param[in] loop_coeff Coefficients of the different loops in the Symanzik action
##  @param[in] num_paths How many contributions from path_length different "staples"
##  @param[in] max_length The maximum number of non-zero of links in any path in the action
##  @param[in] dt The integration step size (for MILC this is dt*beta/3)
##  @param[in] param The parameters of the external fields and the computation settings
##

proc computeGaugeForceQuda*(mom: pointer; sitelink: pointer;
                           input_path_buf: ptr ptr ptr cint; path_length: ptr cint;
                           loop_coeff: ptr cdouble; num_paths: cint;
                           max_length: cint; dt: cdouble;
                           qudaGaugeParam: ptr QudaGaugeParam): cint {.
    importc: "computeGaugeForceQuda", header: "quda.h".}
##
##  Compute the product of gauge links along a path and add to/overwrite the output field
##
##  @param[in,out] out The output field to be updated
##  @param[in] sitelink The gauge field from which we compute the products of gauge links
##  @param[in] input_path_buf[dim][num_paths][path_length]
##  @param[in] path_length One less that the number of links in a loop (e.g., 3 for a staple)
##  @param[in] loop_coeff Coefficients of the different loops in the Symanzik action
##  @param[in] num_paths How many contributions from path_length different "staples"
##  @param[in] max_length The maximum number of non-zero of links in any path in the action
##  @param[in] dt The integration step size (for MILC this is dt*beta/3)
##  @param[in] param The parameters of the external fields and the computation settings
##

proc computeGaugePathQuda*(`out`: pointer; sitelink: pointer;
                          input_path_buf: ptr ptr ptr cint; path_length: ptr cint;
                          loop_coeff: ptr cdouble; num_paths: cint; max_length: cint;
                          dt: cdouble; qudaGaugeParam: ptr QudaGaugeParam): cint {.
    importc: "computeGaugePathQuda", header: "quda.h".}
##
##  Compute the traces of products of gauge links along paths using the resident field
##
##  @param[in,out] traces The computed traces
##  @param[in] sitelink The gauge field from which we compute the products of gauge links
##  @param[in] path_length The number of links in each loop
##  @param[in] loop_coeff Multiplicative coefficients for each loop
##  @param[in] num_paths Total number of loops
##  @param[in] max_length The maximum number of non-zero of links in any path in the action
##  @param[in] factor An overall normalization factor
##

proc computeGaugeLoopTraceQuda*(traces: ptr double_complex;
                               input_path_buf: ptr ptr cint; path_length: ptr cint;
                               loop_coeff: ptr cdouble; num_paths: cint;
                               max_length: cint; factor: cdouble) {.
    importc: "computeGaugeLoopTraceQuda", header: "quda.h".}
##
##  Evolve the gauge field by step size dt, using the momentum field
##  I.e., Evalulate U(t+dt) = e(dt pi) U(t)
##
##  @param gauge The gauge field to be updated
##  @param momentum The momentum field
##  @param dt The integration step size step
##  @param conj_mom Whether to conjugate the momentum matrix
##  @param exact Whether to use an exact exponential or Taylor expand
##  @param param The parameters of the external fields and the computation settings
##

proc updateGaugeFieldQuda*(gauge: pointer; momentum: pointer; dt: cdouble;
                          conj_mom: cint; exact: cint; param: ptr QudaGaugeParam) {.
    importc: "updateGaugeFieldQuda", header: "quda.h".}
##
##  Apply the staggered phase factors to the gauge field.  If the
##  imaginary chemical potential is non-zero then the phase factor
##  exp(imu/T) will be applied to the links in the temporal
##  direction.
##
##  @param gauge_h The gauge field
##  @param param The parameters of the gauge field
##

proc staggeredPhaseQuda*(gauge_h: pointer; param: ptr QudaGaugeParam) {.
    importc: "staggeredPhaseQuda", header: "quda.h".}
##
##  Project the input field on the SU(3) group.  If the target
##  tolerance is not met, this routine will give a runtime error.
##
##  @param gauge_h The gauge field to be updated
##  @param tol The tolerance to which we iterate
##  @param param The parameters of the gauge field
##

proc projectSU3Quda*(gauge_h: pointer; tol: cdouble; param: ptr QudaGaugeParam) {.
    importc: "projectSU3Quda", header: "quda.h".}
##
##  Evaluate the momentum contribution to the Hybrid Monte Carlo
##  action.
##
##  @param momentum The momentum field
##  @param param The parameters of the external fields and the computation settings
##  @return momentum action
##

proc momActionQuda*(momentum: pointer; param: ptr QudaGaugeParam): cdouble {.
    importc: "momActionQuda", header: "quda.h".}
##
##  Allocate a gauge (matrix) field on the device and optionally download a host gauge field.
##
##  @param gauge The host gauge field (optional - if set to 0 then the gauge field zeroed)
##  @param geometry The geometry of the matrix field to create (1 - scalar, 4 - vector, 6 - tensor)
##  @param param The parameters of the external field and the field to be created
##  @return Pointer to the gauge field (cast as a void*)
##

proc createGaugeFieldQuda*(gauge: pointer; geometry: cint; param: ptr QudaGaugeParam): pointer {.
    importc: "createGaugeFieldQuda", header: "quda.h".}
##
##  Copy the QUDA gauge (matrix) field on the device to the CPU
##
##  @param outGauge Pointer to the host gauge field
##  @param inGauge Pointer to the device gauge field (QUDA device field)
##  @param param The parameters of the host and device fields
##

proc saveGaugeFieldQuda*(outGauge: pointer; inGauge: pointer;
                        param: ptr QudaGaugeParam) {.importc: "saveGaugeFieldQuda",
    header: "quda.h".}
##
##  Reinterpret gauge as a pointer to a GaugeField and call destructor.
##
##  @param gauge Gauge field to be freed
##

proc destroyGaugeFieldQuda*(gauge: pointer) {.importc: "destroyGaugeFieldQuda",
    header: "quda.h".}
##
##  Compute the clover field and its inverse from the resident gauge field.
##
##  @param param The parameters of the clover field to create
##

proc createCloverQuda*(param: ptr QudaInvertParam) {.importc: "createCloverQuda",
    header: "quda.h".}
##
##  Compute the clover force contributions from a set of partial
##  fractions stemming from a rational approximation suitable for use
##  within MILC.
##
##  @param mom Force matrix
##  @param dt Integrating step size
##  @param x Array of solution vectors
##  @param p Array of intermediate vectors
##  @param coeff Array of residues for each contribution (multiplied by stepsize)
##  @param kappa2 -kappa*kappa parameter
##  @param ck -clover_coefficient * kappa / 8
##  @param nvec Number of vectors
##  @param multiplicity Number fermions this bilinear reresents
##  @param gauge Gauge Field
##  @param gauge_param Gauge field meta data
##  @param inv_param Dirac and solver meta data
##

proc computeCloverForceQuda*(mom: pointer; dt: cdouble; x: ptr pointer; p: ptr pointer;
                            coeff: ptr cdouble; kappa2: cdouble; ck: cdouble;
                            nvector: cint; multiplicity: cdouble; gauge: pointer;
                            gauge_param: ptr QudaGaugeParam;
                            inv_param: ptr QudaInvertParam) {.
    importc: "computeCloverForceQuda", header: "quda.h".}
##
##  Compute the force from a clover or twisted clover determinant or
##  a set of partial fractions stemming from a rational approximation
##  suitable for use from within tmLQCD.
##
##  @param h_mom Host force matrix
##  @param h_x Array of solution vectors x_i = ( Q^2 + s_i )^{-1} b
##  @param h_x0 Array of source vector necessary to compute the force of a ratio of determinant
##  @param coeff Array of coefficients for the rational approximation or {1.0} for the determinant.
##  @param nvector Number of solution vectors and coefficients
##  @param gauge_param Gauge field meta data
##  @param inv_param Dirac and solver meta data
##  @param detratio if 0 compute the force of a determinant otherwise compute the force from a ratio of determinants
##

proc computeTMCloverForceQuda*(h_mom: pointer; h_x: ptr pointer; h_x0: ptr pointer;
                              coeff: ptr cdouble; nvector: cint;
                              gauge_param: ptr QudaGaugeParam;
                              inv_param: ptr QudaInvertParam; detratio: cint) {.
    importc: "computeTMCloverForceQuda", header: "quda.h".}
##
##  Compute the naive staggered force.  All fields must be in the same precision.
##
##  @param mom Momentum field
##  @param dt Integrating step size
##  @param delta Additional scale factor when updating momentum (mom += delta * [force]_TA
##  @param gauge Gauge field (at present only supports resident gauge field)
##  @param x Array of single-parity solution vectors (at present only supports resident solutions)
##  @param gauge_param Gauge field meta data
##  @param invert_param Dirac and solver meta data
##

proc computeStaggeredForceQuda*(mom: pointer; dt: cdouble; delta: cdouble;
                               gauge: pointer; x: ptr pointer;
                               gauge_param: ptr QudaGaugeParam;
                               invert_param: ptr QudaInvertParam) {.
    importc: "computeStaggeredForceQuda", header: "quda.h".}
##
##  Compute the fermion force for the HISQ quark action and integrate the momentum.
##  @param momentum        The momentum field we are integrating
##  @param dt              The stepsize used to integrate the momentum
##  @param level2_coeff    The coefficients for the second level of smearing in the quark action.
##  @param fat7_coeff      The coefficients for the first level of smearing (fat7) in the quark action.
##  @param w_link          Unitarized link variables obtained by applying fat7 smearing and unitarization to the original links.
##  @param v_link          Fat7 link variables.
##  @param u_link          SU(3) think link variables.
##  @param quark           The input fermion field.
##  @param num             The number of quark fields
##  @param num_naik        The number of naik contributions
##  @param coeff           The coefficient multiplying the fermion fields in the outer product
##  @param param.          The field parameters.
##

proc computeHISQForceQuda*(momentum: pointer; dt: cdouble;
                          level2_coeff: array[6, cdouble];
                          fat7_coeff: array[6, cdouble]; w_link: pointer;
                          v_link: pointer; u_link: pointer; quark: ptr pointer;
                          num: cint; num_naik: cint; coeff: ptr ptr cdouble;
                          param: ptr QudaGaugeParam) {.
    importc: "computeHISQForceQuda", header: "quda.h".}
##
##      @brief Generate Gaussian distributed fields and store in the
##      resident gauge field. We create a Gaussian-distributed su(n)
##      field and exponentiate it, e.g., U = exp(sigma * H), where H is
##      the distributed su(n) field and sigma is the width of the
##      distribution (sigma = 0 results in a free field, and sigma = 1 has
##      maximum disorder).
##
##      @param seed The seed used for the RNG
##      @param sigma Width of Gaussian distrubution
##

proc gaussGaugeQuda*(seed: culonglong; sigma: cdouble) {.importc: "gaussGaugeQuda",
    header: "quda.h".}
##
##  @brief Generate Gaussian distributed fields and store in the
##  resident momentum field. We create a Gaussian-distributed su(n)
##  field, e.g., sigma * H, where H is the distributed su(n) field
##  and sigma is the width of the distribution (sigma = 0 results
##  in a free field, and sigma = 1 has maximum disorder).
##
##  @param seed The seed used for the RNG
##  @param sigma Width of Gaussian distrubution
##

proc gaussMomQuda*(seed: culonglong; sigma: cdouble) {.importc: "gaussMomQuda",
    header: "quda.h".}
##
##  Computes the total, spatial and temporal plaquette averages of the loaded gauge configuration.
##  @param[out] Array for storing the averages (total, spatial, temporal)
##

proc plaqQuda*(plaq: array[3, cdouble]) {.importc: "plaqQuda", header: "quda.h".}
##
##      @brief Computes the trace of the Polyakov loop of the current resident field
##      in a given direction.
##
##      @param[out] ploop Trace of the Polyakov loop in direction dir
##      @param[in] dir Direction of Polyakov loop
##

proc polyakovLoopQuda*(ploop: array[2, cdouble]; dir: cint) {.
    importc: "polyakovLoopQuda", header: "quda.h".}
##
##  Performs a deep copy from the internal extendedGaugeResident field.
##  @param Pointer to externally allocated GaugeField
##

proc copyExtendedResidentGaugeQuda*(resident_gauge: pointer) {.
    importc: "copyExtendedResidentGaugeQuda", header: "quda.h".}
##
##  Performs gaussian/Wuppertal smearing on a given spinor using the gauge field
##  gaugeSmeared, if it exist, or gaugePrecise if no smeared field is present.
##  @param h_in   Input spinor field
##  @param h_out  Output spinor field
##  @param param  Contains all metadata regarding host and device
##                storage and operator which will be applied to the spinor
##  @param n_steps Number of steps to apply.
##  @param coeff  Width of the Gaussian distribution
##  @param smear_type Gaussian/Wuppertal smearing
##

proc performFermionSmearQuda*(h_out: pointer; h_in: pointer;
                             param: ptr QudaInvertParam; n_steps: cint;
                             coeff: cdouble; smear_type: QudaFermionSmearType) {.
    importc: "performFermionSmearQuda", header: "quda.h".}
##
##  LEGACY
##  Performs Wuppertal smearing on a given spinor using the gauge field
##  gaugeSmeared, if it exist, or gaugePrecise if no smeared field is present.
##  @param h_out  Result spinor field
##  @param h_in   Input spinor field
##  @param param  Contains all metadata regarding host and device
##                storage and operator which will be applied to the spinor
##  @param n_steps Number of steps to apply.
##  @param alpha  Alpha coefficient for Wuppertal smearing.
##

proc performWuppertalnStep*(h_out: pointer; h_in: pointer;
                           param: ptr QudaInvertParam; n_steps: cuint; alpha: cdouble) {.
    importc: "performWuppertalnStep", header: "quda.h".}
##
##  LEGACY
##  Performs gaussian smearing on a given spinor using the gauge field
##  gaugeSmeared, if it exist, or gaugePrecise if no smeared field is present.
##  @param h_in   Input spinor field
##  @param h_out  Output spinor field
##  @param param  Contains all metadata regarding host and device
##                storage and operator which will be applied to the spinor
##  @param n_steps Number of steps to apply.
##  @param omega  Width of the Gaussian distribution
##

proc performGaussianSmearNStep*(h_out: pointer; h_in: pointer;
                               param: ptr QudaInvertParam; n_steps: cint;
                               omega: cdouble) {.
    importc: "performGaussianSmearNStep", header: "quda.h".}
##
##  Performs APE, Stout, or Over Imroved STOUT smearing on gaugePrecise and stores it in gaugeSmeared
##  @param[in] smear_param Parameter struct that defines the computation parameters
##  @param[in,out] obs_param Parameter struct that defines which
##  observables we are making and the resulting observables.
##

proc performGaugeSmearQuda*(smear_param: ptr QudaGaugeSmearParam;
                           obs_param: ptr QudaGaugeObservableParam) {.
    importc: "performGaugeSmearQuda", header: "quda.h".}
##
##  Performs Wilson Flow on gaugePrecise and stores it in gaugeSmeared
##  @param[in] smear_param Parameter struct that defines the computation parameters
##  @param[in,out] obs_param Parameter struct that defines which
##  observables we are making and the resulting observables.
##

proc performWFlowQuda*(smear_param: ptr QudaGaugeSmearParam;
                      obs_param: ptr QudaGaugeObservableParam) {.
    importc: "performWFlowQuda", header: "quda.h".}
##
##  Performs Gradient Flow (gauge + fermion) on gaugePrecise and stores it in gaugeSmeared
##  @param[out] h_out Output fermion field set
##  @param[in] h_in Input fermion field set
##  @param[in] inv_param Dirac/Laplacian and solver meta data
##  @param[in] smear_param Parameter struct that defines the computation parameters
##  @param[in,out] obs_param Parameter struct that defines which
##  observables we are making and the resulting observables.
##  @param[in] nSpinors Number of spinors in the input and output fields
##

proc performGFlowQuda*(h_out: ptr pointer; h_in: ptr pointer;
                      inv_param: ptr QudaInvertParam;
                      smear_param: ptr QudaGaugeSmearParam;
                      obs_param: ptr QudaGaugeObservableParam; nSpinors: csize_t) {.
    importc: "performGFlowQuda", header: "quda.h".}
##
##  Performs Adjoint Gradient Flow (gauge + fermion) the "safe" way on gaugePrecise and stores it in gaugeSmeared
##  @param[out] h_out Output fermion field set
##  @param[in] h_in Input fermion field set
##  @param[in] inv_param Dirac/Laplacian and solver meta data
##  @param[in] smear_param Parameter struct that defines the computation parameters
##  @param[in,out] obs_param Parameter struct that defines which
##  observables we are making and the resulting observables.
##  @param[in] nSpinors Number of spinors in the input and output fields
##

proc performAdjGFlowSafe*(h_out: ptr pointer; h_in: ptr pointer;
                         inv_param: ptr QudaInvertParam;
                         smear_param: ptr QudaGaugeSmearParam; nSpinors: csize_t) {.
    importc: "performAdjGFlowSafe", header: "quda.h".}
##
##  Performs Adjoint Gradient Flow (gauge + fermion) the Hierarchical way on gaugePrecise and stores it in gaugeSmeared
##  @param[out] h_out Output fermion field set
##  @param[in] h_in Input fermion field set
##  @param[in] inv_param Dirac/Laplacian and solver meta data
##  @param[in] smear_param Parameter struct that defines the computation parameters
##  @param[in,out] obs_param Parameter struct that defines which
##  observables we are making and the resulting observables.
##  @param[in] nSpinors Number of spinors in the input and output fields
##

proc performAdjGFlowHier*(h_out: ptr pointer; h_in: ptr pointer;
                         inv_param: ptr QudaInvertParam;
                         smear_param: ptr QudaGaugeSmearParam; nSpinors: csize_t) {.
    importc: "performAdjGFlowHier", header: "quda.h".}
##
##  @brief Calculates a variety of gauge-field observables.  If a
##  smeared gauge field is presently loaded (in gaugeSmeared) the
##  observables are computed on this, else the resident gauge field
##  will be used.
##  @param[in,out] param Parameter struct that defines which
##  observables we are making and the resulting observables.
##

proc gaugeObservablesQuda*(param: ptr QudaGaugeObservableParam) {.
    importc: "gaugeObservablesQuda", header: "quda.h".}
##
##  Public function to perform color contractions of the host spinors x and y.
##  @param[in] x pointer to host data
##  @param[in] y pointer to host data
##  @param[out] result pointer to the 16 spin projections per lattice site
##  @param[in] cType Which type of contraction (open, degrand-rossi, etc)
##  @param[in] param meta data for construction of ColorSpinorFields.
##  @param[in] X spacetime data for construction of ColorSpinorFields.
##

proc contractQuda*(x: pointer; y: pointer; result: pointer; cType: QudaContractType;
                  param: ptr QudaInvertParam; X: ptr cint) {.importc: "contractQuda",
    header: "quda.h".}
##
##  @param[in] x pointer to host data array
##  @param[in] y pointer to host data array
##  @param[out] result pointer to the spin*spin projections per lattice slice site
##  @param[in] cType Which type of contraction (open, degrand-rossi, etc)
##  @param[in] param meta data for construction of ColorSpinorFields.
##  @param[in] src_colors color dilution parameter
##  @param[in] X local lattice dimansions
##  @param[in] source_position source position array
##  @param[in] number of momentum modes
##  @param[in] mom_modes momentum modes
##  @param[in] fft_type Fourier phase factor type (cos, sin or exp{ikx})
##

proc contractFTQuda*(x: ptr pointer; y: ptr pointer; result: ptr pointer;
                    cType: QudaContractType; cs_param_ptr: pointer;
                    src_colors: cint; X: ptr cint; source_position: ptr cint;
                    n_mom: cint; mom_modes: ptr cint; fft_type: ptr QudaFFTSymmType) {.
    importc: "contractFTQuda", header: "quda.h".}
##
##  @brief Gauge fixing with overrelaxation with support for single and multi GPU.
##  @param[in,out] gauge, gauge field to be fixed
##  @param[in] gauge_dir, 3 for Coulomb gauge fixing, other for Landau gauge fixing
##  @param[in] Nsteps, maximum number of steps to perform gauge fixing
##  @param[in] verbose_interval, print gauge fixing info when iteration count is a multiple of this
##  @param[in] relax_boost, gauge fixing parameter of the overrelaxation method, most common value is 1.5 or 1.7.
##  @param[in] tolerance, torelance value to stop the method, if this value is zero then the method stops when
##  iteration reachs the maximum number of steps defined by Nsteps
##  @param[in] reunit_interval, reunitarize gauge field when iteration count is a multiple of this
##  @param[in] stopWtheta, 0 for MILC criterion and 1 to use the theta value
##  @param[in] param The parameters of the external fields and the computation settings
##

proc computeGaugeFixingOVRQuda*(gauge: pointer; gauge_dir: cuint; Nsteps: cuint;
                               verbose_interval: cuint; relax_boost: cdouble;
                               tolerance: cdouble; reunit_interval: cuint;
                               stopWtheta: cuint; param: ptr QudaGaugeParam): cint {.
    importc: "computeGaugeFixingOVRQuda", header: "quda.h".}
##
##  @brief Gauge fixing with Steepest descent method with FFTs with support for single GPU only.
##  @param[in,out] gauge, gauge field to be fixed
##  @param[in] gauge_dir, 3 for Coulomb gauge fixing, other for Landau gauge fixing
##  @param[in] Nsteps, maximum number of steps to perform gauge fixing
##  @param[in] verbose_interval, print gauge fixing info when iteration count is a multiple of this
##  @param[in] alpha, gauge fixing parameter of the method, most common value is 0.08
##  @param[in] autotune, 1 to autotune the method, i.e., if the Fg inverts its tendency we decrease the alpha value
##  @param[in] tolerance, torelance value to stop the method, if this value is zero then the method stops when
##  iteration reachs the maximum number of steps defined by Nsteps
##  @param[in] stopWtheta, 0 for MILC criterion and 1 to use the theta value
##  @param[in] param The parameters of the external fields and the computation settings
##

proc computeGaugeFixingFFTQuda*(gauge: pointer; gauge_dir: cuint; Nsteps: cuint;
                               verbose_interval: cuint; alpha: cdouble;
                               autotune: cuint; tolerance: cdouble;
                               stopWtheta: cuint; param: ptr QudaGaugeParam): cint {.
    importc: "computeGaugeFixingFFTQuda", header: "quda.h".}
##
##  @brief Strided Batched GEMM
##  @param[in] arrayA The array containing the A matrix data
##  @param[in] arrayB The array containing the B matrix data
##  @param[in] arrayC The array containing the C matrix data
##  @param[in] native Boolean to use either the native or generic version
##  @param[in] param The data defining the problem execution.
##

proc blasGEMMQuda*(arrayA: pointer; arrayB: pointer; arrayC: pointer;
                  native: QudaBoolean; param: ptr QudaBLASParam) {.
    importc: "blasGEMMQuda", header: "quda.h".}
##
##  @brief Strided Batched in-place matrix inversion via LU
##  @param[in] Ainv The array containing the A inverse matrix data
##  @param[in] A The array containing the A matrix data
##  @param[in] use_native Boolean to use either the native or generic version
##  @param[in] param The data defining the problem execution.
##

proc blasLUInvQuda*(Ainv: pointer; A: pointer; use_native: QudaBoolean;
                   param: ptr QudaBLASParam) {.importc: "blasLUInvQuda",
    header: "quda.h".}
##
##  @brief Flush the chronological history for the given index
##  @param[in] index Index for which we are flushing
##

proc flushChronoQuda*(index: cint) {.importc: "flushChronoQuda", header: "quda.h".}
##
##  Create deflation solver resources.
##
##

proc newDeflationQuda*(param: ptr QudaEigParam): pointer {.
    importc: "newDeflationQuda", header: "quda.h".}
##
##  Free resources allocated by the deflated solver
##

proc destroyDeflationQuda*(df_instance: pointer) {.importc: "destroyDeflationQuda",
    header: "quda.h".}
##
##  @brief Flush the memory pools associated with the supplied type.
##  At present this only supports the options QUDA_MEMORY_DEVICE and
##  QUDA_MEMORY_HOST_PINNED, and any other type will result in an
##  error.
##  @param[in] type The memory type whose pool we wish to flush.
##

proc flushPoolQuda*(`type`: QudaMemoryType) {.importc: "flushPoolQuda",
    header: "quda.h".}
proc setMPICommHandleQuda*(mycomm: pointer) {.importc: "setMPICommHandleQuda",
    header: "quda.h".}
##  Parameter set for quark smearing operations

type
  QudaQuarkSmearParam* {.importc: "QudaQuarkSmearParam", header: "quda.h", bycopy.} = object
    ## -------------------------------------------------
    ##  Used to store information pertinent to the operator *
    inv_param* {.importc: "inv_param".}: ptr QudaInvertParam
    ##  Number of steps to apply *
    n_steps* {.importc: "n_steps".}: cint
    ##  The width of the Gaussian *
    width* {.importc: "width".}: cdouble
    ##  if nonzero then compute two-link, otherwise reuse gaugeSmeared*
    compute_2link* {.importc: "compute_2link".}: cint
    ##  if nonzero then delete two-link, otherwise keep two-link for future use*
    delete_2link* {.importc: "delete_2link".}: cint
    ##  Set if the input spinor is on a time slice *
    t0* {.importc: "t0".}: cint
    ##  Time taken for the smearing operations *
    secs* {.importc: "secs".}: cdouble
    ##  Flops count for the smearing operations *
    gflops* {.importc: "gflops".}: cdouble
    energy* {.importc: "energy".}: cdouble
    ## < The energy consumed by the smearing operations
    power* {.importc: "power".}: cdouble
    ## < The mean power of the smearing operations
    temp* {.importc: "temp".}: cdouble
    ## < The mean temperature of the device for the duration of the smearing operations
    clock* {.importc: "clock".}: cdouble
    ## < The mean clock frequency of the device for the duration of the smearing operations


##
##  Performs two-link Gaussian smearing on a given spinor (for staggered fermions).
##  @param[in,out] h_in Input spinor field to smear
##  @param[in] smear_param   Contains all metadata the operator which will be applied to the spinor
##

proc performTwoLinkGaussianSmearNStep*(h_in: pointer;
                                      smear_param: ptr QudaQuarkSmearParam) {.
    importc: "performTwoLinkGaussianSmearNStep", header: "quda.h".}
##
##  @brief Performs contractions between a set of quark fields and
##  eigenvectors of the 3-d Laplace operator.
##  @param[in,out] host_sinks An array representing the inner
##  products between the quark fields and the eigen-vector fields.
##  Ordered as [nQuark][nEv][Lt][nSpin][complexity].
##  @param[in] host_quark Array of quark fields we are taking the inner over
##  @param[in] n_quark Number of quark fields
##  @param[in] tile_quark Tile size for quark fields (batch size)
##  @param[in] host_evec Array of eigenvectors we are taking the inner over
##  @param[in] n_evec Number of eigenvectors
##  @param[in] tile_evec Tile size for eigenvectors (batch size)
##  @param[in] inv_param Meta-data structure
##  @param[in] X Lattice dimensions
##

## !!!Ignored construct:  void laphSinkProject ( double _Complex * host_sinks , void * * host_quark , int n_quark , int tile_quark , void * * host_evec , int nevec , int tile_evec , QudaInvertParam * inv_param , const int X [ 4 ] ) ;
## Error: token expected: ; but got: (!!!

##  remove NVRTC WAR

##  #include <quda_new_interface.h>
