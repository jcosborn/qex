import enum_quda, quda_constants
import base/stdUtils

## *
##  @file  quda.h
##  @brief Main header file for the QUDA library
##
##  Note to QUDA developers: When adding new members to QudaGaugeParam
##  and QudaInvertParam, be sure to update lib/check_params.h as well
##  as the Fortran interface in lib/quda_fortran.F90.
##

type
  #ConstInt* {.importc:"const int".} = cint
  double_complex* {.importc:"double _Complex".} = object
converter toDoubleComplex*(x: array[2,float]): double_complex =
  var r = cast[ptr array[2,float]](addr result)
  r[] = x
## *
##  Parameters having to do with the gauge field or the
##  interpretation of the gauge field by various Dirac operators
##

type
  QudaGaugeParam* {.importc: "QudaGaugeParam", header: "quda.h", bycopy.} = object
    location* {.importc: "location".}: QudaFieldLocation ## *< The location of the gauge field
    X* {.importc: "X".}: array[4, cint] ## *< The local space-time dimensions (without checkboarding)
    anisotropy* {.importc: "anisotropy".}: cdouble ## *< Used for Wilson and Wilson-clover
    tadpole_coeff* {.importc: "tadpole_coeff".}: cdouble ## *< Used for staggered only
    scale* {.importc: "scale".}: cdouble ## *< Used by staggered long links
    `type`* {.importc: "type".}: QudaLinkType ## *< The link type of the gauge field (e.g., Wilson, fat, long, etc.)
    gauge_order* {.importc: "gauge_order".}: QudaGaugeFieldOrder ## *< The ordering on the input gauge field
    t_boundary* {.importc: "t_boundary".}: QudaTboundary ## *< The temporal boundary condition that will be used for fermion fields
    cpu_prec* {.importc: "cpu_prec".}: QudaPrecision ## *< The precision used by the caller
    cuda_prec* {.importc: "cuda_prec".}: QudaPrecision ## *< The precision of the cuda gauge field
    reconstruct* {.importc: "reconstruct".}: QudaReconstructType ## *< The reconstruction type of the cuda gauge field
    cuda_prec_sloppy* {.importc: "cuda_prec_sloppy".}: QudaPrecision ## *< The precision of the sloppy gauge field
    reconstruct_sloppy* {.importc: "reconstruct_sloppy".}: QudaReconstructType ## *< The recontruction type of the sloppy gauge field
    cuda_prec_refinement_sloppy* {.importc: "cuda_prec_refinement_sloppy".}: QudaPrecision ## *< The precision of the sloppy gauge field for the refinement step in multishift
    reconstruct_refinement_sloppy* {.importc: "reconstruct_refinement_sloppy".}: QudaReconstructType ## *< The recontruction type of the sloppy gauge field for the refinement step in multishift
    cuda_prec_precondition* {.importc: "cuda_prec_precondition".}: QudaPrecision ## *< The precision of the preconditioner gauge field
    reconstruct_precondition* {.importc: "reconstruct_precondition".}: QudaReconstructType ## *< The recontruction type of the preconditioner gauge field
    cuda_prec_eigensolver* {.importc: "cuda_prec_eigensolver".}: QudaPrecision ## *< The precision of the eigensolver gauge field
    reconstruct_eigensolver* {.importc: "reconstruct_eigensolver".}: QudaReconstructType ## *< The recontruction type of the eigensolver gauge field
    gauge_fix* {.importc: "gauge_fix".}: QudaGaugeFixed ## *< Whether the input gauge field is in the axial gauge or not
    ga_pad* {.importc: "ga_pad".}: cint ## *< The pad size that the cudaGaugeField will use (default=0)
    site_ga_pad* {.importc: "site_ga_pad".}: cint ## *< Used by link fattening and the gauge and fermion forces
    staple_pad* {.importc: "staple_pad".}: cint ## *< Used by link fattening
    llfat_ga_pad* {.importc: "llfat_ga_pad".}: cint ## *< Used by link fattening
    mom_ga_pad* {.importc: "mom_ga_pad".}: cint ## *< Used by the gauge and fermion forces
    staggered_phase_type* {.importc: "staggered_phase_type".}: QudaStaggeredPhase ## *< Set the staggered phase type of the links
    staggered_phase_applied* {.importc: "staggered_phase_applied".}: cint ## *< Whether the staggered phase has already been applied to the links
    i_mu* {.importc: "i_mu".}: cdouble ## *< Imaginary chemical potential
    overlap* {.importc: "overlap".}: cint ## *< Width of overlapping domains
    overwrite_mom* {.importc: "overwrite_mom".}: cint ## *< When computing momentum, should we overwrite it or accumulate to to
    use_resident_gauge* {.importc: "use_resident_gauge".}: cint ## *< Use the resident gauge field as input
    use_resident_mom* {.importc: "use_resident_mom".}: cint ## *< Use the resident momentum field as input
    make_resident_gauge* {.importc: "make_resident_gauge".}: cint ## *< Make the result gauge field resident
    make_resident_mom* {.importc: "make_resident_mom".}: cint ## *< Make the result momentum field resident
    return_result_gauge* {.importc: "return_result_gauge".}: cint ## *< Return the result gauge field
    return_result_mom* {.importc: "return_result_mom".}: cint ## *< Return the result momentum field
    gauge_offset* {.importc: "gauge_offset".}: csize_t ## *< Offset into MILC site struct to the gauge field (only if gauge_order=MILC_SITE_GAUGE_ORDER)
    mom_offset* {.importc: "mom_offset".}: csize_t ## *< Offset into MILC site struct to the momentum field (only if gauge_order=MILC_SITE_GAUGE_ORDER)
    site_size* {.importc: "site_size".}: csize_t ## *< Size of MILC site struct (only if gauge_order=MILC_SITE_GAUGE_ORDER)


## *
##  Parameters relating to the solver and the choice of Dirac operator.
##

type
  QudaInvertParam* {.importc: "QudaInvertParam", header: "quda.h", bycopy.} = object
    input_location* {.importc: "input_location".}: QudaFieldLocation ## *< The location of the input field
    output_location* {.importc: "output_location".}: QudaFieldLocation ## *< The location of the output field
    dslash_type* {.importc: "dslash_type".}: QudaDslashType ## *< The Dirac Dslash type that is being used
    inv_type* {.importc: "inv_type".}: QudaInverterType ## *< Which linear solver to use
    mass* {.importc: "mass".}: cdouble ## *< Used for staggered only
    kappa* {.importc: "kappa".}: cdouble ## *< Used for Wilson and Wilson-clover
    m5* {.importc: "m5".}: cdouble ## *< Domain wall height
    Ls* {.importc: "Ls".}: cint  ## *< Extent of the 5th dimension (for domain wall)
    b_5* {.importc: "b_5".}: array[QUDA_MAX_DWF_LS, double_complex] ## *< Mobius coefficients - only real part used if regular Mobius
    c_5* {.importc: "c_5".}: array[QUDA_MAX_DWF_LS, double_complex] ## *< Mobius coefficients - only real part used if regular Mobius
                                                               ## *<
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
    mu* {.importc: "mu".}: cdouble ## *< Twisted mass parameter
    epsilon* {.importc: "epsilon".}: cdouble ## *< Twisted mass parameter
    twist_flavor* {.importc: "twist_flavor".}: QudaTwistFlavorType ## *< Twisted mass flavor
    laplace3D* {.importc: "laplace3D".}: cint ## *< omit this direction from laplace operator: x,y,z,t -> 0,1,2,3 (-1 is full 4D)
    tol* {.importc: "tol".}: cdouble ## *< Solver tolerance in the L2 residual norm
    tol_restart* {.importc: "tol_restart".}: cdouble ## *< Solver tolerance in the L2 residual norm (used to restart InitCG)
    tol_hq* {.importc: "tol_hq".}: cdouble ## *< Solver tolerance in the heavy quark residual norm
    compute_true_res* {.importc: "compute_true_res".}: cint ## * Whether to compute the true residual post solve
    true_res* {.importc: "true_res".}: cdouble ## *< Actual L2 residual norm achieved in solver
    true_res_hq* {.importc: "true_res_hq".}: cdouble ## *< Actual heavy quark residual norm achieved in solver
    maxiter* {.importc: "maxiter".}: cint ## *< Maximum number of iterations in the linear solver
    reliable_delta* {.importc: "reliable_delta".}: cdouble ## *< Reliable update tolerance
    reliable_delta_refinement* {.importc: "reliable_delta_refinement".}: cdouble ## *< Reliable update tolerance used in post multi-shift solver refinement
    use_alternative_reliable* {.importc: "use_alternative_reliable".}: cint ## *< Whether to use alternative reliable updates
    use_sloppy_partial_accumulator* {.importc: "use_sloppy_partial_accumulator".}: cint ## *< Whether to keep the partial solution accumuator in sloppy precision
                                                                                    ## *< This parameter determines how often we accumulate into the
                                                                                    ##        solution vector from the direction vectors in the solver.
                                                                                    ##        E.g., running with solution_accumulator_pipeline = 4, means we
                                                                                    ##        will update the solution vector every four iterations using the
                                                                                    ##        direction vectors from the prior four iterations.  This
                                                                                    ##        increases performance of mixed-precision solvers since it means
                                                                                    ##        less high-precision vector round-trip memory travel, but
                                                                                    ##        requires more low-precision memory allocation.
    solution_accumulator_pipeline* {.importc: "solution_accumulator_pipeline".}: cint ## *< This parameter determines how many consecutive reliable update
                                                                                  ##     residual increases we tolerate before terminating the solver,
                                                                                  ##     i.e., how long do we want to keep trying to converge
    max_res_increase* {.importc: "max_res_increase".}: cint ## *< This parameter determines how many total reliable update
                                                        ##     residual increases we tolerate before terminating the solver,
                                                        ##     i.e., how long do we want to keep trying to converge
    max_res_increase_total* {.importc: "max_res_increase_total".}: cint ## *< This parameter determines how many consecutive heavy-quark
                                                                    ##     residual increases we tolerate before terminating the solver,
                                                                    ##     i.e., how long do we want to keep trying to converge
    max_hq_res_increase* {.importc: "max_hq_res_increase".}: cint ## *< This parameter determines how many total heavy-quark residual
                                                              ##     restarts we tolerate before terminating the solver, i.e., how long
                                                              ##     do we want to keep trying to converge
    max_hq_res_restart_total* {.importc: "max_hq_res_restart_total".}: cint ## *< After how many iterations shall the heavy quark residual be updated
    heavy_quark_check* {.importc: "heavy_quark_check".}: cint
    pipeline* {.importc: "pipeline".}: cint ## *< Whether to use a pipelined solver with less global sums
    num_offset* {.importc: "num_offset".}: cint ## *< Number of offsets in the multi-shift solver
    num_src* {.importc: "num_src".}: cint ## *< Number of sources in the multiple source solver
    num_src_per_sub_partition* {.importc: "num_src_per_sub_partition".}: cint ## *< Number of sources in the multiple source solver, but per sub-partition
                                                                          ## *< The grid of sub-partition according to which the processor grid will be partitioned.
                                                                          ##     Should have:
                                                                          ##       split_grid[0] * split_grid[1] * split_grid[2] * split_grid[3] * num_src_per_sub_partition == num_src. *
    split_grid* {.importc: "split_grid".}: array[QUDA_MAX_DIM, cint]
    overlap* {.importc: "overlap".}: cint ## *< Width of domain overlaps
                                      ## * Offsets for multi-shift solver
    offset* {.importc: "offset".}: array[QUDA_MAX_MULTI_SHIFT, cdouble] ## * Solver tolerance for each offset
    tol_offset* {.importc: "tol_offset".}: array[QUDA_MAX_MULTI_SHIFT, cdouble] ## *
                                                                           ## Solver
                                                                           ## tolerance for each shift when
                                                                           ## refinement is
                                                                           ## applied using the
                                                                           ## heavy-quark
                                                                           ## residual
    tol_hq_offset* {.importc: "tol_hq_offset".}: array[QUDA_MAX_MULTI_SHIFT, cdouble] ##
                                                                                 ## *
                                                                                 ## Actual
                                                                                 ## L2
                                                                                 ## residual
                                                                                 ## norm
                                                                                 ## achieved
                                                                                 ## in
                                                                                 ## solver
                                                                                 ## for
                                                                                 ## each
                                                                                 ## offset
    true_res_offset* {.importc: "true_res_offset".}: array[QUDA_MAX_MULTI_SHIFT,
        cdouble]              ## * Iterated L2 residual norm achieved in multi shift solver for each offset
    iter_res_offset* {.importc: "iter_res_offset".}: array[QUDA_MAX_MULTI_SHIFT,
        cdouble]              ## * Actual heavy quark residual norm achieved in solver for each offset
    true_res_hq_offset* {.importc: "true_res_hq_offset".}: array[
        QUDA_MAX_MULTI_SHIFT, cdouble] ## * Residuals in the partial faction expansion
    residue* {.importc: "residue".}: array[QUDA_MAX_MULTI_SHIFT, cdouble] ## * Whether we should evaluate the action after the linear solver
    compute_action* {.importc: "compute_action".}: cint ## * Computed value of the bilinear action (complex-valued)
                                                    ## 	invert: \phi^\dagger A^{-1} \phi
                                                    ## 	multishift: \phi^\dagger r(x) \phi = \phi^\dagger (sum_k residue[k] * (A + offset[k])^{-1} ) \phi
    action* {.importc: "action".}: array[2, cdouble]
    solution_type* {.importc: "solution_type".}: QudaSolutionType ## *< Type of system to solve
    solve_type* {.importc: "solve_type".}: QudaSolveType ## *< How to solve it
    matpc_type* {.importc: "matpc_type".}: QudaMatPCType ## *< The preconditioned matrix type
    dagger* {.importc: "dagger".}: QudaDagType ## *< Whether we are using the Hermitian conjugate system or not
    mass_normalization* {.importc: "mass_normalization".}: QudaMassNormalizationT ## *< The mass normalization is being used by the caller
    solver_normalization* {.importc: "solver_normalization".}: QudaSolverNormalization ## *< The normalization desired in the solver
    preserve_source* {.importc: "preserve_source".}: QudaPreserveSource ## *< Preserve the source or not in the linear solver (deprecated)
    cpu_prec* {.importc: "cpu_prec".}: QudaPrecision ## *< The precision used by the input fermion fields
    cuda_prec* {.importc: "cuda_prec".}: QudaPrecision ## *< The precision used by the QUDA solver
    cuda_prec_sloppy* {.importc: "cuda_prec_sloppy".}: QudaPrecision ## *< The precision used by the QUDA sloppy operator
    cuda_prec_refinement_sloppy* {.importc: "cuda_prec_refinement_sloppy".}: QudaPrecision ## *< The precision of the sloppy gauge field for the refinement step in multishift
    cuda_prec_precondition* {.importc: "cuda_prec_precondition".}: QudaPrecision ## *< The precision used by the QUDA preconditioner
    cuda_prec_eigensolver* {.importc: "cuda_prec_eigensolver".}: QudaPrecision ## *< The precision used by the QUDA eigensolver
    dirac_order* {.importc: "dirac_order".}: QudaDiracFieldOrder ## *< The order of the input and output fermion fields
    gamma_basis* {.importc: "gamma_basis".}: QudaGammaBasis ## *< Gamma basis of the input and output host fields
    clover_location* {.importc: "clover_location".}: QudaFieldLocation ## *< The location of the clover field
    clover_cpu_prec* {.importc: "clover_cpu_prec".}: QudaPrecision ## *< The precision used for the input clover field
    clover_cuda_prec* {.importc: "clover_cuda_prec".}: QudaPrecision ## *< The precision used for the clover field in the QUDA solver
    clover_cuda_prec_sloppy* {.importc: "clover_cuda_prec_sloppy".}: QudaPrecision ## *< The precision used for the clover field in the QUDA sloppy operator
    clover_cuda_prec_refinement_sloppy* {.
        importc: "clover_cuda_prec_refinement_sloppy".}: QudaPrecision ## *< The precision of the sloppy clover field for the refinement step in multishift
    clover_cuda_prec_precondition* {.importc: "clover_cuda_prec_precondition".}: QudaPrecision ## *< The precision used for the clover field in the QUDA preconditioner
    clover_cuda_prec_eigensolver* {.importc: "clover_cuda_prec_eigensolver".}: QudaPrecision ## *< The precision used for the clover field in the QUDA eigensolver
    clover_order* {.importc: "clover_order".}: QudaCloverFieldOrder ## *< The order of the input clover field
    use_init_guess* {.importc: "use_init_guess".}: QudaUseInitGuess ## *< Whether to use an initial guess in the solver or not
    clover_coeff* {.importc: "clover_coeff".}: cdouble ## *< Coefficient of the clover term
    clover_rho* {.importc: "clover_rho".}: cdouble ## *< Real number added to the clover diagonal (not to inverse)
    compute_clover_trlog* {.importc: "compute_clover_trlog".}: cint ## *< Whether to compute the trace log of the clover term
    trlogA* {.importc: "trlogA".}: array[2, cdouble] ## *< The trace log of the clover term (even/odd computed separately)
    compute_clover* {.importc: "compute_clover".}: cint ## *< Whether to compute the clover field
    compute_clover_inverse* {.importc: "compute_clover_inverse".}: cint ## *< Whether to compute the clover inverse field
    return_clover* {.importc: "return_clover".}: cint ## *< Whether to copy back the clover matrix field
    return_clover_inverse* {.importc: "return_clover_inverse".}: cint ## *< Whether to copy back the inverted clover matrix field
    verbosity* {.importc: "verbosity".}: QudaVerbosity ## *< The verbosity setting to use in the solver
    sp_pad* {.importc: "sp_pad".}: cint ## *< The padding to use for the fermion fields
    cl_pad* {.importc: "cl_pad".}: cint ## *< The padding to use for the clover fields
    iter* {.importc: "iter".}: cint ## *< The number of iterations performed by the solver
    gflops* {.importc: "gflops".}: cdouble ## *< The Gflops rate of the solver
    secs* {.importc: "secs".}: cdouble ## *< The time taken by the solver
    tune* {.importc: "tune".}: QudaTune ## *< Enable auto-tuning? (default = QUDA_TUNE_YES)
                                    ## * Number of steps in s-step algorithms
    Nsteps* {.importc: "Nsteps".}: cint ## * Maximum size of Krylov space used by solver
    gcrNkrylov* {.importc: "gcrNkrylov".}: cint ##
                                            ##  The following parameters are related to the solver
                                            ##  preconditioner, if enabled.
                                            ##
                                            ## *
                                            ##  The inner Krylov solver used in the preconditioner.  Set to
                                            ##  QUDA_INVALID_INVERTER to disable the preconditioner entirely.
                                            ##
    inv_type_precondition* {.importc: "inv_type_precondition".}: QudaInverterType ## *
                                                                              ## Preconditioner
                                                                              ## instance,
                                                                              ## e.g.,
                                                                              ## multigrid
    preconditioner* {.importc: "preconditioner".}: pointer ## * Deflation instance
    deflation_op* {.importc: "deflation_op".}: pointer ## * defines deflation
    eig_param* {.importc: "eig_param".}: pointer ## * If true, deflate the initial guess
    deflate* {.importc: "deflate".}: QudaBoolean ## * Dirac Dslash used in preconditioner
    dslash_type_precondition* {.importc: "dslash_type_precondition".}: QudaDslashType ##
                                                                                  ## *
                                                                                  ## Verbosity
                                                                                  ## of
                                                                                  ## the
                                                                                  ## inner
                                                                                  ## Krylov
                                                                                  ## solver
    verbosity_precondition* {.importc: "verbosity_precondition".}: QudaVerbosity ## *
                                                                             ## Tolerance in the
                                                                             ## inner
                                                                             ## solver
    tol_precondition* {.importc: "tol_precondition".}: cdouble ## * Maximum number of iterations allowed in the inner solver
    maxiter_precondition* {.importc: "maxiter_precondition".}: cint ## * Relaxation parameter used in GCR-DD (default = 1.0)
    omega* {.importc: "omega".}: cdouble ## * Basis for CA algorithms
    ca_basis* {.importc: "ca_basis".}: QudaCABasis ## * Minimum eigenvalue for Chebyshev CA basis
    ca_lambda_min* {.importc: "ca_lambda_min".}: cdouble ## * Maximum eigenvalue for Chebyshev CA basis
    ca_lambda_max* {.importc: "ca_lambda_max".}: cdouble ## * Number of preconditioner cycles to perform per iteration
    precondition_cycle* {.importc: "precondition_cycle".}: cint ## * Whether to use additive or multiplicative Schwarz preconditioning
    schwarz_type* {.importc: "schwarz_type".}: QudaSchwarzType ## *
                                                           ##  Whether to use the L2 relative residual, Fermilab heavy-quark
                                                           ##  residual, or both to determine convergence.  To require that both
                                                           ##  stopping conditions are satisfied, use a bitwise OR as follows:
                                                           ##
                                                           ##  p.residual_type = (QudaResidualType) (QUDA_L2_RELATIVE_RESIDUAL
                                                           ##                                      | QUDA_HEAVY_QUARK_RESIDUAL);
                                                           ##
    residual_type* {.importc: "residual_type".}: QudaResidualType ## *Parameters for deflated solvers
                                                              ## * The precision of the Ritz vectors
    cuda_prec_ritz* {.importc: "cuda_prec_ritz".}: QudaPrecision ## * How many vectors to compute after one solve
                                                             ##   for eigCG recommended values 8 or 16
                                                             ##
    n_ev* {.importc: "n_ev".}: cint ## * EeigCG  : Search space dimension
                                ##   gmresdr : Krylov subspace dimension
                                ##
    max_search_dim* {.importc: "max_search_dim".}: cint ## * For systems with many RHS: current RHS index
    rhs_idx* {.importc: "rhs_idx".}: cint ## * Specifies deflation space volume: total number of eigenvectors is n_ev*deflation_grid
    deflation_grid* {.importc: "deflation_grid".}: cint ## * eigCG: selection criterion for the reduced eigenvector set
    eigenval_tol* {.importc: "eigenval_tol".}: cdouble ## * mixed precision eigCG tuning parameter:  minimum search vector space restarts
    eigcg_max_restarts* {.importc: "eigcg_max_restarts".}: cint ## * initCG tuning parameter:  maximum restarts
    max_restart_num* {.importc: "max_restart_num".}: cint ## * initCG tuning parameter:  tolerance for cg refinement corrections in the deflation stage
    inc_tol* {.importc: "inc_tol".}: cdouble ## * Whether to make the solution vector(s) after the solve
    make_resident_solution* {.importc: "make_resident_solution".}: cint ## * Whether to use the resident solution vector(s)
    use_resident_solution* {.importc: "use_resident_solution".}: cint ## * Whether to use the solution vector to augment the chronological basis
    chrono_make_resident* {.importc: "chrono_make_resident".}: cint ## * Whether the solution should replace the last entry in the chronology
    chrono_replace_last* {.importc: "chrono_replace_last".}: cint ## * Whether to use the resident chronological basis
    chrono_use_resident* {.importc: "chrono_use_resident".}: cint ## * The maximum length of the chronological history to store
    chrono_max_dim* {.importc: "chrono_max_dim".}: cint ## * The index to indicate which chrono history we are augmenting
    chrono_index* {.importc: "chrono_index".}: cint ## * Precision to store the chronological basis in
    chrono_precision* {.importc: "chrono_precision".}: QudaPrecision ## * Which external library to use in the linear solvers (MAGMA or Eigen)
    extlib_type* {.importc: "extlib_type".}: QudaExtLibType ## * Whether to use the platform native or generic BLAS / LAPACK
    native_blas_lapack* {.importc: "native_blas_lapack".}: QudaBoolean


##  Parameter set for solving eigenvalue problems.

type
  QudaEigParam* {.importc: "QudaEigParam", header: "quda.h", bycopy.} = object
    invert_param* {.importc: "invert_param".}: ptr QudaInvertParam ##  EIGENSOLVER PARAMS
                                                              ## -------------------------------------------------
                                                              ## * Used to store information pertinent to the operator *
    ## * Type of eigensolver algorithm to employ *
    eig_type* {.importc: "eig_type".}: QudaEigType ## * Use Polynomial Acceleration *
    use_poly_acc* {.importc: "use_poly_acc".}: QudaBoolean ## * Degree of the Chebysev polynomial *
    poly_deg* {.importc: "poly_deg".}: cint ## * Range used in polynomial acceleration *
    a_min* {.importc: "a_min".}: cdouble
    a_max* {.importc: "a_max".}: cdouble ## * Whether to preserve the deflation space between solves.  If
                                     ##         true, the space will be stored in an instance of the
                                     ##         deflation_space struct, pointed to by preserve_deflation_space
    preserve_deflation* {.importc: "preserve_deflation".}: QudaBoolean ## * This is where we store the deflation space.  This will point
                                                                   ##         to an instance of
                                                                   ## deflation_space. When a deflated solver is enabled, the deflation space will be obtained from this.
    preserve_deflation_space* {.importc: "preserve_deflation_space".}: pointer ## * If we
                                                                           ## restore the
                                                                           ## deflation
                                                                           ## space, this
                                                                           ## boolean
                                                                           ## indicates
                                                                           ##
                                                                           ## whether we are also
                                                                           ## preserving the
                                                                           ## evalues or
                                                                           ## recomputing
                                                                           ##         them.  For
                                                                           ## example if a
                                                                           ## different mass shift is being used
                                                                           ##         than the one used to
                                                                           ## generate the
                                                                           ## space, then this
                                                                           ## should be
                                                                           ##
                                                                           ## false, but
                                                                           ## preserve_deflation would be true
    preserve_evals* {.importc: "preserve_evals".}: QudaBoolean ## * What type of Dirac operator we are using *
                                                           ## * If !(use_norm_op) && !(use_dagger) use M. *
                                                           ## * If use_dagger, use Mdag *
                                                           ## * If use_norm_op, use MdagM *
                                                           ## * If use_norm_op && use_dagger use MMdag. *
    use_dagger* {.importc: "use_dagger".}: QudaBoolean
    use_norm_op* {.importc: "use_norm_op".}: QudaBoolean ## * Use Eigen routines to eigensolve the upper Hessenberg via QR *
    use_eigen_qr* {.importc: "use_eigen_qr".}: QudaBoolean ## * Performs an MdagM solve, then constructs the left and right SVD. *
    compute_svd* {.importc: "compute_svd".}: QudaBoolean ## * Performs the \gamma_5 OP solve by Post multipling the eignvectors with
                                                     ##         \gamma_5 before computing the eigenvalues
    compute_gamma5* {.importc: "compute_gamma5".}: QudaBoolean ## * If true, the solver will error out if the convergence criteria are not met *
    require_convergence* {.importc: "require_convergence".}: QudaBoolean ## * Which part of the spectrum to solve *
    spectrum* {.importc: "spectrum".}: QudaEigSpectrumType ## * Size of the eigenvector search space *
    n_ev* {.importc: "n_ev".}: cint ## * Total size of Krylov space *
    n_kr* {.importc: "n_kr".}: cint ## * Max number of locked eigenpairs (deduced at runtime) *
    nLockedMax* {.importc: "nLockedMax".}: cint ## * Number of requested converged eigenvectors *
    n_conv* {.importc: "n_conv".}: cint ## * Number of requested converged eigenvectors to use in deflation *
    n_ev_deflate* {.importc: "n_ev_deflate".}: cint ## * Tolerance on the least well known eigenvalue's residual *
    tol* {.importc: "tol".}: cdouble ## * Tolerance on the QR iteration *
    qr_tol* {.importc: "qr_tol".}: cdouble ## * For IRLM/IRAM, check every nth restart *
    check_interval* {.importc: "check_interval".}: cint ## * For IRLM/IRAM, quit after n restarts *
    max_restarts* {.importc: "max_restarts".}: cint ## * For the Ritz rotation, the maximal number of extra vectors the solver may allocate *
    batched_rotate* {.importc: "batched_rotate".}: cint ## * For block method solvers, the block size *
    block_size* {.importc: "block_size".}: cint ## * In the test function, cross check the device result against ARPACK *
    arpack_check* {.importc: "arpack_check".}: QudaBoolean ## * For Arpack cross check, name of the Arpack logfile *
    arpack_logfile* {.importc: "arpack_logfile".}: array[512, char] ## * Name of the QUDA logfile (residua, upper
                                                               ## Hessenberg/tridiag matrix updates) *
    QUDA_logfile* {.importc: "QUDA_logfile".}: array[512, char] ## -------------------------------------------------
                                                           ##  EIG-CG PARAMS
                                                           ## -------------------------------------------------
    nk* {.importc: "nk".}: cint
    np* {.importc: "np".}: cint  ## * Whether to load eigenvectors
    import_vectors* {.importc: "import_vectors".}: QudaBoolean ## * The precision of the Ritz vectors
    cuda_prec_ritz* {.importc: "cuda_prec_ritz".}: QudaPrecision ## * The memory type used to keep the Ritz vectors
    mem_type_ritz* {.importc: "mem_type_ritz".}: QudaMemoryType ## * Location where deflation should be done
    location* {.importc: "location".}: QudaFieldLocation ## * Whether to run the verification checks once set up is complete
    run_verify* {.importc: "run_verify".}: QudaBoolean ## * Filename prefix where to load the null-space vectors
    vec_infile* {.importc: "vec_infile".}: array[256, char] ## * Filename prefix for where to save the null-space vectors
    vec_outfile* {.importc: "vec_outfile".}: array[256, char] ## * The precision with which to save the vectors
    save_prec* {.importc: "save_prec".}: QudaPrecision ## * Whether to inflate single-parity eigen-vector I/O to a full
                                                   ##         field (e.g., enabling this is required for compatability with
                                                   ##         MILC I/O)
    io_parity_inflate* {.importc: "io_parity_inflate".}: QudaBoolean ## * The Gflops rate of the eigensolver setup
    gflops* {.importc: "gflops".}: cdouble ## *< The time taken by the eigensolver setup
    secs* {.importc: "secs".}: cdouble ## * Which external library to use in the deflation operations (MAGMA or Eigen)
    extlib_type* {.importc: "extlib_type".}: QudaExtLibType ## -------------------------------------------------

  QudaMultigridParam* {.importc: "QudaMultigridParam", header: "quda.h", bycopy.} = object
    invert_param* {.importc: "invert_param".}: ptr QudaInvertParam
    eig_param* {.importc: "eig_param".}: array[QUDA_MAX_MG_LEVEL, ptr QudaEigParam] ## *
                                                                              ## Number of
                                                                              ## multigrid
                                                                              ## levels
    n_level* {.importc: "n_level".}: cint ## * Geometric block sizes to use on each level
    geo_block_size* {.importc: "geo_block_size".}: array[QUDA_MAX_MG_LEVEL,
        array[QUDA_MAX_DIM, cint]] ## * Spin block sizes to use on each level
    spin_block_size* {.importc: "spin_block_size".}: array[QUDA_MAX_MG_LEVEL, cint] ## *
                                                                               ## Number
                                                                               ## of
                                                                               ## null-space
                                                                               ## vectors
                                                                               ## to
                                                                               ## use
                                                                               ## on
                                                                               ## each
                                                                               ## level
    n_vec* {.importc: "n_vec".}: array[QUDA_MAX_MG_LEVEL, cint] ## * Precision to store the null-space vectors in (post block orthogonalization)
    precision_null* {.importc: "precision_null".}: array[QUDA_MAX_MG_LEVEL,
        QudaPrecision]        ## * Number of times to repeat Gram-Schmidt in block orthogonalization
    n_block_ortho* {.importc: "n_block_ortho".}: array[QUDA_MAX_MG_LEVEL, cint] ## *
                                                                           ## Verbosity on each level of the
                                                                           ## multigrid
    verbosity* {.importc: "verbosity".}: array[QUDA_MAX_MG_LEVEL, QudaVerbosity] ## *
                                                                            ## Inverter to use in the
                                                                            ## setup
                                                                            ## phase
    setup_inv_type* {.importc: "setup_inv_type".}: array[QUDA_MAX_MG_LEVEL,
        QudaInverterType]     ## * Number of setup iterations
    num_setup_iter* {.importc: "num_setup_iter".}: array[QUDA_MAX_MG_LEVEL, cint] ## *
                                                                             ## Tolerance to use in the
                                                                             ## setup
                                                                             ## phase
    setup_tol* {.importc: "setup_tol".}: array[QUDA_MAX_MG_LEVEL, cdouble] ## * Maximum number of iterations for each setup solver
    setup_maxiter* {.importc: "setup_maxiter".}: array[QUDA_MAX_MG_LEVEL, cint] ## *
                                                                           ## Maximum
                                                                           ## number of
                                                                           ## iterations for
                                                                           ## refreshing the
                                                                           ## null-space
                                                                           ## vectors
    setup_maxiter_refresh* {.importc: "setup_maxiter_refresh".}: array[
        QUDA_MAX_MG_LEVEL, cint] ## * Basis to use for CA-CGN(E/R) setup
    setup_ca_basis* {.importc: "setup_ca_basis".}: array[QUDA_MAX_MG_LEVEL,
        QudaCABasis]          ## * Basis size for CACG setup
    setup_ca_basis_size* {.importc: "setup_ca_basis_size".}: array[
        QUDA_MAX_MG_LEVEL, cint] ## * Minimum eigenvalue for Chebyshev CA basis
    setup_ca_lambda_min* {.importc: "setup_ca_lambda_min".}: array[
        QUDA_MAX_MG_LEVEL, cdouble] ## * Maximum eigenvalue for Chebyshev CA basis
    setup_ca_lambda_max* {.importc: "setup_ca_lambda_max".}: array[
        QUDA_MAX_MG_LEVEL, cdouble] ## * Null-space type to use in the setup phase
    setup_type* {.importc: "setup_type".}: QudaSetupType ## * Pre orthonormalize vectors in the setup phase
    pre_orthonormalize* {.importc: "pre_orthonormalize".}: QudaBoolean ## * Post
                                                                   ## orthonormalize vectors in the setup phase
    post_orthonormalize* {.importc: "post_orthonormalize".}: QudaBoolean ## * The solver that wraps around the coarse grid correction and smoother
    coarse_solver* {.importc: "coarse_solver".}: array[QUDA_MAX_MG_LEVEL,
        QudaInverterType]     ## * Tolerance for the solver that wraps around the coarse grid correction and smoother
    coarse_solver_tol* {.importc: "coarse_solver_tol".}: array[QUDA_MAX_MG_LEVEL,
        cdouble]              ## * Maximum number of iterations for the solver that wraps around the coarse grid correction and smoother
    coarse_solver_maxiter* {.importc: "coarse_solver_maxiter".}: array[
        QUDA_MAX_MG_LEVEL, cint] ## * Basis to use for CA-CGN(E/R) coarse solver
    coarse_solver_ca_basis* {.importc: "coarse_solver_ca_basis".}: array[
        QUDA_MAX_MG_LEVEL, QudaCABasis] ## * Basis size for CACG coarse solver
    coarse_solver_ca_basis_size* {.importc: "coarse_solver_ca_basis_size".}: array[
        QUDA_MAX_MG_LEVEL, cint] ## * Minimum eigenvalue for Chebyshev CA basis
    coarse_solver_ca_lambda_min* {.importc: "coarse_solver_ca_lambda_min".}: array[
        QUDA_MAX_MG_LEVEL, cdouble] ## * Maximum eigenvalue for Chebyshev CA basis
    coarse_solver_ca_lambda_max* {.importc: "coarse_solver_ca_lambda_max".}: array[
        QUDA_MAX_MG_LEVEL, cdouble] ## * Smoother to use on each level
    smoother* {.importc: "smoother".}: array[QUDA_MAX_MG_LEVEL, QudaInverterType] ## *
                                                                             ## Tolerance to use for the
                                                                             ## smoother /
                                                                             ## solver on
                                                                             ## each
                                                                             ## level
    smoother_tol* {.importc: "smoother_tol".}: array[QUDA_MAX_MG_LEVEL, cdouble] ## *
                                                                            ## Number of
                                                                            ## pre-smoother
                                                                            ## applications on each
                                                                            ## level
    nu_pre* {.importc: "nu_pre".}: array[QUDA_MAX_MG_LEVEL, cint] ## * Number of post-smoother applications on each level
    nu_post* {.importc: "nu_post".}: array[QUDA_MAX_MG_LEVEL, cint] ## * Over/under relaxation factor for the smoother at each level
    omega* {.importc: "omega".}: array[QUDA_MAX_MG_LEVEL, cdouble] ## * Precision to use for halo communication in the smoother
    smoother_halo_precision* {.importc: "smoother_halo_precision".}: array[
        QUDA_MAX_MG_LEVEL, QudaPrecision] ## * Whether to use additive or multiplicative Schwarz preconditioning in the smoother
    smoother_schwarz_type* {.importc: "smoother_schwarz_type".}: array[
        QUDA_MAX_MG_LEVEL, QudaSchwarzType] ## * Number of Schwarz cycles to apply
    smoother_schwarz_cycle* {.importc: "smoother_schwarz_cycle".}: array[
        QUDA_MAX_MG_LEVEL, cint] ## * The type of residual to send to the next coarse grid, and thus the
                               ## 	type of solution to receive back from this coarse grid
    coarse_grid_solution_type* {.importc: "coarse_grid_solution_type".}: array[
        QUDA_MAX_MG_LEVEL, QudaSolutionType] ## * The type of smoother solve to do on each grid (e/o preconditioning or not)
    smoother_solve_type* {.importc: "smoother_solve_type".}: array[
        QUDA_MAX_MG_LEVEL, QudaSolveType] ## * The type of multigrid cycle to perform at each level
    cycle_type* {.importc: "cycle_type".}: array[QUDA_MAX_MG_LEVEL,
        QudaMultigridCycleType] ## * Whether to use global reductions or not for the smoother / solver at each level
    global_reduction* {.importc: "global_reduction".}: array[QUDA_MAX_MG_LEVEL,
        QudaBoolean]          ## * Location where each level should be done
    location* {.importc: "location".}: array[QUDA_MAX_MG_LEVEL, QudaFieldLocation] ## *
                                                                              ## Location
                                                                              ## where
                                                                              ## the
                                                                              ## coarse-operator
                                                                              ## construction
                                                                              ## will be
                                                                              ## computedn
    setup_location* {.importc: "setup_location".}: array[QUDA_MAX_MG_LEVEL,
        QudaFieldLocation]    ## * Whether to use eigenvectors for the nullspace or, if the coarsest instance deflate
    use_eig_solver* {.importc: "use_eig_solver".}: array[QUDA_MAX_MG_LEVEL,
        QudaBoolean] ## * Minimize device memory allocations during the adaptive setup,
                    ##         placing temporary fields in mapped memory instad of device
                    ##         memory
    setup_minimize_memory* {.importc: "setup_minimize_memory".}: QudaBoolean ## * Whether to compute the null vectors or reload them
    compute_null_vector* {.importc: "compute_null_vector".}: QudaComputeNullVector ## *
                                                                               ## Whether
                                                                               ## to
                                                                               ## generate
                                                                               ## on
                                                                               ## all
                                                                               ## levels
                                                                               ## or
                                                                               ## just
                                                                               ## on
                                                                               ## level 0
    generate_all_levels* {.importc: "generate_all_levels".}: QudaBoolean ## * Whether to run the
                                                                     ## verification checks once set up is complete
    run_verify* {.importc: "run_verify".}: QudaBoolean ## * Whether to run null Vs eigen vector overlap checks once set up is complete
    run_low_mode_check* {.importc: "run_low_mode_check".}: QudaBoolean ## * Whether to run null vector oblique checks once set up is complete
    run_oblique_proj_check* {.importc: "run_oblique_proj_check".}: QudaBoolean ## *
                                                                           ## Whether to load the
                                                                           ## null-space
                                                                           ## vectors to disk
                                                                           ## (requires QIO)
    vec_load* {.importc: "vec_load".}: array[QUDA_MAX_MG_LEVEL, QudaBoolean] ## * Filename prefix where to load the
                                                                        ## null-space vectors
    vec_infile* {.importc: "vec_infile".}: array[QUDA_MAX_MG_LEVEL, array[256, char]] ##
                                                                                ## *
                                                                                ## Whether
                                                                                ## to
                                                                                ## store
                                                                                ## the
                                                                                ## null-space
                                                                                ## vectors
                                                                                ## to
                                                                                ## disk
                                                                                ## (requires
                                                                                ## QIO)
    vec_store* {.importc: "vec_store".}: array[QUDA_MAX_MG_LEVEL, QudaBoolean] ## *
                                                                          ## Filename prefix for where to save the
                                                                          ## null-space
                                                                          ## vectors
    vec_outfile* {.importc: "vec_outfile".}: array[QUDA_MAX_MG_LEVEL,
        array[256, char]]      ## * Whether to use and initial guess during coarse grid deflation
    coarse_guess* {.importc: "coarse_guess".}: QudaBoolean ## * Whether to preserve the deflation space during MG update
    preserve_deflation* {.importc: "preserve_deflation".}: QudaBoolean ## * The Gflops rate of the multigrid solver setup
    gflops* {.importc: "gflops".}: cdouble ## *< The time taken by the multigrid solver setup
    secs* {.importc: "secs".}: cdouble ## * Multiplicative factor for the mu parameter
    mu_factor* {.importc: "mu_factor".}: array[QUDA_MAX_MG_LEVEL, cdouble] ## * Boolean for
                                                                      ## aggregation type, implies staggered or not
    transfer_type* {.importc: "transfer_type".}: array[QUDA_MAX_MG_LEVEL,
        QudaTransferType]     ## * Whether to use tensor cores (if available)
    use_mma* {.importc: "use_mma".}: QudaBoolean ## * Whether to do a full (false) or thin (true) update in the context of updateMultigridQuda
    thin_update_only* {.importc: "thin_update_only".}: QudaBoolean

  QudaGaugeObservableParam* {.importc: "QudaGaugeObservableParam",
                             header: "quda.h", bycopy.} = object
    su_project* {.importc: "su_project".}: QudaBoolean ## *< Whether to porject onto the manifold prior to measurement
    compute_plaquette* {.importc: "compute_plaquette".}: QudaBoolean ## *< Whether to compute the plaquette
    plaquette* {.importc: "plaquette".}: array[3, cdouble] ## *< Total, spatial and temporal field energies, respectively
    compute_qcharge* {.importc: "compute_qcharge".}: QudaBoolean ## *< Whether to compute the topological charge and field energy
    qcharge* {.importc: "qcharge".}: cdouble ## *< Computed topological charge
    energy* {.importc: "energy".}: array[3, cdouble] ## *< Total, spatial and temporal field energies, respectively
    compute_qcharge_density* {.importc: "compute_qcharge_density".}: QudaBoolean ## *< Whether to compute the topological charge density
    qcharge_density* {.importc: "qcharge_density".}: pointer ## *< Pointer to host array of length volume where the q-charge density will be copied

  QudaBLASParam* {.importc: "QudaBLASParam", header: "quda.h", bycopy.} = object
    trans_a* {.importc: "trans_a".}: QudaBLASOperation ## *< operation op(A) that is non- or (conj.) transpose.
    trans_b* {.importc: "trans_b".}: QudaBLASOperation ## *< operation op(B) that is non- or (conj.) transpose.
    m* {.importc: "m".}: cint    ## *< number of rows of matrix op(A) and C.
    n* {.importc: "n".}: cint    ## *< number of columns of matrix op(B) and C.
    k* {.importc: "k".}: cint    ## *< number of columns of op(A) and rows of op(B).
    lda* {.importc: "lda".}: cint ## *< leading dimension of two-dimensional array used to store the matrix A.
    ldb* {.importc: "ldb".}: cint ## *< leading dimension of two-dimensional array used to store matrix B.
    ldc* {.importc: "ldc".}: cint ## *< leading dimension of two-dimensional array used to store matrix C.
    a_offset* {.importc: "a_offset".}: cint ## *< position of the A array from which begin read/write.
    b_offset* {.importc: "b_offset".}: cint ## *< position of the B array from which begin read/write.
    c_offset* {.importc: "c_offset".}: cint ## *< position of the C array from which begin read/write.
    a_stride* {.importc: "a_stride".}: cint ## *< stride of the A array in strided(batched) mode
    b_stride* {.importc: "b_stride".}: cint ## *< stride of the B array in strided(batched) mode
    c_stride* {.importc: "c_stride".}: cint ## *< stride of the C array in strided(batched) mode
    alpha* {.importc: "alpha".}: double_complex ## *< scalar used for multiplication.
    beta* {.importc: "beta".}: double_complex ## *< scalar used for multiplication. If beta==0, C does not have to be a valid input.
    batch_count* {.importc: "batch_count".}: cint ## *< number of pointers contained in arrayA, arrayB and arrayC.
    data_type* {.importc: "data_type".}: QudaBLASDataType ## *< Specifies if using S(C) or D(Z) BLAS type
    data_order* {.importc: "data_order".}: QudaBLASDataOrder ## *< Specifies if using Row or Column major


##
##  Interface functions, found in interface_quda.cpp
##
## *
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
## *
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

## *
##  @param mycomm User provided MPI communicator in place of MPI_COMM_WORLD
##

proc qudaSetCommHandle*(mycomm: pointer) {.importc: "qudaSetCommHandle",
                                        header: "quda.h".}
## *
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

proc initCommsGridQuda*(nDim: cint; dims: ptr ConstInt; `func`: QudaCommsMap;
                       fdata: pointer) {.importc: "initCommsGridQuda",
                                       header: "quda.h".}
## *
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
## *
##  Initialize the library persistant memory allocations (both host
##  and device).  This is a low-level interface that is called by
##  initQuda.  Calling initQudaMemory requires that the user has
##  previously called initQudaDevice.
##

proc initQudaMemory*() {.importc: "initQudaMemory", header: "quda.h".}
## *
##  Initialize the library.  This function is actually a wrapper
##  around calls to initQudaDevice() and initQudaMemory().
##
##  @param device  CUDA device number to use.  In a multi-GPU build,
##                 this parameter may either be set explicitly on a
##                 per-process basis or set to -1 to enable a default
##                 allocation of devices to processes.
##

proc initQuda*(device: cint) {.importc: "initQuda", header: "quda.h".}
## *
##  Finalize the library.
##

proc endQuda*() {.importc: "endQuda", header: "quda.h".}
## *
##  @brief update the radius for halos.
##  @details This should only be needed for automated testing when
##  different partitioning is applied within a single run.
##

proc updateR*() {.importc: "updateR", header: "quda.h".}
## *
##  A new QudaGaugeParam should always be initialized immediately
##  after it's defined (and prior to explicitly setting its members)
##  using this function.  Typical usage is as follows:
##
##    QudaGaugeParam gauge_param = newQudaGaugeParam();
##

proc newQudaGaugeParam*(): QudaGaugeParam {.importc: "newQudaGaugeParam",
    header: "quda.h".}
## *
##  A new QudaInvertParam should always be initialized immediately
##  after it's defined (and prior to explicitly setting its members)
##  using this function.  Typical usage is as follows:
##
##    QudaInvertParam invert_param = newQudaInvertParam();
##

proc newQudaInvertParam*(): QudaInvertParam {.importc: "newQudaInvertParam",
    header: "quda.h".}
## *
##  A new QudaMultigridParam should always be initialized immediately
##  after it's defined (and prior to explicitly setting its members)
##  using this function.  Typical usage is as follows:
##
##    QudaMultigridParam mg_param = newQudaMultigridParam();
##

proc newQudaMultigridParam*(): QudaMultigridParam {.
    importc: "newQudaMultigridParam", header: "quda.h".}
## *
##  A new QudaEigParam should always be initialized immediately
##  after it's defined (and prior to explicitly setting its members)
##  using this function.  Typical usage is as follows:
##
##    QudaEigParam eig_param = newQudaEigParam();
##

proc newQudaEigParam*(): QudaEigParam {.importc: "newQudaEigParam",
                                     header: "quda.h".}
## *
##  A new QudaGaugeObservableParam should always be initialized
##  immediately after it's defined (and prior to explicitly setting
##  its members) using this function.  Typical usage is as follows:
##
##    QudaGaugeParam obs_param = newQudaGaugeObservableParam();
##

proc newQudaGaugeObservableParam*(): QudaGaugeObservableParam {.
    importc: "newQudaGaugeObservableParam", header: "quda.h".}
## *
##  A new QudaBLASParam should always be initialized immediately
##  after it's defined (and prior to explicitly setting its members)
##  using this function.  Typical usage is as follows:
##
##    QudaBLASParam blas_param = newQudaBLASParam();
##

proc newQudaBLASParam*(): QudaBLASParam {.importc: "newQudaBLASParam",
                                       header: "quda.h".}
## *
##  Print the members of QudaGaugeParam.
##  @param param The QudaGaugeParam whose elements we are to print.
##

proc printQudaGaugeParam*(param: ptr QudaGaugeParam) {.
    importc: "printQudaGaugeParam", header: "quda.h".}
## *
##  Print the members of QudaInvertParam.
##  @param param The QudaInvertParam whose elements we are to print.
##

proc printQudaInvertParam*(param: ptr QudaInvertParam) {.
    importc: "printQudaInvertParam", header: "quda.h".}
## *
##  Print the members of QudaMultigridParam.
##  @param param The QudaMultigridParam whose elements we are to print.
##

proc printQudaMultigridParam*(param: ptr QudaMultigridParam) {.
    importc: "printQudaMultigridParam", header: "quda.h".}
## *
##  Print the members of QudaEigParam.
##  @param param The QudaEigParam whose elements we are to print.
##

proc printQudaEigParam*(param: ptr QudaEigParam) {.importc: "printQudaEigParam",
    header: "quda.h".}
## *
##  Print the members of QudaGaugeObservableParam.
##  @param param The QudaGaugeObservableParam whose elements we are to print.
##

proc printQudaGaugeObservableParam*(param: ptr QudaGaugeObservableParam) {.
    importc: "printQudaGaugeObservableParam", header: "quda.h".}
## *
##  Print the members of QudaBLASParam.
##  @param param The QudaBLASParam whose elements we are to print.
##

proc printQudaBLASParam*(param: ptr QudaBLASParam) {.importc: "printQudaBLASParam",
    header: "quda.h".}
## *
##  Load the gauge field from the host.
##  @param h_gauge Base pointer to host gauge field (regardless of dimensionality)
##  @param param   Contains all metadata regarding host and device storage
##

proc loadGaugeQuda*(h_gauge: pointer; param: ptr QudaGaugeParam) {.
    importc: "loadGaugeQuda", header: "quda.h".}
## *
##  Free QUDA's internal copy of the gauge field.
##

proc freeGaugeQuda*() {.importc: "freeGaugeQuda", header: "quda.h".}
## *
##  Save the gauge field to the host.
##  @param h_gauge Base pointer to host gauge field (regardless of dimensionality)
##  @param param   Contains all metadata regarding host and device storage
##

proc saveGaugeQuda*(h_gauge: pointer; param: ptr QudaGaugeParam) {.
    importc: "saveGaugeQuda", header: "quda.h".}
## *
##  Load the clover term and/or the clover inverse from the host.
##  Either h_clover or h_clovinv may be set to NULL.
##  @param h_clover    Base pointer to host clover field
##  @param h_cloverinv Base pointer to host clover inverse field
##  @param inv_param   Contains all metadata regarding host and device storage
##

proc loadCloverQuda*(h_clover: pointer; h_clovinv: pointer;
                    inv_param: ptr QudaInvertParam) {.importc: "loadCloverQuda",
    header: "quda.h".}
## *
##  Free QUDA's internal copy of the clover term and/or clover inverse.
##

proc freeCloverQuda*() {.importc: "freeCloverQuda", header: "quda.h".}
## *
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
## *
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
## *
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
## *
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
##  @param h_gauge     Base pointer to host gauge field (regardless of dimensionality)
##  @param gauge_param Contains all metadata regarding host and device storage for gauge field
##

proc invertMultiSrcQuda*(hp_x: ptr pointer; hp_b: ptr pointer;
                        param: ptr QudaInvertParam; h_gauge: pointer;
                        gauge_param: ptr QudaGaugeParam) {.
    importc: "invertMultiSrcQuda", header: "quda.h".}
## *
##  @brief Really the same with @invertMultiSrcQuda but for staggered-style fermions, by accepting pointers
##  to fat links and long links.
##  @param _hp_x       Array of solution spinor fields
##  @param _hp_b       Array of source spinor fields
##  @param param       Contains all metadata regarding host and device storage and solver parameters
##  @param milc_fatlinks     Base pointer to host **fat** gauge field (regardless of dimensionality)
##  @param milc_longlinks    Base pointer to host **long** gauge field (regardless of dimensionality)
##  @param gauge_param Contains all metadata regarding host and device storage for gauge field
##

proc invertMultiSrcStaggeredQuda*(hp_x: ptr pointer; hp_b: ptr pointer;
                                 param: ptr QudaInvertParam;
                                 milc_fatlinks: pointer; milc_longlinks: pointer;
                                 gauge_param: ptr QudaGaugeParam) {.
    importc: "invertMultiSrcStaggeredQuda", header: "quda.h".}
## *
##  @brief Really the same with @invertMultiSrcQuda but for clover-style fermions, by accepting pointers
##  to direct and inverse clover field pointers.
##  @param _hp_x       Array of solution spinor fields
##  @param _hp_b       Array of source spinor fields
##  @param param       Contains all metadata regarding host and device storage and solver parameters
##  @param h_gauge     Base pointer to host gauge field (regardless of dimensionality)
##  @param gauge_param Contains all metadata regarding host and device storage for gauge field
##  @param h_clover    Base pointer to the direct clover field
##  @param h_clovinv   Base pointer to the inverse clover field
##

proc invertMultiSrcCloverQuda*(hp_x: ptr pointer; hp_b: ptr pointer;
                              param: ptr QudaInvertParam; h_gauge: pointer;
                              gauge_param: ptr QudaGaugeParam; h_clover: pointer;
                              h_clovinv: pointer) {.
    importc: "invertMultiSrcCloverQuda", header: "quda.h".}
## *
##  Solve for multiple shifts (e.g., masses).
##  @param _hp_x    Array of solution spinor fields
##  @param _hp_b    Source spinor fields
##  @param param  Contains all metadata regarding host and device
##                storage and solver parameters
##

proc invertMultiShiftQuda*(hp_x: ptr pointer; hp_b: pointer;
                          param: ptr QudaInvertParam) {.
    importc: "invertMultiShiftQuda", header: "quda.h".}
## *
##  Setup the multigrid solver, according to the parameters set in param.  It
##  is assumed that the gauge field has already been loaded via
##  loadGaugeQuda().
##  @param param  Contains all metadata regarding host and device
##                storage and solver parameters
##

proc newMultigridQuda*(param: ptr QudaMultigridParam): pointer {.
    importc: "newMultigridQuda", header: "quda.h".}
## *
##  @brief Free resources allocated by the multigrid solver
##  @param mg_instance Pointer to instance of multigrid_solver
##  @param param Contains all metadata regarding host and device
##  storage and solver parameters
##

proc destroyMultigridQuda*(mg_instance: pointer) {.importc: "destroyMultigridQuda",
    header: "quda.h".}
## *
##  @brief Updates the multigrid preconditioner for the new gauge / clover field
##  @param mg_instance Pointer to instance of multigrid_solver
##  @param param Contains all metadata regarding host and device
##  storage and solver parameters, of note contains a flag specifying whether
##  to do a full update or a thin update.
##

proc updateMultigridQuda*(mg_instance: pointer; param: ptr QudaMultigridParam) {.
    importc: "updateMultigridQuda", header: "quda.h".}
## *
##  @brief Dump the null-space vectors to disk
##  @param[in] mg_instance Pointer to the instance of multigrid_solver
##  @param[in] param Contains all metadata regarding host and device
##  storage and solver parameters (QudaMultigridParam::vec_outfile
##  sets the output filename prefix).
##

proc dumpMultigridQuda*(mg_instance: pointer; param: ptr QudaMultigridParam) {.
    importc: "dumpMultigridQuda", header: "quda.h".}
## *
##  Apply the Dslash operator (D_{eo} or D_{oe}).
##  @param h_out  Result spinor field
##  @param h_in   Input spinor field
##  @param param  Contains all metadata regarding host and device
##                storage
##  @param parity The destination parity of the field
##

proc dslashQuda*(h_out: pointer; h_in: pointer; inv_param: ptr QudaInvertParam;
                parity: QudaParity) {.importc: "dslashQuda", header: "quda.h".}
## *
##  @brief Perform the solve like @dslashQuda but for multiple rhs by spliting the comm grid into
##  sub-partitions: each sub-partition does one or more rhs'.
##  The QudaInvertParam object specifies how the solve should be performed on each sub-partition.
##  Unlike @invertQuda, the interface also takes the host side gauge as
##  input - gauge field is not required to be loaded beforehand.
##  @param _hp_x       Array of solution spinor fields
##  @param _hp_b       Array of source spinor fields
##  @param param       Contains all metadata regarding host and device storage and solver parameters
##  @param parity      Parity to apply dslash on
##  @param h_gauge     Base pointer to host gauge field (regardless of dimensionality)
##  @param gauge_param Contains all metadata regarding host and device storage for gauge field
##

proc dslashMultiSrcQuda*(hp_x: ptr pointer; hp_b: ptr pointer;
                        param: ptr QudaInvertParam; parity: QudaParity;
                        h_gauge: pointer; gauge_param: ptr QudaGaugeParam) {.
    importc: "dslashMultiSrcQuda", header: "quda.h".}
## *
##  @brief Really the same with @dslashMultiSrcQuda but for staggered-style fermions, by accepting pointers
##  to fat links and long links.
##  @param _hp_x       Array of solution spinor fields
##  @param _hp_b       Array of source spinor fields
##  @param param       Contains all metadata regarding host and device storage and solver parameters
##  @param parity      Parity to apply dslash on
##  @param milc_fatlinks     Base pointer to host **fat** gauge field (regardless of dimensionality)
##  @param milc_longlinks    Base pointer to host **long** gauge field (regardless of dimensionality)
##  @param gauge_param Contains all metadata regarding host and device storage for gauge field
##

proc dslashMultiSrcStaggeredQuda*(hp_x: ptr pointer; hp_b: ptr pointer;
                                 param: ptr QudaInvertParam; parity: QudaParity;
                                 milc_fatlinks: pointer; milc_longlinks: pointer;
                                 gauge_param: ptr QudaGaugeParam) {.
    importc: "dslashMultiSrcStaggeredQuda", header: "quda.h".}
## *
##  @brief Really the same with @dslashMultiSrcQuda but for clover-style fermions, by accepting pointers
##  to direct and inverse clover field pointers.
##  @param _hp_x       Array of solution spinor fields
##  @param _hp_b       Array of source spinor fields
##  @param param       Contains all metadata regarding host and device storage and solver parameters
##  @param parity      Parity to apply dslash on
##  @param h_gauge     Base pointer to host gauge field (regardless of dimensionality)
##  @param gauge_param Contains all metadata regarding host and device storage for gauge field
##  @param h_clover    Base pointer to the direct clover field
##  @param h_clovinv   Base pointer to the inverse clover field
##

proc dslashMultiSrcCloverQuda*(hp_x: ptr pointer; hp_b: ptr pointer;
                              param: ptr QudaInvertParam; parity: QudaParity;
                              h_gauge: pointer; gauge_param: ptr QudaGaugeParam;
                              h_clover: pointer; h_clovinv: pointer) {.
    importc: "dslashMultiSrcCloverQuda", header: "quda.h".}
## *
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
## *
##  Apply the full Dslash matrix, possibly even/odd preconditioned.
##  @param h_out  Result spinor field
##  @param h_in   Input spinor field
##  @param param  Contains all metadata regarding host and device
##                storage
##

proc MatQuda*(h_out: pointer; h_in: pointer; inv_param: ptr QudaInvertParam) {.
    importc: "MatQuda", header: "quda.h".}
## *
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
## *
##  Either downloads and sets the resident momentum field, or uploads
##  and returns the resident momentum field
##
##  @param[in,out] mom The external momentum field
##  @param[in] param The parameters of the external field
##

proc momResidentQuda*(mom: pointer; param: ptr QudaGaugeParam) {.
    importc: "momResidentQuda", header: "quda.h".}
## *
##  Compute the gauge force and update the mometum field
##
##  @param mom The momentum field to be updated
##  @param sitelink The gauge field from which we compute the force
##  @param input_path_buf[dim][num_paths][path_length]
##  @param path_length One less that the number of links in a loop (e.g., 3 for a staple)
##  @param loop_coeff Coefficients of the different loops in the Symanzik action
##  @param num_paths How many contributions from path_length different "staples"
##  @param max_length The maximum number of non-zero of links in any path in the action
##  @param dt The integration step size (for MILC this is dt*beta/3)
##  @param param The parameters of the external fields and the computation settings
##

proc computeGaugeForceQuda*(mom: pointer; sitelink: pointer;
                           input_path_buf: ptr ptr ptr cint; path_length: ptr cint;
                           loop_coeff: ptr cdouble; num_paths: cint;
                           max_length: cint; dt: cdouble;
                           qudaGaugeParam: ptr QudaGaugeParam): cint {.
    importc: "computeGaugeForceQuda", header: "quda.h".}
## *
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
## *
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
## *
##  Project the input field on the SU(3) group.  If the target
##  tolerance is not met, this routine will give a runtime error.
##
##  @param gauge_h The gauge field to be updated
##  @param tol The tolerance to which we iterate
##  @param param The parameters of the gauge field
##

proc projectSU3Quda*(gauge_h: pointer; tol: cdouble; param: ptr QudaGaugeParam) {.
    importc: "projectSU3Quda", header: "quda.h".}
## *
##  Evaluate the momentum contribution to the Hybrid Monte Carlo
##  action.
##
##  @param momentum The momentum field
##  @param param The parameters of the external fields and the computation settings
##  @return momentum action
##

proc momActionQuda*(momentum: pointer; param: ptr QudaGaugeParam): cdouble {.
    importc: "momActionQuda", header: "quda.h".}
## *
##  Allocate a gauge (matrix) field on the device and optionally download a host gauge field.
##
##  @param gauge The host gauge field (optional - if set to 0 then the gauge field zeroed)
##  @param geometry The geometry of the matrix field to create (1 - scalar, 4 - vector, 6 - tensor)
##  @param param The parameters of the external field and the field to be created
##  @return Pointer to the gauge field (cast as a void*)
##

proc createGaugeFieldQuda*(gauge: pointer; geometry: cint; param: ptr QudaGaugeParam): pointer {.
    importc: "createGaugeFieldQuda", header: "quda.h".}
## *
##  Copy the QUDA gauge (matrix) field on the device to the CPU
##
##  @param outGauge Pointer to the host gauge field
##  @param inGauge Pointer to the device gauge field (QUDA device field)
##  @param param The parameters of the host and device fields
##

proc saveGaugeFieldQuda*(outGauge: pointer; inGauge: pointer;
                        param: ptr QudaGaugeParam) {.importc: "saveGaugeFieldQuda",
    header: "quda.h".}
## *
##  Reinterpret gauge as a pointer to cudaGaugeField and call destructor.
##
##  @param gauge Gauge field to be freed
##

proc destroyGaugeFieldQuda*(gauge: pointer) {.importc: "destroyGaugeFieldQuda",
    header: "quda.h".}
## *
##  Compute the clover field and its inverse from the resident gauge field.
##
##  @param param The parameters of the clover field to create
##

proc createCloverQuda*(param: ptr QudaInvertParam) {.importc: "createCloverQuda",
    header: "quda.h".}
## *
##  Compute the clover force contributions in each dimension mu given
##  the array of solution fields, and compute the resulting momentum
##  field.
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
## *
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
## *
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
## *
##      @brief Generate Gaussian distributed fields and store in the
##      resident gauge field.  We create a Gaussian-distributed su(n)
##      field and exponentiate it, e.g., U = exp(sigma * H), where H is
##      the distributed su(n) field and beta is the width of the
##      distribution (beta = 0 results in a free field, and sigma = 1 has
##      maximum disorder).
##
##      @param seed The seed used for the RNG
##      @param sigma Width of Gaussian distrubution
##

proc gaussGaugeQuda*(seed: culonglong; sigma: cdouble) {.importc: "gaussGaugeQuda",
    header: "quda.h".}
## *
##  Computes the total, spatial and temporal plaquette averages of the loaded gauge configuration.
##  @param Array for storing the averages (total, spatial, temporal)
##

proc plaqQuda*(plaq: array[3, cdouble]) {.importc: "plaqQuda", header: "quda.h".}
## *
##  Performs a deep copy from the internal extendedGaugeResident field.
##  @param Pointer to externally allocated GaugeField
##

proc copyExtendedResidentGaugeQuda*(resident_gauge: pointer) {.
    importc: "copyExtendedResidentGaugeQuda", header: "quda.h".}
## *
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
## *
##  Performs APE smearing on gaugePrecise and stores it in gaugeSmeared
##  @param n_steps Number of steps to apply.
##  @param alpha  Alpha coefficient for APE smearing.
##  @param meas_interval Measure the Q charge every Nth step
##

proc performAPEnStep*(n_steps: cuint; alpha: cdouble; meas_interval: cint) {.
    importc: "performAPEnStep", header: "quda.h".}
## *
##  Performs STOUT smearing on gaugePrecise and stores it in gaugeSmeared
##  @param n_steps Number of steps to apply.
##  @param rho    Rho coefficient for STOUT smearing.
##  @param meas_interval Measure the Q charge every Nth step
##

proc performSTOUTnStep*(n_steps: cuint; rho: cdouble; meas_interval: cint) {.
    importc: "performSTOUTnStep", header: "quda.h".}
## *
##  Performs Over Imroved STOUT smearing on gaugePrecise and stores it in gaugeSmeared
##  @param n_steps Number of steps to apply.
##  @param rho    Rho coefficient for STOUT smearing.
##  @param epsilon Epsilon coefficient for Over Improved STOUT smearing.
##  @param meas_interval Measure the Q charge every Nth step
##

proc performOvrImpSTOUTnStep*(n_steps: cuint; rho: cdouble; epsilon: cdouble;
                             meas_interval: cint) {.
    importc: "performOvrImpSTOUTnStep", header: "quda.h".}
## *
##  Performs Wilson Flow on gaugePrecise and stores it in gaugeSmeared
##  @param n_steps Number of steps to apply.
##  @param step_size Size of Wilson Flow step
##  @param meas_interval Measure the Q charge and field energy every Nth step
##  @param wflow_type 1x1 Wilson or 2x1 Symanzik flow type
##

proc performWFlownStep*(n_steps: cuint; step_size: cdouble; meas_interval: cint;
                       wflow_type: QudaWFlowType) {.importc: "performWFlownStep",
    header: "quda.h".}
## *
##  @brief Calculates a variety of gauge-field observables.  If a
##  smeared gauge field is presently loaded (in gaugeSmeared) the
##  observables are computed on this, else the resident gauge field
##  will be used.
##  @param[in,out] param Parameter struct that defines which
##  observables we are making and the resulting observables.
##

proc gaugeObservablesQuda*(param: ptr QudaGaugeObservableParam) {.
    importc: "gaugeObservablesQuda", header: "quda.h".}
## *
##  Public function to perform color contractions of the host spinors x and y.
##  @param[in] x pointer to host data
##  @param[in] y pointer to host data
##  @param[out] result pointer to the 16 spin projections per lattice site
##  @param[in] cType Which type of contraction (open, degrand-rossi, etc)
##  @param[in] param meta data for construction of ColorSpinorFields.
##  @param[in] X spacetime data for construction of ColorSpinorFields.
##

proc contractQuda*(x: pointer; y: pointer; result: pointer; cType: QudaContractType;
                  param: ptr QudaInvertParam; X: ptr ConstInt) {.
    importc: "contractQuda", header: "quda.h".}
## *
##  @brief Gauge fixing with overrelaxation with support for single and multi GPU.
##  @param[in,out] gauge, gauge field to be fixed
##  @param[in] gauge_dir, 3 for Coulomb gauge fixing, other for Landau gauge fixing
##  @param[in] Nsteps, maximum number of steps to perform gauge fixing
##  @param[in] verbose_interval, print gauge fixing info when iteration count is a multiple of this
##  @param[in] relax_boost, gauge fixing parameter of the overrelaxation method, most common value is 1.5 or 1.7.
##  @param[in] tolerance, torelance value to stop the method, if this value is zero then the method stops when iteration reachs the maximum number of steps defined by Nsteps
##  @param[in] reunit_interval, reunitarize gauge field when iteration count is a multiple of this
##  @param[in] stopWtheta, 0 for MILC criterion and 1 to use the theta value
##  @param[in] param The parameters of the external fields and the computation settings
##  @param[out] timeinfo
##

proc computeGaugeFixingOVRQuda*(gauge: pointer; gauge_dir: cuint; Nsteps: cuint;
                               verbose_interval: cuint; relax_boost: cdouble;
                               tolerance: cdouble; reunit_interval: cuint;
                               stopWtheta: cuint; param: ptr QudaGaugeParam;
                               timeinfo: ptr cdouble): cint {.
    importc: "computeGaugeFixingOVRQuda", header: "quda.h".}
## *
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
##  @param[out] timeinfo
##

proc computeGaugeFixingFFTQuda*(gauge: pointer; gauge_dir: cuint; Nsteps: cuint;
                               verbose_interval: cuint; alpha: cdouble;
                               autotune: cuint; tolerance: cdouble;
                               stopWtheta: cuint; param: ptr QudaGaugeParam;
                               timeinfo: ptr cdouble): cint {.
    importc: "computeGaugeFixingFFTQuda", header: "quda.h".}
## *
##  @brief Strided Batched GEMM
##  @param[in] arrayA The array containing the A matrix data
##  @param[in] arrayB The array containing the A matrix data
##  @param[in] arrayC The array containing the A matrix data
##  @param[in] native boolean to use either the native or generic version
##  @param[in] param The data defining the problem execution.
##

proc blasGEMMQuda*(arrayA: pointer; arrayB: pointer; arrayC: pointer;
                  native: QudaBoolean; param: ptr QudaBLASParam) {.
    importc: "blasGEMMQuda", header: "quda.h".}
## *
##  @brief Flush the chronological history for the given index
##  @param[in] index Index for which we are flushing
##

proc flushChronoQuda*(index: cint) {.importc: "flushChronoQuda", header: "quda.h".}
## *
##  Open/Close MAGMA library
##
##

proc openMagma*() {.importc: "openMagma", header: "quda.h".}
proc closeMagma*() {.importc: "closeMagma", header: "quda.h".}
## *
##  Create deflation solver resources.
##
##

proc newDeflationQuda*(param: ptr QudaEigParam): pointer {.
    importc: "newDeflationQuda", header: "quda.h".}
## *
##  Free resources allocated by the deflated solver
##

proc destroyDeflationQuda*(df_instance: pointer) {.importc: "destroyDeflationQuda",
    header: "quda.h".}
proc setMPICommHandleQuda*(mycomm: pointer) {.importc: "setMPICommHandleQuda",
    header: "quda.h".}
##  remove NVRTC WAR

##  #include <quda_new_interface.h>
