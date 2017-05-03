import enum_quda
import quda_constants

## *
##  @file  quda.h
##  @brief Main header file for the QUDA library
## 
##  Note to QUDA developers: When adding new members to QudaGaugeParam
##  and QudaInvertParam, be sure to update lib/check_params.h as well
##  as the Fortran interface in lib/quda_fortran.F90.
## 

## *
##  Parameters having to do with the gauge field or the
##  interpretation of the gauge field by various Dirac operators
## 

type
  QudaGaugeParam* {.importc: "QudaGaugeParam", header: "quda.h".} = object
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
    cuda_prec_precondition* {.importc: "cuda_prec_precondition".}: QudaPrecision ## *< The precision of the preconditioner gauge field
    reconstruct_precondition* {.importc: "reconstruct_precondition".}: QudaReconstructType ## *< The recontruction type of the preconditioner gauge field
    gauge_fix* {.importc: "gauge_fix".}: QudaGaugeFixed ## *< Whether the input gauge field is in the axial gauge or not
    ga_pad* {.importc: "ga_pad".}: cint ## *< The pad size that the cudaGaugeField will use (default=0)
    site_ga_pad* {.importc: "site_ga_pad".}: cint ## *< Used by link fattening and the gauge and fermion forces
    staple_pad* {.importc: "staple_pad".}: cint ## *< Used by link fattening
    llfat_ga_pad* {.importc: "llfat_ga_pad".}: cint ## *< Used by link fattening
    mom_ga_pad* {.importc: "mom_ga_pad".}: cint ## *< Used by the gauge and fermion forces
    gaugeGiB* {.importc: "gaugeGiB".}: cdouble ## *< The storage used by the gauge fields
    preserve_gauge* {.importc: "preserve_gauge".}: cint ## *< Used by link fattening
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
  

## *
##  Parameters relating to the solver and the choice of Dirac operator.
## 

type
  QudaInvertParam* {.importc: "QudaInvertParam", header: "quda.h".} = object
    input_location* {.importc: "input_location".}: QudaFieldLocation ## *< The location of the input field
    output_location* {.importc: "output_location".}: QudaFieldLocation ## *< The location of the output field
    dslash_type* {.importc: "dslash_type".}: QudaDslashType ## *< The Dirac Dslash type that is being used
    inv_type* {.importc: "inv_type".}: QudaInverterType ## *< Which linear solver to use
    mass* {.importc: "mass".}: cdouble ## *< Used for staggered only
    kappa* {.importc: "kappa".}: cdouble ## *< Used for Wilson and Wilson-clover
    m5* {.importc: "m5".}: cdouble ## *< Domain wall height
    Ls* {.importc: "Ls".}: cint  ## *< Extent of the 5th dimension (for domain wall)
    b_5* {.importc: "b_5".}: array[QUDA_MAX_DWF_LS, cdouble] ## *< MDWF coefficients
    c_5* {.importc: "c_5".}: array[QUDA_MAX_DWF_LS, cdouble] ## *< will be used only for the mobius type of Fermion
    mu* {.importc: "mu".}: cdouble ## *< Twisted mass parameter
    epsilon* {.importc: "epsilon".}: cdouble ## *< Twisted mass parameter
    twist_flavor* {.importc: "twist_flavor".}: QudaTwistFlavorType ## *< Twisted mass flavor
    tol* {.importc: "tol".}: cdouble ## *< Solver tolerance in the L2 residual norm
    tol_restart* {.importc: "tol_restart".}: cdouble ## *< Solver tolerance in the L2 residual norm (used to restart InitCG)
    tol_hq* {.importc: "tol_hq".}: cdouble ## *< Solver tolerance in the heavy quark residual norm
    true_res* {.importc: "true_res".}: cdouble ## *< Actual L2 residual norm achieved in solver
    true_res_hq* {.importc: "true_res_hq".}: cdouble ## *< Actual heavy quark residual norm achieved in solver
    maxiter* {.importc: "maxiter".}: cint ## *< Maximum number of iterations in the linear solver
    reliable_delta* {.importc: "reliable_delta".}: cdouble ## *< Reliable update tolerance
    use_sloppy_partial_accumulator* {.importc: "use_sloppy_partial_accumulator".}: cint ## *< Whether to keep the partial solution accumuator in sloppy precision
                                                                                    ## *< This parameter determines how many consective reliable update
                                                                                    ##     residual increases we tolerate before terminating the solver,
                                                                                    ##     i.e., how long do we want to keep trying to converge
    max_res_increase* {.importc: "max_res_increase".}: cint ## *< This parameter determines how many total reliable update
                                                        ##     residual increases we tolerate before terminating the solver,
                                                        ##     i.e., how long do we want to keep trying to converge
    max_res_increase_total* {.importc: "max_res_increase_total".}: cint ## *< After how many iterations shall the heavy quark residual be updated
    heavy_quark_check* {.importc: "heavy_quark_check".}: cint
    pipeline* {.importc: "pipeline".}: cint ## *< Whether to use a pipelined solver with less global sums
    num_offset* {.importc: "num_offset".}: cint ## *< Number of offsets in the multi-shift solver
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
        QUDA_MAX_MULTI_SHIFT, cdouble]
    solution_type* {.importc: "solution_type".}: QudaSolutionType ## *< Type of system to solve
    solve_type* {.importc: "solve_type".}: QudaSolveType ## *< How to solve it
    matpc_type* {.importc: "matpc_type".}: QudaMatPCType ## *< The preconditioned matrix type
    dagger* {.importc: "dagger".}: QudaDagType ## *< Whether we are using the Hermitian conjugate system or not
    mass_normalization* {.importc: "mass_normalization".}: QudaMassNormalization ## *< The mass normalization is being used by the caller
    solver_normalization* {.importc: "solver_normalization".}: QudaSolverNormalization ## *< The normalization desired in the solver
    preserve_source* {.importc: "preserve_source".}: QudaPreserveSource ## *< Preserve the source or not in the linear solver (deprecated)
    cpu_prec* {.importc: "cpu_prec".}: QudaPrecision ## *< The precision used by the input fermion fields
    cuda_prec* {.importc: "cuda_prec".}: QudaPrecision ## *< The precision used by the QUDA solver
    cuda_prec_sloppy* {.importc: "cuda_prec_sloppy".}: QudaPrecision ## *< The precision used by the QUDA sloppy operator
    cuda_prec_precondition* {.importc: "cuda_prec_precondition".}: QudaPrecision ## *< The precision used by the QUDA preconditioner
    dirac_order* {.importc: "dirac_order".}: QudaDiracFieldOrder ## *< The order of the input and output fermion fields
    gamma_basis* {.importc: "gamma_basis".}: QudaGammaBasis ## *< Gamma basis of the input and output host fields
    clover_location* {.importc: "clover_location".}: QudaFieldLocation ## *< The location of the clover field
    clover_cpu_prec* {.importc: "clover_cpu_prec".}: QudaPrecision ## *< The precision used for the input clover field
    clover_cuda_prec* {.importc: "clover_cuda_prec".}: QudaPrecision ## *< The precision used for the clover field in the QUDA solver
    clover_cuda_prec_sloppy* {.importc: "clover_cuda_prec_sloppy".}: QudaPrecision ## *< The precision used for the clover field in the QUDA sloppy operator
    clover_cuda_prec_precondition* {.importc: "clover_cuda_prec_precondition".}: QudaPrecision ## *< The precision used for the clover field in the QUDA preconditioner
    clover_order* {.importc: "clover_order".}: QudaCloverFieldOrder ## *< The order of the input clover field
    use_init_guess* {.importc: "use_init_guess".}: QudaUseInitGuess ## *< Whether to use an initial guess in the solver or not
    clover_coeff* {.importc: "clover_coeff".}: cdouble ## *< Coefficient of the clover term
    compute_clover_trlog* {.importc: "compute_clover_trlog".}: cint ## *< Whether to compute the trace log of the clover term
    trlogA* {.importc: "trlogA".}: array[2, cdouble] ## *< The trace log of the clover term (even/odd computed separately)
    verbosity* {.importc: "verbosity".}: QudaVerbosity ## *< The verbosity setting to use in the solver
    sp_pad* {.importc: "sp_pad".}: cint ## *< The padding to use for the fermion fields
    cl_pad* {.importc: "cl_pad".}: cint ## *< The padding to use for the clover fields
    iter* {.importc: "iter".}: cint ## *< The number of iterations performed by the solver
    spinorGiB* {.importc: "spinorGiB".}: cdouble ## *< The memory footprint of the fermion fields
    cloverGiB* {.importc: "cloverGiB".}: cdouble ## *< The memory footprint of the clover fields
    gflops* {.importc: "gflops".}: cdouble ## *< The Gflops rate of the solver
    secs* {.importc: "secs".}: cdouble ## *< The time taken by the solver
    tune* {.importc: "tune".}: QudaTune ## *< Enable auto-tuning? (default = QUDA_TUNE_YES)
                                    ## * Number of steps in s-step algorithms
    Nsteps* {.importc: "Nsteps".}: cint ## * Maximum size of Krylov space used by solver
    gcrNkrylov* {.importc: "gcrNkrylov".}: cint ## 
                                            ##  The following parameters are related to the domain-decomposed
                                            ##  preconditioner, if enabled.
                                            ## 
                                            ## *
                                            ##  The inner Krylov solver used in the preconditioner.  Set to
                                            ##  QUDA_INVALID_INVERTER to disable the preconditioner entirely.
                                            ## 
    inv_type_precondition* {.importc: "inv_type_precondition".}: QudaInverterType ## *
                                                                              ##       Dirac Dslash used in preconditioner
                                                                              ## 
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
    omega* {.importc: "omega".}: cdouble ## * Number of preconditioner cycles to perform per iteration
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
    nev* {.importc: "nev".}: cint ## * EeigCG  : Search space dimension 
                              ##   gmresdr : Krylov subspace dimension  
                              ## 
    max_search_dim* {.importc: "max_search_dim".}: cint ## for magma library this parameter must be multiple 16?
                                                    ## * For systems with many RHS: current RHS index
    rhs_idx* {.importc: "rhs_idx".}: cint ## * Specifies deflation space volume: total number of eigenvectors is nev*deflation_grid
    deflation_grid* {.importc: "deflation_grid".}: cint ## * eigCG: specifies whether to use reduced eigenvector set
    use_reduced_vector_set* {.importc: "use_reduced_vector_set".}: cint ## * eigCG: selection criterion for the reduced eigenvector set
    eigenval_tol* {.importc: "eigenval_tol".}: cdouble ## * mixed precision eigCG tuning parameter:  whether to use cg refinement corrections in the incremental stage
    use_cg_updates* {.importc: "use_cg_updates".}: cint ## * mixed precision eigCG tuning parameter:  tolerance for cg refinement corrections in the incremental stage
    cg_iterref_tol* {.importc: "cg_iterref_tol".}: cdouble ## * mixed precision eigCG tuning parameter:  minimum search vector space restarts
    eigcg_max_restarts* {.importc: "eigcg_max_restarts".}: cint ## * initCG tuning parameter:  maximum restarts
    max_restart_num* {.importc: "max_restart_num".}: cint ## * initCG tuning parameter:  decrease in absolute value of the residual within each restart cycle
    inc_tol* {.importc: "inc_tol".}: cdouble ## * Whether to make the solution vector(s) after the solve
    make_resident_solution* {.importc: "make_resident_solution".}: cint ## * Whether to use the resident solution vector(s)
    use_resident_solution* {.importc: "use_resident_solution".}: cint


##  Parameter set for solving the eigenvalue problems.
##  Eigen problems are tightly related with Ritz algorithm.
##  And the Lanczos algorithm use the Ritz operator.
##  For Ritz matrix operation, 
##  we need to know about the solution type of dirac operator.
##  For acceleration, we are also using chevisov polynomial method.
##  And nk, np values are needed Implicit Restart Lanczos method
##  which is optimized form of Lanczos algorithm

type
  QudaEigParam* {.importc: "QudaEigParam", header: "quda.h".} = object
    invert_param* {.importc: "invert_param".}: ptr QudaInvertParam
    RitzMat_lanczos* {.importc: "RitzMat_lanczos".}: QudaSolutionType
    RitzMat_Convcheck* {.importc: "RitzMat_Convcheck".}: QudaSolutionType
    eig_type* {.importc: "eig_type".}: QudaEigType
    MatPoly_param* {.importc: "MatPoly_param".}: ptr cdouble
    NPoly* {.importc: "NPoly".}: cint
    Stp_residual* {.importc: "Stp_residual".}: cdouble
    nk* {.importc: "nk".}: cint
    np* {.importc: "np".}: cint
    f_size* {.importc: "f_size".}: cint
    eigen_shift* {.importc: "eigen_shift".}: cdouble


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
  QudaCommsMap* = proc (coords: ptr cint; fdata: pointer): cint

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

proc initCommsGridQuda*(nDim: cint; dims: ptr cint; `func`: QudaCommsMap; fdata: pointer) {.
    importc: "initCommsGridQuda", header: "quda.h".}
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
##  A new QudaEigParam should always be initialized immediately
##  after it's defined (and prior to explicitly setting its members)
##  using this function.  Typical usage is as follows:
## 
##    QudaEigParam eig_param = newQudaEigParam();
## 

proc newQudaEigParam*(): QudaEigParam {.importc: "newQudaEigParam", header: "quda.h".}
## *
##  Print the members of QudaGaugeParam.
##  @param param The QudaGaugeParam whose elements we are to print.
## 

proc printQudaGaugeParam*(param: ptr QudaGaugeParam) {.
    importc: "printQudaGaugeParam", header: "quda.h".}
## *
##  Print the members of QudaGaugeParam.
##  @param param The QudaGaugeParam whose elements we are to print.
## 

proc printQudaInvertParam*(param: ptr QudaInvertParam) {.
    importc: "printQudaInvertParam", header: "quda.h".}
## *
##  Print the members of QudaEigParam.
##  @param param The QudaEigParam whose elements we are to print.
## 

proc printQudaEigParam*(param: ptr QudaEigParam) {.importc: "printQudaEigParam",
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
##  Solve for multiple shifts (e.g., masses).
##  @param _hp_x    Array of solution spinor fields
##  @param _hp_b    Array of source spinor fields
##  @param param  Contains all metadata regarding host and device
##                storage and solver parameters
## 

proc invertMultiShiftQuda*(hp_x: ptr pointer; hp_b: pointer;
                          param: ptr QudaInvertParam) {.
    importc: "invertMultiShiftQuda", header: "quda.h".}
## *
##  Deflated solvers interface (e.g., based on invremental deflation space constructors, like incremental eigCG).
##  @param _h_x    Outnput: array of solution spinor fields (typically O(10))
##  @param _h_b    Input: array of source spinor fields (typically O(10))
##  @param _h_u    Input/Output: array of Ritz spinor fields (typically O(100))
##  @param _h_h    Input/Output: complex projection mutirx (typically O(100))
##  @param param  Contains all metadata regarding host and device
##                storage and solver parameters
## 

proc incrementalEigQuda*(h_x: pointer; h_b: pointer; param: ptr QudaInvertParam;
                        h_u: pointer; inv_eigenvals: ptr cdouble) {.
    importc: "incrementalEigQuda", header: "quda.h".}
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
##  Apply the Dslash operator (D_{eo} or D_{oe}) for 4D EO preconditioned DWF.
##  @param h_out  Result spinor field
##  @param h_in   Input spinor field
##  @param param  Contains all metadata regarding host and device
##                storage
##  @param parity The destination parity of the field
##  @param test_type Choose a type of dslash operators 
## 

proc dslashQuda_4dpc*(h_out: pointer; h_in: pointer; inv_param: ptr QudaInvertParam;
                     parity: QudaParity; test_type: cint) {.
    importc: "dslashQuda_4dpc", header: "quda.h".}
## *
##  Apply the Dslash operator (D_{eo} or D_{oe}) for Mobius DWF.
##  @param h_out  Result spinor field
##  @param h_in   Input spinor field
##  @param param  Contains all metadata regarding host and device
##                storage
##  @param parity The destination parity of the field
##  @param test_type Choose a type of dslash operators 
## 

proc dslashQuda_mdwf*(h_out: pointer; h_in: pointer; inv_param: ptr QudaInvertParam;
                     parity: QudaParity; test_type: cint) {.
    importc: "dslashQuda_mdwf", header: "quda.h".}
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
                parity: ptr QudaParity; inverse: cint) {.importc: "cloverQuda",
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

proc set_dim*(a2: ptr cint) {.importc: "set_dim", header: "quda.h".}
proc pack_ghost*(cpuLink: ptr pointer; cpuGhost: ptr pointer; nFace: cint;
                precision: QudaPrecision) {.importc: "pack_ghost", header: "quda.h".}
proc setFatLinkPadding*(`method`: QudaComputeFatMethod; param: ptr QudaGaugeParam) {.
    importc: "setFatLinkPadding", header: "quda.h".}
proc computeKSLinkQuda*(fatlink: pointer; longlink: pointer; ulink: pointer;
                       inlink: pointer; path_coeff: ptr cdouble;
                       param: ptr QudaGaugeParam; `method`: QudaComputeFatMethod) {.
    importc: "computeKSLinkQuda", header: "quda.h".}
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
##  action.  The momentum field is assumed to be in MILC order.
## 
##  @param momentum The momentum field
##  @param param The parameters of the external fields and the computation settings
##  @return momentum action
## 

proc momActionQuda*(momentum: pointer; param: ptr QudaGaugeParam): cdouble {.
    importc: "momActionQuda", header: "quda.h".}
## *
##  Take a gauge field on the host, load it onto the device and extend it.
##  Return a pointer to the extended gauge field object.
## 
##  @param gauge The CPU gauge field (optional - if set to 0 then the gauge field zeroed)
##  @param geometry The geometry of the matrix field to create (1 - scaler, 4 - vector, 6 - tensor)
##  @param param The parameters of the external field and the field to be created
##  @return Pointer to the gauge field (cast as a void*)
## 

proc createExtendedGaugeFieldQuda*(gauge: pointer; geometry: cint;
                                  param: ptr QudaGaugeParam): pointer {.
    importc: "createExtendedGaugeFieldQuda", header: "quda.h".}
## *
##  Allocate a gauge (matrix) field on the device and optionally download a host gauge field.
## 
##  @param gauge The host gauge field (optional - if set to 0 then the gauge field zeroed)
##  @param geometry The geometry of the matrix field to create (1 - scaler, 4 - vector, 6 - tensor)
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
##  Take a gauge field on the device and copy to the extended gauge
##  field.  The precisions and reconstruct types can differ between
##  the input and output field, but they must be compatible (same volume, geometry).
## 
##  @param outGauge Pointer to the output extended device gauge field (QUDA extended device field)
##  @param inGauge Pointer to the input device gauge field (QUDA gauge field)
## 

proc extendGaugeFieldQuda*(outGauge: pointer; inGauge: pointer) {.
    importc: "extendGaugeFieldQuda", header: "quda.h".}
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
##  Compute the sigma trace field (part of clover force computation).
##  All the pointers here are for QUDA native device objects.  The
##  precisions of all fields must match.  This function requires that
##  there is a persistent clover field.
##  
##  @param out Sigma trace field  (QUDA device field, geometry = 1)
##  @param dummy (not used)
##  @param mu mu direction
##  @param nu nu direction
##  @param dim array of local field dimensions
## 

proc computeCloverTraceQuda*(`out`: pointer; dummy: pointer; mu: cint; nu: cint;
                            dim: array[4, cint]) {.
    importc: "computeCloverTraceQuda", header: "quda.h".}
## *
##  Compute the derivative of the clover term (part of clover force
##  computation).  All the pointers here are for QUDA native device
##  objects.  The precisions of all fields must match.
##  
##  @param out Clover derivative field (QUDA device field, geometry = 1)
##  @param gauge Gauge field (extended QUDA device field, gemoetry = 4)
##  @param oprod Matrix field (outer product) which is multiplied by the derivative
##  @param mu mu direction
##  @param nu nu direction
##  @param coeff Coefficient of the clover derviative (including stepsize and clover coefficient)
##  @param parity Parity for which we are computing
##  @param param Gauge field meta data
##  @param conjugate Whether to make the oprod field anti-hermitian prior to multiplication
## 

proc computeCloverDerivativeQuda*(`out`: pointer; gauge: pointer; oprod: pointer;
                                 mu: cint; nu: cint; coeff: cdouble;
                                 parity: QudaParity; param: ptr QudaGaugeParam;
                                 conjugate: cint) {.
    importc: "computeCloverDerivativeQuda", header: "quda.h".}
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
##  Compute the quark-field outer product needed for gauge generation
##   
##  @param oprod The outer product to be computed.
##  @param quark The input fermion field.
##  @param num The number of quark fields
##  @param coeff The coefficient multiplying the fermion fields in the outer product
##  @param param The parameters of the outer-product field.
## 

proc computeStaggeredOprodQuda*(oprod: ptr pointer; quark: ptr pointer; num: cint;
                               coeff: ptr ptr cdouble; param: ptr QudaGaugeParam) {.
    importc: "computeStaggeredOprodQuda", header: "quda.h".}
## *
##  Compute the naive staggered force (experimental).  All fields are
##  QUDA device fields and must be in the same precision.
## 
##  mom Momentum field (QUDA device field)
##  quark Quark field solution vectors
##  coeff Step-size coefficient
## 

proc computeStaggeredForceQuda*(mom: pointer; quark: pointer; coeff: ptr cdouble) {.
    importc: "computeStaggeredForceQuda", header: "quda.h".}
## *
##  Compute the fermion force for the asqtad quark action. 
##  @param momentum          The momentum contribution from the quark action.
##  @param act_path_coeff    The coefficients that define the asqtad action.
##  @param one_link_src      The quark field outer product corresponding to the one-link term in the action. 
##  @param naik_src          The quark field outer product corresponding to the naik term in the action.
##  @param link              The gauge field.
##  @param param             The field parameters.
## 

proc computeAsqtadForceQuda*(momentum: pointer; flops: ptr clonglong;
                            act_path_coeff: array[6, cdouble];
                            one_link_src: array[4, pointer];
                            naik_src: array[4, pointer]; link: pointer;
                            param: ptr QudaGaugeParam) {.
    importc: "computeAsqtadForceQuda", header: "quda.h".}
## *
##  Compute the fermion force for the HISQ quark action. 
##  @param momentum        The momentum contribution from the quark action.
##  @param level2_coeff    The coefficients for the second level of smearing in the quark action.
##  @param fat7_coeff      The coefficients for the first level of smearing (fat7) in the quark action.
##  @param staple_src      Quark outer-product for the staple.
##  @param one_link_src    Quark outer-product for the one-link term in the action.
##  @param naik_src        Quark outer-product for the three-hop term in the action.
##  @param w_link          Unitarized link variables obtained by applying fat7 smearing and unitarization to the original links.
##  @param v_link          Fat7 link variables. 
##  @param u_link          SU(3) think link variables. 
##  @param param.          The field parameters.
## 

proc computeHISQForceQuda*(momentum: pointer; flops: ptr clonglong;
                          level2_coeff: array[6, cdouble];
                          fat7_coeff: array[6, cdouble];
                          staple_src: array[4, pointer];
                          one_link_src: array[4, pointer];
                          naik_src: array[4, pointer]; w_link: pointer;
                          v_link: pointer; u_link: pointer;
                          param: ptr QudaGaugeParam) {.
    importc: "computeHISQForceQuda", header: "quda.h".}
proc computeHISQForceCompleteQuda*(momentum: pointer;
                                  level2_coeff: array[6, cdouble];
                                  fat7_coeff: array[6, cdouble];
                                  quark_array: ptr pointer; num_terms: cint;
                                  quark_coeff: ptr ptr cdouble; w_link: pointer;
                                  v_link: pointer; u_link: pointer;
                                  param: ptr QudaGaugeParam) {.
    importc: "computeHISQForceCompleteQuda", header: "quda.h".}
## *
##  Computes the total, spatial and temporal plaquette averages of the loaded gauge configuration.
##  @param Array for storing the averages (total, spatial, temporal)
## 

proc plaqQuda*(plaq: array[3, cdouble]) {.importc: "plaqQuda", header: "quda.h".}
## *
##  Performs APE smearing on gaugePrecise and stores it in gaugeSmeared
##  @param nSteps Number of steps to apply.
##  @param alpha  Alpha coefficient for APE smearing.
## 

proc performAPEnStep*(nSteps: cuint; alpha: cdouble) {.importc: "performAPEnStep",
    header: "quda.h".}
## *
##  Calculates the topological charge from gaugeSmeared, if it exist, or from gaugePrecise if no smeared fields are present.
## 

proc qChargeCuda*(): cdouble {.importc: "qChargeCuda", header: "quda.h".}
## *
##  @brief Gauge fixing with overrelaxation with support for single and multi GPU.
##  @param[in,out] gauge, gauge field to be fixed
##  @param[in] gauge_dir, 3 for Coulomb gauge fixing, other for Landau gauge fixing
##  @param[in] Nsteps, maximum number of steps to perform gauge fixing
##  @param[in] verbose_interval, print gauge fixing info when iteration count is a multiple of this
##  @param[in] relax_boost, gauge fixing parameter of the overrelaxation method, most common value is 1.5 or 1.7.
##  @param[in] tolerance, torelance value to stop the method, if this value is zero then the method stops when iteration reachs the maximum number of steps defined by Nsteps
##  @param[in] reunit_interval, reunitarize gauge field when iteration count is a multiple of this
##  @param[in] stopWtheta, 0 for MILC criterium and 1 to use the theta value
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
##  @param[in] tolerance, torelance value to stop the method, if this value is zero then the method stops when iteration reachs the maximum number of steps defined by Nsteps
##  @param[in] stopWtheta, 0 for MILC criterium and 1 to use the theta value
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
##  Open/Close MAGMA library
## 
## 

proc openMagma*() {.importc: "openMagma", header: "quda.h".}
proc closeMagma*() {.importc: "closeMagma", header: "quda.h".}
## *
##  Clean deflation solver resources.
## 
## 

proc destroyDeflationQuda*(param: ptr QudaInvertParam; X: ptr cint; h_u: pointer;
                          inv_eigenvals: ptr cdouble) {.
    importc: "destroyDeflationQuda", header: "quda.h".}
##  #include <quda_new_interface.h>
