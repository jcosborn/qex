import enum_quda, quda

## *
##  @file    quda_milc_interface.h
##
##  @section Description
##
##  The header file defines the milc interface to enable easy
##  interfacing between QUDA and the MILC software packed.
##

when __COMPUTE_CAPABILITY__ >= 600:
  const
    USE_QUDA_MANAGED* = 1
## *
##  Parameters related to MILC site struct
##

type
  QudaMILCSiteArg_t* {.importc: "QudaMILCSiteArg_t",
                      header: "quda_milc_interface.h", bycopy.} = object
    site* {.importc: "site".}: pointer ## * Pointer to beginning of site array
    link* {.importc: "link".}: pointer ## * Pointer to link field (only used if site is not set)
    link_offset* {.importc: "link_offset".}: csize_t ## * Offset to link entry in site struct (bytes)
    mom* {.importc: "mom".}: pointer ## * Pointer to link field (only used if site is not set)
    mom_offset* {.importc: "mom_offset".}: csize_t ## * Offset to mom entry in site struct (bytes)
    size* {.importc: "size".}: csize_t ## * Size of site struct (bytes)


## *
##  Parameters related to linear solvers.
##

type
  QudaInvertArgs_t* {.importc: "QudaInvertArgs_t", header: "quda_milc_interface.h",
                     bycopy.} = object
    max_iter* {.importc: "max_iter".}: cint ## * Maximum number of iterations
    evenodd* {.importc: "evenodd".}: QudaParity ## * Which parity are we working on ? (options are QUDA_EVEN_PARITY, QUDA_ODD_PARITY, QUDA_INVALID_PARITY
    mixed_precision* {.importc: "mixed_precision".}: cint ## * Whether to use mixed precision or not (1 - yes, 0 - no)
    boundary_phase* {.importc: "boundary_phase".}: array[4, cdouble] ## * Boundary conditions
    make_resident_solution* {.importc: "make_resident_solution".}: cint ## * Make the solution resident and don't copy back
    use_resident_solution* {.importc: "use_resident_solution".}: cint ## * Use the resident solution
    solver_type* {.importc: "solver_type".}: QudaInverterType ## * Type of solver to use
    tadpole* {.importc: "tadpole".}: cdouble ## * Tadpole improvement factor - set to 1.0 for
                                         ##                         HISQ fermions since the tadpole factor is
                                         ##                         baked into the links during their construction
    naik_epsilon* {.importc: "naik_epsilon".}: cdouble ## * Naik epsilon parameter (HISQ fermions only).


## *
##  Parameters related to deflated solvers.
##

type
  QudaEigArgs_t* {.importc: "QudaEigArgs_t", header: "quda_milc_interface.h", bycopy.} = object
    prec_ritz* {.importc: "prec_ritz".}: QudaPrecision
    nev* {.importc: "nev".}: cint
    max_search_dim* {.importc: "max_search_dim".}: cint
    deflation_grid* {.importc: "deflation_grid".}: cint
    tol_restart* {.importc: "tol_restart".}: cdouble
    eigcg_max_restarts* {.importc: "eigcg_max_restarts".}: cint
    max_restart_num* {.importc: "max_restart_num".}: cint
    inc_tol* {.importc: "inc_tol".}: cdouble
    eigenval_tol* {.importc: "eigenval_tol".}: cdouble
    solver_ext_lib* {.importc: "solver_ext_lib".}: QudaExtLibType
    deflation_ext_lib* {.importc: "deflation_ext_lib".}: QudaExtLibType
    location_ritz* {.importc: "location_ritz".}: QudaFieldLocation
    mem_type_ritz* {.importc: "mem_type_ritz".}: QudaMemoryType
    vec_infile* {.importc: "vec_infile".}: cstring
    vec_outfile* {.importc: "vec_outfile".}: cstring


## *
##  Parameters related to problem size and machine topology.
##

type
  QudaLayout_t* {.importc: "QudaLayout_t", header: "quda_milc_interface.h", bycopy.} = object
    latsize* {.importc: "latsize".}: ptr cint ## * Local lattice dimensions
    machsize* {.importc: "machsize".}: ptr cint ## * Machine grid size
    device* {.importc: "device".}: cint ## * GPU device  number


## *
##  Parameters used to create a QUDA context.
##

type
  QudaInitArgs_t* {.importc: "QudaInitArgs_t", header: "quda_milc_interface.h",
                   bycopy.} = object
    verbosity* {.importc: "verbosity".}: QudaVerbosity ## * How verbose QUDA should be (QUDA_SILENT, QUDA_VERBOSE or QUDA_SUMMARIZE)
    layout* {.importc: "layout".}: QudaLayout_t ## * Layout for QUDA to use


##  passed to the initialization struct
## *
##  Parameters for defining HISQ calculations
##

type
  QudaHisqParams_t* {.importc: "QudaHisqParams_t", header: "quda_milc_interface.h",
                     bycopy.} = object
    reunit_allow_svd* {.importc: "reunit_allow_svd".}: cint ## * Allow SVD for reuniarization
    reunit_svd_only* {.importc: "reunit_svd_only".}: cint ## * Force use of SVD for reunitarization
    reunit_svd_abs_error* {.importc: "reunit_svd_abs_error".}: cdouble ## * Absolute error bound for SVD to apply
    reunit_svd_rel_error* {.importc: "reunit_svd_rel_error".}: cdouble ## * Relative error bound for SVD to apply
    force_filter* {.importc: "force_filter".}: cdouble ## * UV filter to apply to force


## *
##  Parameters for defining fat-link calculations
##

type
  QudaFatLinkArgs_t* {.importc: "QudaFatLinkArgs_t",
                      header: "quda_milc_interface.h", bycopy.} = object
    su3_source* {.importc: "su3_source".}: cint ## * is the incoming gauge field SU(3)
    use_pinned_memory* {.importc: "use_pinned_memory".}: cint ## * use page-locked memory in QUDA


## *
##  Optional: Set the MPI Comm Handle if it is not MPI_COMM_WORLD
##
##  @param input Pointer to an MPI_Comm handle, static cast as a void *
##

proc qudaSetMPICommHandle*(mycomm: pointer) {.importc: "qudaSetMPICommHandle",
    header: "quda_milc_interface.h".}
## *
##  Initialize the QUDA context.
##
##  @param input Meta data for the QUDA context
##

proc qudaInit*(input: QudaInitArgs_t) {.importc: "qudaInit",
                                     header: "quda_milc_interface.h".}
## *
##  Set set the local dimensions and machine topology for QUDA to use
##
##  @param layout Struct defining local dimensions and machine topology
##

proc qudaSetLayout*(layout: QudaLayout_t) {.importc: "qudaSetLayout",
    header: "quda_milc_interface.h".}
## *
##  Destroy the QUDA context.
##

proc qudaFinalize*() {.importc: "qudaFinalize", header: "quda_milc_interface.h".}
## *
##  Allocate pinned memory suitable for CPU-GPU transfers
##  @param bytes The size of the requested allocation
##  @return Pointer to allocated memory
##

proc qudaAllocatePinned*(bytes: csize_t): pointer {.importc: "qudaAllocatePinned",
    header: "quda_milc_interface.h".}
## *
##  Free pinned memory
##  @param ptr Pointer to memory to be free
##

proc qudaFreePinned*(`ptr`: pointer) {.importc: "qudaFreePinned",
                                    header: "quda_milc_interface.h".}
## *
##  Allocate managed memory to reduce CPU-GPU transfers
##  @param bytes The size of the requested allocation
##  @return Pointer to allocated memory
##

proc qudaAllocateManaged*(bytes: csize_t): pointer {.importc: "qudaAllocateManaged",
    header: "quda_milc_interface.h".}
## *
##  Free managed memory
##  @param ptr Pointer to memory to be free
##

proc qudaFreeManaged*(`ptr`: pointer) {.importc: "qudaFreeManaged",
                                     header: "quda_milc_interface.h".}
## *
##  Set the algorithms to use for HISQ fermion calculations, e.g.,
##  SVD parameters for reunitarization.
##
##  @param hisq_params Meta data desribing the algorithms to use for HISQ fermions
##

proc qudaHisqParamsInit*(hisq_params: QudaHisqParams_t) {.
    importc: "qudaHisqParamsInit", header: "quda_milc_interface.h".}
## *
##  Compute the fat and long links using the input gauge field.  All
##  fields passed here are host fields, that must be preallocated.
##  The precision of all fields must match.
##
##  @param precision The precision of the fields
##  @param fatlink_args Meta data for the algorithms to deploy
##  @param act_path_coeff Array of coefficients for each path in the action
##  @param inlink Host gauge field used for input
##  @param fatlink Host fat-link field that is computed
##  @param longlink Host long-link field that is computed
##

proc qudaLoadKSLink*(precision: cint; fatlink_args: QudaFatLinkArgs_t;
                    act_path_coeff: array[6, cdouble]; inlink: pointer;
                    fatlink: pointer; longlink: pointer) {.
    importc: "qudaLoadKSLink", header: "quda_milc_interface.h".}
## *
##  Compute the fat links and unitzarize using the input gauge field.
##  All fields passed here are host fields, that must be
##  preallocated.  The precision of all fields must match.
##
##  @param precision The precision of the fields
##  @param fatlink_args Meta data for the algorithms to deploy
##  @param path_coeff Array of coefficients for each path in the action
##  @param inlink Host gauge field used for input
##  @param fatlink Host fat-link field that is computed
##  @param ulink Host unitarized field that is computed
##

proc qudaLoadUnitarizedLink*(precision: cint; fatlink_args: QudaFatLinkArgs_t;
                            path_coeff: array[6, cdouble]; inlink: pointer;
                            fatlink: pointer; ulink: pointer) {.
    importc: "qudaLoadUnitarizedLink", header: "quda_milc_interface.h".}
## *
##  Apply the improved staggered operator to a field. All fields
##  passed and returned are host (CPU) field in MILC order.
##
##  @param external_precision Precision of host fields passed to QUDA (2 - double, 1 - single)
##  @param quda_precision Precision for QUDA to use (2 - double, 1 - single)
##  @param inv_args Struct setting some solver metadata
##  @param milc_fatlink Fat-link field on the host
##  @param milc_longlink Long-link field on the host
##  @param source Right-hand side source field
##  @param solution Solution spinor field
##

proc qudaDslash*(external_precision: cint; quda_precision: cint;
                inv_args: QudaInvertArgs_t; milc_fatlink: pointer;
                milc_longlink: pointer; source: pointer; solution: pointer;
                num_iters: ptr cint) {.importc: "qudaDslash",
                                    header: "quda_milc_interface.h".}
## *
##  Solve Ax=b using an improved staggered operator with a
##  domain-decomposition preconditioner.  All fields are fields
##  passed and returned are host (CPU) field in MILC order.  This
##  function requires that persistent gauge and clover fields have
##  been created prior.  This interface is experimental.
##
##  @param external_precision Precision of host fields passed to QUDA (2 - double, 1 - single)
##  @param precision Precision for QUDA to use (2 - double, 1 - single)
##  @param mass Fermion mass parameter
##  @param inv_args Struct setting some solver metadata
##  @param target_residual Target residual
##  @param target_relative_residual Target Fermilab residual
##  @param domain_overlap Array specifying the overlap of the domains in each dimension
##  @param fatlink Fat-link field on the host
##  @param longlink Long-link field on the host
##  @param source Right-hand side source field
##  @param solution Solution spinor field
##  @param final_residual True residual
##  @param final_relative_residual True Fermilab residual
##  @param num_iters Number of iterations taken
##

proc qudaDDInvert*(external_precision: cint; quda_precision: cint; mass: cdouble;
                  inv_args: QudaInvertArgs_t; target_residual: cdouble;
                  target_fermilab_residual: cdouble; domain_overlap: ptr cint;
                  fatlink: pointer; longlink: pointer; source: pointer;
                  solution: pointer; final_residual: ptr cdouble;
                  final_fermilab_residual: ptr cdouble; num_iters: ptr cint) {.
    importc: "qudaDDInvert", header: "quda_milc_interface.h".}
## *
##  Solve Ax=b for an improved staggered operator. All fields are fields
##  passed and returned are host (CPU) field in MILC order.  This
##  function requires that persistent gauge and clover fields have
##  been created prior.  This interface is experimental.
##
##  @param external_precision Precision of host fields passed to QUDA (2 - double, 1 - single)
##  @param quda_precision Precision for QUDA to use (2 - double, 1 - single)
##  @param mass Fermion mass parameter
##  @param inv_args Struct setting some solver metadata
##  @param target_residual Target residual
##  @param target_relative_residual Target Fermilab residual
##  @param milc_fatlink Fat-link field on the host
##  @param milc_longlink Long-link field on the host
##  @param source Right-hand side source field
##  @param solution Solution spinor field
##  @param final_residual True residual
##  @param final_relative_residual True Fermilab residual
##  @param num_iters Number of iterations taken
##

proc qudaInvert*(external_precision: cint; quda_precision: cint; mass: cdouble;
                inv_args: QudaInvertArgs_t; target_residual: cdouble;
                target_fermilab_residual: cdouble; milc_fatlink: pointer;
                milc_longlink: pointer; source: pointer; solution: pointer;
                final_resid: ptr cdouble; final_rel_resid: ptr cdouble;
                num_iters: ptr cint) {.importc: "qudaInvert",
                                    header: "quda_milc_interface.h".}
## *
##  Prepare a staggered/HISQ multigrid solve with given fat and
##  long links. All fields passed are host (CPU) fields
##  in MILC order. This function requires persistent gauge fields.
##  This interface is experimental.
##
##  @param external_precision Precision of host fields passed to QUDA (2 - double, 1 - single)
##  @param quda_precision Precision for QUDA to use (2 - double, 1 - single)
##  @param mass Fermion mass parameter
##  @param inv_args Struct setting some solver metadata; required for tadpole, naik coeff
##  @param milc_fatlink Fat-link field on the host
##  @param milc_longlink Long-link field on the host
##  @param mg_param_file Path to an input text file describing the MG solve, to be documented on QUDA wiki
##  @return Void pointer wrapping a pack of multigrid-related structures
##

proc qudaMultigridCreate*(external_precision: cint; quda_precision: cint;
                         mass: cdouble; inv_args: QudaInvertArgs_t;
                         milc_fatlink: pointer; milc_longlink: pointer;
                         mg_param_file: cstring): pointer {.
    importc: "qudaMultigridCreate", header: "quda_milc_interface.h".}
## *
##  Solve Ax=b for an improved staggered operator using MG.
##  All fields are fields passed and returned are host (CPU)
##  field in MILC order.  This function requires that persistent
##  gauge and clover fields have been created prior. It also
##  requires a multigrid parameter built from qudaSetupMultigrid
##  This interface is experimental.
##
##  @param external_precision Precision of host fields passed to QUDA (2 - double, 1 - single)
##  @param quda_precision Precision for QUDA to use (2 - double, 1 - single)
##  @param mass Fermion mass parameter
##  @param inv_args Struct setting some solver metadata
##  @param target_residual Target residual
##  @param target_relative_residual Target Fermilab residual
##  @param milc_fatlink Fat-link field on the host
##  @param milc_longlink Long-link field on the host
##  @param mg_pack_ptr MG preconditioner structure created by qudaSetupMultigrid
##  @param mg_rebuild_type whether to do a full (1) or thin (0) MG rebuild
##  @param source Right-hand side source field
##  @param solution Solution spinor field
##  @param final_residual True residual
##  @param final_relative_residual True Fermilab residual
##  @param num_iters Number of iterations taken
##

proc qudaInvertMG*(external_precision: cint; quda_precision: cint; mass: cdouble;
                  inv_args: QudaInvertArgs_t; target_residual: cdouble;
                  target_fermilab_residual: cdouble; milc_fatlink: pointer;
                  milc_longlink: pointer; mg_pack_ptr: pointer;
                  mg_rebuild_type: cint; source: pointer; solution: pointer;
                  final_residual: ptr cdouble;
                  final_fermilab_residual: ptr cdouble; num_iters: ptr cint) {.
    importc: "qudaInvertMG", header: "quda_milc_interface.h".}
## *
##  Clean up a staggered/HISQ multigrid object, freeing all internal
##  fields and otherwise allocated memory.
##
##  @param mg_pack_ptr Void pointer mapping to the multigrid structure returned by qudaSetupMultigrid
##

proc qudaMultigridDestroy*(mg_pack_ptr: pointer) {.importc: "qudaMultigridDestroy",
    header: "quda_milc_interface.h".}
## *
##  Solve Ax=b for an improved staggered operator with many right hand sides.
##  All fields are fields passed and returned are host (CPU) field in MILC order.
##  This function requires that persistent gauge and clover fields have
##  been created prior.  This interface is experimental.
##
##  @param external_precision Precision of host fields passed to QUDA (2 - double, 1 - single)
##  @param quda_precision Precision for QUDA to use (2 - double, 1 - single)
##  @param mass Fermion mass parameter
##  @param inv_args Struct setting some solver metadata
##  @param target_residual Target residual
##  @param target_relative_residual Target Fermilab residual
##  @param milc_fatlink Fat-link field on the host
##  @param milc_longlink Long-link field on the host
##  @param source array of right-hand side source fields
##  @param solution array of solution spinor fields
##  @param final_residual True residual
##  @param final_relative_residual True Fermilab residual
##  @param num_iters Number of iterations taken
##  @param num_src Number of source fields
##

proc qudaInvertMsrc*(external_precision: cint; quda_precision: cint; mass: cdouble;
                    inv_args: QudaInvertArgs_t; target_residual: cdouble;
                    target_fermilab_residual: cdouble; fatlink: pointer;
                    longlink: pointer; sourceArray: ptr pointer;
                    solutionArray: ptr pointer; final_residual: ptr cdouble;
                    final_fermilab_residual: ptr cdouble; num_iters: ptr cint;
                    num_src: cint) {.importc: "qudaInvertMsrc",
                                   header: "quda_milc_interface.h".}
## *
##  Solve for multiple shifts (e.g., masses) using an improved
##  staggered operator.  All fields are fields passed and returned
##  are host (CPU) field in MILC order.  This function requires that
##  persistent gauge and clover fields have been created prior.  When
##  a pure double-precision solver is requested no reliable updates
##  are used, else reliable updates are used with a reliable_delta
##  parameter of 0.1.
##
##  @param external_precision Precision of host fields passed to QUDA (2 - double, 1 - single)
##  @param precision Precision for QUDA to use (2 - double, 1 - single)
##  @param num_offsets Number of shifts to solve for
##  @param offset Array of shift offset values
##  @param inv_args Struct setting some solver metadata
##  @param target_residual Array of target residuals per shift
##  @param target_relative_residual Array of target Fermilab residuals per shift
##  @param milc_fatlink Fat-link field on the host
##  @param milc_longlink Long-link field on the host
##  @param source Right-hand side source field
##  @param solutionArray Array of solution spinor fields
##  @param final_residual Array of true residuals
##  @param final_relative_residual Array of true Fermilab residuals
##  @param num_iters Number of iterations taken
##

proc qudaMultishiftInvert*(external_precision: cint; precision: cint;
                          num_offsets: cint; offset: ptr cdouble;
                          inv_args: QudaInvertArgs_t;
                          target_residual: ptr cdouble;
                          target_fermilab_residual: ptr cdouble;
                          milc_fatlink: pointer; milc_longlink: pointer;
                          source: pointer; solutionArray: ptr pointer;
                          final_residual: ptr cdouble;
                          final_fermilab_residual: ptr cdouble; num_iters: ptr cint) {.
    importc: "qudaMultishiftInvert", header: "quda_milc_interface.h".}
## *
##  Solve for a system with many RHS using an improved
##  staggered operator.
##  The solving procedure consists of two computation phases :
##  1) incremental pahse : call eigCG solver to accumulate low eigenmodes
##  2) deflation phase : use computed eigenmodes to deflate a regular CG
##  All fields are fields passed and returned
##  are host (CPU) field in MILC order.  This function requires that
##  persistent gauge and clover fields have been created prior.
##
##  @param external_precision Precision of host fields passed to QUDA (2 - double, 1 - single)
##  @param precision Precision for QUDA to use (2 - double, 1 - single)
##  @param num_offsets Number of shifts to solve for
##  @param offset Array of shift offset values
##  @param inv_args Struct setting some solver metadata
##  @param target_residual Array of target residuals per shift
##  @param target_relative_residual Array of target Fermilab residuals per shift
##  @param milc_fatlink Fat-link field on the host
##  @param milc_longlink Long-link field on the host
##  @param source Right-hand side source field
##  @param solution Array of solution spinor fields
##  @param eig_args contains info about deflation space
##  @param rhs_idx  bookkeep current rhs
##  @param last_rhs_flag  is this the last rhs to solve?
##  @param final_residual Array of true residuals
##  @param final_relative_residual Array of true Fermilab residuals
##  @param num_iters Number of iterations taken
##

proc qudaEigCGInvert*(external_precision: cint; quda_precision: cint; mass: cdouble;
                     inv_args: QudaInvertArgs_t; target_residual: cdouble;
                     target_fermilab_residual: cdouble; fatlink: pointer;
                     longlink: pointer; source: pointer; solution: pointer;
                     eig_args: QudaEigArgs_t; rhs_idx: cint; last_rhs_flag: cint;
                     final_residual: ptr cdouble;
                     final_fermilab_residual: ptr cdouble; num_iters: ptr cint) {.
    importc: "qudaEigCGInvert", header: "quda_milc_interface.h".}
  ## current rhs
  ## is this the last rhs to solve?
## *
##  Solve Ax=b using a Wilson-Clover operator.  All fields are fields
##  passed and returned are host (CPU) field in MILC order.  This
##  function creates the gauge and clover field from the host fields.
##  Reliable updates are used with a reliable_delta parameter of 0.1.
##
##  @param external_precision Precision of host fields passed to QUDA (2 - double, 1 - single)
##  @param quda_precision Precision for QUDA to use (2 - double, 1 - single)
##  @param kappa Kappa value
##  @param clover_coeff Clover coefficient
##  @param inv_args Struct setting some solver metadata
##  @param target_residual Target residual
##  @param milc_link Gauge field on the host
##  @param milc_clover Clover field on the host
##  @param milc_clover_inv Inverse clover on the host
##  @param clover_coeff Clover coefficient
##  @param source Right-hand side source field
##  @param solution Solution spinor field
##  @param final_residual True residual returned by the solver
##  @param final_residual True Fermilab residual returned by the solver
##  @param num_iters Number of iterations taken
##

proc qudaCloverInvert*(external_precision: cint; quda_precision: cint;
                      kappa: cdouble; clover_coeff: cdouble;
                      inv_args: QudaInvertArgs_t; target_residual: cdouble;
                      target_fermilab_residual: cdouble; milc_link: pointer;
                      milc_clover: pointer; milc_clover_inv: pointer;
                      source: pointer; solution: pointer;
                      final_residual: ptr cdouble;
                      final_fermilab_residual: ptr cdouble; num_iters: ptr cint) {.
    importc: "qudaCloverInvert", header: "quda_milc_interface.h".}
## *
##  Solve for a system with many RHS using using a Wilson-Clover operator.
##  The solving procedure consists of two computation phases :
##  1) incremental pahse : call eigCG solver to accumulate low eigenmodes
##  2) deflation phase : use computed eigenmodes to deflate a regular CG
##  All fields are fields passed and returned
##  are host (CPU) field in MILC order.  This function requires that
##  persistent gauge and clover fields have been created prior.
##
##  @param external_precision Precision of host fields passed to QUDA (2 - double, 1 - single)
##  @param quda_precision Precision for QUDA to use (2 - double, 1 - single)
##  @param kappa Kappa value
##  @param clover_coeff Clover coefficient
##  @param inv_args Struct setting some solver metadata
##  @param target_residual Target residual
##  @param milc_link Gauge field on the host
##  @param milc_clover Clover field on the host
##  @param milc_clover_inv Inverse clover on the host
##  @param clover_coeff Clover coefficient
##  @param source Right-hand side source field
##  @param solution Solution spinor field
##  @param eig_args contains info about deflation space
##  @param rhs_idx  bookkeep current rhs
##  @param last_rhs_flag  is this the last rhs to solve?
##  @param final_residual Array of true residuals
##  @param final_relative_residual Array of true Fermilab residuals
##  @param num_iters Number of iterations taken
##

proc qudaEigCGCloverInvert*(external_precision: cint; quda_precision: cint;
                           kappa: cdouble; clover_coeff: cdouble;
                           inv_args: QudaInvertArgs_t; target_residual: cdouble;
                           target_fermilab_residual: cdouble; milc_link: pointer;
                           milc_clover: pointer; milc_clover_inv: pointer;
                           source: pointer; solution: pointer;
                           eig_args: QudaEigArgs_t; rhs_idx: cint;
                           last_rhs_flag: cint; final_residual: ptr cdouble;
                           final_fermilab_residual: ptr cdouble;
                           num_iters: ptr cint) {.importc: "qudaEigCGCloverInvert",
    header: "quda_milc_interface.h".}
  ## current rhs
  ## is this the last rhs to solve?
## *
##  Load the gauge field from the host.
##
##  @param external_precision Precision of host fields passed to QUDA (2 - double, 1 - single)
##  @param quda_precision Precision for QUDA to use (2 - double, 1 - single)
##  @param inv_args Meta data
##  @param milc_link Base pointer to host gauge field (regardless of dimensionality)
##

proc qudaLoadGaugeField*(external_precision: cint; quda_precision: cint;
                        inv_args: QudaInvertArgs_t; milc_link: pointer) {.
    importc: "qudaLoadGaugeField", header: "quda_milc_interface.h".}
## *
##      Free the gauge field allocated in QUDA.
##

proc qudaFreeGaugeField*() {.importc: "qudaFreeGaugeField",
                           header: "quda_milc_interface.h".}
## *
##  Load the clover field and its inverse from the host.  If null
##  pointers are passed, the clover field and / or its inverse will
##  be computed dynamically from the resident gauge field.
##
##  @param external_precision Precision of host fields passed to QUDA (2 - double, 1 - single)
##  @param quda_precision Precision for QUDA to use (2 - double, 1 - single)
##  @param inv_args Meta data
##  @param milc_clover Pointer to host clover field.  If 0 then the
##  clover field is computed dynamically within QUDA.
##  @param milc_clover_inv Pointer to host inverse clover field.  If
##  0 then the inverse if computed dynamically within QUDA.
##  @param solution_type The type of solution required  (mat, matpc)
##  @param solve_type The solve type to use (normal/direct/preconditioning)
##  @param clover_coeff Clover coefficient
##  @param compute_trlog Whether to compute the trlog of the clover field when inverting
##  @param Array for storing the trlog (length two, one for each parity)
##

proc qudaLoadCloverField*(external_precision: cint; quda_precision: cint;
                         inv_args: QudaInvertArgs_t; milc_clover: pointer;
                         milc_clover_inv: pointer;
                         solution_type: QudaSolutionType;
                         solve_type: QudaSolveType; clover_coeff: cdouble;
                         compute_trlog: cint; trlog: ptr cdouble) {.
    importc: "qudaLoadCloverField", header: "quda_milc_interface.h".}
## *
##      Free the clover field allocated in QUDA.
##

proc qudaFreeCloverField*() {.importc: "qudaFreeCloverField",
                            header: "quda_milc_interface.h".}
## *
##  Solve for multiple shifts (e.g., masses) using a Wilson-Clover
##  operator with multi-shift CG.  All fields are fields passed and
##  returned are host (CPU) field in MILC order.  This function
##  requires that persistent gauge and clover fields have been
##  created prior.  When a pure double-precision solver is requested
##  no reliable updates are used, else reliable updates are used with
##  a reliable_delta parameter of 0.1.
##
##  @param external_precision Precision of host fields passed to QUDA (2 - double, 1 - single)
##  @param quda_precision Precision for QUDA to use (2 - double, 1 - single)
##  @param num_offsets Number of shifts to solve for
##  @param offset Array of shift offset values
##  @param kappa Kappa value
##  @param clover_coeff Clover coefficient
##  @param inv_args Struct setting some solver metadata
##  @param target_residual Array of target residuals per shift
##  @param milc_link Ignored
##  @param milc_clover Ignored
##  @param milc_clover_inv Ignored
##  @param clover_coeff Clover coefficient
##  @param source Right-hand side source field
##  @param solutionArray Array of solution spinor fields
##  @param final_residual Array of true residuals
##  @param num_iters Number of iterations taken
##

proc qudaCloverMultishiftInvert*(external_precision: cint; quda_precision: cint;
                                num_offsets: cint; offset: ptr cdouble;
                                kappa: cdouble; clover_coeff: cdouble;
                                inv_args: QudaInvertArgs_t;
                                target_residual: ptr cdouble; milc_link: pointer;
                                milc_clover: pointer; milc_clover_inv: pointer;
                                source: pointer; solutionArray: ptr pointer;
                                final_residual: ptr cdouble; num_iters: ptr cint) {.
    importc: "qudaCloverMultishiftInvert", header: "quda_milc_interface.h".}
## *
##  Compute the fermion force for the HISQ quark action.  All fields
##  are host fields in MILC order, and the precision of these fields
##  must match.
##
##  @param precision       The precision of the fields
##  @param num_terms The number of quark fields
##  @param num_naik_terms The number of naik contributions
##  @param dt Integrating step size
##  @param coeff The coefficients multiplying the fermion fields in the outer product
##  @param quark_field The input fermion field.
##  @param level2_coeff    The coefficients for the second level of smearing in the quark action.
##  @param fat7_coeff      The coefficients for the first level of smearing (fat7) in the quark action.
##  @param w_link          Unitarized link variables obtained by applying fat7 smearing and unitarization to the original links.
##  @param v_link          Fat7 link variables.
##  @param u_link          SU(3) think link variables.
##  @param milc_momentum        The momentum contribution from the quark action.
##

proc qudaHisqForce*(precision: cint; num_terms: cint; num_naik_terms: cint;
                   dt: cdouble; coeff: ptr ptr cdouble; quark_field: ptr pointer;
                   level2_coeff: array[6, cdouble]; fat7_coeff: array[6, cdouble];
                   w_link: pointer; v_link: pointer; u_link: pointer;
                   milc_momentum: pointer) {.importc: "qudaHisqForce",
    header: "quda_milc_interface.h".}
## *
##  Compute the gauge force and update the momentum field.  All fields
##  here are CPU fields in MILC order, and their precisions should
##  match.
##
##  @param precision The precision of the field (2 - double, 1 - single)
##  @param num_loop_types 1, 2 or 3
##  @param milc_loop_coeff Coefficients of the different loops in the Symanzik action
##  @param eb3 The integration step size (for MILC this is dt*beta/3)
##  @param arg Metadata for MILC's internal site struct array
##

proc qudaGaugeForce*(precision: cint; num_loop_types: cint;
                    milc_loop_coeff: array[3, cdouble]; eb3: cdouble;
                    arg: ptr QudaMILCSiteArg_t) {.importc: "qudaGaugeForce",
    header: "quda_milc_interface.h".}
## *
##  Compute the gauge force and update the momentum field.  All fields
##  here are CPU fields in MILC order, and their precisions should
##  match.
##
##  @param precision The precision of the field (2 - double, 1 - single)
##  @param num_loop_types 1, 2 or 3
##  @param milc_loop_coeff Coefficients of the different loops in the Symanzik action
##  @param eb3 The integration step size (for MILC this is dt*beta/3)
##  @param arg Metadata for MILC's internal site struct array
##  @param phase_in whether staggered phases are applied
##

proc qudaGaugeForcePhased*(precision: cint; num_loop_types: cint;
                          milc_loop_coeff: array[3, cdouble]; eb3: cdouble;
                          arg: ptr QudaMILCSiteArg_t; phase_in: cint) {.
    importc: "qudaGaugeForcePhased", header: "quda_milc_interface.h".}
## *
##  Evolve the gauge field by step size dt, using the momentum field
##  I.e., Evalulate U(t+dt) = e(dt pi) U(t).  All fields are CPU fields in MILC order.
##
##  @param precision Precision of the field (2 - double, 1 - single)
##  @param dt The integration step size step
##  @param arg Metadata for MILC's internal site struct array
##

proc qudaUpdateU*(precision: cint; eps: cdouble; arg: ptr QudaMILCSiteArg_t) {.
    importc: "qudaUpdateU", header: "quda_milc_interface.h".}
## *
##  Evolve the gauge field by step size dt, using the momentum field
##  I.e., Evalulate U(t+dt) = e(dt pi) U(t).  All fields are CPU fields in MILC order.
##
##  @param precision Precision of the field (2 - double, 1 - single)
##  @param dt The integration step size step
##  @param arg Metadata for MILC's internal site struct array
##  @param phase_in whether staggered phases are applied
##

proc qudaUpdateUPhased*(precision: cint; eps: cdouble; arg: ptr QudaMILCSiteArg_t;
                       phase_in: cint) {.importc: "qudaUpdateUPhased",
                                       header: "quda_milc_interface.h".}
## *
##  Evolve the gauge field by step size dt, using the momentum field
##  I.e., Evalulate U(t+dt) = e(dt pi) U(t).  All fields are CPU fields in MILC order.
##
##  @param precision Precision of the field (2 - double, 1 - single)
##  @param dt The integration step size step
##  @param arg Metadata for MILC's internal site struct array
##  @param phase_in whether staggered phases are applied
##  @param want_gaugepipe whether to enabled QUDA gaugepipe for HMC
##

proc qudaUpdateUPhasedPipeline*(precision: cint; eps: cdouble;
                               arg: ptr QudaMILCSiteArg_t; phase_in: cint;
                               want_gaugepipe: cint) {.
    importc: "qudaUpdateUPhasedPipeline", header: "quda_milc_interface.h".}
## *
##  Download the momentum from MILC and place into QUDA's resident
##  momentum field.  The source momentum field can either be as part
##  of a MILC site struct (QUDA_MILC_SITE_GAUGE_ORDER) or as a
##  separate field (QUDA_MILC_GAUGE_ORDER).
##
##  @param precision Precision of the field (2 - double, 1 - single)
##  @param arg Metadata for MILC's internal site struct array
##

proc qudaMomLoad*(precision: cint; arg: ptr QudaMILCSiteArg_t) {.
    importc: "qudaMomLoad", header: "quda_milc_interface.h".}
## *
##  Upload the momentum to MILC from QUDA's resident momentum field.
##  The destination momentum field can either be as part of a MILC site
##  struct (QUDA_MILC_SITE_GAUGE_ORDER) or as a separate field
##  (QUDA_MILC_GAUGE_ORDER).
##
##  @param precision Precision of the field (2 - double, 1 - single)
##  @param arg Metadata for MILC's internal site struct array
##

proc qudaMomSave*(precision: cint; arg: ptr QudaMILCSiteArg_t) {.
    importc: "qudaMomSave", header: "quda_milc_interface.h".}
## *
##  Evaluate the momentum contribution to the Hybrid Monte Carlo
##  action.  MILC convention is applied, subtracting 4.0 from each
##  momentum matrix to increase stability.
##
##  @param precision Precision of the field (2 - double, 1 - single)
##  @param arg Metadata for MILC's internal site struct array
##  @return momentum action
##

proc qudaMomAction*(precision: cint; arg: ptr QudaMILCSiteArg_t): cdouble {.
    importc: "qudaMomAction", header: "quda_milc_interface.h".}
## *
##  Apply the staggered phase factors to the gauge field.  If the
##  imaginary chemical potential is non-zero then the phase factor
##  exp(imu/T) will be applied to the links in the temporal
##  direction.
##
##  @param prec Precision of the gauge field
##  @param gauge_h The gauge field
##  @param flag Whether to apply to remove the staggered phase
##  @param i_mu Imaginary chemical potential
##

proc qudaRephase*(prec: cint; gauge: pointer; flag: cint; i_mu: cdouble) {.
    importc: "qudaRephase", header: "quda_milc_interface.h".}
## *
##  Project the input field on the SU(3) group.  If the target
##  tolerance is not met, this routine will give a runtime error.
##
##  @param prec Precision of the gauge field
##  @param tol The tolerance to which we iterate
##  @param arg Metadata for MILC's internal site struct array
##

proc qudaUnitarizeSU3*(prec: cint; tol: cdouble; arg: ptr QudaMILCSiteArg_t) {.
    importc: "qudaUnitarizeSU3", header: "quda_milc_interface.h".}
## *
##  Project the input field on the SU(3) group.  If the target
##  tolerance is not met, this routine will give a runtime error.
##
##  @param prec Precision of the gauge field
##  @param tol The tolerance to which we iterate
##  @param arg Metadata for MILC's internal site struct array
##  @param phase_in whether staggered phases are applied
##

proc qudaUnitarizeSU3Phased*(prec: cint; tol: cdouble; arg: ptr QudaMILCSiteArg_t;
                            phase_in: cint) {.importc: "qudaUnitarizeSU3Phased",
    header: "quda_milc_interface.h".}
## *
##  Compute the clover force contributions in each dimension mu given
##  the array solution fields, and compute the resulting momentum
##  field.
##
##  @param mom Momentum matrix
##  @param dt Integrating step size
##  @param x Array of solution vectors
##  @param p Array of intermediate vectors
##  @param coeff Array of residues for each contribution
##  @param kappa kappa parameter
##  @param ck -clover_coefficient * kappa / 8
##  @param nvec Number of vectors
##  @param multiplicity Number of fermions represented by this bilinear
##  @param gauge Gauge Field
##  @param precision Precision of the fields
##  @param inv_args Struct setting some solver metadata
##

proc qudaCloverForce*(mom: pointer; dt: cdouble; x: ptr pointer; p: ptr pointer;
                     coeff: ptr cdouble; kappa: cdouble; ck: cdouble; nvec: cint;
                     multiplicity: cdouble; gauge: pointer; precision: cint;
                     inv_args: QudaInvertArgs_t) {.importc: "qudaCloverForce",
    header: "quda_milc_interface.h".}
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
##

proc qudaCloverTrace*(`out`: pointer; dummy: pointer; mu: cint; nu: cint) {.
    importc: "qudaCloverTrace", header: "quda_milc_interface.h".}
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
##  @param precision Precision of the fields (2 = double, 1 = single)
##  @param parity Parity for which we are computing
##  @param conjugate Whether to make the oprod field anti-hermitian prior to multiplication
##

proc qudaCloverDerivative*(`out`: pointer; gauge: pointer; oprod: pointer; mu: cint;
                          nu: cint; coeff: cdouble; precision: cint; parity: cint;
                          conjugate: cint) {.importc: "qudaCloverDerivative",
    header: "quda_milc_interface.h".}
## *
##  Take a gauge field on the host, load it onto the device and extend it.
##  Return a pointer to the extended gauge field object.
##
##  @param gauge The CPU gauge field (optional - if set to 0 then the gauge field zeroed)
##  @param geometry The geometry of the matrix field to create (1 - scaler, 4 - vector, 6 - tensor)
##  @param precision The precision of the fields (2 - double, 1 - single)
##  @return Pointer to the gauge field (cast as a void*)
##

proc qudaCreateExtendedGaugeField*(gauge: pointer; geometry: cint; precision: cint): pointer {.
    importc: "qudaCreateExtendedGaugeField", header: "quda_milc_interface.h".}
## *
##  Take the QUDA resident gauge field and extend it.
##  Return a pointer to the extended gauge field object.
##
##  @param gauge The CPU gauge field (optional - if set to 0 then the gauge field zeroed)
##  @param geometry The geometry of the matrix field to create (1 - scaler, 4 - vector, 6 - tensor)
##  @param precision The precision of the fields (2 - double, 1 - single)
##  @return Pointer to the gauge field (cast as a void*)
##

proc qudaResidentExtendedGaugeField*(gauge: pointer; geometry: cint; precision: cint): pointer {.
    importc: "qudaResidentExtendedGaugeField", header: "quda_milc_interface.h".}
## *
##  Allocate a gauge (matrix) field on the device and optionally download a host gauge field.
##
##  @param gauge The host gauge field (optional - if set to 0 then the gauge field zeroed)
##  @param geometry The geometry of the matrix field to create (1 - scaler, 4 - vector, 6 - tensor)
##  @param precision The precision of the field to be created (2 - double, 1 - single)
##  @return Pointer to the gauge field (cast as a void*)
##

proc qudaCreateGaugeField*(gauge: pointer; geometry: cint; precision: cint): pointer {.
    importc: "qudaCreateGaugeField", header: "quda_milc_interface.h".}
## *
##  Copy the QUDA gauge (matrix) field on the device to the CPU
##
##  @param outGauge Pointer to the host gauge field
##  @param inGauge Pointer to the device gauge field (QUDA device field)
##

proc qudaSaveGaugeField*(gauge: pointer; inGauge: pointer) {.
    importc: "qudaSaveGaugeField", header: "quda_milc_interface.h".}
## *
##  Reinterpret gauge as a pointer to cudaGaugeField and call destructor.
##
##  @param gauge Gauge field to be freed
##

proc qudaDestroyGaugeField*(gauge: pointer) {.importc: "qudaDestroyGaugeField",
    header: "quda_milc_interface.h".}
## *
##  @brief Gauge fixing with overrelaxation with support for single and multi GPU.
##  @param[in] precision, 1 for single precision else for double precision
##  @param[in] gauge_dir, 3 for Coulomb gauge fixing, other for Landau gauge fixing
##  @param[in] Nsteps, maximum number of steps to perform gauge fixing
##  @param[in] verbose_interval, print gauge fixing info when iteration count is a multiple of this
##  @param[in] relax_boost, gauge fixing parameter of the overrelaxation method, most common value is 1.5 or 1.7.
##  @param[in] tolerance, torelance value to stop the method, if this value is zero then the method stops when iteration reachs the maximum number of steps defined by Nsteps
##  @param[in] reunit_interval, reunitarize gauge field when iteration count is a multiple of this
##  @param[in] stopWtheta, 0 for MILC criterium and 1 to use the theta value
##  @param[in,out] milc_sitelink, MILC gauge field to be fixed
##

proc qudaGaugeFixingOVR*(precision: cint; gauge_dir: cuint; Nsteps: cint;
                        verbose_interval: cint; relax_boost: cdouble;
                        tolerance: cdouble; reunit_interval: cuint;
                        stopWtheta: cuint; milc_sitelink: pointer) {.
    importc: "qudaGaugeFixingOVR", header: "quda_milc_interface.h".}
## *
##  @brief Gauge fixing with Steepest descent method with FFTs with support for single GPU only.
##  @param[in] precision, 1 for single precision else for double precision
##  @param[in] gauge_dir, 3 for Coulomb gauge fixing, other for Landau gauge fixing
##  @param[in] Nsteps, maximum number of steps to perform gauge fixing
##  @param[in] verbose_interval, print gauge fixing info when iteration count is a multiple of this
##  @param[in] alpha, gauge fixing parameter of the method, most common value is 0.08
##  @param[in] autotune, 1 to autotune the method, i.e., if the Fg inverts its tendency we decrease the alpha value
##  @param[in] tolerance, torelance value to stop the method, if this value is zero then the method stops when iteration reachs the maximum number of steps defined by Nsteps
##  @param[in] stopWtheta, 0 for MILC criterium and 1 to use the theta value
##  @param[in,out] milc_sitelink, MILC gauge field to be fixed
##

proc qudaGaugeFixingFFT*(precision: cint; gauge_dir: cuint; Nsteps: cint;
                        verbose_interval: cint; alpha: cdouble; autotune: cuint;
                        tolerance: cdouble; stopWtheta: cuint;
                        milc_sitelink: pointer) {.importc: "qudaGaugeFixingFFT",
    header: "quda_milc_interface.h".}
##  The below declarations are for removed functions from prior versions of QUDA.
## *
##  Note this interface function has been removed.  This stub remains
##  for compatibility only.
##

proc qudaAsqtadForce*(precision: cint; act_path_coeff: array[6, cdouble];
                     one_link_src: array[4, pointer]; naik_src: array[4, pointer];
                     link: pointer; milc_momentum: pointer) {.
    importc: "qudaAsqtadForce", header: "quda_milc_interface.h".}
## *
##  Note this interface function has been removed.  This stub remains
##  for compatibility only.
##

proc qudaComputeOprod*(precision: cint; num_terms: cint; num_naik_terms: cint;
                      coeff: ptr ptr cdouble; scale: cdouble;
                      quark_field: ptr pointer; oprod: array[3, pointer]) {.
    importc: "qudaComputeOprod", header: "quda_milc_interface.h".}