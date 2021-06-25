import qex
import quda, enum_quda, quda_constants
import bitops
import std/decls
#{.experimental: "views".}
proc cstrcpy(dest: ptr char, src: ptr char): ptr char {.discardable,importc:"strcpy".}

template printfQuda(args: varargs[untyped]) =
  cprintf(args)
template errorQuda(args: varargs[untyped]) =
  cprintf(args)

var xdim*,ydim*,zdim*,tdim*: int32
var Lsdim* = 16
var verbosity* = QUDA_SUMMARIZE
#var verbosity* = QUDA_DEBUG_VERBOSE
var cpu_prec* = QUDA_DOUBLE_PRECISION
var cuda_prec* = QUDA_DOUBLE_PRECISION
var cuda_prec_sloppy* = QUDA_SINGLE_PRECISION
var cuda_prec_precondition* = QUDA_SINGLE_PRECISION
#var cuda_prec_eigensolver* = QUDA_DOUBLE_PRECISION
var cuda_prec_eigensolver* = QUDA_SINGLE_PRECISION
var cuda_prec_refinement_sloppy* = QUDA_SINGLE_PRECISION
var prec_sloppy* = QUDA_SINGLE_PRECISION
var prec_refinement_sloppy* = QUDA_SINGLE_PRECISION
var prec_precondition* = QUDA_SINGLE_PRECISION
var prec_eigensolver* = QUDA_DOUBLE_PRECISION
var prec_null* = QUDA_SINGLE_PRECISION
var link_recon* = QUDA_RECONSTRUCT_NO
var link_recon_sloppy* = QUDA_RECONSTRUCT_NO
var link_recon_precondition* = QUDA_RECONSTRUCT_NO
var link_recon_eigensolver* = QUDA_RECONSTRUCT_NO
var fermion_t_boundary = QUDA_ANTI_PERIODIC_T
var dslash_type = QUDA_WILSON_DSLASH
var laplace3D = 4'i32
var niter = 100
var maxiter_precondition = 10
var gcrNkrylov = 10
var ca_basis = QUDA_POWER_BASIS
var ca_lambda_min = 0.0
var ca_lambda_max = -1.0
var pipeline = 0
#var solution_accumulator_pipeline = 0
var solution_accumulator_pipeline = 1
var precon_type = QUDA_INVALID_INVERTER
var precon_schwarz_type = QUDA_INVALID_SCHWARZ
var precon_schwarz_cycle = 1
var multishift = 1
var mass* = 0.1
var kappa* = -1.0
var mu = 0.1
var epsilon = 0.01
var m5 = -1.5
var b5 = 1.5
var c5 = 0.5
var anisotropy = 1.0
var eps_naik = 0.0
var clover_coeff = 0.1
var tol = 1e-7
var tol_restart = 5e3 * tol
var tol_precondition = 1e-1
var tol_hq = 0.0
var reliable_delta = 0.1
var alternative_reliable = false
var twist_flavor = QUDA_TWIST_SINGLET
var eofa_pm = 1
var eofa_shift = -1.2345
var eofa_mq1 = 1.0
var eofa_mq2 = 0.085
var eofa_mq3 = 1.0

var matpc_type = QUDA_MATPC_EVEN_EVEN
#var solve_type* = QUDA_NORMOP_PC_SOLVE
var solve_type* = QUDA_DIRECT_SOLVE

var mg_levels* = 2

type mgarray[T] = array[QUDA_MAX_MG_LEVEL,T]

var nvec: mgarray[int]
var mg_vec_infile: mgarray[string]
var mg_vec_outfile: mgarray[string]
var solver_location: mgarray[QudaFieldLocation]
var setup_location: mgarray[QudaFieldLocation]
var nu_pre: mgarray[int]
var nu_post: mgarray[int]
var n_block_ortho: mgarray[int]
var mu_factor: mgarray[float]
var mg_verbosity: mgarray[QudaVerbosity]
var setup_inv: mgarray[QudaInverterType]
var coarse_solve_type: mgarray[QudaSolveType]
var smoother_solve_type: mgarray[QudaSolveType]
var num_setup_iter: mgarray[int]
var setup_tol: mgarray[float]
var setup_maxiter: mgarray[int]
var setup_maxiter_refresh: mgarray[int]
var setup_ca_basis: mgarray[QudaCABasis]
var setup_ca_basis_size: mgarray[int]
var setup_ca_lambda_min: mgarray[float]
var setup_ca_lambda_max: mgarray[float]
var setup_type = QUDA_NULL_VECTOR_SETUP
var pre_orthonormalize = false
var post_orthonormalize = true
var omega = 0.85
var coarse_solver: mgarray[QudaInverterType]
var coarse_solver_tol: mgarray[float]
var smoother_type: mgarray[QudaInverterType]
var smoother_halo_prec* = QUDA_SINGLE_PRECISION
var smoother_tol: mgarray[float]
var coarse_solver_maxiter: mgarray[int]
var coarse_solver_ca_basis: mgarray[QudaCABasis]
var coarse_solver_ca_basis_size: mgarray[int]
var coarse_solver_ca_lambda_min: mgarray[float]
var coarse_solver_ca_lambda_max: mgarray[float]
var generate_nullspace = true
var generate_all_levels = true
var mg_schwarz_type: mgarray[QudaSchwarzType]
var mg_schwarz_cycle: mgarray[int]
var mg_evolve_thin_updates = false
var verify_results = true
var low_mode_check = false
#var low_mode_check = true
var oblique_proj_check = false
#var oblique_proj_check = true
var mg_use_mma = false
var mg_eig_coarse_guess = false
var geo_block_size: mgarray[array[4,int32]]

proc setLat*(x: openArray[int]) =
  xdim = int32 x[0]
  ydim = int32 x[1]
  zdim = int32 x[2]
  tdim = int32 x[3]

proc setGaugeParam*(gauge_param: var QudaGaugeParam) =
  gauge_param.`type` = QUDA_SU3_LINKS
  gauge_param.X[0] = xdim
  gauge_param.X[1] = ydim
  gauge_param.X[2] = zdim
  gauge_param.X[3] = tdim
  gauge_param.cpu_prec = cpu_prec
  gauge_param.cuda_prec = cuda_prec
  gauge_param.cuda_prec_sloppy = cuda_prec
  gauge_param.cuda_prec_precondition = cuda_prec
  gauge_param.cuda_prec_eigensolver = cuda_prec
  gauge_param.reconstruct = link_recon
  gauge_param.reconstruct_sloppy = link_recon
  gauge_param.reconstruct_precondition = link_recon
  gauge_param.reconstruct_eigensolver = link_recon
  gauge_param.reconstruct_refinement_sloppy = link_recon
  gauge_param.anisotropy = 1.0
  gauge_param.tadpole_coeff = 1.0
  gauge_param.ga_pad = 0
  gauge_param.mom_ga_pad = 0
  gauge_param.gauge_fix = QUDA_GAUGE_FIXED_NO

proc setWilsonGaugeParam*(gauge_param: var QudaGaugeParam) =
  setGaugeParam(gauge_param)
  gauge_param.anisotropy = anisotropy
  gauge_param.`type` = QUDA_WILSON_LINKS
  gauge_param.gauge_order = QUDA_QDP_GAUGE_ORDER
  gauge_param.t_boundary = fermion_t_boundary
  gauge_param.cuda_prec_sloppy = cuda_prec_sloppy
  gauge_param.cuda_prec_precondition = cuda_prec_precondition
  gauge_param.cuda_prec_eigensolver = cuda_prec_eigensolver
  gauge_param.cuda_prec_refinement_sloppy = cuda_prec_refinement_sloppy
  gauge_param.reconstruct_sloppy = link_recon_sloppy
  gauge_param.reconstruct_precondition = link_recon_precondition
  gauge_param.reconstruct_eigensolver = link_recon_eigensolver
  gauge_param.reconstruct_refinement_sloppy = link_recon_sloppy
  var pad_size: cint = 0
  ##  For multi-GPU, ga_pad must be large enough to store a time-slice
  if nRanks > 1:
    var x_face_size = gauge_param.X[1] * gauge_param.X[2] * gauge_param.X[3] div 2
    var y_face_size = gauge_param.X[0] * gauge_param.X[2] * gauge_param.X[3] div 2
    var z_face_size = gauge_param.X[0] * gauge_param.X[1] * gauge_param.X[3] div 2
    var t_face_size = gauge_param.X[0] * gauge_param.X[1] * gauge_param.X[2] div 2
    pad_size = max(max(x_face_size, y_face_size), max(z_face_size, t_face_size))
  gauge_param.ga_pad = pad_size

proc setStaggeredGaugeParam*(gauge_param: var QudaGaugeParam) =
  setGaugeParam(gauge_param)
  gauge_param.cuda_prec_sloppy = prec_sloppy
  gauge_param.cuda_prec_refinement_sloppy = prec_refinement_sloppy
  gauge_param.cuda_prec_precondition = prec_precondition
  gauge_param.cuda_prec_eigensolver = prec_eigensolver
  gauge_param.reconstruct_sloppy = link_recon_sloppy
  gauge_param.reconstruct_precondition = link_recon_precondition
  gauge_param.reconstruct_eigensolver = link_recon_eigensolver
  gauge_param.reconstruct_refinement_sloppy = link_recon_sloppy
  ##  For HISQ, this must always be set to 1.0, since the tadpole
  ##  correction is baked into the coefficients for the first fattening.
  ##  The tadpole doesn't mean anything for the second fattening
  ##  since the input fields are unitarized.
  gauge_param.tadpole_coeff = 1.0
  if dslash_type == QUDA_ASQTAD_DSLASH:
    gauge_param.scale = -1.0/24.0
    if eps_naik != 0:
      gauge_param.scale = gauge_param.scale * (1.0 + eps_naik)
  else:
    gauge_param.scale = 1.0
  gauge_param.gauge_order = QUDA_MILC_GAUGE_ORDER
  gauge_param.t_boundary = fermion_t_boundary
  gauge_param.staggered_phase_type = QUDA_STAGGERED_PHASE_MILC
  gauge_param.`type` = QUDA_WILSON_LINKS
  var pad_size: cint = 0
  when defined(MULTI_GPU):
    var x_face_size = gauge_param.X[1] * gauge_param.X[2] * gauge_param.X[3] div 2
    var y_face_size = gauge_param.X[0] * gauge_param.X[2] * gauge_param.X[3] div 2
    var z_face_size = gauge_param.X[0] * gauge_param.X[1] * gauge_param.X[3] div 2
    var t_face_size = gauge_param.X[0] * gauge_param.X[1] * gauge_param.X[2] div 2
    pad_size = max(max(x_face_size, y_face_size), max(z_face_size, t_face_size))
  gauge_param.ga_pad = pad_size

proc setInvertParam*(inv_param: var QudaInvertParam) =
  ##  Set dslash type
  inv_param.dslash_type = dslash_type
  ##  Use kappa or mass normalisation
  if kappa == -1.0:
    inv_param.mass = mass
    inv_param.kappa = 1.0/(2.0 * (1 + 3/anisotropy + mass))
    if dslash_type == QUDA_LAPLACE_DSLASH:
      inv_param.kappa = 1.0/(8 + mass)
  else:
    inv_param.kappa = kappa
    inv_param.mass = 0.5/kappa - (1.0 + 3.0/anisotropy)
    if dslash_type == QUDA_LAPLACE_DSLASH:
      inv_param.mass = 1.0/kappa - 8.0
  printfQuda("Kappa = %.8f Mass = %.8f\n", inv_param.kappa, inv_param.mass)
  ##  Use 3D or 4D laplace
  inv_param.laplace3D = laplace3D
  ##  Some fermion specific parameters
  if dslash_type == QUDA_TWISTED_MASS_DSLASH or
      dslash_type == QUDA_TWISTED_CLOVER_DSLASH:
    inv_param.mu = mu
    inv_param.epsilon = epsilon
    inv_param.twist_flavor = twist_flavor
    inv_param.Ls = if (inv_param.twist_flavor == QUDA_TWIST_NONDEG_DOUBLET): 2 else: 1
  elif dslash_type == QUDA_DOMAIN_WALL_DSLASH or
      dslash_type == QUDA_DOMAIN_WALL_4D_DSLASH or
      dslash_type == QUDA_MOBIUS_DWF_DSLASH or
      dslash_type == QUDA_MOBIUS_DWF_EOFA_DSLASH:
    inv_param.m5 = m5
    #kappa5 = 0.5 div (5 + inv_param.m5)
    inv_param.Ls = int32 Lsdim
    for k in 0..<Lsdim:
      ##  for mobius only
      ##  b5[k], c[k] values are chosen for arbitrary values,
      ##  but the difference of them are same as 1.0
      inv_param.b_5[k] = [b5,0]
      inv_param.c_5[k] = [c5,0]
    inv_param.eofa_pm = int32 eofa_pm
    inv_param.eofa_shift = eofa_shift
    inv_param.mq1 = eofa_mq1
    inv_param.mq2 = eofa_mq2
    inv_param.mq3 = eofa_mq3
  else:
    inv_param.Ls = 1
  ##  Set clover specific parameters
  if dslash_type == QUDA_CLOVER_WILSON_DSLASH or
      dslash_type == QUDA_TWISTED_CLOVER_DSLASH:
    inv_param.clover_cpu_prec = cpu_prec
    inv_param.clover_cuda_prec = cuda_prec
    inv_param.clover_cuda_prec_sloppy = cuda_prec_sloppy
    inv_param.clover_cuda_prec_precondition = cuda_prec_precondition
    inv_param.clover_cuda_prec_eigensolver = cuda_prec_eigensolver
    inv_param.clover_cuda_prec_refinement_sloppy = cuda_prec_sloppy
    inv_param.clover_order = QUDA_PACKED_CLOVER_ORDER
    inv_param.clover_coeff = clover_coeff
  #inv_param.inv_type = inv_type
  inv_param.inv_type = QUDA_BICGSTAB_INVERTER
  #inv_param.solution_type = solution_type
  inv_param.solution_type = QUDA_MAT_SOLUTION
  inv_param.solve_type = solve_type
  inv_param.matpc_type = matpc_type
  inv_param.dagger = QUDA_DAG_NO
  #inv_param.mass_normalization = normalization
  inv_param.mass_normalization = QUDA_KAPPA_NORMALIZATION
  inv_param.solver_normalization = QUDA_DEFAULT_NORMALIZATION
  inv_param.pipeline = int32 pipeline
  inv_param.Nsteps = 2
  inv_param.gcrNkrylov = int32 gcrNkrylov
  inv_param.ca_basis = ca_basis
  inv_param.ca_lambda_min = ca_lambda_min
  inv_param.ca_lambda_max = ca_lambda_max
  inv_param.tol = tol
  inv_param.tol_restart = tol_restart
  if tol_hq == 0 and tol == 0:
    errorQuda("qudaInvert: requesting zero residual\n")
    quit(1)
  inv_param.residual_type = QudaResidualType 0
  if tol != 0:
    inv_param.residual_type = QudaResidualType(
      bitor(ord(inv_param.residual_type), ord(QUDA_L2_RELATIVE_RESIDUAL)) )
  if tol_hq != 0:
    inv_param.residual_type = QudaResidualType(
      bitor(ord(inv_param.residual_type), ord(QUDA_HEAVY_QUARK_RESIDUAL)) )
  inv_param.tol_hq = tol_hq
  ##  specify a tolerance for the residual for heavy quark residual
  ##  Offsets used only by multi-shift solver
  ##  These should be set in the application code. We set the them here by way of
  ##  example
  inv_param.num_offset = int32 multishift
  for i in 0..<inv_param.num_offset:
    inv_param.offset[i] = 0.06 + i * i * 0.1
  ##  these can be set individually
  for i in 0..<inv_param.num_offset:
    inv_param.tol_offset[i] = inv_param.tol
    inv_param.tol_hq_offset[i] = inv_param.tol_hq
  inv_param.maxiter = int32 niter
  inv_param.reliable_delta = reliable_delta
  inv_param.reliable_delta_refinement = reliable_delta
  inv_param.use_alternative_reliable = int32 ord alternative_reliable
  inv_param.use_sloppy_partial_accumulator = 0
  inv_param.solution_accumulator_pipeline = int32 solution_accumulator_pipeline
  inv_param.max_res_increase = 1
  ##  domain decomposition preconditioner parameters
  inv_param.inv_type_precondition = precon_type
  inv_param.schwarz_type = precon_schwarz_type
  inv_param.precondition_cycle = int32 precon_schwarz_cycle
  inv_param.tol_precondition = tol_precondition
  inv_param.maxiter_precondition = int32 maxiter_precondition
  inv_param.verbosity_precondition = mg_verbosity[0]
  inv_param.cuda_prec_precondition = cuda_prec_precondition
  inv_param.cuda_prec_eigensolver = cuda_prec_eigensolver
  inv_param.omega = 1.0
  inv_param.cpu_prec = cpu_prec
  inv_param.cuda_prec = cuda_prec
  inv_param.cuda_prec_sloppy = cuda_prec_sloppy
  inv_param.cuda_prec_refinement_sloppy = cuda_prec_refinement_sloppy
  inv_param.preserve_source = QUDA_PRESERVE_SOURCE_YES
  inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS
  inv_param.dirac_order = QUDA_DIRAC_ORDER
  inv_param.input_location = QUDA_CPU_FIELD_LOCATION
  inv_param.output_location = QUDA_CPU_FIELD_LOCATION
  inv_param.sp_pad = 0
  inv_param.cl_pad = 0
  inv_param.verbosity = verbosity
  inv_param.extlib_type = QUDA_EIGEN_EXTLIB
  inv_param.native_blas_lapack = QUDA_BOOLEAN_TRUE
  #inv_param.native_blas_lapack = QUDA_BOOLEAN_FALSE

proc setQudaMgSolveTypes*() =
  for i in 0..<QUDA_MAX_MG_LEVEL:
    if coarse_solve_type[i] == QUDA_INVALID_SOLVE:
      coarse_solve_type[i] = solve_type
    if smoother_solve_type[i] == QUDA_INVALID_SOLVE:
      smoother_solve_type[i] = QUDA_DIRECT_PC_SOLVE

proc setQudaDefaultMgTestParams*() =
  # We give here some default values
  for i in 0..<QUDA_MAX_MG_LEVEL:
    mg_verbosity[i] = QUDA_SUMMARIZE
    #mg_verbosity[i] = QUDA_DEBUG_VERBOSE
    setup_inv[i] = QUDA_BICGSTAB_INVERTER
    num_setup_iter[i] = 1
    setup_tol[i] = 5e-6
    setup_maxiter[i] = 500
    setup_maxiter_refresh[i] = 20
    mu_factor[i] = 1.0
    coarse_solve_type[i] = QUDA_INVALID_SOLVE
    smoother_solve_type[i] = QUDA_INVALID_SOLVE
    mg_schwarz_type[i] = QUDA_INVALID_SCHWARZ
    mg_schwarz_cycle[i] = 1
    smoother_type[i] = QUDA_GCR_INVERTER
    #smoother_type[i] = QUDA_CA_GCR_INVERTER
    smoother_tol[i] = 0.25
    coarse_solver[i] = QUDA_GCR_INVERTER
    coarse_solver_tol[i] = 0.25
    coarse_solver_maxiter[i] = 100
    solver_location[i] = QUDA_CUDA_FIELD_LOCATION
    setup_location[i] = QUDA_CUDA_FIELD_LOCATION
    nu_pre[i] = 0
    nu_post[i] = 8
    n_block_ortho[i] = 1
    mg_vec_infile[i] = ""
    mg_vec_outfile[i] = ""

    # Default eigensolver params
    #[
    mg_eig[i] = false
    mg_eig_tol[i] = 1e-3
    mg_eig_n_ev[i] = nvec[i]
    mg_eig_n_kr[i] = 3 * nvc[i]
    mg_eig_require_convergence[i] = QUDA_BOOLEAN_TRUE
    mg_eig_type[i] = QUDA_EIG_TR_LANCZOS
    mg_eig_spectrum[i] = QUDA_SPECTRUM_SR_EIG
    mg_eig_check_interval[i] = 5
    mg_eig_max_restarts[i] = 100
    mg_eig_use_normop[i] = QUDA_BOOLEAN_FALSE
    mg_eig_use_dagger[i] = QUDA_BOOLEAN_FALSE
    mg_eig_use_poly_acc[i] = QUDA_BOOLEAN_TRUE
    mg_eig_poly_deg[i] = 100
    mg_eig_amin[i] = 1.0
    mg_eig_amax[i] = -1.0  # use power iterations
    mg_eig_save_prec[i] = QUDA_DOUBLE_PRECISION
    ]#

    setup_ca_basis[i] = QUDA_POWER_BASIS
    setup_ca_basis_size[i] = 4
    setup_ca_lambda_min[i] = 0.0
    setup_ca_lambda_max[i] = -1.0  # use power iterations

    coarse_solver_ca_basis[i] = QUDA_POWER_BASIS
    coarse_solver_ca_basis_size[i] = 4
    coarse_solver_ca_lambda_min[i] = 0.0
    coarse_solver_ca_lambda_max[i] = -1.0

proc setMultigridInvertParam*(inv_param: var QudaInvertParam) =
  inv_param.Ls = 1
  inv_param.sp_pad = 0
  inv_param.cl_pad = 0
  inv_param.cpu_prec = cpu_prec
  inv_param.cuda_prec = cuda_prec
  inv_param.cuda_prec_sloppy = cuda_prec_sloppy
  inv_param.cuda_prec_precondition = cuda_prec_precondition
  inv_param.cuda_prec_eigensolver = cuda_prec_eigensolver
  inv_param.preserve_source = QUDA_PRESERVE_SOURCE_NO
  inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS
  inv_param.dirac_order = QUDA_DIRAC_ORDER
  if dslash_type == QUDA_CLOVER_WILSON_DSLASH or
      dslash_type == QUDA_TWISTED_CLOVER_DSLASH:
    inv_param.clover_cpu_prec = cpu_prec
    inv_param.clover_cuda_prec = cuda_prec
    inv_param.clover_cuda_prec_sloppy = cuda_prec_sloppy
    inv_param.clover_cuda_prec_precondition = cuda_prec_precondition
    inv_param.clover_cuda_prec_eigensolver = cuda_prec_eigensolver
    inv_param.clover_cuda_prec_refinement_sloppy = cuda_prec_sloppy
    inv_param.clover_order = QUDA_PACKED_CLOVER_ORDER
  inv_param.input_location = QUDA_CPU_FIELD_LOCATION
  inv_param.output_location = QUDA_CPU_FIELD_LOCATION
  inv_param.dslash_type = dslash_type
  if kappa == -1.0:
    inv_param.mass = mass
    inv_param.kappa = 1.0/(2.0 * (1 + 3/anisotropy + mass))
  else:
    inv_param.kappa = kappa
    inv_param.mass = 0.5/kappa - (1 + 3/anisotropy)
  if dslash_type == QUDA_TWISTED_MASS_DSLASH or
      dslash_type == QUDA_TWISTED_CLOVER_DSLASH:
    inv_param.mu = mu
    inv_param.epsilon = epsilon
    inv_param.twist_flavor = twist_flavor
    inv_param.Ls = if (inv_param.twist_flavor == QUDA_TWIST_NONDEG_DOUBLET): 2 else: 1
    if twist_flavor == QUDA_TWIST_NONDEG_DOUBLET:
      printfQuda("Twisted-mass doublet non supported (yet)\n")
      quit(0)
  inv_param.clover_coeff = clover_coeff
  inv_param.dagger = QUDA_DAG_NO
  inv_param.mass_normalization = QUDA_KAPPA_NORMALIZATION
  ##  do we want full solution or single-parity solution
  inv_param.solution_type = QUDA_MAT_SOLUTION
  ##  do we want to use an even-odd preconditioned solve or not
  inv_param.solve_type = solve_type
  inv_param.matpc_type = matpc_type
  inv_param.inv_type = QUDA_GCR_INVERTER
  inv_param.verbosity = verbosity
  inv_param.verbosity_precondition = mg_verbosity[0]
  inv_param.inv_type_precondition = QUDA_MG_INVERTER
  inv_param.pipeline = int32 pipeline
  inv_param.gcrNkrylov = int32 gcrNkrylov
  inv_param.tol = tol
  ##  require both L2 relative and heavy quark residual to determine convergence
  inv_param.residual_type = QUDA_L2_RELATIVE_RESIDUAL
  inv_param.tol_hq = tol_hq
  ##  specify a tolerance for the residual for heavy quark residual
  ##  Offsets used only by multi-shift solver
  ##  should be set in application
  inv_param.num_offset = int32 multishift
  for i in 0..<inv_param.num_offset:
    inv_param.offset[i] = 0.06 + i * i * 0.1
  ##  these can be set individually
  for i in 0..<inv_param.num_offset:
    inv_param.tol_offset[i] = inv_param.tol
    inv_param.tol_hq_offset[i] = inv_param.tol_hq
  inv_param.maxiter = int32 niter
  inv_param.reliable_delta = reliable_delta
  ##  domain decomposition preconditioner is disabled when using MG
  inv_param.schwarz_type = QUDA_INVALID_SCHWARZ
  inv_param.precondition_cycle = 1
  inv_param.tol_precondition = 1e-1
  inv_param.maxiter_precondition = 1
  inv_param.omega = 1.0
  inv_param.native_blas_lapack = QUDA_BOOLEAN_TRUE
  #inv_param.native_blas_lapack = QUDA_BOOLEAN_FALSE

proc setMultigridParam*(mg_param: var QudaMultigridParam) =
  #var inv_param {.byaddr.} = mg_param.invert_param
  template inv_param: untyped = mg_param.invert_param[]
  ##  this will be used to setup SolverParam parent in MGParam class
  ##  Whether or not to use native BLAS LAPACK
  inv_param.native_blas_lapack = QUDA_BOOLEAN_TRUE
  #inv_param.native_blas_lapack = QUDA_BOOLEAN_FALSE
  inv_param.Ls = 1
  inv_param.sp_pad = 0
  inv_param.cl_pad = 0
  inv_param.cpu_prec = cpu_prec
  inv_param.cuda_prec = cuda_prec
  inv_param.cuda_prec_sloppy = cuda_prec_sloppy
  inv_param.cuda_prec_precondition = cuda_prec_precondition
  inv_param.cuda_prec_eigensolver = cuda_prec_eigensolver
  inv_param.preserve_source = QUDA_PRESERVE_SOURCE_NO
  inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS
  inv_param.dirac_order = QUDA_DIRAC_ORDER
  if dslash_type == QUDA_CLOVER_WILSON_DSLASH or
      dslash_type == QUDA_TWISTED_CLOVER_DSLASH:
    inv_param.clover_cpu_prec = cpu_prec
    inv_param.clover_cuda_prec = cuda_prec
    inv_param.clover_cuda_prec_sloppy = cuda_prec_sloppy
    inv_param.clover_cuda_prec_precondition = cuda_prec_precondition
    inv_param.clover_cuda_prec_eigensolver = cuda_prec_eigensolver
    inv_param.clover_cuda_prec_refinement_sloppy = cuda_prec_sloppy
    inv_param.clover_order = QUDA_PACKED_CLOVER_ORDER
    inv_param.clover_coeff = clover_coeff
  inv_param.input_location = QUDA_CPU_FIELD_LOCATION
  inv_param.output_location = QUDA_CPU_FIELD_LOCATION
  inv_param.dslash_type = dslash_type
  if kappa == -1.0:
    inv_param.mass = mass
    inv_param.kappa = 1.0/(2.0 * (1 + 3/anisotropy + mass))
  else:
    inv_param.kappa = kappa
    inv_param.mass = 0.5/kappa - (1 + 3/anisotropy)
  if dslash_type == QUDA_TWISTED_MASS_DSLASH or
      dslash_type == QUDA_TWISTED_CLOVER_DSLASH:
    inv_param.mu = mu
    inv_param.epsilon = epsilon
    inv_param.twist_flavor = twist_flavor
    inv_param.Ls = if (inv_param.twist_flavor == QUDA_TWIST_NONDEG_DOUBLET): 2 else: 1
    if twist_flavor == QUDA_TWIST_NONDEG_DOUBLET:
      printfQuda("Twisted-mass doublet non supported (yet)\n")
      quit(0)
  inv_param.dagger = QUDA_DAG_NO
  inv_param.mass_normalization = QUDA_KAPPA_NORMALIZATION
  inv_param.matpc_type = matpc_type
  inv_param.solution_type = QUDA_MAT_SOLUTION
  inv_param.solve_type = QUDA_DIRECT_SOLVE
  #mg_param.invert_param = inv_param
  mg_param.n_level = int32 mg_levels
  for i in 0..<mg_param.n_level:
    for j in 0..<4:
      ##  if not defined use 4
      mg_param.geo_block_size[i][j] = if geo_block_size[i][j]>0:
                                        geo_block_size[i][j] else: 4
    for j in 4..<QUDA_MAX_DIM:
      mg_param.geo_block_size[i][j] = 1
    #mg_param.use_eig_solver[i] = if mg_eig[i]: QUDA_BOOLEAN_TRUE else:QUDA_BOOLEAN_FALSE
    mg_param.use_eig_solver[i] = QUDA_BOOLEAN_FALSE
    mg_param.verbosity[i] = mg_verbosity[i]
    mg_param.setup_inv_type[i] = setup_inv[i]
    mg_param.num_setup_iter[i] = int32 num_setup_iter[i]
    mg_param.setup_tol[i] = setup_tol[i]
    mg_param.setup_maxiter[i] = int32 setup_maxiter[i]
    mg_param.setup_maxiter_refresh[i] = int32 setup_maxiter_refresh[i]
    ##  Basis to use for CA-CGN(E/R) setup
    mg_param.setup_ca_basis[i] = setup_ca_basis[i]
    ##  Basis size for CACG setup
    mg_param.setup_ca_basis_size[i] = int32 setup_ca_basis_size[i]
    ##  Minimum and maximum eigenvalue for Chebyshev CA basis setup
    mg_param.setup_ca_lambda_min[i] = setup_ca_lambda_min[i]
    mg_param.setup_ca_lambda_max[i] = setup_ca_lambda_max[i]
    mg_param.spin_block_size[i] = 1
    mg_param.n_vec[i] = if nvec[i] == 0: 24 else: nvec[i]
    ##  default to 24 vectors if not set
    mg_param.n_block_ortho[i] = int32 n_block_ortho[i]
    ##  number of times to Gram-Schmidt
    mg_param.precision_null[i] = prec_null
    ##  precision to store the null-space basis
    mg_param.smoother_halo_precision[i] = smoother_halo_prec
    ##  precision of the halo exchange in the smoother
    mg_param.nu_pre[i] = int32 nu_pre[i]
    mg_param.nu_post[i] = int32 nu_post[i]
    mg_param.mu_factor[i] = mu_factor[i]
    mg_param.cycle_type[i] = QUDA_MG_CYCLE_RECURSIVE
    ##  Is not a staggered solve, always aggregate
    mg_param.transfer_type[i] = QUDA_TRANSFER_AGGREGATE
    ##  set the coarse solver wrappers including bottom solver
    mg_param.coarse_solver[i] = coarse_solver[i]
    mg_param.coarse_solver_tol[i] = coarse_solver_tol[i]
    mg_param.coarse_solver_maxiter[i] = int32 coarse_solver_maxiter[i]
    ##  Basis to use for CA-CGN(E/R) coarse solver
    mg_param.coarse_solver_ca_basis[i] = coarse_solver_ca_basis[i]
    ##  Basis size for CACG coarse solver/
    mg_param.coarse_solver_ca_basis_size[i] = int32 coarse_solver_ca_basis_size[i]
    ##  Minimum and maximum eigenvalue for Chebyshev CA basis
    mg_param.coarse_solver_ca_lambda_min[i] = coarse_solver_ca_lambda_min[i]
    mg_param.coarse_solver_ca_lambda_max[i] = coarse_solver_ca_lambda_max[i]
    mg_param.smoother[i] = smoother_type[i]
    ##  set the smoother / bottom solver tolerance (for MR smoothing this will be ignored)
    mg_param.smoother_tol[i] = smoother_tol[i]
    ##  set to QUDA_DIRECT_SOLVE for no even/odd preconditioning on the smoother
    ##  set to QUDA_DIRECT_PC_SOLVE for to enable even/odd preconditioning on the smoother
    mg_param.smoother_solve_type[i] = smoother_solve_type[i]
    ##  set to QUDA_ADDITIVE_SCHWARZ for Additive Schwarz precondioned smoother (presently only impelemented for MR)
    mg_param.smoother_schwarz_type[i] = mg_schwarz_type[i]
    ##  if using Schwarz preconditioning then use local reductions only
    mg_param.global_reduction[i] = if (mg_schwarz_type[i] == QUDA_INVALID_SCHWARZ): QUDA_BOOLEAN_TRUE else: QUDA_BOOLEAN_FALSE
    ##  set number of Schwarz cycles to apply
    mg_param.smoother_schwarz_cycle[i] = int32 mg_schwarz_cycle[i]
    ##  Set set coarse_grid_solution_type: this defines which linear
    ##  system we are solving on a given level
    ##  * QUDA_MAT_SOLUTION - we are solving the full system and inject
    ##    a full field into coarse grid
    ##  * QUDA_MATPC_SOLUTION - we are solving the e/o-preconditioned
    ##    system, and only inject single parity field into coarse grid
    ##
    ##  Multiple possible scenarios here
    ##
    ##  1. **Direct outer solver and direct smoother**: here we use
    ##  full-field residual coarsening, and everything involves the
    ##  full system so coarse_grid_solution_type = QUDA_MAT_SOLUTION
    ##
    ##  2. **Direct outer solver and preconditioned smoother**: here,
    ##  only the smoothing uses e/o preconditioning, so
    ##  coarse_grid_solution_type = QUDA_MAT_SOLUTION_TYPE.
    ##  We reconstruct the full residual prior to coarsening after the
    ##  pre-smoother, and then need to project the solution for post
    ##  smoothing.
    ##
    ##  3. **Preconditioned outer solver and preconditioned smoother**:
    ##  here we use single-parity residual coarsening throughout, so
    ##  coarse_grid_solution_type = QUDA_MATPC_SOLUTION.  This is a bit
    ##  questionable from a theoretical point of view, since we don't
    ##  coarsen the preconditioned operator directly, rather we coarsen
    ##  the full operator and preconditioned that, but it just works.
    ##  This is the optimal combination in general for Wilson-type
    ##  operators: although there is an occasional increase in
    ##  iteration or two), by working completely in the preconditioned
    ##  space, we save the cost of reconstructing the full residual
    ##  from the preconditioned smoother, and re-projecting for the
    ##  subsequent smoother, as well as reducing the cost of the
    ##  ancillary blas operations in the coarse-grid solve.
    ##
    ##  Note, we cannot use preconditioned outer solve with direct
    ##  smoother
    ##
    ##  Finally, we have to treat the top level carefully: for all
    ##  other levels the entry into and out of the grid will be a
    ##  full-field, which we can then work in Schur complement space or
    ##  not (e.g., freedom to choose coarse_grid_solution_type).  For
    ##  the top level, if the outer solver is for the preconditioned
    ##  system, then we must use preconditoning, e.g., option 3.) above.
    if i == 0:
      ##  top-level treatment
      if coarse_solve_type[0] != solve_type:
        errorQuda("Mismatch between top-level MG solve type %d and outer solve type %d",
                  coarse_solve_type[0], solve_type)
      if solve_type == QUDA_DIRECT_SOLVE:
        mg_param.coarse_grid_solution_type[i] = QUDA_MAT_SOLUTION
      elif solve_type == QUDA_DIRECT_PC_SOLVE:
        mg_param.coarse_grid_solution_type[i] = QUDA_MATPC_SOLUTION
      else:
        errorQuda("Unexpected solve_type = %d\n", solve_type)
    else:
      if coarse_solve_type[i] == QUDA_DIRECT_SOLVE:
        mg_param.coarse_grid_solution_type[i] = QUDA_MAT_SOLUTION
      elif coarse_solve_type[i] == QUDA_DIRECT_PC_SOLVE:
        mg_param.coarse_grid_solution_type[i] = QUDA_MATPC_SOLUTION
      else:
        errorQuda("Unexpected solve_type = %d\n", coarse_solve_type[i])
    mg_param.omega[i] = omega
    ##  over/under relaxation factor
    mg_param.location[i] = solver_location[i]
    mg_param.setup_location[i] = setup_location[i]
  ##  whether to run GPU setup but putting temporaries into mapped (slow CPU) memory
  mg_param.setup_minimize_memory = QUDA_BOOLEAN_FALSE
  ##  only coarsen the spin on the first restriction
  mg_param.spin_block_size[0] = 2
  mg_param.setup_type = setup_type
  mg_param.pre_orthonormalize = if pre_orthonormalize: QUDA_BOOLEAN_TRUE else: QUDA_BOOLEAN_FALSE
  mg_param.post_orthonormalize = if post_orthonormalize: QUDA_BOOLEAN_TRUE else: QUDA_BOOLEAN_FALSE
  mg_param.compute_null_vector = if generate_nullspace: QUDA_COMPUTE_NULL_VECTOR_YES else: QUDA_COMPUTE_NULL_VECTOR_NO
  mg_param.generate_all_levels = if generate_all_levels: QUDA_BOOLEAN_TRUE else: QUDA_BOOLEAN_FALSE
  mg_param.run_verify = if verify_results: QUDA_BOOLEAN_TRUE else: QUDA_BOOLEAN_FALSE
  mg_param.run_low_mode_check = if low_mode_check: QUDA_BOOLEAN_TRUE else: QUDA_BOOLEAN_FALSE
  mg_param.run_oblique_proj_check = if oblique_proj_check: QUDA_BOOLEAN_TRUE else: QUDA_BOOLEAN_FALSE
  mg_param.use_mma = if mg_use_mma: QUDA_BOOLEAN_TRUE else: QUDA_BOOLEAN_FALSE
  ##  Whether or not to use thin restarts in the evolve tests
  mg_param.thin_update_only = if mg_evolve_thin_updates: QUDA_BOOLEAN_TRUE else: QUDA_BOOLEAN_FALSE
  ##  set file i/o parameters
  for i in 0..<mg_param.n_level:
    if mg_vec_infile[i] != "":
      cstrcpy(addr mg_param.vec_infile[i][0], addr mg_vec_infile[i][0])
      mg_param.vec_load[i] = QUDA_BOOLEAN_TRUE
    if mg_vec_outfile[i] != "":
      cstrcpy(addr mg_param.vec_outfile[i][0], addr mg_vec_outfile[i][0])
      mg_param.vec_store[i] = QUDA_BOOLEAN_TRUE
  mg_param.coarse_guess = if mg_eig_coarse_guess: QUDA_BOOLEAN_TRUE else: QUDA_BOOLEAN_FALSE
  ##  these need to tbe set for now but are actually ignored by the MG setup
  ##  needed to make it pass the initialization test
  inv_param.inv_type = QUDA_GCR_INVERTER
  inv_param.tol = 1e-10
  inv_param.maxiter = 1000
  inv_param.reliable_delta = reliable_delta
  inv_param.gcrNkrylov = 10
  inv_param.verbosity = verbosity
  inv_param.verbosity_precondition = verbosity
