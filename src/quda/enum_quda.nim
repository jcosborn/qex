const INT_MIN = -2147483648

const
  QUDA_INVALID_ENUM* = INT_MIN

type
  QudaMemoryType* {.size: sizeof(cint).} = enum
    QUDA_MEMORY_INVALID = QUDA_INVALID_ENUM,
    QUDA_MEMORY_DEVICE = 0, QUDA_MEMORY_PINNED, QUDA_MEMORY_MAPPED


##
##  Types used in QudaGaugeParam
##

type
  QudaLinkType* {.size: sizeof(cint).} = enum
    QUDA_INVALID_LINKS = QUDA_INVALID_ENUM,
    QUDA_SU3_LINKS = 0, QUDA_GENERAL_LINKS, QUDA_THREE_LINKS, QUDA_MOMENTUM_LINKS, QUDA_COARSE_LINKS, ##  used for coarse-gauge field with multigrid
    QUDA_SMEARED_LINKS       ##  used for loading and saving gaugeSmeared in the interface
  QudaGaugeFieldOrder* {.size: sizeof(cint).} = enum
    QUDA_INVALID_GAUGE_ORDER = QUDA_INVALID_ENUM,
    QUDA_FLOAT_GAUGE_ORDER = 1, QUDA_FLOAT2_GAUGE_ORDER = 2, ##  no reconstruct and double precision
    QUDA_FLOAT4_GAUGE_ORDER = 4, ##  8 and 12 reconstruct half and single
    QUDA_QDP_GAUGE_ORDER,     ##  expect *gauge[mu], even-odd, spacetime, row-column color
    QUDA_QDPJIT_GAUGE_ORDER,  ##  expect *gauge[mu], even-odd, complex-column-row-spacetime
    QUDA_CPS_WILSON_GAUGE_ORDER, ##  expect *gauge, even-odd, mu, spacetime, column-row color
    QUDA_MILC_GAUGE_ORDER,    ##  expect *gauge, even-odd, mu, spacetime, row-column order
    QUDA_MILC_SITE_GAUGE_ORDER, ##  packed into MILC site AoS [even-odd][spacetime] array, and [dir][row][col] inside
    QUDA_BQCD_GAUGE_ORDER,    ##  expect *gauge, mu, even-odd, spacetime+halos, column-row order
    QUDA_TIFR_GAUGE_ORDER,    ##  expect *gauge, mu, even-odd, spacetime, column-row order
    QUDA_TIFR_PADDED_GAUGE_ORDER ##  expect *gauge, mu, parity, t, z+halo, y, x/2, column-row order
  QudaTboundary* {.size: sizeof(cint).} = enum
    QUDA_INVALID_T_BOUNDARY = QUDA_INVALID_ENUM,
    QUDA_ANTI_PERIODIC_T = -1, QUDA_PERIODIC_T = 1
  QudaPrecision* {.size: sizeof(cint).} = enum
    QUDA_INVALID_PRECISION = QUDA_INVALID_ENUM,
    QUDA_QUARTER_PRECISION = 1, QUDA_HALF_PRECISION = 2, QUDA_SINGLE_PRECISION = 4,
    QUDA_DOUBLE_PRECISION = 8
  QudaReconstructType* {.size: sizeof(cint).} = enum
    QUDA_RECONSTRUCT_INVALID = QUDA_INVALID_ENUM,
    QUDA_RECONSTRUCT_8 = 8,     ##  reconstruct from 8 real numbers
    QUDA_RECONSTRUCT_9 = 9,     ##  used for storing HISQ long-link variables
    QUDA_RECONSTRUCT_10 = 10,   ##  10-number parameterization used for storing the momentum field
    QUDA_RECONSTRUCT_12 = 12, ##  reconstruct from 12 real numbers
    QUDA_RECONSTRUCT_13 = 13,   ##  used for storing HISQ long-link variables
    QUDA_RECONSTRUCT_NO = 18    ##  store all 18 real numbers explicitly
  QudaGaugeFixed* {.size: sizeof(cint).} = enum
    QUDA_GAUGE_FIXED_INVALID = QUDA_INVALID_ENUM,
    QUDA_GAUGE_FIXED_NO = 0,      ##  no gauge fixing
    QUDA_GAUGE_FIXED_YES     ##  gauge field stored in temporal gauge

const
  QUDA_WILSON_LINKS = QUDA_SU3_LINKS
  QUDA_ASQTAD_FAT_LINKS = QUDA_GENERAL_LINKS
  QUDA_ASQTAD_LONG_LINKS = QUDA_THREE_LINKS
  QUDA_ASQTAD_MOM_LINKS = QUDA_MOMENTUM_LINKS
  QUDA_ASQTAD_GENERAL_LINKS = QUDA_GENERAL_LINKS






##
##  Types used in QudaInvertParam
##

type
  QudaDslashType* {.size: sizeof(cint).} = enum
    QUDA_INVALID_DSLASH = QUDA_INVALID_ENUM,
    QUDA_WILSON_DSLASH = 0, QUDA_CLOVER_WILSON_DSLASH, QUDA_DOMAIN_WALL_DSLASH,
    QUDA_DOMAIN_WALL_4D_DSLASH, QUDA_MOBIUS_DWF_DSLASH, QUDA_STAGGERED_DSLASH,
    QUDA_ASQTAD_DSLASH, QUDA_TWISTED_MASS_DSLASH, QUDA_TWISTED_CLOVER_DSLASH,
    QUDA_LAPLACE_DSLASH, QUDA_COVDEV_DSLASH
  QudaInverterType* {.size: sizeof(cint).} = enum
    QUDA_INVALID_INVERTER = QUDA_INVALID_ENUM,
    QUDA_CG_INVERTER = 0, QUDA_BICGSTAB_INVERTER, QUDA_GCR_INVERTER, QUDA_MR_INVERTER,
    QUDA_MPBICGSTAB_INVERTER, QUDA_SD_INVERTER, QUDA_XSD_INVERTER,
    QUDA_PCG_INVERTER, QUDA_MPCG_INVERTER, QUDA_EIGCG_INVERTER,
    QUDA_INC_EIGCG_INVERTER, QUDA_GMRESDR_INVERTER, QUDA_GMRESDR_PROJ_INVERTER,
    QUDA_GMRESDR_SH_INVERTER, QUDA_FGMRESDR_INVERTER, QUDA_MG_INVERTER,
    QUDA_BICGSTABL_INVERTER, QUDA_CGNE_INVERTER, QUDA_CGNR_INVERTER,
    QUDA_CG3_INVERTER, QUDA_CG3NE_INVERTER, QUDA_CG3NR_INVERTER,
    QUDA_CA_CG_INVERTER, QUDA_CA_GCR_INVERTER
  QudaEigType* {.size: sizeof(cint).} = enum
    QUDA_INVALID_TYPE = QUDA_INVALID_ENUM,
    QUDA_LANCZOS = 0,             ## Normal Lanczos eigen solver
    QUDA_IMP_RST_LANCZOS     ## implicit restarted lanczos solver
  QudaSolutionType* {.size: sizeof(cint).} = enum
    QUDA_INVALID_SOLUTION = QUDA_INVALID_ENUM,
    QUDA_MAT_SOLUTION = 0, QUDA_MATDAG_MAT_SOLUTION, QUDA_MATPC_SOLUTION,
    QUDA_MATPC_DAG_SOLUTION, QUDA_MATPCDAG_MATPC_SOLUTION,
    QUDA_MATPCDAG_MATPC_SHIFT_SOLUTION
  QudaSolveType* {.size: sizeof(cint).} = enum
    QUDA_INVALID_SOLVE = QUDA_INVALID_ENUM,
    QUDA_DIRECT_SOLVE = 0, QUDA_NORMOP_SOLVE, QUDA_DIRECT_PC_SOLVE,
    QUDA_NORMOP_PC_SOLVE, QUDA_NORMERR_SOLVE, QUDA_NORMERR_PC_SOLVE
  QudaMultigridCycleType* {.size: sizeof(cint).} = enum
    QUDA_MG_CYCLE_INVALID = QUDA_INVALID_ENUM,
    QUDA_MG_CYCLE_VCYCLE = 0, QUDA_MG_CYCLE_FCYCLE, QUDA_MG_CYCLE_WCYCLE,
    QUDA_MG_CYCLE_RECURSIVE
  QudaSchwarzType* {.size: sizeof(cint).} = enum
    QUDA_INVALID_SCHWARZ = QUDA_INVALID_ENUM,
    QUDA_ADDITIVE_SCHWARZ = 0, QUDA_MULTIPLICATIVE_SCHWARZ
  QudaResidualType* {.size: sizeof(cint).} = enum
    QUDA_INVALID_RESIDUAL = QUDA_INVALID_ENUM,
    QUDA_L2_RELATIVE_RESIDUAL = 1, ##  L2 relative residual (default)
    QUDA_L2_ABSOLUTE_RESIDUAL = 2, ##  L2 absolute residual
    QUDA_HEAVY_QUARK_RESIDUAL = 4 ##  Fermilab heavy quark residual





const
  QUDA_NORMEQ_SOLVE = QUDA_NORMOP_SOLVE
  QUDA_NORMEQ_PC_SOLVE = QUDA_NORMOP_PC_SOLVE




##  Whether the preconditioned matrix is (1-k^2 Deo Doe) or (1-k^2 Doe Deo)
##
##  For the clover-improved Wilson Dirac operator, QUDA_MATPC_EVEN_EVEN
##  defaults to the "symmetric" form, (1 - k^2 A_ee^-1 D_eo A_oo^-1 D_oe),
##  and likewise for QUDA_MATPC_ODD_ODD.
##
##  For the "asymmetric" form, (A_ee - k^2 D_eo A_oo^-1 D_oe), select
##  QUDA_MATPC_EVEN_EVEN_ASYMMETRIC.
##

type
  QudaMatPCType* {.size: sizeof(cint).} = enum
    QUDA_MATPC_INVALID = QUDA_INVALID_ENUM,
    QUDA_MATPC_EVEN_EVEN = 0, QUDA_MATPC_ODD_ODD, QUDA_MATPC_EVEN_EVEN_ASYMMETRIC,
    QUDA_MATPC_ODD_ODD_ASYMMETRIC
  QudaDagType* {.size: sizeof(cint).} = enum
    QUDA_DAG_INVALID = QUDA_INVALID_ENUM,
    QUDA_DAG_NO = 0, QUDA_DAG_YES
  QudaMassNormalization* {.size: sizeof(cint).} = enum
    QUDA_INVALID_NORMALIZATION = QUDA_INVALID_ENUM,
    QUDA_KAPPA_NORMALIZATION = 0, QUDA_M_NORMALIZATION, # avoid redef in nim
    QUDA_ASYMMETRIC_MASS_NORMALIZATION
  QudaSolverNormalization* {.size: sizeof(cint).} = enum
    QUDA_DEFAULT_NORMALIZATION, ##  leave source and solution untouched
    QUDA_SOURCE_NORMALIZATION ##  normalize such that || src || = 1
  QudaPreserveSource* {.size: sizeof(cint).} = enum
    QUDA_PRESERVE_SOURCE_INVALID = QUDA_INVALID_ENUM,
    QUDA_PRESERVE_SOURCE_NO = 0,  ##  use the source for the residual
    QUDA_PRESERVE_SOURCE_YES ##  keep the source intact
  QudaDiracFieldOrder* {.size: sizeof(cint).} = enum
    QUDA_INVALID_DIRAC_ORDER = QUDA_INVALID_ENUM,
    QUDA_INTERNAL_DIRAC_ORDER = 0, ##  internal dirac order used, varies on precision and dslash type
    QUDA_DIRAC_ORDER,         ##  even-odd, color inside spin
    QUDA_QDP_DIRAC_ORDER,     ##  even-odd, spin inside color
    QUDA_QDPJIT_DIRAC_ORDER,  ##  even-odd, complex-color-spin-spacetime
    QUDA_CPS_WILSON_DIRAC_ORDER, ##  odd-even, color inside spin
    QUDA_LEX_DIRAC_ORDER,     ##  lexicographical order, color inside spin
    QUDA_TIFR_PADDED_DIRAC_ORDER ##  padded z dimension for TIFR RHMC code
  QudaCloverFieldOrder* {.size: sizeof(cint).} = enum
    QUDA_INVALID_CLOVER_ORDER = QUDA_INVALID_ENUM,
    QUDA_FLOAT_CLOVER_ORDER = 1, ##  even-odd float ordering
    QUDA_FLOAT2_CLOVER_ORDER = 2, ##  even-odd float2 ordering
    QUDA_FLOAT4_CLOVER_ORDER = 4, ##  even-odd float4 ordering
    QUDA_PACKED_CLOVER_ORDER, ##  even-odd, QDP packed
    QUDA_QDPJIT_CLOVER_ORDER, ##  (diagonal / off-diagonal)-chirality-spacetime
    QUDA_BQCD_CLOVER_ORDER   ##  even-odd, super-diagonal packed and reordered
  QudaVerbosity* {.size: sizeof(cint).} = enum
    QUDA_INVALID_VERBOSITY = QUDA_INVALID_ENUM,
    QUDA_SILENT = 0, QUDA_SUMMARIZE, QUDA_VERBOSE, QUDA_DEBUG_VERBOSE
  QudaTune* {.size: sizeof(cint).} = enum
    QUDA_TUNE_INVALID = QUDA_INVALID_ENUM,
    QUDA_TUNE_NO = 0, QUDA_TUNE_YES
  QudaPreserveDirac* {.size: sizeof(cint).} = enum
    QUDA_PRESERVE_DIRAC_INVALID = QUDA_INVALID_ENUM,
    QUDA_PRESERVE_DIRAC_NO = 0, QUDA_PRESERVE_DIRAC_YES











##
##  Type used for "parity" argument to dslashQuda()
##

type
  QudaParity* {.size: sizeof(cint).} = enum
    QUDA_INVALID_PARITY = QUDA_INVALID_ENUM,
    QUDA_EVEN_PARITY = 0, QUDA_ODD_PARITY


##
##  Types used only internally
##

type
  QudaDiracType* {.size: sizeof(cint).} = enum
    QUDA_INVALID_DIRAC = QUDA_INVALID_ENUM,
    QUDA_WILSON_DIRAC = 0, QUDA_WILSONPC_DIRAC, QUDA_CLOVER_DIRAC, QUDA_CLOVERPC_DIRAC,
    QUDA_DOMAIN_WALL_DIRAC, QUDA_DOMAIN_WALLPC_DIRAC, QUDA_DOMAIN_WALL_4DPC_DIRAC,
    QUDA_MOBIUS_DOMAIN_WALL_DIRAC, QUDA_MOBIUS_DOMAIN_WALLPC_DIRAC,
    QUDA_STAGGERED_DIRAC, QUDA_STAGGEREDPC_DIRAC, QUDA_ASQTAD_DIRAC,
    QUDA_ASQTADPC_DIRAC, QUDA_TWISTED_MASS_DIRAC, QUDA_TWISTED_MASSPC_DIRAC,
    QUDA_TWISTED_CLOVER_DIRAC, QUDA_TWISTED_CLOVERPC_DIRAC, QUDA_COARSE_DIRAC,
    QUDA_COARSEPC_DIRAC, QUDA_GAUGE_LAPLACE_DIRAC, QUDA_GAUGE_LAPLACEPC_DIRAC,
    QUDA_GAUGE_COVDEV_DIRAC


##  Where the field is stored

type
  QudaFieldLocation* {.size: sizeof(cint).} = enum
    QUDA_INVALID_FIELD_LOCATION = QUDA_INVALID_ENUM,
    QUDA_CPU_FIELD_LOCATION = 1, QUDA_CUDA_FIELD_LOCATION = 2


##  Which sites are included

type
  QudaSiteSubset* {.size: sizeof(cint).} = enum
    QUDA_INVALID_SITE_SUBSET = QUDA_INVALID_ENUM,
    QUDA_PARITY_SITE_SUBSET = 1, QUDA_FULL_SITE_SUBSET = 2


##  Site ordering (always t-z-y-x, with rightmost varying fastest)

type
  QudaSiteOrder* {.size: sizeof(cint).} = enum
    QUDA_INVALID_SITE_ORDER = QUDA_INVALID_ENUM,
    QUDA_LEXICOGRAPHIC_SITE_ORDER = 0, ##  lexicographic ordering
    QUDA_EVEN_ODD_SITE_ORDER, ##  QUDA and QDP use this
    QUDA_ODD_EVEN_SITE_ORDER ##  CPS uses this


##  Degree of freedom ordering

type
  QudaFieldOrder* {.size: sizeof(cint).} = enum
    QUDA_INVALID_FIELD_ORDER = QUDA_INVALID_ENUM,
    QUDA_FLOAT_FIELD_ORDER = 1, ##  spin-color-complex-space
    QUDA_FLOAT2_FIELD_ORDER = 2, ##  (spin-color-complex)/2-space-(spin-color-complex)%2
    QUDA_FLOAT4_FIELD_ORDER = 4, ##  (spin-color-complex)/4-space-(spin-color-complex)%4
    QUDA_SPACE_SPIN_COLOR_FIELD_ORDER, ##  CPS/QDP++ ordering
    QUDA_SPACE_COLOR_SPIN_FIELD_ORDER, ##  QLA ordering (spin inside color)
    QUDA_QDPJIT_FIELD_ORDER,  ##  QDP field ordering (complex-color-spin-spacetime)
    QUDA_QOP_DOMAIN_WALL_FIELD_ORDER, ##  QOP domain-wall ordering
    QUDA_PADDED_SPACE_SPIN_COLOR_FIELD_ORDER ##  TIFR RHMC ordering
  QudaFieldCreate* {.size: sizeof(cint).} = enum
    QUDA_INVALID_FIELD_CREATE = QUDA_INVALID_ENUM,
    QUDA_NULL_FIELD_CREATE = 0,   ##  create new field
    QUDA_ZERO_FIELD_CREATE,   ##  create new field and zero it
    QUDA_COPY_FIELD_CREATE,   ##  create copy to field
    QUDA_REFERENCE_FIELD_CREATE ##  create reference to field
  QudaGammaBasis* {.size: sizeof(cint).} = enum
    QUDA_INVALID_GAMMA_BASIS = QUDA_INVALID_ENUM,
    QUDA_DEGRAND_ROSSI_GAMMA_BASIS = 0, QUDA_UKQCD_GAMMA_BASIS,
    QUDA_CHIRAL_GAMMA_BASIS
  QudaSourceType* {.size: sizeof(cint).} = enum
    QUDA_INVALID_SOURCE = QUDA_INVALID_ENUM,
    QUDA_POINT_SOURCE = 0, QUDA_RANDOM_SOURCE, QUDA_CONSTANT_SOURCE,
    QUDA_SINUSOIDAL_SOURCE, QUDA_CORNER_SOURCE
  QudaNoiseType* {.size: sizeof(cint).} = enum
    QUDA_NOISE_INVALID = QUDA_INVALID_ENUM,
    QUDA_NOISE_GAUSS = 0, QUDA_NOISE_UNIFORM






##  used to select projection method for deflated solvers

type
  QudaProjectionType* {.size: sizeof(cint).} = enum
    QUDA_INVALID_PROJECTION = QUDA_INVALID_ENUM,
    QUDA_MINRES_PROJECTION = 0, QUDA_GALERKIN_PROJECTION


##  used to select preconditioning method in domain-wall fermion

type
  QudaDWFPCType* {.size: sizeof(cint).} = enum
    QUDA_PC_INVALID = QUDA_INVALID_ENUM,
    QUDA_5D_PC = 0, QUDA_4D_PC
  QudaTwistFlavorType* {.size: sizeof(cint).} = enum
    QUDA_TWIST_INVALID = QUDA_INVALID_ENUM,
    QUDA_TWIST_DEG_DOUBLET = -2, QUDA_TWIST_NO = 0,
    QUDA_TWIST_SINGLET = 1,
    QUDA_TWIST_NONDEG_DOUBLET = +2
  QudaTwistDslashType* {.size: sizeof(cint).} = enum
    QUDA_DSLASH_INVALID = QUDA_INVALID_ENUM,
    QUDA_DEG_TWIST_INV_DSLASH = 0, QUDA_DEG_DSLASH_TWIST_INV,
    QUDA_DEG_DSLASH_TWIST_XPAY, QUDA_NONDEG_DSLASH
  QudaTwistCloverDslashType* {.size: sizeof(cint).} = enum
    QUDA_TC_DSLASH_INVALID = QUDA_INVALID_ENUM,
    QUDA_DEG_CLOVER_TWIST_INV_DSLASH = 0, QUDA_DEG_DSLASH_CLOVER_TWIST_INV,
    QUDA_DEG_DSLASH_CLOVER_TWIST_XPAY
  QudaTwistGamma5Type* {.size: sizeof(cint).} = enum
    QUDA_TWIST_GAMMA5_INVALID = QUDA_INVALID_ENUM,
    QUDA_TWIST_GAMMA5_DIRECT = 0, QUDA_TWIST_GAMMA5_INVERSE
  QudaUseInitGuess* {.size: sizeof(cint).} = enum
    QUDA_USE_INIT_GUESS_INVALID = QUDA_INVALID_ENUM,
    QUDA_USE_INIT_GUESS_NO = 0, QUDA_USE_INIT_GUESS_YES
  QudaComputeNullVector* {.size: sizeof(cint).} = enum
    QUDA_COMPUTE_NULL_VECTOR_INVALID = QUDA_INVALID_ENUM,
    QUDA_COMPUTE_NULL_VECTOR_NO = 0, QUDA_COMPUTE_NULL_VECTOR_YES
  QudaSetupType* {.size: sizeof(cint).} = enum
    QUDA_INVALID_SETUP_TYPE = QUDA_INVALID_ENUM,
    QUDA_NULL_VECTOR_SETUP = 0, QUDA_TEST_VECTOR_SETUP
  QudaBoolean* {.size: sizeof(cint).} = enum
    QUDA_BOOLEAN_INVALID = QUDA_INVALID_ENUM,
    QUDA_BOOLEAN_NO = 0, QUDA_BOOLEAN_YES = 1
  QudaDirection* {.size: sizeof(cint).} = enum
    QUDA_BACKWARDS = -1, QUDA_FORWARDS = +1, QUDA_BOTH_DIRS = 2
  QudaLinkDirection* {.size: sizeof(cint).} = enum
    QUDA_LINK_BACKWARDS, QUDA_LINK_FORWARDS, QUDA_LINK_BIDIRECTIONAL
  QudaFieldGeometry* {.size: sizeof(cint).} = enum
    QUDA_INVALID_GEOMETRY = QUDA_INVALID_ENUM,
    QUDA_SCALAR_GEOMETRY = 1, QUDA_VECTOR_GEOMETRY = 4, QUDA_TENSOR_GEOMETRY = 6,
    QUDA_COARSE_GEOMETRY = 8
  QudaGhostExchange* {.size: sizeof(cint).} = enum
    QUDA_GHOST_EXCHANGE_INVALID = QUDA_INVALID_ENUM,
    QUDA_GHOST_EXCHANGE_NO = 0, QUDA_GHOST_EXCHANGE_PAD, QUDA_GHOST_EXCHANGE_EXTENDED
  QudaStaggeredPhase* {.size: sizeof(cint).} = enum
    QUDA_STAGGERED_PHASE_INVALID = QUDA_INVALID_ENUM,
    QUDA_STAGGERED_PHASE_NO = 0, QUDA_STAGGERED_PHASE_MILC = 1,
    QUDA_STAGGERED_PHASE_CPS = 2, QUDA_STAGGERED_PHASE_TIFR = 3
  QudaContractType* {.size: sizeof(cint).} = enum
    QUDA_CONTRACT_INVALID = QUDA_INVALID_ENUM,
    QUDA_CONTRACT = 0, QUDA_CONTRACT_PLUS, QUDA_CONTRACT_MINUS, QUDA_CONTRACT_GAMMA5,
    QUDA_CONTRACT_GAMMA5_PLUS, QUDA_CONTRACT_GAMMA5_MINUS, QUDA_CONTRACT_TSLICE,
    QUDA_CONTRACT_TSLICE_PLUS, QUDA_CONTRACT_TSLICE_MINUS
















## Allows to choose an appropriate external library

type
  QudaExtLibType* {.size: sizeof(cint).} = enum
    QUDA_EXTLIB_INVALID = QUDA_INVALID_ENUM,
    QUDA_CUSOLVE_EXTLIB = 0, QUDA_EIGEN_EXTLIB, QUDA_MAGMA_EXTLIB

