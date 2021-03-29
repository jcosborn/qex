const
  QUDA_VERSION_MAJOR* = 1
  QUDA_VERSION_MINOR* = 0
  QUDA_VERSION_SUBMINOR* = 0

## *
##  @def   QUDA_VERSION
##  @brief This macro is deprecated.  Use QUDA_VERSION_MAJOR, etc., instead.
##

const
  QUDA_VERSION* = ((QUDA_VERSION_MAJOR shl 16) or (QUDA_VERSION_MINOR shl 8) or
      QUDA_VERSION_SUBMINOR)

## *
##  @def   QUDA_MAX_DIM
##  @brief Maximum number of dimensions supported by QUDA.  In practice, no
##         routines make use of more than 5.
##

const
  QUDA_MAX_DIM* = 6

## *
##  @def   QUDA_MAX_GEOMETRY
##  @brief Maximum geometry supported by a field.  This essentially is
##  the maximum number of dimensions supported per lattice site.
##

const
  QUDA_MAX_GEOMETRY* = 8

## *
##  @def QUDA_MAX_MULTI_SHIFT
##  @brief Maximum number of shifts supported by the multi-shift solver.
##         This number may be changed if need be.
##

const
  QUDA_MAX_MULTI_SHIFT* = 32

## *
##  @def QUDA_MAX_BLOCK_SRC
##  @brief Maximum number of sources that can be supported by the block solver
##

const
  QUDA_MAX_BLOCK_SRC* = 64

## *
##  @def QUDA_MAX_ARRAY
##  @brief Maximum array length used in QudaInvertParam arrays
##

const
  QUDA_MAX_ARRAY_SIZE* = (if QUDA_MAX_MULTI_SHIFT > QUDA_MAX_BLOCK_SRC: QUDA_MAX_MULTI_SHIFT else: QUDA_MAX_BLOCK_SRC)

## *
##  @def   QUDA_MAX_DWF_LS
##  @brief Maximum length of the Ls dimension for domain-wall fermions
##

const
  QUDA_MAX_DWF_LS* = 32

## *
##  @def QUDA_MAX_MG_LEVEL
##  @brief Maximum number of multi-grid levels.  This number may be
##  increased if needed.
##

const
  QUDA_MAX_MG_LEVEL* = 4

## *
##  @def QUDA_MAX_MULTI_REDUCE
##  @brief Maximum number of simultaneous reductions that can take
##  place.  This number may be increased if needed.
##

const
  QUDA_MAX_MULTI_REDUCE* = 1024
