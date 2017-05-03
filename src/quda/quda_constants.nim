const
  QUDA_VERSION_MAJOR* = 0
  QUDA_VERSION_MINOR* = 8
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
  QUDA_MAX_DIM* = 5

## *
##  @def QUDA_MAX_MULTI_SHIFT
##  @brief Maximum number of shifts supported by the multi-shift solver.
##         This number may be changed if need be.
## 

const
  QUDA_MAX_MULTI_SHIFT* = 32

## *
##  @def   QUDA_MAX_DWF_LS
##  @brief Maximum length of the Ls dimension for domain-wall fermions
## 

const
  QUDA_MAX_DWF_LS* = 128
