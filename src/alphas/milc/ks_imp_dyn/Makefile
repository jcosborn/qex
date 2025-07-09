#  MIMD version 7
#  Standard application makefile.
#  Replaces Make_vanilla, Make_linux_mpi, Make_linux_qdp, Make_qcdoc_gcc
#  Do not use for making libraries
#  Copy this file into the application directory and run make there.
#

# This file
MAKEFILE = Makefile

#----------------------------------------------------------------------
#  User choices - edit to suit 
#----------------------------------------------------------------------
# 1. Host and accelerator architecture.  Controls optimization flags here and in libraries.
#    Can control BINEXT below, a suffix appended to the name of the executable.

ARCH ?= # epyc hsw skx clx icx spr knl pow8 pow9
#GPU_ARCH ?= # nvidia amd intel

#----------------------------------------------------------------------
# 2. Compiler family

COMPILER ?= intel # intel, ibm, cray-intel, rocm
OFFLOAD ?= # cuda hip sycl openmp

#----------------------------------------------------------------------
# 3. MPP vs Scalar

# Compiling with MPI?  false for a scalar machine
MPP = true

#----------------------------------------------------------------------
# 4. Generic Precision 

# 1 = single precision; 2 = double
PRECISION = 2

#----------------------------------------------------------------------
# 5. Set compiler.
# Choices include mpicc cc gcc pgcc g++

ifeq ($(strip ${COMPILER}),intel)

  ifeq ($(strip ${MPP}),true)
    MY_CC ?= mpicc
    MY_CXX ?= mpicxx
  else
    MY_CC  ?= icx
    MY_CXX ?= icpx
  endif

else ifeq ($(strip ${COMPILER}),cray-intel)

  ifeq ($(strip ${MPP}),true)
    MY_CC ?= cc
    MY_CXX ?= CC
  else
    MY_CC  ?= icc
    MY_CXX ?= icpc
  endif

else ifeq ($(strip ${COMPILER}),gnu)

  ifeq ($(strip ${MPP}),true)
    MY_CC ?= mpicc
    MY_CXX ?= mpiCC
  else
    MY_CC  ?= gcc-8
    MY_CXX ?= g++-8
  endif

else ifeq ($(strip ${COMPILER}),ibm)

  ifeq ($(strip ${MPP}),true)
    MY_CC ?= mpixlc_r
    MY_CXX ?= mpixlcxx_r
  else
    MY_CC ?= xlc_r
    MY_CXX ?= xlc++_r
  endif

else ifeq ($(strip ${COMPILER}),amdclang)

  ifeq ($(strip ${MPP}),true)
    MY_CC ?= mpicc
    MY_CXX ?= mpicxx
  else
    MY_CC ?= amdclang
    MY_CXX ?= amdclang++
  endif

endif

# Accelerator

ifeq ($(strip ${OFFLOAD}),sycl)

  MY_CC += -cc=icx
  MY_CXX += -cxx=icpx -fsycl

endif

CC = ${MY_CC}
CXX = ${MY_CXX}

# If the above construction doesn't work, override the definitions here

# CC =
# CXX =

#----------------------------------------------------------------------
# 6. Compiler optimization level
# Choices include -g -O, etc
# Power9 recommendations are -Ofast

OPT              ?= -O3 -g

# OpenMP?

OMP = true #true

#----------------------------------------------------------------------
# 7. Other compiler optimization flags.  Uncomment stanza to suit.

#-------------- Gnu C -------------------------------------

ifeq ($(strip ${COMPILER}),gnu)

  OCFLAGS += -std=c99
  OCXXFLAGS += -std=gnu++17

  ifeq ($(strip ${ARCH}),pow8)
    ARCH_FLAG = -mcpu=power8
  endif

  ifeq ($(strip ${ARCH}),pow9)
    ARCH_FLAG = -mcpu=power9 -mtune=power9
  endif

  ifeq ($(strip ${OMP}),true)
    OCFLAGS += -fopenmp
    OCXXFLAGS += -fopenmp
    LDFLAGS += -fopenmp
  endif

# Other Gnu options
#OCFLAGS += -mavx # depends on architecture
# enable all warnings with exceptions
OCFLAGS += -Wall -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unknown-pragmas -Wno-unused-function

endif

#------------------------ BlueGene -----------------------------------

ifeq ($(strip ${COMPILER}),ibm)

  OCFLAGS = -std=gnu99
  OCXXFLAGS += -std=c++17

  # consider running with XLSMPOPTS=stack=10M or higher
  ifeq ($(strip ${OMP}),true)
    OCFLAGS += -qsmp=omp
    OCXXFLAGS += -qsmp=omp
    LDFLAGS += -qsmp=omp
  endif

endif

#-------------- Intel icc/ecc -----------------------------------

ifeq ($(strip ${COMPILER}),intel)

  OCFLAGS += -std=c99
  OCXXFLAGS += -std=c++17

  ifeq ($(strip ${ARCH}),knl)
  ARCH_FLAG = -xMIC-AVX512
  BINEXT=.knl
  else ifeq ($(strip ${ARCH}),knc)
  ARCH_FLAG = -mmic
  BINEXT=.knc
  else ifeq ($(strip ${ARCH}),skx)
  ARCH_FLAG = -xCORE-AVX512 -qopt-zmm-usage=high
  BINEXT=.skx
  else ifeq ($(strip ${ARCH}),clx)
  ARCH_FLAG = -xCORE-AVX512 -qopt-zmm-usage=high
  else ifeq ($(strip ${ARCH}),spr)
  ARCH_FLAG = ""
  else ifeq ($(strip ${ARCH}),icx)
  ARCH_FLAG = -xCORE-AVX512 -qopt-zmm-usage=high
  else ifeq ($(strip ${ARCH}),hsw)
  ARCH_FLAG = -xCORE-AVX2
  BINEXT=.hsw
  else
  ARCH_FLAG = -mavx
  BINEXT=
  endif

  OCFLAGS += ${ARCH_FLAG}
  OCXXFLAGS += ${ARCH_FLAG}
  LDFLAGS += ${ARCH_FLAG}
  OCFLAGS += -debug
  OCXXFLAGS += -debug

  ifeq ($(strip ${OMP}),true)
    OCFLAGS += -fopenmp
    OCXXFLAGS += -fopenmp
    LDFLAGS += -fopenmp
  endif

endif

#-------------- Cray-Intel  -----------------------------------

ifeq ($(strip ${COMPILER}),cray-intel)

  OCFLAGS += -std=c99
  OCXXFLAGS += -std=c++17

  ifeq ($(strip ${ARCH}),knl)
  ARCH_FLAG = -xMIC-AVX512
  BINEXT=.knl
  else ifeq ($(strip ${ARCH}),knc)
  ARCH_FLAG = -mmic
  BINEXT=.knc
  else ifeq ($(strip ${ARCH}),hsw)
  ARCH_FLAG = -xCORE-AVX2
  BINEXT=.hsw
  else
  ARCH_FLAG = -mavx
  BINEXT=
  endif

  OCFLAGS += ${ARCH_FLAG}
  OCXXFLAGS += ${ARCH_FLAG}
  LDFLAGS += ${ARCH_FLAG}
  OCFLAGS += # -parallel-source-info=2 -debug inline-debug-info -qopt-report=5
  OCXXFLAGS += # -parallel-source-info=2 -debug inline-debug-info -qopt-report=5

  ifeq ($(strip ${OMP}),true)
    OCFLAGS += -qopenmp
    OCXXFLAGS += -qopenmp
    LDFLAGS += -qopenmp
  endif

endif

#-------------- Portland Group ----------------------------
#OCFLAGS = -tp p6 -Munroll=c:4,n:4
#OCFLAGS= -mpentiumpro -march=pentiumpro -funroll-all-loops -malign-double -D_REENTRANT  # Pentium pro

#-------------- AMD Clang ----------------------------
ifeq ($(strip ${COMPILER}),amdclang)

  ifeq ($(strip ${OMP}),true)
    OCFLAGS += -fopenmp
    OCXXFLAGS += -fopenmp
    LDFLAGS += -fopenmp
  endif

endif

#----------------------------------------------------------------------
# 8. Choose large file support.

ifeq ($(strip ${COMPILER}),ibm)
  CLFS = -D_LARGE_FILES   # AIX
else
  CLFS = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE # Large files gcc only
endif

#----------------------------------------------------------------------
# 9. Installation-specific MPI includes and libraries
#    Not needed if using mpicc or on single processor

#----------------- Intel MPI Custom -------------------------------------------

# IMPI = -I/opt/intel/impi/5.1.3.181/compilers_and_libraries_2016.2.181/linux/mpi/intel64/include

# ifeq ($(strip ${ARCH}),knl)
#   LMPI = -L/opt/intel/compiler/latest/compilers_and_libraries_2016/linux/mpi/mic/lib
# else ifeq ($(strip ${ARCH}),knc) 
#   LMPI = -L/opt/intel/compiler/latest/compilers_and_libraries_2016/linux/mpi/mic/lib
# else
#   LMPI = -L/opt/intel/compiler/latest/compilers_and_libraries_2016/linux/mpi/intel64/lib
# endif

#----------------- MVICH ----------------------------------------------
#IMPI = -I/uufs/icebox/sys/src/mpich/1.2.0-via/include  # MVICH
#LMPI = -L/uufs/icebox/sys/pkg/mpich/1.2.0-via/lib -lmpi -lvipl -lpthread # MVICH

#----------------------------------------------------------------------
# 10. I/O routines
# Both io_nonansi and io_ansi should work on a scalar machine
# Solaris 2.6 gave "bad file number" errors with io_ansi.  CD

MACHINE_DEP_IO   = io_ansi.o # (io_ansi.o io_nonansi.o io_dcap.o)

# Uncomment if you have and want to support dcache I/O
# (Forces use of io_dcap.o)

# WANTDCAP = true

# The location of the installed dcap libraries. Uncomment and enter if
# it is not already defined as a system environment variable.

# DCAP_DIR = /usr/local/develop/dcache

# Choose the appropriate library path according to the addressing
# size of the machine

# DCAPLIB  = lib64 # (lib64 lib)

#----------------------------------------------------------------------
# 11. SciDAC package options

# Edit these "wants"

WANTQOP ?= # true # or blank. Implies HAVEQDP, HAVEQOP, HAVEQMP.

WANTQIO ?= # true # or blank.  Implies HAVEQMP.

WANTQMP ?= # true or blank.

# QMP_MPI or QMP_SPI
QMP_BACKEND = QMP_MPI

# Edit these locations for the installed SciDAC packages
# It is assumed that these are the parents of "include" and "lib"

SCIDAC = ${HOME}/scidac/install
TAG=
# Parallel versions
QMPPAR ?= ${SCIDAC}/qmp${TAG}
QIOPAR ?= $(SCIDAC)/qio${TAG}
# Single processor versions
QMPSNG ?= ${SCIDAC}/qmp-single${TAG}
QIOSNG ?= $(SCIDAC)/qio-single${TAG}
QLA ?= ${SCIDAC}/qla${TAG}
# Either version
QDP ?= ${SCIDAC}/qdp${TAG}
QOPQDP ?= ${SCIDAC}/qopqdp${TAG}

QOP = ${QOPQDP}

# Make_template_scidac defines these macros:
# HAVEQOP HAVEQDP HAVEQIO HAVEQMP (Says what we are compiling with)
# LIBSCIDAC INCSCIDAC (The -L and -I compiler and linker lists)
# SCIDAC_LIBRARIES SCIDAC_HEADERS  (Lists for make dependencies)
# CSCIDAC (List of compiler macros for SciDAC modules)

include ../Make_template_scidac

#----------------------------------------------------------------------
# 12. APE link I/O
WANT_APE_IO ?= # true

#----------------------------------------------------------------------
# 12. Intel MKL for FFTW and LAPACK

ifeq ($(strip ${COMPILER}),intel)
  INCFFTW = -mkl
  LIBFFTW = -mkl
else ifeq ($(strip ${COMPILER}),cray-intel)
  INCFFTW = -mkl
  LIBFFTW = -mkl
endif

#----------------------------------------------------------------------
# 12. FFTW3 Options

WANTFFTW ?= #true

ifeq ($(strip ${WANTFFTW}),true)
  FFTW ?= ${HOME}/fftw/build-gcc

  FFTW_HEADERS = ${FFTW}/include
  INCFFTW = -I${FFTW_HEADERS}
  LIBFFTW = -L${FFTW}/lib
  LIBFFTW += -lfftw3 -lfftw3f
  PACKAGE_HEADERS += ${FFTW_HEADERS}
endif

#----------------------------------------------------------------------
# 13. LAPACK Options (for qopqdp-lapack and arb_overlap )

#LIBLAPACK = -L/opt/ibmcmp/xlf/bg/11.1/lib /soft/apps/LAPACK/liblapack_bgp.a /soft/apps/LIBGOTO/libgoto.a -lxlf90 -lxlsmp # LAPACK on BG/P

# Utah physics and math Redhat-linux
# LIBLAPACK = -L/usr/uumath/lib64  -llapack-3.6.0 -lblas-3.6.0 -L/usr/lib/gcc/x86_64-redhat-linux/4.8.2 -lgfortran

# Utah physics and math Centos-linux.  Must link with gfortran. Incompatible with Grid!
# LIBLAPACK = -L/usr/local/lib64 -llapack -lblas
# LDLAPACK = gfortran

# FNAL cluster (Jim's installation of ATLAS)
# LDFLAGS = -Wl,-rpath,"/usr/local/atlas-3.10-lapack-3.4.2/lib" -L/usr/local/atlas-3.10-lapack-3.4.2/lib
# LIBS = $(LDFLAGS) -lprimme -lm  -llapack -lptf77blas -lptcblas -latlas -lgfortran -lpthread

# NERSC Cori Haswell
# LIBLAPACK = -L${LIBSCI_BASE_DIR}/INTEL/15.0/haswell/lib -lsci_intel

# NERSC Cori KNL
# LIBLAPACK = -L${LIBSCI_BASE_DIR}/INTEL/15.0/mic_knl/lib -lsci_intel

# NERSC Edison
# LIBLAPACK = -L${LIBSCI_BASE_DIR}/INTEL/15.0/ivybridge/lib -lsci_intel

# falco Centos-linux.  Must link with gfortran. 
# LIBLAPACK = -L/usr/uumath/lib64 -llapack -lblas
# LDLAPACK = gfortran

#----------------------------------------------------------------------
# 14. PRIMME Options (for arb_overlap and ks_eigen).  REQUIRES LAPACK AS WELL.

WANTPRIMME = #true

# PRIMME version 2.0

ifeq ($(strip ${WANTPRIMME}),true)
  PRIMME_HEADERS = ${HOME}/PRIMME/include
  INCPRIMME = -I${PRIMME_HEADERS}
  PACKAGE_HEADERS += ${PRIMME_HEADERS}
  LIBPRIMME = -L${HOME}/PRIMME/lib -lprimme

  CEIG ?= # -DPRIMME_PRECOND -DPOLY_EIGEN

endif

#----------------------------------------------------------------------
# 14. ARPACK Options (for ks_eigen).  REQUIRES LAPACK AS WELL.

WANTARPACK = #true

ifeq ($(strip ${WANTARPACK}),true)
#  LIBARPACK = -L/usr/lib64 -lparpack  -larpack -lifcore -llapack -lblas
  LIBARPACK = -L/usr/lib64 -larpack
endif

#----------------------------------------------------------------------
# 15. GPU/QUDA Options

WANTQUDA    ?= false

ifeq ($(strip ${WANTQUDA}),true)

  HAVE_QUDA = true
  HAVE_GPU = true
  CGPU += -DHAVE_QUDA

  WANT_CL_BCG_GPU ?= #true
  WANT_FN_CG_GPU ?= #true
  WANT_FL_GPU ?= #true
  WANT_FF_GPU ?= #true
  WANT_GF_GPU ?= #true
  WANT_EIG_GPU ?= #true
  WANT_GSMEAR_GPU ?= #true
  WANT_KS_CONT_GPU ?= #true
  WANT_SHIFT_GPU ?= #true
  WANT_SPIN_TASTE_GPU ?= #true
  WANT_GAUGEFIX_OVR_GPU ?= #true
  WANT_MULTIGRID ?= #true

endif

ifeq ($(strip ${WANTQUDA}),true)
  ifeq ($(strip ${OFFLOAD}),)
    OFFLOAD = CUDA
  endif

  QUDA_HOME ?= ${HOME}/quda

  INCQUDA = -I${QUDA_HOME}/include -I${QUDA_HOME}/tests
  PACKAGE_HEADERS += ${QUDA_HOME}/include
  LIBQUDA ?= -Wl,-rpath ${QUDA_HOME}/lib -L${QUDA_HOME}/lib -lquda
  QUDA_LIBRARIES = ${QUDA_HOME}/lib
  QUDA_HEADERS = ${QUDA_HOME}/include

  ifeq ($(strip ${OFFLOAD}),CUDA)
    CUDA_HOME ?= /usr/local/cuda
    INCQUDA += -I${CUDA_HOME}/include
    PACKAGE_HEADERS += ${CUDA_HOME}/include
    LIBQUDA += -L${CUDA_HOME}/lib64 -L${CUDA_MATH}/lib64 -L${CUDA_COMP}/lib -lcudart -lcuda -lcublas -lcufft -ldl
  endif

# Verbosity choices:
# SET_QUDA_SILENT, SET_QUDA_SUMMARIZE, SET_QUDA_VERBOSE, SET_QUDA_DEBUG_VERBOSE

  QUDA_VERBOSITY ?= SUMMARIZE
  ifeq ($(strip ${QUDA_VERBOSITY}),SILENT)
    CGPU += -DSET_QUDA_SILENT # silent output
  else ifeq ($(strip ${QUDA_VERBOSITY}),SUMMARIZE)
    CGPU += -DSET_QUDA_SUMMARIZE # summary output
  else ifeq ($(strip ${QUDA_VERBOSITY}),VERBOSE)
    CGPU += -DSET_QUDA_VERBOSE # verbose output, outputs autotuning information, residual history of solvers
  else ifeq ($(strip ${QUDA_VERBOSITY}),DEBUG_VERBOSE)
    CGPU += -DSET_QUDA_DEBUG_VERBOSE # debug-level output
  endif

endif

#----------------------------------------------------------------------
# 16. QPhiX Options

WANTQPHIX = #true
WANT_FN_CG_QPHIX = true
WANT_GF_QPHIX = true

QPHIX_HOME = ../QPhiX_MILC/milc-qphix

ifeq ($(strip ${WANTQPHIX}), true)

  INCQPHIX = -I${QPHIX_HOME}
  PACKAGE_HEADERS += ${QPHIX_HOME}
  HAVE_QPHIX = true
  CPHI = -DHAVE_QPHIX
  QPHIX_HEADERS = ${QPHIX_HOME}

  ifeq ($(strip ${MPP}),true)

  # MPI versions of QPHIX

  ifeq ($(strip ${ARCH}),knl)
    LIBQPHIX = -L${QPHIX_HOME} -lqphixmilc_avx512 -lrt
  else ifeq ($(strip ${ARCH}),knc)
    LIBQPHIX = -L${QPHIX_HOME} -lqphixmilc_mic -lrt
  else ifeq ($(strip ${ARCH}),hsw)
    LIBQPHIX = -L${QPHIX_HOME} -lqphixmilc_avx2 -lrt
  else ifeq ($(strip ${ARCH}),skx)
    LIBQPHIX = -L${QPHIX_HOME} -lqphixmilc_skx -lrt
  endif

  else

  # Non-MPI versions

  ifeq ($(strip ${ARCH}),knl)
    LIBQPHIX = -L${QPHIX_HOME} -lqphixmilc_avx512_single -lrt
  else ifeq ($(strip ${ARCH}),knc)
    LIBQPHIX = -L${QPHIX_HOME} -lqphixmilc_mic_single -lrt
  else ifeq ($(strip ${ARCH}),hsw)
    LIBQPHIX = -L${QPHIX_HOME} -lqphixmilc_avx2_single -lrt
  else ifeq ($(strip ${ARCH}),skx)
    LIBQPHIX = -L${QPHIX_HOME} -lqphixmilc_skx_single -lrt
  endif

  endif

  QPHIX_HEADERS   = ${QPHIX_HOME}
  PACKAGE_HEADERS += ${QPHIX_HEADERS}
  QPHIX_LIBRARIES = ${QPHIX_HOME}

  ifeq ($(strip ${WANT_FN_CG_QPHIX}),true)
    HAVE_FN_CG_QPHIX = true
    CPHI += -DUSE_CG_QPHIX
  endif

  ifeq ($(strip ${WANT_GF_QPHIX}),true)
    HAVE_GF_QPHIX = true
    CPHI += -DUSE_GF_QPHIX
  endif

endif

#----------------------------------------------------------------------
# 16. Hadrons Options

WANTHADRONS ?= false # true implies WANTGRID = true

ifeq ($(strip ${WANTHADRONS}), true)

  HAVE_HADRONS = true
  CGPU += -DHAVE_HADRONS

  ifeq ($(strip ${MPP}),true)
    ifeq ($(strip ${ARCH}),knl)
      HADRONS_ARCH = avx512
    else ifeq ($(strip ${ARCH}),skx)
      HADRONS_ARCH = avx512
    else ifeq ($(strip ${ARCH}),hsw)
      HADRONS_ARCH = avx2
    endif
  else
    # Scalar version                                                                

    HADRONS_ARCH = scalar

  endif

  HADRONS_HOME = ../Grid/install-hadrons-${HADRONS_ARCH}
  HADRONS_LIBRARIES = ${HADRONS_HOME}/lib
  LIBHADRONS = -L${HADRONS_LIBRARIES} -lHadrons -ldl
  HADRONS_HEADERS = ${HADRONS_HOME}/include
  INCHADRONS = -I${HADRONS_HEADERS}

  PACKAGE_HEADERS += ${HADRONS_HEADERS}/Hadrons
  PACKAGE_DEPS += Hadrons

  LDFLAGS += -fopenmp

endif

#----------------------------------------------------------------------
# 16. Grid Options

WANTGRID ?= false

ifeq ($(strip ${WANTHADRONS}), true)
  WANTGRID = true
endif

ifeq ($(strip ${WANTGRID}), true)

  HAVE_GRID = true
  HAVE_GPU = true
  CGPU += -DHAVE_GRID

  WANT_FN_CG_GPU ?= false   
  WANT_FL_GPU ?= false       # Under development
  WANT_FF_GPU ?= false       # Future
  WANT_GF_GPU ?= false       # Future
  WANT_EIG_GPU ?= false      # Automatic for now

  GRID_SHMEM_MAX ?= 2048        # Megabytes
  GRID_DEVICE_MEM_MAX ?= 32768  # Megabytes
  GRID_SHMEM_MPI ?= 1
  GRID_ACCELERATOR_THREADS ?= 8
  GRID_MULTI_CG  ?= GRID_5DCG # GRID_5DCG GRID_BLOCKCG GRID_MRHSCG

  CPHI += -DGRID_SHMEM_MAX=${GRID_SHMEM_MAX}
  CPHI += -DGRID_SHMEM_MPI=${GRID_SHMEM_MPI}
  CPHI += -DGRID_DEVICE_MEM_MAX=${GRID_DEVICE_MEM_MAX}
  CPHI += -DGRID_ACCELERATOR_THREADS=${GRID_ACCELERATOR_THREADS}
  CPHI += -DGRID_COMMS_OVERLAP=${GRID_COMMS_OVERLAP}
  CPHI += -DGRID_MULTI_CG=${GRID_MULTI_CG}

endif

ifeq ($(strip ${WANTGRID}), true)

  # Accelerator offloads
  ifeq ($(strip ${OFFLOAD}),sycl)
    GRID_ARCH = gpu-sycl
  else ifeq ($(strip ${OFFLOAD}),cuda)
    GRID_ARCH = gpu-cuda
  else ifeq ($(strip ${OFFLOAD}),hip)
    GRID_ARCH = gpu-hip
  else
    # CPU only
    ifeq ($(strip ${ARCH}),knl)
      GRID_ARCH = avx512
    else ifeq ($(strip ${ARCH}),skx)
      GRID_ARCH = avx512
    else ifeq ($(strip ${ARCH}),clx)
      GRID_ARCH = avx512
    else ifeq ($(strip ${ARCH}),icx)
      GRID_ARCH = avx512
    else ifeq ($(strip ${ARCH}),hsw)
      GRID_ARCH = avx2
    else ifeq ($(strip ${ARCH}),epyc)
      GRID_ARCH = avx2
    else
      # Scalar version                                                                
      GRID_ARCH = scalar
    endif
  endif

  GRID_HOME = ../Grid/install-grid-${GRID_ARCH}
  GRID_LIBRARIES = ${GRID_HOME}/lib
  LIBGRID = -L${GRID_LIBRARIES} -lGrid -lcrypto -lz
  GRID_HEADERS = ${GRID_HOME}/include
  INCGRID = -I${GRID_HEADERS}

  PACKAGE_HEADERS += ${GRID_HEADERS}/Grid
  PACKAGE_DEPS += Grid

endif

#----------------------------------------------------------------------
# 16. QPhiXJ (JLab) Options

WANTQPHIXJ = #true

# Choose vectorization parameters.
# Choices 4, 8 (or 1 for scalar)
QPHIXJ_SOALEN=4

ifeq ($(strip ${WANTQPHIXJ}), true)

  HAVE_QPHIXJ = true
  CPHI += -DHAVE_QPHIXJ

  ifeq ($(strip ${MPP}),true)

    # QMP versions of QPHIXJ

    QPHIXJ_HOME = ../QPhiX_JLab/install

    ifeq ($(strip ${ARCH}),knl)
      QPHIXJ_ARCH = avx512
    else ifeq ($(strip ${ARCH}),hsw)
      QPHIXJ_ARCH = avx2
    endif
  else
    # Scalar version
    QPHIXJ_ARCH = scalar
    QPHIXJ_SOALEN = 1
  endif

  # NOTE: These are QMP versions of QPHIXJ so requires QMP

  QPHIXJ_HOME = ../QPhiX_JLab/install/dslash-${QPHIXJ_ARCH}-s${QPHIXJ_SOALEN}
  QPHIXJ_LIBRARIES = ${QPHIXJ_HOME}/lib
  LIBQPHIXJ = -L${QPHIXJ_LIBRARIES} -lqphix_solver -lqphix_codegen
  QPHIXJ_HEADERS = ${QPHIXJ_HOME}/include
  INCQPHIXJ = -I${QPHIXJ_HEADERS}

  PACKAGE_HEADERS += ${QPHIXJ_HEADERS}
  PACKAGE_DEPS += QPhiX_JLab

endif

#----------------------------------------------------------------------
# 17. Linker (need the C++ linker for QUDA, QPHIX, GRID)

ifeq ($(strip ${LDLAPACK}),)
  ifeq ($(strip ${WANTQUDA}),true)
    LD  = ${CXX}
  else ifeq ($(strip ${WANTQPHIX}),true)
    LD  = ${CXX}
  else ifeq ($(strip ${WANTQPHIXJ}),true)
    LD  = ${CXX}
  else ifeq ($(strip ${WANTGRID}),true)
    LD  = ${CXX}
  else
    LD  = ${CC}
  endif
else
  LD = ${LDLAPACK}
endif

#----------------------------------------------------------------------
# 18. Extra ld flags

# VTune

# VTUNE_VERSION = 2016u2
# INCVTUNE = -I/opt/intel/vtune/${VTUNE_VERSION}/vtune_amplifier_xe_2016/include
# LIBVTUNE = -L/opt/intel/vtune/${VTUNE_VERSION}/vtune_amplifier_xe_2016/lib64 -l ittnotify
# 
# OCFLAGS += -DVTUNE

#----------------------------------------------------------------------
# 19. Inlining choices

# USE INLINE SSE WITH EXTREME CAUTION!  IT MAY GIVE WRONG RESULTS.

# SSE ASM and explicit C inlining is available for some of the library
# functions.

# Use SSE for P3 or P4 or Athlon Thunderbird and compilers, such
# as gcc that accept ASM macros

# Both SSE and C inline macros can be invoked selectively by defining
# SSE_INLINE and C_INLINE and changing the function call to the macro
# name.

# To invoke the inline macros globally (where available) without
# changing the code, define, instead, SSE_GLOBAL_INLINE or
# C_GLOBAL_INLINE.

# You may use SSE and C inlining at the same time.  The SSE version
# takes precedence when both are available.

# See also the libraries Make_SSE_nasm for building non-inline SSE
# Some compilers don't like -DSSE_INLINE with the debugging -g option.

# Choose nothing or
#  [ C_INLINE | C_GLOBAL_INLINE ] [ SSE_INLINE | SSE_GLOBAL_INLINE ]
INLINEOPT = -DC_GLOBAL_INLINE # -DSSE_GLOBAL_INLINE #-DC_INLINE

# There are special single-precision macros for the AMD Opteron
# To get them, uncomment the next line

#INLINEOPT += -DSSEOPTERON

#----------------------------------------------------------------------
# 20. Miscellaneous macros for performance control and metric

#     Define them with a -D prefix.

#------------------------------
# git code version

GIT_VERSION := "$(shell git describe --abbrev=4 --dirty --always --tags)"
CGITVER = -DMILC_CODE_VERSION=\"$(GIT_VERSION)\"

#------------------------------
# Print timing statistics.
# Applications: many

# Use any combination of these
# CGTIME CG Solver
# FFTIME Fermion force
# FLTIME Link fattening
# GFTIME Gauge force
# IOTIME I/O timing
# PRTIME print time (clover_invert2)

# REMAP  report remapping time for QDP, QOP in conjunction with above

CTIME ?= -DNERSC_TIME -DCGTIME -DFFTIME -DFLTIME -DGFTIME -DREMAP -DPRTIME -DIOTIME

#------------------------------
# Profiling
# Applications:  QDP

# QDP_PROFILE         Generates a report for all QDP routines

CPROF =#

#------------------------------
# Debugging and diagnostics
# Applications:  all

# COM_CRC             Message passing test.  Checksums on all gathers.
# CHECK_MALLOC        Report malloc/free activity.
#                     (Then process stdout using check_malloc.pl)
# CG_DEBUG            Print debugging information for the inverters.
# CG_OK               Print inverter convergence information even when OK
# REMAP_STDIO_APPEND  All nodes append to stdout.
# NO_FREOPEN          Don't use freopen to remap stdin and stdout
#
# HISQ_SVD_VALUES_INFO Print HISQ SVD diagnostics
# HISQ_SVD_COUNTER    Print summary count of SVD uses
# HISQ_FORCE_FILTER_COUNTER Print summary count of force filter applications.

CDEBUG = -DCG_OK -DREMAP_STDIO_APPEND -DHISQ_SVD_COUNTER -DHISQ_FORCE_FILTER_COUNTER # -DCHECK_MALLOC 

#------------------------------
# Backward compatibility

# As of version qopqdp 0.9.0 the normalization convention for the
# staggered inverter changed.  If you are using a version of QOPQDP
# with the old convention, define this macro:
CCOMPAT += #-DOLD_QOPQDP_NORM

# Prior to version 7.7.2 the conversion from staggeredd to naive was peculiar.
CCOMPAT += #-DOLD_STAGGERED2NAIVE
CCOMPAT += -DOLD_GAUSSRAND

#------------------------------
# Layout
# Applications: all

# These are currently selected only by editing Make_template.
#  Choices
#    layout_hyper_prime.o          Standard hypercubic
#    layout_timeslices.o           Puts full timeslices on each node.
#    layout_squares.o              Puts 2D slices on each node.
#    layout_hyper_tstretch.o       Rarely used.  Fewer timeslices on last node.
#                                  In case node no is incommensurate with nt
#    layout_timeslices_2.o         Untested.  Supposedly more forgiving.

#------------------------------
# Compute node grid layout

#     layout_hyper_prime views the machine as a grid and selects a layout
#     that distributes the lattice uniformly across the nodes and
#     minimizes the surface to volume ratio.  On fixed grid machines
#     such as the IBM Bluegene it may be better to control the
#     geometry.  In that case some of our applications support
#     the compiler macro FIX_NODE_GEOM.  You must then provide
#     and extra list of dimensions in the parameter input file.
#     See e.g. ks_imp_rhmc.

CGEOM ?=#-DFIX_NODE_GEOM

#------------------------------
# I/O node grid layout

#     For some applications and architectures it is more efficient to
#     split large files and have a subset of processors handle the I/O.
#     We define the I/O node partitions by hypercubes in the same manner 
#     as the compute node geometry above.  The I/O node geometry must
#     be commensurate with the compute node geometry.
#     For applications that support it, the I/O node geometry is specified
#     by the macro FIX_IONODE_GEOM.  Then the parameter input file
#     includes a list of dimensions.

CGEOM +=# -DFIX_IONODE_GEOM

#------------------------------
# Improved staggered CG inverter and Dslash
# Applications: ks_imp_dyn ks_imp_rhmc ks_imp_invert_multi ks_hl_spectrum

#  Note, some options still require editing Make_template
#  Choices 
#   dslash_fn.o                   Overlaps computation with backward gathers.
#   dslash_fn2.o                  Does all gathers before computation
#                                   but has fewer FORALLSITES loops
#   dslash_fn_dblstore.o          Double store version of dslash_fn.o
#                                   Supports GATHER13

# DBLSTORE_FN    Copies backward links.  Requires more memory.
#                You must also call for dslash_fn_dblstore.o for now.
# D_FN_GATHER13  Combine third and next neighbor gathers in Dslash.
#                For now, works only with dslash_fn_dblstore.o
# FEWSUMS        Fewer CG reduction calls

# If we are using QUDA, the backward links are unused, so we should
# avoid unecessary overhead and use the standard dslash.  Note that
# dslash_fn also has hooks in place to offload any dslash_fn_field
# calls to QUDA
ifeq ($(strip ${WANTQUDA}),true)
  KSCGSTORE = -DFEWSUMS
else
  KSCGSTORE = -DDBLSTORE_FN -DFEWSUMS -DD_FN_GATHER13
endif

#------------------------------
# Staggered fermion force routines
# Applications: ks_imp_dyn ks_imp_rhmc

# These are currently selected only by editing Make_template

# Choices
#  fermion_force_asqtad.o    Optimized for the Asqtad action
#  fermion_force_general.o   Takes any quark action

#------------------------------
# Prefetching
# Applications: all

# PREFETCH (Not working yet)

CPREFETCH = #

#------------------------------
# Multimass improved KS CG solvers
# Applications: ks_imp_rhmc ks_measure ks_spectrum

# Choices
# KS_MULTICG=OFFSET  The basic multicg solver.
# KS_MULTICG=HYBRID  Solve with multicg and polish off with single mass CG.
# KS_MULTICG=FAKE    Iterate the single mass solver.
# KS_MULTICG=REVERSE Iterate in reverse order
# KS_MULTICG=REVHYB  Same as HYBRID but with vectors in reverse order.
# NO_REFINE          No refinements except for masses with nonzero Naik eps
# CPU_REFINE         Refine on CPU only (if at all), not GPU
# PRIMME_PRECOND
# POLY_EIGEN
# MATVEC_PRECOND
# CHEBYSHEV_EIGEN

KSCGMULTI ?= -DKS_MULTICG=HYBRID # -DNO_REFINE # -DHALF_MIXED

#------------------------------
# Mixed precision
# enabled mixed-precision solvers for QUDA and GRID  (if set, overrides HALF_MIXED and MAX_MIXED macros below)
WANT_MIXED_PRECISION ?= 2
WANT_MIXED_PRECISION_GPU ?= ${WANT_MIXED_PRECISION}

# HALF_MIXED         (not QUDA or GRID) If PRECISION=2, do multimass solve in single precision
#                    and single-mass refinements in double
# HALF_MIXED         (QUDA) If PRECISION=2, use double-single mixed-precision solvers
# MAX_MIXED          (QUDA) Use double-half or single-half mixed-precision solvers
#                    (for multi-shift, behavior is as HALF_MIXED)

KSCGMIXED ?= #

ifeq ($(strip ${WANT_MIXED_PRECISION_GPU}),1)
  CGPU += -DHALF_MIXED # use single precision where appropriate
else ifeq ($(strip ${WANT_MIXED_PRECISION_GPU}),2)
  CGPU += -DMAX_MIXED # use half precision where appropriate
endif

#------------------------------
# Multi source types
# MULTISOURCE (deprecated)
# MULTICOLORSOURCE (deprecated)

KSSOURCE ?= #

#------------------------------
# Multifermion force routines
# Applications: ks_imp_rhmc

# Choices
# KS_MULTIFF=FNMAT    Construct matrix parallel transporters
#                     and use with sum of outer products of sources 
# KS_MULTIFF=FNMATREV Older version of FNMAT.  Traverses path
#                     in reverse order
# KS_MULTIFF=ASVEC    Use improved Asqtad, parallel transporting
#                     groups of source vectors.  See VECLENGTH.

# Additional options
# VECLENGTH=n        Number of source vectors to process in one group.
#                    Applies only to the ASVEC option

KSFFMULTI = -DKS_MULTIFF=FNMAT

#------------------------------
# KS gaussian smearing
GAUSS_SMEAR_KS_TWOLINK ?= true   # Compute two-link field to use in smearing

#------------------------------
# RHMC molecular dynamics algorithm
# Applications: ks_imp_rhmc

# Choices
# (Not needed if the make target has a preset value in Make_template)

# INT_ALG=INT_LEAPFROG
# INT_ALG=INT_OMELYAN
# INT_ALG=INT_2EPS_3TO1
# INT_ALG=INT_2EPS_2TO1
# INT_ALG=INT_2G1F
# INT_ALG=INT_3G1F
# INT_ALG=INT_4MN4FP
# INT_ALG=INT_4MN5FV
# INT_ALG=INT_FOURSTEP
# INT_ALG=INT_PLAY

KSRHMCINT =#

#------------------------------

# Staggered spin-taste operator shift for applications that construct
# nonlocal interpolating operators.  By default the shift operator is
# symmetric (forward and backward).  This macro makes it one-sided
# (forward).

# The Fat-Naik variants of the nonlocal operators are currently
# constructed only from the forward fat links.

KSSHIFT = # -DONE_SIDED_SHIFT

#------------------------------
# Clover inverter choice
# Applications: clover_invert, clover_invert2
#

# CL_CG=BICG  Biconjugate gradient
# CL_CG=CG    Standard CG
# CL_CG=MR    Minimum residue
# CL_CG=HOP   Hopping

# HALF_MIXED  Do double-precision inversion with single, or single with half (if supported)
# MAX_MIXED   Do double-precision inversion with half-precision (if supported)
# SCALE_PROP  Do rescaling for the clover propagator

CLCG = -DCL_CG=BICG 

#------------------------------
# Propagator storage
# Applications: clover_invert2
#

# CLOV_LEAN   Write intermediate propagators to disk, saving memory
#             but increasing run time somewhat

CLMEM = #-DCLOV_LEAN

#----------------------------------------------------------------------
# Extra include paths

INCADD = ${INCFFTW} ${INCPRIMME} ${INCQUDA} ${INCQPHIX} ${INCQPHIXJ} ${INCHADRONS} ${INCGRID} ${INCVTUNE}

#----------------------------------------------------------------------
#  Extra libraries

LIBADD = ${LIBFFTW} ${LIBPRIMME} ${LIBARPACK} ${LIBLAPACK} ${LIBQUDA} ${LIBQPHIX} \
  ${LIBQPHIXJ} ${LIBHADRONS} ${LIBGRID} ${LIBVTUNE}

#------------------------------
# Summary

CODETYPE = ${CTIME} ${CPROF} ${CDEBUG} ${CGEOM} ${CEIG} ${KSCGSTORE} ${CPREFETCH} \
 ${KSCGMULTI} ${KSCGMIXED} ${KSSOURCE} ${KSFFMULTI} ${KSRHMCINT} ${KSSHIFT} \
 ${CLCG} ${CLMEM} ${CQOP} ${CCOMPAT} ${CGITVER}

#----------------------------------------------------------------------
# MILC library make file in libraries directory.  
#    CHECK IT FOR FURTHER OPTIONS!

MAKELIBRARIES = Make_vanilla

#----------------------------------------------------------------------
# End of user choices.  Please, also, check choices in include/config.h.
#----------------------------------------------------------------------

# Definitions of compiler macros -- don't change.

ifeq ($(strip ${WANT_CL_BCG_GPU}),true)
  HAVE_CL_GPU = true
  CGPU += -DUSE_CL_GPU
endif

ifeq ($(strip ${WANT_FN_CG_GPU}),true)
  HAVE_FN_CG_GPU = true
  CGPU += -DUSE_CG_GPU
endif

ifeq ($(strip ${WANT_GA_GPU}),true)
  HAVE_GA_GPU = true
  CGPU += -DUSE_GA_GPU
endif

ifeq ($(strip ${WANT_GF_GPU}),true)
  HAVE_GF_GPU = true
  CGPU += -DUSE_GF_GPU
endif

ifeq ($(strip ${WANT_FL_GPU}),true)
  HAVE_FL_GPU = true
  CGPU += -DUSE_FL_GPU
endif

ifeq ($(strip ${WANT_FF_GPU}),true)
  HAVE_FF_GPU = true
  CGPU += -DUSE_FF_GPU
endif

ifeq ($(strip ${WANT_EIG_GPU}),true)
  HAVE_EIG_GPU = true
  CGPU += -DUSE_EIG_GPU
endif

ifeq ($(strip ${WANT_GSMEAR_GPU}),true)
  HAVE_GSMEAR_GPU = true
  CGPU += -DUSE_GSMEAR_GPU
endif

ifeq ($(strip ${WANT_KS_CONT_GPU}),true)
  HAVE_KS_CONT_GPU = true
  CGPU += -DUSE_KS_CONT_GPU
endif

ifeq ($(strip ${WANT_SHIFT_GPU}),true)
  HAVE_SHIFT_GPU = true
  CGPU += -DUSE_SHIFT_GPU
endif

ifeq ($(strip ${WANT_SPIN_TASTE_GPU}),true)
  HAVE_SPIN_TASTE_GPU = true
  CGPU += -DUSE_SPIN_TASTE_GPU
endif

ifeq ($(strip ${WANT_GAUGEFIX_OVR_GPU}),true)
  HAVE_GAUGEFIX_OVR_GPU = true
  CGPU += -DUSE_GAUGEFIX_OVR_GPU
endif

ifeq ($(strip ${WANT_MULTIGRID}),true)
  CGPU += -DMULTIGRID
endif

ifeq ($(strip ${OMP}),true)
  OCFLAGS += -DOMP
  OCXXFLAGS += -DOMP
endif

ifeq ($(strip ${MPP}),true)
  ifeq ($(strip ${HAVEQMP}),true)
     COMMTYPE = QMP
     COMMPKG = com_qmp.o
  else
     COMMTYPE = MPI_COMMS
     COMMPKG = com_mpi.o
  endif
else
  ifeq ($(strip ${HAVEQMP}),true)
     COMMTYPE = SINGLE
     COMMPKG = com_qmp.o
  else
     COMMTYPE = SINGLE
     COMMPKG = com_vanilla.o
  endif
endif

ifeq ($(strip ${HAVEQDP}),true)
  QDPPREC = -DQDP_PrecisionInt=${PRECISION}
endif

ifeq ($(strip ${HAVEQOP}),true)
  QOPPREC = -DQOP_PrecisionInt=${PRECISION}
endif

ifeq ($(strip ${WANTQPHIX}),true)
  QPHIXPREC = -DQPHIX_PrecisionInt=${PRECISION}
endif

ifeq ($(strip ${WANTQPHIXJ}),true)
  QPHIXPREC = -DQPHIXJ_PrecisionInt=${PRECISION}
endif

ifeq ($(strip ${WANTGRID}),true)
  GRIDPREC = -DGRID_PrecisionInt=${PRECISION}
endif

ifeq ($(strip ${WANTDCAP}),true)
   MACHINE_DEP_IO = io_dcap.o
   OCFLAGS += -I${DCAP_DIR}/include
   LDFLAGS += -L${DCAP_DIR}/${DCAPLIB} -Wl,--rpath,${DCAP_DIR}/${DCAPLIB} -ldcap
endif

ifeq ($(strip ${WANT_APE_IO}),true)
  HAVE_APE_IO = true
  OCFLAGS += -DAPE_LINKS_FILE
endif

ifeq ($(strip ${WANTFFTW}),true)
  HAVEFFTW = true
endif

ifeq ($(strip ${WANTPRIMME}),true)
  HAVE_PRIMME = true
  OCFLAGS += -DHAVE_PRIMME
endif

ifeq ($(strip ${WANTARPACK}),true)
  HAVE_ARPACK = true
  OCFLAGS += -DHAVE_ARPACK
endif

ifeq ($(strip ${GAUSS_SMEAR_KS_TWOLINK}),true)
  OCFLAGS += -DGAUSS_SMEAR_KS_TWOLINK
endif

# Make_template_combos defines convenience macros for interdependent
# groups of compilation units.  They are used to specify build lists.

include ../Make_template_combos

CPREC = -DMILC_PRECISION=${PRECISION} ${QDPPREC} ${QOPPREC} ${QPHIXPREC} ${GRIDPREC}
DARCH = ${CSCIDAC} ${CGPU} ${CPHI}

# Complete set of compiler flags - do not change
CFLAGS = ${OPT} ${OCFLAGS} -D${COMMTYPE} ${CODETYPE} ${INLINEOPT} \
	${CPREC} ${CLFS} ${INCSCIDAC} -I${MYINCLUDEDIR} ${DARCH} \
	${DEFINES} ${ADDDEFINES} ${IMPI} ${INCADD}
CXXFLAGS = ${OPT} ${OCXXFLAGS} -D${COMMTYPE} ${CODETYPE} ${INLINEOPT} \
        ${CPREC} ${CLFS} ${INCSCIDAC} -I${MYINCLUDEDIR} ${DARCH} \
	${DEFINES} ${ADDDEFINES} ${IMPI} ${INCADD}

ILIB = ${LIBSCIDAC} ${LMPI} ${LIBADD}

.PHONY: time check test_clean
time:
	make -f Make_time time

check: test_clean
	cd test ; perl ../../check.pl ${EXEC} ${CASE} ${PREC} < checklist

test_clean:
	cd test ; make test_clean

include Make_template
