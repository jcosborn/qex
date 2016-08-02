BGQ = $(shell [ -e /bgsys ] && echo "bgq")
ifeq ($(BGQ),bgq)
  QEXDIR = $(HOME)/nim/qex
  QMPDIR = $(HOME)/lqcd/install/qmp
  QIODIR = $(HOME)/lqcd/install/qio
  #CC = mpicc
  CC = $(QEXDIR)/mpixlc2
  #CC = mpixlc_r
  LD = $(CC)
  #CC_TYPE = clang
  CC_TYPE = gcc
  #CFLAGS_ALWAYS = "-Wall -std=gnu99"
  #CFLAGS_DEBUG = "-g3 -O0"
  #CFLAGS_SPEED = "-O3 -march=native"
  CFLAGS_ALWAYS = ""
  CFLAGS_DEBUG = "-g3 -O0"
  CFLAGS_SPEED = "-O3"
  VERBOSITY = 3
  SIMD = QPX
else
  QEXDIR = /media/data/weinbe2/nim/qex
  #QEXDIR = $(HOME)/lqcd/src/qex-201602290051
  QMPDIR = /media/data/weinbe2/lqcd/install/qmp
  QIODIR = /media/data/weinbe2/lqcd/install/qio
  CC = gcc
  LD = $(CC)
  CC_TYPE = gcc
  CFLAGS_ALWAYS = "-Wall -std=gnu99"
  CFLAGS_DEBUG = "-g3 -O0"
  CFLAGS_SPEED = "-O3 -march=native"
  VERBOSITY = 1
  #SIMD = SSE
  #SIMD = AVX
  SIMD = SSE,AVX
  VLEN = 8
endif

#### edit above

ifneq ($(origin QEXDIR), undefined)
  ENVS += QEXDIR=$(QEXDIR)
endif
ifneq ($(origin QMPDIR), undefined)
  ENVS += QMPDIR=$(QMPDIR)
endif
ifneq ($(origin QIODIR), undefined)
  ENVS += QIODIR=$(QIODIR)
endif
ifneq ($(origin CC_TYPE), undefined)
  ENVS += CC_TYPE=$(CC_TYPE)
endif
ifneq ($(origin CC), undefined)
  ENVS += CC=$(CC)
endif
ifneq ($(origin LD), undefined)
  ENVS += LD=$(LD)
endif
ifneq ($(origin CFLAGS_ALWAYS), undefined)
  ENVS += CFLAGS_ALWAYS=$(CFLAGS_ALWAYS)
endif
ifneq ($(origin CFLAGS_DEBUG), undefined)
  ENVS += CFLAGS_DEBUG=$(CFLAGS_DEBUG)
endif
ifneq ($(origin CFLAGS_SPEED), undefined)
  ENVS += CFLAGS_SPEED=$(CFLAGS_SPEED)
endif
ifneq ($(origin VERBOSITY), undefined)
  ENVS += VERBOSITY=$(VERBOSITY)
endif
ifneq ($(origin SIMD), undefined)
  ENVS += SIMD=$(SIMD)
endif
ifneq ($(origin VLEN), undefined)
  ENVS += VLEN=$(VLEN)
endif

DBG = "-d:release"
ifeq ($(debug),1)
  DBG = "-d:debug"
endif
ifeq ($(run),1)
  DBG += "-r"
endif
ifeq ($(c),1)
  DBG += "-c"
endif
DBG += "--warning[SmallLshouldNotBeUsed]:off"
DBG += "--hint[XDeclaredButNotUsed]:off"
#DBG += "--verbosity:2"
#DBG += "--listCmd"
DBG += "--implicitStatic:on"
#DBG += "--parallelBuild:4"

TARGETS = $(MAKECMDGOALS)

ifeq ($(TARGETS),clean)
  TARGETS = dummy
endif
clean:
	\rm -rf nimcache
ifeq ($(TARGETS),all)
  T1 = $(wildcard $(QEXDIR)/src/*.nim)
  T2 = $(notdir $(T1))
  T3 = $(basename $(T2))
  T4 = $(filter-out simdQpx,$(T3))
  TARGETS = $(T4)
endif
all: $(TARGETS)

$(TARGETS): _force
	$(ENVS) nim c $(DBG) --nimcache:`pwd`/nimcache -o:`pwd`/$(@F) $(QEXDIR)/src/$@

.PHONY: _force
