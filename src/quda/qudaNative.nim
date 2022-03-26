import quda, enum_quda

import os

# TODO select cpp when import this file?

when not defined(qudaDir):
  {.fatal:"Must define qudaDir to use QUDA.".}
when not defined(cudaLibDir):
  {.fatal:"Must define cudaLibDir to use QUDA.".}
const qudaDir {.strdefine.} = ""
const qudaTestUtilsDir {.strdefine.} = "../../quda/tests/utils"
const cudaLibDir {.strdefine.} = ""
const cudaLib = "-L" & cudaLibDir & " -lcudart -lcublas -lcufft -Wl,-rpath," & cudaLibDir & " -L" & cudaLibDir & "/stubs -lcuda"
{.passC: "-I" & qudaDir & "/include -I" & cudaLibDir & "/../include".}

const qmpDir {.strdefine.} = getEnv("QMPDIR")
const qioDir {.strdefine.} = getEnv("QIODIR")

when qioDir.len > 0:
  when qmpDir.len > 0:
    # Assume quda is built with QIO and QMP.
    {.passL: " -L" & qudaDir & "/lib -lquda -Wl,-rpath," & qudaDir & "/lib " & cudaLib & " -L" & qioDir & "/lib -lqio -llime -L" & qmpDir & "/lib -lqmp".}
  else:
    # Assume QUDA is built with QIO.
    {.passL: " -L" & qudaDir & "/lib -lquda -Wl,-rpath," & qudaDir & "/lib " & cudaLib & " -L" & qioDir & "/lib -lqio -llime".}
else:
  {.passL: " -L" & qudaDir & "/lib -lquda -Wl,-rpath," & qudaDir & "/lib " & cudaLib.}

proc delete[T](x:ptr T) {.importcpp: "delete @", nodecl.}

########
const lfh = "lattice_field.h"
type LatticeFieldParam {.importcpp:"quda::LatticeFieldParam",header:lfh,bycopy.} = object
proc constructLatticeFieldParam(param:QudaGaugeParam):LatticeFieldParam {.importcpp:"quda::LatticeFieldParam(@)",constructor,header:lfh.}

########
const csfh = "color_spinor_field.h"
type
  ColorSpinorParam {.importcpp:"quda::ColorSpinorParam",header:csfh,bycopy.} = object
  ColorSpinorField {.importcpp:"quda::ColorSpinorField",header:csfh,bycopy.} = object
proc ColorSpinorFieldCreate(param:ColorSpinorParam):ptr ColorSpinorField {.importc:"quda::ColorSpinorField::Create",header:csfh.}
proc Source(this:ColorSpinorField, sourceType:QudaSourceType) {.importcpp:"#.Source(@)",header:csfh.}
proc Vptr(this:ColorSpinorField):pointer {.importcpp:"#.V(@)",header:csfh.}

########
const rqh = "random_quda.h"
type RNG {.importcpp:"quda::RNG",header:rqh,bycopy.} = object
proc newRNG(param:LatticeFieldParam, seedin:uint64):ptr RNG {.importcpp:"new quda::RNG(@)",header:rqh.}
proc Init(this:RNG) {.importcpp:"#.Init(@)",header:rqh.}
proc Release(this:RNG) {.importcpp:"#.Release(@)",header:rqh.}

########
const clph = qudaTestUtilsDir & "/command_line_params.h"
var
  mass{.importc,header:clph.}:cdouble
  kappa{.importc,header:clph.}:cdouble
  mu{.importc,header:clph.}:cdouble
  epsilon{.importc,header:clph.}:cdouble
  m5{.importc,header:clph.}:cdouble
  b5{.importc,header:clph.}:cdouble
  c5{.importc,header:clph.}:cdouble
  dim{.importc,header:clph.}:array[4,cint]  # actually a std::array<int,4>
  Lsdim{.importc,header:clph.}:cint
  dslash_type{.importc,header:clph.}:QudaDslashType

########
{.passL: " -L" & qudaDir & "/tests -lquda_test -Wl,-rpath," & qudaDir & "/lib".}
const huh = qudaTestUtilsDir & "/host_utils.h"
var
  V{.importc,header:huh.}:cint
  cpu_prec{.importc,header:huh.}:QudaPrecision
  gaugeSiteSize{.importc:"gauge_site_size",header:huh.}:cint
  cloverSiteSize{.importc:"clover_site_size",header:huh.}:cint
  hostGaugeDataTypeSize{.importc:"host_gauge_data_type_size",header:huh.}:csize_t
  hostCloverDataTypeSize{.importc:"host_clover_data_type_size",header:huh.}:csize_t
  hostSpinorDataTypeSize{.importc:"host_spinor_data_type_size",header:huh.}:csize_t
proc setQudaPrecisions {.importc,header:huh.}
proc setQudaDefaultMgTestParams {.importc,header:huh.}
proc setWilsonGaugeParam(gauge_param:QudaGaugeParam) {.importc,header:huh.}
proc setInvertParam(inv_param:QudaInvertParam) {.importc,header:huh.}
proc setDims(X:ptr cint) {.importc,header:huh.}
proc dw_setDims(X:ptr cint, L5:int) {.importc,header:huh.}
proc constructHostGaugeField(gauge:ptr pointer, gauge_param:QudaGaugeParam, argc:cint, argv:ptr ptr cchar) {.importc,header:huh.}
proc constructHostCloverField(clover:pointer, clover_inv:pointer, inv_param:QudaInvertParam) {.importc,header:huh.}
proc constructWilsonTestSpinorParam(csParam:ptr ColorSpinorParam, inv_param:ptr QudaInvertParam, gauge_param:ptr QudaGaugeParam) {.importc,header:huh.}

#[ using host_utils for now
proc setGaugeParam(p:var QudaGaugeParam, lo:any) =
  p.type = QUDA_SU3_LINKS
  for i in 0..<4:
    p.X[i] = lo.physGeom[i].cint
  p.anisotropy = 1.0
  p.tadpoleCoeff = 1.0
  p.gaugeFix = QUDA_GAUGE_FIXED_NO

proc setWilsonGaugeParam(p:var QudaGaugeParam, lo:any) =
  setGaugeParam(p,lo);
  p.type = QUDA_WILSON_LINKS
  p.gaugeOrder = QUDA_QDP_GAUGE_ORDER
  var
    xf = p.X[1]*p.X[2]*p.X[3]/2
    yf = p.X[0]*p.X[2]*p.X[3]/2
    zf = p.X[0]*p.X[1]*p.X[3]/2
    tf = p.X[0]*p.X[1]*p.X[2]/2
  p.gaPad = max(xf,yf,zf,tf)

proc setInvertParam(p:var QudaInvertParam, d:QudaDslashType, mass=0.1, kappa=-1.0, cloverCoeff=0.1) =
  let anisotropy = 1.0
  p.dslashType = d
  if kappa == -1.0:
    p.mass = mass
    p.kappa = 1.0 / (2.0 * (1.0 + 3.0 / anisotropy + mass))
    if d == QUDA_LAPLACE_DSLASH:
      p.kappa = 1.0 / (8.0 + mass)
  else:
    p.kappa = kappa;
    p.mass = 0.5 / kappa - (1.0 + 3.0 / anisotropy)
    if d == QUDA_LAPLACE_DSLASH:
      p.mass = 1.0 / kappa - 8.0
  echo "Kappa = ",p.kappa," Mass = ",p.mass
  p.laplace3D = 4
  if d == QUDA_DOMAIN_WALL_DSLASH or d == QUDA_DOMAIN_WALL_4D_DSLASH or d == QUDA_MOBIUS_DWF_DSLASH or d == QUDA_MOBIUS_DWF_EOFA_DSLASH:
    echo "DWF Unimplemented"
    # TODO
    #[
    p.m5 = m5
    kappa5 = 0.5 / (5.0 + p.m5)
    p.Ls = Lsdim
    for k in 0..<Lsdim:  # for mobius only
      # b5[k], c[k] values are chosen for arbitrary values,
      # but the difference of them are same as 1.0
      p.b_5[k] = b5
      p.c_5[k] = c5
    p.eofa_pm = eofa_pm
    p.eofa_shift = eofa_shift
    p.mq1 = eofa_mq1
    p.mq2 = eofa_mq2
    p.mq3 = eofa_mq3
    ]#
  else:
    p.Ls = 1
  if d == QUDA_CLOVER_WILSON_DSLASH or d == QUDA_TWISTED_CLOVER_DSLASH:
    p.clover_cpu_prec = QUDA_DOUBLE_PRECISION;
    p.clover_cuda_prec = QUDA_DOUBLE_PRECISION;
    p.clover_cuda_prec_sloppy = QUDA_HALF_PRECISION;
    p.clover_cuda_prec_precondition = QUDA_DOUBLE_PRECISION;
    p.clover_cuda_prec_eigensolver = QUDA_DOUBLE_PRECISION;
    p.clover_cuda_prec_refinement_sloppy = QUDA_SINGLE_PRECISION;
    p.clover_order = QUDA_PACKED_CLOVER_ORDER;
    p.clover_coeff = clover_coeff;
]#

when isMainModule:
  import qex
  letParam:
    gaugefile = ""
    lat =
      if fileExists(gaugefile):
        getFileLattice gaugefile
      else:
        if gaugefile.len > 0:
          qexWarn "Nonexistent gauge file: ", gaugefile
        @[8,8,8,8]
    dslash = ""

  # from quda, default values in tests/utls/command_line_params.cpp
  mass = floatParam("mass", mass)
  kappa = floatParam("kappa", kappa)
  mu = floatParam("mu", mu)
  epsilon = floatParam("epsilon", epsilon)
  m5 = floatParam("m5", m5)
  b5 = floatParam("b5", b5)
  c5 = floatParam("c5", c5)
  Lsdim = intParam("Lsdim", 4).cint
  for i in 0..<4: dim[i] = lat[i].cint
  if dslash.len>0:
    dslash_type = case dslash
      of "wilson": QUDA_WILSON_DSLASH
      of "clover": QUDA_CLOVER_WILSON_DSLASH
      of "twisted-mass": QUDA_TWISTED_MASS_DSLASH
      of "twisted-clover": QUDA_TWISTED_CLOVER_DSLASH
      of "clover-hasenbusch-twist": QUDA_CLOVER_HASENBUSCH_TWIST_DSLASH
      of "staggered": QUDA_STAGGERED_DSLASH
      of "asqtad": QUDA_ASQTAD_DSLASH
      of "domain-wall": QUDA_DOMAIN_WALL_DSLASH
      of "domain-wall-4d": QUDA_DOMAIN_WALL_4D_DSLASH
      of "mobius": QUDA_MOBIUS_DWF_DSLASH
      of "mobius-eofa": QUDA_MOBIUS_DWF_EOFA_DSLASH
      of "laplace": QUDA_LAPLACE_DSLASH
      else:
        qexWarn "Unrecognized dslash: '" & dslash & "'.  Using default wilson."
        dslash_type

  qexInit()
  echoParams()
  echo "dslash_type: ",dslash_type
  echo "rank ", myRank, "/", nRanks
  threads: echo "thread ", threadNum, "/", numThreads

  let
    lo = lat.newLayout
  var
    rgC: array[4,cint]

  proc qudaCommsMap(coords0: ptr cint; fdata: pointer): cint {.cdecl.} =
    let pl = cast[ptr type(lo)](fdata)
    let coords = cast[ptr UncheckedArray[cint]](coords0)
    let r = pl[].rankFromRankCoords(coords)
    r.cint
  for i in 0..<4: rgC[i] = lo.rankGeom[i].cint

  # --------

  setQudaDefaultMgTestParams()

  setQudaPrecisions()

  initCommsGridQuda(rgC.len.cint, rgC[0].addr, qudaCommsMap, unsafeAddr(lo))

  var
    gaugeParam = newQudaGaugeParam()
    invParam = newQudaInvertParam()

  setWilsonGaugeParam(gaugeParam)
  setInvertParam(invParam)

  initQuda(0)

  if dslash_type == QUDA_DOMAIN_WALL_DSLASH or dslash_type == QUDA_DOMAIN_WALL_4D_DSLASH or
      dslash_type == QUDA_MOBIUS_DWF_DSLASH or dslash_type == QUDA_MOBIUS_DWF_EOFA_DSLASH:
    dw_setDims(gaugeParam.X[0].addr, invParam.Ls)
  else:
    setDims(gaugeParam.X[0].addr)

  var gauge: array[4,pointer]
  for i in 0..<4:
    gauge[i] = alloc(csize_t(V * gaugeSiteSize) * hostGaugeDataTypeSize)
  constructHostGaugeField(gauge[0].addr, gaugeParam, 0, nil)
  loadGaugeQuda(gauge[0].addr, gaugeParam.addr)

  var
    clover = alloc(csize_t(V * cloverSiteSize) * hostCloverDataTypeSize)
    cloverInv = alloc(csize_t(V * cloverSiteSize) * hostSpinorDataTypeSize)
  if dslash_type == QUDA_CLOVER_WILSON_DSLASH or dslash_type == QUDA_TWISTED_CLOVER_DSLASH:
    constructHostCloverField(clover, cloverInv, invParam)
    loadCloverQuda(clover, cloverInv, invParam.addr)

  var plaq: array[3,cdouble]
  plaqQuda(plaq)
  echo "Computed plaquette is ",plaq[0],", (spatial = ",plaq[1],", temporal = ",plaq[2],")"

  var csParam:ColorSpinorParam
  constructWilsonTestSpinorParam(csParam.addr, invParam.addr, gaugeParam.addr)
  # var check = ColorSpinorFieldCreate(csParam)

  var rng = newRNG(constructLatticeFieldParam(gaugeParam), 1234)
  Init(rng[])

  var
    fin = ColorSpinorFieldCreate(csParam)
    fout = ColorSpinorFieldCreate(csParam)
  Source(fin[], QUDA_RANDOM_SOURCE)

  invertQuda(Vptr(fout[]), Vptr(fin[]), invParam.addr)
  echo "Done: ",invParam.iter," iter / ",invParam.secs," secs = ",invParam.gflops/invParam.secs," Gflops"

  Release(rng[])
  delete(rng)

  # TODO verify inversion

  freeGaugeQuda()
  for i in 0..<4:
    dealloc(gauge[i])

  if dslash_type == QUDA_CLOVER_WILSON_DSLASH or dslash_type == QUDA_TWISTED_CLOVER_DSLASH:
    freeCloverQuda()
    dealloc(clover)
    dealloc(cloverInv)

  # --------

  endQuda()
  # finalizeComms()
  qexFinalize()
  echoTimers()
