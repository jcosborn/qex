const hdr = currentSourcePath()[0..^11] & "lapack.h"
const lapackLib {.strdefine.} = "-llapack -lblas"
#const lapackLib {.strdefine.} = "/usr/lib/lapack/liblapack.a -lblas -lgfortran"
#const lapackLib {.strdefine.} = "-L/usr/lib/lapack -llapack"
{.passL: lapackLib.}
#{.pragma: lapack, header: hdr.}
{.pragma: lapack.}

type
  fint* = cint
  #fint* = int64
  doublereal* = float64
  doublecomplex* = dcomplex
  scomplex* = object
    re*,im*: float32
  dcomplex* = object
    re*,im*: float64

proc cgemm*(transa: cstring; transb: cstring; m: ptr fint; n: ptr fint;
            k: ptr fint; alpha: ptr scomplex; a: ptr scomplex;
            lda: ptr fint; b: ptr scomplex; ldb: ptr fint;
            beta: ptr scomplex; c: ptr scomplex; ldc: ptr fint) {.
              lapack, importc: "cgemm_".}

proc zgemm*(transa: cstring; transb: cstring; m: ptr fint; n: ptr fint;
            k: ptr fint; alpha: ptr dcomplex; a: ptr dcomplex;
            lda: ptr fint; b: ptr dcomplex; ldb: ptr fint;
            beta: ptr dcomplex; c: ptr dcomplex; ldc: ptr fint) {.
              lapack, importc: "zgemm_".}

proc dsterf*(n: ptr fint, d: ptr float64, e: ptr float64,
             info: ptr fint) {.lapack, importC:"dsterf_".}

proc dstebz*(rnge: cstring; order: cstring; n: ptr fint;
             vl: ptr doublereal; vu: ptr doublereal; il: ptr fint;
             iu: ptr fint; abstol: ptr doublereal; d: ptr doublereal;
             e: ptr doublereal; m: ptr fint; nsplit: ptr fint;
             w: ptr doublereal; iblock: ptr fint; isplit: ptr fint;
             work: ptr doublereal; iwork: ptr fint; info: ptr fint) {.
               lapack, importc: "dstebz_".}

proc dgetrf*(m: ptr fint, n: ptr fint, a: ptr float64, lda: ptr fint,
             ipiv: ptr fint, info: ptr fint) {.lapack, importc:"dgetrf_".}

#[
proc zstegr*(JOBZ: ptr char,
             RANGE: ptr char,
	     n: ptr fint,
	     d: ptr float64, #
	     e: ptr float64,
	     vl: ptr float64,
	     vu: ptr float64,
	     il: ptr fint,
	     iu: ptr fint,
	     abstol: ptr float64,
	     m: ptr fint,
	     w: ptr float64,
	     z: ptr dcomplex, # dimension( ldz, * )
	     ldz: fint,
	     isuppz: ptr fint,
	     work: ptr float64,
	     lwork: fint,
	     iwork: ptr fint,
	     liwork: fint,
	     info: fint
)
]#

proc zheev*(jobz: cstring; uplo: cstring; n: ptr fint; a: ptr dcomplex;
            lda: ptr fint; w: ptr float64; work: ptr dcomplex; lwork: ptr fint;
            rwork: ptr float64; info: ptr fint) {.lapack, importc: "zheev_".}

proc zhegv*(itype: ptr fint; jobz: cstring; uplo: cstring; n: ptr fint;
            a: ptr doublecomplex; lda: ptr fint; b: ptr doublecomplex;
            ldb: ptr fint; w: ptr doublereal; work: ptr doublecomplex;
            lwork: ptr fint; rwork: ptr doublereal;
            info: ptr fint) {.lapack, importc: "zhegv_".}

proc zgeev*(jobvl: cstring; jobvr: cstring; n: ptr fint; a: ptr doublecomplex;
            lda: ptr fint; w: ptr doublecomplex; vl: ptr doublecomplex;
            ldvl: ptr fint; vr: ptr doublecomplex; ldvr: ptr fint;
            work: ptr doublecomplex; lwork: ptr fint; rwork: ptr doublereal;
            info: ptr fint) {.lapack, importc: "zgeev_".}

proc dbdsqr*(uplo: cstring; n: ptr fint; ncvt: ptr fint; nru: ptr fint;
             ncc: ptr fint; d: ptr float64; e: ptr float64;
             vt: ptr float64; ldvt: ptr fint; u: ptr float64;
             ldu: ptr fint; c: ptr float64; ldc: ptr fint; work: ptr float64;
             info: ptr fint) {.lapack, importc: "dbdsqr_".}

proc dbdsdc*(uplo: cstring; compq: cstring; n: ptr fint;
             d: ptr float64; e: ptr float64; u: ptr float64;
             ldu: ptr fint; vt: ptr float64; ldvt: ptr fint;
             q: ptr float64; iq: ptr fint; work: ptr float64;
             iwork: ptr fint; info: ptr fint) {.lapack, importc: "dbdsdc_"}

proc dbdsvdx*(uplo: cstring, jobz: cstring, range: cstring, n: ptr fint,
              d: ptr float64, e: ptr float64, vl: ptr float64, vu: ptr float64,
              il: ptr fint, iu: ptr fint, ns: ptr fint, s: ptr float64,
              z: ptr float64, ldz: ptr fint, work: ptr float64,
              iwork: ptr fint, info: ptr fint) {.lapack, importc: "dbdsvdx_".}

proc dlasq1*(n: ptr fint; d: ptr float64; e: ptr float64;
             work: ptr float64; info: ptr fint) {.lapack, importc: "dlasq1_".}


when isMainModule:
  template toPtrInt32(x: int): ptr int32 =
    var t = x.int32
    addr t
  template toPtrScomplex(x: int): ptr scomplex =
    var t = scomplex(re: x.float32, im: 0'f32)
    addr t
  template `&`(x: int): untyped = toPtrInt32(x)
  template `&&`(x: int): untyped = toPtrScomplex(x)
  var c,a,b: ptr scomplex
  var cr,cc,bc: int

  cgemm("C","N", &cc,&cr,&bc, &&1, b,&bc, a,&bc, &&0, c,&cc)
