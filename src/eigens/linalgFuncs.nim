import lapack
import base
import maths

type
  dvec* = ref object
    dat*: ptr carray[float64]
    `len`*: int
    stride: int
    isSub: bool
  dmat* = ref object
    dat*: ptr carray[float64]
    nrows*, ncols*: int
    rowStride, colStride: int
    isSub: bool
  zmat* = ref object
    dat*: ptr carray[ComplexType[float64]]
    nrows*, ncols*: int
    rowStride, colStride: int
    isSub: bool

proc dvec_free(x: dvec) =
  if not x.isSub:
    dealloc(x.dat)
template dvec_alloc*(x: var dvec, n: int) =
  x.new(dvec_free)
  x.dat = cast[type(x.dat)](alloc(n*sizeof(float64)))
  x.`len` = n
  x.stride = 1
  x.isSub = false
proc newDvec*(n: int): dvec =
  result.dvec_alloc(n)
template `&`*(x: dvec): untyped = cast[ptr carray[float]](addr x[0])
template vec_size(x: dvec): int = x.len
template `[]`*(x: dvec, i: int): untyped = x.dat[i]
template `[]=`*(x: dvec, i: int, y: untyped): untyped = x.dat[i] = y
template dvec_get*(x: dvec, i: int): float64 = x[i]
template dvec_set*(x: dvec, i: int, y: float64) = x[i] = y
template dsubvec(x: dvec, y: dvec, i,n: int) =
  x.new
  x.dat = cast[type(x.dat)](addr(y[i]))
  x.`len` = n
  x.stride = 1
  x.isSub = true
template v_eq_zero(x: dvec) =
  for i in 0..<x.len: x[i] = 0.0
template v_eq_v(x,y: dvec) =
  for i in 0..<x.len: x[i] = y[i]
template norm2_v(x: dvec): float =
  var r = 0.0
  for i in 0..<x.len: r += x[i]*x[i]
  r
template dot*(x,y: dvec): float =
  var r = 0.0
  for i in 0..<x.len: r += x[i]*y[i]
  r
template v_eq_r_times_v(x: dvec, r: float, y: dvec) =
  for i in 0..<x.len: x[i] = r*y[i]
template v_peq_r_times_v(x: dvec, r: float, y: dvec) =
  for i in 0..<x.len: x[i] += r*y[i]
template v_meq_r_times_v(x: dvec, r: float, y: dvec) =
  for i in 0..<x.len: x[i] -= r*y[i]
template daxpy*(a: float, x: dvec, y: dvec) = v_peq_r_times_v(y, a, x)
proc normalize*(x: dvec) =
  let n2 = norm2_v(x)
  let s = 1.0/sqrt(n2)
  for i in 0..<x.len: x[i] *= s

proc dmat_free(x: dmat) =
  if not x.isSub:
    dealloc(x.dat)
template dmat_alloc*(x: dmat, nr,nc: int) =
  x.new(dmat_free)
  x.dat = cast[type(x.dat)](alloc(nr*nc*sizeof(float64)))
  x.nrows = nr
  x.ncols = nc
  x.rowStride = 1
  x.colStride = nr
  x.isSub = false
proc newDmat*(nr,nc: int): dmat =
  result.dmat_alloc(nr, nc)

proc zmat_free(x: zmat) =
  if not x.isSub:
    dealloc(x.dat)
template zmat_alloc*(x: zmat, nr,nc: int) =
  x.new(zmat_free)
  x.dat = cast[type(x.dat)](alloc(nr*nc*2*sizeof(float64)))
  x.nrows = nr
  x.ncols = nc
  x.rowStride = 1
  x.colStride = nr
  x.isSub = false
proc newZmat*(nr,nc: int): zmat =
  result.zmat_alloc(nr, nc)


template mat_nrows*(x: dmat): int = x.nrows
template mat_ncols*(x: dmat): int = x.ncols
template `[]`*(x: dmat|zmat, i,j: int): untyped = x.dat[i+j*x.nrows]
template `[]=`*(x: dmat|zmat, i,j: int, y: untyped): untyped =
  x.dat[i+j*x.nrows] = y
template dmat_get*(x: dmat, i,j: int): float64 = x[i,j]
template dmat_set(x: dmat, i,j: int, y: float64) = x[i,j] = y
template norm2*(x: dmat): float =
  var r = 0.0
  for i in 0..<x.nrows:
    for j in 0..<x.ncols:
      r += x[i,j]*x[i,j]
  r
proc `+=`*(r,x: dmat) =
  for i in 0..<r.nrows:
    for j in 0..<r.ncols:
      r[i,j] += x[i,j]
proc `-=`*(r,x: dmat) =
  for i in 0..<r.nrows:
    for j in 0..<r.ncols:
      r[i,j] -= x[i,j]

template dcolvec*(v: dvec, m: dmat, k: int) =
  v.new
  v.`len` = m.nrows
  v.stride = m.rowStride
  v.isSub = true
  v.dat = cast[type(v.dat)](addr(m[0,k]))
template colvec*(v: dvec, m: dmat, k: int) = dcolvec(v, m, k)

type
  zheevTmp* = object
    work*: ptr dcomplex
    rwork*: ptr cdouble
    lwork*: fint
    n*: fint

proc zewrk*(n: int): ptr zheevTmp =
  var w {.global.}: zheevTmp
  var s {.global.}: int = 0
  if n > s:
    if s > 0:
      dealloc(w.work)
      dealloc(w.rwork)
    s = n
    w.n = n.fint
    w.lwork = 3 * n.fint
    w.work = cast[ptr dcomplex](alloc(w.lwork*sizeof(dcomplex)))
    w.rwork = cast[ptr cdouble](alloc(3*n*sizeof(cdouble)))
  return addr(w)

proc zeigs*(m: ptr float64; d: ptr float64; n: int) =
  var info: fint
  var t = zewrk(n)
  var nn = fint(n)
  var an = addr nn
  var mm = cast[ptr dcomplex](m)
  zheev("V", "L", an, mm, an, d, t.work, addr(t.lwork), t.rwork, addr(info))

type
  zgeev_tmp* = object
    work*: ptr dcomplex
    rwork*: ptr cdouble
    lwork*: fint
    n*: fint

proc zgewrk*(n: int): ptr zgeev_tmp =
  var w {.global.}: zgeev_tmp
  var s {.global.}: int = 0
  if n > s:
    if s > 0:
      dealloc(w.work)
      dealloc(w.rwork)
    s = n
    w.n = fint(n)
    w.lwork = fint(64 * n + 1)
    w.work = cast[ptr dcomplex](alloc(w.lwork*sizeof(dcomplex)))
    w.rwork = cast[ptr cdouble](alloc((2*n+1)*sizeof(cdouble)))
  return addr(w)

proc zgeigs*(m: ptr float64; e: ptr float64; n: int) =
  var info: fint
  var t = zgewrk(n)
  var nn = fint(n)
  var an = addr nn
  var mm = cast[ptr dcomplex](m)
  var ee = cast[ptr dcomplex](e)
  zgeev("N", "N", an, mm, an, ee, nil, an, nil, an, t.work,
        addr(t.lwork), t.rwork, addr(info))

#[
proc zgeigsv*(m: ptr zmat; d: ptr zvec; vl: ptr zmat; vr: ptr zmat; n: cint) = 
  var info: cint
  var t: ptr zgeev_tmp = zgewrk(mat_nrows(m[]))
  if mat_row_stride(m[]) == 1: 
    zgeev("V", "V", addr(mat_nrows(m[])), mat_data(m[]), addr(mat_ltd(m[])), 
          vec_data(d[]), mat_data(vl[]), addr(mat_ltd(vl[])), mat_data(vr[]), 
          addr(mat_ltd(vr[])), t.work, addr(t.lwork), t.rwork, addr(info))
    if mat_conjugate(m[]) != mat_conjugate(vl[]): zmat_conjugate(vl)
    if mat_conjugate(m[]) != mat_conjugate(vr[]): zmat_conjugate(vr)
  else: 
    zgeev("V", "V", addr(mat_nrows(m[])), mat_data(m[]), addr(mat_ltd(m[])), 
          vec_data(d[]), mat_data(vr[]), addr(mat_ltd(vr[])), mat_data(vl[]), 
          addr(mat_ltd(vl[])), t.work, addr(t.lwork), t.rwork, addr(info))
    if mat_conjugate(m[]) == mat_conjugate(vl[]): zmat_conjugate(vl)
    if mat_conjugate(m[]) == mat_conjugate(vr[]): zmat_conjugate(vr)
  if mat_row_stride(vl[]) != 1: zmat_transpose(vl)
  if mat_row_stride(vr[]) != 1: zmat_transpose(vr)
  if mat_conjugate(m[]) != vec_conjugate(d[]): zvec_conjugate(d)
]#

# Hermetian generalized eigenvalues
# A x = lambda B x
proc zeigsgv*(a: ptr float64; b: ptr float64; e: ptr float64; n: int) =
  var itype = 1.fint
  var jobz = "V"
  var uplo = "L"
  var nn = fint(n)
  var an = addr nn
  #var aa = cast[ptr dcomplex](a)
  #var bb = cast[ptr dcomplex](b)
  var sz = n*n*sizeof(dcomplex)
  var aa = cast[ptr dcomplex](alloc(sz))
  var bb = cast[ptr dcomplex](alloc(sz))
  copyMem(aa, a, sz)
  copyMem(bb, b, sz)
  var lwork = (65*n).fint
  var work = cast[ptr dcomplex](alloc(lwork*sizeof(dcomplex)))
  var rwork = cast[ptr float64](alloc(3*n*sizeof(float64)))
  var info = fint(0)
  var eps = 1e-16
  while true:
    zhegv(addr itype, jobz, uplo, an, aa, an, bb, an, e,
          work, addr lwork, rwork, addr info)
    if info==0: break
    echo "error: zhegv info: ", info, "  eps: ", eps
    copyMem(aa, a, sz)
    copyMem(bb, b, sz)
    let bx = cast[ptr carray[dcomplex]](bb)
    for i in 0..<n:
      bx[i*(n+1)].re += eps
    eps *= 10.0
  copyMem(a, aa, sz)
  copyMem(b, bb, sz)
  dealloc aa
  dealloc bb
  dealloc work
  dealloc rwork

proc svdbi1*(ev: ptr carray[float64], a: ptr carray[float64],
             b: ptr carray[float64], n: int) =
  var n2 = fint(2*n)
  var d = cast[ptr carray[float64]](alloc0(n2*sizeof(float64)))
  var e = cast[ptr carray[float64]](alloc(n2*sizeof(float64)))
  e[0] = a[0]
  for i in 0..(n-2):
    e[2*i+1] = b[i]
    e[2*i+2] = a[i+1]
  var info = fint(0)
  dsterf(addr(n2), addr(d[0]), addr(e[0]), addr(info))
  for i in 0..<n:
    ev[i] = d[n+i]
  dealloc(d)
  dealloc(e)

proc svdbi*(ev: ptr carray[float64], a: ptr carray[float64],
            b: ptr carray[float64], n: int) =
  var nn = fint(n)
  var d = cast[ptr carray[float64]](alloc(n*sizeof(float64)))
  var e = cast[ptr carray[float64]](alloc(n*sizeof(float64)))
  for i in 0..(n-2):
    d[i] = a[i]
    e[i] = b[i]
  d[n-1] = a[n-1]
  var work = cast[ptr carray[float64]](alloc(4*n*sizeof(float64)))
  var info = fint(0)
  dlasq1(addr(nn), addr(d[0]), addr(e[0]), addr(work[0]), addr(info))
  for i in 0..<n:
    ev[i] = d[n-1-i]
  dealloc(work)
  dealloc(e)
  dealloc(d)


proc svdBidiag1*(d: ptr carray[float64], e: ptr carray[float64],
                 v: ptr carray[float64], u: ptr carray[float64], n,k: int) =
  let uplo = "U"
  var nn = fint(n)
  var an = addr nn
  #var nv = ffint(k)
  var nv = fint(n)
  var ak = addr nv
  var nc = fint(0)
  var vv = cast[ptr carray[float64]](alloc(n*nv*sizeof(float64)))
  var vu = cast[ptr carray[float64]](alloc(n*nv*sizeof(float64)))
  var work = cast[ptr float64](alloc(4*n*sizeof(float64)))
  var info = fint(0)

  echo nn
  echo nv, "  ", ak[]
  # dbdsqr computes all singular vectors
  dbdsqr(uplo, an, ak, ak, addr nc, &d, &e, &vv, an, &vu, ak,
         nil, an, work, addr info)
  echo "info: ", info
  var sv = 0.0
  var su = 0.0
  for i in 0..<(n*nv):
    sv += vv[i]*vv[i]
    su += vu[i]*vu[i]
  echo "v2: ", sv, "  u2: ", su

  if not u.isNil:
    for i in 0..<n:
      for j in 0..<k:
        #v[j*n+i] = vv[j*n+i]
        #u[j*n+i] = vu[i*k+j]
        v[j*n+i] = vv[i*k+j]
        u[j*n+i] = vu[j*n+i]
        #echo i, " ", j, "  ", v[j*n+i], "  ", u[j*n+i]
  dealloc vv
  dealloc vu
  dealloc work

proc svdBidiag*(d: ptr carray[float64], e: ptr carray[float64],
                v: ptr carray[float64], u: ptr carray[float64], n,k: int) =
  let uplo = "U"
  let compq = "I"
  var nn = fint(n)
  var work = cast[ptr float64](alloc((3*n*n+4*n)*sizeof(float64)))
  var iwork = cast[ptr fint](alloc(2*8*n*sizeof(fint)))
  var info = fint(0)
  var vu = cast[ptr carray[float64]](alloc(n*n*sizeof(float64)))
  var vv = cast[ptr carray[float64]](alloc(n*n*sizeof(float64)))

  # dbdsdc calculates all singular vectors
  dbdsdc(uplo, compq, addr nn, &d, &e, &vu, addr nn, &vv, addr nn,
         nil, nil, work, iwork, addr info)
  if info!=0:
    echo "dbdsdc info: ", info

  for i in 0..<k:
    let j = n-1-i
    if i<j:
      swap(d[i], d[j])

  if not u.isNil:
    for i in 0..<k:
      let ii = n-1-i
      for j in 0..<n:
        v[i*n+j] = vv[j*n+ii]
        u[i*n+j] = vu[ii*n+j]

  dealloc vv
  dealloc vu
  dealloc work
  dealloc iwork

proc svdBidiag3*(d: ptr carray[float64], e: ptr carray[float64],
                 v: ptr carray[float64], u: ptr carray[float64], n,k: int) =
  let nv = fint(max(k,1))
  let uplo = "U"
  let jobz = if k>0: "V" else: "N"
  let rnge = if k>0: "I" else: "A"
  var nn = fint(n)
  var vl = 0.0
  var vu = 1e99
  var il = fint(n-nv+1)
  var iu = fint(n)
  var ns = fint(0)
  var s = cast[ptr carray[float64]](alloc(n*sizeof(float64)))
  #echo cast[int](s)
  var z = cast[ptr carray[float64]](alloc(2*n*(nv+1)*sizeof(float64)))
  #echo cast[int](z)
  var ldz = fint(2*n)
  var work = cast[ptr float64](alloc(14*n*sizeof(float64)))
  #echo cast[int](work)
  var iwork = cast[ptr fint](alloc(2*12*n*sizeof(fint)))  ## check
  #echo cast[int](iwork)
  var info = fint(0)
  #echo "here1"

  # #[
  dbdsvdx(uplo, jobz, rnge, addr nn, &d, &e, addr vl, addr vu,
          addr il, addr iu, addr ns, &s, &z, addr ldz,
          work, iwork, addr info)
  #        ]#
  #echo "here1"
  echo "info: ", info

  for i in 0..<k:
    d[i] = s[k-1-i]
    #echo i, "  ", d[i]

  if not u.isNil:
    for i in 0..<k:
      #let ii = 2*n*(k-1-i)
      let ii = 2*n*i
      for j in 0..<n:
        u[i*n+j] = z[ii+j]
        v[i*n+j] = z[ii+n+j]
  #echo "here1"

  #echo cast[int](s)
  dealloc s
  #echo cast[int](z)
  dealloc z
  #echo cast[int](work)
  dealloc work
  #echo cast[int](iwork)
  dealloc iwork
  #echo "here2"

proc dsvdbi2*(a: dvec, b: dvec, n: int) =
  var d = cast[ptr carray[float64]](alloc(n*sizeof(float64)))
  var e = cast[ptr carray[float64]](alloc(n*sizeof(float64)))
  for i in 0..<n:
    d[i] = a[i]
    if i<n-1: e[i] = b[i]
  svdBidiag(d, e, nil, nil, n, 0)
  for i in 0..<n:
    a[i] = d[i]
    if i<n-1: b[i] = e[i]
  dealloc(d)
  dealloc(e)

import eigens/svdBi4
export svdBi4

when isMainModule:
  import base/basicOps
  import maths/complexType
  template adj*(x: SomeNumber): untyped = x
  import times
  type
    Dcmplx = Complex[float64,float64]
    Zmat = object
      nrows,ncols: int
      data: seq[Dcmplx]

  proc newDcmplx(x,y: float64): Dcmplx =
    newComplex(x,y)
  proc newZmat(nrows,ncols: int): Zmat =
    result.nrows = nrows
    result.ncols = ncols
    result.data.newSeq(nrows*ncols)

  template `[]`(x: Zmat): untyped = addr(x.data[0])
  template `[]`(x: Zmat, i,j: untyped): untyped =
    x.data[i*x.ncols+j]
  template `[]=`(x: Zmat, i,j,y: untyped): untyped =
    x.data[i*x.ncols+j] = y

  var nr = 4
  var nc = nr
  var m1 = newZmat(nr,nc)
  var m2 = newZmat(nr,nc)
  for i in 0..<nr:
    for j in 0..<nc:
      m1[i,j] = newDcmplx(i.float+j.float,i.float-j.float)
      m2[i,j] = m1[i,j]
    #m[i,i] = newDcmplx(i.float+0.1,0.0)
  var ev = newSeq[float64](nr)
  template `[]`*(x: seq): untyped = addr(x[0])

  proc testZeigs =
    zeigs(cast[ptr float64](m2[]), ev[], nr)

    for i in 0..<nr:
      echo ev[i]

    var x = newSeq[Dcmplx](nr)
    for i in 0..<nr:
      for j in 0..<nr:
        x[j] := 0
        for k in 0..<nr:
          #x[j] += m1[k,j]*m2[i,k]
          x[j] += m1[j,k]*m2[i,k].adj
        echo x[j] - ev[i]*m2[i,j].adj
  #testZeigs()

  proc testSvdbd(n: int) =
    var v = newSeq[float](n)
    var d = newSeq[float](n)
    var e = newSeq[float](n)
    template `&`(x: seq): untyped = cast[ptr carray[float]](addr x[0])
    for i in 0..<n:
      d[i] = (i+1).float
      e[i] = -(i+1).float
    svdbi(&v, &d, &e, n)
    svdbi(&v, &d, &e, n)
    let t0 = epochTime()
    svdbi(&v, &d, &e, n)
    let t1 = epochTime()
    #for i in 0..<n:
    #  echo i, ": ", v[i]
    echo n, " time: ", t1-t0
  testSvdbd(100)
  testSvdbd(200)
  testSvdbd(400)
  testSvdbd(800)
  testSvdbd(1600)
  testSvdbd(3200)
  testSvdbd(6400)
  testSvdbd(12800)

  proc testSvdbdv =
    var d = newSeq[float](nr)
    var e = newSeq[float](nr)
    var v = newSeq[float](nr*nr)
    var u = newSeq[float](nr*nr)
    template `&`(x: seq): untyped = cast[ptr carray[float]](addr x[0])
    for i in 0..<nr:
      d[i] = (i+1).float
      e[i] = -(i+1).float
    svdBidiag(&d, &e, &v, &u, nr, nr)
    #for j in 0..<nr:
    #    v[i+nr*j]
  #testSvdbdv()
