import strutils
import primme, primme/complex
import base, field, comms/qmp

export complex

type
  PP = primme_params or primme_svds_params

proc sumReal*[P:PP](sendBuf: pointer; recvBuf: pointer; count: ptr cint;
                    primme: ptr P; ierr: ptr cint) {.noconv.} =
  for i in 0..<count[]:
    asarray[float](recvBuf)[i] = asarray[float](sendBuf)[i]
  QMP_sum_double_array(cast[ptr cdouble](recvBuf), count[])
  ierr[] = 0

# WARNING: low level implementation details follow.
template convPrimmeArray(ff:Field, aa:ptr ccomplex[float], ss:int, body:untyped) {.dirty.} =
  const
    nc = ff[0].len
    vl = ff.V
    cl = 2*vl
  let
    n = nc*f.l.nEven div vl
    skip = nc*ss div vl
  var
    a = asarray[float]aa
    f = asarray[float]ff.s.data
  tfor i, 0..<n:
    let s = i + skip
    forO j, 0, vl.pred:
      body
proc toPrimmeArray*(f:Field, a:ptr ccomplex[float], skip:int = 0) =
  ## `skip` is in units of sites.
  convPrimmeArray(f,a,skip):
    a[cl*i+2*j] = f[cl*s+j]
    a[cl*i+2*j+1] = f[cl*s+j+vl]
proc fromPrimmeArray*(f:Field, a:ptr ccomplex[float], skip:int = 0) =
  ## `skip` is in units of sites.
  convPrimmeArray(f,a,skip):
    f[cl*s+j] = a[cl*i+2*j]
    f[cl*s+j+vl] = a[cl*i+2*j+1]

# WARNING: low level implementation details follow.
template convPrimmeArrayGauge(ff:seq[Field], aa:ptr PRIMME_COMPLEX_DOUBLE, ss:int, body:untyped) {.dirty.} =
  const
    nc = ff[0][0].nrows
    vl = ff[0].V
    cl = 2*vl
  let
    nd = ff.len
    n = nc*nc*f[0].l.nSites div vl
    skip = nc*nc*ss div vl
  for mu in 0..<nd:
    var
      a = asarray[float](asarray[float](aa)[mu*cl*n].addr)
      f = asarray[float]ff[mu].s.data
    tfor i, 0..<n:
      let s = i + skip
      forO j, 0, vl.pred:
        body
proc toPrimmeArrayGauge*(f:seq[Field], a:ptr PRIMME_COMPLEX_DOUBLE, skip:int = 0) =
  ## `skip` is in units of sites.
  convPrimmeArrayGauge(f,a,skip):
    a[cl*i+2*j] = f[cl*s+j]
    a[cl*i+2*j+1] = f[cl*s+j+vl]
proc fromPrimmeArrayGauge*(f:seq[Field], a:ptr PRIMME_COMPLEX_DOUBLE, skip:int = 0) =
  ## `skip` is in units of sites.
  convPrimmeArrayGauge(f,a,skip):
    f[cl*s+j] = a[cl*i+2*j]
    f[cl*s+j+vl] = a[cl*i+2*j+1]

template ff*(x:untyped):auto = formatFloat(x,ffScientific,17)
