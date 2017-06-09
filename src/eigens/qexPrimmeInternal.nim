import strutils
import primme
import base, field, comms/qmp

type
  OpInfo*[S,F] = object
    s*:ptr S
    m*:float
    x*,y*:F
proc sumReal*[P](sendBuf: pointer; recvBuf: pointer; count: ptr cint;
               primme: ptr P; ierr: ptr cint) {.noconv.} =
  for i in 0..<count[]:
    asarray[float](recvBuf)[i] = asarray[float](sendBuf)[i]
  QMP_sum_double_array(cast[ptr cdouble](recvBuf), count[])
  ierr[] = 0

# WARNING: low level implementation details follow.
template convPrimmeArray(ff:Field, aa:ptr float, ss:int, body:untyped) {.dirty.} =
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
    forO j, 0, <vl:
      body
proc toPrimmeArray*(f:Field, a:ptr float, skip:int = 0) =
  ## `skip` is in units of sites.
  convPrimmeArray(f,a,skip):
    a[cl*i+2*j] = f[cl*s+j]
    a[cl*i+2*j+1] = f[cl*s+j+vl]
proc fromPrimmeArray*(f:Field, a:ptr float, skip:int = 0) =
  ## `skip` is in units of sites.
  convPrimmeArray(f,a,skip):
    f[cl*s+j] = a[cl*i+2*j]
    f[cl*s+j+vl] = a[cl*i+2*j+1]

template ff*(x:untyped):auto = formatFloat(x,ffScientific,17)
