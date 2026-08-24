## Geometry report: counts, weights, and the measured residual of every identity of
## doc/03-targets.md T1.1.  Rows are printed for L = 1, 2, 4, ... up to `lev`.

import base
import std/[algorithm, math, strformat]
import core/lattice

qexInit()
freezeTimers()

letParam:
  lev = 8
  nt = 32
  at = 0.2
  tilt = defaultTilt

installStandardParams()
echoParams()
processHelpParam()

func mdiff(m, n: Mat2): float =
  for i in 0..1:
    for j in 0..1: result = max(result, abs(m[i][j] - n[i][j]))

func dualExact(l, ls: float): float =
  ## Area of the two right spherical triangles with legs l/2 and l*:
  ## tan(E/2) = tan(a/2) tan(b/2), so E_pair = 4 arctan(tan(l/4) tan(l*/2)).
  ## Flat limit l l*/2.
  4.0*arctan(tan(0.25*l)*tan(0.5*ls))

echo "L    N_V    N_E    N_F  V-E+F        abar         a/at        maxM"
var levs: seq[int]
var lv = 1
while lv <= lev:
  levs.add lv
  lv *= 2
for l in levs:
  let
    s = newSphere(l, tilt)
    lat = newLat(s, nt, at)
  echo &"{l:<3d}{s.nv:6d} {s.ne:6d} {s.nf:6d} {s.nv-s.ne+s.nf:6d}  {s.abar:.10f}  {lat.asOverAt:11.6f} {lat.maxM:11.6f}"

echo ""
echo "L        min A_y      max A_y     mean A_y        min l        max l       min l*  #A_y vals  #l vals"
for l in levs:
  let s = newSphere(l, tilt)
  var
    amin = s.area[0]
    amax = s.area[0]
    asum = 0.0
    lmin = s.edges[0].len
    lmax = s.edges[0].len
    dmin = s.edges[0].dl[0]
  for a in s.area:
    amin = min(amin, a)
    amax = max(amax, a)
    asum += a
  for e in s.edges:
    lmin = min(lmin, e.len)
    lmax = max(lmax, e.len)
    dmin = min(dmin, min(e.dl[0], e.dl[1]))
  var av = s.area
  av.sort
  var lvv = newSeq[float](s.ne)
  for i in 0..<s.ne: lvv[i] = s.edges[i].len
  lvv.sort
  var na = 1
  for i in 1..<av.len:
    if av[i] - av[i-1] > 1e-12: inc na
  var nl = 1
  for i in 1..<lvv.len:
    if lvv[i] - lvv[i-1] > 1e-12: inc nl
  echo &"{l:<3d}{amin:13.9f}{amax:13.9f}{asum/float(s.nv):13.9f}{lmin:13.9f}{lmax:13.9f}{dmin:13.9f}{na:9d}{nl:9d}"

echo ""
echo "identities (max residual over the whole lattice)"
echo "L    |sumAtri-4pi| |sumAy-4pi|   dual-exact   sum(subA)    subA->Ay     holonomy    global-hol   omega-oracle   Z2  defect/2pi  nf"
for l in levs:
  let s = newSphere(l, tilt)
  var
    sa = 0.0
    sy = 0.0
    rex = 0.0
    rsub = 0.0
    rhol = 0.0
    rome = 0.0
  var gl = id2
  for f in 0..<s.nf:
    sa += s.faces[f].area
    var de = 0.0
    var ss = 0.0
    for i in 0..2:
      let e = s.faces[f].e[i]
      de += dualExact(s.edges[e].len, s.edges[e].dl[if s.faces[f].s[i] > 0: 0 else: 1])
      ss += s.faces[f].sub[i]
    rex = max(rex, abs(de - s.faces[f].area))
    rsub = max(rsub, abs(ss - s.faces[f].area))
    let h = s.holonomy f
    rhol = max(rhol, mdiff(h, expIsig3(s.faces[f].area)))
    gl = gl*h
  var ray = 0.0
  for y in 0..<s.nv:
    sy += s.area[y]
    var t = 0.0
    for k, f in s.nbf[y]:
      for i in 0..2:
        if s.faces[f].v[i] == y: t += s.faces[f].sub[i]
    ray = max(ray, abs(t - s.area[y]))
  var
    z2 = 0.0
    oc = newSeq[float](s.ne)
  for e in 0..<s.ne:
    oc[e] = s.omegaChart e
    let
      d = oc[e] - s.edges[e].omega
      k = round(d/(2.0*PI))
      sc = if (int(k) and 1) == 1: -1.0 else: 1.0
    rome = max(rome, abs(d - 2.0*PI*k))
    z2 = max(z2, abs(sc*s.edges[e].sgn - 1.0))    # 0 or 2; must be a vertex coboundary
  # (B.4): the chart lift alone is short by 2 pi on each face containing a chart pole.
  # sum_f sum_i s_i omega_i = 0 while sum_f A_tri = 4 pi, so the total defect is -4 pi
  # spread over exactly chi(S^2) = 2 faces.
  var
    dtot = 0.0
    ndef = 0
  for f in 0..<s.nf:
    var d = -s.faces[f].area
    for i in 0..2: d += float(s.faces[f].s[i])*oc[s.faces[f].e[i]]
    dtot += d
    if abs(d) > PI: inc ndef
  echo &"{l:<3d}{abs(sa-4.0*PI):13.3e}{abs(sy-4.0*PI):13.3e}{rex:13.3e}{rsub:13.3e}{ray:13.3e}{rhol:13.3e}{mdiff(gl,id2):13.3e}{rome:15.3e}{z2:7.1f}{dtot/(2.0*PI):10.6f}{ndef:5d}"

echo ""
echo "O(abar^2) residuals, relative to A_tri (must fall by ~4 per level doubling)"
echo "L     closure off  closure diag   flat sum(l l*/2)   poleGap site/link"
for l in levs:
  let
    s = newSphere(l, tilt)
    cl = s.checkClosure
    pg = s.poleGap
    at3 = 4.0*PI/float(s.nf)
  var rfl = 0.0
  for f in 0..<s.nf:
    var d = 0.0
    for i in 0..2:
      let e = s.faces[f].e[i]
      d += 0.5*s.edges[e].len*s.edges[e].dl[if s.faces[f].s[i] > 0: 0 else: 1]
    rfl = max(rfl, abs(d - s.faces[f].area))
  echo &"{l:<3d}{cl.offDiag/at3:13.3e}{cl.diag/at3:14.3e}{rfl/at3:19.3e}          {pg.site:.4f} {pg.link:.4f}"

echo ""
echo "flat equilateral limit (kappa -> 1/sqrt3 = 0.57735, kappa'/(abar/at) -> sqrt3/2 = 0.86603)"
echo "L    mean kappa    -1/sqrt3   mean kappa'/(a/at)   -sqrt3/2    Var(l)/abar^2   max |Ay/((sqrt3/2) ay^2)-1|    kappa range"
for l in levs:
  let
    s = newSphere(l, tilt)
    lat = newLat(s, nt, at)
  var
    km = 0.0
    kmin = lat.kap[0]
    kmax = lat.kap[0]
  for k in lat.kap:
    km += k
    kmin = min(kmin, k)
    kmax = max(kmax, k)
  km /= float(s.ne)
  var kt = 0.0
  for k in lat.kapT: kt += k
  kt /= float(s.nv)
  var vl = 0.0
  for e in s.edges: vl += (e.len - s.abar)*(e.len - s.abar)
  vl /= float(s.ne)*s.abar*s.abar
  var ry = 0.0
  for y in 0..<s.nv:
    var ay = 0.0
    for e in s.nbe[y]: ay += s.edges[e].len
    ay /= float(s.nbe[y].len)
    ry = max(ry, abs(s.area[y]/(0.5*sqrt(3.0)*ay*ay) - 1.0))
  echo &"{l:<3d}{km:12.6f}{km-1.0/sqrt(3.0):12.3e}{kt/lat.asOverAt:21.6f}{kt/lat.asOverAt-0.5*sqrt(3.0):12.3e}{vl:16.3e}{ry:26.3e}   {kmin:.5f}-{kmax:.5f}"

echo ""
echo "flat cross-check on a regular hexagonal patch of polar radius d (errors must go like d^2)"
echo "d           l*/a-1/(2sqrt3)  A_tri/a^2-sqrt3/4   A_y/a^2-sqrt3/2   kappa-1/sqrt3   kappa'/(a/at)-sqrt3/2"
for d in [1e-1, 1e-2, 1e-3]:
  let o: Vec3 = [0.0, 0.0, 1.0]
  var p, cc: array[6, Vec3]
  for k in 0..5:
    let ph = float(k)*PI/3.0
    p[k] = [sin(d)*cos(ph), sin(d)*sin(ph), cos(d)]
  for k in 0..5: cc[k] = circum(o, p[k], p[(k+1) mod 6])
  var a = 0.0
  for k in 0..5: a += geodesic(o, p[k]) + geodesic(p[k], p[(k+1) mod 6])
  a /= 12.0
  var ay = 0.0
  for k in 0..5: ay += sphArea(o, cc[k], cc[(k+1) mod 6])
  let
    r3 = sqrt 3.0
    ld = dualLen(cc[0], o, p[0]) + dualLen(cc[5], p[0], o)
  echo &"{d:<8.0e}{0.5*ld/a-0.5/r3:17.3e}{sphArea(o,p[0],p[1])/(a*a)-0.25*r3:19.3e}{ay/(a*a)-0.5*r3:18.3e}{ld/a-1.0/r3:16.3e}{ay/(a*a)-0.5*r3:24.3e}"

echo ""
echo "site coordination and abar"
for l in levs:
  let s = newSphere(l, tilt)
  var n5, n6, no = 0
  for y in 0..<s.nv:
    case s.nbr[y].len
    of 5: inc n5
    of 6: inc n6
    else: inc no
  echo &"L={l:<3d} 5-fold {n5:5d}  6-fold {n6:5d}  other {no:5d}   abar = {s.abar:.6f}"

processSaveParams()
writeParamFile()
qexFinalize()
