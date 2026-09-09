## Shared by the honeycomb test suites: PASS/FAIL reporting, scalar N x N
## complex algebra and site indexing for brute-force references, weak Abelian
## plane-wave and constant-flux backgrounds, decay-rate extraction.
## Importing this module sets the unittest console formatter.

import std/[math, complex, unittest]
import qex except epsilon
import physics/qcdTypes
import ../hcgeom, ../hcgauge, ../hcaction

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))

const nc* = getDefaultNc()

template ok*(msg: string, cond: bool) =
  ## a PASS/FAIL line carrying the measured numbers; a failure marks the
  ## enclosing unittest `test` failed
  if cond:
    echo "PASS: ", msg
  else:
    echo "FAIL: ", msg
    fail()

# ---------------------------------------------------------------------------
# single-site access (QEX element accessors return lane proxies, not floats)
# ---------------------------------------------------------------------------

proc cdot*(x, y: auto): Complex64 =
  ## <x|y>, first argument conjugated (QEX convention)
  let z = dot(x, y)
  complex64(z.re, z.im)

proc siteNorm2*(f: auto, i: int): float =
  ## |psi(i)|^2 over spin and colour
  for s in 0..<4:
    for c in 0..<nc:
      result += toF(f{i}[s][c].re)^2 + toF(f{i}[s][c].im)^2

proc spinor*(f: auto, i, c: int): array[4, Complex64] =
  ## the 4 spin components of colour c at site i
  for s in 0..<4:
    result[s] = complex64(toF f{i}[s][c].re, toF f{i}[s][c].im)

# ---------------------------------------------------------------------------
# scalar colour matrices, independent of the QEX matrix code
# ---------------------------------------------------------------------------

type Mat* = array[nc, array[nc, array[2, float]]]   ## [row][col][re/im]

proc mid*(): Mat =
  for i in 0..<nc: result[i][i][0] = 1.0

proc mmul*(a, b: Mat): Mat =
  for i in 0..<nc:
    for j in 0..<nc:
      var re = 0.0
      var im = 0.0
      for k in 0..<nc:
        re += a[i][k][0]*b[k][j][0] - a[i][k][1]*b[k][j][1]
        im += a[i][k][0]*b[k][j][1] + a[i][k][1]*b[k][j][0]
      result[i][j][0] = re
      result[i][j][1] = im

proc mdag*(a: Mat): Mat =
  for i in 0..<nc:
    for j in 0..<nc:
      result[i][j][0] = a[j][i][0]
      result[i][j][1] = -a[j][i][1]

proc madd*(a, b: Mat): Mat =
  for i in 0..<nc:
    for j in 0..<nc:
      result[i][j][0] = a[i][j][0] + b[i][j][0]
      result[i][j][1] = a[i][j][1] + b[i][j][1]

proc mscale*(a: Mat, s: float): Mat =
  for i in 0..<nc:
    for j in 0..<nc:
      result[i][j][0] = s*a[i][j][0]
      result[i][j][1] = s*a[i][j][1]

proc reTr*(a: Mat): float =
  for i in 0..<nc: result += a[i][i][0]

proc imTr*(a: Mat): float =
  for i in 0..<nc: result += a[i][i][1]

proc reTrMulAdj*(a, b: Mat): float =
  ## Re Tr(a b^dag)
  for i in 0..<nc:
    for j in 0..<nc:
      result += a[i][j][0]*b[i][j][0] + a[i][j][1]*b[i][j][1]

proc mTAH*(a: Mat): Mat =
  ## (a - a^dag)/2 - tr/nc
  let d = mdag(a)
  for i in 0..<nc:
    for j in 0..<nc:
      result[i][j][0] = 0.5*(a[i][j][0] - d[i][j][0])
      result[i][j][1] = 0.5*(a[i][j][1] - d[i][j][1])
  let tre = reTr(result)/nc.float
  let tim = imTr(result)/nc.float
  for i in 0..<nc:
    result[i][i][0] -= tre
    result[i][i][1] -= tim

proc mmaxdiff*(a, b: Mat): float =
  for i in 0..<nc:
    for j in 0..<nc:
      result = max(result, abs(a[i][j][0] - b[i][j][0]))
      result = max(result, abs(a[i][j][1] - b[i][j][1]))

proc getMat*(f: auto, idx: int): Mat =
  for a in 0..<nc:
    for b in 0..<nc:
      result[a][b][0] = toF f{idx}[a, b].re
      result[a][b][1] = toF f{idx}[a, b].im

proc setMat*(f: auto, idx: int, m: Mat) =
  for a in 0..<nc:
    for b in 0..<nc:
      f{idx}[a, b].re := m[a][b][0]
      f{idx}[a, b].im := m[a][b][1]

# ---------------------------------------------------------------------------
# site indexing on the cell torus
# ---------------------------------------------------------------------------

proc wrapCoord*(c: Cell, geom: openArray[int]): Cell =
  for mu in 0..<nDim:
    result[mu] = ((c[mu] mod geom[mu]) + geom[mu]) mod geom[mu]

proc siteIndex*(lo: Layout, c: Cell, geom: openArray[int]): int =
  ## local site index of the wrapped cell c; single-rank layouts
  let ri = lo.rankIndex(wrapCoord(c, geom))
  doAssert ri.rank == lo.myRank, "single-rank tests only"
  ri.index

iterator lexCells*(geom: openArray[int]): Cell =
  for x3 in 0..<geom[3]:
    for x2 in 0..<geom[2]:
      for x1 in 0..<geom[1]:
        for x0 in 0..<geom[0]:
          yield [x0, x1, x2, x3]

# ---------------------------------------------------------------------------
# gauge-field arithmetic on all 24 link fields
# ---------------------------------------------------------------------------

proc wilsonAction*(gact: GaugeActionCoeffs, g: auto): float =
  ## S = beta sum_{x,mu<nu} (1 - Re Tr P/N); gaugeAction1 omits the constant
  gact.gaugeAction1(g) + gact.plaq*6.0*float(g[0].l.physVol)

# ---------------------------------------------------------------------------
# Abelian backgrounds embedded via T = diag(1,-1,0), U = exp(i phi T) with
# phi the exact straight-line integral of A along the link
# ---------------------------------------------------------------------------

template setPhase*(m: untyped, th: float) =
  ## m := exp(i th T); m must be the identity
  m[0, 0].re := cos(th)
  m[0, 0].im := sin(th)
  m[1, 1].re := cos(th)
  m[1, 1].im := -sin(th)

proc phaseIntegral*(a, b: float): float =
  ## int_0^1 cos(a + t b) dt
  if abs(b) < 1e-12: cos(a)
  else: (sin(a+b) - sin(a))/b

proc linkPhase*(x, n, eps, pv: array[4, float]): float =
  ## int_x^{x+n} A.dl  for  A_mu(x) = eps_mu cos(p.x)
  var a = 0.0
  var b = 0.0
  var en = 0.0
  for mu in 0..<4:
    a += pv[mu]*x[mu]
    b += pv[mu]*n[mu]
    en += eps[mu]*n[mu]
  en*phaseIntegral(a, b)

proc fluxPhase*(x, n: array[4, float]; f1, f2: float; l0, l2: int): float =
  ## int_x^{x+n} A.dl for A_1 = f1 x0, A_3 = f2 x2 in the fundamental domain,
  ## with the transition-function corrections -f1 L0 x1(end), -f2 L2 x3(end)
  ## for links crossing the x0 = L0 / x2 = L2 seams (diagonal links cross
  ## with nonzero transverse displacement, hence the endpoint rule)
  result = n[1]*f1*(x[0] + 0.5*n[0]) + n[3]*f2*(x[2] + 0.5*n[2])
  if x[0] + n[0] >= float(l0) - 1e-9:
    result -= f1*float(l0)*(x[1] + n[1])
  if x[2] + n[2] >= float(l2) - 1e-9:
    result -= f2*float(l2)*(x[3] + n[3])

proc setHcPhases(hg: HcGauge, phase: proc(x, n: array[4, float]): float) =
  ## every link of hg from `phase(x, n)`, x the start point, n the unit step
  hg.unit
  let lo = hg.lo
  for i in lo.sites:
    var y, yb: array[4, float]
    for mu in 0..<4:
      y[mu] = lo.coords[mu][i].float    # A(y) at integer coords
      yb[mu] = y[mu] + 0.5              # B(y) at half-integer coords
    for mu in 0..<4:
      var n: array[4, float]
      n[mu] = 1.0
      setPhase(hg.uA[mu]{i}, phase(y, n))
      setPhase(hg.uB[mu]{i}, phase(yb, n))
    for d in 0..<nDiag:
      var n: array[4, float]
      for mu in 0..<4:
        n[mu] = float((d shr mu) and 1) - 0.5
      setPhase(hg.uD[d]{i}, phase(yb, n))

proc setAbelianHc*(hg: HcGauge, eps, pv: array[4, float]) =
  ## plane wave A_mu(x) = eps_mu cos(p.x) on the honeycomb
  setHcPhases(hg, proc(x, n: array[4, float]): float = linkPhase(x, n, eps, pv))

proc setAbelianCubic*(g: auto, lo: auto, eps, pv: array[4, float]) =
  ## the same plane wave on the cubic lattice
  g.unit
  for i in lo.sites:
    var y: array[4, float]
    for mu in 0..<4:
      y[mu] = lo.coords[mu][i].float
    for mu in 0..<4:
      var n: array[4, float]
      n[mu] = 1.0
      setPhase(g[mu]{i}, linkPhase(y, n, eps, pv))

proc setFluxHc*(hg: HcGauge, n1, n2: int) =
  ## constant field strength F_01 = 2 pi n1/(L0 L1), F_23 = 2 pi n2/(L2 L3)
  let
    geom = hg.lo.physGeom
    f1 = 2.0*PI*float(n1)/float(geom[0]*geom[1])
    f2 = 2.0*PI*float(n2)/float(geom[2]*geom[3])
  setHcPhases(hg, proc(x, n: array[4, float]): float =
    fluxPhase(x, n, f1, f2, geom[0], geom[2]))

proc setFluxCubic*(g: auto, lo: auto, n1, n2: int) =
  ## the same fluxes on the cubic lattice (refCubicMeas -abeliantest)
  let lat = lo.physGeom
  let
    f1 = 2.0*PI*float(n1)/float(lat[0]*lat[1])
    f2 = 2.0*PI*float(n2)/float(lat[2]*lat[3])
  g.unit
  for i in lo.sites:
    let
      x0 = lo.coords[0][i]
      x1 = lo.coords[1][i]
      x2 = lo.coords[2][i]
      x3 = lo.coords[3][i]
    setPhase(g[1]{i}, f1*x0.float)
    if x0 == lat[0]-1:
      setPhase(g[0]{i}, -f1*float(lat[0]*x1))
    setPhase(g[3]{i}, f2*x2.float)
    if x2 == lat[2]-1:
      setPhase(g[2]{i}, -f2*float(lat[2]*x3))

const
  epsAmp* = 3.0e-3
    ## weak-field amplitude for the heat-kernel rate measurements: S(t) in the
    ## late windows stays >~1e5 above the cancellation floor of the action sum
    ## (at 1e-3 the rate moved by ~5e-5 rel between thread counts), while the
    ## nonlinear corrections O(epsAmp^2) ~ 1e-5 sit far below the 1e-3 gates
  ehat* = [0.0, 1.0, -0.6, 0.3]
    ## transverse polarisation for p along direction 0

# ---------------------------------------------------------------------------
# decay rates and spectra
# ---------------------------------------------------------------------------

proc plateauRate*(ss: seq[float], dt: float, frac = 0.25):
    tuple[rate, drift: float] =
  ## rate_i = ln(S_{i-1}/S_i)/(2 dt) averaged over the last `frac` of the
  ## series; `drift` is the max deviation inside that window
  var r: seq[float]
  for i in 1..<ss.len:
    r.add ln(ss[i-1]/ss[i])/(2.0*dt)
  let m = max(2, int(frac*r.len.float))
  var s = 0.0
  for i in r.len-m..<r.len: s += r[i]
  let avg = s/m.float
  var d = 0.0
  for i in r.len-m..<r.len: d = max(d, abs(r[i]-avg))
  (avg, d)

proc fit3*(x, y: array[3, float]): array[3, float] =
  ## exact y = c0 + c1 x + c2 x^2 through 3 points
  let
    d = (x[1]-x[0])*(x[2]-x[0])*(x[2]-x[1])
    c2 = ((y[2]-y[0])*(x[1]-x[0]) - (y[1]-y[0])*(x[2]-x[0]))/d
    c1 = (y[1]-y[0])/(x[1]-x[0]) - c2*(x[1]+x[0])
    c0 = y[0] - c1*x[0] - c2*x[0]*x[0]
  [c0, c1, c2]

proc worstNearest*(vals, exact: openArray[Complex64]): float =
  ## worst over vals of the distance to the nearest exact value; exact non-empty
  for v in vals:
    var d = abs(v - exact[0])
    for e in exact: d = min(d, abs(v - e))
    result = max(result, d)
