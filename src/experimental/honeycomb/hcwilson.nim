## Wilson-Dirac operator on the 16-cell honeycomb, with optional tree-level
## clover term (a = 1, 24 unit neighbour vectors n_i, (gamma.n_i)^2 = 1):
##
##   D psi(x) = (m + 4 r) psi(x) + (1/6) sum_{i=1}^{24} (gamma.n_i - r) U_i(x) psi(x+n_i)
##            - (c_SW r/2) sum_{a>b} gamma_a gamma_b Fhat_ab(x) psi(x)
##
## Free: D(p) = M + i gamma.K, M = (r/6) sum_i (1 - cos p.n_i) -> r p^2/2,
## K_mu = (1/6) sum_i n_i,mu sin p.n_i -> p_mu.  The Wilson term is -(r/2) D^2
## as on the cubic lattice, so c_SW = 1 is the standard tree-level value;
## Fhat_ab = +i a^2 F_ab T from hctopo (sigma_ab F_ab = gamma_a gamma_b Fhat_ab).
## D^dag is the same hopping sum with the spin factor of the opposite
## direction; gamma5 D gamma5 = D^dag is tested, not built in.
##
## Fields: one Dirac field per sublattice on the cell layout, psiA(y) at A(y),
## psiB(y) at B(y+1/2).  Hops per cell y:
##   A row: uA[mu](y) psiA(y+e_mu); uA[mu](y-e_mu)^dag psiA(y-e_mu);
##          uD[dbar](y-dbar)^dag psiB(y-dbar), dbar = delta xor 15, direction d(delta)
##   B row: same with uB; uD[delta](y) psiA(y+delta), direction d(delta)
## Plane waves are phased with the integer cell coordinate on both sublattices
## (no half-site phase); antiperiodic time shifts k_3 by pi/N_t.
## Gammas: QEX DeGrand-Rossi gamma1..gamma4 for mu = 0..3, gamma5 = diag(1,1,-1,-1).

import base, layout, field, maths
import physics/qcdTypes
import rng
import hcgeom, hcgauge, hctopo

export hcgeom, hcgauge, hctopo

type
  HcFermion*[F] = object
    a*, b*: F

  SpinMat* = typeof(gamma0)

proc newDiracField*(lo: Layout, nc: static[int] = getDefaultNc()): auto =
  ## one Dirac field on the cell layout, any SIMD length
  type C = typeof(lo.newDComplexV)
  type DF = Spin[VectorArray[4, Color[VectorArray[nc, C]]]]
  lo.newField(DF)

proc newHcFermion*(lo: Layout): auto =
  ## zeroed; allocates
  type F = typeof(newDiracField(lo))
  var r: HcFermion[F]
  r.a = newDiracField(lo)
  r.b = newDiracField(lo)
  r.a := 0
  r.b := 0
  r

proc newOneOf*[F](x: HcFermion[F]): HcFermion[F] =
  result.a = x.a.newOneOf
  result.b = x.b.newOneOf
  result.a := 0
  result.b := 0

# allocation free from here on

proc `:=`*[F](r: HcFermion[F], x: HcFermion[F]) =
  r.a := x.a
  r.b := x.b

proc `:=`*[F](r: HcFermion[F], v: SomeNumber) =
  r.a := v
  r.b := v

proc gaussian*(x: HcFermion, r: var RNGField) =
  x.a.gaussian r
  x.b.gaussian r

proc norm2*(x: HcFermion): float =
  norm2(x.a) + norm2(x.b)

proc dot*(x, y: HcFermion): auto =
  ## QEX convention: the first argument is conjugated
  dot(x.a, y.a) + dot(x.b, y.b)

proc norm2diff*(x, y: HcFermion): float =
  norm2diff(x.a, y.a) + norm2diff(x.b, y.b)

proc applyGamma5*(r: HcFermion, x: HcFermion) =
  for e in r.a:
    r.a[e] := gamma5 * x.a[e]
  for e in r.b:
    r.b[e] := gamma5 * x.b[e]

template mul*(r: var Spin, x: Color, y: Spin2) =
  ## colour matrix times Dirac fermion; needed by Transporter, absent in QEX
  mixin mul
  mul(r[], x, y[])

# ---------------------------------------------------------------------------
# spin matrices
# ---------------------------------------------------------------------------

proc gammaDotDir*(dir: int): SpinMat =
  ## gamma.n for direction index `dir` of hcgeom
  let n = toFloat(dirVec(dir))
  var t: SpinMat
  t := n[0]*gamma1 + n[1]*gamma2 + n[2]*gamma3 + n[3]*gamma4
  t

proc hopMat*(dir: int, rw: float): SpinMat =
  ## (gamma.n_dir - rw)/6
  var t: SpinMat
  t := gammaDotDir(dir) - rw*gamma0
  var s: SpinMat
  s := (1.0/6.0)*t
  s

proc buildCloverGam(): array[6, SpinMat] =
  var t: SpinMat
  template setp(a, b: int, ga, gb: typed) =
    t := ga*gb
    result[pairIndex(a, b)] = t
  setp(1, 0, gamma2, gamma1)
  setp(2, 0, gamma3, gamma1)
  setp(2, 1, gamma3, gamma2)
  setp(3, 0, gamma4, gamma1)
  setp(3, 1, gamma4, gamma2)
  setp(3, 2, gamma4, gamma3)

let cloverGam* = buildCloverGam()
  ## cloverGam[pairIndex(a,b)] = gamma_a gamma_b

# ---------------------------------------------------------------------------
# the operator
# ---------------------------------------------------------------------------

type
  HcWilson*[MF, FF, TR, SHF, SHM, TW] = ref object
    ## gauge ref, pre-shifted diagonal links, shifters, hop matrices and the
    ## clover field.  Call `gaugeRefresh` after the links change in place.
    g*: HcGauge[MF]
    uDsh*: array[nDiag, MF]      ## uDsh[d](y) = uD[d](y-d)
    tfA, tbA, tfB, tbB: array[nDim, TR]
    shAf: HcShift16[FF, SHF]     ## f[d](y) = psiA(y+d)
    shBb: HcShift16[FF, SHF]     ## f[d](y) = psiB(y-d)
    shDb: array[nDim, SHM]
    gm: array[nDirs, SpinMat]    ## (gamma.n_dir - rw)/6
    rwCur: float
    rwValid: bool
    cSW*: float
    tw*: TW                      ## TopoWork with Fhat (nil if cSW = 0)

proc setHopMats(w: HcWilson, rw: float) =
  if w.rwValid and w.rwCur == rw: return
  for dir in 0..<nDirs:
    w.gm[dir] = hopMat(dir, rw)
  w.rwCur = rw
  w.rwValid = true

proc gaugeRefresh*(w: HcWilson) =
  ## rebuild uDsh from w.g (32 single-axis shifts) and the clover field;
  ## outside `threads:`
  threads:
    w.uDsh[0] := w.g.uD[0]
    for d in 1..<nDiag:
      var cur = w.g.uD[d]
      for mu in 0..<nDim:
        if ((d shr mu) and 1) == 1:
          cur = w.shDb[mu] ^* cur
      w.uDsh[d] := cur
  if w.cSW != 0.0:
    fmunu(w.tw, w.g)

proc newHcWilson*[MF](g: HcGauge[MF], cSW = 0.0): auto =
  ## allocates; outside `threads:`
  var proto = newDiracField(g.lo)
  type FF = typeof(proto)
  type TR = typeof(newTransporter(g.uA[0], proto, 0, 1))
  type SHF = typeof(newShifter(proto, 0, 1))
  type SHM = typeof(newShifter(g.uA[0], 0, 1))
  type TW = typeof(newTopoWork(g))
  var w = HcWilson[MF, FF, TR, SHF, SHM, TW]()
  w.g = g
  w.cSW = cSW
  for d in 0..<nDiag:
    w.uDsh[d] = g.uD[d].newOneOf
  for mu in 0..<nDim:
    w.tfA[mu] = newTransporter(g.uA[mu], proto, mu, 1)
    w.tbA[mu] = newTransporter(g.uA[mu], proto, mu, -1)
    w.tfB[mu] = newTransporter(g.uB[mu], proto, mu, 1)
    w.tbB[mu] = newTransporter(g.uB[mu], proto, mu, -1)
    w.shDb[mu] = newShifter(g.uA[0], mu, -1)
  w.shAf = newHcShift16(proto, 1)
  w.shBb = newHcShift16(proto, -1)
  if cSW != 0.0:
    w.tw = newTopoWork(g)
  w.gaugeRefresh
  w

proc applyDirac(w: HcWilson, r: HcFermion, x: HcFermion,
                m, rw: float, dag: bool) =
  ## r = D x (dag = false) or D^dag x; opens `threads:`; r must not alias x
  doAssert not (r.a == x.a or r.b == x.b), "r must not alias x"
  w.setHopMats(rw)
  var cf: array[nDirs, int]
  for i in 0..<nDirs:
    cf[i] = if dag: opposite(i) else: i
  w.shAf.setSrc x.a
  w.shBb.setSrc x.b
  let mass = m + 4.0*rw
  var gc {.noinit.}: array[6, SpinMat]
  let useClover = w.cSW != 0.0
  if useClover:
    for p in 0..<6:
      gc[p] := (-0.5*w.cSW*rw)*cloverGam[p]
  threads:
    w.shAf.run
    w.shBb.run
    r.a := mass*x.a
    r.b := mass*x.b
    for mu in 0..<nDim:
      block:
        let gmm = w.gm[cf[axisIndex(mu, false)]]
        let h = w.tfA[mu] ^* x.a
        for e in r.a:
          r.a[e] += gmm * h[e]
      block:
        let gmm = w.gm[cf[axisIndex(mu, true)]]
        let h = w.tbA[mu] ^* x.a
        for e in r.a:
          r.a[e] += gmm * h[e]
      block:
        let gmm = w.gm[cf[axisIndex(mu, false)]]
        let h = w.tfB[mu] ^* x.b
        for e in r.b:
          r.b[e] += gmm * h[e]
      block:
        let gmm = w.gm[cf[axisIndex(mu, true)]]
        let h = w.tbB[mu] ^* x.b
        for e in r.b:
          r.b[e] += gmm * h[e]
    for d in 0..<nDiag:
      let gmm = w.gm[cf[diagIndex(d)]]
      for e in r.b:
        r.b[e] += gmm * (w.g.uD[d][e] * w.shAf.f[d][e])
      let dl = d xor 15
      for e in r.a:
        r.a[e] += gmm * (w.uDsh[dl][e].adj * w.shBb.f[dl][e])
    if useClover:
      for p in 0..<6:
        for e in r.a:
          r.a[e] += gc[p] * (w.tw.f[0][p][e] * x.a[e])
        for e in r.b:
          r.b[e] += gc[p] * (w.tw.f[1][p][e] * x.b[e])

proc D*(w: HcWilson, r: var HcFermion, x: HcFermion,
        m: float, rw: float = 1.0) =
  applyDirac(w, r, x, m, rw, false)

proc Ddag*(w: HcWilson, r: var HcFermion, x: HcFermion,
           m: float, rw: float = 1.0) =
  ## the clover term is Hermitian and commutes with gamma5: unchanged in D^dag
  applyDirac(w, r, x, m, rw, true)

proc setBC*(g: HcGauge) =
  ## antiperiodic time: flip uA[3], uB[3] and the 8 uD with bit 3 set at
  ## cell time N_t-1; inside `threads:`.  Refresh any HcWilson afterwards.
  template flip(u: untyped) =
    let nt1 = u.l.physGeom[3] - 1
    for i in u.l.sites:
      if u.l.coords[3][i] == nt1:
        u{i} *= -1
  flip g.uA[3]
  flip g.uB[3]
  for d in 0..<nDiag:
    if ((d shr 3) and 1) == 1:
      flip g.uD[d]
