## Lattice realizations of the analytic maps and their graph actions.

import base/globals
setDefaultNc(1)
setVLENmax(4)

import qex
import math

import maps
export maps

import ../graph/core
import ../graph/scalar
import ../graph/gauge
import ../graph/gauge/shared
import ../graph/gauge/basic_ops
import ../graph/gauge/action/ops
import ../graph/support/op as opsupport
import layout, physics/qcdTypes

proc mapDelta(m: CircleMap; p: auto): auto =
  result = p
  for k in 0..<simdLength(p):
    result[k] = circleEval(m, p[k]).y-p[k]

proc mapLnGp(m: CircleMap; p: auto): auto =
  result = p
  for k in 0..<simdLength(p):
    result[k] = ln(circleEval(m, p[k]).dy)

proc mapDlogj(m: CircleMap; p: auto): auto =
  result = p
  for k in 0..<simdLength(p):
    let e = circleEval(m, p[k])
    result[k] = e.ddy/e.dy

proc mapDeltaGpm1(m: CircleMap; p: auto): tuple[delta, gpm1: typeof(p)] =
  result.delta = p
  result.gpm1 = p
  for k in 0..<simdLength(p):
    let e = circleEval(m, p[k])
    result.delta[k] = e.y-p[k]
    result.gpm1[k] = e.dy-1.0

type PortalMask* = DLatticeColorMatrixV

proc portalMasksEvenOdd*(proto: PortalMask): seq[PortalMask] =
  let
    lo = proto.l
    nd = lo.nDim
    me = proto.newOneOf
    mo = proto.newOneOf
  doAssert lo.physGeom[0] mod 2 == 0 and lo.physGeom[1] mod 2 == 0, "even/odd portals require even lattice extents"
  threads:
    var co = newSeq[cint](nd)
    for j in lo.sites:
      lo.coord(co, (lo.myRank, j))
      let even = ((co[0].int+co[1].int) and 1) == 0
      me{j} := (if even: 1.0 else: 0.0)
      mo{j} := (if even: 0.0 else: 1.0)
  @[me, mo]

proc portalMasksStride*(proto: PortalMask; stride, offx, offy: int): seq[PortalMask] =
  let
    lo = proto.l
    lat = lo.physGeom
    nd = lo.nDim
  doAssert stride >= 2, "portal stride must be at least two"
  doAssert lat[0] mod stride == 0 and lat[1] mod stride == 0, "portal stride must divide each lattice extent"
  let
    ox = ((offx mod stride)+stride) mod stride
    oy = ((offy mod stride)+stride) mod stride
    ms = proto.newOneOf
  threads:
    var co = newSeq[cint](nd)
    for j in lo.sites:
      lo.coord(co, (lo.myRank, j))
      let active = co[0].int mod stride == ox and co[1].int mod stride == oy
      ms{j} := (if active: 1.0 else: 0.0)
  @[ms]

template fillPlaqAngle(pa, sf, g: untyped) =
  block:
    let p = (g[0]*(sf[0] ^* g[1]))*(g[1]*(sf[1] ^* g[0])).adj
    for x in pa:
      let z = p[x][0, 0]
      pa[x][0, 0].re := atan2(z.im, z.re)

proc plaqAngleField*(g: shared.Gauge): auto =
  let
    sf = newShifters(g[0], 1)
    p = g[0].newOneOf
  threads:
    fillPlaqAngle(p, sf, g)
  p

proc portalSmearGrad(V, cU: Ggauge; mask: PortalMask; m: CircleMap): Ggauge =
  let
    pa = V.gval[0].newOneOf
    D = V.gval[0].newOneOf
    GP = V.gval[0].newOneOf
    U0 = V.gval[0].newOneOf
    U1 = V.gval[0].newOneOf
    R0 = V.gval[0].newOneOf
    R1 = V.gval[0].newOneOf
    B = V.gval[0].newOneOf
    sf = newShifters(V.gval[0], 1)
    sbD = newShifters(V.gval[0], -1)
    sbB = newShifters(V.gval[0], -1)
  proc kf(v: Gvalue) =
    let
      V = Ggauge(v.inputs[0])
      cU = Ggauge(v.inputs[1])
      cV = Ggauge(v)
    threads:
      fillPlaqAngle(pa, sf, V.gval)
      threadBarrier()
      for x in D:
        let
          p = pa[x][0, 0].re
          mk = mask[x][0, 0].re
          e = mapDeltaGpm1(m, p)
        D[x][0, 0].re := mk*e.delta
        GP[x][0, 0].re := mk*e.gpm1
      threadBarrier()
      let
        Dm0 = sbD[0] ^* D
        Dm1 = sbD[1] ^* D
      threadBarrier()
      for z in U0:
        let
          dt0 = 0.25*(D[z][0, 0].re-Dm1[z][0, 0].re)
          dt1 = 0.25*(Dm0[z][0, 0].re-D[z][0, 0].re)
          c0 = cos(dt0)
          s0 = sin(dt0)
          c1 = cos(dt1)
          s1 = sin(dt1)
          v0r = V.gval[0][z][0, 0].re
          v0i = V.gval[0][z][0, 0].im
          v1r = V.gval[1][z][0, 0].re
          v1i = V.gval[1][z][0, 0].im
          u0r = c0*v0r-s0*v0i
          u0i = c0*v0i+s0*v0r
          u1r = c1*v1r-s1*v1i
          u1i = c1*v1i+s1*v1r
        U0[z][0, 0].re := u0r
        U0[z][0, 0].im := u0i
        U1[z][0, 0].re := u1r
        U1[z][0, 0].im := u1i
        let
          cu0r = cU.gval[0][z][0, 0].re
          cu0i = cU.gval[0][z][0, 0].im
          cu1r = cU.gval[1][z][0, 0].re
          cu1i = cU.gval[1][z][0, 0].im
        R0[z][0, 0].re := cu0i*u0r-cu0r*u0i
        R1[z][0, 0].re := cu1i*u1r-cu1r*u1i
      threadBarrier()
      let
        R1p0 = sf[0] ^* R1
        R0p1 = sf[1] ^* R0
      threadBarrier()
      for x in B:
        let a = 0.25*(R0[x][0, 0].re+R1p0[x][0, 0].re-
          R0p1[x][0, 0].re-R1[x][0, 0].re)
        B[x][0, 0].re := a*GP[x][0, 0].re
      threadBarrier()
      let
        Bm0 = sbB[0] ^* B
        Bm1 = sbB[1] ^* B
      threadBarrier()
      for z in cV.gval[0]:
        let
          v0r = V.gval[0][z][0, 0].re
          v0i = V.gval[0][z][0, 0].im
          v1r = V.gval[1][z][0, 0].re
          v1i = V.gval[1][z][0, 0].im
          u0r = U0[z][0, 0].re
          u0i = U0[z][0, 0].im
          u1r = U1[z][0, 0].re
          u1i = U1[z][0, 0].im
          cu0r = cU.gval[0][z][0, 0].re
          cu0i = cU.gval[0][z][0, 0].im
          cu1r = cU.gval[1][z][0, 0].re
          cu1i = cU.gval[1][z][0, 0].im
          p0r = u0r*v0r+u0i*v0i
          p0i = u0r*v0i-u0i*v0r
          p1r = u1r*v1r+u1i*v1i
          p1i = u1r*v1i-u1i*v1r
          k0 = B[z][0, 0].re-Bm1[z][0, 0].re
          k1 = Bm0[z][0, 0].re-B[z][0, 0].re
        cV.gval[0][z][0, 0].re := p0r*cu0r-p0i*cu0i-k0*v0i
        cV.gval[0][z][0, 0].im := p0r*cu0i+p0i*cu0r+k0*v0r
        cV.gval[1][z][0, 0].re := p1r*cu1r-p1i*cu1i-k1*v1i
        cV.gval[1][z][0, 0].im := p1r*cu1i+p1i*cu1r+k1*v1r
  proc kb(zb, z: Gvalue; i: int; input: Gvalue): Gvalue =
    raiseUnsupportedPath("portalSmearGrad backward", "portal field-map gradient is not differentiated")
  let f = Gfunc(forward: kf, backward: kb, name: "portalSmearGrad")
  graphNode(V.gaugeNodeLike, @[Gvalue(V), Gvalue(cU)], f, "portalSmearGrad")

proc portalLogDetJ(V: Ggauge; mask: PortalMask; m: CircleMap): Gscalar

proc portalSmear*(V: Ggauge; mask: PortalMask; m: CircleMap): Ggauge =
  let
    pa = V.gval[0].newOneOf
    D = V.gval[0].newOneOf
    sf = newShifters(V.gval[0], 1)
    sb = newShifters(V.gval[0], -1)
  proc forward(v: Gvalue) =
    let
      V = Ggauge(v.inputs[0])
      U = Ggauge(v)
    threads:
      fillPlaqAngle(pa, sf, V.gval)
      threadBarrier()
      for x in D:
        D[x][0, 0].re := mask[x][0, 0].re*mapDelta(m, pa[x][0, 0].re)
      threadBarrier()
      let
        Dm0 = sb[0] ^* D
        Dm1 = sb[1] ^* D
      threadBarrier()
      for z in U.gval[0]:
        let
          dt0 = 0.25*(D[z][0, 0].re-Dm1[z][0, 0].re)
          dt1 = 0.25*(Dm0[z][0, 0].re-D[z][0, 0].re)
          c0 = cos(dt0)
          s0 = sin(dt0)
          c1 = cos(dt1)
          s1 = sin(dt1)
          v0r = V.gval[0][z][0, 0].re
          v0i = V.gval[0][z][0, 0].im
          v1r = V.gval[1][z][0, 0].re
          v1i = V.gval[1][z][0, 0].im
        U.gval[0][z][0, 0].re := c0*v0r-s0*v0i
        U.gval[0][z][0, 0].im := c0*v0i+s0*v0r
        U.gval[1][z][0, 0].re := c1*v1r-s1*v1i
        U.gval[1][z][0, 0].im := c1*v1i+s1*v1r
  proc backward(zb, z: Gvalue; i: int; input: Gvalue): Gvalue =
    let cU = requireUpstream(zb, "portalSmear backward", Ggauge)
    Gvalue(portalSmearGrad(Ggauge(z.inputs[0]), cU, mask, m))
  proc ldj(z: Gvalue): tuple[ld, via: Gvalue] =
    (Gvalue(portalLogDetJ(Ggauge(z.inputs[0]), mask, m)), z.inputs[0])
  let f = Gfunc(forward: forward, backward: backward, logdet: ldj, name: "portalSmear")
  graphNode(V.gaugeNodeLike, @[Gvalue(V)], f, "portalSmear")

proc portalLogDetJ(V: Ggauge; mask: PortalMask; m: CircleMap): Gscalar =
  let
    sf = newShifters(V.gval[0], 1)
    pf = V.gval[0].newOneOf
  proc forward(v: Gvalue) =
    let
      V = Ggauge(v.inputs[0])
      z = Gscalar(v)
    var res = 0.0
    threads:
      fillPlaqAngle(pf, sf, V.gval)
      threadBarrier()
      var s = 0.0
      for x in pf:
        s += simdSum(mask[x][0, 0].re*mapLnGp(m, pf[x][0, 0].re))
      s.threadRankSum
      threadSingle: res = s
    z.sval = res
  proc gradKernel(V: Ggauge): Ggauge =
    let
      pa = V.gval[0].newOneOf
      D = V.gval[0].newOneOf
      sf = newShifters(V.gval[0], 1)
      sb = newShifters(V.gval[0], -1)
    proc kf(v: Gvalue) =
      let
        V = Ggauge(v.inputs[0])
        cV = Ggauge(v)
      threads:
        fillPlaqAngle(pa, sf, V.gval)
        threadBarrier()
        for x in D:
          D[x][0, 0].re := mask[x][0, 0].re*mapDlogj(m, pa[x][0, 0].re)
        threadBarrier()
        let
          Dm0 = sb[0] ^* D
          Dm1 = sb[1] ^* D
        threadBarrier()
        for z in cV.gval[0]:
          let
            v0r = V.gval[0][z][0, 0].re
            v0i = V.gval[0][z][0, 0].im
            v1r = V.gval[1][z][0, 0].re
            v1i = V.gval[1][z][0, 0].im
            k0 = D[z][0, 0].re-Dm1[z][0, 0].re
            k1 = Dm0[z][0, 0].re-D[z][0, 0].re
          cV.gval[0][z][0, 0].re := -k0*v0i
          cV.gval[0][z][0, 0].im := k0*v0r
          cV.gval[1][z][0, 0].re := -k1*v1i
          cV.gval[1][z][0, 0].im := k1*v1r
    proc kb(zb, z: Gvalue; i: int; input: Gvalue): Gvalue =
      raiseUnsupportedPath("portalLogDetJ gradient backward", "portal log-Jacobian gradient is not differentiated")
    let f = Gfunc(forward: kf, backward: kb, name: "portalLogDetJgrad")
    graphNode(V.gaugeNodeLike, @[Gvalue(V)], f, "portalLogDetJgrad")
  proc backward(zb, z: Gvalue; i: int; input: Gvalue): Gvalue =
    scaledUpstreamOr(zb, Gscalar, gradKernel(Ggauge(z.inputs[0])))
  let f = Gfunc(forward: forward, backward: backward, name: "portalLogDetJ")
  graphNode(scalarNodeLike(V), @[Gvalue(V)], f, "portalLogDetJ")

proc smearFlow*(V: Ggauge; m: CircleMap; masks: seq[PortalMask]; sweeps: int): Ggauge =
  doAssert sweeps >= 0, "portal sweep count must be nonnegative"
  result = V
  for _ in 0..<sweeps:
    for mask in masks:
      result = portalSmear(result, mask, m)

proc portalSmearHost*(V: shared.Gauge; mask: PortalMask; m: CircleMap): tuple[u: shared.Gauge, lndet: float] =
  let
    pa = plaqAngleField(V)
    D = V[0].newOneOf
    u0 = V[0].newOneOf
    u1 = V[1].newOneOf
    sb = newShifters(V[0], -1)
  var ld = 0.0
  threads:
    var s = 0.0
    for x in D:
      var d = pa[x][0, 0].re
      for k in 0..<simdLength(d):
        if mask[x][0, 0].re[k] != 0.0:
          let e = circleEval(m, pa[x][0, 0].re[k])
          d[k] = e.y-pa[x][0, 0].re[k]
          s += ln(e.dy)
        else:
          d[k] = 0.0
      D[x][0, 0].re := d
    s.threadRankSum
    threadSingle: ld = s
    threadBarrier()
    let
      Dm0 = sb[0] ^* D
      Dm1 = sb[1] ^* D
    threadBarrier()
    for z in u0:
      let
        dt0 = 0.25*(D[z][0, 0].re-Dm1[z][0, 0].re)
        dt1 = 0.25*(Dm0[z][0, 0].re-D[z][0, 0].re)
        c0 = cos(dt0)
        s0 = sin(dt0)
        c1 = cos(dt1)
        s1 = sin(dt1)
        v0r = V[0][z][0, 0].re
        v0i = V[0][z][0, 0].im
        v1r = V[1][z][0, 0].re
        v1i = V[1][z][0, 0].im
      u0[z][0, 0].re := c0*v0r-s0*v0i
      u0[z][0, 0].im := c0*v0i+s0*v0r
      u1[z][0, 0].re := c1*v1r-s1*v1i
      u1[z][0, 0].im := c1*v1i+s1*v1r
  (@[u0, u1], ld)

proc smearFlowHost*(V: shared.Gauge; m: CircleMap; masks: seq[PortalMask]; sweeps: int): tuple[u: shared.Gauge, lndet: float] =
  doAssert sweeps >= 0, "portal sweep count must be nonnegative"
  result.u = V
  for _ in 0..<sweeps:
    for mask in masks:
      let e = portalSmearHost(result.u, mask, m)
      result.u = e.u
      result.lndet += e.lndet

proc invertPortalSmear*(U: auto; mask: PortalMask; m: CircleMap;
                        tol = 2e-14; maxIter = 80) =
  let
    pa = plaqAngleField(U)
    D = U[0].newOneOf
    sb = newShifters(U[0], -1)
  threads:
    for x in D:
      var d = pa[x][0, 0].re
      for k in 0..<simdLength(d):
        if mask[x][0, 0].re[k] != 0.0:
          let q = pa[x][0, 0].re[k]
          d[k] = q-circleInv(m, q, tol, maxIter)
        else:
          d[k] = 0.0
      D[x][0, 0].re := d
    threadBarrier()
    let
      Dm0 = sb[0] ^* D
      Dm1 = sb[1] ^* D
    threadBarrier()
    for z in U[0]:
      let
        dt0 = 0.25*(D[z][0, 0].re-Dm1[z][0, 0].re)
        dt1 = 0.25*(Dm0[z][0, 0].re-D[z][0, 0].re)
        c0 = cos(dt0)
        s0 = sin(dt0)
        c1 = cos(dt1)
        s1 = sin(dt1)
        u0r = U[0][z][0, 0].re
        u0i = U[0][z][0, 0].im
        u1r = U[1][z][0, 0].re
        u1i = U[1][z][0, 0].im
      U[0][z][0, 0].re := c0*u0r+s0*u0i
      U[0][z][0, 0].im := c0*u0i-s0*u0r
      U[1][z][0, 0].re := c1*u1r+s1*u1i
      U[1][z][0, 0].im := c1*u1i-s1*u1r

proc invertPortalFlow*(U: auto; m: CircleMap; masks: seq[PortalMask]; sweeps: int;
                       tol = 2e-14; maxIter = 80) =
  doAssert sweeps >= 0, "portal sweep count must be nonnegative"
  if sweeps == 0: return
  for _ in countdown(sweeps-1, 0):
    for j in countdown(masks.len-1, 0):
      invertPortalSmear(U, masks[j], m, tol, maxIter)

# Local action correction for the block5 coupling layer.
proc blockNeighborAngles(pa: DLatticeColorMatrixV): array[4, DLatticeColorMatrixV] =
  let
    n0 = pa.newOneOf
    n1 = pa.newOneOf
    n2 = pa.newOneOf
    n3 = pa.newOneOf
    sf = newShifters(pa, 1)
    sb = newShifters(pa, -1)
  threads:
    let
      s0 = sb[1] ^* pa
      s1 = sf[0] ^* pa
      s2 = sf[1] ^* pa
      s3 = sb[0] ^* pa
    threadBarrier()
    for x in n0:
      n0[x][0, 0].re := s0[x][0, 0].re
      n1[x][0, 0].re := s1[x][0, 0].re
      n2[x][0, 0].re := s2[x][0, 0].re
      n3[x][0, 0].re := s3[x][0, 0].re
  [n0, n1, n2, n3]

type BlockForce = tuple[effMinusAux, kc: float; kn: array[4, float]]

proc blockCorrectionOp(V: Ggauge; mask: PortalMask; beta: float; evalBlock: proc(p: float; n: array[4, float]): BlockForce): Gscalar =
  proc forward(v: Gvalue) =
    let
      V = Ggauge(v.inputs[0])
      z = Gscalar(v)
      pa = plaqAngleField(V.gval)
      nb = blockNeighborAngles(pa)
    var acc = 0.0
    threads:
      var s = 0.0
      for x in pa:
        for k in 0..<simdLength(pa[x][0, 0].re):
          if mask[x][0, 0].re[k] != 0.0:
            let n = [nb[0][x][0, 0].re[k], nb[1][x][0, 0].re[k],
                     nb[2][x][0, 0].re[k], nb[3][x][0, 0].re[k]]
            s += evalBlock(pa[x][0, 0].re[k], n).effMinusAux
      s.threadRankSum
      threadSingle: acc = s
    z.sval = acc
  proc gradKernel(V: Ggauge): Ggauge =
    proc kf(v: Gvalue) =
      let
        V = Ggauge(v.inputs[0])
        cV = Ggauge(v)
        pa = plaqAngleField(V.gval)
        nb = blockNeighborAngles(pa)
        kc = V.gval[0].newOneOf
        k0 = V.gval[0].newOneOf
        k1 = V.gval[0].newOneOf
        k2 = V.gval[0].newOneOf
        k3 = V.gval[0].newOneOf
        dK = V.gval[0].newOneOf
      threads:
        for x in pa:
          var
            vc = pa[x][0, 0].re
            v0 = vc
            v1 = vc
            v2 = vc
            v3 = vc
          for k in 0..<simdLength(vc):
            if mask[x][0, 0].re[k] != 0.0:
              let
                p = pa[x][0, 0].re[k]
                n = [nb[0][x][0, 0].re[k], nb[1][x][0, 0].re[k],
                     nb[2][x][0, 0].re[k], nb[3][x][0, 0].re[k]]
                e = evalBlock(p, n)
              vc[k] = e.kc-beta*sin(p)
              v0[k] = e.kn[0]-beta*sin(n[0])
              v1[k] = e.kn[1]-beta*sin(n[1])
              v2[k] = e.kn[2]-beta*sin(n[2])
              v3[k] = e.kn[3]-beta*sin(n[3])
            else:
              vc[k] = 0.0
              v0[k] = 0.0
              v1[k] = 0.0
              v2[k] = 0.0
              v3[k] = 0.0
          kc[x][0, 0].re := vc
          k0[x][0, 0].re := v0
          k1[x][0, 0].re := v1
          k2[x][0, 0].re := v2
          k3[x][0, 0].re := v3
      let
        sf = newShifters(V.gval[0], 1)
        sb = newShifters(V.gval[0], -1)
      threads:
        let
          q0 = sf[1] ^* k0
          q1 = sb[0] ^* k1
          q2 = sb[1] ^* k2
          q3 = sf[0] ^* k3
        threadBarrier()
        for x in dK:
          dK[x][0, 0].re := kc[x][0, 0].re+q0[x][0, 0].re+
            q1[x][0, 0].re+q2[x][0, 0].re+q3[x][0, 0].re
        threadBarrier()
        let
          dm0 = sb[0] ^* dK
          dm1 = sb[1] ^* dK
        threadBarrier()
        for x in cV.gval[0]:
          let
            a0 = dK[x][0, 0].re-dm1[x][0, 0].re
            a1 = dm0[x][0, 0].re-dK[x][0, 0].re
            v0r = V.gval[0][x][0, 0].re
            v0i = V.gval[0][x][0, 0].im
            v1r = V.gval[1][x][0, 0].re
            v1i = V.gval[1][x][0, 0].im
          cV.gval[0][x][0, 0].re := -a0*v0i
          cV.gval[0][x][0, 0].im := a0*v0r
          cV.gval[1][x][0, 0].re := -a1*v1i
          cV.gval[1][x][0, 0].im := a1*v1r
    proc kb(zb, z: Gvalue; i: int; input: Gvalue): Gvalue =
      raiseUnsupportedPath("blockCorrection gradient backward", "block correction gradient is not differentiated")
    let f = Gfunc(forward: kf, backward: kb, name: "blockCorrectionGrad")
    graphNode(V.gaugeNodeLike, @[Gvalue(V)], f, "blockCorrectionGrad")
  proc backward(zb, z: Gvalue; i: int; input: Gvalue): Gvalue =
    scaledUpstreamOr(zb, Gscalar, gradKernel(Ggauge(z.inputs[0])))
  let f = Gfunc(forward: forward, backward: backward, name: "blockCorrection")
  graphNode(scalarNodeLike(V), @[Gvalue(V)], f, "blockCorrection")

# Link2: transform one link through its two adjacent plaquette angles.
proc pairAngles(pa: DLatticeColorMatrixV; dir: int): tuple[plus, minus: DLatticeColorMatrixV] =
  doAssert dir in 0..1
  let
    pp = pa.newOneOf
    pm = pa.newOneOf
    sb = newShifters(pa, -1)
  threads:
    let
      p0 = sb[0] ^* pa
      p1 = sb[1] ^* pa
    threadBarrier()
    if dir == 0:
      for x in pp:
        pp[x][0, 0].re := pa[x][0, 0].re
        pm[x][0, 0].re := p1[x][0, 0].re
    else:
      for x in pp:
        pp[x][0, 0].re := p0[x][0, 0].re
        pm[x][0, 0].re := pa[x][0, 0].re
  (pp, pm)

proc pairLogDetJ(V: Ggauge; mask: PortalMask; p: PairMap; dir: int): Gscalar =
  proc forward(v: Gvalue) =
    let
      V = Ggauge(v.inputs[0])
      z = Gscalar(v)
      pa = plaqAngleField(V.gval)
      ps = pairAngles(pa, dir)
    var acc = 0.0
    threads:
      var s = 0.0
      for x in pa:
        for k in 0..<simdLength(pa[x][0, 0].re):
          if mask[x][0, 0].re[k] != 0.0:
            let
              pp = ps.plus[x][0, 0].re[k]
              pm = ps.minus[x][0, 0].re[k]
            s += evalPair(p, pp, pm).logdet
      s.threadRankSum
      threadSingle: acc = s
    z.sval = acc
  proc gradKernel(V: Ggauge): Ggauge =
    proc kf(v: Gvalue) =
      let
        V = Ggauge(v.inputs[0])
        cV = Ggauge(v)
        pa = plaqAngleField(V.gval)
        ps = pairAngles(pa, dir)
        kp = V.gval[0].newOneOf
        km = V.gval[0].newOneOf
        dK = V.gval[0].newOneOf
      threads:
        for x in pa:
          var
            vp = pa[x][0, 0].re
            vm = vp
          for k in 0..<simdLength(vp):
            if mask[x][0, 0].re[k] != 0.0:
              let
                pp = ps.plus[x][0, 0].re[k]
                pm = ps.minus[x][0, 0].re[k]
                e = evalPair(p, pp, pm)
              vp[k] = e.logdetPlus
              vm[k] = e.logdetMinus
            else:
              vp[k] = 0.0
              vm[k] = 0.0
          kp[x][0, 0].re := vp
          km[x][0, 0].re := vm
      let
        sf = newShifters(V.gval[0], 1)
        sb = newShifters(V.gval[0], -1)
      threads:
        if dir == 0:
          let q = sf[1] ^* km
          threadBarrier()
          for x in dK: dK[x][0, 0].re := kp[x][0, 0].re+q[x][0, 0].re
        else:
          let q = sf[0] ^* kp
          threadBarrier()
          for x in dK: dK[x][0, 0].re := q[x][0, 0].re+km[x][0, 0].re
        threadBarrier()
        let
          dm0 = sb[0] ^* dK
          dm1 = sb[1] ^* dK
        threadBarrier()
        for x in cV.gval[0]:
          let
            a0 = dK[x][0, 0].re-dm1[x][0, 0].re
            a1 = dm0[x][0, 0].re-dK[x][0, 0].re
            v0r = V.gval[0][x][0, 0].re
            v0i = V.gval[0][x][0, 0].im
            v1r = V.gval[1][x][0, 0].re
            v1i = V.gval[1][x][0, 0].im
          cV.gval[0][x][0, 0].re := -a0*v0i
          cV.gval[0][x][0, 0].im := a0*v0r
          cV.gval[1][x][0, 0].re := -a1*v1i
          cV.gval[1][x][0, 0].im := a1*v1r
    proc kb(zb, z: Gvalue; i: int; input: Gvalue): Gvalue =
      raiseUnsupportedPath("pairLogDetJ gradient backward", "pair log-Jacobian gradient is not differentiated")
    let f = Gfunc(forward: kf, backward: kb, name: "pairLogDetJgrad")
    graphNode(V.gaugeNodeLike, @[Gvalue(V)], f, "pairLogDetJgrad")
  proc backward(zb, z: Gvalue; i: int; input: Gvalue): Gvalue =
    scaledUpstreamOr(zb, Gscalar, gradKernel(Ggauge(z.inputs[0])))
  let f = Gfunc(forward: forward, backward: backward, name: "pairLogDetJ")
  graphNode(scalarNodeLike(V), @[Gvalue(V)], f, "pairLogDetJ")

proc pairSmearGrad(V, cU: Ggauge; mask: PortalMask; p: PairMap; dir: int): Ggauge =
  proc kf(v: Gvalue) =
    let
      V = Ggauge(v.inputs[0])
      cU = Ggauge(v.inputs[1])
      cV = Ggauge(v)
      pa = plaqAngleField(V.gval)
      ps = pairAngles(pa, dir)
      D = V.gval[0].newOneOf
      DP = V.gval[0].newOneOf
      DM = V.gval[0].newOneOf
    threads:
      for x in pa:
        var
          d = pa[x][0, 0].re
          dp = d
          dm = d
        for k in 0..<simdLength(d):
          if mask[x][0, 0].re[k] != 0.0:
            let e = evalPair(p, ps.plus[x][0, 0].re[k], ps.minus[x][0, 0].re[k])
            d[k] = e.delta
            dp[k] = e.deltaPlus
            dm[k] = e.deltaMinus
          else:
            d[k] = 0.0
            dp[k] = 0.0
            dm[k] = 0.0
        D[x][0, 0].re := d
        DP[x][0, 0].re := dp
        DM[x][0, 0].re := dm
    let
      U0 = V.gval[0].newOneOf
      U1 = V.gval[0].newOneOf
      R = V.gval[0].newOneOf
      SP = V.gval[0].newOneOf
      SM = V.gval[0].newOneOf
      dP = V.gval[0].newOneOf
      sf = newShifters(V.gval[0], 1)
      sb = newShifters(V.gval[0], -1)
    threads:
      for x in U0:
        let
          d = D[x][0, 0].re
          c = cos(d)
          s = sin(d)
          v0r = V.gval[0][x][0, 0].re
          v0i = V.gval[0][x][0, 0].im
          v1r = V.gval[1][x][0, 0].re
          v1i = V.gval[1][x][0, 0].im
        if dir == 0:
          let
            u0r = c*v0r-s*v0i
            u0i = c*v0i+s*v0r
          U0[x][0, 0].re := u0r
          U0[x][0, 0].im := u0i
          U1[x][0, 0].re := v1r
          U1[x][0, 0].im := v1i
          let
            cr = cU.gval[0][x][0, 0].re
            ci = cU.gval[0][x][0, 0].im
          R[x][0, 0].re := ci*u0r-cr*u0i
        else:
          let
            u1r = c*v1r-s*v1i
            u1i = c*v1i+s*v1r
          U1[x][0, 0].re := u1r
          U1[x][0, 0].im := u1i
          U0[x][0, 0].re := v0r
          U0[x][0, 0].im := v0i
          let
            cr = cU.gval[1][x][0, 0].re
            ci = cU.gval[1][x][0, 0].im
          R[x][0, 0].re := ci*u1r-cr*u1i
      threadBarrier()
      for x in SP:
        SP[x][0, 0].re := R[x][0, 0].re*DP[x][0, 0].re
        SM[x][0, 0].re := R[x][0, 0].re*DM[x][0, 0].re
      threadBarrier()
      if dir == 0:
        let q = sf[1] ^* SM
        threadBarrier()
        for x in dP: dP[x][0, 0].re := SP[x][0, 0].re+q[x][0, 0].re
      else:
        let q = sf[0] ^* SP
        threadBarrier()
        for x in dP: dP[x][0, 0].re := q[x][0, 0].re+SM[x][0, 0].re
      threadBarrier()
      let
        dm0 = sb[0] ^* dP
        dm1 = sb[1] ^* dP
      threadBarrier()
      for x in cV.gval[0]:
        let
          v0r = V.gval[0][x][0, 0].re
          v0i = V.gval[0][x][0, 0].im
          v1r = V.gval[1][x][0, 0].re
          v1i = V.gval[1][x][0, 0].im
          u0r = U0[x][0, 0].re
          u0i = U0[x][0, 0].im
          u1r = U1[x][0, 0].re
          u1i = U1[x][0, 0].im
          c0r = cU.gval[0][x][0, 0].re
          c0i = cU.gval[0][x][0, 0].im
          c1r = cU.gval[1][x][0, 0].re
          c1i = cU.gval[1][x][0, 0].im
          p0r = u0r*v0r+u0i*v0i
          p0i = u0r*v0i-u0i*v0r
          p1r = u1r*v1r+u1i*v1i
          p1i = u1r*v1i-u1i*v1r
          a0 = dP[x][0, 0].re-dm1[x][0, 0].re
          a1 = dm0[x][0, 0].re-dP[x][0, 0].re
        cV.gval[0][x][0, 0].re := p0r*c0r-p0i*c0i-a0*v0i
        cV.gval[0][x][0, 0].im := p0r*c0i+p0i*c0r+a0*v0r
        cV.gval[1][x][0, 0].re := p1r*c1r-p1i*c1i-a1*v1i
        cV.gval[1][x][0, 0].im := p1r*c1i+p1i*c1r+a1*v1r
  proc kb(zb, z: Gvalue; i: int; input: Gvalue): Gvalue =
    raiseUnsupportedPath("pairSmearGrad backward", "pair field-map gradient is not differentiated")
  let f = Gfunc(forward: kf, backward: kb, name: "pairSmearGrad")
  graphNode(V.gaugeNodeLike, @[Gvalue(V), Gvalue(cU)], f, "pairSmearGrad")

proc pairSmear*(V: Ggauge; mask: PortalMask; p: PairMap; dir: int): Ggauge =
  proc forward(v: Gvalue) =
    let
      V = Ggauge(v.inputs[0])
      U = Ggauge(v)
      pa = plaqAngleField(V.gval)
      ps = pairAngles(pa, dir)
      D = V.gval[0].newOneOf
    threads:
      for x in pa:
        var d = pa[x][0, 0].re
        for k in 0..<simdLength(d):
          if mask[x][0, 0].re[k] != 0.0:
            d[k] = evalPair(p, ps.plus[x][0, 0].re[k], ps.minus[x][0, 0].re[k]).delta
          else:
            d[k] = 0.0
        D[x][0, 0].re := d
      threadBarrier()
      for x in U.gval[0]:
        let
          d = D[x][0, 0].re
          c = cos(d)
          s = sin(d)
          v0r = V.gval[0][x][0, 0].re
          v0i = V.gval[0][x][0, 0].im
          v1r = V.gval[1][x][0, 0].re
          v1i = V.gval[1][x][0, 0].im
        if dir == 0:
          U.gval[0][x][0, 0].re := c*v0r-s*v0i
          U.gval[0][x][0, 0].im := c*v0i+s*v0r
          U.gval[1][x][0, 0].re := v1r
          U.gval[1][x][0, 0].im := v1i
        else:
          U.gval[1][x][0, 0].re := c*v1r-s*v1i
          U.gval[1][x][0, 0].im := c*v1i+s*v1r
          U.gval[0][x][0, 0].re := v0r
          U.gval[0][x][0, 0].im := v0i
  proc backward(zb, z: Gvalue; i: int; input: Gvalue): Gvalue =
    let cU = requireUpstream(zb, "pairSmear backward", Ggauge)
    Gvalue(pairSmearGrad(Ggauge(z.inputs[0]), cU, mask, p, dir))
  proc ldj(z: Gvalue): tuple[ld, via: Gvalue] =
    (Gvalue(pairLogDetJ(Ggauge(z.inputs[0]), mask, p, dir)), z.inputs[0])
  let f = Gfunc(forward: forward, backward: backward, logdet: ldj, name: "pairSmear")
  graphNode(V.gaugeNodeLike, @[Gvalue(V)], f, "pairSmear")

proc pairSmearHost*(V: shared.Gauge; mask: PortalMask; p: PairMap; dir: int): tuple[u: shared.Gauge, lndet: float] =
  let
    pa = plaqAngleField(V)
    ps = pairAngles(pa, dir)
    D = V[0].newOneOf
    u0 = V[0].newOneOf
    u1 = V[1].newOneOf
  var ld = 0.0
  threads:
    var sld = 0.0
    for x in pa:
      var d = pa[x][0, 0].re
      for k in 0..<simdLength(d):
        if mask[x][0, 0].re[k] != 0.0:
          let e = evalPair(p, ps.plus[x][0, 0].re[k], ps.minus[x][0, 0].re[k])
          d[k] = e.delta
          sld += e.logdet
        else:
          d[k] = 0.0
      D[x][0, 0].re := d
    sld.threadRankSum
    threadSingle: ld = sld
    threadBarrier()
    for x in u0:
      let
        d = D[x][0, 0].re
        c = cos(d)
        s = sin(d)
        v0r = V[0][x][0, 0].re
        v0i = V[0][x][0, 0].im
        v1r = V[1][x][0, 0].re
        v1i = V[1][x][0, 0].im
      if dir == 0:
        u0[x][0, 0].re := c*v0r-s*v0i
        u0[x][0, 0].im := c*v0i+s*v0r
        u1[x][0, 0].re := v1r
        u1[x][0, 0].im := v1i
      else:
        u1[x][0, 0].re := c*v1r-s*v1i
        u1[x][0, 0].im := c*v1i+s*v1r
        u0[x][0, 0].re := v0r
        u0[x][0, 0].im := v0i
  (@[u0, u1], ld)

proc invertPairSmear*(U: auto; mask: PortalMask; p: PairMap; dir: int) =
  let
    pa = plaqAngleField(U)
    ps = pairAngles(pa, dir)
    D = U[0].newOneOf
  threads:
    for x in pa:
      var d = pa[x][0, 0].re
      for k in 0..<simdLength(d):
        if mask[x][0, 0].re[k] != 0.0:
          let
            pp = ps.plus[x][0, 0].re[k]
            pm = ps.minus[x][0, 0].re[k]
            a = invertPair(p, pp, pm)
          d[k] = a[0]-pp
        else:
          d[k] = 0.0
      D[x][0, 0].re := d
    threadBarrier()
    for x in U[0]:
      let
        d = D[x][0, 0].re
        c = cos(d)
        s = sin(d)
      if dir == 0:
        let
          ur = U[0][x][0, 0].re
          ui = U[0][x][0, 0].im
        U[0][x][0, 0].re := c*ur-s*ui
        U[0][x][0, 0].im := c*ui+s*ur
      else:
        let
          ur = U[1][x][0, 0].re
          ui = U[1][x][0, 0].im
        U[1][x][0, 0].re := c*ur-s*ui
        U[1][x][0, 0].im := c*ui+s*ur

type PairLayer* = object
  mask*: PortalMask
  dir*: int

proc pairFlowHost*(V: shared.Gauge; p: PairMap; layers: seq[PairLayer]; rounds: int): tuple[u: shared.Gauge, lndet: float] =
  if rounds < 1: mapFail("pair-flow rounds must be positive")
  result.u = V
  for _ in 0..<rounds:
    for layer in layers:
      let e = pairSmearHost(result.u, layer.mask, p, layer.dir)
      result.u = e.u
      result.lndet += e.lndet

proc invertPairFlow*(U: auto; p: PairMap; layers: seq[PairLayer]; rounds: int) =
  if rounds < 1: mapFail("pair-flow rounds must be positive")
  for _ in countdown(rounds-1, 0):
    for j in countdown(layers.len-1, 0):
      invertPairSmear(U, layers[j].mask, p, layers[j].dir)

proc pairFlow*(V: Ggauge; p: PairMap; layers: seq[PairLayer]; rounds: int): Ggauge =
  if rounds < 1: mapFail("pair-flow rounds must be positive")
  result = V
  for _ in 0..<rounds:
    for layer in layers:
      result = pairSmear(result, layer.mask, p, layer.dir)

# Block5 coupling: four coordinates at fixed five-plaquette flux.
type
  Vec4* = array[4, float]
  Vec5* = array[5, float]
  Mat4* = array[4, Vec4]
  Mat5* = array[5, Vec5]

  BlockMap* = object
    beta*: float
    order*: array[4, int]
    stages*: array[4, seq[ContextMap]]
    invTol*: float
    invIter*: int

  BlockEval* = object
    auxiliary*, physical*: Vec5
    physicalZ*: Vec4
    jacobian*: Mat4
    derivatives*: Mat5
    det*: float
    detGradient*: Vec4
    detFluxGradient*: float
    logdet*: float
    logdetGradient*: Vec4
    logdetFluxGradient*: float
    action*: float
    gradient*: Vec4
    fluxForce*: float
    deltas*: Vec4

proc det4*(a: Mat4): float =
  var x = a
  result = 1.0
  for j in 0..3:
    var p = j
    for i in j+1..3:
      if abs(x[i][j]) > abs(x[p][j]): p = i
    if x[p][j] == 0.0: return 0.0
    if p != j:
      swap(x[p], x[j])
      result = -result
    let d = x[j][j]
    result *= d
    for i in j+1..3:
      let q = x[i][j]/d
      for k in j+1..3: x[i][k] -= q*x[j][k]

proc blockContext(y: Vec4; flux: float; active: int): tuple[c: array[maxMapContext, float]; index: array[maxMapContext, int]] =
  var q = 0
  for j in 0..3:
    if j != active:
      result.c[q] = y[j]
      result.index[q] = j
      inc q
  result.c[3] = flux
  result.index[3] = -1

proc evalBlockMap*(m: BlockMap; z: Vec4; flux: float): BlockEval =
  # p4=flux-Σ_i z_i; each scalar coupling updates one z_i at fixed flux.
  var y = z
  var jac: Mat4
  var jf: Vec4
  result.det = 1.0
  for i in 0..3: jac[i][i] = 1.0
  for i in 0..3: result.auxiliary[i] = z[i]
  result.auxiliary[4] = flux-z[0]-z[1]-z[2]-z[3]

  let depth = m.stages[0].len
  for round in 0..<depth:
    for active in m.order:
      let
        bc = blockContext(y, flux, active)
        e = evalContext(m.stages[active][round], y[active], bc.c)
        lx = e.dxx/e.dx
      var lc: array[maxMapContext, float]
      for k in 0..<maxMapContext: lc[k] = e.dxc[k]/e.dx
      for col in 0..3:
        result.logdetGradient[col] += lx*jac[active][col]
        for k in 0..<maxMapContext:
          if bc.index[k] >= 0:
            result.logdetGradient[col] += lc[k]*jac[bc.index[k]][col]
      result.logdetFluxGradient += lx*jf[active]
      for k in 0..<maxMapContext:
        result.logdetFluxGradient += lc[k]*(if bc.index[k] < 0: 1.0 else: jf[bc.index[k]])

      var row: Vec4
      for col in 0..3:
        row[col] = e.dx*jac[active][col]
        for k in 0..<maxMapContext:
          if bc.index[k] >= 0:
            row[col] += e.dc[k]*jac[bc.index[k]][col]
      var rowf = e.dx*jf[active]
      for k in 0..<maxMapContext:
        rowf += e.dc[k]*(if bc.index[k] < 0: 1.0 else: jf[bc.index[k]])
      jac[active] = row
      jf[active] = rowf
      y[active] = e.y
      result.det *= e.dx
      result.logdet += ln(e.dx)

  result.physicalZ = y
  for i in 0..3:
    result.physical[i] = y[i]
    result.jacobian[i] = jac[i]
    for j in 0..3: result.derivatives[i][j] = jac[i][j]
    result.derivatives[i][4] = jf[i]
  result.physical[4] = flux-y[0]-y[1]-y[2]-y[3]
  for j in 0..3:
    for i in 0..3: result.derivatives[4][j] -= jac[i][j]
    result.derivatives[4][4] -= jf[j]
  result.derivatives[4][4] += 1.0

  for j in 0..3:
    result.detGradient[j] = result.det*result.logdetGradient[j]
  result.detFluxGradient = result.det*result.logdetFluxGradient

  for i in 0..4: result.action -= m.beta*cos(result.physical[i])
  result.action -= result.logdet
  for col in 0..3:
    for i in 0..4:
      result.gradient[col] += result.derivatives[i][col]*m.beta*sin(result.physical[i])
    result.gradient[col] -= result.logdetGradient[col]
  for i in 0..4:
    result.fluxForce += result.derivatives[i][4]*m.beta*sin(result.physical[i])
  result.fluxForce -= result.logdetFluxGradient
  for i in 0..3:
    result.deltas[i] = result.auxiliary[i+1]-result.physical[i+1]

proc invertBlockMap*(m: BlockMap; physicalZ: Vec4; flux: float): Vec4 =
  result = physicalZ
  let depth = m.stages[0].len
  for round in countdown(depth-1, 0):
    for oi in countdown(3, 0):
      let
        active = m.order[oi]
        bc = blockContext(result, flux, active)
      result[active] = invertContext(m.stages[active][round], result[active],
        bc.c, m.invTol, m.invIter)

proc blockForce(m: BlockMap; p: float; n: array[4, float]): BlockForce =
  let
    z: Vec4 = [p, n[0], n[1], n[2]]
    flux = p+n[0]+n[1]+n[2]+n[3]
    e = evalBlockMap(m, z, flux)
  var aux = 0.0
  for x in [p, n[0], n[1], n[2], n[3]]: aux -= m.beta*cos(x)
  result.effMinusAux = e.action-aux
  result.kc = e.gradient[0]+e.fluxForce
  for j in 0..2: result.kn[j] = e.gradient[j+1]+e.fluxForce
  result.kn[3] = e.fluxForce

proc blockCorrection*(V: Ggauge; mask: PortalMask; m: BlockMap): Gscalar =
  blockCorrectionOp(V, mask, m.beta,
    proc(p: float; n: array[4, float]): BlockForce = blockForce(m, p, n))

proc blockSmearHost*(V: shared.Gauge; mask: PortalMask; m: BlockMap): tuple[u: shared.Gauge, lndet: float] =
  let
    pa = plaqAngleField(V)
    nb = blockNeighborAngles(pa)
    a0 = V[0].newOneOf
    a1 = V[0].newOneOf
    a2 = V[0].newOneOf
    a3 = V[0].newOneOf
  var ld = 0.0
  threads:
    var sld = 0.0
    for x in pa:
      var
        d0 = pa[x][0, 0].re
        d1 = d0
        d2 = d0
        d3 = d0
      for k in 0..<simdLength(d0):
        if mask[x][0, 0].re[k] != 0.0:
          let
            p = pa[x][0, 0].re[k]
            n = [nb[0][x][0, 0].re[k], nb[1][x][0, 0].re[k],
                 nb[2][x][0, 0].re[k], nb[3][x][0, 0].re[k]]
            z: Vec4 = [p, n[0], n[1], n[2]]
            e = evalBlockMap(m, z, p+n[0]+n[1]+n[2]+n[3])
          d0[k] = e.deltas[0]
          d1[k] = e.deltas[1]
          d2[k] = e.deltas[2]
          d3[k] = e.deltas[3]
          sld += e.logdet
        else:
          d0[k] = 0.0
          d1[k] = 0.0
          d2[k] = 0.0
          d3[k] = 0.0
      a0[x][0, 0].re := d0
      a1[x][0, 0].re := d1
      a2[x][0, 0].re := d2
      a3[x][0, 0].re := d3
    sld.threadRankSum
    threadSingle: ld = sld
  let
    u0 = V[0].newOneOf
    u1 = V[1].newOneOf
    sb = newShifters(V[0], -1)
  threads:
    let
      a2m1 = sb[1] ^* a2
      a1m0 = sb[0] ^* a1
    threadBarrier()
    for x in u0:
      let
        d0 = a0[x][0, 0].re-a2m1[x][0, 0].re
        d1 = -a3[x][0, 0].re+a1m0[x][0, 0].re
        c0 = cos(d0)
        s0 = sin(d0)
        c1 = cos(d1)
        s1 = sin(d1)
        v0r = V[0][x][0, 0].re
        v0i = V[0][x][0, 0].im
        v1r = V[1][x][0, 0].re
        v1i = V[1][x][0, 0].im
      u0[x][0, 0].re := c0*v0r-s0*v0i
      u0[x][0, 0].im := c0*v0i+s0*v0r
      u1[x][0, 0].re := c1*v1r-s1*v1i
      u1[x][0, 0].im := c1*v1i+s1*v1r
  (@[u0, u1], ld)

proc invertBlockSmear*(U: auto; mask: PortalMask; m: BlockMap) =
  let
    pa = plaqAngleField(U)
    nb = blockNeighborAngles(pa)
    a0 = U[0].newOneOf
    a1 = U[0].newOneOf
    a2 = U[0].newOneOf
    a3 = U[0].newOneOf
  threads:
    for x in pa:
      var
        d0 = pa[x][0, 0].re
        d1 = d0
        d2 = d0
        d3 = d0
      for k in 0..<simdLength(d0):
        if mask[x][0, 0].re[k] != 0.0:
          let
            p = pa[x][0, 0].re[k]
            n = [nb[0][x][0, 0].re[k], nb[1][x][0, 0].re[k],
                 nb[2][x][0, 0].re[k], nb[3][x][0, 0].re[k]]
            flux = p+n[0]+n[1]+n[2]+n[3]
            original = invertBlockMap(m, [p, n[0], n[1], n[2]], flux)
            original4 = flux-original[0]-original[1]-original[2]-original[3]
          d0[k] = n[0]-original[1]
          d1[k] = n[1]-original[2]
          d2[k] = n[2]-original[3]
          d3[k] = n[3]-original4
        else:
          d0[k] = 0.0
          d1[k] = 0.0
          d2[k] = 0.0
          d3[k] = 0.0
      a0[x][0, 0].re := d0
      a1[x][0, 0].re := d1
      a2[x][0, 0].re := d2
      a3[x][0, 0].re := d3
  let sb = newShifters(U[0], -1)
  threads:
    let
      a2m1 = sb[1] ^* a2
      a1m0 = sb[0] ^* a1
    threadBarrier()
    for x in U[0]:
      let
        d0 = a0[x][0, 0].re-a2m1[x][0, 0].re
        d1 = -a3[x][0, 0].re+a1m0[x][0, 0].re
        c0 = cos(d0)
        s0 = sin(d0)
        c1 = cos(d1)
        s1 = sin(d1)
        u0r = U[0][x][0, 0].re
        u0i = U[0][x][0, 0].im
        u1r = U[1][x][0, 0].re
        u1i = U[1][x][0, 0].im
      U[0][x][0, 0].re := c0*u0r-s0*u0i
      U[0][x][0, 0].im := c0*u0i+s0*u0r
      U[1][x][0, 0].re := c1*u1r-s1*u1i
      U[1][x][0, 0].im := c1*u1i+s1*u1r

proc blockAction*(gc: Gactcoeff; mask: PortalMask; m: BlockMap): proc(V: Ggauge): Gscalar =
  result = proc(V: Ggauge): Gscalar = gaugeAction(gc, V)+blockCorrection(V, mask, m)

# Runtime selectors, tuning parameters, and geometry dispatch.
# plaq4/scalar: identity, sine, sqfourier, cspline, bspline, fejer.
# link2/scalar: all above except cspline, plus stout.
# block5/chain: all link2 bases except stout.
# block5/coupling: identity, sine, sqfourier, bspline.
