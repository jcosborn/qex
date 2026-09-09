## Triangle gauge action and force on the 16-cell honeycomb.
##
##   S(U) = (beta/2) sum_x sum_{i=1}^{32} (1 - Re Tr P_i(x)/N)     (apexTris: each triangle once)
##
## Each link sits in 8 triangles; with P = U_l V_k^dag (QEX staple convention)
##   D_l = (beta/2N) sum_{k=1}^{8} V_k                  actionDeriv (not projected)
##   f_l = projectTAH(U_l D_l^dag)                      force
##   d/ds S(exp(sP) U)|_0 = sum_l redot(P_l, f_l)       (as QEX gaugeForce)
## so p -= dt f, U := exp(dt p) U conserves S + (1/2) sum_l redot(p_l, p_l).
##
## Staples per apexTris entry (mu, delta): deltaP = delta or 2^mu,
## db = delta xor 15, dbp = deltaP xor 15.
##   apex-B  P = uD[delta](z) uA[mu](z+delta) uD[deltaP](z)^dag
##     V(uD[delta])  = uD[deltaP] uA[mu](z+delta)^dag
##     V(uD[deltaP]) = uD[delta] uA[mu](z+delta)
##     V(uA[mu])(y)  = [uD[delta]^dag uD[deltaP]](y-delta)
##   apex-A  P = uD[db](z)^dag uB[mu](z) uD[dbp](z+e_mu)      (re-based at z = apex-db)
##     V(uB[mu])  = uD[db] uD[dbp](z+e_mu)^dag
##     V(uD[db])  = uB[mu] uD[dbp](z+e_mu)
##     V(uD[dbp])(y) = [uB[mu]^dag uD[db]](y-e_mu)
##
## ActionWork holds every shifter and scratch field: build it once outside
## `threads:`, then action/actionDeriv/force never allocate.

import base, layout, field, maths, rng
import physics/qcdTypes
import gauge
import hcgeom, hcgauge

export hcgauge

type
  ActionWork*[F, SH, SS] = ref object
    shA*: array[nDim, SH]           ## forward 16-trees: uA[mu](x+delta)
    sF*: array[nDim, SS]            ## +e_mu shifters: uD[dbp](x+e_mu)
    sB*: array[3, array[nDim, SS]]  ## -e_mu shifter chains (3 levels)
    t*: F                           ## 2-link product scratch

proc newActionWork*[F](g: HcGauge[F]): auto =
  ## allocates; outside `threads:`.  Serves any gauge of the same shape.
  type SH = type(newHcShift16(g.uA[0], 1))
  type SS = type(newShifter(g.uD[0], 0, 1))
  var w = ActionWork[F, SH, SS]()
  for mu in 0..<nDim:
    w.shA[mu] = newHcShift16(g.uA[mu], 1)
    w.sF[mu] = newShifter(g.uD[0], mu, 1)
    for l in 0..<3:
      w.sB[l][mu] = newShifter(g.uD[0], mu, -1)
  w.t = g.uA[0].newOneOf
  w

proc rebind[F, SH, SS](w: ActionWork[F, SH, SS], g: HcGauge[F]) =
  ## ref assignments; outside `threads:`
  for mu in 0..<nDim:
    w.shA[mu].setSrc g.uA[mu]

proc action*[F, SH, SS](w: ActionWork[F, SH, SS], beta: float,
                        g: HcGauge[F]): float =
  ## S(U) through redot(staple, link), one triangle per apex label
  const nc = g.uA[0][0].nrows
  rebind(w, g)
  var tsum = 0.0
  threads:
    for mu in 0..<nDim:
      w.shA[mu].run
    var a = 0.0
    for t in apexTris:
      let
        mu = t.mu
        delta = t.delta
        deltaP = t.deltaP
        db = delta xor 15
        dbp = deltaP xor 15
      w.t := g.uD[deltaP] * w.shA[mu].f[delta].adj
      a += redot(w.t, g.uD[delta])
      let s2 = w.sF[mu] ^* g.uD[dbp]
      w.t := g.uD[db] * s2.adj
      a += redot(w.t, g.uB[mu])
    threadMaster: tsum = a
  let nTri = float(nTriPerSite*2*g.lo.physVol)
  0.5*beta*(nTri - tsum/float(nc))

proc actionDeriv*[F, SH, SS](w: ActionWork[F, SH, SS], beta: float,
                             g: HcGauge[F], f: var HcGauge[F]) =
  ## f_l := (beta/2N) sum_{k=1}^{8} V_k(l); f overwritten, distinct from g
  const nc = g.uA[0][0].nrows
  let ff = f
  rebind(w, g)
  let cf = beta/(2.0*float(nc))
  threads:
    for u in ff.links: u := 0
    for mu in 0..<nDim:
      w.shA[mu].run
    for t in apexTris:
      let
        mu = t.mu
        delta = t.delta
        deltaP = t.deltaP
        db = delta xor 15
        dbp = deltaP xor 15
      ff.uD[delta]  += g.uD[deltaP] * w.shA[mu].f[delta].adj
      ff.uD[deltaP] += g.uD[delta] * w.shA[mu].f[delta]
      w.t := g.uD[delta].adj * g.uD[deltaP]
      block:
        var cur = w.t
        var lev = 0
        for b in 0..<nDim:
          if ((delta shr b) and 1) != 0:
            cur = w.sB[lev][b] ^* cur
            inc lev
        ff.uA[mu] += cur
      let s2 = w.sF[mu] ^* g.uD[dbp]
      ff.uB[mu] += g.uD[db] * s2.adj
      ff.uD[db] += g.uB[mu] * s2
      w.t := g.uB[mu].adj * g.uD[db]
      let rs = w.sB[0][mu] ^* w.t
      ff.uD[dbp] += rs
    if cf != 1.0:
      for u in ff.links: u := cf*u

proc force*[F, SH, SS](w: ActionWork[F, SH, SS], beta: float,
                       g: HcGauge[F], f: var HcGauge[F]) =
  ## f_l = projectTAH(U_l D_l^dag); outside `threads:`
  actionDeriv(w, beta, g, f)
  contractProjectTAH(g.links, f.links)

proc triSum*[F, SH, SS](w: ActionWork[F, SH, SS], beta: float,
                        g: HcGauge[F]): float =
  ## triangleSum through the action: 1 - S/((beta/2) 32 N_sites)
  1.0 - action(w, beta, g)/(0.5*beta*float(nTriPerSite*2*g.lo.physVol))
