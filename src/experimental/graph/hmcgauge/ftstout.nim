## Stout field-transformation HMC:
##
##   U = f(V),  S_eff(V) = S(U) - log det f'(V)
##   ds = gaugeActionDeriv(actWilson(1),W)
##   F = projectTAH(W ds†),  W' = exp(alpha F) W.
##
## Each sweep applies the substep to every (parity,dir).

import std/tables
import ../[core, scalar, gauge]
import integrator
import qex
from ../gauge/action/domain import evalGaugeForceValue

type
  StoutAction* = object
    action*: GaugeAction
    flow*: proc(V: Ggauge): Ggauge {.closure.}

proc smearFlow*(V: Ggauge, c1: Gactcoeff, alpha: Gscalar, sweeps: int): Ggauge =
  ## U = f(V); ln det f'(V) is logDetJ(result, V). c1 and alpha must be
  ## graph-independent of V (the stout logdet hooks reject them otherwise).
  if sweeps < 0:
    raiseValueError("stout sweeps must be >= 0, got " & $sweeps)
  if sweeps > 0:
    for n in V.gval[0].l.physGeom:
      if (n and 1) != 0:
        raiseValueError("stout flow requires even lattice extents")
  let nd = V.gval.len
  result = V
  for s in 0..<sweeps:
    for parity in 0..1:
      for dir in 0..<nd:
        result = stoutUpdateLogDetJ(result, c1, alpha, parity, dir).Wnew

proc smearedField*(V: Ggauge, rho: float, sweeps: int): Ggauge =
  ## U = f(V); ln det f'(V) is logDetJ(result, V).
  smearFlow(V, actWilson(scalar.toGvalue(V.runtime, 1.0)), scalar.toGvalue(V.runtime, rho), sweeps)

proc stoutAction*(gc: Gactcoeff, rho: float, sweeps: int): StoutAction =
  ## S_eff(V) = S(f(V)) - log det f'(V).
  if sweeps < 0:
    raiseValueError("stout sweeps must be >= 0, got " & $sweeps)
  let
    c1 = actWilson(scalar.toGvalue(gc.runtime, 1.0))
    alpha = scalar.toGvalue(gc.runtime, rho)
  var flows = initTable[NodeKey, Ggauge]()
  proc getFlow(V: Ggauge): Ggauge =
    let key = V.nodeKey
    if key notin flows:
      flows[key] = smearFlow(V, c1, alpha, sweeps)
    flows[key]
  result.flow = getFlow
  result.action = proc(V: Ggauge): Gscalar =
    let u = getFlow(V)
    gaugeAction(gc, u) - logDetJ(u, V)

proc invertStoutFlow*(U: auto, rho: float, sweeps: int; rdf2req = 1e-24, maxiter = 1000, verbose = false): tuple[iter: int, rdf2: float] =
  ## Invert U = f(V) in place, reversing the subset order:
  ##   V_S^(k+1) = exp(-rho*projectTAH(V_S^k ds_S†))*U_S.
  ## ds_S depends only on the frozen complement, so compute it once per substep.
  ## Return total iterations and max sum|delta F|^2/sum|F|^2.
  if sweeps < 0:
    raiseValueError("stout sweeps must be >= 0, got " & $sweeps)
  if rdf2req <= 0.0:
    raiseValueError("stout inverse rdf2req must be > 0, got " & $rdf2req)
  if maxiter < 1:
    raiseValueError("stout inverse maxiter must be >= 1, got " & $maxiter)
  if sweeps > 0:
    for n in U[0].l.physGeom:
      if (n and 1) != 0:
        raiseValueError("stout inverse requires even lattice extents")
  let
    c1 = GaugeActionCoeffs(plaq: 1.0)
    nd = U.len
    lo = U[0].l
    ds = lo.newGauge       # staple-derivative (same as the forward op)
    Upost = lo.newGauge    # the post-substep links: fixed RHS of the fixed point
    fprev = lo.newGauge    # previous-iteration projTAH direction (for the residual)
  result = (iter: 0, rdf2: 0.0)
  for sweep in countdown(sweeps-1, 0):
    for parity in countdown(1, 0):
      for dir in countdown(nd-1, 0):
        let sub = lo.getSubset(if parity == 0: "even" else: "odd")
        threads:
          Upost[dir] := U[dir]
          fprev[dir] := 0
        evalGaugeForceValue(c1, U, ds)   # ds_S frozen ⇒ constant over the solve
        var
          iter = 0
          rdf2 = 0.0
          df2o = -1.0
        while iter < maxiter:
          inc iter
          threads:
            var df2 = 0.0
            var f2 = 0.0
            for x in sub:
              let s = U[dir][x] * ds[dir][x].adj
              var t {.noinit.}: evalType(U[dir][x])
              t.projectTAH s
              df2 += norm2(t - fprev[dir][x]).simdSum
              f2 += t.norm2.simdSum
              fprev[dir][x] := t
              U[dir][x] := exp(-rho*t) * Upost[dir][x]
            threadBarrier()
            threadRankSum df2
            threadRankSum f2
            threadMaster:
              rdf2 = df2/f2
              if verbose:
                echo "stout inverse s", sweep, " p", parity, " d", dir,
                  " it", iter, " rdf2 ", rdf2
              if df2o >= 0.0 and df2o < df2:
                qexWarn "stout inverse df^2 increased (subset parity ", parity,
                  " dir ", dir, "): iter ", iter, " ", df2, " > ", df2o
              df2o = df2
          if rdf2 < rdf2req: break
        if rdf2 >= rdf2req:
          raiseValueError("stout inverse failed at sweep " & $sweep &
            ", parity " & $parity & ", dir " & $dir & " after " &
            $maxiter & " iterations: rdf2 " & $rdf2)
        result.iter += iter
        if rdf2 > result.rdf2: result.rdf2 = rdf2
