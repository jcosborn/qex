## Measurements and jackknife-analysis helpers shared by the 2D U(1) HMC apps,
## including the full `obstat` jackknife report. Field-transformation apps
## pass their per-trajectory ln det f'(V) as `Jvals`; `statsByQ` also accepts each
## trajectory's MD-force aggregates as `mdvals`.

import qex
import math
import maths/special
import utils/resample
import integrator
import ./stats

proc tunnelingRate*(xs: Ensemble[seq[float]]): float =
  var r = 0.0
  for i in 1..<xs.len:
    let d = xs[i]-xs[i-1]
    r += (d*d-r)/float(i)
  sqrt(r)

proc topo2DU1*(g: array or seq): float =
  ## Integer topological charge of a 2D U(1) configuration: the sum over sites of
  ## the plaquette angle, divided by 2*pi. `g` is a gauge snapshot (seq of 1x1
  ## color-matrix fields).
  let
    lo = g[0].l
    nd = lo.nDim
    t = newTransporters(g, g[0], 1)
  var p = 0.0
  threads:
    var tp: type(atan2(g[0][0][0,0].im, g[0][0][0,0].re))
    for mu in 1..<nd:
      for nu in 0..<mu:
        let tpl = (t[mu]^*g[nu]) * (t[nu]^*g[mu]).adj
        for i in tpl:
          tp += atan2(tpl[i][0,0].im, tpl[i][0,0].re)
    var v = tp.simdSum
    v.threadRankSum
    threadSingle: p += v
  p/TAU

proc topoMaxP2DU1*(g: array or seq): tuple[topo, maxP: float] =
  ## Topological charge (as `topo2DU1`) and the maximum |plaquette angle| over the
  ## lattice, in one pass. Used for proposed-configuration monitoring.
  let
    lo = g[0].l
    nd = lo.nDim
    t = newTransporters(g, g[0], 1)
  var
    psum = 0.0
    pmax = 0.0
  threads:
    var tp: type(atan2(g[0][0][0,0].im, g[0][0][0,0].re))
    var mx = 0.0
    for mu in 1..<nd:
      for nu in 0..<mu:
        let tpl = (t[mu]^*g[nu]) * (t[nu]^*g[mu]).adj
        for i in tpl:
          let a = atan2(tpl[i][0,0].im, tpl[i][0,0].re)
          tp += a
          let am = simdMax(abs(a))
          if am > mx: mx = am
    var v = tp.simdSum
    v.threadRankSum
    mx.threadRankMax
    threadSingle:
      psum += v
      pmax = mx
  (topo: psum/TAU, maxP: pmax)

proc infVolPlaq*(beta: float): float = besselI1(beta)/besselI0(beta)

proc infVolChiQ*(beta: float, tolerance = 1.0e-15): float =
  let i0beta = besselI0(beta)
  var
    term = PI * PI / 12.0
    sumVal = term
    n = 1
  while abs(term) > tolerance*abs(sumVal):
    let inbeta = besselIn(n, beta)
    # term = ((-1)^n / n^2) * (I_n(beta) / I_0(beta))
    term = inbeta / (i0beta * float(n*n))
    if n mod 2 == 0:
      sumVal += term
    else:
      sumVal -= term
    inc n
  sumVal / (PI * PI)

proc statsByQ*(Hvals, Avals: seq[float]; dQchanged: seq[bool]; jkBlockSize: int;
    Jvals: seq[float] = @[]; mdvals: seq[MdForceStats] = @[]) =
  ## Per-trajectory statistics split by whether the MD proposal changed the topological
  ## charge (`dQchanged[i]`, from the proposal monitor). `Hvals` contains signed dH.
  ## Each quantity is a grouped-jackknife mean ± standard error. Blocks are deleted
  ## from the original trajectory sequence before selecting the Q class: Pacc
  ## (⟨min(1,exp(−dH))⟩), exp(−dH), dHrms, and — when supplied — lnDet/lnDetrms (`Jvals`,
  ## ln det f'(V)) and the per-trajectory MD-force aggregates (`mdvals`). These
  ## mirror the overall obstat quantities; the overall Pacc is already printed by
  ## obstat, so only the two Q classes are shown here.
  var APvals = newSeq[float](Avals.len)
  for i in 0..<Avals.len:
    APvals[i] = min(1.0, Avals[i])
  proc report(tag: string; changed: bool) =
    var
      h: seq[float]
      selected = newSeq[bool](Avals.len)
      weights = newSeq[float](Avals.len)
    for i in 0..<Avals.len:
      if dQchanged[i] == changed:
        selected[i] = true
        weights[i] = 1.0
        h.add Hvals[i]
    echo tag, " ntraj: ", h.len
    if h.len == 0:
      return
    proc je(xs: seq[float]; isRms: bool): string =
      let s = xs.weightedJackknife(weights, jkBlockSize, isRms)
      $s.mean & " ± " & (if s.hasStdev: $s.stdev else: "n/a")
    let he = extrema(h)
    echo tag, " Pacc: ", je(APvals, false)
    echo tag, " exp(-dH): ", je(Avals, false)
    echo tag, " dHrms: ", je(Hvals, true)
    echo tag, " dH min/max: ", he.lo, " / ", he.hi
    if Jvals.len > 0:
      echo tag, " lnDet: ", je(Jvals, false)
      echo tag, " lnDetrms: ", je(Jvals, true)
    if mdvals.len > 0:
      echoMdStats(mdvals, jkBlockSize, tag, selected)
  report("dQ=0", false)
  report("dQ!=0", true)

proc obstat*(Hvals, Avals, Pvals, Qvals: seq[float]; beta: float; vol, ntraj, jkBlockSize: int;
    Jvals: seq[float] = @[]; mdvals: seq[MdForceStats] = @[]) =
  ## Jackknife analysis of the per-trajectory observables shared by the 2D U(1) HMC
  ## apps. `Hvals` contains signed dH. When `Jvals` (per-trajectory ln det f'(V))
  ## is non-empty, the field-transformation app reports lnDet/lnDetrms first; the
  ## plain sampler leaves `Jvals` empty and the report starts at dHrms.
  let dHrms = Hvals.jackknife(jkBlockSize, rms)
  let expmdh = Avals.jackknife(jkBlockSize, mean)
  var APvals = newseq[float](Avals.len)
  for i in 0..<Avals.len:
    let a = Avals[i]
    APvals[i] = if a > 1.0: 1.0 else: a
  let pacc = APvals.jackknife(jkBlockSize, mean)
  let Pmean = Pvals.jackknife(jkBlockSize, mean)
  let Pac = Pvals.jackknife(jkBlockSize, intAutocorr)
  let Qmean = Qvals.jackknife(jkBlockSize, mean)
  let Qac = Qvals.jackknife(jkBlockSize, intAutocorr, 0.0)
  let dQrms = Qvals.jackknife(jkBlockSize, tunnelingRate)

  var qmax = 0
  for qv in Qvals:
    let q = int(round(qv))
    if q > qmax:
      qmax = q
    elif -q > qmax:
      qmax = -q
  var qdist = newseq[JackknifeStat[float]]()
  proc count(xs: Ensemble[seq[float]], x: int): float =
    for i in 0..<xs.len:
      if x == int(abs(round(xs[i]))):
        result += 1.0
    if x != 0:
      result *= 0.5
  for i in 0..qmax:
    qdist.add Qvals.jackknife(jkBlockSize, count, i)

  var Q2vals = newseq[float](Qvals.len)
  for i in 0..<Qvals.len:
    let a = Qvals[i]
    Q2vals[i] = a*a
  let Q2 = Q2vals.jackknife(jkBlockSize, mean)
  let Q2ac = Q2vals.jackknife(jkBlockSize, intAutocorr)
  proc err(x: JackknifeStat[float]; scale = 1.0): string =
    if x.hasStdev: $(scale*x.stdev) else: "n/a"

  echo "all ntraj: ", Hvals.len
  if Jvals.len > 0:
    let lnDet = Jvals.jackknife(jkBlockSize, mean)
    let lnDetrms = Jvals.jackknife(jkBlockSize, rms)
    echo "all lnDet: ", lnDet.mean, " ± ", err(lnDet)
    echo "all lnDetrms: ", lnDetrms.mean, " ± ", err(lnDetrms)
  let he = extrema(Hvals)
  echo "all dHrms: ", dHrms.mean, " ± ", err(dHrms)
  echo "all dH min/max: ", he.lo, " / ", he.hi
  echo "all exp(-dH): ", expmdh.mean, " ± ", err(expmdh)
  echo "all Pacc: ", pacc.mean, " ± ", err(pacc)
  if mdvals.len > 0:
    echoMdStats(mdvals, jkBlockSize, "all")
  echo "all Pmean: ", Pmean.mean, " ± ", err(Pmean), "  dP: ", Pmean.mean - infVolPlaq(beta)
  echo "all Tau_P: ", Pac.mean, " ± ", err(Pac)
  echo "all Qmean: ", Qmean.mean, " ± ", err(Qmean)
  echo "all Tau_Q: ", Qac.mean, " ± ", err(Qac)
  echo "all dQrms: ", dQrms.mean, " ± ", err(dQrms)
  echo "all Q2/V: ", Q2.mean/float(vol), " ± ", err(Q2, 1.0/float(vol)), "  dQ2/V: ", Q2.mean/float(vol) - infVolChiQ(beta)
  echo "all Tau_Q2: ", Q2ac.mean, " ± ", err(Q2ac)
  for i in 0..qmax:
    echo "all P(Q=", i, "): ", qdist[i].mean/float(ntraj), " ± ", err(qdist[i], 1.0/float(ntraj))
