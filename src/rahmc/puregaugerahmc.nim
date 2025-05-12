# Implementation of repelling-attracting HMC (arXiv:2403.04607)
# GitHub: https://github.com/sidv23/ra-hmc/tree/4cae4af6873867f9a307a8a0b40f83d76ac56cc3

import qex, gauge, physics/qcdTypes, maths/groupOps, gauge/stoutsmear
import math, os, sequtils, strformat, strutils, times, typetraits, json, tables
import utilityProcs

type 
  GaugeActType = enum 
    ActWilson, 
    ActAdjoint, 
    ActRect, 
    ActSymanzik, 
    ActIwasaki, 
    ActDBW2

# ---------------------

qexInit()

let gact = ActSymanzik
var gamma: float

letParam:
  lat = @[16,16,16,16]
  beta = 4.7
  adjFac = -0.25
  nAdjFac = -0.25
  rectFac = -1.0/12.0
  tau = 1.0
  hmcSteps = 100
  rahmcSteps = 100
  rahmcTrajLenScaleFac = 20.0
  trajs = 500
  hmcStepSize = tau/hmcSteps
  rahmcStepSize = tau/rahmcSteps/rahmcTrajLenScaleFac
  seed: uint64 = 987654321
  noMetropolisUntil = 20
  hmcUntil = 40

let
  gc = case gact
    of ActWilson: GaugeActionCoeffs(plaq: beta)
    of ActAdjoint: GaugeActionCoeffs(plaq: beta, adjplaq: beta*adjFac)
    of ActRect: gaugeActRect(beta, rectFac)
    of ActSymanzik: Symanzik(beta)
    of ActIwasaki: Iwasaki(beta)
    of ActDBW2: DBW2(beta)

var
  lo = lat.newLayout
  vol = lo.physVol

var 
  r = lo.newRNGField(MRG32k3a, seed)
  g = lo.newGauge()
  bg = lo.newGauge()
  f = lo.newGauge()
  p = lo.newGauge()
  R: MRG32k3a 

R.seed(seed, 987654321)
g.unit

var
  maxFlowTime = (0.5*lat[^1])*(0.5*lat[^1])/8.0
  flowInfo = %* {
    "C0p0": {
      "action": "Wilson",
      "path": "./",
      "time-increments": [0.02,0.1],
      "maximum-flow-times": [5.0,maxFlowTime]
    }
  }

# ---------------------

proc `/`(x,y: int): int = x div y

proc flipDampingFactor() = 
  gamma = -gamma

proc heatbath() = 
  threads: p.randomTAH(r)

proc momentumFlip() =
  var q = lo.newGauge()
  threads:
    for mu in 0..<p.len:
      q[mu] := (-1.0)*p[mu]
      threadBarrier()
      p[mu] := q[mu]

proc action(): float =
  result = case gact
    of ActAdjoint: gc.actionA(g)
    else: gc.gaugeAction1(g)

proc hamiltonian(): float = 0.5*p.norm2 + action() - float(16*vol)

proc force() =
  case gact:
    of ActAdjoint: gc.forceA(g,f)
    else: gc.gaugeForce(g,f)

proc updatePConformal(b,eps: float) =
  var q = lo.newGauge()
  force()
  threads:
    for mu in 0..<p.len:
      q[mu] := p[mu] 
      threadBarrier()  
      p[mu] := b*q[mu] - eps*f[mu]

proc updateP(eps: float) =
  force()
  threads:
    for mu in 0..<p.len: p[mu] -= eps*f[mu]

proc updateU(eps: float) =
  threads:
    for mu in 0..<g.len:
      for s in g[mu]: g[mu][s] := exp(eps*p[mu][s])*g[mu][s]

proc evolveRAHMC() =
  let
    e = rahmcStepSize
    eby2 = 0.5*rahmcStepSize
  var
    b = exp(gamma*eby2)
    b2 = exp(gamma*e)

  # Steps with gamma < 0
  updatePConformal(b,eby2)
  for step in 1..<rahmcSteps/2:
    updateU(e)
    updatePConformal(b2,(1.0+b2)*eby2)
  updateU(e)
  updatePConformal(b,b*eby2)

  # Steps with gamma > 0
  g.reunit
  (b,b2) = (1/b,1/b2)
  updatePConformal(b,eby2)
  for step in 1..<rahmcSteps/2:
    updateU(e)
    updatePConformal(b2,(1.0+b2)*eby2)
  updateU(e)
  updatePConformal(b,b*eby2)

proc evolveHMC() =
  updateP(0.5*hmcStepSize)
  for step in 1..<hmcSteps:
    updateU(hmcStepSize)
    if step != hmcSteps-1: updateP(hmcStepSize)
    else: updateP(0.5*hmcStepSize)

# ---------------------

for traj in 0..<trajs:
  var 
    hi,hf,dh: float
    alg = case (traj < hmcUntil)
      of true: "HMC"
      of false: "raHMC"

  gamma = case (traj < hmcUntil)
    of true: 0.0
    of false: R.uniform
  
  bg.setGauge(g)
  if traj < hmcUntil: heatbath()
  else: momentumFlip()
  hi = hamiltonian()

  if traj < hmcUntil: evolveHMC()
  else: evolveRAHMC()
  
  hf = hamiltonian()
  dh = hf-hi
  if traj > noMetropolisUntil:
    if R.uniform <= exp(-dh):
      echo alg," (",traj,")"," ACCEPT: dH: ", dH
      g.reunit
    else:
      echo alg," (",traj,")"," REJECT: dH: ", dH
      g.setGauge(bg)
  else:
    echo alg," (",traj,")"," ACCEPT (NO METROPOLIS): dH: ", dH
    g.reunit

  g.mplaq
  g.ploop

  let filename = "C0p0_" & $(traj) & ".log"
  flowInfo["C0p0"].setFilename(filename)
  g.gradientFlow(flowInfo):
    let output = measurements.formatMeasurements(style="KS_nHYP_FA")
    f.write(output & "\n")


qexFinalize()

# ---------------------

#[
  # Steps with gamma < 0
  force()
  threads:
    for mu in 0..<p.len: p[mu] := b*p[mu] - eby2*f[mu]
  for step in 1..<steps:
    threads:
      for mu in 0..<g.len:
        for s in g[mu]: g[mu][s] := exp(e*p[mu][s])*g[mu][s]
    force()
    threads:
      for mu in 0..<p.len: p[mu] := b2*p[mu] - (1.0 + b2)*eby2*f[mu]
  threads:
    for mu in 0..<g.len:
      for s in g[mu]: g[mu][s] := exp(e*p[mu][s])*g[mu][s]
  force()
  threads:
    for mu in 0..<p.len: p[mu] := b*(p[mu] - eby2*f[mu])
]#