## Calculates tree-level normalization factor for gradient flow
## 
## This code is both a very useful tool and a demonstration of what's possible with
## QEX and the Nim programming language. 
## 
## Based on: JHEP09(2014)018
## Validated against: https://github.com/akhilc2/gradientflow/tree/main
## 
## Author: Curtis Taylor Peterson

import std/[strutils]
import std/[strformat]
import std/[math]

import qex

const nd = 4

type
  GaugeActionKind* = enum
    Symanzik,
    Clover

  GaugeAction* = object
    case kind: GaugeActionKind
    of Symanzik:
      cr: float
    of Clover: discard
    
  GaugeFlow*[V: static[int]] = object
    l*: Layout[V]
    g*: GaugeAction
    f*: GaugeAction
    e*: GaugeAction
    alpha*: float

type 
  Coordinate[V: static[int]] = array[V, float]
  CoordinateV[V: static[int]] = array[nd, Coordinate[V]]

#[ gauge action & gauge flow ]#

proc newGaugeAction*(kind: GaugeActionKind; cr: float = 0.0): GaugeAction =
  case kind
  of Symanzik:
    result.kind = kind
    result.cr = cr
  of Clover: result.kind = kind

proc newGaugeFlow*[V: static[int]](
  lo: Layout[V]; 
  g, f, e: GaugeAction;
  alpha: float = 1.0
): GaugeFlow[V] =
  result.l = lo
  result.g = g
  result.f = f
  result.e = e
  result.alpha = alpha

proc summarize(action: GaugeAction): string =
  result = case action.kind:
    of Symanzik: fmt"Symanzik (rectangle coefficient: {action.cr})"
    of Clover: "Clover"

proc summarize*(flow: GaugeFlow): void =
  echo "Gauge action: ", flow.g.summarize()
  echo "Flow action: ", flow.f.summarize()
  echo "Operator action: ", flow.e.summarize()
  echo "Gauge fixing parameter: ", fmt"{flow.alpha}"

#[ coordinate & coordinate arithmetic ]#

proc newCoordinateV[V: static[int]](lo: Layout[V]; s: int): CoordinateV[V] =
  for mu in 0..<nd:
    for lane in 0..<V: result[mu][lane] = lo.coords[mu][s*V + lane].float

proc `+=`[V: static[int]](a: var Coordinate[V]; b: Coordinate[V]) =
  for lane in 0..<V: a[lane] += b[lane]

proc `-`[V: static[int]](a: Coordinate[V]): Coordinate[V] =
  for lane in 0..<V: result[lane] = -a[lane]

proc `-`[V: static[int]](a, b: Coordinate[V]): Coordinate[V] =
  for lane in 0..<V: result[lane] = a[lane] - b[lane]

proc `-`[V: static[int]](a: float; b: Coordinate[V]): Coordinate[V] =
  for lane in 0..<V: result[lane] = a - b[lane]

proc `*`[V: static[int]](a: Coordinate[V]; b: Coordinate[V]): Coordinate[V] =
  for lane in 0..<V: result[lane] = a[lane]*b[lane]

proc `*`[V: static[int]](a: float; b: Coordinate[V]): Coordinate[V] = 
  for lane in 0..<V: result[lane] = a*b[lane]

proc `/`[V: static[int]](a: Coordinate[V]; b: float): Coordinate[V] =
  for lane in 0..<V: result[lane] = a[lane]/b

proc sin[V: static[int]](a: Coordinate[V]): Coordinate[V] =
  for lane in 0..<V: result[lane] = sin(a[lane])

proc cos[V: static[int]](a: Coordinate[V]): Coordinate[V] =
  for lane in 0..<V: result[lane] = cos(a[lane])

proc `-=`[V: static[int]](a: var DComplexV; b: Coordinate[V]) =
  for lane in 0..<V:
    a.re[asSimd(lane)] -= b[lane]
    a.im[lane] = 0.0

proc `:=`[V: static[int]](a: var DComplexV; b: Coordinate[V]) =
  for lane in 0..<V:
    a.re[lane] = b[lane]
    a.im[lane] = 0.0

#[ tree-level normalization factor ]#

template action(action: GaugeAction; s: int): MatrixArray[nd, nd, DComplexV] =
  block:
    var result: MatrixArray[nd, nd, DComplexV]
    var p = lo.newCoordinateV(s)
    var p2: Coordinate[V]
    case action.kind
    of Symanzik:
      let cr = action.cr
      var p4: Coordinate[V]
      for mu in 0..<nd:
        var p2mu {.noinit.}: Coordinate[V]
        p[mu] = 2.0*sin(PI*p[mu]/lo.physGeom[mu].float)
        p2mu = p[mu]*p[mu]
        p2 += p2mu
        p4 += p2mu*p2mu
      for mu in 0..<nd:
        for nu in 0..<nd:
          if mu == nu: result[mu, nu] := p2 - cr*p4 - cr*p2*p[mu]*p[mu]
          else: result[mu, nu] := 0.0
          result[mu, nu] -= p[mu]*p[nu]*(1.0 - cr*p[mu]*p[mu] - cr*p[nu]*p[nu])
    of Clover:
      var q: CoordinateV[V]
      for mu in 0..<nd:
        q[mu] = cos(PI*p[mu]/lo.physGeom[mu].float)
        p[mu] = sin(2.0*PI*p[mu]/lo.physGeom[mu].float)
        p2 += p[mu]*p[mu]
      for mu in 0..<nd:
        for nu in 0..<nd:
          var kmunu: Coordinate[V] = -p[mu]*p[nu]
          if mu == nu: kmunu += p2
          result[mu, nu] := kmunu*q[mu]*q[nu]
    result

template gaugeFix[V](flow: GaugeFlow[V]; s: int): MatrixArray[nd, nd, DComplexV] =
  block:
    var result: MatrixArray[nd, nd, DComplexV]
    var p = lo.newCoordinateV(s)
    let alpha = flow.alpha
    for mu in 0..<nd: p[mu] = 2.0*sin(PI*p[mu]/lo.physGeom[mu].float)
    for mu in 0..<nd:
      for nu in 0..<nd: result[mu, nu] := p[mu]*p[nu]/alpha
    result

template inverse(a: MatrixArray[nd, nd, DComplexV]): MatrixArray[nd, nd, DComplexV] =
  block:
    var result: MatrixArray[nd, nd, DComplexV]
    result.inverse(a)
    result

proc treeLevelNormalization*[V: static[int]](flow: GaugeFlow[V]; t: float): float =
  assert flow.l.nDim == nd
  let lo = flow.l
  var sf, sg, se: Field[V, MatrixArray[nd, nd, DComplexV]]
  var sum: float = 0.0

  sf.new flow.l
  sg.new flow.l
  se.new flow.l

  threads:
    for s in sf:
      var gf: MatrixArray[nd, nd, DComplexV] = flow.gaugeFix(s)
      sf[s] := exp(-t*(flow.f.action(s) + gf))
      sg[s] := inverse(flow.g.action(s) + gf)
      se[s] := flow.e.action(s)
    threadBarrier()
    threadMaster:
      if myRank == 0: sg{0} := 0.0
    threadBarrier()
    var tSum: float = 0.0
    for s in sf: tSum += trace(sf[s]*sg[s]*sf[s]*se[s]).simdSum().re
    threadBarrier()
    threadSum(tSum)
    threadMaster: sum = tSum
  rankSum(sum)

  return 64.0*PI*PI*t*t*(2.0 + sum)/3.0/lo.physVol.float

qexInit()

letParam:
  lattice = @[24, 24, 24, 48]
  
  gauge = "Symanzik"  # default: Luschwer-Weiss gauge action
  flow = "Symanzik"   # default: Wilson flow
  operator = "Clover" # default: clover E(t) discretization

  gaugeRectangleCoeff = -1.0/12.0
  flowRectangleCoeff = 0.0
  operatorRectangleCoeff = -1.0/12.0

  gaugeFixCoeff = 1.0

  minimumFlowTime = 0.0
  maximumFlowTime = (0.5*lattice[0].float)*(0.5*lattice[0].float) / 8.0
  flowTimeEpsilon = 0.02

let lo = lattice.newLayout()
let gaugeFlow = lo.newGaugeFlow(
  parseEnum[GaugeActionKind](gauge).newGaugeAction(gaugeRectangleCoeff),
  parseEnum[GaugeActionKind](flow).newGaugeAction(flowRectangleCoeff),
  parseEnum[GaugeActionKind](operator).newGaugeAction(operatorRectangleCoeff),
  alpha = gaugeFixCoeff
)

gaugeFlow.summarize()

echo ""
echo "    t/a^2    C(a^2/t, sqrt{8t}/L)"
echo "------------ --------------------"
var flowTime = minimumFlowTime
while round(flowTime, 10) <= round(maximumFlowTime, 10):
  let normalization = gaugeFlow.treeLevelNormalization(flowTime)
  echo fmt"{flowTime:.10f}", " ", fmt"{normalization:.18f}"
  flowTime += flowTimeEpsilon
echo ""

qexFinalize()
