import experimental/symbolic/graph
import ../../gauge, physics/qcdTypes

#
# gauge support
#

# TODO: needs more thought
type Gauge = seq[DLatticeColorMatrixV]

type SymNodeValueGauge = ref object of SymNodeValueConcrete
  gaugeValue: Gauge

method getGauge*(v: SymNodeValue): Gauge {.base.} =
  raiseValueError("Custom method required for value: " & $v)

method getGauge*(v:SymNodeValueGauge): Gauge = v.gaugeValue

method `$`*(v: SymNodeValueGauge): string = "gaugeValue"

proc assign*(z: SymNode, v: Gauge) =
  z.assign SymNodeValueGauge(gaugeValue: v)

method copySymNodeValue*(v: SymNodeValueGauge): SymNodeValueGauge =
  # TODO: we don't need this, if we don't take gradient after assign
  let z = v.gaugeValue[0].l.newGauge
  threads:
    for mu in 0..<z.len:
      z[mu] := v.gaugeValue[mu]
  SymNodeValueGauge(gaugeValue: z)

method identSymNodeValue*(z: SymNodeValueGauge, x: SymNodeValueGauge) =
  threads:
    for mu in 0..<z.gaugeValue.len:
      z.gaugeValue[mu] := x.gaugeValue[mu]

method identAllocateSymNodeValue*(z: SymNode, x: SymNodeValueGauge) =
  # TODO: leave it uninitialized?
  z.value = SymNodeValueGauge(gaugeValue: x.gaugeValue[0].l.newGauge)

method iaddSymNodeValue*(z: SymNodeValueGauge, x: SymNodeValueGauge, y: SymNodeValueGauge) =
  threads:
    for mu in 0..<z.gaugeValue.len:
      z.gaugeValue[mu] := x.gaugeValue[mu] + y.gaugeValue[mu]

method iaddAllocateSymNodeValue*(z: SymNode, x: SymNodeValueGauge, y: SymNodeValueGauge) =
  z.value = SymNodeValueGauge(gaugeValue: x.gaugeValue[0].l.newGauge)

method norm2SymNodeValue*(z: SymNodeValue, x: SymNodeValue) {.base.} =
  raiseErrorBaseMethod("args:\n  " & z.repr & "\n  " & x.repr)

method norm2AllocateSymNodeValue*(z: SymNode, x: SymNodeValue) {.base.} =
  raiseErrorBaseMethod("args:\n  " & z.nodeRepr & "\n  " & x.repr)

method imulSymNodeValue*(z: SymNodeValue, x: SymNodeValue, y: SymNodeValue) {.base.} =
  raiseErrorBaseMethod("args:\n  " & z.repr & "\n  " & x.repr & "\n  " & y.repr)

method imulAllocateSymNodeValue*(z: SymNode, x: SymNodeValue, y: SymNodeValue) {.base.} =
  raiseErrorBaseMethod("args:\n  " & z.nodeRepr & "\n  " & x.repr & "\n  " & y.repr)

method norm2SymNodeValue*(z: SymNodeValueFloat, x: SymNodeValueGauge) =
  threads:
    var t = 0.0
    for mu in 0..<x.gaugeValue.len:
      t += x.gaugeValue[mu].norm2
    threadMaster: z.floatValue = t

method norm2AllocateSymNodeValue*(z: SymNode, x: SymNodeValueGauge) =
  z.value = SymNodeValueFloat()

method imulSymNodeValue*(z: SymNodeValueFloat, x: SymNodeValueFloat, y: SymNodeValueFloat) =
  z.floatValue = x.floatValue * y.floatValue

method imulAllocateSymNodeValue*(z: SymNode, x: SymNodeValueFloat, y: SymNodeValueFloat) =
  z.value = SymNodeValueFloat()

method imulSymNodeValue*(z: SymNodeValueGauge, x: SymNodeValueFloat, y: SymNodeValueGauge) =
  threads:
    for mu in 0..<z.gaugeValue.len:
      z.gaugeValue[mu] := x.floatValue * y.gaugeValue[mu]

method imulAllocateSymNodeValue*(z: SymNode, x: SymNodeValueFloat, y: SymNodeValueGauge) =
  z.value = SymNodeValueGauge(gaugeValue: y.gaugeValue[0].l.newGauge)

method imulSymNodeValue*(z: SymNodeValueGauge, x: SymNodeValueGauge, y: SymNodeValueFloat) =
  threads:
    for mu in 0..<z.gaugeValue.len:
      z.gaugeValue[mu] := x.gaugeValue[mu] * y.floatValue

method imulAllocateSymNodeValue*(z: SymNode, x: SymNodeValueGauge, y: SymNodeValueFloat) =
  z.value = SymNodeValueGauge(gaugeValue: x.gaugeValue[0].l.newGauge)

method imulSymNodeValue*(z: SymNodeValueGauge, x: SymNodeValueGauge, y: SymNodeValueGauge) =
  threads:
    for mu in 0..<z.gaugeValue.len:
      z.gaugeValue[mu] := x.gaugeValue[mu] * y.gaugeValue[mu]

method imulAllocateSymNodeValue*(z: SymNode, x: SymNodeValueGauge, y: SymNodeValueGauge) =
  z.value = SymNodeValueGauge(gaugeValue: x.gaugeValue[0].l.newGauge)

#
# more functions
#

proc mul*(x: SymNode, y: SymNode): SymNode

proc mulForward(z: SymNode) =
  imulSymNodeValue(z.value, z.inputs[0].value, z.inputs[1].value)

proc mulAllocate(z: SymNode) =
  imulAllocateSymNodeValue(z, z.inputs[0].value, z.inputs[1].value)

proc mulBackward(z: SymNode, i: int, dep: SymNode): SymNode =
  if i != 0 and i != 1:
    raiseError("mul has 2 operands, got i = " & $i)
  let g = z.gradientDependentOrNil dep
  if g == nil:
    if i == 0:
      return z.inputs[1]
    else:
      return z.inputs[0]
  else:
    # assume noncommuting operands here
    if i == 0:
      return g.mul z.inputs[1]
    else:
      return z.inputs[0].mul g

proc mul*(x: SymNode, y: SymNode): SymNode =
  newSymNode(
    inputs = @[x, y],
    forward = mulForward,
    allocateValue = mulAllocate,
    backward = mulBackward,
    name = "mul")

proc norm2Forward(z: SymNode) =
  norm2SymNodeValue(z.value, z.inputs[0].value)

proc norm2Allocate(z: SymNode) =
  norm2AllocateSymNodeValue(z, z.inputs[0].value)

proc norm2Backward(z: SymNode, i: int, dep: SymNode): SymNode =
  if i != 0:
    raiseError("norm2 has 1 operand, got i = " & $i)
  let g = z.gradientDependentOrNil dep
  let two = newSymNode(value = SymNodeValueFloat(floatValue: 2.0), name = "Two[norm2]")
  if g == nil:
    return two.mul z.inputs[0]
  else:
    let m = two.mul g
    m.tagUpdate
    return m.mul z.inputs[0]

proc norm2*(x: SymNode): SymNode =
  newSymNode(
    inputs = @[x],
    forward = norm2Forward,
    allocateValue = norm2Allocate,
    backward = norm2Backward,
    name = "norm2")

#
# SymNode and float
#

# TODO perhaps make a generic converter: anyValue -> SymNode

proc `*`(x: float, y: SymNode): SymNode =
  return newSymNode(value = SymNodeValueFloat(floatValue: x)).mul y

proc `*`(x: SymNode, y: float): SymNode =
  return x.mul newSymNode(value = SymNodeValueFloat(floatValue: y))

when isMainModule:
  let p = newSym("p")
  let k = 0.5 * p.norm2
  let dkdp = k.gradient p

  import qex
  import os, strutils

  qexInit()
  type GaugeActType = enum ActWilson, ActAdjoint, ActRect, ActSymanzik, ActIwasaki, ActDBW2
  converter toGaugeActType(s:string):GaugeActType = parseEnum[GaugeActType](s)
  letParam:
    gaugefile = ""
    savefile = "config"
    savefreq = 10
    lat =
      if fileExists(gaugefile):
        getFileLattice gaugefile
      else:
        if gaugefile.len > 0:
          qexWarn "Nonexistent gauge file: ", gaugefile
        @[4,4,4,8]
    gact:GaugeActType = "ActWilson"
    beta = 6.0
    adjFac = -0.25
    rectFac = -1.4088
    seed:uint64 = 4321
  echoParams()
  echo "rank ", myRank, "/", nRanks
  threads: echo "thread ", threadNum, "/", numThreads

  let
    gc = case gact
      of ActWilson: GaugeActionCoeffs(plaq: beta)
      of ActAdjoint: GaugeActionCoeffs(plaq: beta, adjplaq: beta*adjFac)
      of ActRect: gaugeActRect(beta, rectFac)
      of ActSymanzik: Symanzik(beta)
      of ActIwasaki: Iwasaki(beta)
      of ActDBW2: DBW2(beta)
    lo = lat.newLayout
    vol = lo.physVol

  echo gc
  var r = lo.newRNGField(MRG32k3a, seed)
  var R:MRG32k3a  # global RNG
  R.seed(seed, 987654321)

  var mom = lo.newgauge
  threads:
    mom.randomTAH r

  p.assign mom
  optimize(output = @[k,dkdp], variables = @[p])

  k.allocate
  dkdp.allocate
  echo k.treerepr
  echo dkdp.treerepr

  k.eval
  dkdp.eval

  echo k.value
  let dp = dkdp.value.getGauge
  var d2 = 0.0
  threads:
    var d2t = 0.0
    for mu in 0..<mom.len:
      mom[mu] -= dp[mu]
      d2t += mom[mu].norm2
    threadMaster: d2 = d2t
  echo d2
