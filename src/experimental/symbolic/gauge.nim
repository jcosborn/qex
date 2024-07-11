import experimental/symbolic/graph
import layout, ../../gauge, physics/qcdTypes

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

method copySymNodeValue*(v: SymNodeValueGauge): SymNodeValue =
  # TODO: we don't need this, if we don't take gradient after assign
  let z = v.gaugeValue[0].l.newGauge
  threads:
    for mu in 0..<z.len:
      z[mu] := v.gaugeValue[mu]
  SymNodeValueGauge(gaugeValue: z)

method identAllocateSymNodeValue*(x: SymNodeValueGauge): SymNodeValue =
  # TODO: leave it uninitialized?
  SymNodeValueGauge(gaugeValue: x.gaugeValue[0].l.newGauge)

method identSymNodeValue*(z: SymNodeValueGauge, x: SymNodeValueGauge) =
  threads:
    for mu in 0..<z.gaugeValue.len:
      z.gaugeValue[mu] := x.gaugeValue[mu]

method zeroAllocateSymNodeValue*(x: SymNodeValueGauge): SymNodeValue =
  identAllocateSymNodeValue(x)

method zeroSymNodeValue*(z: SymNodeValueGauge, x: SymNodeValueGauge) =
  threads:
    for mu in 0..<z.gaugeValue.len:
      z.gaugeValue[mu] := 0.0

method addAllocateSymNodeValue*(x: SymNodeValueGauge, y: SymNodeValueGauge): SymNodeValue =
  SymNodeValueGauge(gaugeValue: x.gaugeValue[0].l.newGauge)

method addSymNodeValue*(z: SymNodeValueGauge, x: SymNodeValueGauge, y: SymNodeValueGauge) =
  threads:
    for mu in 0..<z.gaugeValue.len:
      z.gaugeValue[mu] := x.gaugeValue[mu] + y.gaugeValue[mu]

method mulAllocateSymNodeValue*(x: SymNodeValueGauge, y: SymNodeValueGauge): SymNodeValue =
  SymNodeValueGauge(gaugeValue: x.gaugeValue[0].l.newGauge)

method mulSymNodeValue*(z: SymNodeValueGauge, x: SymNodeValueGauge, y: SymNodeValueGauge) =
  threads:
    for mu in 0..<z.gaugeValue.len:
      z.gaugeValue[mu] := x.gaugeValue[mu] * y.gaugeValue[mu]

method mulAllocateSymNodeValue*(x: SymNodeValueFloat, y: SymNodeValueGauge): SymNodeValue =
  SymNodeValueGauge(gaugeValue: y.gaugeValue[0].l.newGauge)

method mulSymNodeValue*(z: SymNodeValueGauge, x: SymNodeValueFloat, y: SymNodeValueGauge) =
  threads:
    for mu in 0..<z.gaugeValue.len:
      z.gaugeValue[mu] := x.floatValue * y.gaugeValue[mu]

method mulAllocateSymNodeValue*(x: SymNodeValueGauge, y: SymNodeValueFloat): SymNodeValue =
  SymNodeValueGauge(gaugeValue: x.gaugeValue[0].l.newGauge)

method mulSymNodeValue*(z: SymNodeValueGauge, x: SymNodeValueGauge, y: SymNodeValueFloat) =
  threads:
    for mu in 0..<z.gaugeValue.len:
      z.gaugeValue[mu] := x.gaugeValue[mu] * y.floatValue

#[
method subAllocateSymNodeValue*(x: SymNodeValueFloat, y: SymNodeValueGauge): SymNodeValue =
  SymNodeValueGauge(gaugeValue: y.gaugeValue[0].l.newGauge)

method subSymNodeValue*(z: SymNodeValueGauge, x: SymNodeValueFloat, y: SymNodeValueGauge) =
  threads:
    for mu in 0..<z.gaugeValue.len:
      z.gaugeValue[mu] := x.floatValue - y.gaugeValue[mu]
]#

method subAllocateSymNodeValue*(x: SymNodeValueGauge, y: SymNodeValueFloat): SymNodeValue =
  SymNodeValueGauge(gaugeValue: x.gaugeValue[0].l.newGauge)

method subSymNodeValue*(z: SymNodeValueGauge, x: SymNodeValueGauge, y: SymNodeValueFloat) =
  threads:
    for mu in 0..<z.gaugeValue.len:
      z.gaugeValue[mu] := x.gaugeValue[mu] - y.floatValue

method subAllocateSymNodeValue*(x: SymNodeValueGauge, y: SymNodeValueGauge): SymNodeValue =
  SymNodeValueGauge(gaugeValue: x.gaugeValue[0].l.newGauge)

method subSymNodeValue*(z: SymNodeValueGauge, x: SymNodeValueGauge, y: SymNodeValueGauge) =
  threads:
    for mu in 0..<z.gaugeValue.len:
      z.gaugeValue[mu] := x.gaugeValue[mu] - y.gaugeValue[mu]

#
# more functions
#

proc norm2*(x: SymNode): SymNode
proc adjoint*(x: SymNode): SymNode

method norm2AllocateSymNodeValue*(x: SymNodeValue): SymNodeValue {.base.} =
  raiseErrorBaseMethod("args:\n  " & x.repr)

method norm2SymNodeValue*(z: SymNodeValue, x: SymNodeValue) {.base.} =
  raiseErrorBaseMethod("args:\n  " & z.repr & "\n  " & x.repr)

method norm2AllocateSymNodeValue*(x: SymNodeValueGauge): SymNodeValue =
  SymNodeValueFloat()

method norm2SymNodeValue*(z: SymNodeValueFloat, x: SymNodeValueGauge) =
  threads:
    var t = 0.0
    for mu in 0..<x.gaugeValue.len:
      t += x.gaugeValue[mu].norm2
    threadMaster: z.floatValue = t

proc norm2Allocate(z: SymNode) =
  z.value = norm2AllocateSymNodeValue(z.inputs[0].value)

proc norm2Forward(z: SymNode) =
  norm2SymNodeValue(z.value, z.inputs[0].value)

proc norm2Backward(z: SymNode, i: int, dep: SymNode): SymNode =
  if i != 0:
    raiseError("norm2 has 1 operand, got i = " & $i)
  let g = z.gradientDependentOrNil dep
  let two = newSymNodeFloat(floatValue = 2.0, name = "Two[norm2]")
  if g == nil:
    return two.mul z.inputs[0]
  else:
    let m = two.mul g
    return m.mul z.inputs[0]

proc norm2*(x: SymNode): SymNode =
  newSymNode(
    inputs = @[x],
    forward = norm2Forward,
    allocateValue = norm2Allocate,
    backward = norm2Backward,
    name = "norm2")

method adjointAllocateSymNodeValue*(x: SymNodeValue): SymNodeValue {.base.} =
  raiseErrorBaseMethod("args:\n  " & x.repr)

method adjointSymNodeValue*(z: SymNodeValue, x: SymNodeValue) {.base.} =
  raiseErrorBaseMethod("args:\n  " & z.repr & "\n  " & x.repr)

method adjointAllocateSymNodeValue*(x: SymNodeValueGauge): SymNodeValue =
  SymNodeValueGauge(gaugeValue: x.gaugeValue[0].l.newGauge)

method adjointSymNodeValue*(z: SymNodeValueGauge, x: SymNodeValueGauge) =
  threads:
    for mu in 0..<x.gaugeValue.len:
      z.gaugeValue[mu] := x.gaugeValue[mu].adj

proc adjointAllocate(z: SymNode) =
  z.value = adjointAllocateSymNodeValue(z.inputs[0].value)

proc adjointForward(z: SymNode) =
  adjointSymNodeValue(z.value, z.inputs[0].value)

proc adjointBackward(z: SymNode, i: int, dep: SymNode): SymNode =
  if i != 0:
    raiseError("adjoint has 1 operand, got i = " & $i)
  let g = z.gradientDependentOrNil dep
  if g == nil:
    raiseValueError("unimplemented")
  else:
    return g.adjoint

proc adjoint*(x: SymNode): SymNode =
  newSymNode(
    inputs = @[x],
    forward = adjointForward,
    allocateValue = adjointAllocate,
    backward = adjointBackward,
    name = "adjoint")

#
# gauge action
#

type SymNodeValueGaugeActionCoeffs* = ref object of SymNodeValueConcrete
  gaugeActionCoeffsValue*: GaugeActionCoeffs

method getGaugeActionCoeffs*(v: SymNodeValue): GaugeActionCoeffs {.base.} =
  raiseValueError("Custom method required for value: " & $v)

method getGaugeActionCoeffs*(v: SymNodeValueGaugeActionCoeffs): GaugeActionCoeffs =
  v.gaugeActionCoeffsValue

method setGaugeActionCoeffs*(v: SymNodeValue, c: GaugeActionCoeffs) {.base.} =
  raiseValueError("Custom method required for value: " & $v)

method setGaugeActionCoeffs*(v: SymNodeValueGaugeActionCoeffs, c: GaugeActionCoeffs) =
  v.gaugeActionCoeffsValue = c

proc gaugeAction*(c: SymNodeValue, beta: SymNode, g: SymNode): SymNode
proc gaugeActionDeriv*(c: SymNodeValue, beta: SymNode, g: SymNode): SymNode
proc gaugeActionDeriv2*(c: SymNodeValue, beta: SymNode, g: SymNode, f: SymNode): SymNode
proc projectTAH*(x: SymNode): SymNode
proc contractProjectTAH*(x: SymNode, y: SymNode): SymNode

proc gaugeForce*(c: SymNodeValue, beta: SymNode, g: SymNode): SymNode =
  contractProjectTAH(g, gaugeActionDeriv(c, beta, g))

method gaugeActionAllocateSymNodeValue*(x: SymNodeValue, y: SymNodeValue): SymNodeValue {.base.} =
  raiseErrorBaseMethod("args:\n  " & x.repr & "\n  " & y.repr)

method gaugeActionSymNodeValue*(z: SymNodeValue, c: SymNodeValue, x: SymNodeValue, y: SymNodeValue) {.base.} =
  raiseErrorBaseMethod("args:\n  " & z.repr & "\n  " & x.repr)

method gaugeActionSymNodeValue*(z: SymNodeValueFloat, c: SymNodeValueGaugeActionCoeffs, beta: SymNodeValueFloat, g: SymNodeValueGauge) =
  let gc = beta.floatValue * c.getGaugeActionCoeffs
  if gc.adjplaq == 0:
    z.floatValue = gc.gaugeAction1 g.gaugeValue
  elif gc.rect == 0 and gc.pgm == 0:
    z.floatValue = gc.actionA g.gaugeValue
  else:
    raiseValueError("Gauge coefficient unsupported: " & $c.getGaugeActionCoeffs)

method gaugeActionAllocateSymNodeValue*(beta: SymNodeValueFloat, g: SymNodeValueGauge): SymNodeValue =
  SymNodeValueFloat()

proc gaugeActionForward(z: SymNode) =
  gaugeActionSymNodeValue(z.value, z.getArg, z.inputs[0].value, z.inputs[1].value)

proc gaugeActionAllocate(z: SymNode) =
  z.value = gaugeActionAllocateSymNodeValue(z.inputs[0].value, z.inputs[1].value)

proc gaugeActionBackward(z: SymNode, i: int, dep: SymNode): SymNode =
  ## Always treat beta as a multiplicative factor
  let g = z.gradientDependentOrNil dep
  case i
  of 0:
    if g == nil:
      return z / z.inputs[0]
    else:
      return (g * z) / z.inputs[0]
  of 1:
    let f = gaugeForce(z.getArg, z.inputs[0], z.inputs[1])
    if g == nil:
      return f
    else:
      return g.mul f
  else:
    raiseError("gaugeAction has 2 operand, got i = " & $i)

proc gaugeAction*(c: SymNodeValue, beta: SymNode, g: SymNode): SymNode =
  newSymNode(
    inputs = @[beta, g],
    arg = c,
    forward = gaugeActionForward,
    allocateValue = gaugeActionAllocate,
    backward = gaugeActionBackward,
    name = "gaugeAction")

method gaugeActionDerivAllocateSymNodeValue*(x: SymNodeValue, y: SymNodeValue): SymNodeValue {.base.} =
  raiseErrorBaseMethod("args:\n  " & x.repr & "\n  " & y.repr)

method gaugeActionDerivSymNodeValue*(z: SymNodeValue, c: SymNodeValue, x: SymNodeValue, y: SymNodeValue) {.base.} =
  raiseErrorBaseMethod("args:\n  " & z.repr & "\n  " & x.repr)

method gaugeActionDerivAllocateSymNodeValue*(beta: SymNodeValueFloat, g: SymNodeValueGauge): SymNodeValue =
  SymNodeValueGauge(gaugeValue: g.gaugeValue[0].l.newGauge)

method gaugeActionDerivSymNodeValue*(z: SymNodeValueGauge, c: SymNodeValueGaugeActionCoeffs, beta: SymNodeValueFloat, g: SymNodeValueGauge) =
  let gc = beta.floatValue * c.getGaugeActionCoeffs
  if gc.adjplaq == 0:
    gc.gaugeActionDeriv(g.gaugeValue, z.gaugeValue)
  elif gc.rect == 0 and gc.pgm == 0:
    gc.gaugeADeriv(g.gaugeValue, z.gaugeValue)
  else:
    raiseValueError("Gauge coefficient unsupported: " & $c.getGaugeActionCoeffs)

proc gaugeActionDerivForward(z: SymNode) =
  gaugeActionDerivSymNodeValue(z.value, z.getArg, z.inputs[0].value, z.inputs[1].value)

proc gaugeActionDerivAllocate(z: SymNode) =
  z.value = gaugeActionDerivAllocateSymNodeValue(z.inputs[0].value, z.inputs[1].value)

proc gaugeActionDerivBackward(z: SymNode, i: int, dep: SymNode): SymNode =
  ## Always treat beta as a multiplicative factor
  let g = z.gradientDependentOrNil dep
  if g == nil:
    raiseValueError("gradient of " & $dep & " with respect to " & $z & " does not exists")
  case i
  of 0:
    # beta derivative
    #return g.dot(z) / z.inputs[0]
    raiseValueError("unimplemented")
  of 1:
    return g.mul gaugeActionDeriv2(z.getArg, z.inputs[0], z.inputs[1], g)
  else:
    raiseError("gaugeForce has 2 operand, got i = " & $i)

proc gaugeActionDeriv*(c: SymNodeValue, beta: SymNode, g: SymNode): SymNode =
  newSymNode(
    inputs = @[beta, g],
    arg = c,
    forward = gaugeActionDerivForward,
    allocateValue = gaugeActionDerivAllocate,
    backward = gaugeActionDerivBackward,
    name = "gaugeActionDeriv")

method gaugeActionDeriv2AllocateSymNodeValue*(x: SymNodeValue, y: SymNodeValue, z: SymNodeValue): SymNodeValue {.base.} =
  raiseErrorBaseMethod("args:\n  " & x.repr & "\n  " & y.repr & "\n  " & z.repr)

method gaugeActionDeriv2SymNodeValue*(z: SymNodeValue, c:SymNodeValue, x: SymNodeValue, y: SymNodeValue, h: SymNodeValue) {.base.} =
  raiseErrorBaseMethod("args:\n  " & z.repr & "\n  " & c.repr & "\n  " & x.repr & "\n  " & y.repr & "\n  " & h.repr)

method gaugeActionDeriv2AllocateSymNodeValue*(beta: SymNodeValueFloat, g: SymNodeValueGauge, c: SymNodeValueGauge): SymNodeValue  =
  SymNodeValueGauge(gaugeValue: g.gaugeValue[0].l.newGauge)

method gaugeActionDeriv2SymNodeValue*(z: SymNodeValueGauge, c: SymNodeValueGaugeActionCoeffs, beta: SymNodeValueFloat, g: SymNodeValueGauge, h: SymNodeValueGauge) =
  let gc = beta.floatValue * c.getGaugeActionCoeffs
  if gc.adjplaq == 0 and gc.rect == 0 and gc.pgm == 0:
    gc.gaugeDerivDeriv2(g.gaugeValue, h.gaugeValue, z.gaugeValue)
  else:
    raiseValueError("Gauge coefficient unsupported: " & $c.getGaugeActionCoeffs)

proc gaugeActionDeriv2Forward(z: SymNode) =
  gaugeActionDeriv2SymNodeValue(z.value, z.getArg, z.inputs[0].value, z.inputs[1].value, z.inputs[2].value)

proc gaugeActionDeriv2Allocate(z: SymNode) =
  z.value = gaugeActionDeriv2AllocateSymNodeValue(z.inputs[0].value, z.inputs[1].value, z.inputs[2].value)

proc gaugeActionDeriv2Backward(z: SymNode, i: int, dep: SymNode): SymNode =
  raiseValueError("unimplemented")

proc gaugeActionDeriv2*(c: SymNodeValue, beta: SymNode, g: SymNode, f: SymNode): SymNode =
  ## Remember to set arg for GaugeActionCoeffs with beta factored out
  newSymNode(
    inputs = @[beta, g, f],
    arg = c,
    forward = gaugeActionDeriv2Forward,
    allocateValue = gaugeActionDeriv2Allocate,
    backward = gaugeActionDeriv2Backward,
    name = "gaugeActionDeriv2")

method contractProjectTAHAllocateSymNodeValue*(x: SymNodeValue, y: SymNodeValue): SymNodeValue {.base.} =
  raiseErrorBaseMethod("args:\n  " & x.repr & "\n  " & y.repr)

method contractProjectTAHSymNodeValue*(z: SymNodeValue, x: SymNodeValue, y: SymNodeValue) {.base.} =
  raiseErrorBaseMethod("args:\n  " & z.repr & "\n  " & x.repr)

method contractProjectTAHAllocateSymNodeValue*(x: SymNodeValueGauge, y: SymNodeValueGauge): SymNodeValue  =
  SymNodeValueGauge(gaugeValue: x.gaugeValue[0].l.newGauge)

method contractProjectTAHSymNodeValue*(z: SymNodeValueGauge, x: SymNodeValueGauge, y: SymNodeValueGauge) =
  let nd = z.gaugeValue.len
  type T = type(z.gaugeValue[0])
  let f = cast[ptr cArray[T]](unsafeAddr(z.gaugeValue[0]))
  let u = cast[ptr cArray[T]](unsafeAddr(x.gaugeValue[0]))
  let o = cast[ptr cArray[T]](unsafeAddr(y.gaugeValue[0]))
  threads:
    for mu in 0..<nd:
      for e in o[mu]:
        let s = u[mu][e]*o[mu][e].adj
        f[mu][e].projectTAH s

proc contractProjectTAHForward(z: SymNode) =
  contractProjectTAHSymNodeValue(z.value, z.inputs[0].value, z.inputs[1].value)

proc contractProjectTAHAllocate(z: SymNode) =
  z.value = contractProjectTAHAllocateSymNodeValue(z.inputs[0].value, z.inputs[1].value)

proc contractProjectTAHBackward(z: SymNode, i: int, dep: SymNode): SymNode =
  let g = z.gradientDependentOrNil dep
  if g == nil:
    raiseValueError("gradient of " & $dep & " with respect to " & $z & " does not exists")
  # TODO: how can we optimize it so these are fused in a single loop
  let gTAH = g.projectTAH
  case i
  of 0:
    return gTAH * z.inputs[1]
  of 1:
    return gTAH.adjoint * z.inputs[0]
  else:
    raiseError("contractProjectTAH has 2 operand, got i = " & $i)

proc contractProjectTAH*(x: SymNode, y: SymNode): SymNode =
  newSymNode(
    inputs = @[x, y],
    forward = contractProjectTAHForward,
    allocateValue = contractProjectTAHAllocate,
    backward = contractProjectTAHBackward,
    name = "contractProjectTAH")

method projectTAHAllocateSymNodeValue*(x: SymNodeValue): SymNodeValue {.base.} =
  raiseErrorBaseMethod("args:\n  " & x.repr)

method projectTAHSymNodeValue*(z: SymNodeValue, x: SymNodeValue) {.base.} =
  raiseErrorBaseMethod("args:\n  " & z.repr & "\n  " & x.repr)

method projectTAHAllocateSymNodeValue*(x: SymNodeValueGauge): SymNodeValue  =
  SymNodeValueGauge(gaugeValue: x.gaugeValue[0].l.newGauge)

method projectTAHSymNodeValue*(z: SymNodeValueGauge, x: SymNodeValueGauge) =
  let nd = z.gaugeValue.len
  type T = type(z.gaugeValue[0])
  let f = cast[ptr cArray[T]](unsafeAddr(z.gaugeValue[0]))
  let u = cast[ptr cArray[T]](unsafeAddr(x.gaugeValue[0]))
  threads:
    for mu in 0..<nd:
      for e in u[mu]:
        f[mu][e].projectTAH u[mu][e]

proc projectTAHForward(z: SymNode) =
  projectTAHSymNodeValue(z.value, z.inputs[0].value)

proc projectTAHAllocate(z: SymNode) =
  z.value = projectTAHAllocateSymNodeValue(z.inputs[0].value)

proc projectTAHBackward(z: SymNode, i: int, dep: SymNode): SymNode =
  let g = z.gradientDependentOrNil dep
  if g == nil:
    raiseValueError("gradient of " & $dep & " with respect to " & $z & " does not exists")
  case i
  of 0:
    return g.projectTAH
  else:
    raiseError("projectTAH has 1 operand, got i = " & $i)

proc projectTAH*(x: SymNode): SymNode =
  newSymNode(
    inputs = @[x],
    forward = projectTAHForward,
    allocateValue = projectTAHAllocate,
    backward = projectTAHBackward,
    name = "projectTAH")

#
# SymNode and float
#

when isMainModule:
  proc hamiltonian(gc: SymNodeValue, b, x, p: SymNode): SymNode =
    gaugeAction(gc, b, x) + 0.5 * p.norm2  # TODO: -16*volume(p)

  proc intABA(eps, x, p: SymNode, h: proc): (SymNode, SymNode) =
    var x = exp(0.5 * eps * h(x,p).gradient(p)) * x
    var p = p - eps * h(x,p).gradient(x)
    x = exp(0.5 * eps * h(x,p).gradient(p)) * x
    return (x, p)

  let
    b = newSym("beta")
    gc = SymNodeValueGaugeActionCoeffs()
    u = newSym("u")
    p = newSym("p")
    eps = newSym("eps")

  proc ham(x, p: SymNode): SymNode = hamiltonian(gc, b, x, p)

  let
    (uu, pp) = intABA(eps, u, p, ham)
    h = ham(u, p)
    hh = ham(uu, pp)
    dhdp = h.gradient p
    dhdu = h.gradient u
    s = h.inputs[0]
    k = h.inputs[1]

  let checkdhdp = norm2(dhdp - p)

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
    dt = 0.1
  echoParams()
  echo "rank ", myRank, "/", nRanks
  threads: echo "thread ", threadNum, "/", numThreads

  gc.setGaugeActionCoeffs(case gact
    of ActWilson: GaugeActionCoeffs(plaq: 1.0)
    of ActAdjoint: GaugeActionCoeffs(plaq: 1.0, adjplaq: adjFac)
    of ActRect: gaugeActRect(1.0, rectFac)
    of ActSymanzik: Symanzik(1.0)
    of ActIwasaki: Iwasaki(1.0)
    of ActDBW2: DBW2(1.0))
  let
    lo = lat.newLayout
    vol = lo.physVol

  echo gc
  var r = lo.newRNGField(MRG32k3a, seed)
  var R:MRG32k3a  # global RNG
  R.seed(seed, 987654321)

  var g = lo.newgauge
  var mom = lo.newgauge
  threads:
    g.unit
    mom.randomTAH r

  b.assign beta
  u.assign g
  p.assign mom
  eps.assign dt
  optimize(output = @[hh,h,dhdp,checkdhdp], variables = @[p,u,eps,b])

  hh.allocate
  h.allocate
  k.allocate
  dhdp.allocate
  checkdhdp.allocate
  echo "k:\n",k.treerepr
  echo "dhdp:\n",dhdp.treerepr
  echo "norm2(dhdp-p):\n",checkdhdp.treerepr
  echo "hh:\n",hh.treerepr

  hh.eval
  h.eval
  dhdp.eval
  checkdhdp.eval

  echo h.value.getFloat, "  ", s.value.getFloat, "  ", k.value.getFloat
  echo hh.value.getFloat, "  ", hh.inputs[0].value.getFloat, "  ", hh.inputs[1].value.getFloat
  if checkdhdp.value.getFloat < 1e-30 * k.value.getFloat:
    echo "Pass: check dhdp"
  else:
    echo "Failed: check dhdp, norm2(dhdp-p) = ", checkdhdp.value.getFloat

  echo "k:\n",k.treerepr
  echo "dhdp:\n",dhdp.treerepr
  echo "norm2(dhdp-p):\n",checkdhdp.treerepr

  threads:
    mom.randomTAH r
  p.updated

  h.eval
  dhdp.eval
  checkdhdp.eval

  echo h.value.getFloat, "  ", s.value.getFloat, "  ", k.value.getFloat
  if checkdhdp.value.getFloat < 1e-30 * k.value.getFloat:
    echo "Pass: check dhdp"
  else:
    echo "Failed: check dhdp, norm2(dhdp-p) = ", checkdhdp.value.getFloat

  echo "k:\n",k.treerepr
  echo "dhdp:\n",dhdp.treerepr
  echo "norm2(dhdp-p):\n",checkdhdp.treerepr
