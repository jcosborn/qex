import qex
import alphas

import std/[json, strutils, sequtils]

const
  banner = """
|-----------------------------------------------------------------------|
 Quantum EXpressions (QEX) gradient flow

 QEX GitHub: https://github.com/jcosborn/qex

 QEX authors: James Osborn & Xiao-Yong Jin (Argonne National Laboratory)
 QEX gradient flow authors: 
   - James Osborn (Argonne National Laboratory)
   - Xiao-Yong Jin (Argonne National Laboratory)
 Modified and extended by:  
   - Curtis Taylor Peterson (Michigan State University)
     * Contact: curtistaylorpetersonwork@gmail.com
 Cite: Proceedings of Science (PoS) LATTICE2016 (2017) 271
|-----------------------------------------------------------------------|
"""

import alphas

import std/[strutils]
from std/math import round

let # defaults for testing/regression
  defaultLat = @[8, 8, 8, 16]
  lastFlowTime = (0.5*defaultLat[^1].float)*(0.5*defaultLat[^1].float)/8.0
  defaultMeasurements = ["plaquette", "rectangle", "clover", "topology", "polyakov"]
  defaultInputs = %* {
    "lattice-geometry": defaultLat,
    #"rank-geometry": [2, 2, 1, 2],
    #"simd-geometry": [2, 2, 1, 2],
    "flows": {
      "wilson": {
        "step-sizes": [0.02, 0.1],
        "step-size-transition-flow-times": [5.0, lastFlowTime],
        "measurements": defaultMeasurements
      }, 
      "symanzik": {
        "step-sizes": [0.02, 0.1],
        "step-size-transition-flow-times": [5.0, lastFlowTime],
        "measurements": defaultMeasurements
      },
      "zeuthen": {
        "step-sizes": [0.02, 0.1],
        "step-size-transition-flow-times": [5.0, lastFlowTime],
        "measurements": defaultMeasurements
      }
    }
  }

const CrSymanzik = -1.0/12.0

type
  GradientFlowKind = enum 
    WilsonFlow, 
    RectangleFlow, 
    IwasakiFlow, 
    DBW2Flow, 
    SymanzikFlow,
    ZeuthenFlow
  MeasurementKind = enum 
    Plaquette, 
    Clover, 
    Rectangle,
    Topology, 
    Polyakov

type
  Laplacian[U,U0] = object
    sf, sb: seq[Transporter[U,U,U0]]

type
  GradientFlow[U,U0] = object
    case kind: GradientFlowKind
      of ZeuthenFlow: 
        lap*: Laplacian[U,U0]
      else: discard
    gc: GaugeActionCoeffs
    meas: seq[MeasurementKind]
    u*: seq[U]
    step*: int
    flowTime*: float
    logFile*: File

#[ constructors ]#

converter strToGradientFlowKind(str: string): GradientFlowKind =
  case str: # I'm only dealing with these flows at the moment
    of "wilsonflow", "WilsonFlow", "wilson", "Wilson": WilsonFlow
    of "symanzikflow", "SymanzikFlow", "symanzik", "Symanzik": SymanzikFlow
    of "zeuthenflow", "ZeuthenFlow", "zeuthen", "Zeuthen": ZeuthenFlow
    else:
      qexError "unsupported flow kind: " & str
      WilsonFlow

converter strToMeasurementKind(str: string): MeasurementKind =
  case str:
    of "plaquette", "Plaquette", "plaq", "Plaq": Plaquette
    of "clover", "Clover", "clov", "Clov": Clover
    of "rectangle", "Rectangle", "rect", "Rect": Rectangle
    of "topology", "Topology", "topo", "Topo": Topology
    of "polyakov", "Polyakov", "poly", "Poly": Polyakov
    else:
      qexError "unsupported measurement kind: " & str
      Plaquette

proc newLaplacian[U](u: openArray[U]): auto =
  type U0 = evalType(u[0][0])
  result = Laplacian[U,U0](
    sf: u.newTransporters(u[0], 1),
    sb: u.newTransporters(u[0], -1)
  )

proc newGradientFlow[U](
  u: seq[U]; 
  kind: GradientFlowKind;
  plaq: float = 1.0;
  rect: float = CrSymanzik;
  meas: seq[MeasurementKind] = @[
    Plaquette,
    Rectangle, 
    Clover, 
    Topology, 
    Polyakov
  ]
): auto =
  type U0 = evalType(u[0][0])
  result = GradientFlow[U,U0](kind: kind, meas: meas, u: u.newOneOf())
  result.gc = case result.kind
    of WilsonFlow: GaugeActionCoeffs(plaq: plaq)
    of RectangleFlow: gaugeActRect(plaq, rect)
    of IwasakiFlow: Iwasaki(plaq)
    of DBW2Flow: DBW2(plaq)
    of SymanzikFlow: Symanzik(plaq)
    of ZeuthenFlow: 
      result.lap = u.newLaplacian()
      Symanzik(plaq)

#[ main gradient flow templates/procedures ]#

template ioBlock(flow: var GradientFlow; filename: string; body: untyped): untyped =
  if flow.logFile.isNil: flow.logFile = filename.open(fmWrite)
  else: qexError "log file already opened"
  body
  if not flow.logFile.isNil: flow.logFile.close()

template gradient[U](flow: var GradientFlow; f: var seq[U]; u: seq[U]) =
  let nd = u[0].l.nDim
  let ft = f.newOneOf()
  flow.gc.gaugeForce(u, f)
  case flow.kind:
    of ZeuthenFlow: # arXiv:1411.6706; Eqn. (6.2)
      flow.lap.sf.setLinks(u)
      flow.lap.sb.setLinks(u)
      let coeff = -0.5*CrSymanzik
      var
        sf = flow.lap.sf
        sb = flow.lap.sb
      var ft = f.newOneOf()
      threads:
        ft := f
        threadBarrier()
        for mu in 0..<nd:
          f[mu] += coeff*(sf[mu]^*ft[mu] + sb[mu]^*ft[mu] - 2.0*ft[mu])
    else: discard

proc isIn(meas: MeasurementKind, ms: seq[MeasurementKind]): bool =
  for m in ms:
    if meas == m: return true
  return false

template measurements[U,U0](flow: var GradientFlow[U,U0]; u: seq[U]): untyped =
  let prec = 18
  let 
    nc = u[0][0].nrows
    nd = u[0].l.nDim
    physVol = u[0].l.physVol
  let meas = flow.meas
  var output = newSeq[string]()
  var 
    fmunu: seq[seq[U]]
    es, et: U
  if Clover.isIn(meas) or Topology.isIn(meas): fmunu = u.fmunu(1)
  output.add "FLOW " & flow.flowTime.formatFloat(ffDecimal, 3)
  for meas in flow.meas:
    case meas:
      of Plaquette:
        let 
          plq = u.plaq
          nl = plq.len div 2
        let 
          ss = plq[0..<nl].sum
          st = plq[nl..^1].sum
        output.add (ss + st).formatFloat(ffDecimal, prec)
      of Rectangle: 
        let 
          frct = 1.0/physVol.float/nd.float/(nd.float - 1.0)
          rc = GaugeActionCoeffs(plaq: 0.0, rect: -frct)
          rect = rc.gaugeAction1(u)
        output.add rect.formatFloat(ffDecimal, prec)
      of Clover:
        let (es, et) = fmunu.densityE()
        let e = (nc.float - es - et)/nc.float
        output.add e.formatFloat(ffDecimal, prec)
      of Polyakov:
        type P = typeOf(u.wline @[1])
        let geom = u[0].l.physGeom
        var poly = newSeq[P](nd)
        for l in 0..<nd: poly[l] = u.wline repeat(l+1, geom[l])
        let
          pls = poly[0..^(nd-2)].sum() / (nd.float - 1.0)
          plt = poly[^1]
        output.add pls.re().formatFloat(ffDecimal, prec)
        output.add pls.im().formatFloat(ffDecimal, prec)
        output.add plt.re().formatFloat(ffDecimal, prec)
        output.add plt.im().formatFloat(ffDecimal, prec)
      of Topology: output.add fmunu.topoQ().formatFloat(ffDecimal, prec)
  let outputStr = output.join(" ") & "\n"
  flow.logFile.write(outputStr)
  echo outputStr.replace("FLOW", $flow.kind)

template gradientFlow[U,U0](flow: var GradientFlow[U,U0]; steps: int; eps: float) = 
  ## Gradient flow
  ## Originally written by James Osborn & Xiaoyong Jin.
  ## d/dt Vt = Z(Vt) Vt
  ## Runge-Kutta:
  ## W0 <- Vt
  ## W1 <- exp(1/4 Z0) W0
  ## W2 <- exp(8/9 Z1 - 17/36 Z0) W1
  ## V(t+eps) <- exp(3/4 Z2 - 8/9 Z1 + 17/36 Z0) W2
  ## where Zi = eps Z(Wi)
  let 
    nc = flow.u[0][0].nrows.float
    nd = flow.u[0].l.nDim
  let epsnc = eps * nc  # compensate force normalization
  var
    p = flow.u[0].l.newGauge  # mom
    f = flow.u[0].l.newGauge  # force
  var n = 0
  var u = flow.u.newOneOf()
  
  # gradient flow evolution and measurements
  threads: u := flow.u
  while true: # arXiv:1006.4518
    flow.gradient(f, u)
    threads:
      for mu in 0..<f.len:
        for e in u[mu]:
          var v {.noinit.}: flow.U0
          v := (-1.0/4.0)*epsnc*f[mu][e]
          let t = exp(v)*u[mu][e]
          p[mu][e] := v
          u[mu][e] := t
    flow.gradient(f, u)
    threads:
      for mu in 0..<f.len:
        for e in u[mu]:
          var v {.noinit.}: flow.U0
          v := (-8.0/9.0)*epsnc*f[mu][e] + (-17.0/9.0)*p[mu][e]
          let t = exp(v)*u[mu][e]
          p[mu][e] := v
          u[mu][e] := t
    flow.gradient(f, u)
    threads:
      for mu in 0..<f.len:
        for e in u[mu]:
          var v {.noinit.}: flow.U0
          v := (-3.0/4.0)*epsnc*f[mu][e] - p[mu][e]
          let t = exp(v)*u[mu][e]
          u[mu][e] := t
    inc n
    flow.step += 1
    flow.flowTime += eps
    flow.measurements(u)
    if n > (steps - 1): 
      threads: flow.u := u
      break

when isMainModule:
  echo banner

  proc getMeasurements(self: JsonNode): seq[MeasurementKind] =
    result = newSeq[MeasurementKind]()
    if self.hasKey("measurements"):
      for m in self["measurements"].items: 
        result.add m.getStr().strToMeasurementKind()
    else: qexError "no measurements specified for flow"
  
  qexInit()

  # command line information
  let 
    cmd = readCMD()
    flowInfo = case cmd.hasKey("json")
      of true: readJSON(cmd["json"].getStr())
      of false:
        echo "using default lattice parameters"
        defaultInputs
    cfg = case cmd.hasKey("configuration")
      of true: $cmd["configuration"].getInt()
      of false:
        echo "no configuration specified, using default"
        "0"
    filename = case cmd.hasKey("base-filename")
      of true: cmd["base-filename"].getStr()
      of false: "__REGRESSIONMODE__"
    outputDir = case cmd.hasKey("output-dir")
      of true: cmd["output-dir"].getStr()
      of false: ""
    latLayout = case flowInfo.hasKey("lattice-geometry")
      of true: flowInfo["lattice-geometry"].getIntSeq()
      of false: 
        echo "using default lattice geometry"
        defaultLat
    lo = case flowInfo.hasKey("rank-geometry")
      of true:
        let rnkLayout = flowInfo["rank-geometry"].getIntSeq()
        case flowInfo.hasKey("simd-geometry"):
          of true:
            let vecLayout = flowInfo["simd-geometry"].getIntSeq()
            newLayout(latLayout, VLEN, rnkLayout, vecLayout)
          of false: newLayout(latLayout, rnkLayout)
      of false: newLayout(latLayout)
  var u = lo.newGauge()

  # echo information
  echo "command line info: ", $cmd
  echo "gradient flow info: ", $flowInfo

  # read in gauge field
  if filename != "__REGRESSIONMODE__":
    u.readGauge(filename & "_" & cfg & ".lat")
  else: alphasrng.random(u) # for testing/regression

  # cycle through flows
  for flow in flowInfo["flows"].pairs:
    let 
      flowKind = flow.key.strToGradientFlowKind()
      measKinds = flow.val.getMeasurements()
      stepSizes = flow.val["step-sizes"].getFloatSeq()
      maxFlowTimes = flow.val["step-size-transition-flow-times"].getFloatSeq()
      transitions = maxFlowTimes.len
    let outputFilename = outputDir & flow.key & "_" & cfg & ".log"
    var flowObject = u.newGradientFlow(flowKind, meas = measKinds)

    # run flow for specified step sizes
    if stepSizes.len != transitions:
      qexError "number of step sizes must match number of transition flow times"
    flowObject.ioBlock(outputFilename):
      threads: flowObject.u := u
      flowObject.measurements(u) # t/a^2 = 0.0 measurement
      for transition in 0..<transitions:
        let 
          eps = stepSizes[transition]
          steps = round((maxFlowTimes[transition] - flowObject.flowTime)/eps)
        flowObject.gradientFlow(steps.int, eps)

  qexFinalize()