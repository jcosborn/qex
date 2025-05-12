import qex
import sequtils, json, strutils, io, sets
import gauge
import gauge/gaugeUtils

proc norm2*(p:auto): float =
  var output: float
  threads:
    var temp = 0.0
    for mu in 0..<p.len: temp += p[mu].norm2
    threadBarrier()
    threadMaster: output = temp
  output

proc reunit*(g:auto) =
  tic()
  threads:
    let d = g.checkSU
    threadBarrier()
    echo "unitary deviation avg: ",d.avg," max: ",d.max
    g.projectSU
    threadBarrier()
    let dd = g.checkSU
    echo "new unitary deviation avg: ",dd.avg," max: ",dd.max
  toc("reunit")

proc mplaq*(g:auto) =
  tic()
  let
    pl = g.plaq
    nl = pl.len div 2
    ps = pl[0..<nl].sum * 2.0
    pt = pl[nl..^1].sum * 2.0
  echo "MEASplaq ss: ",ps,"  st: ",pt,"  tot: ",0.5*(ps+pt)
  toc("plaq")

proc ploop*(g:auto) =
  tic()
  let pg = g[0].l.physGeom
  var pl = newseq[typeof(g.wline @[1])](pg.len)
  for i in 0..<pg.len:
    pl[i] = g.wline repeat(i+1, pg[i])
  let
    pls = pl[0..^2].sum / float(pl.len-1)
    plt = pl[^1]
  echo "MEASploop spatial: ",pls.re," ",pls.im," temporal: ",plt.re," ",plt.im
  toc("ploop")

proc setGauge*(u: auto; v: auto) =
  threads:
    for mu in 0..<u.len:
      for s in u[mu]: u[mu][s] := v[mu][s]


template gradientFlow(
    gc: GaugeActionCoeffs;
    g: array|seq;
    steps: int;
    eps: float;
    measure: untyped
  ): untyped =
  #[ Gradient flow w/ Wilson or rectangle action
     Originally written by James Osborn & Xiaoyong Jin.
     d/dt Vt = Z(Vt) Vt
     Runge-Kutta:
     W0 <- Vt
     W1 <- exp(1/4 Z0) W0
     W2 <- exp(8/9 Z1 - 17/36 Z0) W1
     V(t+eps) <- exp(3/4 Z2 - 8/9 Z1 + 17/36 Z0) W2
     where
     Zi = eps Z(Wi)
  ]#
  proc flowProc {.gensym.} =
    tic("flowProc")
    const nc = g[0][0].nrows.float
    var
      p = g[0].l.newGauge  # mom
      f = g[0].l.newGauge  # force
      n = 1
    while true:
      let t = n * eps
      let epsnc = eps * nc  # compensate force normalization
      gc.gaugeForce(g,f)
      threads:
        for mu in 0..<f.len:
          for e in g[mu]:
            var v {.noinit.}:type(load1(f[0][0]))
            v := (-1.0/4.0)*epsnc*f[mu][e]
            let t = exp(v)*g[mu][e]
            p[mu][e] := v
            g[mu][e] := t
      gc.gaugeForce(g,f)
      threads:
        for mu in 0..<f.len:
          for e in g[mu]:
            var v {.noinit.}:type(load1(f[0][0]))
            v := (-8.0/9.0)*epsnc*f[mu][e] + (-17.0/9.0)*p[mu][e]
            let t = exp(v)*g[mu][e]
            p[mu][e] := v
            g[mu][e] := t
      gc.gaugeForce(g,f)
      threads:
        for mu in 0..<f.len:
          for e in g[mu]:
            var v {.noinit.}:type(load1(f[0][0]))
            v := (-3.0/4.0)*epsnc*f[mu][e] - p[mu][e]
            let t = exp(v)*g[mu][e]
            g[mu][e] := t
      let wflowT {.inject.} = t
      measure
      inc n
      if steps>0 and n>steps: break
    toc("end")
  flowProc()

template gradientFlow(
    gc: GaugeActionCoeffs;
    g: array|seq;
    eps: float;
    measure: untyped
  ): untyped =
  gc.gradientFlow(g,0,eps):
    let flowTime {.inject.} = wflowT
    measure

proc flowMeasurements(u: auto; loop: int; tau: float): JsonNode =
  var
    pls,plt: ComplexProxy[ComplexObj[float64,float64]]
    poly: seq[ComplexProxy[ComplexObj[float64,float64]]]
    t2Ess,t2Est,t2Ees,t2Eet: float
  let
    f = u.fmunu(loop)
    (es, et) = f.densityE
    q = f.topoQ
    pl = u.plaq
    nl = pl.len div 2
    ss = 6.0*pl[0..<nl].sum
    st = 6.0*pl[nl..^1].sum
    pg = u[0].l.physGeom
  poly = newSeq[ComplexProxy[ComplexObj[float64,float64]]](pg.len)
  for i in 0..<pg.len: poly[i] = u.wline repeat(i+1,pg[i])
  pls = poly[0..^2].sum/float(poly.len-1)
  plt = poly[^1]
  (t2Ess,t2Est) = (6.0*tau*tau*(3.0-ss),6.0*tau*tau*(3.0-st))
  (t2Ees,t2Eet) = (tau*tau*es,tau*tau*et)
  result = %* {
    "flow-time": tau,
    "plaquette":0.5*ss+0.5*st,
    "clover":es+et,
    "t2E-plaquette":t2Ess+t2Est,
    "t2E-clover":t2Ees+t2Eet,
    "t2E-spacelike-plaquette":t2Ess,
    "t2E-timelike-plaquette":t2Est,
    "t2E-spacelike-clover":t2Ees,
    "t2E-timelike-clover":t2Eet,
    "topological-charge":q,
    "Re(spacelike-Polyakov-loop)":3.0*pls.re,
    "Im(spacelike-Polyakov-loop)":3.0*pls.im,
    "Re(timelike-Polyakov-loop)":3.0*plt.re,
    "Im(timelike-Polyakov-loop)":3.0*plt.im,
  }

proc get(info:JsonNode;key:string): seq[float] =
  result = newSeq[float]()
  for el in info[key].getElems(): result.add getFloat(el)

template gradientFlow*(u: auto; info: JsonNode; body: untyped) =
  var 
    v = u[0].l.newGauge
    f {.inject,used.}: File
    tau {.inject.}: float
    measurements {.inject.}: JsonNode
  for flow,flowInfo in info:
    let
      loops = case info[flow].hasKey("loops")
        of true: info[flow]["loops"].getInt()
        of false: 1
      beta = case info[flow].hasKey("beta")
        of true: info[flow]["beta"].getFloat()
        of false: 1.0
      cr = case info[flow]["action"].getStr()
        of "Rectangle": info[flow]["cr"].getFloat()
        else: 0.0
      gc = case info[flow]["action"].getStr()
        of "Rectangle": gaugeActRect(beta,cr)
        else: GaugeActionCoeffs(plaq:beta)
      fn = info[flow]["path"].getStr() & info[flow]["filename"].getStr()
      dts = info[flow].get("time-increments")
      maxFlts = info[flow].get("maximum-flow-times")
    var lastMaxFlt = 0.0
    threads: v := u
    v.reunit
    f = fn.open(fmWrite)
    for (dt,maxFlt) in zip(dts,maxFlts):
      gc.gradientFlow(v,dt):
        tau = flowTime + lastMaxFlt
        measurements = v.flowMeasurements(loops,tau)
        body
        if flowTime>=maxFlt: break
      lastMaxFlt = maxFlt
    f.close()
      
proc formatMeasurements*(
    measurements: JsonNode;
    style: string = "default"
  ): string =
  result = ""
  let tau = measurements["flow-time"].getFloat()
  case style:
    of "default":
      var concatenation = @["FLOW",tau.formatFloat(ffDecimal,3)]
      result = concatenation.join(" ")
      for key,measurement in measurements:
        let O = measurement.getFloat()
        concatenation = @[result,O.formatFloat(ffDecimal,13)] 
        result = concatenation.join(" ")
    of "KS_nHYP_FA": # Style of https://github.com/daschaich/KS_nHYP_FA
      let
        standIn = %* {"stand-in": 0.0}
        observables = @[
          measurements["plaquette"],
          measurements["clover"],
          measurements["t2E-clover"],
          standIn["stand-in"],
          measurements["t2E-plaquette"],
          measurements["topological-charge"],
          measurements["t2E-spacelike-clover"],
          measurements["t2E-timelike-clover"],
          measurements["Re(timelike-Polyakov-loop)"],
          measurements["Im(timelike-Polyakov-loop)"],
          measurements["Re(spacelike-Polyakov-loop)"],
          measurements["Im(spacelike-Polyakov-loop)"]
        ]
      var concatenation = @["FLOW",tau.formatFloat(ffDecimal,2)]
      result = concatenation.join(" ")
      for observable in observables:
        let O = observable.getFloat()
        concatenation = @[result,O.formatFloat(ffDecimal,13)] 
        result = concatenation.join(" ")
    else: discard

proc setFilename*(info: auto; fn: string) =
  info["filename"] = %* fn 