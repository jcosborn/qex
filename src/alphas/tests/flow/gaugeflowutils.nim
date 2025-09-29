import qex
import gauge
import gauge/[gaugeUtils]
import json,sequtils,strutils

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

proc flowMeasurements*(u: auto; loop: int; tau: float): JsonNode =
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
    #v.reunit
    f = fn.open(fmWrite)
    let m = u.flowMeasurements(1, 0.0)
    f.write(m.formatMeasurements(style = logStyle) & "\n")
    for flowTimeInfo in zip(dts,maxFlts):
      let (dtau,tmax) = flowTimeInfo
      gc.gradientFlow(v,dtau):
        tau = flowTime + lastMaxFlt
        if tau > tmax: break
        measurements = v.flowMeasurements(loops,tau)
        body
      lastMaxFlt = tmax
    f.close()

# old
# FLOW 4.980 0.9905002120655 NORECT 0.9082318792468 1.1954405117255
# FLOW 5.000 0.9905658406963 NORECT 0.9087972026386 1.2217308887443
# FLOW 5.100 0.9909487862018 NORECT 0.9120949058955 1.3686707539197

# new
# FLOW 4.980 0.990500212065502095 0.966141443113515686 0.908231879246808593 1.195440511725532673 0.006748192248644567 0.004765314947830970 0.032650452152740060 0.049245689401411948
# FLOW 5.000 0.990565840696340905 0.966361698075400621 0.908797202638550905 1.221730888744295473 0.006743259215785240 0.004550885130445288 0.032005542691391359 0.049815500426037426
# FLOW 5.100 0.990886226394149139 0.967436850021683115 0.911556254425717860 1.345716794336233368
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
          standIn["stand-in"], # i.e., nothing; where rectangle could be
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
    of "regress":
      var concatenation = @["FLOW",tau.formatFloat(ffDecimal,3)]
      concatenation.add (measurements["plaquette"].getFloat()/3.0).formatFloat(ffDecimal,13)
      concatenation.add "NORECT" # no rectangle measurement
      concatenation.add ((3.0 - measurements["clover"].getFloat())/3.0).formatFloat(ffDecimal,13)
      concatenation.add measurements["topological-charge"].getFloat().formatFloat(ffDecimal,13)
      concatenation.add measurements["Re(spacelike-Polyakov-loop)"].getFloat().formatFloat(ffDecimal,13)
      concatenation.add measurements["Im(spacelike-Polyakov-loop)"].getFloat().formatFloat(ffDecimal,13)
      concatenation.add measurements["Re(timelike-Polyakov-loop)"].getFloat().formatFloat(ffDecimal,13)
      concatenation.add measurements["Im(timelike-Polyakov-loop)"].getFloat().formatFloat(ffDecimal,13)
      result = concatenation.join(" ")
    else: discard