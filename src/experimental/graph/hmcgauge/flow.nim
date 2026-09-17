import ../[core, scalar, gauge]
import integrator
import std/tables

type
  GaugeFlow* = proc(V: Ggauge): Ggauge {.closure.}
  FlowAction* = object
    action*: GaugeAction
    flow*: GaugeFlow

proc flowAction*(gc: Gactcoeff, map: GaugeFlow): FlowAction =
  ## S_eff(V) = S(f(V)) - log det f'(V).
  ## Build map(V) once per input node so action, force and flow consumers share
  ## f(V): plans clone their graphs at construction, so nothing downstream can
  ## merge two equivalent chains, and logDetJ reuses its cached chain sum.
  ## Either returned closure retains every cached flow for its lifetime.
  var flows = initTable[NodeKey, Ggauge]()
  proc getFlow(V: Ggauge): Ggauge =
    let key = V.nodeKey
    if key notin flows:
      flows[key] = map(V)
    flows[key]
  result.flow = getFlow
  result.action = proc(V: Ggauge): Gscalar =
    let u = getFlow(V)
    gaugeAction(gc, u) - logDetJ(u, V)
