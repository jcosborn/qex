import gauge/shared, gauge/basic_ops, gauge/fused_ops
import gauge/action/ops

export
  Gauge, Ggauge,
  reunitGauge, gaugeSnapshot,
  update, mutateGauge, toGvalue,
  retr, adj, norm2, redot, exp, expDeriv, projTAH,
  `-`, `+`, `*`,
  adjmul, muladj, contractProjTAH, axexp, axexpmuly,
  Gactcoeff,
  actWilson, actSymanzik, actIwasaki, actDBW2, actAdj,
  gaugeAction, gaugeActionDeriv, gaugeActionDeriv2, gaugeForce
