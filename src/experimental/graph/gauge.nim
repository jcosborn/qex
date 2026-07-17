import gauge/shared, gauge/basic_ops, gauge/fused_ops
import gauge/action/ops
import gauge/stout

export
  Gauge, Ggauge,
  reunitGauge, gaugeSnapshot,
  update, mutateGauge, toGvalue,
  retr, adj, norm2, redot, exp, expDeriv, projTAH,
  `-`, `+`, `*`,
  axpy, adjmul, muladj, contractProjTAH, axexp, axexpmuly,
  Gactcoeff,
  actWilson, actSymanzik, actIwasaki, actDBW2, actAdj,
  gaugeAction, gaugeActionDeriv, gaugeActionDeriv2, gaugeForce,
  blendSubset, stoutUpdate, stoutLogDetJ, stoutUpdateLogDetJ
