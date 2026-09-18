import gauge/types, gauge/basic_ops, gauge/matfun, gauge/fused_ops
import gauge/matrix
import gauge/stencil
import gauge/action/ops
import gauge/stout

export matrix

export
  Gauge, Ggauge,
  gaugeSnapshot,
  update, mutateGauge, toGvalue,
  retr, adj, norm2, redot, exp, expDeriv, projTAH,
  `-`, `+`, `*`,
  axpy, adjmul, muladj, contractProjTAH, axexp, axexpmuly,
  Gactcoeff,
  actWilson, actSymanzik, actIwasaki, actDBW2, actAdj,
  gaugeAction, gaugeActionDeriv, gaugeActionDeriv2, gaugeForce,
  gaugeActionGraph, adjPlaqAction,
  plaqSum, stapleSum,
  blendSubset, stoutUpdate, stoutLogDetJ, stoutLogDetJGraph, stoutUpdateLogDetJ
