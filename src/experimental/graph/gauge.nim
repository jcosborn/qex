import gauge/shared, gauge/basic_ops, gauge/fused_ops, gauge/action

export
  Gauge, Ggauge,
  requireGauge, requireSameGaugeShape,
  reunitGauge, gaugeSnapshot, zeroGaugeStorage,
  update, mutateGauge, toGvalue, gaugeNodeLike,
  valCopy, copyCompatible,
  retr, adj, norm2, redot, exp, expDeriv, projTAH,
  `-`, `+`, `*`,
  adjmul, muladj, contractProjTAH, axexp, axexpmuly,
  Gactcoeff, requireActCoeff,
  actWilson, actSymanzik, actIwasaki, actDBW2, actAdj,
  gaugeAction, gaugeActionDeriv, gaugeActionDeriv2, gaugeForce
