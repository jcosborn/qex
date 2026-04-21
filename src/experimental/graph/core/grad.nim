import base
import grad_cache, grad_engine, cond

export
  GradCacheStats, clearGradCache, resetGradCacheStats,
  resetGradCache, dumpGradientList,
  gradOrZero,
  cond

proc grad*(dep: Gvalue, x: Gvalue): Gvalue =
  grad_engine.grad(dep, x)

proc gradIsolated*(dep: Gvalue, x: Gvalue): Gvalue =
  grad_engine.gradIsolated(dep, x)
