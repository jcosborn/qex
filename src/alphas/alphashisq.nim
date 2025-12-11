when defined(gridDir):
  import extra/[gridsmear]#, gridgauge]
  import extra/[gridgauge]
  export gridsmear
  export gridgauge
  #import extra/[qexgauge]
  #export qexgauge
else:
  import extra/[qexsmear, qexgauge]
  export qexsmear
  export qexgauge

template gridBackend*(work: untyped): untyped = 
  when defined(gridDir): work

template qexBackend*(work: untyped): untyped = 
  when not defined(gridDir): work

proc newHISQ*(
    lepage: float = 0.0; 
    naik: float = 1.0,
    reunitMethod: ProjectionMethod = CayleyHamilton,
    reunitEps: float = 1e-16,
    reunitMaxiters: int = 10,
    delta: float = 1e-5
  ): HisqCoefs =
  result = HisqCoefs(naik: -naik/24.0)
  result.fat7first.setHisqFat7(lepage,0.0)
  result.fat7second.setHisqFat7(2.0-lepage,naik)
  result.projection = newUnitaryProjection(reunitMethod, reunitEps, reunitMaxiters)
  case reunitMethod:
    of CayleyHamilton: result.projection.delta = delta
    else: discard