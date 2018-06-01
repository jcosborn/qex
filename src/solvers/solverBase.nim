import base, layout

type
  SolverParams* = object
    r2req*:float
    maxits*:int
    verbosity*:int
    finalIterations*:int
    seconds*: float
    subset*:Subset
    subsetName*:string

proc initSolverParams*():SolverParams =
  result.r2req = 1e-6
  result.maxits = 2000
  result.verbosity = 1
  result.subsetName = "all"
