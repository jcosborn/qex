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

