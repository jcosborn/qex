import base, layout

type
  SolverParams* = object
    r2req*:float
    maxits*:int
    sloppySolve*:SloppyType
    verbosity*:int
    finalIterations*:int
    seconds*: float
    subset*:Subset
    subsetName*:string
  SloppyType* = enum
    SloppyNone, SloppySingle, SloppyHalf

proc initSolverParams*():SolverParams =
  result.r2req = 1e-6
  result.maxits = 50000
  result.sloppySolve = SloppyNone
  result.verbosity = 1
  result.subsetName = "all"
