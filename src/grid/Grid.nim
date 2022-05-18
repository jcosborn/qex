template haveGrid*():bool = defined(gridDir)

when haveGrid():

  import grid/gridImpl
  export gridImpl

else:  # put stubs here

  import base
  import layout
  import physics/stagD

  proc gridSolveEE*(s:Staggered; r,t:Field; m:SomeNumber; sp: var SolverParams) =
    qexError "Grid not compiled in (define gridDir)"
