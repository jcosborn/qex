template haveQuda*():bool = defined(qudaDir)

when haveQuda():

  import qudaWrapperImpl
  export qudaWrapperImpl

else:  # put stubs here

  import base
  import layout
  import physics/stagD
  import gauge

  proc qudaSolveEE*(s:Staggered; r,t:Field; m:SomeNumber; sp: var SolverParams) =
    qexError "QUDA not compiled in (define qudaDir)"

  proc qudaGaugeForce*(c: GaugeActionCoeffs, g: openArray, f: openArray) =
    qexError "QUDA not compiled in (define qudaDir)"
