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

  proc qudaSolveOO*(s:Staggered; r,t:Field; m:SomeNumber; sp: var SolverParams) =
    qexError "QUDA not compiled in (define qudaDir)"

  proc qudaGaugeForce*[G,F](c: GaugeActionCoeffs, g: openArray[G], f: openArray[F]) =
    qexError "QUDA not compiled in (define qudaDir)"
