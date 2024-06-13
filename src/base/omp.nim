import os

when defined(noOpenmp):
  static: echo "OpenMP disabled"
  template omp_set_num_threads*(x: cint) = discard
  template omp_get_num_threads*(): cint = cint 1
  template omp_get_max_threads*(): cint = cint 1
  template omp_get_thread_num*(): cint = cint 0
  template ompPragma(p:string) = discard
  template ompBlock*(p:string; body:untyped) =
    block:
      body
else:
  static: echo "Using OpenMP"
  when existsEnv("OMPFLAG"):
    const ompFlag = getEnv("OMPFLAG")
  else:
    const ompFlag = "-fopenmp"
  {. passC: ompFlag .}
  {. passL: ompFlag .}
  {. pragma: omp, header:"omp.h" .}
  proc omp_set_num_threads*(x: cint) {.omp.}
  proc omp_get_num_threads*(): cint {.omp.}
  proc omp_get_max_threads*(): cint {.omp.}
  proc omp_get_thread_num*(): cint {.omp.}
  #proc forceOmpOn() {.omp.}
  template ompPragma(p:string) =
    #forceOmpOn()
    #{. emit:["#pragma omp ", p] .}
    {. emit:["_Pragma(\"omp ", p, "\")"] .}
  template ompBlock*(p:string; body:untyped) =
    #{. emit:"#pragma omp " & p .}
    #{. emit:"{ /* Inserted by ompBlock " & p & " */".}
    #{. emit:["#pragma omp ", p] .}
    ompPragma(p)
    block:
      body
    #{. emit:"} /* End ompBlock " & p & " */".}

template ompBarrier* = ompPragma("barrier")

template ompParallel*(body:untyped) =
  ompBlock("parallel"):
    if(omp_get_thread_num()!=0):
      setupForeignThreadGc()
    body
template ompMaster*(body:untyped) = ompBlock("master", body)
template ompSingle*(body:untyped) = ompBlock("single", body)
template ompCritical*(body:untyped) = ompBlock("critical", body)

when isMainModule:
  proc test =
    echo "main: ", ompGetThreadNum(), "/", ompGetNumThreads()
    ompParallel:
      echo "parallel: ", ompGetThreadNum(), "/", ompGetNumThreads()
      ompBarrier()
      ompMaster:
        echo "master: ", ompGetThreadNum(), "/", ompGetNumThreads()
      echo "parallel: ", ompGetThreadNum(), "/", ompGetNumThreads()
      ompSingle:
        echo "single: ", ompGetThreadNum(), "/", ompGetNumThreads()
      echo "parallel: ", ompGetThreadNum(), "/", ompGetNumThreads()
      ompCritical:
        echo "critical: ", ompGetThreadNum(), "/", ompGetNumThreads()
      echo "parallel: ", ompGetThreadNum(), "/", ompGetNumThreads()
      ompBarrier()
  test()
