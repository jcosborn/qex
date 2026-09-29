import os, macros

when defined(noOpenmp):
  static: echo "OpenMP disabled"
  template omp_set_num_threads*(x: cint) = discard
  template omp_get_num_threads*(): cint = cint 1
  template omp_get_max_threads*(): cint = cint 1
  template omp_get_thread_num*(): cint = cint 0
  template ompPragma(p:string) = discard
  template ompPragmaNv(nv, other: string) = discard
  template ompBlock*(p:string; body:untyped) =
    block:
      body
else:
  static: echo "Using OpenMP"
  when existsEnv("OMPFLAG"):
    const ompFlag {.strDefine.} = getEnv("OMPFLAG")
  else:
    const ompFlag {.strDefine.} = "-fopenmp"
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
  template ompPragma(p:string,body:typed) =
    {. push stackTrace:off, lineTrace:off, line_dir:off .}
    {. emit:["_Pragma(\"omp ", p, "\")"] .}
    body
    {. pop .}
  template ompPragmaNv(nv, other: string) =
    ## the pragmas nv for NVHPC, other for the other compilers
    {. emit:["\n#ifdef __NVCOMPILER\n", nv, "\n#else\n", other, "\n#endif\n"] .}
  template ompBlock*(p:string; body:untyped) =
    #{. emit:"#pragma omp " & p .}
    #{. emit:"{ /* Inserted by ompBlock " & p & " */".}
    #{. emit:["#pragma omp ", p] .}
    ompPragma(p)
    block:
      body
    #{. emit:"} /* End ompBlock " & p & " */".}

  macro ompPragma2*(p: varargs[untyped]): auto =
    var b = newNimNode(nnkBracket)
    b.add newLit "_Pragma(\"omp "
    for i in 0..<p.len:
      if p[i].kind == nnkTupleConstr:
        for j in 0..<p[i].len:
          b.add p[i][j]
      else:
        b.add p[i]
    b.add newLit "\")"
    #echo p.treerepr
    result = quote do:
      {. emit:`b` .}
  #template ompBlock2*(p: varargs[untyped]) =
  #  ompPragma2(p[0..^2])
  #  block:
  #    p[^1]
  macro ompBlock2*(p: varargs[untyped]): auto =
    {. push stackTrace:off, lineTrace:off, line_dir:off .}
    #echo p.treerepr
    let body = p[^1]
    var p2 = newNimNode(nnkCall).add bindSym"ompPragma2"
    for i in 0..(p.len-2):
      p2.add p[i]
    #echo body.treerepr
    result = quote do:
      `p2`
      block:
        `body`
    #echo result.treerepr
    {. pop .}

template ompBarrier* = ompPragma("barrier")
template ompFlush* = ompPragma("flush")
# NVHPC (nvc, nvc++ 25.11) rejects the memory order clauses of OpenMP 5.0:
# there a flush without clauses, a full fence, orders the atomic accesses
template ompFlushAcquire* = ompPragmaNv("_Pragma(\"omp flush\")", "_Pragma(\"omp flush acquire\")")
template ompFlushRelease* = ompPragmaNv("_Pragma(\"omp flush\")", "_Pragma(\"omp flush release\")")
template ompFlushSeqCst* = ompPragmaNv("_Pragma(\"omp flush\")", "_Pragma(\"omp flush seq_cst\")")
template ompAtomicRead*(body) =
  {. push stackTrace:off, lineTrace:off, line_dir:off .}
  ompPragmaNv("_Pragma(\"omp atomic read\")", "_Pragma(\"omp atomic read acquire\")")
  body
  ompPragmaNv("_Pragma(\"omp flush\")", "")
  {. pop .}
template ompAtomicWrite*(body) =
  {. push stackTrace:off, lineTrace:off, line_dir:off .}
  ompPragmaNv("_Pragma(\"omp flush\")\n_Pragma(\"omp atomic write\")", "_Pragma(\"omp atomic write release\")")
  body
  {. pop .}

template ompParallel*(body:untyped) =
  ompBlock("parallel"):
    when(declared(setupForeignThreadGc)):
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
