import macros
import base/metaUtils
import expr

proc addChildrenFrom*(dst,src: NimNode): NimNode =
  for c in src: dst.add(c)
  result = dst
macro procInst*(p: typed): auto =
  #echo "begin procInst:"
  #echo p.treerepr
  result = p[0]
macro makeCall*(p: proc, x: tuple): NimNode =
  result = newCall(p).addChildrenFrom(x)

type
  CudaDim3* {.importc:"dim3",header:"cuda_runtime.h".} = object
    x*, y*, z*: cint
  cudaError_t* {.importc,header:"cuda_runtime.h".} = object
  cudaMemcpyKind* {.importc,header:"cuda_runtime.h".} = object
var
  cudaSuccess*{.importC,header:"cuda_runtime.h".}: cudaError_t
  cudaErrorNotSupported*{.importC,header:"cuda_runtime.h".}: cudaError_t
  cudaMemcpyHostToDevice*{.importC,header:"cuda_runtime.h".}: cudaMemcpyKind
  cudaMemcpyDeviceToHost*{.importC,header:"cuda_runtime.h".}: cudaMemcpyKind

#template toPointer*(x: pointer): pointer = x
#template toPointer*[T](x: ptr T): pointer = pointer(x)
#template toPointer*(x: seq): pointer = toPointer(x[0])
#template toPointer*(x: not (pointer|seq)): pointer = pointer(unsafeAddr(x))
template toPointer*(x: typed): pointer =
  #dumpType: x
  when x is pointer: x
  elif x is ptr: x
  elif x is seq: toPointer(x[0])
  else: pointer(unsafeAddr(x))
template dataAddr*(x: typed): pointer =
  #dumpType: x
  when x is seq: dataAddr(x[0])
  elif x is array: dataAddr(x[0])
  #elif x is ptr: x
  else: pointer(unsafeAddr(x))
  #else: x

proc cudaGetLastError*(): cudaError_t
  {.importC,header:"cuda_runtime.h".}
proc cudaGetErrorStringX*(error: cudaError_t): ptr char
  {.importC:"cudaGetErrorString",header:"cuda_runtime.h".}
proc cudaGetErrorString*(error: cudaError_t): cstring =
  var s {.codegendecl:"const $# $#".} = cudaGetErrorStringX(error)
  result = s
proc `$`*(error: cudaError_t): string =
  let s = cudaGetErrorString(error)
  result = $s
converter toBool*(e: cudaError_t): bool =
  cast[cint](e) != cast[cint](cudaSuccess)

proc cudaMalloc*(p:ptr pointer, size: csize): cudaError_t
  {.importC,header:"cuda_runtime.h".}
template cudaMalloc*(p:pointer, size: csize): cudaError_t =
  cudaMalloc((ptr pointer)(p.addr), size)
proc cudaFree*(p: pointer): cudaError_t
  {.importC,header:"cuda_runtime.h".}
proc cudaMallocManaged*(p: ptr pointer, size: csize): cudaError_t
  {.importC,header:"cuda_runtime.h".}

proc cudaMemcpyX*(dst,src: pointer, count: csize, kind: cudaMemcpyKind):
  cudaError_t {.importC:"cudaMemcpy",header:"cuda_runtime.h".}
template cudaMemcpy*(dst,src: typed, count: csize,
                     kind: cudaMemcpyKind): cudaError_t =
  let pdst = toPointer(dst)
  let psrc = toPointer(src)
  cudaMemcpyX(pdst, psrc, count, kind)

proc cudaLaunchKernel(p:pointer, gd,bd: CudaDim3, args: ptr pointer):
  cudaError_t {.importC,header:"cuda_runtime.h".}

proc cudaDeviceReset*(): cudaError_t
  {.importC,header:"cuda_runtime.h".}
proc cudaDeviceSynchronize*(): cudaError_t
  {.importC,header:"cuda_runtime.h".}

#proc printf*(fmt:cstring):cint {.importc,varargs,header:"<stdio.h>",discardable.}
#proc fprintf*(stream:ptr FILE,fmt:cstring):cint {.importc,varargs,header:"<stdio.h>".}
#proc malloc*(size: csize):pointer {.importc,header:"<stdlib.h>".}

template cudaDefs(body: untyped): untyped {.dirty.} =
  var gridDim{.global,importC,noDecl.}: CudaDim3
  var blockIdx{.global,importC,noDecl.}: CudaDim3
  var blockDim{.global,importC,noDecl.}: CudaDim3
  var threadIdx{.global,importC,noDecl.}: CudaDim3
  template getGridDim: untyped {.used.} = gridDim
  template getBlockIdx: untyped {.used.} = blockIdx
  template getBlockDim: untyped {.used.} = blockDim
  template getThreadIdx: untyped {.used.} = threadIdx
  template getThreadNum: untyped {.used.} = blockDim.x * blockIdx.x + threadIdx.x
  template getNumThreads: untyped {.used.} = gridDim.x * blockDim.x
  bind inlineProcs
  inlineProcs:
    body

template cudaLaunch*(p: proc; blocksPerGrid,threadsPerBlock: SomeInteger;
                     arg: varargs[pointer,dataAddr]) =
  var pp: proc = p
  var gridDim, blockDim: CudaDim3
  gridDim.x = blocksPerGrid
  gridDim.y = 1
  gridDim.z = 1
  blockDim.x = threadsPerBlock
  blockDim.y = 1
  blockDim.z = 1
  var args: array[arg.len, pointer]
  for i in 0..<arg.len: args[i] = arg[i]
  #echo "really launching kernel"
  let err = cudaLaunchKernel(pp, gridDim, blockDim, addr args[0])
  if err:
    echo err
    quit cast[cint](err)

template `<<`*(p: proc, x: tuple): untyped = (p,x)
template getInst*(p: untyped): untyped =
  #when compiles((var t=p; t)): p
  #else:
  procInst(p)
    #var t =
    #t
macro `>>`*(px: tuple, y: any): auto =
  #echo "begin >>:"
  #echo px.treerepr
  #echo "kernel type:"
  #echo px[0].getTypeImpl.treerepr
  #echo "kernel args:"
  #echo y.treerepr
  #var a = y
  #if y.kind != nnkPar: a = newNimNode(nnkPar).addChildrenFrom(y)
  result = newCall(ident("cudaLaunch"))
  let krnl = newCall(px[0]).addChildrenFrom(y)
  #echo "kernel inst call:"
  #echo krnl.treerepr
  result.add getAst(getInst(krnl))[0]
  result.add px[1][0]
  result.add px[1][1]
  for c in y:
    result.add c
  #echo "kernel launch body:"
  #echo result.treerepr

proc cudaproc(s:string, p:NimNode):NimNode =
  #echo "begin cuda:"
  #echo s
  #let ss = s.strVal
  #echo "proc:"
  #echo p.treerepr
  p.expectKind nnkProcDef
  result = p
  # if p.kind == nnkProcDef:
  #   result = p
  # else:
  #   result = p[0]
  result.addPragma parseExpr("{.codegenDecl:\""&s&" $# $#$#\".}")[0]
  result.body = getAst(cudaDefs(result.body))
  var sl = newStmtList()
  sl.add( quote do:
    {.push checks: off.}
    {.push stacktrace: off.} )
  sl.add result
  result = sl
  #echo "end cuda:"
  #echo result.treerepr
macro cudaGlobal*(p: untyped): untyped = cudaproc("__global__",p)

template onGpu*(nn,tpb: untyped, body: untyped): untyped =
  block:
    var v = packVars(body, getGpuPtr)
    type ByCopy[T] {.bycopy.} = object
      d: T
    proc kern(xx: ByCopy[type(v)]) {.cudaGlobal.} =
      template deref(k: int): untyped = xx.d[k]
      substVars(body, deref)
    let ni = nn.int32
    let threadsPerBlock = tpb.int32
    let blocksPerGrid = (ni+threadsPerBlock-1) div threadsPerBlock
    #echo "launching kernel"
    cudaLaunch(kern, blocksPerGrid, threadsPerBlock, v)
    discard cudaDeviceSynchronize()
template onGpu*(nn: untyped, body: untyped): untyped = onGpu(nn, 64, body)
template onGpu*(body: untyped): untyped = onGpu(512*64, 64, body)

when isMainModule:
  type FltArr = UncheckedArray[float32]

  proc vectorAdd(A: FltArr; B: FltArr; C: var FltArr; n: int32) {.cudaGlobal.} =
    var i = blockDim.x * blockIdx.x + threadIdx.x
    if i < n:
      C[i] = A[i] + B[i]

  proc test =
    var n = 50000.cint
    var
      a = newSeq[float32](n)
      b = newSeq[float32](n)
      c = newSeq[float32](n)
    var threadsPerBlock: cint = 256
    var blocksPerGrid: cint = (n + threadsPerBlock - 1) div threadsPerBlock

    cudaLaunch(vectorAdd, blocksPerGrid, threadsPerBlock, a, b, c, n)

    template getGpuPtr(x: int): untyped = x
    template getGpuPtr[T](x: seq[T]): untyped = addr(x[0])
    template `[]`(x: ptr SomeNumber, i: SomeInteger): untyped {.used.} =
      cast[ptr UncheckedArray[type(x[])]](x)[][i]
    template `[]=`(x: ptr SomeNumber, i: SomeInteger, y:untyped): untyped {.used.} =
      cast[ptr UncheckedArray[type(x[])]](x)[][i] = y

    onGpu(n):
      let i = getBlockDim().x * getBlockIdx().x + getThreadIdx().x
      if i < n:
        c[i] = a[i] + b[i]

  test()
