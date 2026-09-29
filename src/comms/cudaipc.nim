## CUDA IPC: device memory of one process mapped into another on the same
## node, so kernels can store straight into the memory of a peer, as
## comms/zeipc does with Level Zero.  ipcExport describes device memory of
## this process, a peer maps it with ipcOpen and unmaps it with ipcClose.
import qex

{.passL: "-lcuda".}
{.emit: "/*INCLUDESECTION*/\n#include <cuda.h>".}

type
  IpcMemHandle {.importc: "cudaIpcMemHandle_t", header: "cuda_runtime.h", completeStruct, bycopy.} = object
    reserved: array[64, char]
  GpuIpc* = object
    ## the handle of the allocation holding the memory and the offset of the
    ## memory in the allocation
    off: int
    h: IpcMemHandle

proc cudaIpcGetMemHandle(h: ptr IpcMemHandle; p: pointer): cint {.importc, header: "cuda_runtime.h".}
proc cudaIpcOpenMemHandle(p: ptr pointer; h: IpcMemHandle; flags: cuint): cint {.importc, header: "cuda_runtime.h".}
proc cudaIpcCloseMemHandle(p: pointer): cint {.importc, header: "cuda_runtime.h".}

proc ipcExport*(p: pointer): GpuIpc =
  ## device memory at p for a peer; the handle covers the whole allocation
  var base: pointer
  var err: cint
  {.emit: ["{ CUdeviceptr b; size_t s; ", err, " = cuMemGetAddressRange(&b, &s, (CUdeviceptr)", p, "); ", base, " = (void*)b; }"].}
  if err != 0 or cudaIpcGetMemHandle(addr result.h, base) != 0:
    qexError("cudaIpcGetMemHandle failed")
  result.off = cast[int](p) - cast[int](base)

proc ipcOpen*(e: GpuIpc): tuple[base, p: pointer] =
  ## maps the memory e of a peer: p is its address here, base goes to ipcClose
  if cudaIpcOpenMemHandle(addr result.base, e.h, 1) != 0:  # cudaIpcMemLazyEnablePeerAccess
    qexError("cudaIpcOpenMemHandle failed")
  result.p = cast[pointer](cast[int](result.base) + e.off)

proc ipcClose*(base: pointer) =
  discard cudaIpcCloseMemHandle(base)
