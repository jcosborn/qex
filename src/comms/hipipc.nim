## HIP IPC: device memory of one process mapped into another on the same
## node, so kernels can store straight into the memory of a peer, as
## comms/cudaipc does with CUDA.  ipcExport describes device memory of this
## process, a peer maps it with ipcOpen and unmaps it with ipcClose.
import qex

type
  IpcMemHandle {.importc: "hipIpcMemHandle_t", header: "hip/hip_runtime.h", completeStruct, bycopy.} = object
    reserved: array[64, char]
  GpuIpc* = object
    ## the handle of the allocation holding the memory and the offset of the
    ## memory in the allocation
    off: int
    h: IpcMemHandle

proc hipIpcGetMemHandle(h: ptr IpcMemHandle; p: pointer): cint {.importc, header: "hip/hip_runtime.h".}
proc hipIpcOpenMemHandle(p: ptr pointer; h: IpcMemHandle; flags: cuint): cint {.importc, header: "hip/hip_runtime.h".}
proc hipIpcCloseMemHandle(p: pointer): cint {.importc, header: "hip/hip_runtime.h".}
proc hipMemGetAddressRange(b: ptr pointer; s: ptr csize_t; p: pointer): cint {.importc, header: "hip/hip_runtime.h".}

proc ipcExport*(p: pointer): GpuIpc =
  ## device memory at p for a peer; the handle covers the whole allocation
  var base: pointer
  var bytes: csize_t
  if hipMemGetAddressRange(addr base, addr bytes, p) != 0 or hipIpcGetMemHandle(addr result.h, base) != 0:
    qexError("hipIpcGetMemHandle failed")
  result.off = cast[int](p) - cast[int](base)

proc ipcOpen*(e: GpuIpc): tuple[base, p: pointer] =
  ## maps the memory e of a peer: p is its address here, base goes to ipcClose
  if hipIpcOpenMemHandle(addr result.base, e.h, 1) != 0:  # hipIpcMemLazyEnablePeerAccess
    qexError("hipIpcOpenMemHandle failed")
  result.p = cast[pointer](cast[int](result.base) + e.off)

proc ipcClose*(base: pointer) =
  discard hipIpcCloseMemHandle(base)
