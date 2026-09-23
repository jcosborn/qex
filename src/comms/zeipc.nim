## Level Zero IPC: device memory of one process mapped into another on the
## same node, as on Aurora and Sunspot, so kernels can store straight into
## the memory of a peer.  zeExport describes device memory of this process,
## a peer maps it with zeOpen and unmaps it with zeClose.
import qex
import backend/accel

{.passL: "-lze_loader".}

type
  Ctx {.importc: "ze_context_handle_t", header: "level_zero/ze_api.h".} = pointer
  Dev {.importc: "ze_device_handle_t", header: "level_zero/ze_api.h".} = pointer
  ZeIpcMemHandle {.importc: "ze_ipc_mem_handle_t", header: "level_zero/ze_api.h".} = object
    data: array[64, char]
  ZeIpc* = object
    ## the handle of the allocation holding the memory, with the file
    ## descriptor of the owner process pid in its first bytes, and the
    ## offset of the memory in the allocation
    pid, fd, off: int
    h: ZeIpcMemHandle

proc zeMemGetIpcHandle(ctx: Ctx; p: pointer; h: var ZeIpcMemHandle): cint {.importc, header: "level_zero/ze_api.h".}
proc zeMemOpenIpcHandle(ctx: Ctx; dev: Dev; h: ZeIpcMemHandle; flags: uint32; p: var pointer): cint {.importc, header: "level_zero/ze_api.h".}
proc zeMemCloseIpcHandle(ctx: Ctx; p: pointer): cint {.importc, header: "level_zero/ze_api.h".}
proc zeMemGetAddressRange(ctx: Ctx; p: pointer; base: var pointer; bytes: var csize_t): cint {.importc, header: "level_zero/ze_api.h".}
proc zeHandles(): (Ctx, Dev) =
  let (c, d) = gpuZeContext()
  (cast[Ctx](c), cast[Dev](d))
proc syscall(n: clong): clong {.importc, header: "<unistd.h>", varargs.}
proc getpid(): cint {.importc, header: "<unistd.h>".}
proc close(fd: cint): cint {.importc, header: "<unistd.h>".}
var SYS_pidfd_open {.importc, header: "<sys/syscall.h>".}: clong
var SYS_pidfd_getfd {.importc, header: "<sys/syscall.h>".}: clong

proc zeExport*(p: pointer): ZeIpc =
  ## device memory at p for a peer; the handle covers the whole allocation,
  ## e.g. a block of the memory pool of the runtime
  let (ctx, _) = zeHandles()
  var base: pointer
  var bytes: csize_t
  if zeMemGetAddressRange(ctx, p, base, bytes) != 0 or zeMemGetIpcHandle(ctx, base, result.h) != 0:
    qexError("zeMemGetIpcHandle failed")
  result.pid = getpid()
  result.off = cast[int](p) - cast[int](base)
  copyMem(addr result.fd, addr result.h.data[0], sizeof(cint))

proc zeOpen*(e: ZeIpc): tuple[base, p: pointer] =
  ## maps the memory e of a peer: p is its address here, base goes to zeClose
  let (ctx, dev) = zeHandles()
  let pfd = cint syscall(SYS_pidfd_open, e.pid, 0)
  let fd = cint syscall(SYS_pidfd_getfd, pfd, e.fd, 0)
  if pfd < 0 or fd < 0: qexError("pidfd_getfd failed")
  discard close(pfd)
  var h = e.h
  copyMem(addr h.data[0], unsafeAddr fd, sizeof(cint))
  if zeMemOpenIpcHandle(ctx, dev, h, 0, result.base) != 0:
    qexError("zeMemOpenIpcHandle failed")
  result.p = cast[pointer](cast[int](result.base) + e.off)

proc zeClose*(base: pointer) =
  let (ctx, _) = zeHandles()
  discard zeMemCloseIpcHandle(ctx, base)
