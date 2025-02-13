include system/ansi_c

{. pragma: execinfo, header:"execinfo.h" .}

proc backtrace(buffer: ptr pointer, size: cint): cint {.importc,execinfo.}
proc backtrace_symbols(buffer: ptr pointer, size: cint): ptr UncheckedArray[cstring]
    {.importc,execinfo.}

proc print_trace =
  #void *array[10];
  const nmax = 20
  var arr: array[nmax, pointer]
  #char **strings;

  let size = backtrace(addr arr[0], nmax)
  let strings = backtrace_symbols(addr arr[0], size);
  if not isNil strings:
    echo "Stack trace of size: ", size
    for i in 0..<size:
      echo strings[i]
  cfree strings

proc sigtrace(sig: cint) {.noconv.} =
  print_trace()
  quit()

proc setTrace* =
  c_signal(SIGSEGV, sigtrace)

when isMainModule:
  print_trace()
  setTrace()
  discard c_raise(SIGSEGV)
