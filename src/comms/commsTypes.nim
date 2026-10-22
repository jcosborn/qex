type
  Comm* = ref object of RootObj

# globals

var defaultComm*: Comm
template getDefaultComm*(): Comm = defaultComm
template getComm*(): Comm = getDefaultComm()  # temporary alias
var myRank* = 0
var nRanks* = 1

when defined(qexBorrowCheck) and not defined(gcDestructors):
  # Debug check for memory handed to a transfer that uses it after the call
  # returns or across allocations (async sends/recvs, scatter, pointer gathers).
  # Under refc a GC object with refcount 0 stays alive only while some stack
  # slot holds its address, which the C compiler may drop before the transfer
  # ends (makeGatherMap in comms/gather).  The check asks for a counted owner
  # even where a later use keeps the object on the stack.  -d:nimFulldebug
  # exports the lookup.
  import system/ansi_c
  when not declared(interiorAllocatedPtr):
    {.error: "-d:qexBorrowCheck needs -d:nimFulldebug".}
  proc checkBorrow*(p: pointer, what: cstring) =
    let c = interiorAllocatedPtr(p)
    # refc Cell: refcount is the first word, the count starts at bit 3
    if c != nil and (cast[ptr int](c)[] shr 3) == 0:
      discard c_fprintf(cstderr, "[BORROW] %s: %p lies in a GC object with refcount 0; keep its owner referenced (GC_ref) until the transfer ends\n", what, p)
      c_abort()
else:
  template checkBorrow*(p: pointer, what: cstring) = discard
