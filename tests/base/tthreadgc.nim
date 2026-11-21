import testutils
import qex

# refc: a worker thread of a threads block collects its own heap, and what
# only the locals of the block refer to must stay.  tthreadgc.nims gives a
# collection at every allocation and freed memory filled with 0xFF.
#
# ompParallel moves the stack bottom of a worker thread to the end of its
# stack (gcWholeStack in base/omp.nim).  Without that call the bottom stays
# where Nim puts it, in the frame of setupForeignThreadGc below the first
# block the thread ran: the collector then misses the locals of that block,
# and of a later block those above that frame, which depends on the
# compiler.  So the block runs first in the program, before those of
# qexInit, and once more after them.

proc changed(): seq[int] =
  ## Per thread: of the values of 24 seqs kept in a local array of a threads
  ## block, those that changed while the thread allocated 100 more seqs.
  var c = newSeq[int](getMaxThreads())
  let pc = cast[ptr UncheckedArray[int]](addr c[0])
  threads:
    var a: array[24, seq[int]]
    for k in 0..<a.len:
      a[k] = newSeq[int](8)
      for j in 0..<8: a[k][j] = 1000*threadNum + 10*k + j
    var b = 0
    for rep in 0..<100:
      let t = newSeq[int](16)
      b += t.len
    var nb = 0
    for k in 0..<a.len:
      let p = cast[ptr UncheckedArray[int]](addr a[k][0])  # past the length of a freed seq
      for j in 0..<8:
        if p[j] != 1000*threadNum + 10*k + j: inc nb
    pc[threadNum] = nb + 0*b
  c

let c = changed()
qexInit()

suite "Collections on the threads of a threads block":
  test "what only locals of the first block of the program refer to is kept":
    check c == newSeq[int](c.len)

  test "what only locals of a later block refer to is kept":
    let d = changed()
    check d == newSeq[int](d.len)

qexFinalize()
