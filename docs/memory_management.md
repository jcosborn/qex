# Memory management under ARC/ORC

QEX currently by default uses `--mm:refc` though that may change in a
future version.  Passing `--mm:` in `nimargs` selects another
mode. Code must work under `orc`, the other destructor modes (`arc`,
`atomicArc`, `yrc`), which all define `gcDestructors`, and `refc`. In
the destructor modes every copy of a `ref`, seq or string runs hooks
(`=copy`, `=sink`, `=destroy`), and under `arc` and `orc` the
reference counts are not atomic.

## Inside `threads` blocks

All threads of a region share the fields, layouts, communicators and shifters
the body uses. A reference count changed from several threads at once corrupts
the heap, so code in a region copies and destroys no shared reference.

```nim
threads:
  let n = x.l.nSitesOuter                # read through the owner
  let lo {.cursor.} = x.l                # or borrow it under a name
  let comm {.cursor.} = getDefaultComm()
  y := f ^* x                            # `^*` returns `lent F`
```

- Name a shared reference with `{.cursor.}`. Cursor inference does not cover
  globals, values captured by the threads proc, or variables whose fields are
  passed by address.
- Procs used in regions return borrowed values (`lent T`; the value must come
  from the first parameter), or are templates that set a `{.cursor.}` variable
  (the shift helpers of `gaugeProd` in `src/gauge/gaugeUtils.nim`).
- Objects that view other objects mark those fields `{.cursor.}` (the
  expression types in `src/field/fieldET.nim`). The viewed objects must outlive
  the view.
- Allocate outside the region, or in one thread between barriers (the lazy
  `{.global.}` field in `stagD2xx`).
- Seqs and strings local to a thread are fine. Assign a shared seq or string
  from one thread only.
- Do not `return` from a region or catch its exceptions outside it. A worker's
  exception prints the error and exits; the master's propagates after the block.

`{.cursor.}` is a raw pointer without runtime checks: the owner keeps the object
alive while the cursor is in use. Under refc, `{.cursor.}` and `lent` change
nothing.

## Everywhere

- `{.noInit.}` only on types without managed memory: numbers, SIMD values,
  plain arrays and objects. Under ARC/ORC the first assignment, `newSeq` or
  `setLen` of an uninitialized seq, string, ref, closure or object with a
  destructor reads the garbage it holds.
- Free resources in a `=destroy` hook on the object type under
  `when defined(gcDestructors)`; `new(x, finalizer)` is for refc only
  (`GatherBuffers` in `src/comms/commsTypes.nim`).
- Test `defined(gcDestructors)` for code that differs between the modes. Keep
  refc APIs behind `when not defined(gcDestructors)`: `shallowCopy`,
  `GC_ref`/`GC_unref` on seqs, `-d:nimFulldebug`, `interiorAllocatedPtr`.
- `p = f(p)` frees the old `p` before `f` reads it when `p`'s type is a ref
  type with a `static` generic parameter, such as `Field[V,T]`, `Layout[V]` or
  the graph's `Grmat[n]` (nim-lang/Nim#26680, seen in Nim 2.0, 2.2 and 2.4).
  Compute into a new name first: `let q = f(p); p = q` (`stoutLogDetJGraph` in
  `src/experimental/graph/gauge/stout.nim`).
- Leave a buffer of size 0 `nil` when it is later freed (`lbuf` in
  `prepareShiftBuf`, `src/layout/shiftX.nim`). Under ARC/ORC every `alloc(0)`
  returns the same block, so each owner would free the same pointer.
- ARC/ORC destroy global variables when the program ends, after `qexFinalize()`
  has finalized comms. Destructors and finalizers that call QMP or MPI check
  that the library is still initialized (`freeShiftBufsQ`).
- Raw pointers into managed memory (`unsafeAddr x[0]`, `x.data`): under refc the
  caller pins the owner with a `GC_ref`/`GC_unref` pair inside
  `when not defined(gcDestructors)`; under ARC/ORC the owner is a named variable
  whose scope covers every use, since temporaries die at the end of their
  statement.
- Macros that rewrite typed code (`inlineProcs` in `src/base/metaUtils.nim`)
  check what each call refers to before reading its definition: `getImpl` of a
  proc variable, such as the allocator's `outOfMemHook`, is a variable
  declaration. Under ARC/ORC the walked code also reaches system module procs
  that do not inline, such as `rawWrite`. Leave both kinds of call in place.
- Test configs (`.nims`) cannot see the memory manager: `defined(gcRefc)` and
  `compileOption("gc", ...)` are false there. Read the last `--mm:` argument
  (`tests/base/tpiogc.nims`).

## Checking

- `--expandArc:procName` prints a proc with the copies and destroys Nim
  inserted.
- `-d:useMalloc --passC:-fsanitize=address --passL:-fsanitize=address` stops at
  the first invalid heap access.
- Run the tests with several threads, in the default build and in one with
  `nimargs = @["--mm:refc"]` in `qexconfig.nims`. Races on reference counts
  crash in the allocator or in `rememberCycle`.
