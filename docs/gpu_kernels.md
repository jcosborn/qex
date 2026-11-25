# Writing GPU kernels

A kernel is written once and compiled for the backend chosen with
`-d:Backend=...`: CUDA, HIP, SYCL, OpenMP offload, or the CPU (default). A
binary holds one backend. [GPU backends](gpu_backends.md) covers builds,
devices and the staggered solver; this page covers kernel code. Every
kernel of QEX launches through the construct described here.

## The construct

```nim
import backend/[accel, cgfield]

proc scale(r, x: Field; a: float) =
  onGpu(r.l):
    for s in gpuSites(r):
      r[s] := a*x[s]
```

| Form | Does |
|---|---|
| `onGpu(n): body` | launches the kernel, waits, copies results back |
| `let fin = onGpuNowait(n): body` | launches; `fin()` waits and copies back |
| `onGpuQueue(n): body` | launches a kernel whose variables stay on the device; `gpuWaitAsync()` waits |

Each form with a count `n` also takes a sub-group size, `onGpu(n, 16):
body`: sub-groups of 16 threads on Intel GPUs, where the compiler chooses
otherwise; the other backends ignore it.

`n` is the thread count of the launch on a device: a number, or a `Layout`
for its number of sites. Without `n` a device runs a default count. The CPU
backend runs the body in the host threads whatever `n` is.

The body is the code of one thread. It runs once in each thread, with
`gpuThreadNum` in `0..<gpuNumThreads`. Loops take their share of the work
through iterators:

- `gpuRange(m)`: the indices of `0..<m` of this thread.
- `for e in f`, `f.even`, `f.odd`, `f.all`: the sites of the field `f` that
  this thread owns (see [Field code](#field-code)).

A device thread gets every `gpuNumThreads`-th index, so neighboring threads
touch neighboring memory, and with as many threads as indices each thread
has one: a loop `for i in gpuRange(n)` over the launch count `n` is
compiled for a device as its one index, without the loop. A host thread
gets a contiguous share. The result of a kernel must not depend on the
thread count.

A device kernel declares its thread number and count at its start, so
templates, iterators and the procs that a device kernel copies may use
`gpuThreadNum`, `gpuNumThreads` and the iterators; a proc marked
`gpuInline` or `gpuCall` is a function of its own and may not. With the CPU
backend they are `threadNum` and `numThreads` of the host thread.

## Where a launch may stand

- Inside a proc. At the top level of a module, and directly in a `test`
  body of `std/unittest`, Nim emits the pragmas of a device launch at file
  scope of the generated C, where a device compiler stops.
- Inside a `threads:` block every thread of the team calls `onGpu`, and so
  every proc that launches a kernel. A device launches the kernel once; the
  CPU backend runs the body in the team, which meets at a barrier after a
  queued kernel.
- The sums of `backend/accel` (`sumFixed`, `gpuSum`, `SumX`) keep host
  state and stop with an error inside a `threads:` block, as the hop and the
  update of StagGpu do.
- Kernels that name a `Field`, a seq, a `Halo`, a `HaloLayout` or a `GpuSum`
  launch outside `threads:` blocks with a GPU backend: their device copies
  live in the table of [`backend/gpumem`](../src/backend/gpumem.nim), which
  belongs to the main thread.
- With OpenMP offload a kernel launched inside a `threads:` block completes
  at the launch.

## Sites and numbers

A site index has type `SiteV[L]`: one lane of an outer site on a device
(`L` = 1), a whole SIMD site on the CPU (`L` = `V`). `f[s]` of a `GpuField`
takes either and gives a value whose number type is the scalar type for one
lane and the SIMD type for a whole site. Kernel math written with the QEX
types (`Color`, vectors, matrices) and their operations compiles for both.
`gpuType(T)` is the device form of the element type `T`.

Kernels on arrays of reals, as those of the staggered operator, take an
element type instead. A kernel index stands for `W = gpuLanes(T, V)` lanes
of an outer site: one on a device, a hardware vector with the CPU backend.
With `E = gpuElem(T, W)` and `VE = V div W`, real `c` of index `k` of a
field with `NC` reals per site is at

    (k div VE)*(NC*VE) + c*VE + k mod VE

of the field memory read as an array of `E`: the layout of the lanes with
`(E, VE)` for `(T, V)`. `gpuTo(E, x)` gives a number, or a value of the
other precision, in `E`. A kernel written over `E` and `VE` is lane code on
a device and SIMD code on the CPU. Neighbors are the part that differs:
a table per outer site names the outer site of the neighbors and the swap
of its lanes (`perm`), see `dslash` of
[`physics/stagGpu`](../src/physics/stagGpu.nim).

## Field code

A kernel can be written as field code in a `threads:` block:

```nim
proc update(r, a, b: Field; m: Field2; c: float) =
  onGpu(r.l):
    r := c*a + b
    r += m*a
    r.odd := 0
    for e in r.even:
      r[e] := r[e] - m[e].adj*b[e]
```

- `for e in f` visits the sites of `f` that this thread owns: one lane per
  index on a device, a whole SIMD site with the CPU backend. `f.even`,
  `f.odd` and `f.all` visit the owned sites of those subsets.
- A thread owns the same sites in every loop over a field or its subsets:
  every `gpuNumThreads`-th lane on a device, a block of the even and a block
  of the odd outer sites with the CPU backend. Statements that read and
  write the same sites follow each other in one kernel without a barrier.
  A statement that reads other sites written in the same kernel, such as
  neighbors, needs a launch of its own: a device has no barrier across a
  launch.
- `:=`, `+=`, `-=` and `*=` take a field, a subset, or an expression of
  fields and numbers with `+`, `-`, `*`, unary `-` and `adj`. The
  expression is evaluated at each owned site with the operators of the
  site types.
- `f[e]` is the site value of the QEX types. On devices the products,
  sums, adjoints, numbers, `trace`, `redot`, `norm2`, `projectTAH` and
  `expAH` of QEX work in kernels, copied by `inlineProcs`; `exp` exceeds the
  iteration limit of the compiler's VM. A proc that `inlineProcs` copies
  must not declare types in its body: the copy drops type sections.
- A field captured by the body has its device copy in the GPU memory table;
  `gpuFieldView(lo, T, p)` is a field of the layout `lo` with element type
  `T` over device memory `p` from `gpuMalloc`, without a table entry.
- Outside kernels and `threads:` blocks, `x.norm2` and `redot(x, y)` of a
  field view or subset return the sum over the sites and the ranks; each
  launches a kernel and adds the values of its threads in a fixed order.

With the CPU backend such a kernel gives the bits of the same code in a
`threads:` block, apart from sums (their order differs) and functions whose
SIMD form treats the lanes together (`expAH` scales by the largest norm of
the lanes).

Rules that the CPU backend does not check:

- On a device `f[e]` is a view of the field: `let m = f[e]` refers to the
  field, where the CPU copies. Copy with `var m {.noInit.}: evalType(f[e])`
  and `m := f[e]` when the site changes before `m` is used. `adj` is a view
  everywhere: `r := r.adj` reads what it overwrites, as in `threads:` code.
- `gpuThreadNum`, `gpuNumThreads` and the iterators exist only in kernels
  and in the procs that a device kernel copies; elsewhere a device build
  stops in the C compile (`qexTid` undeclared).
- Loops over fields and loops over `gpuRange(m)` give a thread different
  shares; statements on the same data in one kernel use one kind.
- The launch copies `a[i]` of a loop `for e in a[i]` in the body (Nim passes
  it to the iterator by address, and the captures of a SYCL, CUDA or HIP
  kernel are constant); a template that loops over a field argument binds
  it to a `let` first.
- A kernel over a subset launches with the sites of the subset, for
  example `onGpu(lo.nEven): x.even := 0`; with `onGpu(lo)` half the threads
  wait.
- In a generic proc Nim looks up the `items` of `for e in f` where the proc
  is instantiated: generic code of a library writes `items(f)`, or every
  module that instantiates it imports `backend/cgfield`.

A CPU build with `-d:gpuSimd=false` runs field kernels lane by lane through
the views of lanes, as a device does.

## What a body may name

The launch captures every outer variable that the body names:

| Variable | In the kernel |
|---|---|
| has `toGpu(x)` (numbers, `Field`, seq, `Halo`, `GpuSum`, your own types) | `getGpu(x, toGpu(x))`; the finalizer calls `fromGpu(x, g)` |
| other plain data (pointers, objects and arrays without managed fields) | a copy |
| anything else (a string, a ref without `toGpu`) | the compile stops |

Types, constants, static parameters, routines, thread variables and
`importc` variables are used as they are.

- The captured values are immutable: a kernel writes through pointers and
  the views that `toGpu` returns, never to a captured variable.
- A name the body declares is its own. Give a local of the body another
  name than a host variable that the same body uses. A body may declare
  templates for its own use; they see the kernel values.
- A helper template defined outside the body sees the host variables it
  names, behind the back of the launch. Pass kernel data to such templates
  as arguments. Nim warns when a body and an outside template name the same
  variable ("a new symbol ... has been injected").
- `getGpu(x, g)` receives a stand-in of the type of `x`; it selects on the
  type and returns the device value.

To give a type a device form, define the three procs:

```nim
type
  Arr = object      # host: n reals at p
    n: int
    p: ptr UncheckedArray[float]
  Dev = object      # what the kernel gets
    p: ptr UncheckedArray[float]

proc toGpu(x: Arr): Dev = ...                 # allocate or look up, copy in
template getGpu(x: Arr, g: Dev): auto = g
proc fromGpu(x: Arr, g: Dev) = ...            # copy back
```

An object that owns device memory returns a view from `toGpu`, copies
nothing in `fromGpu`, and declares `gpuResident`:

```nim
template gpuResident*(T: typedesc[MyDeviceObject]): bool = true
```

`onGpuQueue` accepts numbers, plain data and such types. A variable that is
copied back, such as a seq or a `Field`, stops its compile.

## What a body may call

| Callee | With a device backend |
|---|---|
| template | expands in the body |
| proc marked `{.gpuInline.}` | called; always inlined; for small arithmetic helpers |
| proc marked `{.gpuCall.}` | called; the compiler decides on inlining; for larger procs such as RNG draws |
| `importc`, `importcpp` and magic procs | called |
| any other proc | its body is copied into the kernel (`inlineProcs`) |

The marks come from `getOptimPragmas()` of
[`base/basicOps`](../src/base/basicOps.nim). A marked proc is compiled for
host and device under CUDA and HIP (`QEX_HD`), so its body must be device
code: no host globals, no CPU SIMD types. An unmarked proc needs no
attribute, since its body lands in the kernel; generic QEX operations work
this way.

`inlineProcs` stops the compile for two forms it cannot copy: a proc that
returns a `var`, and a proc whose `emit` names a parameter or `result` in
its C text. Mark such procs.

With the CPU backend every proc is called.

## Sums

```nim
var gs = newGpuSum[float](n)        # room for kernels of up to n threads
onGpu(n):
  var r = 0.0
  for i in gpuRange(m): r += x[i]*x[i]
  gs.reduce r                       # once in every thread
echo gs.value
```

`T` of a `GpuSum[T]` is `float` or an array of floats. On a device the value
of each thread is stored and `value` adds them in an order fixed by the
thread count (`sumFixed`), so a rank's sum has the same bits on every device
backend. The CPU backend adds over the host threads. `gs.free` releases the
device block.

`SumX` adds exactly through integer words, with equal bits at any rank
count; `gpuSum` adds the values of a kernel loop in the fixed order. The
[staggered section](gpu_backends.md#cg-modes-and-communication) describes
both.

## Memory

- `gpuMalloc`, `gpuFree`: device memory. `gpuMallocHost`, `gpuFreeHost`:
  pinned host memory for fast copies. With the CPU backend both are host
  memory aligned to 64 bytes.
- `gpuMemCpyToGpu`, `gpuMemCpyToCpu` copy and wait; `gpuMemCpyToCpuAsync`
  queues a copy to pinned memory behind the kernels already launched.
- Host objects keep their device copies in the table of `backend/gpumem`,
  keyed by host address. `gpuFlagsIncl` and `gpuFlagsExcl` skip copies in
  and out. An object whose flags skip copies stays alive and in place, or
  its entry is freed with `freeGpuMem`.
- With the CPU backend a `Field`, a seq of plain values and an `alignedMem`
  are their own device copies.
- Device memory named by a queued kernel stays allocated until the kernel
  has completed.

## The queue

Kernels launch in order and run in order. `onGpuQueue` returns before its
kernel completes; `gpuWaitAsync()` waits for all queued work, and the reads
of the device sums wait as well.

With OpenMP offload:

- An `onGpu` kernel waits for the queue where it holds work and runs as a
  target without `nowait`, which completes at its launch. `gpuWaitAsync`
  skips its task wait outside threads blocks when nothing was queued since
  the last wait.
- A wait inside a `threads:` block covers the kernels queued in that block
  only: call `gpuWaitAsync()` before a block whose kernels read what kernels
  queued before it write.
- An `onGpuQueue` or `onGpuNowait` kernel is queued when the kernel reads
  numbers and pointers from outside only. A kernel that reads an object (a
  view with several fields, an array) completes at its launch, as an
  `onGpu` kernel: icx reads a captured object from the stack of the launch
  when the kernel runs. The launch decides on the typed kernel, so a name
  that a macro of the body only turns into other names (`pickB(so, j, C)`
  of stagGpu reads `so0`, `so1`, ...) counts as nothing. Read pointers and
  numbers in the kernels that should overlap with host work.

## The CPU backend

The CPU backend is the development platform and a target of its own.

- An `onGpu` kernel runs in a `threads:` block, or in the calling team, on
  whole SIMD sites. A launch outside a `threads:` block opens one, which
  costs 2 to 5 us of host time with 4 to 16 threads; inside one the team
  meets at a barrier.
- The staggered operator with links of 18 or 14 reals and the update of its
  CG compute on hardware vectors, for one system and for several at once. `-d:gpuSimd=false` runs every kernel lane by
  lane: the reference of the SIMD forms, equal bit for bit where the
  compiler keeps the order of the arithmetic and `fmadd` rounds once.
- `fmadd` of vectors is the fused instruction where the C target has it
  (`-march=native` on a host with FMA). Without the instruction the vectors
  round twice and the numbers, through libm, once.
- Device memory is host memory, the queue is empty, and every halo message
  uses QMP.
- A kernel that compiles here can still fail on a device: the rules of this
  page on captures and callees are checked by every backend, but host
  pointers reached through plain data are not.

## Speed on devices

Rules from measurements on H100, MI300X and PVC.

- Keep the code of a stencil direction in templates over small `gpuInline`
  procs. With the hops of the staggered operator as procs, clang merged the
  local and halo branches of a kernel and used four times the registers on
  MI300X.
- Pass loop bounds to kernel templates as expressions. A template that
  captured a copy of a bound changed the index arithmetic of a device
  kernel.
- Use `gpuCall` for large procs. Forced inlining of the RNG draws made the
  clang kernels 12 to 33% slower.
- Write a product that two code paths must round alike with `fmadd`. A
  compiler that contracts one path and not the other lets halo copies drift
  from their owners.
- Write `a*y + b*x` as `fmadd(a, y, b*x)` where results should not depend
  on the compiler. clang, nvcc and icpx chose which product to fuse by the
  code around the kernel: moving the staggered hop from an older launch to
  `onGpu` changed the CG residuals in the last bits.
- A kernel reads its captured values in place, as locals of the launch.
  Copies at the start of a kernel made nvcc hold them in ordinary registers
  and the staggered hop 2.5 to 4% slower on H100 and GH200.
- A copy to the host that waits costs 7 to 11 us after a kernel; queue it
  behind the kernels with `gpuMemCpyToCpuAsync` and wait once.
- At 16^4 on one device the host time of each launch shows in the time of
  the solver. Launches have profile timers only when built with
  `-d:gpuLaunchTimers`; on PVC and GB300 the timers took 1 to 2% of the
  time of the staggered operator.
- One atomic add per term to few slots makes short kernels wait; give
  `SumX` enough slots or use the fixed order there.
- QEX's `expAH` in a kernel took 12 to 14 times the time of the
  polynomial that gaugeGpu writes out on arrays (`mexp`) with SYCL on PVC:
  its matrix temporaries went to private memory, 92 KB per hardware thread.
  With OpenMP offload it took 1.01 to 1.15 times (`befield2`). gaugeGpu
  keeps `mexp` for the link update.

## Validation

1. Run the kernel with the CPU backend at 1 and 2 ranks; `tests/base/tongpu`
   and `tests/base/tgpumem` cover the construct.
2. Run the Nim phase for each device backend on any machine:
   `nimargs:'@["-d:Backend=CUDA","--compileOnly"]'` in the configuration.
3. Build and run on a device of each backend. The Nim phase does not see
   what only the device compiler checks.
4. Compare with the results of before bit for bit where the order of the
   arithmetic is unchanged, and time against the bandwidth or flop limit
   of the device.
