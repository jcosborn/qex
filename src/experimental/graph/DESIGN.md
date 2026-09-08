# Experimental Graph Design

This document records the contracts that are easiest to miss when reading one
module at a time. It is intentionally about subtree-level design; local formulas,
operator-specific details, and test-only conveniences should stay in code and
nearby tests.

## 1. Layers And Scope

The reusable graph layers are:

- `core.nim`: graph values, runtimes, traversal, evaluation, and reverse-mode
  gradient construction.
- `scalar.nim`: scalar/int graph values and scalar operators.
- `functional.nim`: first-class structural lambdas, `apply`, and structural VJP
  construction.
- `gauge.nim`: gauge-field graph values and gauge operators.

`multi.nim` is support plumbing for multi-output operators, not a second
product-value language.

`pghmc.nim` is an application layer. It owns trajectory construction,
training, and logging, but consumes the graph contracts instead of redefining
runtime identity, cache policy, or dependency semantics.

## 2. Core Graph Model

The graph is value-centric: each node is a value plus a `Gfunc` contract for how
other graph machinery may see it.

Raw `inputs` are only the obvious edges. Nontrivial nodes must answer three
separate dependency questions:

1. **Eval:** what must evaluation visit to compute the current concrete value?
2. **Grad plan:** what differentiable structure can contribute to a target?
3. **Reachable:** what can affect this value or any downstream derivative?

Those surfaces are deliberately separate. Most subtle bugs in this subtree are
wrong-dependency-surface bugs.

### Input views

`Gfunc.inputView` is the single dependency hook. It receives an
`InputWalkMode`:

- `iwmEval`: values needed to compute the current concrete value;
- `iwmReachable`: values that should appear in reachability/debug traversals;
- `iwmBackward`: differentiable deps used by reverse-mode planning and
  propagation.

A nil hook means "walk raw `inputs`" for every mode. Only install a custom hook
when raw inputs are not the right dependency surface.

Custom input views must emit ordinary graph values without mutating `inputs`.
The `backward` hook receives both the backward dep index and the actual `input`,
because an `iwmBackward` surface need not line up with raw input positions.

### Topology-stable evaluation

Evaluation is topology-stable. A `forward` hook may update only the current
node's cached concrete value or runtime-local value storage. It must not change
`inputs`, `gfunc`, unrelated nodes, or lambda producer bindings. `eval` does
not audit this on every run; module-owned hooks are covered by tests, and custom
hooks are trusted to obey the construction contract. A hook that rewires
topology during evaluation is invalid and may fail later near the bad mutation.

### Gfunc identity and metadata

`Gfunc` is a compact operator record, not an abstraction layer. Core and higher
layers read its hooks directly.

When semantics depend on static metadata not represented by raw inputs, close
over that metadata in the `Gfunc` hook that owns it. Graph topology and operator
metadata are fixed after construction; user updates change concrete values, not
which operation a node represents.

Layer-specific node categories stay in the owning layer. Real structural nodes
use hidden `Gfunc` subtypes owned by their modules, such as `Gcond` in core
conditionals and `Gapply` in the functional apply layer; callers should use the
exported predicates rather than inspect public names. Core `Gfunc` has no
apply/lambda kind field or lambda-result permission bit. The public `Gfunc.name`
string is diagnostic text only.

## 3. Runtime Ownership And Freshness

Every `Gvalue` belongs to exactly one non-nil `GraphRuntime` at construction
time. There is no default runtime. Root constructors such as `toGvalue(grt, ...)`
are explicit, and internal literal materialization must be anchored to an
existing runtime-bearing prototype.

The runtime owns mutable graph identity and generic reuse state:

- freshness epochs
- stable node ids
- symbolic graph revision
- gradient cache
- log-Jacobian chain cache
- functional state: apply instantiation cache
- cache stats and debug knobs
- runtime-local run counters

There is no string-keyed runtime extension registry. State used by more than one
graph layer, or needed to make cache invalidation visible, lives as an ordinary
`GraphRuntime` field. Functional state is grouped under `GraphRuntime.functional`
so the core runtime has one clear field for the functional module's cache and
configuration. The functional layer owns that state’s semantics, but the storage
is direct typed runtime data rather than an erased extension object.

Stable ids and cache entries are meaningful only inside one runtime. Mixing
runtimes is a construction error.

`graphNode` does not assign runtimes. It verifies that the result and checked
inputs already agree on one runtime, then installs `inputs` and `gfunc`.
Inputful graph nodes require a non-nil `Gfunc`; zero-input prototypes may have
no function.

Leaf mutation is explicit. Callers mark semantic changes with `updated`, and
derived nodes become current relative to their eval-visible dependencies.
`updated` clears the `staticZeroLeaf` marker because a formerly static zero leaf
has become ordinary mutable storage.

Static zero leaves are a construction invariant. `markStaticZeroLeaf` is only
valid for inputless leaves whose concrete value is zero. `isStaticZeroLeaf` is a
cheap marker-and-value query, not a structural validator; computed graph nodes
that currently evaluate to zero are not static zero leaves.

Caches are structural, not value-based:

- ordinary leaf/capture value updates should usually preserve cached symbolic
  work;
- changing graph topology, selected lambda identity, or lambda-boundary
  structure advances the runtime symbolic revision and invalidates cached
  reductions and gradients;
- cache keys are tied to runtime-local stable ids and the runtime symbolic
  revision, not current numeric values.

If a change makes ordinary captured-value updates routinely invalidate symbolic
caches, it is probably fighting the design.

Arbitrary in-place mutation of structural lambda bodies after cached use is not a
public freshness model. Build a fresh lambda or update captured leaves instead of
rewiring a normalized lambda body and expecting existing structural caches to
track that mutation. Lambda-ref binding/copying is a construction-time mechanism;
rebinding a lambda placeholder after `grad`, `vjpOf`, or apply instantiation has
built symbolic work is the same kind of unsupported topology mutation.

The package is not graph-runtime thread-safe. Treat a runtime as exclusive
mutable state; globals used to connect functional submodules are not a concurrency
contract.

## 4. Construction, Shape, And Nil Contracts

This package prefers construction-time checks over late coercion.

`copyCompatible` is a semantic contract. Operators that choose between
alternatives, such as `cond`, require compatible result shape up front.

Public graph operators should preserve concrete node types whenever the result
type is known. `Gvalue` is the erased storage type for core hooks, input arrays,
`apply` results, and genuinely opaque boundaries; it is not a generic operator
dispatch surface.

Mixed numeric literals are supported only through runtime-anchored helpers.
They must be anchored to the prototype that owns the literal's meaning. For
`apply(fun, 1)`, that means the lambda parameter prototype, not the lambda-valued
function expression itself.

Ordinary graph values are non-nil after construction. Keep nil validation at
boundaries that cross from unconstructed or erased data into graph storage:
runtime attachment, graph-node input checking, multi-output slot/input checking,
and public erased APIs that restore concrete value types. Inside graph internals,
check nil only when nil is part of the local protocol: optional hooks, cache
misses, unresolved lambda bindings, optional construction fields, or root/upstream
conventions.

Fixed-shape internal nodes should index the inputs they constructed. Do not add
per-use input-count wrappers for shapes owned by the same module; malformed
manual node mutation is user error and should crash or fail close to the bad
mutation. Keep checks that prevent silent wrong graphs, such as branch shape
compatibility, mixed runtimes, nil dynamic dependencies, and erased type
recovery.

Module-internal consistency should usually be a unit-test contract, not a
runtime `require...` guard repeated at every use. Add a runtime check only when a
bad state would otherwise fail silently, corrupt shared state, or cross a public
or erased boundary with misleading semantics.

Do not use nil as a second error channel for malformed graph structure. Invalid
ordinary inputs should be rejected at construction; hooks required to return a
value should fail near the hook result if they return nil.

## 5. Gradient Cache And Reverse Mode

`grad(dep, x)` builds a backward graph. It does not evaluate `dep`.

The core gradient engine:

1. walks each node's `iwmBackward` input view to mark nodes that can contribute
   to `x`;
2. visits marked nodes in reverse dependency order;
3. asks each node's ordinary `backward` hook for graph expressions for
   backward-dep adjoints;
4. accumulates those expressions with the target value's graph algebra.

Gradient cache entries are per output node and runtime symbolic revision. They
store complete adjoints only:

- a direct `grad` hit is valid only for a target whose adjoint is marked
  complete;
- intermediate nodes with complete cached adjoints can reuse those adjoints
  across later target builds for the same output/revision;
- pending adjoints are committed only after a successful build, so failed builds
  do not poison later builds.

Public cache lookup follows the same revision contract. `findGrad(input, output)`
returns only adjoints from the current runtime symbolic revision; stale entries
are treated as cache misses.

The symbolic revision changes only when internal symbolic lambda/VJP metadata is
registered. Ordinary `update` calls change concrete values and freshness epochs,
but do not invalidate symbolic gradient cache structure. Rewriting topology,
lambda bodies, bindings, or operator metadata after construction is outside the
contract; build a fresh graph instead.

Structural functional VJP bodies use uncached seeded builds: the public
output-gradient cache owns `grad(dep, target)` reuse, while seed-specific VJP
bodies avoid storing contributions whose meaning depends on a particular
upstream adjoint.

Conditional upstream gradients are split before calling a node's `backward`.
Static zero branches are skipped there, so inactive branches can produce guarded
VJP graphs without constructing or evaluating the inactive apply VJP.

## 5.1 Log-Jacobian Factorization (`logDetJ`)

Determinants do not localize across merging paths, so unlike reverse mode the
log-Jacobian engine works only on declared composition chains:

```text
v = w0 -> w1 -> ... -> wn = u
ln|det(du/dv)| = sum_k ln|det(d w_k / d w_{k-1})|
```

`Gfunc.logdet` is the per-op contract. For a step node `z` it returns
`(ld, via)` meaning:

- every backward path from `z` to values upstream of `via` passes through
  `via`, so `d z/d base = (d z/d via)·(d via/d base)` for any base below it;
- `ld` is a scalar graph expression for `ln|det(d z/d via)|`, the total
  derivative through every edge between `z` and `via` (a re-entrant context
  edge such as a subset-frozen staple built from `via` is inside that factor);
- `z` and `via` are shape-compatible (the step is square).

Structural verification checks squareness. For graph dependence it sees only
`iwmBackward`: `via` must be a backward dependency of the step, and no backward
path to the requested base may bypass it. This trusts the op's input view at the
same level as reverse mode.

An op whose local determinant formula holds graph-valued auxiliaries fixed must
itself reject any that depend on `via` (for stout: `ds`, `alpha`, or
coefficients). A re-entrant auxiliary is valid only when the owning op
internalizes it and proves a frozen/triangular decomposition, as the action-aware
subset stout step does. The local formula is op-owned math, like `backward`, and
must be pinned by tests.

Whether the map is globally invertible is the application's concern; the engine
needs only the factorization, and a singular step shows up as `-inf`/`nan` in
the ordinary output.

`logDetJ(u, v)` builds the chain sum as an ordinary scalar graph:

- it asks each chain node's hook once and memoizes hook results and per-base
  sums in `GraphRuntime.ldjCacheByNode`, keyed by the runtime symbolic
  revision, so repeated calls return the identical node and shared prefixes are
  reused; `resetLdjCache` drops all entries;
- `cond` distributes: each branch is its own chain, the selector contributes
  nothing (the same a.e. convention as cond's backward), and eval follows only
  the selected branch's sum;
- a static-zero `ld` is elided from the sum, so a chain of only static-zero
  steps remains a static zero; the empty chain (`u` is `v`) is a fresh static
  zero scalar;
- a chain node without a hook, or a chain that dead-ends on a leaf off the
  base, fails at construction.

`apply`/lambda values are not chains the engine can walk; `logDetJ` fails at the
apply node. Generic cloning reuses the same `Gfunc`, so its hook may capture
immutable configuration but must derive graph values, including `ld` and `via`,
from `z` and its inputs.

Actions on a graph-level flow u = f(v) take the generic form

```text
S_eff(v) = S(u) - logDetJ(u, v)
```

- The fused stout pair carries its logdet view as a structural input of its
  update view, so cloning preserves the pair without a side cache.
- A fused local correction without a graph-level flow, such as block5 coupling,
  remains valid and declares no factorization.
- Versus a telescoped local correction, the generic action costs one extra step
  forward and one extra pullback per force (2 of 3K kernels for K steps).
  Recover them by fusing each step's smear and logdet kernels, not by telescoping.

## 6. Exceptional Node Checklist

If a node is not "raw inputs plus ordinary backward", document or make obvious:

- eval dependencies;
- grad-plan dependencies;
- reachable dependencies;
- how ordinary `backward` maps raw inputs to symbolic contributions;
- whether the op declares a `logdet` factorization, and that its hook derives
  `ld` and `via` from `z` when the `Gfunc` can be shared or cloned;
- kernel capture and cloning: a `Gfunc`'s closures may capture only immutable
  configuration and layout-only scratch (`newOneOf` buffers and shifters).
  Anything holding input values (`Transporter.link`, saved field references, or
  a captured input node used instead of `v.inputs[i]`) must be rebound from
  `v.inputs` inside the kernel, because generic cloning (`newOneOf` + cloned
  inputs + the same `Gfunc`) shares those closures;
- where erased raw inputs/upstreams are restored to concrete types;
- whether `newOneOf + cloned inputs + gfunc` preserves the node, or whether it
  needs clone-aware handling.

Nodes with additional structural state must expose that state through inputs,
increment the runtime symbolic revision when internal symbolic metadata changes,
or implement clone-aware handling. Functional symbolic VJPs are ordinary graph
nodes, so their function and target dependencies must be visible through inputs.

## 7. `cond`

`cond` expresses data-dependent choice without turning branch selection into
structural churn.

Its contract is:

- selectors are scalar or int graph values;
- branches must have compatible result shape;
- same-type public calls preserve the branch concrete type;
- eval walks only the selector and selected branch;
- reachable and grad-plan walks see the selector and both branches;
- reverse mode gives the selector a zero adjoint and gives each branch the
  upstream adjoint guarded by `cond`.

This split keeps value evaluation lazy while letting `grad` build branch
backward graphs directly. Selector flips reuse symbolic gradients because branch
choice remains a runtime value.

Literal branch conveniences are only for scalar and int branches. Other value
types must pass explicit graph values on both branches. If a branch comes from an
erased boundary such as `apply`, cast it before using it as a typed branch.

Use `cond` for choosing between compatible values, not as a back door for shape
changes.

## 8. Functional Layer: Structural Lambdas And `apply`

The functional layer is a structural lambda calculus embedded in the graph. It
supports higher-order AD by rewriting lambda/application/VJP expressions into
ordinary graph structure.

A lambda value has type:

```text
Lambda<A, B> = paramProto A -> resultProto B
```

Both parameter and result prototypes are required. Result-only placeholders are
malformed because VJP shape depends on every argument level and the final result
cotangent.

### Lambda storage and normalization

`Glambda` owns `param`, `body`, and `captureParams` directly. Captured values
live in `Glambda.inputs`.

`GlambdaRef` owns `kind`, `paramProto`, `resultProto`, and optional `binding`.
`lrkLocal` refs are local placeholders; `lrkProduced` refs are symbolic
lambda-valued results.

Lambdas are normalized into capture form before cached reuse matters. Free
values move into `Glambda.inputs`, and `captureParams` stores fresh body
parameters standing for those inputs. This makes substitution-based
instantiation stable and gives ordinary traversal concrete edges for captures.
Resolved lambdas maintain `captureParams.len == inputs.len`; each pair means
"substitute this capture parameter with this captured value." Normalization,
cloning, and generated-lambda construction own that invariant, and downstream
code may index the two sequences directly.

Resolved `Glambda` values are structural graph values, not fake graph operations.
When captures exist, they install the internal `lambda captures` function only
to make capture freshness visible to ordinary evaluation. Captures still live
directly in `inputs`; the function is not a public lambda-producing operation.
Whole-lambda gradients are rejected; scalar capture gradients are built by
apply's structural VJP path.

`lambdaParam(paramProto, resultProto)` is the public way to build higher-order
lambda parameters. Produced lambda refs are internal placeholders used as result
shapes for structural lambda forms owned by the functional layer, such as real
`apply` and real `cond`.

### Apply

`apply(fun, arg)` is an ordinary graph node:

```text
apply : Lambda<A, B> * A -> B
```

Graph construction does not eagerly instantiate lambda bodies. Bodies are
instantiated when the apply node is evaluated, or when the backward graph builder
structurally transforms a direct lambda body.

`apply` returns erased `Gvalue` because function-valued expressions may resolve
only through structural lambda graph forms. Code that knows the result type
should cast before passing the result to typed graph operators.

For scalar/int literal overloads, literals are materialized from the lambda
parameter prototype. This is required for structural, non-direct lambdas such as
`cond(k, f1, f2)`: the function expression itself is lambda-valued and is not the
argument prototype.

`apply` does not perform construction-time argument compatibility checks. The
node records a symbolic lambda application from the function's result prototype;
body instantiation decides what is actually demanded. If the body never uses an
incompatible argument, evaluation and gradients may still be well-defined. If a
demanded path treats a non-lambda as a lambda or copies an incompatible result
shape, that path fails when it is evaluated or differentiated.

### Apply traversal and caches

Apply traversal is apply-owned. Apply nodes always store the function and
argument. During structural VJP construction they may also store active value
targets as explicit inputs; those targets are construction-time graph edges, not
a later core refresh.

Eval, reachable, and backward input views compute lambda-visible surfaces
dynamically. Apply VJP construction passes active structural targets through an
explicit build context. Any dependency needed by a generated apply must become
an ordinary graph input on that apply before the generated graph is used. Lambda
clone remains substitution; when generated VJP bodies need active value targets
on cloned apply nodes, the VJP builder supplies an explicit source-node to
extra-input map. The clone module appends only those data-listed inputs.

Some function-side structure is needed only to invalidate stale gradient caches,
not to receive a raw adjoint. Apply cache entries key on the selected nominal
lambda closure and the runtime symbolic revision. Whole lambda values are not
differentiable raw inputs; function-side structural identity invalidates stale
pullbacks without introducing first-class lambda cotangents.

Apply instantiation and VJP construction remain local to the apply node code, and
plain `apply` construction must not prepare VJPs eagerly. Backward/VJP builders
own that work. Keep direct runtime apply-cache table operations local to
`ensureInstantiation` unless another caller needs the exact same multi-step
operation and failure semantics.

## 9. Structural VJP Transformation

`vjpOf(fun)` is a symbolic lambda-valued graph node only when the function cannot
be lowered immediately. Direct lambdas lower to ordinary generated lambdas.
Unresolved placeholders keep an explicit `vjpOf`/`vjpOfResult` node whose inputs
are the function expression and, for capture/result targets, the differentiated
graph value. A single apply-owned reducer lowers those delayed nodes once lambda
substitution or evaluation makes the source lambda available.

```text
vjpOf(fun)                 == callVjpOf(fun)
internal vjpOf(fun,target) == captureVjpOf(fun, target)
internal vjpOfResult(...)  == resultVjpOf(fun, target)
```

The final graph may contain only ordinary constructs plus unresolved symbolic VJP
nodes:

```text
lambda(...)
apply(...)
cond(...)
scalar ops
regular custom scalar/multi/gauge graph nodes
vjpOf(...)
vjpOfResult(...)
```

These graph-visible constructs must not exist:

```text
lambdaVjp
GlambdaCotangent
lambdaCapture
lambdaCotangent
lowerLambdaCotangent
structuralVjps metadata on lambdas
```

The builder uses three structural VJP expressions:

```text
callVjpOf(fun)
captureVjpOf(fun, target)
resultVjpOf(fun, target)
```

These names are documentation notation for the type rules below; in code they are
a single `VjpSpec` record whose `kind`/`target` fields select the form:
`callVjpOf` is `VjpSpec(kind = lvkCall, target = lvtkArgument)`, `captureVjpOf` is
`VjpSpec(kind = lvkCall, target = lvtkValue)`, and `resultVjpOf` is
`VjpSpec(kind = lvkResult, ...)`. On-graph they reduce to the `vjpOf` /
`vjpOfResult` nodes (`gvjpOfCall` / `gvjpOfResult`).

Targets are represented as one value, `LambdaVjpTarget`: either the call
argument or a concrete graph value target. Keeping the target kind and value
together avoids parallel fields that can drift out of sync.

For `fun : A -> B`:

```text
callVjpOf(fun) : A -> Cotangent<B> -> Cotangent<A>

callVjpOf(fun) =
  lambda(arg,
    lambda(seed,
      dArg))
```

For a scalar capture or other scalar graph target:

```text
captureVjpOf(fun, target) =
  lambda(arg,
    lambda(seed,
      dTarget))
```

For a function returning another function:

```text
fun : A -> (B -> C)

resultVjpOf(fun, target) =
  lambda(a,
    lambda(b,
      lambda(seed,
        dTarget)))
```

Deeper returned lambdas add one argument lambda per returned function level
before the final seed lambda.

For a direct scalar-result lambda:

```text
fun = lambda(x, body)

callVjpOf(fun) =
  lambda(x2,
    lambda(seed,
      gradSeeded(instantiate(body, x -> x2), x2, seed)))
```

The capture VJP is the same transform with a different target:

```text
captureVjpOf(fun, target) =
  lambda(x2,
    lambda(seed,
      gradSeeded(instantiate(body, x -> x2), target, seed)))
```

For unresolved lambda refs, this type rule is prototype-driven: the seed
parameter is `Cotangent<B>`, derived from the result prototype, not from
`Cotangent<A>`.

For lambda-valued bodies, VJP construction descends into the returned lambda:

```text
callVjpOf(lambda(x, lambda(y, body))) =
  lambda(x2,
    lambda(y2,
      lambda(seed,
        gradSeeded(
          instantiate(body, x -> x2, y -> y2),
          x2,
          seed))))
```

For scalar-result application:

```text
z = apply(fun, arg)

dArg = apply(apply(callVjpOf(fun), arg), seed)
```

For scalar targets captured by `fun`:

```text
dTarget = apply(apply(captureVjpOf(fun, target), arg), seed)
```

For an application whose function position is itself lambda-valued:

```text
fun2 = apply(hof, hofArg)
z = apply(fun2, arg)

captureVjpOf(fun2, target) =
  apply(resultVjpOf(hof, target), hofArg)
```

If `target` is inside an ordinary scalar `hofArg`, the cotangent for `hofArg` is
propagated by scalar reverse mode. If `hofArg` is lambda-valued, there is no
first-class lambda cotangent; unresolved structural VJPs stay as symbolic
`vjpOf` nodes until lambda substitution replaces the formal function with the
actual lambda expression.

### Symbolic VJP substitution

Instantiating a lambda with a lambda-valued argument substitutes ordinary graph
nodes. A symbolic VJP is just another graph node, so the function input is
substituted by normal clone/memo rules:

```text
hof = lambda(f, body)
actual = lambda(v, actualBody)
instantiate(hof, actual)

f                         -> actual
vjpOf(f)                  -> vjpOf(actual)
vjpOf(f, target)          -> vjpOf(actual, target)
vjpOfResult(f, target)    -> vjpOfResult(actual, target)
```

The same rule applies to captured lambda values. If normalization turns a free
lambda value into a capture parameter, symbolic VJP nodes that reference that
value are cloned to reference the capture parameter. Later instantiation replaces
the capture parameter with the actual structural lambda value.

No lambda value owns persistent VJP metadata. The graph expression carries the
function, target, and result prototype directly; cloning remaps those ordinary
inputs with the rest of the body.

### Conditionals and recursion

Conditional lambdas remain symbolic:

```text
f = cond(k, f1, f2)
z = apply(f, x)
```

Graph construction must not rewrite that to:

```text
cond(k, apply(f1, x), apply(f2, x))
```

VJP expressions distribute over the conditional:

```text
callVjpOf(cond(k, f1, f2)) =
  cond(k, callVjpOf(f1), callVjpOf(f2))

captureVjpOf(cond(k, f1, f2), target) =
  cond(k, captureVjpOf(f1, target), captureVjpOf(f2, target))

resultVjpOf(cond(k, f1, f2), target) =
  cond(k, resultVjpOf(f1, target), resultVjpOf(f2, target))
```

Only the selected branch is evaluated at eval time.

Self-application and other recursive lambda expressions use shared generated VJP
shells. The shell is registered before its body is transformed; recursive uses
that reduce to the same nominal lambda closure, VJP kind, target, and compatible
prototypes reuse the active shell. Recursive VJP correctness does not depend on
matching lambda body structure. Generated lambda normalization preserves active
VJP shells so recursive references stay shared instead of cloning open-endedly.

The active recursive VJP stack lives in an `ApplyVjpBuildCtx` passed through
apply/VJP construction. It is not process-global state; independent VJP builds
must not inherit active targets or memo entries from each other.

Persistent cache invalidation is handled by the runtime symbolic revision.
Active recursive VJP sharing is nominal and build-local: selected lambda
identity, VJP kind, target, and compatible prototypes identify the shell while a
VJP transform is running. Completed memo entries are keyed by selected lambda and
capture identity so same-shaped captured lambdas cannot reuse stale branch
pullbacks.

## 10. Custom Graph Functions And Lambda Values

Custom `Gfunc` nodes may appear inside lambda bodies when they produce ordinary
scalar, multi, or gauge values and expose normal forward/backward behavior.

Custom `Gfunc` nodes should produce ordinary values. Lambda-valued behavior is
represented with structural lambda forms:

```text
lambda(...)
apply(hof, arg)
cond(k, f1, f2)
lambdaParam(paramProto, resultProto)
```

Produced lambda refs and symbolic VJP nodes are internal functional-layer data. A
custom `Gfunc(name: "apply")` or `Gfunc(name: "cond")` is still an
opaque function; names are not structural authority.

Representative acceptance checks are:

```text
vjpOf(lambda(x, x * x))
vjpOf(vjpOf(lambda(x, x * x)))

hof = lambda(f, lambda(x, Gscalar(apply(f, x)) * Gscalar(apply(f, x))))
g = lambda(v, a * v + 1.0)
z = apply(apply(hof, g), x)
dzda = grad(z, a)
d2zda2 = grad(dzda, a)
```

Those should evaluate as ordinary higher-order graphs. Direct lambdas should
lower to generated lambdas, while unresolved lambda placeholders may keep
symbolic `vjpOf` nodes until apply/eval selects a concrete lambda expression.

## 11. `Gmulti`

`Gmulti` is operator plumbing for multi-output nodes. It is not a general product
value abstraction.

The key separation is:

- forward slot storage holds concrete evaluated values;
- symbolic slot selection builds a graph expression for one fixed slot.

Slot indices are static operator metadata. Dynamic graph-valued indexing is out
of scope. That keeps selection local, keeps slot shape resolution simple, and
avoids turning multi-output carriers into another value language.

Multi-output carriers may be heterogeneous. Slotwise combination and gradient
accumulation use the contribution slot's algebra, not a single erased `Gvalue`
operator. Indexing a carrier returns an erased fixed-slot expression; callers
that know the slot type should cast it.

If control flow must choose between slots, express that choice outside indexing.
`storedSlot` exposes stored forward slot state, not a symbolic graph node; treat
it as current only after evaluating the carrier or a consumer.

Use `Gmulti` inside fused operators to share real work across related outputs or
input gradients. The packed carrier should have an operator-specific slot
contract documented at the fused operator. Its backward should build shared
subexpressions once, then return slot gradients as a single `Gmulti`.

## 12. Gauge Layer

The gauge layer extends the graph model to gauge-field values and related
operators. Its primary differentiation contract is with respect to gauge fields.

Gauge graph construction should stay in concrete gauge/scalar/coefficient types.
Backward builders recover erased raw inputs or upstreams by direct cast, because
each backward builder is tied to one forward operator and knows the operand
types it stored. A wrong internal cast is a developer error. Let it fail instead
of adding wrappers that only restate the type contract.

`toGvalue(grt, gauge)` copies caller gauge storage into graph-owned storage.
After construction, graph state changes should happen through the returned graph
value and `updated`; mutating the original caller gauge must not silently change a
graph leaf.

`gaugeSnapshot` returns a copy, not a live view. Code that intentionally mutates
graph-owned gauge storage must use `mutateGauge` so freshness is marked even
though storage changes in place. The raw `Ggauge.gval` storage field is exported
for gauge implementation modules that import `gauge/shared`; the top-level gauge
module does not re-export it, and public writers should use `update`/`mutateGauge`.

Zero-valued gauges are ordinary zeroed storage, not a privileged semantic flag.
There is no separate zero-state fast path to keep in sync with the payload.

Gauge-action coefficients are graph values and may participate in coefficient
subgraphs, but the action layer does not promise differentiation of
`gaugeAction` or `gaugeActionDeriv` with respect to those coefficients.
Unsupported coefficient gradients should fail explicitly rather than degrade into
ambiguous behavior.

## 13. `hmcgauge`

`hmcgauge` composes graph primitives into trajectory evaluation, acceptance,
integration, and training.

Its design role is operational:

- own trajectory-level graph construction;
- keep mutation points explicit;
- keep learned parameters paired with their gradient expressions rather than
  spreading that invariant across parallel arrays.

Trajectory snapshot accessors return detached raw gauge storage. Accepted
trajectory commit copies that snapshot back into the graph-owned initial gauge
and marks the mutation there.

Integrator coefficient completion is intentionally narrow. For force-gradient
families, callers either accept the default tuple or provide the full explicit
tuple; partial positional completion is out of scope.

`IntegratorCoeffs` is parsed configuration data. Validation belongs where the
coefficients are turned into a concrete integrator run spec; the data object does
not need variant-object protection or accessor lambda refs just to represent
parsed text.
