# Experimental Graph Design

This note records the subtree-level contracts that are easy to miss when
reading one module at a time. Keep operator formulas, local mechanics, and
test-only conveniences in the code and tests.

## Scope

The reusable graph layers are `core.nim`, `scalar.nim`, `functional.nim`, and
`gauge.nim`.

`multi.nim` is support plumbing for multi-output operators rather than a second
general-purpose value language.

`hmcgauge.nim` is an application layer built on those surfaces. It owns
trajectory construction, training, and logging, but it should consume graph
contracts rather than redefine them.

## Core Model

The graph is value-centric: each node is a value plus a contract for how the
rest of the system is allowed to see it.

Raw `inputs` are only the obvious edges. They are not the whole dependency
surface. Every nontrivial node effectively answers three different questions:

- What must evaluation see to compute the current value?
- What structure is allowed to invalidate a cached symbolic gradient?
- What can affect the value or any derivative downstream?

Those are different questions, so the package keeps separate traversal modes for
them. Most correctness bugs in this subtree come from exposing the wrong set of
dependencies for one mode.

`backwardTarget` exists for the cases where reverse propagation is not well
described as "differentiate each raw input independently". Use it when laziness,
hidden dependencies, or operator-specific semantics make ordinary input-wise
backprop misleading.

`backwardTarget` is not additive with ordinary `backward`. `newGfunc` rejects
nodes that define both, so a node must choose either ordinary raw-input
propagation or target-directed propagation.

Custom dependency walks affect reachability, freshness, and cache signatures.
They do not turn hidden dependencies into raw backward inputs. If a hidden
dependency is differentiable, the node should use `backwardTarget`.

In `GdepWalks`, a nil hook means "walk raw inputs" for that traversal mode.
Only store a custom walk hook when raw inputs are not the right dependency
surface.

When node semantics depend on static metadata that is not present in raw inputs,
that metadata must be exposed structurally or encoded in the node's
`gfunc.signatureKey`. `multiSelect[i]` is an example: the selected slot index is
part of the operator identity.

Layer-specific node categories should stay in the layer that owns them. For
example, deferred-apply nodes are recognized in `functional/apply` by comparing
the node's `Gfunc` with that layer's singleton function objects; core `Gfunc`
does not carry an apply-specific kind field.

`Gfunc` is a compact operator record, not an abstraction layer. Its hooks and
metadata are read directly by graph internals. Keep computed accessors only when
they add behavior, such as preserving the signature-key fallback for function
records whose explicit key is zero.

Gradient-signature traversal can hide subtrees, but each visited node still
records the runtime-local ids of its raw input roots. This lets nodes such as
`cond` and deferred apply-partials keep selected internals lazy while still
making root identity changes visible to cache invalidation.

## Runtime Ownership

Every `Gvalue` belongs to exactly one `GraphRuntime`, and it must have that
runtime at construction time.

The runtime reference is non-nil by construction. Core code assumes this
invariant and uses `value.runtime` directly instead of rechecking it at each
call site.

There is no default runtime in this subtree. Root constructors such as
`toGvalue(grt, ...)` are explicit, and internal literal materialization must be
anchored to an existing runtime-bearing node.

The runtime owns the mutable state that gives a graph its identity and reuse
behavior:

- freshness epochs
- stable node ids
- gradient and apply caches
- cache stats and debug knobs
- runtime-local evaluation counters
- depth guards for callable resolution and apply-partial preparation

Those fields are mutable runtime state, not a separate service layer. Counters,
debug flags, depth limits, and cache stats should be updated directly at the use
site. A runtime helper is justified only when it validates runtime ownership,
keeps a multi-field update atomic, or gives a repeated cache operation one clear
failure model.

Cache reset helpers are runtime/core operations. Higher layers may use the
runtime fields they own, but should not re-export cache controls only to rename
where callers import them from.

Stable ids and cache entries are meaningful only inside one runtime. Mixing
nodes from different runtimes is a construction error, not a valid advanced
case.

`graphNode` does not assign runtimes. Its job is only to verify that the result
node and all inputs already agree on one runtime. Inputful graph nodes must also
have a graph function at construction time; zero-input prototypes may leave it
nil.

`newMultiOutputNode` is the special multi-output constructor that anchors the
carrier runtime from slot prototypes, then verifies input compatibility. Like
`graphNode`, an inputful carrier requires a graph function. Keep that rule local
to `multi.nim`.

The package is not thread-safe at the graph-runtime level. A runtime should be
treated as exclusive mutable state.

## Nil Checks

Ordinary graph values are non-nil by construction. After a value has crossed a
graph-construction boundary, code should use it directly rather than repeating
defensive `nil` checks at every traversal, validator, copy, or accessor site.
That includes public graph operators whose parameters are graph values: passing
`nil` is outside their contract unless that parameter explicitly documents `nil`
as meaningful.

Keep `nil` validation at the few places that deliberately cross from
unconstructed or erased data into graph storage, such as runtime attachment,
`graphNode` input checking, `newMultiOutputNode` slot/input checking, and
explicitly erased public APIs that must restore a concrete value type.

Inside graph internals, check for `nil` only when `nil` is part of the local
protocol. Current examples include optional `Gfunc` hooks, nil dependency-walk
hooks meaning "raw inputs", cache misses, unresolved callable reductions,
unbound callable wrappers, optional lambda construction fields before closure
resolution, root/upstream gradient conventions where `zb == nil` means "use
one-like", and HMC acceptance state where a missing accepted gauge has semantic
meaning.

Do not use `nil` as a second error channel for malformed graph structure. If a
node can contain an invalid ordinary input, catch that at construction. If a
hook is required to return a graph value, returning `nil` is a protocol error
and should fail near the hook result.

Do not add local prechecks merely to turn caller-supplied `nil` into a nicer
operator-specific error. That spreads a non-invariant through normal graph code
and makes the simple path harder to read.

## Freshness And Reuse

Leaf mutation is explicit. Callers mark semantic change with `updated`, and
derived nodes become current relative to the dependencies visible to evaluation.

The important design choice is that caches are structural, not value-based.

That means:

- ordinary leaf-value updates should usually preserve cached symbolic work
- changing callable identity or callable-boundary structure should invalidate
  cached reductions and gradients
- cache keys are tied to runtime-local stable ids, not object addresses or
  current numeric values

If a change makes plain captured-value updates routinely invalidate structural
caches, that change is usually fighting the design rather than refining it.

Gradient cache entries are per output node and signature. They store complete
adjoints separately from lower-level expansion records:

- a direct `findGrad` or `grad` cache hit is valid only for a target whose
  adjoint has been marked complete
- ordinary raw-input expansions may reuse a cached input contribution instead
  of rerunning the node's `backward` hook
- target-directed expansions remember which `(node, target)` pairs have already
  been expanded through `backwardTarget`

Those records are committed only after a gradient build succeeds, so failed
ordinary or target-directed expansions do not poison later builds. When a build
uses any `backwardTarget`, only the requested target adjoint is considered
complete; intermediate adjoints from that build are internal scratch unless a
later ordinary build proves them complete.

`gradIsolated` remains available for callers that intentionally want temporary
cache isolation. Core lazy operators such as `cond` and deferred `apply` should
normally call `grad` inside their target hooks so branch and partial gradients
can share the runtime's structural cache.

## Shape And Construction Contracts

This package prefers construction-time checks over late coercion.

Public graph-building operators should preserve concrete node types whenever the
result type is known. `Gvalue` is the erased storage type used by core hooks,
node input arrays, apply results, and other genuinely opaque boundaries; it is
not a general operator-dispatch surface.

Two consequences matter across modules:

- `copyCompatible` is a real semantic contract, not a convenience predicate
- operators that choose between alternatives, such as `cond`, require result
  shape agreement up front

Mixed numeric literals are supported only through explicit, runtime-anchored
helpers and overloads. They should never smuggle in a global runtime choice.

The remaining dynamic methods are value-prototype contracts used after type
erasure has already happened: allocation/copy compatibility, zero/one creation,
gradient accumulation, upstream scaling, hidden dependency walks, signature
tokens, and value display. Concrete graph math should be ordinary procs, and
erased values should be cast or validated at the boundary that knows the
expected type.

Reverse-mode root and upstream conventions are type-specific. Do not recover
generic erased arithmetic as a fallback; add the needed `addLike`, `scaleLike`,
`zeroLike`, or `oneLike` behavior to the value type that owns the convention.

Closure-body cloning assumes ordinary graph nodes can be reconstructed from
`newOneOf`, cloned `inputs`, and the original `gfunc`. Nodes with additional
structural state must either expose that state through inputs/signature metadata
or receive clone-aware handling.

## Adding An Exceptional Node

If a node is not "raw inputs plus ordinary backward", state all of these before
or near the implementation:

- eval dependencies
- grad-signature dependencies
- depend-mode dependencies
- whether it uses ordinary `backward` or target-directed `backwardTarget`
- which static metadata must enter `signatureKey`
- where erased raw inputs or upstream values are restored to concrete types
- whether `newOneOf + cloned inputs + gfunc` preserves the node

## `cond`

`cond` exists to express data-dependent choice without turning branch selection
into structural churn.

Its contract is:

- public selectors are scalar or int graph values
- branches must have compatible result shape, and public same-type calls keep
  that concrete branch type
- evaluation is lazy and follows only the selector plus the chosen branch
- dependency traversal sees the selector and both branches
- gradient-signature traversal only recurses through the selector, while raw
  branch-root identities still contribute to the cache signature

That split is deliberate. It keeps branch evaluation lazy while allowing
selector flips to reuse symbolic gradients when the branch roots themselves have
not changed.

Literal branch conveniences are only for scalar and int branches. Other value
types must pass explicit graph values on both branches. If a branch comes from a
boundary that returns erased `Gvalue`, such as `apply`, cast it before using it
as a typed branch.

Use `cond` for choosing between compatible values, not as a back door for shape
changes.

## Higher-Order Functions

`apply(fun, x)` is designed to stay lazy even when the callable is produced by
the graph itself.

The higher-order path therefore distinguishes between:

- symbolic callable boundaries, which decide structural identity
- concrete reduction, which instantiates a body only when needed
- deferred partials, which preserve laziness during gradient construction

Closures are normalized before cached reuse matters. That keeps capture
structure explicit and makes substitution-based instantiation the stable model.

Callable storage is deliberately flat. `Glambda` owns `param`, `body`, and
captured bindings directly. `Gwrapper` owns its result prototype and optional
bound value directly, with `kind` distinguishing local wrappers from callable
producer wrappers. Do not reintroduce subtype-only wrappers or getter/setter
pairs unless they enforce a new invariant at construction.

Callable placeholders must carry an explicit result prototype. A wrapper without
one is malformed rather than a wildcard callable.

`apply` returns an erased `Gvalue` because function-valued expressions can be
resolved only as the callable graph is reduced. Code that knows the callable's
result type should cast the apply result before feeding it into typed graph
operators.

Callable freshness is intentionally asymmetric:

- local wrappers expose the current binding symbolically
- callable wrappers represent the last produced callable until their producer is
  reevaluated

Binding a wrapper is a semantic update for both cases. A callable producer's
forward hook may install a binding whose hidden captures are newer than the raw
producer inputs, so the wrapper epoch must remain at least as new as both its
visible producer dependencies and the captures reachable through the bound
callable. Evaluation may raise a node epoch to match visible eval dependencies,
but it must not move an epoch backwards after a forward hook mutates the node.

That asymmetry is part of the public behavior. It is what lets the system reuse
structural work aggressively without pretending that all callables are live
views of their producers.

Deferred-apply evaluation traversal is apply-owned. The apply-partial eval walk
skips through deferred apply nodes using their depend-mode view, while concrete
nodes still expose depend-mode hidden deps so closure boundaries remain
reachable.

Deferred apply-partial nodes carry one target as a raw input, but traversal and
gradient-signature exposure deliberately follow only the base expression.
Target-specific partials are materialized from cached reduced expressions on
demand.

Apply cache lookup and partial materialization are local to the apply node code.
Keep the direct table operations there unless another caller needs the same
multi-step operation with the same error semantics.

## `Gmulti`

`Gmulti` is operator plumbing for multi-output nodes. It is not a general
product-value abstraction.

The key separation is between:

- forward slot storage, which holds concrete evaluated values
- symbolic slot selection, which builds a graph expression for one fixed slot

Slot indices are static operator metadata. Dynamic graph-valued indexing is out
of scope by design. That keeps selection local, keeps slot shape resolution
simple, and avoids turning multi-output carriers into a second general-purpose
data language.

Multi-output carriers may be heterogeneous. Slotwise combination and gradient
accumulation therefore use the contribution slot's algebra, not a single erased
`Gvalue` operator, and accumulation does not recheck the original carrier
prototype. Indexing a multi-output carrier still returns an erased fixed-slot
expression; callers that know the slot type should cast it.

If control flow must choose between slots, express that choice outside the
indexing operation.

`storedSlot` exposes stored forward slot state, not a symbolic graph node.
Callers should treat it as current only after evaluating the carrier or a
consumer.

### Shared-Compute `Gmulti` Patterns

Use `Gmulti` inside fused operators to share real work across related outputs or
input gradients. Keep the public API typed; the operator implementation owns any
`Gmulti` input/output carrier.

The packed carrier should have an operator-specific slot contract documented at
the fused operator, not a pile of one-call slot wrappers. Its backward should
build shared subexpressions once, then return slot gradients as a single
`Gmulti`. See `axexpmuly(a, x, y)` in `gauge/fused_ops.nim`: it packs
`[a, x, y]`, computes `[exp(a*x), exp(a*x)*y]`, reuses the saved exponential
for `y_bar`, and shares the exponential derivative for `a_bar` and `x_bar`.

## Gauge Layer

The gauge layer extends the graph model to gauge-field values and related
operators. Its primary differentiation contract is with respect to gauge fields.

Gauge graph construction should stay in concrete gauge/scalar/coefficient types.
Backward builders are the expected place to recover erased raw inputs or
upstreams, because each backward builder is tied to one forward operator and
knows the concrete operand types it stored.

`toGvalue(grt, gauge)` copies caller gauge storage into graph-owned storage.
After construction, graph state changes should happen through the returned graph
value and `updated`; mutating the original caller gauge must not silently change a
graph leaf.

`gaugeSnapshot` returns a copy, not a live view of graph-owned storage. Code
that intentionally mutates graph-owned gauge storage must use `mutateGauge` so
the freshness epoch is marked even though the storage update happens in place.
The top-level `gauge` module does not re-export `unsafeGaugeStorage`; direct
storage access is for gauge implementation modules that import `gauge/shared`.

Zero-valued gauges are ordinary zeroed storage, not a privileged semantic flag.
There is no separate zero-state fast path to keep in sync with the payload.

Gauge-action coefficients are graph values and can participate in coefficient
subgraphs, but the action layer does not promise differentiation of
`gaugeAction` or `gaugeActionDeriv` with respect to those coefficients.
Unsupported coefficient gradients should fail explicitly rather than degrade
into ambiguous behavior.

## `hmcgauge`

`hmcgauge` composes graph primitives into trajectory evaluation, acceptance,
integration, and training.

Its design role is operational:

- own trajectory-level graph construction
- keep mutation points explicit
- keep learned parameters paired with their gradient expressions rather than
  spreading that invariant across parallel arrays

Trajectory snapshot accessors return detached raw gauge storage. Accepted
trajectory commit copies that snapshot back into the graph-owned initial gauge
and marks the mutation there.

Integrator coefficient completion is intentionally narrow. For the force-gradient
families, callers either accept the default tuple or provide the full explicit
tuple; partial positional completion is out of scope.

`IntegratorCoeffs` is simple parsed configuration data. Validation belongs where
the coefficients are turned into a concrete integrator run spec; the data object
does not need variant-object protection or accessor wrappers just to represent
the parsed text.

This layer should consume graph-core contracts, not redefine them. In
particular, runtime identity, cache policy, and dependency semantics belong
below this layer.
