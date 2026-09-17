# Experimental Graph Design

Contracts shared across modules. Local kernel details belong in code and tests.

## 1. Layers And Scope

| Module | Owns |
| --- | --- |
| `core.nim` | Values, runtimes, traversal, evaluation, reverse mode |
| `scalar.nim` | Scalar/int values and operators |
| `functional.nim` | Structural lambdas, `apply`, VJP construction |
| `field/types.nim` | Shared field and field-collection value storage |
| `gauge.nim` | Gauge values and operators |
| `nn.nim` | Real channel fields, parameter arrays and neural operators |
| `plan.nim` | Private execution graphs, joint scheduling, reusable buffers; import explicitly |
| `multi.nim` | Fused-operator slots and result views |
| `hmcgauge` | Trajectories, sampling and training |

Drivers such as `pghmc.nim` supply configuration, actions and callbacks. `Gmulti`
is operator plumbing; it does not define a general product-value language.

## 2. Core Graph Model

A node is a `Gvalue` plus a `Gfunc` and declared dependencies. Three views serve
three different operations:

| View | Question |
| --- | --- |
| Eval | Which values compute the current result? |
| Backward | Which differentiable paths contribute to the target? |
| Reachable | Which values can affect this result or a later derivative? |

### Input views

`Gfunc.inputView(v, mode, visit)` defines these walks; nil means raw `inputs`.

| `InputWalkMode` | Surface |
| --- | --- |
| `iwmEval` | Current numerical dependencies, preserving lazy branch order |
| `iwmReachable` | Reachability, diagnostics and plan source discovery |
| `iwmBackward` | Reverse planning and propagation |

Write $E(v),R(v),B(v)$ for the corresponding dependency collections. Plans require

$$
E(v)\subseteq\operatorname{inputs}(v)\cup R(v).
$$

Hooks emit ordinary graph values without mutating inputs. `backward` receives
both the backward index and actual input: backward positions need not equal raw
positions. `cond` evaluates one branch; its other views include both.

`graphValues` follows raw/reachable edges and symbolic lambda bindings. At a
resolved lambda it exposes captures, without entering the bound body. This audit
does not alter $B(v)$. There is no stop-gradient operator: partial derivatives
use distinct slots (section 5), while live dependencies remain visible for later
orders.

### Topology-stable evaluation

A forward may change its cached value or runtime-local storage. It must not
change `inputs`, `gfunc`, unrelated nodes or lambda producer bindings. Construction
and tests enforce this contract; ordinary eval does not repeatedly audit topology.

### Gfunc identity and metadata

`Gfunc` is a directly read operator record. Hooks may capture immutable operator
metadata; topology and metadata are fixed after construction. Numerical updates
do not change the represented operation.

Owning modules identify structural nodes through private subtypes/predicates
(`Gcond`, `Gapply`, etc.). `Gfunc.name` is diagnostic only. Core has no lambda/apply
kind field, lambda-result permission bit or equivalent name-based dispatch.

## 3. Runtime Ownership And Freshness

Each value has one non-nil runtime at construction:

$$
\operatorname{runtime}(z)=\operatorname{runtime}(x_i)\quad\text{for all checked inputs }x_i.
$$

There is no default runtime. `toGvalue(rt, ...)` is explicit; literals inherit the
runtime of their meaning-bearing prototype. `graphNode` checks this equality and
installs `inputs`/`gfunc`; it does not assign runtimes. Inputful nodes require a
function; inputless prototypes may have none.

`GraphRuntime` directly owns epochs, stable IDs, symbolic/boundary revisions,
gradient/logdet/apply caches, execution/apply/workspace frames, statistics and
debug state. Functional caches live in its typed `functional` field. IDs and
cache entries are runtime-local; mixed runtimes fail construction. No erased
extension registry is needed. A runtime is exclusive mutable state, not a
concurrent execution interface.

### Value storage and ordinary evaluation

| State | Meaning |
| --- | --- |
| `valueReady` | Successful evaluation |
| `epoch` | Numerical freshness |
| `hasStorage` | Resident payload |
| `valueOverride` | Explicit replacement of a computed value |
| `stale` | Invalid published plan result; ordinary eval must refuse it |

For the maximum epoch $e_{in}$ visited by the eval walk, a non-stale value needs
restoration/evaluation when

$$
\neg\mathrm{valueReady}\;\lor\;\neg\mathrm{hasStorage}\;\lor\;\mathrm{epoch}<e_{in}.
$$

Computed field/gauge outputs start as descriptors: `graphNode` releases storage;
eval acquires it before forwarding. `toGvalue` copies authoritative inputs into
resident leaves. Ordinary evaluation retains completed payloads and work.

| Event | Freshness/storage effect |
| --- | --- |
| First eval of an inputless producer | Run its forward |
| Restore unchanged computed storage | Run the producer; preserve its epoch |
| Explicit update | Mark ready/overridden, advance epoch, clear `stale` and constant-restoration markers |
| Override expires through a newer input or missing storage | Run original producer, clear override, advance epoch |
| Forward fails | Leave unready; permit retry |

`ensureStorage` restores generated zeros/identities from shape and their rule.
`update`, `mutateGauge` and `mutateField` remove that rule; mutation evaluates
before exposing storage and marks freshness afterwards. Restore a discarded
authoritative leaf with `update` before use.

`releaseStorage` replaces descriptors/sequences; escaped raw references retain
the old allocations. `releaseWork` drops communication scratch. The next forward
recreates/rebinds work. Stout `expa` and `m` belong to the payload and participate
in complete copy/alias/storage operations. Shift, hop and halo work always binds
the current descriptors.

`multiValues` owns numerical slot wrappers whose evaluated payloads alias inputs;
symbolic `bundle[i]` selects that original input. Structural slots retain their
input references, including through clone/forward/release. Mixed nested bundles
own separate numerical wrappers; this introduces no function-valued cotangents.

### Operator buffer contracts

| Mode | Forward contract |
| --- | --- |
| `bmFull` | Assign every numerical component; borrow declared inputs synchronously |
| `bmZero` | Destination buffers are cleared before every forward, including retry/override replacement; assign scalar slots explicitly |
| `bmAlias` | Retain raw `aliasInputs`; an empty list means all inputs |
| `bmOpaque` | Ownership unspecified; retained inputs and their ancestry need dedicated storage |

Value families implement `bufferProto`, `bufferCompatible`, `bindBuffer`,
`clearBuffer` and `bufferBytes`. Prototypes include shape/layout/type; packed
producer payloads include all cached components. Nil `bufferProto` declines
pooling. `bmZero` clears dedicated buffers too.

An in-place transfer requires all of:

- `bmFull` and a proof covering every simultaneously aliased operand;
- the input's last use;
- compatible concrete storage;
- no active alias or external owner.

`updated` clears `staticZeroLeaf`: a written constant becomes ordinary mutable
storage. `markStaticZeroLeaf` is valid only for an inputless concrete zero.
`isStaticZeroLeaf` checks marker and value, not arbitrary graph structure.

| Change | Cache consequence |
| --- | --- |
| Ordinary leaf/capture value update | Advance numerical freshness; retain symbolic work |
| Supported lambda binding/internal symbolic metadata change | Advance symbolic revision |
| Conditional lambda selection | Select a different instantiation key |
| New/expired computed override or generated-constant mutation | Advance boundary revision for source audits |

`GlambdaRef.valCopy` accepts compatible produced refs and changes the symbolic
revision. It rejects copying a resolved `Glambda` body. Rebuild derivative
expressions after rebinding. Arbitrary rewrites of cached topology/bodies are
unsupported; build a fresh graph instead.

`tgvalue`, `tgvalueu1`, `tgstorage` and `tglifetime` exercise these contracts.
Graph test helpers reserve command-line arguments for QEX, disabling unittest
name filtering; parameterized validation must still report executed cases.

### Shared execution plans

Construct all derivatives before planning; include every result needed together:

```nim
import graph/[core, scalar, gauge, plan]
let rt = initGraphRuntime()
let x = rt.toGvalue(g)
let f = norm2(x*x)
let p = plan(f, grad(f, x))
discard p.eval()
let value = Gscalar(p[0]).sval
let gradient = Ggauge(p[1])
```

`p.eval()` publishes numerical roots in order as detached leaves. Build further
derivatives from original expressions. Apply function-valued results, including
function slots inside carriers, to numerical arguments before publication.

The plan owns private clones, arena and caches; ordinary numerical sources keep
their storage. Ordinary eval uses the original graph. Active plan frames redirect
nested eval into the plan schedule.

For published wrappers, validity uses both fields:

| Transition | `stale` | `valueReady` | Result identity |
| --- | --- | --- | --- |
| Construction | true | false | New wrapper |
| Execution begins | true | false | Current wrapper |
| Publication via `updated` | false | true | Current wrapper; epoch advances |
| Validation, execution or publication fails | true | false | Every current wrapper, including partially published results |
| Rebuild retires outputs | true | false | Old wrappers will never be republished |
| Unchanged `p.eval()` | preserved | preserved | No forwards |
| `p.clear()` | preserved | preserved | Next eval rebuilds |
| Rejected nested plan call | preserved | preserved | Entry rejected before mutation |

Readiness invalidates dependent plans' cached feeds; `stale` prevents ordinary
eval from silently readying an inputless leaf. Scalar publications obey the same
rule. Ordinary computed nodes may be unready without being stale publications.

Published storage is borrowed until the next execution, including failure. Retry
through `p.eval()`: same-generation retry repairs existing wrappers/consumers;
rebuild requires reacquiring `p[i]` and rebuilding or rebinding consumers, captures
and dependent plans. Reacquisition does not redirect old edges.

Raw `.sval`/`.gval`/`.fval` reads bypass this guard. So does `storedSlot`: even eval
of the extracted slot has no dependency on the published carrier. Use symbolic
`pub[i]`, or evaluate the carrier before immediate stored-slot reads. Copy scalars
or take an owned `gaugeSnapshot` while valid to survive another execution. Raw
references preserve allocations, not their contents against arena reuse.

| Planner invariant | Requirement |
| --- | --- |
| Dependency lifetime | All inputs live through the forward; alias inputs live through their last consumer |
| Release | Detach the private descriptor before reusing its slot; cancel unselected edges lazily |
| Pool | Homogeneous concrete allocations; compatibility checked in both directions, without assuming transitivity |
| Rejected allocation | Requesting node stays dedicated for this private generation |
| In-place ownership | Allocation keeps its original pool; only the current owner returns it |
| Opaque consumer | May dedicate a live allocation |

Computed numerical overrides are external boundaries. Ordinary eval expires them;
a boundary change rebuilds the plan. Abstract lambda bodies ignore concrete
override cuts: each body gets fresh binders and an expression memo; captures use
the outer context. Generated constants and bound refs are isolated by private
cloning. Mutating an original constant or rebinding a supported source ref
invalidates that representation.

Boundary revision gates complete source scans. Ordinary numerical updates need
feed/result epoch and residency checks only. Validation checks the revision again
after evaluating external overrides, which may expire during that evaluation.

Nested `apply` retains outer arguments through the result copy. Private variants
use `(applyId, selectedLambdaId, symbolicRevision)` and are never mirrored into
ordinary apply caches. Plans execute serially within a runtime.

Nested `eval` registers the supplied object; it does not clone captured original
expressions. A captured computed node eligible for pooling can therefore lose
its storage when released. Initialized ordinary numerical leaves are external;
opaque/protected writers retain dedicated storage. Forwards must derive graph
work from their cloned declared inputs and expose its dependencies.

Communication work rebinds current inputs/destinations and completes exchanges
before returning. Plaquette/staple workspace keys include operation, layout,
field count and derivative order. Sharing requires both eval and workspace
frames; external-boundary eval uses its own work. Shifts/hops complete sends and
receives before input reuse and declare no in-place permission.

`clear` releases private caches and arena ownership, preserving published
allocations and validity. The next eval retires wrappers and rebuilds. Failures
restore runtime frames; result leaves cannot reconstruct source expressions.

`GraphPlan.stats` separates arena bytes, workspace bytes, source audits and
forward counts. Chain fixtures require two arena buffers, or one with proven
in-place permission, excluding preserved input storage. General derivatives may
need more. Metadata/communication ownership remains separate from the arena.

## 4. Construction, Shape, And Nil Contracts

| Operation/boundary | Contract |
| --- | --- |
| `copyCompatible` | Semantic result-shape compatibility, checked before constructing choices such as `cond` |
| `newOneOf` | Preserve producer subtype and kernel metadata |
| `valueLike` | Preserve mathematical value family, without producer caches |
| `cond`, `apply`, `slotVar`, numerical slots/selections | Use `valueLike`; e.g. copying a stout result yields an ordinary gauge |
| Owned cache slots | Use `newOneOf`; copy/alias every required payload component |
| `Gvalue` | Erased hook/input/apply boundary, not generic operator dispatch |
| Literal materialization | Anchor to the relevant prototype; for `apply(fun,1)`, use the parameter prototype |

Ordinary values are non-nil after construction. Validate nil/runtime/type/shape
at public or erased boundaries: runtime attachment, graph inputs, multi slots,
branch compatibility and erased type recovery. Internal nil checks belong only
to real protocols: optional hooks, cache misses, unresolved refs, construction
fields or root/upstream conventions.

Fixed-shape internals index what they constructed. Malformed manual mutation is
user error; do not add wrappers or repeated guards merely to restate internal
shape/type invariants. Keep checks that prevent silent wrong graphs/corruption;
otherwise prefer tests. Nil is not a second failure channel: invalid inputs and
required hook results must fail near their boundary.

## 5. Gradient Cache And Reverse Mode

For $z=f(x_1,\ldots,x_n)$ and output cotangent $u$,

$$
\bar x_i=(D_i f)^*u,\qquad
\bar x=\sum_{i:x_i=x}\bar x_i.
$$

`grad(dep,x)` constructs this graph without evaluating `dep`: walk `iwmBackward`,
mark paths to the target, traverse them in reverse order, call backward hooks,
and accumulate with each target's algebra.

The root cotangent is implicit: without a seed a root uses `oneLike`, the
multiplicative unit of its algebra, so `grad(y,y)` is ones per channel for
`Greal`/`Garray` and the site identity for matrix fields, where a root then
differentiates $\sum_x\operatorname{Re}\operatorname{tr}$. Implicit seeding
reaches only `grad(y,y)` and hooks written with `rootedUpstream` (`cond`,
`multi`, `apply`, `slotVar`, the NN field ops); hooks written with
`requireUpstream` raise at the root, and `Ggauge` has no default seed.

Cache key: `(output node, runtime symbolic revision)`. Only complete adjoints
are reusable, including complete intermediate adjoints. Pending contributions
commit after successful construction; failed builds cannot poison the cache.
`findGrad(input,output)` treats old revisions as misses.

Numerical updates change epochs, not symbolic gradient structure. Supported
lambda binding/internal VJP changes advance the symbolic revision; arbitrary
metadata/topology edits require a new graph. Structural VJP bodies use uncached
seeded builds; seed-dependent adjoints do not enter the public output cache.

### Per-slot partials and slot variables

A backward hook owns one input slot. Its value may depend on siblings, but its
partial must not traverse their paths again. Let $s_i=\operatorname{slotVar}(x_i)$:

$$
\operatorname{value}(s_i)=\operatorname{value}(x_i),\quad
D_{x_i}s_i=I,\quad
\bar x_i=\left.D_{s_i}^{*}f(x_1,\ldots,s_i,\ldots,x_n)[u]\right|_{s_i=x_i}.
$$

`s_i` has a distinct identity and a transparent input edge. This isolates the
inner partial while keeping dependence on `x_i` live for later derivatives.
For example, two aliased inputs to multiplication contribute separately:

$$
f(x,x)=x^2,\qquad \bar x_1=u x^\dagger,\quad \bar x_2=x^\dagger u.
$$

Replicas spell the partial with `gradSeeded(replica(slot), slot, seed)`. There is
no `secondPullback` helper: action Hessian backwards use `actionJet`; stout
backwards differentiate scalar replica scores over independent live slots
(section 12). Conditional upstreams split before backward dispatch; static-zero
branches skip inactive VJP construction and evaluation.

`slotVar` uses `valueLike`. A resolved lambda aliases as `GlambdaRef`; resolution
and structural VJPs follow its input edge, without whole-lambda cotangents.
Structural `Gmulti` aliases keep shape prototypes and refresh through that edge
without copying numerical slots.

## 5.1 Log-Jacobian Factorization (`logDetJ`)

For a declared square composition chain,

$$
v=w_0\longrightarrow w_1\longrightarrow\cdots\longrightarrow w_n=u,
\qquad
\log|\det D_vu|=\sum_{k=1}^{n}\log|\det D_{w_{k-1}}w_k|.
$$

`Gfunc.logdet(z) = (ld,via)` declares

$$
D_bz=D_{via}z\,D_b via,\qquad ld=\log|\det D_{via}z|.
$$

Every backward path to a base below `via` must pass through it; `z` and `via`
have compatible shapes. Verification uses only `iwmBackward`: it checks
squareness, dependency and absence of bypass paths. The owning op must supply the
**total** local derivative, including any re-entrant context edges.

A formula holding auxiliaries fixed must reject auxiliaries depending on `via`
(e.g. explicit stout staples, alpha or coefficients). Such dependence is allowed
only when the op proves and internalizes a frozen/triangular decomposition, as
in the action-aware subset stout step. Tests establish that local formula.
Global invertibility is an application concern; singular output remains visible.

| Chain case | Behavior |
| --- | --- |
| Repeated hook/base query | Memoize hook and sum by symbolic revision; reuse shared prefixes |
| `resetLdjCache` | Drop all entries |
| `cond` | Branchwise chain sums; zero selector adjoint; evaluate selected branch only |
| Static-zero local `ld` | Elide from sum |
| `u` is `v` | Fresh static-zero scalar |
| Missing hook, off-base leaf, `apply`/lambda node | Fail construction |

Cloning shares `Gfunc`: hooks may capture immutable configuration, but derive
`ld`, `via` and all graph values from `z.inputs`.

$$
u=f(v),\qquad S_{\rm eff}(v)=S(u)-\operatorname{logDetJ}(u,v).
$$

The fused stout update structurally carries its logdet view, so cloning preserves
the pair. A local correction without a graph flow (block5 coupling) declares no
factorization. Compared with telescoping, the generic action adds one forward
and one pullback per force: $2$ of $3K$ kernels for $K$ steps. Fuse each step's
smear/logdet kernels to recover them.

## 6. Exceptional Node Checklist

For a node beyond raw inputs plus ordinary backward, record:

- eval/backward/reachable dependencies and backward-to-input mapping;
- any square logdet factorization, derived from the current node's inputs;
- concrete recovery at erased boundaries;
- subtype/metadata cloning requirements;
- buffer mode, payload/work ownership, aliases and in-place permissions.

Shared `Gfunc` closures may capture immutable configuration and layout-only
scratch. Rebind every retained input value (`Transporter.link`, fields, etc.)
from `v.inputs` inside the forward; a captured original node is not a cloned
input. Additional structural state must be in inputs, participate in symbolic
revision invalidation, or have explicit clone handling. Symbolic VJPs expose
function/target dependencies as ordinary graph edges.

## 7. `cond`

For compatible numerical branches and a scalar/int selector, branchwise AD uses

$$
z=\begin{cases}a&c\ne0\\b&c=0,\end{cases}
\qquad
\bar c=0,\quad
\bar a=\operatorname{cond}(c,\bar z,0),\quad
\bar b=\operatorname{cond}(c,0,\bar z).
$$

Eval visits the selector and chosen branch; reachable/backward visit both.
The zero selector adjoint is the a.e. convention, not a derivative of the switch.
Selector flips reuse symbolic gradients. Results use the branch's `valueLike`
family, without producer caches. Literal branches exist only for scalar/int
values; other types require graph values. Cast erased results before typed use.
Branch selection cannot change result shape or runtime ownership.

## 8. Functional Layer: Structural Lambdas And `apply`

$$
\operatorname{Lambda}\langle A,B\rangle=A\to B.
$$

Both parameter and result prototypes are required: VJP shape depends on every
argument level and the final result cotangent. Result-only placeholders are
invalid. Higher-order AD rewrites lambda/application/VJP structure into ordinary
graph expressions.

### Lambda storage and normalization

| Value | Stored structure |
| --- | --- |
| `Glambda` | `param`, `body`, `captureParams`; captured values in `inputs` |
| `GlambdaRef` | `kind`, parameter/result prototypes, optional `binding` |
| `lrkLocal` | Local placeholder |
| `lrkProduced` | Internal symbolic function result |

Normalization moves free values into inputs and substitutes fresh capture
parameters in the body:

$$
|\mathrm{captureParams}|=|\mathrm{inputs}|,\qquad
\mathrm{captureParams}_i\mapsto\mathrm{inputs}_i.
$$

Normalization, cloning and generated-lambda construction own this invariant;
consumers may index directly. Resolved lambdas are structural values. Their
internal `lambda captures` hook exposes capture freshness only. Whole-lambda
gradients are rejected; apply builds scalar-capture gradients structurally.
`lambdaParam(paramProto,resultProto)` constructs public function parameters;
produced refs belong to internal `apply`/`cond` result shapes.

### Apply

$$
\operatorname{apply}:(A\to B)\times A\to B,\qquad
\operatorname{apply}(\lambda x.b,a)=b[x\mapsto a].
$$

Construction records the application without instantiating its body or VJP.
Instantiation occurs at eval or structural backward construction. The public
result is erased `Gvalue`; cast before typed operations. Scalar/int literals use
the **parameter** prototype, including for conditional function expressions.

Argument compatibility is demand-driven, not checked at construction. An unused
incompatible argument can be harmless; a demanded non-lambda function or
incompatible result copy fails on evaluation/differentiation of that path.

### Apply traversal and caches

Inputs always contain function and argument. VJP construction adds active value
targets as explicit edges through its build context. Generated dependencies must
exist before use; cloning appends only the extra inputs listed in the supplied
source-node map.

Eval/reachable/backward views are dynamic. Reachable declares every direct eval
dependency across both branches; backward keeps its maximal target frontier.
`cloneValues` shares a context across roots but starts a fresh expression memo
per lambda body. Only active recursive identities cross that boundary. Private
clones isolate generated constants/bound refs and ignore numerical override cuts
inside abstract bodies.

| Cache | Key and ownership |
| --- | --- |
| Ordinary apply | Selected nominal closure and symbolic revision; runtime-owned |
| Planned apply | `(applyId, selectedLambdaId, revision)`; requires eval and apply frames |
| External-boundary eval | Ordinary runtime cache |
| Plan clear | Drop private variants |

Function-side structure invalidates stale pullbacks without receiving raw lambda
cotangents. `ensureInstantiation` owns direct apply-cache operations and their
failure semantics. Apply construction must not prepare VJPs eagerly.

## 9. Structural VJP Transformation

Direct lambdas lower to generated ordinary lambdas. Unresolved function refs keep
`GvjpOf` nodes with function/target inputs and result prototypes; one apply-owned
reducer lowers them after substitution or evaluation resolves the function.

| Notation | `VjpSpec` | Result type for numerical $B$ in $f:A\to B$ |
| --- | --- | --- |
| `callVjpOf(f)` = public `vjpOf(f)` | depth 0, argument target | $A\to B^*\to A^*$ |
| `captureVjpOf(f,t)` | depth 0, value target | $A\to B^*\to T^*$ |
| `resultVjpOf(f,t)` | positive depth | Preserve returned-function arguments before the final cotangent |

$T$ is the differentiated value's type; $B^*$ comes from the **result** prototype.
`LambdaVjpTarget`
keeps argument/value target kind and concrete value together. Depth counts
arguments preserved before differentiation: applications increase it; lambda
shells decrease it; delayed VJPs and active-shell matching retain it.

For $f=\lambda x.b$, fresh $x'$ and cotangent $u$,

$$
\begin{aligned}
\operatorname{callVjpOf}(f)&=\lambda x'.\lambda u.\ D_{x'}^*b[x\mapsto x'][u],\\
\operatorname{captureVjpOf}(f,t)&=\lambda x'.\lambda u.\ D_t^*b[x\mapsto x'][u].
\end{aligned}
$$

The code spells each $D^*[u]$ as `gradSeeded(body,target,seed)`. Returned lambdas
preserve their arguments before the seed:

$$
\operatorname{callVjpOf}(\lambda x.\lambda y.b)
=\lambda x'.\lambda y'.\lambda u.\ D_{x'}^*b[x\mapsto x',y\mapsto y'][u].
$$

For $z=\operatorname{apply}(f,a)$, the argument-slot and capture contributions are

$$
\bar a=\operatorname{apply}(\operatorname{apply}(\operatorname{callVjpOf}(f),a),u),
\qquad
\bar t_{\rm capture}=\operatorname{apply}(\operatorname{apply}(\operatorname{captureVjpOf}(f,t),a),u).
$$

For $f_2=\operatorname{apply}(h,a_0)$,

$$
\operatorname{captureVjpOf}(f_2,t)
=\operatorname{apply}(\operatorname{resultVjpOf}(h,t),a_0).
$$

In $f(a_0)\cdots(a_n)$, a value target receives its direct capture contribution
and chains through inner arguments $a_0,\ldots,a_{n-1}$. The outer apply owns
$a_n$'s slot. Keep original argument nodes live for higher derivatives. Ordinary
scalar arguments use scalar reverse mode. Lambda-valued arguments have no
first-class cotangent; their delayed VJPs wait for structural substitution.

Allowed graph forms are ordinary `lambda`, `apply`, `cond`, scalar/custom
scalar/multi/gauge nodes, plus unresolved `vjpOf`/`vjpOfResult`. No `lambdaVjp`,
`GlambdaCotangent`, `lambdaCapture`, `lambdaCotangent`, `lowerLambdaCotangent`,
or persistent `structuralVjps` metadata belongs on graph values.

### Symbolic VJP substitution

For substitution $\sigma(f)=g$, ordinary clone rules give

$$
\begin{aligned}
f&\mapsto g,\\
\operatorname{vjpOf}(f)&\mapsto\operatorname{vjpOf}(g),\\
\operatorname{vjpOf}(f,t)&\mapsto\operatorname{vjpOf}(g,\sigma(t)),\\
\operatorname{vjpOfResult}(f,t)&\mapsto\operatorname{vjpOfResult}(g,\sigma(t)).
\end{aligned}
$$

Capture normalization first substitutes its capture parameter, then instantiation
substitutes the actual value. Function and target inputs remap with the rest of
the body; lambda values own no persistent VJP metadata.

### Conditionals and recursion

Keep `apply(cond(c,f,g),x)` symbolic at construction; do not eagerly expand it
into branch applications. Each VJP form $V$ distributes structurally:

$$
V(\operatorname{cond}(c,f,g))=\operatorname{cond}(c,V(f),V(g)).
$$

Capture/result forms carry the same target on each branch. Eval visits only the
selected branch.

Register a recursive VJP shell **before** transforming its body. Reuse an active
shell by nominal selected-closure identity, VJP kind, target, depth and compatible
prototypes; body-shape equality is insufficient. Normalization preserves these
shared shells rather than cloning recursion indefinitely.

The active stack lives in a passed `ApplyVjpBuildCtx`, never process-global
state. Independent builds share no active targets/memos. Completed memos include
selected lambda and capture identity; persistent caches obey symbolic revision.

## 10. Custom Graph Functions And Lambda Values

Custom `Gfunc` nodes may produce ordinary scalar/multi/gauge values inside lambda
bodies when they expose normal dependencies and forward/backward rules. Function
values use structural `lambda`, `apply`, `cond` or `lambdaParam` forms. Produced
refs and symbolic VJPs are functional-layer internals. Diagnostic names such as
`"apply"` or `"cond"` confer no structural authority.

Representative acceptance cases:

```text
vjpOf(lambda(x, x*x))
vjpOf(vjpOf(lambda(x, x*x)))
hof = lambda(f, lambda(x, Gscalar(apply(f,x))*Gscalar(apply(f,x))))
g = lambda(v, a*v+1.0)
z = apply(apply(hof,g),x)
grad(z,a)
grad(grad(z,a),a)
```

Direct lambdas reduce to generated lambdas; unresolved refs retain symbolic VJPs
until substitution or eval selects a concrete function.

## 11. `Gmulti`

| Surface | Meaning |
| --- | --- |
| Numerical slots | Concrete evaluated values |
| Lambda slots | Retained structural input values |
| `x[i]` | Symbolic selection for a fixed metadata index; erased result |
| `storedSlot(i)` | Last stored forward value; no carrier evaluation dependency |
| Shape-only carrier | Prototypes only; rejects storedSlot and generic indexing; owner supplies views |

Heterogeneous slots accumulate with each contribution's algebra. Cast known slot
types at erased boundaries. Dynamic graph indices are unsupported; put control
flow outside indexing. Read stored slots only after evaluating the carrier or a
consumer.

Fused operators document their slot contract, share real work once, and return
slot gradients as one `Gmulti`. This is not a general product-value abstraction.

## 12. Gauge Layer

Numerical kernels and conventions: [kernel reference](../../../docs/gauge_kernels.md).

| Graph value | Payload |
| --- | --- |
| `Ggauge` | One `DLatticeColorMatrixV` per direction |
| `Gfield` | One `DLatticeColorMatrixV` |
| `Gcfield` | `DLatticeComplexMatrixV[1]`; equals `Gfield` for Nc=1 |
| `Grfield` | `DLatticeRealMatrixV[1]` |
| `Grmat8` | `DLatticeRealMatrixV[8]` |

`Grmat[n]` graph storage/accumulation supports $n\in\{1,8\}$; numerical
constructors support static square sizes. Scalar fields are explicit $1\times1$
matrix sites.

Use concrete gauge/scalar/coefficient types in construction and direct casts in
operator-specific backwards. `toGvalue` copies caller gauge storage;
`gaugeSnapshot` returns an owned copy. Public mutation uses `update`/`mutateGauge`
to mark freshness. Generated zeros allocate cleared storage on demand and lose
`staticZeroLeaf` on update.

### Grad-complete basic tier

The pairing is real and unnormalized:

$$
\langle A,B\rangle=\sum_{x,\mu}\Re\operatorname{tr}(A_{x,\mu}^{\dagger}B_{x,\mu}),
\qquad df=\langle\nabla f,dA\rangle.
$$

For real matrix fields use $A^T$, with no direction sum for a single field.
A tier is grad-complete if every pullback stays inside that tier plus scalar
operations/constants; repeated AD then needs no per-order operator code.

| $z$ | Pullback for output cotangent $u$ |
| --- | --- |
| $A+B$ | $\bar A=u,\ \bar B=u$ |
| $AB$ | $\bar A=uB^\dagger,\ \bar B=A^\dagger u$ |
| $A^\dagger$ | $\bar A=u^\dagger$ |
| $\lVert A\rVert^2$ | $\bar A=2uA$ |
| $\Re\operatorname{tr}A$ | $\bar A=uI$ |
| $\Pi(A)$, $\Pi=\operatorname{projTAH}$ | $\bar A=\Pi(u)$ |
| $M_S A+(I-M_S)B$ | $\bar A=M_Su,\ \bar B=(I-M_S)u$ |

The closed tier contains the site algebra, `blendSubset`, `linkField`,
`injectLink`, `shift` and `hop`. Field/transport helpers are internal modules,
not gauge-facade exports. Site kernels are shared across storage families;
bundle subset operations remain separate.

For sitewise complex bridges (`gauge/cfield`),

| $z$ | Pullback |
| --- | --- |
| $\operatorname{tr}x$ | $\bar x=uI$ |
| $\operatorname{scale}(c,x)=cx$ | $\bar c=\operatorname{tr}(x^\dagger u),\ \bar x=c^*u$ |
| $\operatorname{dot}(x,y)=\operatorname{tr}(x^\dagger y)$ | $\bar x=u^*y,\ \bar y=ux$ |

$$
\operatorname{norm2}(\operatorname{trace}P)=\sum_x|\operatorname{tr}P_x|^2.
$$

Here $c^*$ is complex conjugation; bars denote cotangents.

For real matrices (`gauge/matrix`, exported by the facade),

$$
Y=A^{-1}B,\quad dY=A^{-1}(dB-dA\,Y),\quad V=A^{-T}u,
\qquad\bar B=V,\quad\bar A=-VY^T.
$$

| Operation | Pullback |
| --- | --- |
| $Y=A^{-1}$ | $\bar A=-Y^TuY^T$ |
| $z=\log\det A$ | $\bar A=uA^{-T}$ |

The LU kernels require nonzero leading pivots; logdet also requires
$\det A>0$, allowing negative individual pivots. Per-site real weights use
`scale`; `sum` reduces over physical sites. Transpose, contractions, scalar
functions and SU(3) bridges remain closed under further differentiation.

$$
\begin{aligned}
z=x/y &: &\bar x&=u/y,&\bar y&=-uz/y,\\
z=e^{ix} &: &\bar x&=\Im(z^*u),\\
z=\arg x &: &\bar x&=iux/|x|^2.
\end{aligned}
$$

Real/complex embedding adjoints are `re` ↔ `complex` and `im` ↔ `imaginary`.
Division requires a nonzero denominator; `ln` and differentiated `sqrt` require
positive inputs. `arg` excludes zero and the principal branch cut.

For halo offset $s$, define $G_s b(x)=b(x-s)$ and $S_s=G_s^*$. With no adjoint
flags, the stencil product and pullbacks are

$$
z=a\,G_s b,\qquad \bar a=u(G_s b)^\dagger,\qquad \bar b=S_s(a^\dagger u).
$$

`gp` adds optional operand adjoints. `lineProducts` memoizes QEX's path plan by
path key; common subproducts and derivatives are built once. Each step uses one
output field and a boundary halo, without an interior copy. Scatter clears its
output/shell each eval, places each site at its unique target, then completes
reverse exchange. Nodes rebind halo inputs each forward. General field transport
uses hop chains:

$$
H_{\mu,+}f(x)=U_\mu(x)f(x+\hat\mu),\qquad
H_{\mu,-}f(x)=U_\mu(x-\hat\mu)^\dagger f(x-\hat\mu).
$$

For the positive, unnormalized plaquette trace sum $P(U)$,

$$
J_m^P(U;h_1,\ldots,h_m)=D^m\nabla P(U)[h_1,\ldots,h_m],
\qquad \operatorname{stapleSum}=J_m^P,
$$

$$
\bar U=J_{m+1}^P(U;h_1,\ldots,h_m,u),\qquad
\bar h_i=J_m^P(U;h_1,\ldots,u,\ldots,h_m).
$$

$P$ is quartic, so $J_m^P=0$ for $m>3$. Fused orders $0\ldots3$ retain repeated
and dependent seeds; higher orders validate every seed's shape/runtime before
returning an ordinary static zero. `plaqSum` has pullback $uJ_0^P$.

Ordinary fused nodes own `PlaqWork`; clones start empty and release drops the
reference. Planned work is keyed by operation/layout/field count/order. Every
call binds current inputs/seeds and completes exchanges. `tgplaqstencil` compares
against hop chains; its U(1) wrapper covers the scalar specialization.
`tstencilmpi` checks faces and, with two split axes, corners.

| Module | Protocol |
| --- | --- |
| `gauge/types` | Storage, ownership, shapes |
| `gauge/basic_ops` | Site algebra and bundle blend/mask |
| `gauge/matfun` | Exponential jets/replicas |
| `gauge/field_ops`, `gauge/transport` | Direction fields, shifts, hops and paths |
| `gauge/cfield`, `gauge/matrix` | Scalar/matrix and SU(3) bridges |
| `gauge/stencil` | Halo moves and fused plaquette jets |
| `gauge/fused_ops` | Packed site kernels |
| `gauge/action` | Coefficients, kernel dispatch, actions and derivative bases |
| `gauge/stout` | Updates, finite logdet and grouped pullbacks |

For each fundamental loop family $k\in\{p,r,g\}$, let $L_k$ be its real trace
sum and $S_k=-L_k/N_c$. Current action jets use

$$
S(c,U)=\sum_kc_kS_k(U),\quad
B_{k,m}=D^m\nabla S_k(U)[h_1,\ldots,h_m],\quad
A_m=\sum_kc_kB_{k,m}.
$$

$$
\bar c_k=\langle u,B_{k,m}\rangle,\qquad
\bar U=A_{m+1}(c,U;h_1,\ldots,h_m,u),\qquad
\bar h_i=A_m(c,U;h_1,\ldots,u,\ldots,h_m).
$$

`actionJet` uses $B_{p,m}=-J_m^P/N_c$; rectangle/parallelogram bases use
`gradSeeded` on reference loops. Hidden bases enter eval/reachable views;
backward sees only coefficients, field and live seeds. A zero coefficient skips
basis evaluation, not its coefficient derivative.

The graph `gaugeActionDeriv` returns $\nabla S$ (the numerical wrapper negates
coefficients for the underlying negative-gradient kernel). For $H=\nabla^2S$,

$$
f_S=M_S\nabla S,\qquad (Df_S)^*b=H(M_Sb).
$$

Thus a subset Hessian pullback reaches all affected links. For $z=Hb$,
$\bar b=Hu$ and $\bar U=A_2(c,U;b,u)$; coefficient partials use fresh slots.

| Coefficient value | Selected family/tangents |
| --- | --- |
| `adjplaq == 0` | Fundamental: plaq, rect, pgm |
| `adjplaq != 0`, `rect == pgm == 0` | Adjoint plaquette: plaq, adjplaq |
| Other combinations | Unsupported |

`gaugeAction`, `gaugeActionDeriv`, `gaugeForce` and fundamental full/subset
Hessians differentiate coefficient subgraphs within these families.

Coefficient AD differentiates the active family only, without differentiating
the selector. In particular,

$$
\left.D_\beta^{\rm AD}\operatorname{gaugeAction}(\operatorname{actAdj}(\beta,a),U)
\right|_{\beta=0}=-\frac1{N_c}\sum_p\Re\operatorname{tr}P_p.
$$

`gaugeActionGraph` supplies fundamental references; `adjPlaqAction` supplies
higher adjoint derivatives. The optimized adjoint field Hessian rejects further
AD. Coefficient partials follow the selected family on every evaluation.

Stout replicas use independent slots for $W,d_s,\alpha$, coefficients and
upstreams. For update replica $R$ and finite logdet replica $\ell$,

$$
q=\left\langle\sum_j u_j,R(W,d_s,\alpha)\right\rangle+u_\ell\ell(W,d_s,\alpha),
\qquad g_i=\nabla_{s_i}q.
$$

For a gradient kernel with outputs $g_i$ and cotangents $v_i$, its backward uses

$$
r(s;v)=\sum_i\langle v_i,g_i(s)\rangle,\qquad\bar s_i=\nabla_{s_i}r.
$$

`replicaInputGrads`/`replicaSlotGrads` differentiate these paired scores over
fresh live slots. Cache slots are eval-only, with zero backward contributions.
In the action-aware step,
$d_s=M_S\nabla S(W)$, so the total field contribution is

$$
\bar W=\bar W_{\rm direct}+H(M_S\bar d_s).
$$

The coefficient overload checks Wilson-only coefficients on every evaluation:
rectangle/parallelogram staples couple active links and invalidate the independent
parity/direction factorization. Its coefficient tangents are restricted to plaq.

`tests/gauge/higher` differentiates one order beyond each built replica, with
aliased and dependent cotangents. `tgloops` covers loop families, zero coefficients
and mixed coefficient/field derivatives; `tgtoweru1`
covers U(1). These checks pin live-slot dependence rather than frozen seeds.

### Exponential and stout contracts

Let $E_m(Y;d)=D^mE(Y)[d_1,\ldots,d_m]$. `expJet` uses fused orders $1\ldots3$;
its shared pullback rule is

$$
\bar Y=E_{m+1}(Y^\dagger;u,d_1^\dagger,\ldots,d_m^\dagger),\qquad
\bar d_i=E_m(Y^\dagger;u,d_1^\dagger,\ldots,\widehat{d_i^\dagger},\ldots,d_m^\dagger).
$$

The hat omits that direction. Above order three, `expTopReplica` uses the same
slot rule recursively. For Nc=1,

$$
E_m(y;d)=e^y\prod_{i=1}^m d_i\quad(m\ge0).
$$

| Operation | Finite numerical map |
| --- | --- |
| Ordinary Nc>1 exponential/tower | Degree 4, scale 20 from `newExpParam`; static kind/order check |
| SU(3) `axexpmuly`/stout primal | Adaptive degree-12 `expAH`, scaled to $\lVert F\rVert_F^2\le1/16$ |
| SU(3) field/staple pullbacks and logdet | Degree 13 at $X/32$, five recoveries |
| SU(3) alpha pullbacks/higher update replicas | Ordinary `expDeriv`/graph exponential |
| U(1) exponential tower/stout differentials | Exact scalar formulas |

The ordinary replica spells the polynomial and squaring as

$$
Y=X/2^{20},\quad q=Y+\tfrac12Y^2+\tfrac16Y^3+\tfrac1{24}Y^4,
\qquad q\leftarrow q(2I+q)\ \text{20 times},\quad E(X)=I+q.
$$

For stout, $M=\alpha Wd_s^\dagger$, $X=-\operatorname{ad}(\Pi(M))$ and
$D=\operatorname{su3ProjectDeriv}(M)$:

$$
P_0=\sum_{k=0}^{13}\frac{(X/32)^k}{(k+1)!},\quad
P_{j+1}=P_j+2^{j-6}XP_j^2\ (j=0,\ldots,4),\quad
K=I+P_5D,\quad\ell=\log\det K.
$$

`expProjectTAHOrder=13` is shared by replica/kernel calls; scale is `expProjectTAHScale=5`.
The seed lives in the SU(3) adjoint image; cotangents may be arbitrary real
matrices. Require $\det K>0$ and nonzero leading LU pivots. U(1) uses
$\ell=\log(1+\Re M)$.

`stoutLogDetJGraph` differentiates **this finite expression** at every order.
The update primal, its first field/staple pullbacks and its higher replicas use
different finite maps (`expAH`, scaled Phi with cached `expAH`, ordinary $E$).
Higher update pullbacks therefore approximate derivatives of the fused update
and first pullback. Accuracy depends on generator norm and Jacobian conditioning;
there is no enforced norm cap or accuracy bound outside tested fixtures.

`tscaledexp`/`tgjac` use dense polynomial, converged series and differentiated
adaptive-exponential references from `scaledexpRef`. Coverage includes repeated
spectra, scalar/SIMD scaling thresholds and $\|F\|_F\le8$. The update approximation
has nine field/staple/mixed-alpha directional cases at norms $0.12,2,8$, each with
steps $10^{-3}$ and $5\cdot10^{-4}$, separate from finite-logdet derivative checks.

## 12.1 Neural Fields

`GfieldOf[F]` in `field/types` supplies the common storage contracts for a field
or sequence of fields. The neural value types support `float32` and `float64`:

| Value | Storage |
| --- | --- |
| `Greal[T]` | `GfieldOf[seq[RealField[T]]]`, one scalar QEX field per channel |
| `Gmask` | A float32 field; zero excludes a site, nonzero selects it |
| `Garray[T]` | Replicated parameter data and shape |

Computed fields allocate on evaluation and support plan buffers. Parameter
arrays allocate directly and do not use field buffer pooling. `gauge/rfield`
converts one channel to/from `Grfield` and composes the matrix operators.
NN kernels and wrappers have no checkpoint or file I/O.

The [numerical kernels](../../nn.nim) own their thread regions. Convolution is
periodic, stride-one cross-correlation; `convParams` constructs dense odd kernels
with weights ordered `[output,input,spatial...]`, last spatial axis fastest.
Forward subset/mask operations preserve the destination complement. Pointwise
kernels permit in-place use; convolution requires disjoint source/destination
fields. `convVjp` clears its destination and accumulates reverse halos;
pointwise pullbacks zero complements unless `passthrough` is requested.

NN parameter pullbacks are graph expressions. For convolution $C(x,w)$ and
$B(x,b)_{oit}=\sum_s b_o(s)x_i(s+t)$,

$$
\partial_x\langle b,C(x,w)\rangle=C^T(b,w),\qquad
\partial_w\langle b,C(x,w)\rangle=B(x,b),
$$
$$
\partial_x\langle a,B(x,b)\rangle=C^T(b,a),\qquad
\partial_b\langle a,B(x,b)\rangle=C(x,a).
$$

`channelSum` and `broadcast` are adjoints for bias and channel-scale pullbacks.
Parameter gradients sum owned sites and reduce across spatial ranks once;
array dot products act on the replicated result without another rank reduction.
These rules permit one parameter derivative before or after repeated input
pullbacks. Native tests cover input orders through four and one parameter
derivative combined with input orders zero through four. Masks are discrete;
`clipMin` has slopes zero, one-half and one below, at and above its threshold.

## 13. `hmcgauge`

`hmcgauge` owns trajectory construction, sampling, integration and training.
`flowAction(gc, map)` caches each transformed graph by its input node and exposes
`flow` and `action` callables, with $S_{\rm eff}(V)=S(f(V))-\log\det J_f(V)$.
Stout flows use this constructor. `measure2du1.runFlowHmc` owns
their physical-field plans, proposal measurements and sampling statistics.
Parameters remain paired with gradient expressions. `stoutAction` and
`smearedField` accept graph rho; `StoutAction.rho` exposes it for mixed derivatives.
Float overloads create a scalar in the input runtime.

$$
H=S+T,\qquad\Delta H=H_f-H_0,\qquad a=e^{-\Delta H},\qquad
L=-\min(1,a)\,(n_{\rm steps}\,dt)^2.
$$

`runHmc` uses one active proposal plan per phase/root set. Roots include initial
and final $H,S,T$, each force's scalar RMS/min/max, final gauge, and requested
reverse momentum, view, loss and gradients. Early scalar reductions avoid
retaining force history solely for diagnostics.

| Callback/data | Ownership/order |
| --- | --- |
| `Proposal(gauge,view,dH,acc,loss,gradients)` | Gauge/view borrowed through callback; scalars/gradients copied |
| `trainStep` | Consume copied gradients; own no execution plan |
| Accepted state | Snapshot final gauge before reverse checks or parameter callbacks |
| Reverse plan | Use saved final values; `finally` restores both original leaves and clears it |
| Commit/measurement | Commit earlier accepted snapshot, then evaluate committed measurements |

Take `gaugeSnapshot` for values retained beyond the callback. `pghmc` omits
training expressions when `trajsTrain == 0`; physical/logdet measurements use a separate
joint plan. Force-gradient integrators accept the default coefficient tuple or
a complete explicit tuple, not partial positional completion. Validate
`IntegratorCoeffs` when constructing the run spec.

`2MNp` is the momentum-first 2MN schedule; one procedure integrates both
orderings with the same minimal-norm default lambda `0.1931833275037836`, and
kicks apply the negated coefficient. Built with `trace = true`, every
integrator records its force, kick, drift and force-gradient shift events as
graph nodes; `tintegratoru1` replays them against the analytic U(1) dynamics,
and plans root the event values to collect a complete trace.

## 14. Validation And Benchmarks

See the [validation and benchmark protocol](graph_validation.md)
for test execution, counters/settings and required direct/planned fingerprint
comparisons.
