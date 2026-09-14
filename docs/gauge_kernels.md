# Numerical gauge kernels

## Implemented

### Fields and storage

Numerical gauge bundles contain one field per direction, with a common layout.
Lattice extents are even. Inputs, seeds and outputs use compatible layouts.

The numerical matrix types are defined in `physics/qcdTypes.nim`:

| Site values | Lattice field |
| --- | --- |
| `DRealMatrixV[n]` | `DLatticeRealMatrixV[n]` |
| `DComplexMatrixV[n]` | `DLatticeComplexMatrixV[n]` |

Both families contain square `n x n` matrices. `[1]` explicitly denotes scalar
matrix sites. `RealMatrix(l,n)` and `ColorMatrix(l,n)` construct fields for the
layout's SIMD width. `toScalar` and `toMatrix` convert between scalar sites and
1x1 matrix sites.

`newShape(l,T)` and `newShape(f)` create field descriptors containing layout,
element size and allocation shape. They allocate no numerical storage. This lets
execution planners assign an existing buffer to a field without first allocating
another one. A descriptor made from a `FieldArray` view describes one field.
Bind or allocate storage before indexing a descriptor. Replacing a descriptor
preserves other references to its former field and storage.

`newShifter(f,dir,len,dest=d)` and `newTransporter(u,f,dir,len,dest=d)` use the
supplied output field. Its layout matches `f`, and its storage is disjoint from
application operands. Rebinding `field` preserves these requirements. The
transporter owns its communication buffers. Test: `tfieldstorage`.

### Matrix-field operations

`field/matrixFields.nim` implements local matrix and scalar operations. Call
these kernels inside `threads`; each traverses its output's site partition.

`siteTrace`, `siteNorm2`, `siteRedot` and `siteDot` return 1x1 fields. `sum` of a
real 1x1 field reduces physical sites and ranks without volume normalization;
every thread participates. `scale` supports real or complex scalar matrix sites.

`solve`, `inverse` and real `logDet` use unpivoted LU. Factors and right-hand
sides are matching square matrices, with nonzero leading pivots. `solve` and
`inverse` permit input/output aliases. For a positive determinant,

```text
A = L U,   ln det A = sum_k ln |U_kk|
```

Individual pivots can be negative. The logarithmic sum avoids constructing a
potentially overflowing or underflowing determinant.

`blendSubset` selects its candidate on a contiguous outer-site interval and
its other input elsewhere. `maskSubset` selects its input on that interval and
zero elsewhere. Selection operates on whole SIMD sites and permits input/output
aliasing.

Scalar functions follow their scalar domains: positive inputs for `ln`,
nonnegative inputs for `sqrt`, and nonzero divisors. `arg` uses the principal
angle in `[-pi,pi]`; its derivatives exclude zero and the branch cut.
Test: `tmatrixFields`.

### Paths, halos and HYP

Directed links use signed, one-based directions. `gaugeUtils.plan` records shared
products in dependency order and one output per requested path. `optimalPairs`
shares adjacent products and their adjoints using

```text
(A B)^dag = B^dag A^dag
```

`gaugeProd(...,origin=true)` shifts outputs to their starting sites. Unshifted,
unadjointed outputs can alias inputs or shared products. Copy a result when it
needs independently mutable storage. Product entries are released after their
last planned use. Test: `tpathplan`.

Halo layouts preserve physical field indices and group halo sites by parity.
Neighbor tables and gather maps translate these indices to periodic coordinates.
Geometry constructors cache immutable layouts and maps and are called outside
`threads`. Halo storage is defined only at entries populated by its gather map.
Reverse exchange accumulates mapped halo entries into physical fields.
Test: `thalo`.

The HYP reverse pass projects intermediates within their forward construction
domains. Neighbor validity and physical cotangent reach determine those domains.
`thypsmearhalo` poisons padding and compares physical results with the independent
shift implementation, including checks on stored floating-point components.

### Gauge actions and workspaces

`GaugeActionCoeffs` defines two families:

| Family | Coefficients | Entry points |
| --- | --- | --- |
| Fundamental loops | `plaq`, `rect`, `pgm` | `gaugeAction1/2/3`, `gaugeActionDeriv`, `gaugeDeriv2`, `gaugeForce*`, fundamental Hessians |
| Plaquette and adjoint plaquette | `plaq`, `adjplaq` | `actionA`, `gaugeADeriv`, `forceA` |

`action(c,g)` and `force(c,g,f)` select the family by nonzero `adjplaq`.
`actionA` uses the additive constant `V d(d-1)/2 * (plaq+adjplaq)` relative to
its unshifted loop sum. At zero `adjplaq`, the dispatcher selects the fundamental
action.

For `sg < nu < mu`, the loop paths and fundamental action are

```text
plaq: [mu,nu,-mu,-nu]
rect: [mu,nu,nu,-mu,-nu,-nu], [mu,mu,nu,-mu,-mu,-nu]
pgm:  [mu,nu,sg,-mu,-nu,-sg], [mu,sg,nu,-mu,-sg,-nu],
      [nu,mu,sg,-nu,-mu,-sg], [mu,-nu,sg,-mu,nu,-sg]
S = -(1/Nc) sum (c_plaq Re tr P + c_rect Re tr R + c_pgm Re tr G)
```

At identity links,
`S = -V (c_plaq d(d-1)/2 + c_rect d(d-1) + c_pgm 4 binom(d,3))`.
Every occurrence of a link contributes, including repeated physical links after
periodic wrapping. Ambient derivatives use the real Frobenius pairing:

```text
gaugeActionDeriv(c,g,f)      : f = -grad S, or accumulation when requested
gaugeDeriv2(c,g,f)          : f += grad S
gaugeDerivDeriv2(c,g,h,f)   : f += Hessian(S)[h]
```

Call complete gauge operations outside `threads`. `newLoopWork(g[0])` creates
an owner for one layout and field type. Passing `work=w` retains loop halo
buffers, a temporary gauge for improved subset derivatives, and the transporter
bundles used by full Hessians. Each evaluation refreshes inputs and releases
input references on return. A workspace serves one evaluation at a time;
independent workspaces have independent storage. Communicator sequencing and
geometry-initialization requirements also apply.

A subset gradient selects one parity and direction. Its Hessian pullback applies
the full Hessian to a masked seed and reaches every affected link. `Subset`
overwrites, `SubsetAdd` accumulates, `SubsetAddBase` adds on selected links and
sets `base+H` elsewhere, and `SubsetSum` sums seeds on the selected parity.
`gaugeDeriv2SubsetWork` with `clear=false` preserves the complement; `clear=true`
zeros it. Rectangle/parallelogram subsets use the workspace's temporary gauge.

Warm full and improved subset Hessians reuse the buffers managed by QEX's raw
allocator. Action kernels allocate their staple scratch, and generic halo gather
allocates temporary message sequences. Tests: `tgaugeloops`, `tgaugeaction`.

### Fused loop kernels

`loopAction(c,g,work=w)` evaluates the combined fundamental action in one stencil
pass. `loopDeriv` computes its gradient and mixed directional derivatives:

```nim
let w = newLoopWork(g[0])
let s = c.loopAction(g,work=w)
c.loopDeriv(g,f,work=w)              # f = grad S
c.loopDeriv(g,[h],f,work=w)          # f = Hessian(S)[h]
c.loopDeriv(g,[h,k],f,work=w)        # f = D^2 grad S[h,k]
```

Seed arrays have a static length. The derivative overwrites `f`. Its degree is
three for plaquette gradients and five when rectangle or parallelogram terms
are active; derivatives above that degree vanish. Outputs are disjoint from
seeds and from the gauge whenever the gauge values are read. At the highest
nonzero derivative, output storage can reuse the gauge.

The general kernels share partial matrix products and their adjoints across
active loops. Each matrix symbol identifies a link direction and relative site
offset. `gaugeUtils.plan(shifts=false)` factors these symbols without introducing
field shifts. Mixed derivative coefficients obey

```text
(A B)_S = sum_{T subset S} A_T B_(S \ T)
```

Only coefficients reachable from the requested outputs are evaluated. Zero
coefficient families are excluded from the paths, halo maps and product plans;
plaquettes are also excluded above derivative order three. `LoopWork` retains
aligned scratch per thread, reused at each site, with no intermediate lattice
fields. The highest derivative uses only seed loads.

`plaqSum` and `stapleSum` provide unnormalized plaquette operations:

```text
P(g) = sum_{x,mu>nu} Re tr U_mu(x) U_nu(x+mu) U_mu(x+nu)^dag U_nu(x)^dag
plaqSum(g) = P(g)
stapleSum(g,ds,f) : f = D^(ds.len) grad P(g)
```

`PlaqWork` owns buffers for a fixed layout, derivative order and action mode.
Every evaluation rebinds inputs; the workspace retains those references until
rebound or released. Its outputs are disjoint from seeds and, below order three,
the gauge. `newOneOf` creates independent workspace storage.

The dedicated staple kernels and production six-link helpers use
`maths/matrixFunctions.productJet`, which propagates
coefficients of a product of affine factors:

```text
P_S <- P_S A_0 + sum_{d in S} P_{S \ {d}} A_d
```

Static factor count and derivative order prune coefficients that cannot reach
the requested derivative. Descending subset updates reuse prefix storage.
Three factors require 2, 5, 9, 9 matrix multiplications at orders 0, 1, 2, 3.
Tests: `tproductjet`, `tplaqstencil`, `tloopfused`.

### SU(3) projected-exponential differentials

For `tr(T_a T_b) = -delta_ab/2`, `groupOps` provides the real 8x8 maps

```text
X_ab(M) = -2 Re tr([T_a,T_b] M) = -ad(projectTAH(M))_ab     su3AdNeg
D_ab(M) = -2 Re tr(T_a T_b M)                                 su3ProjectDeriv
```

Under `sum_ab H_ab X_ab` and `Re tr(G^dag dM)`, their adjoints are

```text
X*(H) =  2 sum_ab H_ab [T_a,T_b]     su3AdNegAdj
D*(H) = -2 sum_ab H_ab T_b T_a       su3ProjectDerivAdj
```

For `Z=exp(F) U`, `F=projectTAH(M)` and `R=ad(F)`, the Lie differential gives

```text
phi1(z) = sum_{k>=0} z^k/(k+1)!
A = exp(R) + phi1(R) D = exp(R) K
K = I + phi1(-R) D,   det A = det K
```

Polynomial evaluation handles the zero modes of `R`. With static scale `s`,

```text
X = -R, P_0 = P_order(X/2^s)
P_order(Y) = sum_{k=0..order} Y^k/(k+1)!
P_{j+1} = P_j + 2^(j-s-1) X P_j^2, j < s
```

`scale` defaults to zero; `expProjectTAHScale=5` selects five recoveries. The
seed degree defaults to 13. The matrix and scale-zero vector kernels share
even/odd coefficient construction. The SU(3) adjoint identity is

```text
X^8 = -s1 X^6 - s1^2 X^4/4 - s3 X^2
s1 = 3 |F|^2,   s3 = |F|^6/2 - 27 (Im det F)^2
```

Structured perturbations stay in the SU(3) adjoint image; cotangents can be
arbitrary real matrices. Recovery values and squares used by the reverse pass
are stored on the stack. For Frobenius cotangent `U` and `Q=P^2`, one reverse
step is

```text
X_bar += c U Q^T, B = c X^T U, P_bar = U + B P^T + P^T B
```

The pullback applies the selected polynomial differential to its cotangent.
Logdet gradients differentiate `ln det(I+Phi D)`. These are approximations to
the analytic Lie differential of the adaptive degree-12 `expAH` update. Their
accuracy depends on the input norm and Jacobian conditioning. Logdet evaluation
requires positive determinant and nonzero leading LU pivots. U(1) uses exact
scalar formulas.

`tscaledexp` and the graph Jacobian tests use `tests/base/scaledexpRef.nim`
to construct references at runtime: dense polynomial powers
and forward jets, a converged analytic exp/Phi series with scaling, pivoted
elimination, and a differentiated degree-12 exponential with its primal branch
held fixed. Cases cover zero and tiny inputs, repeated and rotated spectra,
and scalar/SIMD scaling thresholds through `||F||_F=8`. `texpcontract` compares
with numerical differentiation of the update; `tsu3` checks generator conventions
and structured operations.

### Tests and benchmarks

Test lattice parameters use `letParam`. Defaults based on
`latticeFromLocalLattice` preserve the local test volume as the rank count grows;
`-lat:` supplies an explicit global lattice. Run numerical tests with `make tests`.
Experimental test generation is available through `make tests experimental`.
The generated scripts support the repository's MPI test runner.

`benchGaugeActions`, `benchPlaq` and `benchExpProject` use the standard `tic`/`toc`
profiler and `getElapsedTime`. Parameters select lattice, repetitions and trials;
site benchmarks also select differential operation, scale, precision and norm.
`benchPlaq` compares production and fused actions, projected forces and ambient
Hessian applications with Symanzik plaquette/rectangle coefficients and `pgm`
set to zero or 0.07. Workspaces are warmed before timing; result comparisons
run outside the timed regions.

## Future work

- Compare the general loop kernels on production machines before selecting
  additional production dispatches.
- Add reusable message buffers to the generic halo gather workspace.
- Add vector or rectangular LU right-hand sides when an operator needs them.
- Extend adjoint-plaquette Hessian kernels and validate their coefficient family.
