# Numerical gauge kernels

## Implemented

### Fields and storage

A gauge is one field per direction, with a common layout and even lattice
extents. Inputs, seeds and outputs have compatible layouts.

| Site type (`physics/qcdTypes.nim`) | Lattice field | Shape |
| --- | --- | --- |
| `DRealMatrixV[n]` | `DLatticeRealMatrixV[n]` | $n\times n$ real |
| `DComplexMatrixV[n]` | `DLatticeComplexMatrixV[n]` | $n\times n$ complex |

`[1]` means a scalar matrix site. `RealMatrix(l,n)`/`ColorMatrix(l,n)` use the
layout's SIMD width; `toScalar`/`toMatrix` convert scalar-site representations.

| Storage operation | Contract |
| --- | --- |
| `newShape(l,T)`, `newShape(f)` | Layout/type/allocation descriptors only; bind or allocate before indexing |
| Shape from `FieldArray` view | Describes one field |
| `release(f)`, `release(array)` | Free owned storage immediately; aliases become invalid. Release a `FieldArray` as a whole. |
| Replace descriptor | Other references retain the former field/storage |
| `newShifter(...,dest=d)`, `newTransporter(...,dest=d)` | Matching layout; destination disjoint from operands; transporter owns communication buffers |
| Rebind `field` | Preserve the same layout/disjointness rules |

Test: `tfieldstorage`.

### Matrix-field operations

`field/matrixFields.nim` kernels run inside `threads`, on each thread's site
partition. Local contractions return $1\times1$ fields:

$$
\begin{aligned}
\operatorname{siteTrace}(A)_x&=\operatorname{tr}A_x,&
\operatorname{siteDot}(A,B)_x&=\operatorname{tr}(A_x^\dagger B_x),\\
\operatorname{siteRedot}(A,B)_x&=\Re\operatorname{tr}(A_x^\dagger B_x),&
\operatorname{siteNorm2}(A)_x&=\operatorname{tr}(A_x^\dagger A_x),\\
\operatorname{sum}(r)&=\sum_{\text{physical }x}r_x,&
\operatorname{scale}(c,A)_x&=c_xA_x.
\end{aligned}
$$

`sum` is global across ranks, unnormalized, and requires every thread.
`scale` accepts real or complex scalar matrices.

$$
A=LU,\qquad\operatorname{solve}(A,B)=A^{-1}B,\qquad
\log\det A=\sum_k\log|U_{kk}|\quad(\det A>0).
$$

LU is unpivoted: factors/RHS are matching square matrices and leading pivots
must be nonzero. Individual pivots may be negative. The log sum avoids an
intermediate determinant product. Solve/inverse permit input/output aliases.

For the mask $M_S$ of a contiguous outer-site interval,

$$
\operatorname{blendSubset}(A,B)=M_SA+(I-M_S)B,\qquad
\operatorname{maskSubset}(A)=M_SA.
$$

Selection uses whole SIMD sites and permits aliases. Scalar domains:
$\ln x:x>0$; $\sqrt x:x\ge0$; division: denominator nonzero;
$\arg z\in[-\pi,\pi]$, with derivatives excluding zero and the branch cut.
Test: `tmatrixFields`.

### Fused gauge updates

`field/matrixFields.nim` kernels run inside `threads` on each thread's site
partition, on one field or on every direction of a gauge:

$$
\operatorname{axexp}(a,X)_x=e^{aX_x},\qquad
\operatorname{axexpmuly}(a,X,Y)_x=e^{aX_x}Y_x,\qquad
\operatorname{contractProjectTAH}(X,Y)_x=\operatorname{projectTAH}(X_xY_x^\dagger).
$$

`axexp` and `axexpmuly` accept any $X$ and use the general `exp` of
`maths/matrixFunctions`: the degree 4 polynomial of $e^{2^{-20}aX}-1$ and 20
doublings. `axexpah` and `axexpahmuly` compute the same products with `expAH`,
which requires $X$ anti-Hermitian (traceless for SU(N)): scale to
$\lVert aX\rVert_F^2\le 1/16$, the degree 12 Taylor polynomial of
$r=e^{aX/2^k}-1$, then $k$ squarings in expm1 form, $r\leftarrow r(r+2)$,
$E=1+r$, so each squaring rounds relative to $\lVert r\rVert$ instead of
$\lVert E\rVert=1$. On traceless anti-Hermitian input `expAH` is at least as
accurate as `exp` at every norm tested (`texpaccuracy`, against a
double-double reference: at most 0.4, 1.0 and 1.8 eps for $\lVert X\rVert_F$
in 0.25..0.5, 0.5..1 and 1..2 against 1.3, 2.9 and 5.1 eps for `exp`, N = 3,
float64) and takes 0.4 to 0.5 of its time. `expAH` evaluates a matrix whose
entries are wider than one SIMD register of the build's intrinsics (8 float64
lanes with AVX512, 4 with AVX, 2 with SSE; VLEN lanes without intrinsics) in
slices of one register with the scale of the whole matrix,
about twice as fast at 16 lanes on AVX; the slices may round differently from
one evaluation of the whole matrix.
Real entries are sliced the same way. On 1x1 complex fields the `axexpah`
kernels evaluate $\cos t+i\sin t$ for $aX=it$; a real 1x1 anti-Hermitian
matrix is 0 and takes `expAH` = exp(0) = 1. Lattice code whose argument is
anti-Hermitian by construction uses it: the HMC momentum updates through
`axexpahmuly`, the gradient flow (`wflow`), gauge fixing (`gaugefix`), stout
smearing and its inverse (`stoutsmear`) and `warmSU`.
The optional `expax` receives E from the same site computation; omitting it
removes its store. E-only has no Y argument. The destination may alias an
operand. An optional E destination must be distinct from the product
destination and Y. Subset overloads leave the complement untouched.
`contractProjectTAH(g,f)` with two arguments overwrites `f` in its own threads
region. `stoutsmear.stoutStepKernel` composes them for the subset stout step.
The MD gauge updates of the HMC programs call `axexpahmuly` and the force
projections call `contractProjectTAH`. Tests: `tmatrixFields`, `texpaccuracy`, `taxexpwide`.

Whole-field and parity loops can assign a site to different threads. Join or
insert a barrier between dependent phases that change partitions, including
whole-field initialization before a subset write. `tmatrixFields` checks output
masks, aliases, complements and allocation reuse at one and six threads. A local
E product and a stored/reloaded E product may differ by ordinary FMA roundoff;
the product comparison uses a normalized norm bound while E is checked exactly.

### Paths, halos and HYP

Paths use signed directions $\pm(\mu+1)$. `gaugeUtils.plan` orders shared
products by dependency; `optimalPairs` uses

$$
(AB)^\dagger=B^\dagger A^\dagger.
$$

`gaugeProd(...,origin=true)` shifts outputs to their starting sites. Unshifted,
unadjointed results may alias inputs/shared products; copy before independent
mutation. Each product is freed after its last consumer unless an output
aliases it, and the shift fields on return. Test: `tpathplan`.

Halo layouts preserve physical indices and parity; maps encode periodic
coordinates. Construct/cache geometry outside `threads`. The halo buffer is
allocated uninitialized: only gathered entries are defined, and the reverse
exchange reads every cell its map covers before accumulating into physical
fields, so zero or fill the halo first, as `hypsmear2` does. Test: `thalo`.

Retained gather send/receive buffers belong to the communicator. Gathers through
one communicator are serialized; each call completes its messages before the
buffers can be reused. The buffers grow on demand and are released with their
owner. Finalizing the default communicator releases them explicitly. These
message buffers remain outside `LoopWork.bytes`.

HYP pullbacks project intermediates only within their forward construction
domains, respecting neighbor validity and cotangent reach. `thypsmearhalo`
compares stored floating-point components against independent shifts with
poisoned padding.

### Gauge actions and workspaces

`newLoopWork(f, keep=keepStaples)` opts into tagged staple reuse. The default
frees temporary staple storage before returning; derivatives and rectangle
actions reuse a bounded batch of 12 fields, or 18 for a derivative with a
rectangle term, in four dimensions.

| Family | Coefficients | Numerical entry points |
| --- | --- | --- |
| Fundamental | `plaq`, `rect`, `pgm` | `gaugeAction1/2/3`, `gaugeActionDeriv`, `gaugeDeriv2`, `gaugeForce*`, Hessians |
| Plaquette + adjoint | `plaq`, `adjplaq` | `actionA`, `gaugeADeriv`, `forceA` |

`action`/`force` dispatch on `adjplaq`: zero selects the fundamental family,
nonzero the adjoint family. Fundamental kernels reject `adjplaq != 0` and the
adjoint kernels reject `rect`/`pgm` with `ValueError` (`requireFundamental`,
`requireAdjoint`), so mixed coefficient sets have no entry point.
For $\sigma<\nu<\mu$, paths are

```text
plaq: [mu,nu,-mu,-nu]
rect: [mu,nu,nu,-mu,-nu,-nu], [mu,mu,nu,-mu,-mu,-nu]
pgm:  [mu,nu,sg,-mu,-nu,-sg], [mu,sg,nu,-mu,-sg,-nu],
      [nu,mu,sg,-nu,-mu,-sg], [mu,-nu,sg,-mu,nu,-sg]
```

With volume $V$, dimension $d$, and color count $N_c$,

$$
S=-\frac1{N_c}\left(c_p\sum_P\Re\operatorname{tr}P+
 c_r\sum_R\Re\operatorname{tr}R+c_g\sum_G\Re\operatorname{tr}G\right),
$$

$$
S(I)=-V\left[c_p\binom d2+c_r d(d-1)+4c_g\binom d3\right].
$$

`actionA` adds $V\binom d2(c_p+c_a)$ to its unshifted trace sum.
Every link occurrence contributes, including periodic repeats. Under
$\langle A,B\rangle=\sum\Re\operatorname{tr}(A^\dagger B)$:

| Numerical kernel | Output |
| --- | --- |
| `gaugeActionDeriv(c,g,f)` | $f=-\nabla S$, or accumulation when requested |
| `gaugeDeriv2(c,g,f)` | $f\mathrel{+}=\nabla S$ |
| `gaugeDerivDeriv2(c,g,h,f)` | $f\mathrel{+}=\nabla^2S\,h$ |
| `gaugeForce(c,g,f)`, `force` | $f_\mu=\operatorname{projectTAH}\big(U_\mu(-\nabla S)_\mu^\dagger\big)$: `contractProjectTAH` after the derivative |

Call complete gauge operations outside `threads`. `newLoopWork(g[0])` owns work
for one layout/type: halos, improved-subset temporary gauge and full-Hessian
transporters. Calls refresh/release input references; a workspace serves one
evaluation at a time. Distinct workspaces own distinct storage. Communicator
sequencing and geometry initialization requirements still apply.

For a parity/direction mask $M_S$,

$$
g_S=M_S\nabla S,\qquad (Dg_S)^*b=\nabla^2S(M_Sb).
$$

The Hessian reaches every affected link. `gaugeDerivDeriv2Subset` overwrites;
`gaugeDerivDeriv2SubsetAdd` accumulates; `gaugeDerivDeriv2SubsetAddBase` adds on
selected links and sets base plus Hessian elsewhere; `gaugeDerivDeriv2SubsetSum`
combines selected-parity seeds. The first derivative `gaugeDeriv2SubsetWork`
takes caller-built `sd`, `sf`, `sb` (a direction shifter and parity shift
buffers, as `gaugeDeriv2Subset` builds them), preserves the complement for
`clear=false` and zeros it for `clear=true`. Rectangle/parallelogram subsets
take the full-derivative route through the workspace temporary gauge and ignore
those buffers.

Warm Hessians reuse QEX raw buffers. Action staple scratch and generic halo
message sequences still allocate. Tests: `tgaugeloops`, `tgaugeaction`.

### Fused loop kernels

`loopAction(c,g,work=w)` evaluates the fundamental action. For a seed list $h$,

$$
\operatorname{loopDeriv}(c,g,h)=D^{|h|}\nabla S(g)[h_1,\ldots,h_{|h|}].
$$

```nim
let w = newLoopWork(g[0])
let s = c.loopAction(g,work=w)
c.loopDeriv(g,f,work=w)        # grad S
c.loopDeriv(g,[h],f,work=w)    # Hessian(S) h
c.loopDeriv(g,[h,k],f,work=w)  # D^2 grad S[h,k]
```

Seed count is static. The output overwrites `f`; it is disjoint from seeds and
from any gauge values still read. At the highest nonzero order, only seeds are
loaded and output may reuse gauge storage.

$$
\deg(\nabla S)=
\begin{cases}3&\text{plaquette only},\\5&\text{rectangle/parallelogram active},\end{cases}
\qquad D^m\nabla S=0\quad(m>\deg\nabla S).
$$

Mixed product coefficients satisfy

$$
(AB)_S=\sum_{T\subseteq S}A_TB_{S\setminus T}.
$$

`gaugeUtils.plan(shifts=false)` factors matrix symbols by link direction/offset,
without field shifts. Only coefficients reaching requested outputs are evaluated.
Zero-weight families are absent from paths/halos/plans; plaquettes vanish above
order three. `loopAction`/`loopDeriv` reuse aligned per-thread site scratch and
no intermediate lattice fields; the rectangle/parallelogram subset routes, first
derivative and subset Hessian, use the workspace's scratch gauge, and the full
Hessian uses its transporters.

For plaquettes,

$$
P(U)=\sum_{x,\mu>\nu}\Re\operatorname{tr}
 [U_\mu(x)U_\nu(x+\hat\mu)U_\mu(x+\hat\nu)^\dagger U_\nu(x)^\dagger],
$$

$$
\operatorname{plaqSum}(U)=P(U),\qquad
\operatorname{stapleSum}(U,h)=D^{|h|}\nabla P(U)[h].
$$

`PlaqWork` fixes layout/order/action mode and rebinds inputs each call. It keeps
those references until the next call rebinds them; there is no release, so drop
the workspace to free them. Outputs are disjoint from seeds and, below order
three, the gauge. `newOneOf` allocates independent work.

`productJet` propagates products of affine factors:

$$
P_S\leftarrow P_S A_0+\sum_{d\in S}P_{S\setminus\{d\}}A_d.
$$

Descending subsets reuse prefix storage; static factor count/order prune
unreachable coefficients. Three factors cost $2,5,9,9$ matrix multiplications
at orders $0,1,2,3$. Tests: `tproductjet`, `tplaqstencil`, `tloopfused`.

### SU(3) projected-exponential differentials

Use anti-Hermitian generators with $\operatorname{tr}(T_aT_b)=-\delta_{ab}/2$.
`groupOps` defines

$$
X_{ab}(M)=-2\Re\operatorname{tr}([T_a,T_b]M)
 =-\operatorname{ad}(\operatorname{projectTAH}M)_{ab},
\qquad D_{ab}(M)=-2\Re\operatorname{tr}(T_aT_bM).
$$

Under $\sum_{ab}H_{ab}\,dX_{ab}=\Re\operatorname{tr}(G^\dagger dM)$,

$$
X^*(H)=2\sum_{ab}H_{ab}[T_a,T_b],\qquad
D^*(H)=-2\sum_{ab}H_{ab}T_bT_a.
$$

API: `su3AdNeg`, `su3ProjectDeriv`, `su3AdNegAdj`, `su3ProjectDerivAdj`.
Cotangents $H$ may be arbitrary real matrices; structured perturbations remain
in the SU(3) adjoint image.

For $M=UC$ with fixed $C$, let $Z=e^F U$, $F=\operatorname{projectTAH}(M)$
and $R=\operatorname{ad}F$. The Jacobian in left Lie coordinates is

$$
\phi_1(z)=\sum_{k\ge0}\frac{z^k}{(k+1)!},\qquad
A=e^R+\phi_1(R)D=e^R K,\qquad K=I+\phi_1(-R)D,\qquad\det A=\det K.
$$

Polynomial evaluation handles zero modes. With $X=-R$, static scale $s$ and
seed degree $n$,

$$
P_0=\sum_{k=0}^{n}\frac{(X/2^s)^k}{(k+1)!},\qquad
P_{j+1}=P_j+c_jXP_j^2,\qquad c_j=2^{j-s-1},\quad0\le j<s.
$$

Kernel defaults: $s=0$, $n=13$; `expProjectTAHScale=5` selects five recoveries.
SU(3) logdet-gradient kernels require positive odd $n$. Matrix/scale-zero vector
kernels share even/odd coefficients. SU(3) reduction:

$$
X^8=-s_1X^6-\frac{s_1^2}{4}X^4-s_3X^2,\qquad
s_1=3\|F\|_F^2,\quad s_3=\frac{\|F\|_F^6}{2}-27(\Im\det F)^2.
$$

Recovery values/squares are stack-resident. For cotangent $U$, $Q=P^2$ and
$P'=P+cXP^2$, one reverse step is

$$
\bar X\mathrel{+}=cUQ^T,\qquad B=cX^TU,\qquad
\bar P=U+BP^T+P^TB.
$$

The logdet pullback differentiates its selected finite expression:

$$
\ell=\log\det(I+PD),\quad d\ell=\operatorname{tr}(K^{-1}dK),\quad
\bar P=K^{-T}D^T,\quad\bar D=P^TK^{-T},\quad K=I+PD.
$$

Require positive determinant and nonzero leading LU pivots. Phi/logdet kernels
approximate the analytic Lie differential of adaptive degree-12 `expAH`;
accuracy depends on generator norm and conditioning. U(1) uses exact scalar
formulas.

`scaledexpRef` builds dense polynomial/forward-jet, converged exp/Phi-series,
pivoted-elimination and differentiated degree-12 references (primal branch fixed).
`tscaledexp` and graph Jacobian tests cover zero/tiny inputs, repeated/rotated
spectra and scalar/SIMD scaling thresholds through $\|F\|_F=8$.
`texpcontract` checks numerical update derivatives; `tsu3` checks generators and
structured operations.

### Tests and benchmarks

`letParam` supplies lattice parameters. `latticeFromLocalLattice` derives global
volume from the local lattice/rank count; `-lat:` is explicitly global.
Test generation and the MPI launchers are described in the
[graph validation protocol](../src/experimental/graph/graph_validation.md).

`benchGaugeActions`, `benchPlaq` and `benchExpProject` use `tic`/`toc` and
`getElapsedTime`; parameters select lattice/repetitions/trials and, for site
benchmarks, operation/scale/precision/norm. `benchPlaq` compares production/fused
actions, forces and Hessians for Symanzik coefficients with `pgm=0` or `0.07`.
Warm workspaces first; validate results outside timing.

## Future work

- Measure general loop kernels on production hardware before further dispatch changes.
- Add vector/rectangular LU right-hand sides when needed.
- Extend and validate adjoint-plaquette Hessian kernels.
