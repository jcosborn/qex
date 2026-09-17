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

### Paths, halos and HYP

Paths use signed directions $\pm(\mu+1)$. `gaugeUtils.plan` orders shared
products by dependency; `optimalPairs` uses

$$
(AB)^\dagger=B^\dagger A^\dagger.
$$

`gaugeProd(...,origin=true)` shifts outputs to their starting sites. Unshifted,
unadjointed results may alias inputs/shared products; copy before independent
mutation. Release product entries after their last consumer. Test: `tpathplan`.

Halo layouts preserve physical indices and parity; maps encode periodic
coordinates. Construct/cache geometry outside `threads`. Only gathered halo
entries are defined. Reverse exchange accumulates into physical fields.
Test: `thalo`.

HYP pullbacks project intermediates only within their forward construction
domains, respecting neighbor validity and cotangent reach. `thypsmearhalo`
compares stored floating-point components against independent shifts with
poisoned padding.

### Gauge actions and workspaces

| Family | Coefficients | Numerical entry points |
| --- | --- | --- |
| Fundamental | `plaq`, `rect`, `pgm` | `gaugeAction1/2/3`, `gaugeActionDeriv`, `gaugeDeriv2`, `gaugeForce*`, Hessians |
| Plaquette + adjoint | `plaq`, `adjplaq` | `actionA`, `gaugeADeriv`, `forceA` |

`action`/`force` select the adjoint family iff `adjplaq != 0`.
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

Call complete gauge operations outside `threads`. `newLoopWork(g[0])` owns work
for one layout/type: halos, improved-subset temporary gauge and full-Hessian
transporters. Calls refresh/release input references; a workspace serves one
evaluation at a time. Distinct workspaces own distinct storage. Communicator
sequencing and geometry initialization requirements still apply.

For a parity/direction mask $M_S$,

$$
g_S=M_S\nabla S,\qquad (Dg_S)^*b=\nabla^2S(M_Sb).
$$

The Hessian reaches every affected link. `Subset` overwrites; `SubsetAdd`
accumulates; `SubsetAddBase` adds on selected links and sets base plus Hessian
elsewhere; `SubsetSum` combines selected-parity seeds. `gaugeDeriv2SubsetWork`
preserves the complement for `clear=false` and zeros it for `clear=true`.
Rectangle/parallelogram subsets use the workspace temporary gauge.

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
order three. `LoopWork` reuses aligned per-thread site scratch, with no
intermediate lattice fields.

For plaquettes,

$$
P(U)=\sum_{x,\mu>\nu}\Re\operatorname{tr}
 [U_\mu(x)U_\nu(x+\hat\mu)U_\mu(x+\hat\nu)^\dagger U_\nu(x)^\dagger],
$$

$$
\operatorname{plaqSum}(U)=P(U),\qquad
\operatorname{stapleSum}(U,h)=D^{|h|}\nabla P(U)[h].
$$

`PlaqWork` fixes layout/order/action mode and rebinds inputs each call. It retains
those references until rebinding/release. Outputs are disjoint from seeds and,
below order three, the gauge. `newOneOf` allocates independent work.

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
`make tests` builds numerical suites; `make tests experimental` generates optional
suites. Generated scripts support the repository MPI runner. See the
[graph validation protocol](../src/experimental/graph/graph_validation.md).

`benchGaugeActions`, `benchPlaq` and `benchExpProject` use `tic`/`toc` and
`getElapsedTime`; parameters select lattice/repetitions/trials and, for site
benchmarks, operation/scale/precision/norm. `benchPlaq` compares production/fused
actions, forces and Hessians for Symanzik coefficients with `pgm=0` or `0.07`.
Warm workspaces first; validate results outside timing.

## Future work

- Measure general loop kernels on production hardware before further dispatch changes.
- Reuse generic halo message buffers.
- Add vector/rectangular LU right-hand sides when needed.
- Extend and validate adjoint-plaquette Hessian kernels.
