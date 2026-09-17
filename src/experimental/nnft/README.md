# Learned U(1) field transformations

The application evolves a latent gauge field $V$ and measures the physical field
$U=f(V)$. Its effective action is

$$
S_{\rm eff}(V)=-\beta\sum_x\cos P_x(f(V)) - \ell(V),
$$

where $P_x$ is the plaquette angle and $\ell$ is the clipped angular log
determinant defined below. The learned flow requires a two-dimensional lattice
with even extents of at least four.

## Modules

| Module | Responsibility |
| --- | --- |
| [Numerical NN](../../nn.nim) | Real field channels, convolution, pointwise operations, masks and pullbacks |
| [NN graph](../graph/nn.nim) | Separate graph functions over shared field values and parameter arrays |
| [Array I/O](../../io/arrays.nim), [field I/O](../../io/arrayfields.nim) | NumPy/JAX array interchange and mapping through global lattice coordinates |
| [flow](flow.nim) | Numerical stage tapes and their field pullbacks |
| [model](model.nim) | Parameter leaves per stage and the network as graph functions |
| [expr](expr.nim), [graph](graph.nim) | Composed and fused graph implementations |
| [io](io.nim) | Checkpoint and latent-angle adapters |
| [Driver](nnfthmc.nim) | Model initialization, HMC and physical-field measurements |
| [Shared update](../../gauge/stoutsmear.nim) | Numerical exponential update and joint pullback |
| [Shared flow action](../graph/hmcgauge/flow.nim) | Cached transformation and effective action |

Numerical and graph NN operations consume in-memory values. They have no
checkpoint or application I/O dependency. Channels use existing QEX fields;
convolutions use halo exchange. Graph storage and derivative contracts are in
the [graph design](../graph/DESIGN.md#121-neural-fields).

The [subset-stout driver](../graph/pg2du1ftstouthmc.nim) constructs its own
smearing flow and inverse initialization. Both drivers use the shared numerical
update, transformed-action constructor and
[2D U(1) measurement loop](../graph/hmcgauge/measure2du1.nim).

## Network and flow

Each stage has an independent network of circular, stride-one convolutions
with odd kernels, GELU between them, from the six features to the twelve
coefficients (the reference checkpoint: two $3\times3$ layers, $6\to12\to12$):

```nim
var h = x
for l in 0..<p.weights.len:
  h = bias(conv(h, p.weights[l]), p.biases[l])
  if l < p.weights.high: h = gelu(h)
let z = scale(h, p.scale)
divide(divide(arctan(z), T(PI)), T(3))
```

GELU uses the exact erfc expression. The two divisions execute separately in
the selected neural precision. Features are ordered as
`[sin(P), cos(P), sin(R0), sin(R1), cos(R0), cos(R1)]`, where $P$ is the
plaquette angle and $R_\mu$ is the rectangle extended in direction $\mu$.
Masked sine/cosine features take values zero/one.

For stage $s$, with the eight classes repeating for $s\ge8$, the direction and
active coordinate parities are

$$
d=\lfloor (s\bmod8)/4\rfloor,\qquad
r=\lfloor(s\bmod4)/2\rfloor,\qquad c=s\bmod2.
$$

Only direction $d$ at sites with those row/column parities changes. Plaquette
features exclude the active rows for $d=0$, or columns for $d=1$. The transverse
rectangle feature excludes active sites; the other rectangle feature is masked
out. These masks make the network coefficients independent of active links.

Six network outputs weight open plaquette and rectangle staples $D_k$. On an
active link, the update and stage log determinant are

$$
d_s=-\sum_{k=0}^{5}c_kD_k,\qquad M=V_d d_s^\dagger,\qquad
V'_d=\exp(\operatorname{projTAH}(M))V_d,
$$
$$
\ell_s=\sum_{x\in A_s}\log\max(1+\Re M_x,10^{-8}).
$$

Inactive links are copied. The full flow applies all eight stages and sums
their log determinants; the clip follows the
[`clipMin` convention](../graph/DESIGN.md#121-neural-fields). Graph
log-determinant factorization requires model parameters independent of the
flow input.

`newNnftStage` creates one stage's numerical tape; `evalStage` fills its output
and returns the stage log determinant, and `stageVjp` pulls back through the
tape evaluated at the matching input and parameters. The fused graph stage owns
these calls. `toNnftModel` copies parameters into graph leaves; `model.update`
updates those leaves. `learnedAction` supplies the `flow` and `action` callables
used by the HMC driver.

Fused graph stages allocate their private tapes on evaluation. Field pullbacks
use the tape; parameter and higher input pullbacks use the composed stage with
distinct live slots for primal arguments and cotangents. Plans reuse surrounding
field buffers; the retained stage tape does not use gauge-only buffer pooling.

## Precision and input format

| Data | Precision |
| --- | --- |
| Neural parameters, features, activations and pullbacks | `float32` or `float64` |
| Gauge, momentum, graph scalars and log determinants | `float64` |
| Site masks | `float32`, zero/nonzero selection |
| Checkpoint files | `float32`, promoted when loading a double model |

The checkpoint directory contains `manifest.json` with `version: 1` and an
`arrays` object naming `param_0`, `param_1`, ... in stage order. Each entry is
an [array descriptor](../../io/arrays.nim) with `dtype: "float32"`. Within a
stage, every convolution contributes its bias `[out]` followed by its weights
`[out,in,k0,k1]`, and a scale `[12,1,1]` closes the stage; the shapes define
the channel widths, kernels and stage count. The reference checkpoint has eight
stages of five leaves:

| Array | Value | Shape |
| --- | --- | --- |
| `param_(5*s)` | First bias | `[12]` |
| `param_(5*s+1)` | First weights | `[12,6,3,3]` |
| `param_(5*s+2)` | Second bias | `[12]` |
| `param_(5*s+3)` | Second weights | `[12,12,3,3]` |
| `param_(5*s+4)` | Output scale | `[12,1,1]` |

Weights are stored output-channel first, input-channel second, with the last
spatial dimension fastest. Graph scale leaves have shape `[12]`.
Latent angles are an optional floating array of shape `[2,L0,L1]` in the same
manifest; `loadLatent` maps them to $V_d=\exp(i\theta_d)$ by global coordinates.

## Running

From the repository root with a configured `build_mac`:

```sh
make -C build_mac nnft/nnfthmc
OMP_NUM_THREADS=1 build_mac/bin/nnfthmc \
  -checkpoint:MODEL_DIR -precision:double
```

`MODEL_DIR` is the manifest directory described above. The driver has these defaults:

| Option | Learned default |
| --- | --- |
| `precision` | `single` |
| `lat`, `beta` | `8,8`, `3.0` |
| `seed`, `rng` | `1029`, Philox4x64 |
| `gintalg`, `gintcoeffs` | `2MNp`, the default lambda |
| `dt`, `gsteps` | `0.35`, `10` |
| `trajsThermo`, `trajs` | `0`, `1` |

The driver starts from the cold latent field. `-latent:NAME` selects an angle
array in the checkpoint manifest; `-lat` must match its spatial dimensions.
Physical gauge-file restarts are rejected because there is no learned inverse.
Measurements and optional saves (`-savefreq:N`, default `0`) use $U=f(V)$;
the driver does not save the evolving latent field, so these files cannot
resume a run. The driver samples with fixed parameters and does not construct
a parameter-training objective.

## Tests and comparison

The [test guide](tests/README.md) lists the native kernel, graph, I/O and
application units with their commands, and the optional JAX comparison of the
application in both neural precisions.
