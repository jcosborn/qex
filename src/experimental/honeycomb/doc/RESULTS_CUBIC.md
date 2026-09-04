# Cubic reference pipeline — results

**Task C / PLAN.md task R step 5.** Everything here runs on the ordinary
hypercubic lattice with committed QEX code only; no honeycomb geometry is
involved. The purpose is (a) to build and validate the measurement/analysis
chain that the honeycomb tasks will reuse, and (b) to produce the *cubic*
reference curve of slide 10.

Files: `hcanalysis.nim`, `refCubicGen.nim`, `refCubicMeas.nim`,
`tests/tanalysis.nim`, `doc/plots/cubic/`.

---

## 0. Headline: the `topoQ` normalisation question is settled

> **QEX's `topoQ` (`src/gauge/gaugeUtils.nim:1274`, prefactor `-1/(4 π²`) is
> CORRECT. There is no factor of 2. Use `topoQ(f)` directly.**
>
> `doc/FORMULATION.md` §4.3 and `doc/PLAN.md` task W2 both say the prefactor is
> "suspect by a factor 2". **They are wrong** and should be corrected by their
> owners. The honeycomb clover charge must be normalised the same way as QEX's,
> i.e. `Q = -(1/(4π²)) Σ_x [Re tr(F₁₀F₃₂) − Re tr(F₂₀F₃₁) + Re tr(F₂₁F₃₀)]`
> with the traceless-**antihermitian** `F` (the one for which the plaquette is
> `P ≈ exp(a²F)`), equivalently `q(x) = −(1/32π²) ε_{μνρσ} tr(F̂_μν F̂_ρσ)`.
> Note the **minus sign** relative to FORMULATION §4.3 — that is just the
> hermitian/antihermitian convention; the magnitude `1/32π²` there is right.

### Where the "factor 2" suspicion came from

Combining the frequently-quoted
`Q = (1/32π²) ∫ ε_{μνρσ} F^a_{μν} F^a_{ρσ}` with the (correct) identity
`F^a_{μν}F^a_{ρσ} = −2 tr(F_{μν}F_{ρσ})` gives `Q = −(1/2π²)(a−b+c)`, twice
QEX's value. **The first of those two inputs is a mis-remembered formula.**
The textbook identity is

```
  Q = (1/32π²) ∫ F^a_{μν} F̃^a_{μν}          F̃_{μν} = ½ ε_{μνρσ} F_{ρσ}
    = (1/64π²) ∫ ε_{μνρσ} F^a_{μν} F^a_{ρσ}
```

— the `εFF` form carries `1/64π²`, not `1/32π²`. Check on BPST: a self-dual
instanton has `∫F^aF^a = 32π²` (with `S = ¼∫F^aF^a = 8π²/g²`), so
`(1/32π²)∫F^aF̃^a = 1` ✅, whereas `(1/32π²)∫εF^aF^a` would give 2.
With the correct `1/64π²` and `F^aF^a = −2tr(FF)` one lands exactly on
`Q = −(1/4π²)(a−b+c)`: QEX's expression.

### The numerical proof — exact, not statistical

`refCubicMeas -abeliantest` builds a **constant-field-strength Abelian SU(3)**
configuration in the Cartan direction `T = diag(1,−1,0)`:

```
  U₁(x) = exp(i φ₁ x₀ T)                 φ₁ = 2π n₁/(L₀L₁)
  U₀(x) = 1  except at x₀=L₀−1: exp(−i φ₁ L₀ x₁ T)
  U₃(x) = exp(i φ₂ x₂ T)                 φ₂ = 2π n₂/(L₂L₃)
  U₂(x) = 1  except at x₂=L₂−1: exp(−i φ₂ L₂ x₃ T)
```

Every `(0,1)` plaquette is then exactly `exp(iφ₁T)`, every `(2,3)` plaquette
exactly `exp(iφ₂T)`, and all other plaquettes are 1. The bundle is a direct sum
of three U(1) line bundles with charges `q = (1,−1,0)` and fluxes `q_i n_k`, so
by Atiyah–Singer the **exact** topological charge is

```
  Q_exact = Σ_i q_i² n₁ n₂ = 2 n₁ n₂ .
```

For this configuration the clover is exact: `fmunu` returns
`F₀₁ = i sin(φ₁) T`, `F₂₃ = i sin(φ₂) T`, so QEX must return
`topoQ = (V/2π²) sin φ₁ sin φ₂ → 2 n₁ n₂`. Measured
(`doc/plots/cubic/abeliantest.dat`):

| L | n₁ | n₂ | Q_exact | QEX `topoQ` | ratio |
|---|---|---|---|---|---|
| 8  | 1 | 1 |  2 |  1.993582728090 | 0.9967914 |
| 12 | 1 | 1 |  2 |  1.998731082901 | 0.9993655 |
| 16 | 1 | 1 |  2 |  1.999598437023 | 0.9997992 |
| 24 | 1 | 1 |  2 |  1.999920673805 | 0.9999603 |
| 32 | 1 | 1 |  2 |  1.999974900425 | 0.9999875 |
| 16 | 1 | 2 |  4 |  3.997992390718 | 0.9994981 |
| 16 | 2 | 3 | 12 | 11.984345658302 | 0.9986955 |
| 16 | −1| 2 | −4 | −3.997992390718 | 0.9994981 |
| 24 | 3 | 3 | 18 | 17.993575393716 | 0.9996431 |

The residual is the expected `O(a²)` clover artefact
`1 − (φ₁²+φ₂²)/6` to all printed digits (e.g. L=16, n=1: predicted deviation
`2×(2π/256)²/6 = 2.008e−4`, measured `2.008e−4`), and it scales as `1/L⁴`
exactly as it must. The closed-form lattice prediction
`(V/2π²) sin φ₁ sin φ₂` is reproduced to `1.6e−13`. `densityE` is likewise
verified: `E = 2sin²φ₁ + 2sin²φ₂` to `2e−16`.

**A factor 2 is excluded at the 10⁻³ level**: `2×topoQ` would be 3.999… where
the exact answer is 2, in every row. (`refCubicMeas -abeliantest` prints
`2×topoQ` for comparison.)

### Supporting evidence from real configurations

`doc/plots/cubic/qflow.dat` is the flow history of one β=5.8, 12⁴
configuration with the Symanzik-improved clover (`fmunuloop:5`). `Q` reaches a
clean plateau (figure below) at **3.9936** at `t = 5 ≈ 4 t₀` (and is 3.9906 already at
`t = 3`), i.e. the integer 4. Note this alone cannot distinguish `Q` from
`2Q` (4 and 8 are both integers) — which is exactly why the Abelian test above
was needed.

![qflow](plots/cubic/qflow.png)

**Warning for the honeycomb topology task:** with the *plain* clover
(`fmunuloop:1`) at `a²/t₀ ≳ 0.5` the charge is **not** close to an integer at
`t ≈ t₀`; deviations of 0.2–0.45 are normal, and at β=5.6 (`a²/t₀ ≈ 1.5`) the
charge does not even plateau under flow — instantons of size `~a` dissolve. Do
not use "Q is near an integer" as a normalisation test at coarse spacing.

---

## 1. `hcanalysis.nim` — the reusable analysis module

Pure Nim (only `std/math`, `std/strutils`); no QEX dependency, so it compiles
and runs in under 2 s and the honeycomb tasks can use it unchanged.

| symbol | what it does |
|---|---|
| `findT0(t, t2E, target=0.3, order=1)` | flow time where `t²⟨E⟩ = target`; `order=1` linear (the usual convention), `order=3` cubic-through-4-points solved by bisection |
| `findW0(t, tdt2E, target=0.3)` | **returns `w0²`** (the flow time at the crossing), the direct analogue of `findT0`; take `sqrt` for `w0` |
| `findCrossing(x, y, target, order)` | the generic first-upward-crossing finder behind both |
| `derivT2E(t, t2E)` | `W(t) = t d/dt[t²⟨E⟩]` by centred differences on a possibly non-uniform grid |
| `mean`, `variance`, `stddev`, `stderrMean` | elementary statistics |
| `autocorr`, `autocorrTime`, `autocorrTimeW` | `ρ(t)` and `τ_int` with the Madras–Sokal automatic window (`W ≥ 5 τ_int`) plus its error |
| `binned(x, b)` | block averages |
| `jackknife(x, f, bin=1)`, `jackknifeMean` | delete-`bin` jackknife of an arbitrary estimator |
| `fitPoly(x, y, dy, powers)` / `fitPolyCov` | weighted linear least squares in arbitrary integer powers, with covariance and `χ²/dof` |
| `evalPoly`, `readColumns` | helpers |
| `noCrossing = -1.0` | sentinel for "target never reached" — **not** `NaN`, because QEX builds with `-Ofast -ffast-math` under which `x != x` is optimised away |

### `tests/tanalysis.nim` — 23/23 PASS

```
export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
cd build_mac && make src/experimental/honeycomb/tests/tanalysis.nim && ./bin/tanalysis
```

* `findT0` exact on a linear curve; on `t²E = A t³` the linear interpolation
  error scales as `dt²` (9.8e−5 at dt=0.05) and the cubic one is `2e−16`;
* `findW0` + `derivT2E` recover the analytic crossing of `W(t)=0.3`;
* jackknife of the mean reproduces `sqrt(var/n)` to `1e−14` **exactly**;
  jackknife of the variance matches `sqrt(2/n)`; binning inflates the error of
  a block-correlated series by ≈`sqrt(bin)`;
* `autocorrTime` gives 0.51 for white noise (exact 0.5) and 4.45 for AR(1) with
  `a=0.8` (exact 4.5);
* `fitPoly` recovers exact coefficients for `[0,1]` and `[0,1,2]`, gives
  `⟨χ²/dof⟩ = 1.014` over 400 Gaussian-noise replicas with 66.5 % 1σ coverage
  of the intercept, and correctly ignores a point with a huge error bar.

---

## 2. `refCubicGen.nim` — quenched SU(3) HMC

Wilson plaquette action, `GaugeActionCoeffs(plaq: β)` + `gaugeForce`, momenta
from `randomTAH`, Omelyan integrators from `algorithms/integrator`
(`-gintalg:`, default `2MN`; **`4MN5FV` with `nsteps=4` is what was used** —
it gives >92 % acceptance at 20 force evaluations per unit of MD time, better
than `2MN` with `nsteps=12` at 24 force evaluations).

```
OMP_NUM_THREADS=6 ./bin/refCubicGen -lat:12,12,12,12 -beta:5.8 \
   -ntraj:4000 -nwarm:300 -tau:1.0 -nsteps:4 -gintalg:4MN5FV \
   -savefreq:20 -outdir:$TMPDIR/qexref/b58L12 -revCheckFreq:2000
```

### HMC correctness

| ensemble | trajs | acceptance | `⟨exp(−ΔH)⟩` | `⟨plaq⟩` | literature `⟨plaq⟩` |
|---|---|---|---|---|---|
| β=5.6, 8⁴  | 6000 | 98.7 % | 0.99997(48)  | 0.524163(42) | 0.5242 |
| β=5.7, 12⁴ | 4000 | 96.6 % | 1.00053(140) | 0.549085(21) | 0.5493 |
| β=5.8, 12⁴ | 4000 | 96.1 % | 0.99826(152) | 0.567602(17) | 0.5670 |
| β=5.9, 16⁴ | 1379 | 92.6 % | 1.00182(476) | 0.581855(15) | 0.5822 |
| β=6.0, 16⁴ | 1537 | 92.6 % | 1.00175(499) | 0.593703(14) | 0.5937 |

(errors are naive standard errors of the mean, i.e. they ignore
autocorrelation; `⟨exp(−ΔH)⟩` is over **all** proposals, accepted and
rejected, which is what the identity `⟨exp(−ΔH)⟩ = 1` refers to.)

`⟨exp(−ΔH)⟩ = 1` within errors everywhere, acceptance ≫ 70 %, and reversibility
`|ΔH_fwd + ΔH_bwd| ≤ 5×10⁻¹⁰` on every check. Plaquettes agree with the known
Wilson-action values.

---

## 3. `refCubicMeas.nim` — flow, `t₀`, `Q`

Uses `gauge/wflow.gaugeFlow` (Lüscher RK3) and `gaugeUtils.fmunu/densityE/topoQ`.
For each configuration it records `(t, t²⟨E⟩, W(t), Q(t))`, extracts
`t₀`, `w₀²` and `Q(t₀)`, and then jackknifes the ensemble.

Modes: default measurement, `-flowdump:` (full history of the first config),
`-abeliantest:1` (§0), `-chifit:<file>` (continuum extrapolation of a
three-column table, using `hcanalysis.fitPoly`).

### Practical lessons that the honeycomb flow task will hit too

1. **`t₀` does not exist if the box is too small.** On 8⁴ at β=6.0
   (`L/√t₀ ≈ 4.5`) `t²⟨E⟩` rises only to **0.127** at `t ≈ 2` and then *falls*
   — the flow radius `√(8t)` reaches `L/2` and the field is smoothed to a
   constant. `t²⟨E⟩ = 0.3` is never reached. A rule of thumb of `L/√t₀ ≳ 4` is far too
   loose; **`L/√t₀ ≳ 9–10` is needed** (`√(8t₀)/L ≲ 0.3`). All (β, L) used
   below have `L/√t₀ = 8.8 … 12.1`.
2. **QEX's SIMD layout rejects many lattice sizes.** With the standard
   `vlen = 4` build, `newLayout` splits two directions by 2 and then requires
   the resulting outer geometry to be even, so `L = 10` and `L = 14` fail with
   `can't lay out inner geom`. Use `L ∈ {8, 12, 16, 20, 24}`.
3. **Thread oversubscription is catastrophic, not gradual.** On this 16-logical-core
   machine, one `refCubicGen` at 8⁴: 1 thread 1.48 s, 4 threads 0.53 s,
   8 threads 0.54 s, **16 threads 17.5 s** (33× slower than 4). QEX's spin-wait
   barriers thrash. Use 4–6 threads per process, and do not run more concurrent
   jobs than there are physical cores.
4. RK3 flow stability requires `eps ≲ 0.15` (`eps·p̂²_max ≲ 2.5` with
   `p̂²_max = 16`); `eps ≈ t₀/25` was used, which is comfortably inside that and
   accurate.

---

## 4. Physics run — the cubic curve of slide 10

Five ensembles, `(β, L)` chosen so that `L/√t₀ ≈ 9.8–12.1` (≈ 1.6–2.0 fm with
`√(8t₀) = 0.47 fm`), i.e. roughly fixed physical volume, and `L ∈ {8,12,16}`
because of the SIMD-layout restriction. `t₀/a²` was **measured**, not assumed.

`Q ≡ Q(t₀)` from the plain clover, `χ = ⟨Q²⟩/V`, `V = L⁴a⁴`; `t₀` and `⟨Q²⟩`
jackknifed **together** over the same configurations (they are correlated).

| β | L | n_cfg | t₀/a² | a²/t₀ | L/√t₀ | ⟨Q⟩ | ⟨Q²⟩ | **10⁴ t₀²χ** | τ_int(Q) | n_eff |
|---|---|---|---|---|---|---|---|---|---|---|
| 5.6 | 8  | 197 | 0.66886(266) | 1.4951(59) | 9.78  | −0.013 | 1.741(187) | **1.902(201)** | 17 | 170 |
| 5.7 | 12 | 133 | 0.98604(399) | 1.0142(41) | 12.09 | −0.132 | 6.425(820) | **3.013(384)** | 22 | 60 |
| 5.8 | 12 | 159 | 1.51760(1032)| 0.6589(45) | 9.74  | +0.049 | 3.672(446) | **4.079(496)** | 22 | 74 |
| 5.9 | 16 | 55  | 2.30886(2266)| 0.4331(43) | 10.53 | −0.483 | 6.852(1319)| **5.574(1112)**| 24 | 14 |
| 6.0 | 16 | 59  | 3.31134(4853)| 0.3020(44) | 8.79  | +0.121 | 2.136(577) | **3.573(935)** | ≳31 | ≲11 |

`τ_int(Q)` is in **trajectories** (τ_int measured in units of the save interval
× the save interval, which is 30/20/20/12/12 trajectories of τ=1);
`n_eff = n_cfg/(2 τ_int)` in units of the save interval.
Data: `doc/plots/cubic/chitop.dat`.

`t₀/a²` agrees with the Necco–Sommer `r₀/a` parametrisation combined with
`t₀/r₀² = 0.1107` to 2–6 % (0.953 / 1.491 / 2.226 / 3.190 expected for
β = 5.7…6.0 vs 0.986 / 1.518 / 2.309 / 3.311 measured), which is an independent
check that the flow, `E` and the `t²⟨E⟩ = 0.3` convention are all right.

### Continuum extrapolation

![chitop](plots/cubic/chitop.png)

| fit | points | `10⁴t₀²χ|_{a=0}` | slope `c₁` | χ²/dof |
|---|---|---|---|---|
| **O(a²)**, β = 5.6…5.9 (**primary**) | 4 | **5.91 ± 0.68** | −2.70 ± 0.51 | 0.38 |
| O(a⁴), β = 5.6…5.9 | 4 | 7.65 ± 2.41 | −6.29 ± 4.78 (c₂ = 1.64 ± 2.17) | 0.19 |
| O(a²), all 5 points | 5 | 5.43 ± 0.59 | −2.36 ± 0.45 | 0.92 |
| O(a²), β = 5.6…5.8 only | 3 | 5.65 ± 0.75 | −2.51 ± 0.56 | 0.08 |
| reference (Cè et al., 1506.06052) | — | **6.67 ± 0.07** | −1.6 (slide 10) | — |

**How close it lands, and why.**

* The O(a²) extrapolation gives **5.91 ± 0.68**, i.e. **1.1 σ below** the
  literature value 6.67(7); the O(a⁴) extrapolation gives **7.65 ± 2.41**,
  1.4 σ high but with a 31 % error. The result is *consistent* with 6.7 but is
  not a precision determination — our error is 12 %, versus 1 % in Cè et al.
* The dominant limitation is the **lattice-spacing range**. Our points cover
  `a²/t₀ = 0.30 … 1.50`; slide 10's red points cover `0.15 … 0.35`. Over our
  range the cut-off effect is enormous — `10⁴t₀²χ` falls from 6.7 to 1.9, a
  factor 3.5 — so a *pure* `O(a²)` form is being asked to work far outside its
  domain, and the resulting intercept is biased low by exactly the amount seen
  (the O(a⁴) fit, which has the freedom to curve, recovers 7.65). Note that our
  fitted O(a²) slope −2.70(51) is considerably steeper than the −1.6 of slide
  10's red line, which is what one expects when the fit is contaminated by
  higher orders. **This is the same physics the paper is about**: the cubic
  lattice needs `a²/t₀ ≲ 0.35` for a controlled extrapolation, whereas the
  16-cell data on slide 10 are flat out to `a²/t₀ ≈ 1`.
* The second limitation is **topological freezing**. At β = 6.0
  (`a ≈ 0.09 fm`, 16⁴) the HMC charge performs a slow random walk with
  `|Q| ≲ 2` over the whole ~1500-trajectory stream: `⟨Q²⟩` was still creeping
  upwards (1.08 → 1.63 → 2.21 → 2.48 → 2.14 as configurations accumulated,
  against an expected ≈ 3.7) and the Madras–Sokal window for `Q²` reached a
  quarter of the sample. Its error is therefore **not trustworthy** and it is
  excluded from the primary fit (open grey symbol on the plot). β = 5.9 is
  marginal — about 14 independent samples, `⟨Q⟩ = −0.48` still 1.5 σ from 0.
  β ≤ 5.8 is fine (`⟨Q⟩` consistent with 0, frequent sign changes,
  `n_eff ≈ 60–170`).
* Definition systematic: measuring `Q` at a longer flow time raises `χ`.
  On 100 β=5.6 configurations, `10⁴t₀²χ` from `Q(t₀)` is **1.86(29)** and from
  `Q(4t₀)` is **2.61(39)** — a 40 % shift, i.e. comparable to the statistical
  error and a genuine O(a²) ambiguity of the definition at this very coarse
  spacing (`t₀ = 0.67a²`, so the flow radius at `t₀` is only 2.3 a). We quote
  `Q(t₀)`, matching the paper's `Q_flow = Q_flow(t=t₀)` (slide 17); a longer
  reference flow time would move all points up somewhat and the continuum
  value not at all. Data: `doc/plots/cubic/data/cfgmeas.b56L8.longflow.log`.

### Topological charge distribution (paper Fig. 1 analogue)

![qhist](plots/cubic/qhist.png)

197 configurations at β=5.6 on 8⁴: `⟨Q⟩ = −0.013`, `⟨Q²⟩ = 1.74`, symmetric and
close to Gaussian, no sign of a stuck sector.

### Cost actually spent

~70 minutes of wall clock on Monte Carlo (5 generation streams + up to 6
measurement jobs, 10–19 threads total on a 16-logical-core machine), which
overshoots the ~40 min guidance. Two reasons, both recorded here so the
honeycomb runs can avoid them: (i) the runs were sized from the single-job
6-thread throughput and then lost a factor ~5 to oversubscription once five
jobs ran at once; (ii) background jobs cannot be killed in this sandbox
(`kill`, `pkill` and `killall` are all denied), so an over-long run cannot be
cut short. **Size the job list for the number of physical cores first.**

Throughput (single job, 6 threads, 4MN5FV `nsteps=4`, τ=1):
0.03 s/traj at 8⁴, 0.15 s at 12⁴, 0.47 s at 16⁴; flow+`t₀`+`Q` measurement
1.2 s/config at 8⁴ (33 flow steps, 16 `fmunu`), and `fmunu` — not the flow —
dominates on small volumes.

---

## 5. Reproduction

```bash
export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
cd build_mac
make src/experimental/honeycomb/tests/tanalysis.nim && ./bin/tanalysis
make src/experimental/honeycomb/refCubicGen.nim
make src/experimental/honeycomb/refCubicMeas.nim

# the topoQ verdict (seconds)
./bin/refCubicMeas -abeliantest:1 -lat:16,16,16,16 -n1:1 -n2:1

# ensembles and measurement
P=../src/experimental/honeycomb/doc/plots/cubic
sh $P/run_gen.sh  $TMPDIR/qexref     # 5 HMC streams
sh $P/run_meas.sh $TMPDIR/qexref 6   # flow + t0 + Q, writes the ENSEMBLE lines

# chitop.dat + fits + PNG.  harvest.py re-derives the ENSEMBLE numbers from the
# per-configuration CFG lines (digit-for-digit identical to refCubicMeas's own
# output, verified on the 200-configuration beta=5.6 ensemble); it exists so
# that a measurement run that had to be truncated can still be analysed, and so
# that several partial runs over interleaved subsets of one stream can be
# merged in Monte Carlo time order.
python3 $P/harvest.py $TMPDIR/qexref --skip 3 --frozen 6.0 > $P/chitop.dat
sh $P/mkplot.sh
gnuplot $P/qflow.gp $P/qhist.gp
```
