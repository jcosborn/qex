# NNFT tests

Model, precision, checkpoint and driver contracts are in the
[application guide](../README.md).

Routine tests are native Nim units with synthetic inputs. Kernel checks use
analytic derivatives, explicit periodic-coordinate sums, adjoint identities and
finite differences. Graph checks cover storage, aliases, plans and cloning.
I/O tests create their own temporary files. None requires Python or JAX.

| Layer | Entry points |
| --- | --- |
| Numerical NN kernels | `tests/nn/tnn.nim` |
| NN graph operations and derivatives | `graph/tests/tgnn.nim` |
| Real/matrix field bridges | `graph/tests/trfield.nim`, `trfieldu1.nim` |
| Integrator | `graph/tests/tintegratoru1.nim` |
| Learned stage and composed graph | `nnft/tests/tnnftunit.nim`, `tnnftexpr.nim` |
| Generic and application I/O | `tests/io/tarrays.nim`, `tarrayfields.nim`, `nnft/tests/tnnftio.nim` |

For focused runs from the repository root with a configured `build_mac`:

```sh
OMP_NUM_THREADS=1 make -C build_mac ARGS=--assertions:on run nn/tnn io/tarrays io/tarrayfields
OMP_NUM_THREADS=1 make -C build_mac ARGS=--assertions:on run \
  graph/tests/tgnn graph/tests/trfield graph/tests/trfieldu1 \
  graph/tests/tintegratoru1 nnft/tests/tnnftunit nnft/tests/tnnftexpr nnft/tests/tnnftio
```

`make tests` discovers the generic tests. Graph units use the optional
`experimental/graph` group; application units use `experimental/nnft`. For
example, `make run tests experimental/nnft` builds and runs the application units.
NN and learned-flow units exercise both neural precisions.

## Optional JAX comparison

`compare.py` and `compare.nim` form one application comparison outside unit-test
discovery. The Python reference contains the flow and trajectory equations and
creates its own model and inputs. It checks the transformed field, log determinant,
action, force and one short trajectory in both neural precisions on a synthetic
$8\times12$ case with two steps of `dt=0.08`; these are comparison settings, not
driver defaults. Gauge, momentum and scalar storage are double precision in both
implementations.

```sh
make -C build_mac ARGS=--assertions:on nnft/tests/compare
python3 src/experimental/nnft/tests/compare.py build_mac/bin/compare
```

This command requires JAX and NumPy. Temporary comparison files are removed by
default; `--output DIRECTORY` keeps inputs, outputs and the report in a new
directory. `--checkpoint MODEL.npz` optionally uses the application's forty
float32 parameter arrays. Retained `DIRECTORY/input` is also usable with the
application's `-checkpoint` option with `-lat:8,12 -latent:theta`.
