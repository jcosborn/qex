#!/usr/bin/env python3
"""Optional NNFT/JAX comparison; takes a configured nnft/tests/compare binary.

The default model and 8x12 input are synthetic. Both references use double angles and
momenta, with float32 or float64 neural operations, as in the QEX application.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile

os.environ.setdefault("JAX_PLATFORMS", "cpu")
import jax
import jax.numpy as jnp
import numpy as np

jax.config.update("jax_enable_x64", True)


def loops(theta):
    u, v = theta
    p = u - v - jnp.roll(u, -1, 1) + jnp.roll(v, -1, 0)
    # Rectangles are the sum of two adjacent oriented plaquettes.
    return p, p + jnp.roll(p, -1, 0), p + jnp.roll(p, -1, 1)


def conv(x, w, b):
    x = jnp.pad(x[None], ((0, 0), (0, 0), (1, 1), (1, 1)), mode="wrap")
    y = jax.lax.conv_general_dilated(
        x, w, (1, 1), "VALID", dimension_numbers=("NCHW", "OIHW", "NCHW"))
    return y[0] + b[:, None, None]


def flow(theta, params):
    row, col = jnp.indices(theta.shape[1:])
    logdet = jnp.float64(0)
    for s, (b1, w1, b2, w2, scale) in enumerate(params):
        d, r, c = s // 4, s % 4 // 2, s % 2
        active = (row % 2 == r) & (col % 2 == c)
        pmask = row % 2 != r if d == 0 else col % 2 != c
        p, r0, r1 = loops(theta)
        fm = (p * pmask, r0 * (~active if d == 1 else False),
              r1 * (~active if d == 0 else False))
        x = jnp.stack([jnp.sin(fm[0]), jnp.cos(fm[0]),
                       jnp.sin(fm[1]), jnp.sin(fm[2]),
                       jnp.cos(fm[1]), jnp.cos(fm[2])]).astype(w1.dtype)
        x = jax.nn.gelu(conv(x, w1, b1), approximate=False)
        k = (jnp.arctan(conv(x, w2, b2) * scale) / jnp.asarray(np.pi, w1.dtype)
             / jnp.asarray(3, w1.dtype)).astype(jnp.float64)
        if d == 0:
            angles = jnp.stack([p, jnp.roll(p, 1, 1), jnp.roll(r0, 1, 0),
                               jnp.roll(r0, (1, 1), (0, 1)), r0, jnp.roll(r0, 1, 1)])
            k = k[jnp.array([0, 1, 4, 5, 6, 7])]
            signs = jnp.array([-1, 1, -1, 1, -1, 1])
        else:
            angles = jnp.stack([p, jnp.roll(p, 1, 0), jnp.roll(r1, 1, 1),
                               jnp.roll(r1, (1, 1), (0, 1)), r1, jnp.roll(r1, 1, 0)])
            k = k[jnp.array([2, 3, 8, 9, 10, 11])]
            signs = jnp.array([1, -1, 1, -1, 1, -1])
        shift = jnp.sum(k * jnp.sin(angles) * signs[:, None, None], axis=0)
        diagonal = 1 - jnp.sum(k * jnp.cos(angles), axis=0)
        logdet += jnp.sum(jnp.where(active, jnp.log(jnp.maximum(diagonal, 1e-8)), 0))
        theta = theta.at[d].add(jnp.where(active, shift, 0))
    return theta, logdet


def reference(theta, momentum, params, cfg):
    def action(x):
        y, ld = flow(x, params)
        return -cfg["beta"] * jnp.cos(loops(y)[0]).sum() - ld

    original = theta
    value_grad = jax.jit(jax.value_and_grad(action))
    y, ld = jax.jit(lambda x: flow(x, params))(theta)
    initial, force = value_grad(theta)
    h0 = initial + (momentum * momentum).sum() / 2
    dt, lam = cfg["dt"], cfg["lambda"]
    p = momentum - lam * dt * force
    for i in range(cfg["steps"]):
        theta = theta + dt / 2 * p
        p = p - (1 - 2 * lam) * dt * value_grad(theta)[1]
        theta = theta + dt / 2 * p
        p = p - (lam if i + 1 == cfg["steps"] else 2 * lam) * dt * value_grad(theta)[1]
    final = value_grad(theta)[0]
    h1 = final + (p * p).sum() / 2
    dh = h1 - h0
    probability = jnp.minimum(1, jnp.exp(-dh))
    accepted = cfg["uniform"] < probability
    committed = jnp.where(accepted, theta, original)
    return {name: np.asarray(value, dtype=np.float64) for name, value in {
        "flow_re": jnp.cos(y), "flow_im": jnp.sin(y), "logdet": ld,
        "action": initial, "force": force, "proposal_re": jnp.cos(theta),
        "proposal_im": jnp.sin(theta), "momentum": p, "H0": h0, "H1": h1,
        "deltaH": dh, "probability": probability,
        "accepted": accepted, "committed_re": jnp.cos(committed), "committed_im": jnp.sin(committed),
    }.items()}


def spec(name, value):
    return {"file": name + ".bin", "shape": list(value.shape),
            "dtype": str(value.dtype), "byte_order": "little", "order": "C"}


def write(path, value):
    np.asarray(value, dtype=value.dtype.newbyteorder("<")).tofile(path)


def check(actual, expected, precision):
    atol, rtol = (2e-6, 2e-6) if precision == "single" else (5e-11, 5e-12)
    errors = {}
    for name, value in expected.items():
        got = np.fromfile(actual / (name + ".bin"), dtype="<f8").reshape(value.shape)
        np.testing.assert_allclose(got, value, atol=atol, rtol=rtol, equal_nan=False,
                                   err_msg=precision + " " + name)
        errors[name] = float(np.max(np.abs(got - value)))
    return {"arrays": len(errors), "nn_dtype": "float32" if precision == "single" else "float64",
            "atol": atol, "rtol": rtol, "max_errors": errors}


def run(binary, directory, checkpoint=None):
    inputs = directory / "input"
    inputs.mkdir()
    cfg = {"beta": 3.0, "dt": 0.08, "lambda": 0.1931833275037836, "steps": 2, "uniform": 0.73}
    manifest = {"version": 1, "config": cfg, "arrays": {}, "results": {}}
    stored = None
    if checkpoint:
        with np.load(checkpoint, allow_pickle=False) as archive:
            stored = {"param_" + str(i): archive["param_" + str(i)] for i in range(40)}
    params = []
    for s in range(8):
        layer = []
        for j, (shape, amp) in enumerate([((12,), .12), ((12, 6, 3, 3), .05),
                                         ((12,), .17), ((12, 12, 3, 3), .04), ((12, 1, 1), .65)]):
            name = "param_" + str(5 * s + j)
            if stored is None:
                a = np.arange(np.prod(shape), dtype=np.float64).reshape(shape)
                value = (amp * np.sin(.37 * a + .61 * s + .23 * j)).astype(np.float32)
                if j == 4:
                    value += np.float32(.8)
            else:
                value = stored[name]
                if value.shape != shape or value.dtype != np.dtype("float32"):
                    raise ValueError(name + " must be a float32 array of shape " + str(shape))
            manifest["arrays"][name] = spec(name, value)
            write(inputs / (name + ".bin"), value)
            layer.append(value)
        params.append(layer)
    d, row, col = np.indices((2, 8, 12))
    theta = .7 * np.sin(.43 * row + .29 * col + .83 * d) + .05 * (row - col)
    momentum = .6 * np.cos(.31 * row - .37 * col + .71 * d)
    for name, value in [("theta", theta), ("initial_momentum", momentum)]:
        manifest["arrays"][name] = spec(name, value)
        write(inputs / (name + ".bin"), value)
    report = {"lattice": [8, 12], "gauge_momentum_scalars": "float64",
              "parameters": "float32 storage, promoted for double NN evaluation"}
    for precision, dtype in [("single", jnp.float32), ("double", jnp.float64)]:
        print("Comparing " + precision, flush=True)
        p = [[jnp.asarray(x, dtype) for x in layer] for layer in params]
        expected = reference(jnp.asarray(theta), jnp.asarray(momentum), p, cfg)
        ref = directory / (precision + "-reference")
        actual = directory / precision
        ref.mkdir()
        actual.mkdir()
        for name, value in expected.items():
            manifest["results"][name] = spec(name, value)
            write(ref / (name + ".bin"), value)
        (inputs / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        command = [str(binary), "-input:" + str(inputs), "-output:" + str(actual), "-precision:" + precision]
        env = dict(os.environ)
        env.setdefault("OMP_NUM_THREADS", "1")
        logpath = directory / (precision + ".log")
        with logpath.open("w") as log:
            result = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError("QEX comparison failed:\n" + logpath.read_text())
        report[precision] = check(actual, expected, precision)
        print(precision + ": " + str(len(expected)) + " arrays passed", flush=True)
    report["passed"] = True
    (directory / "summary.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path, help="compiled nnft/tests/compare")
    parser.add_argument("--output", type=Path, help="new directory in which to retain inputs, outputs and logs")
    parser.add_argument("--checkpoint", type=Path, help="optional NPZ containing float32 param_0 through param_39")
    args = parser.parse_args()
    binary = args.binary.resolve(strict=True)
    if args.output:
        directory = args.output.resolve()
        directory.mkdir(parents=True)
        run(binary, directory, args.checkpoint)
    else:
        with tempfile.TemporaryDirectory(prefix="qex-nnft-") as tmp:
            run(binary, Path(tmp), args.checkpoint)
