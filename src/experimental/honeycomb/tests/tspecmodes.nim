## Candidate processing and report fields without an eigensolve.

import std/[complex, strutils, unittest]
import qex except epsilon
import physics/qcdTypes
import ../hcSpectrum

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))
qexInit()

proc checkModes[F](vecs: seq[F]) =
  # D = 2 I. The last candidate has lambda = 5 + 4i, hence residual 5.
  let sig = -0.25
  let mu = complex64(1.0/(2.0 - sig), 0.0)
  let bad = complex64(5.0, 4.0)
  let mus = @[mu, mu, mu, complex64(1.0, 0.0)/(bad - complex64(sig, 0.0))]
  let vs = @[vecs[0], vecs[1], vecs[2], vecs[2]]
  let n2 = vnorm2(vecs[2])
  var calls = 0
  let applyD = proc(r: var F; x: F) =
    inc calls
    copyP(r, x)
    scaleP(r, 2.0)
  let modes = measureModes(mus, vs, sig, applyD)
  check modes.len == 4
  for i in 0..2:
    check abs(modes[i].lam - complex64(2.0, 0.0)) < 1e-13
    check modes[i].resid < 1e-13
  check abs(modes[0].chi - 1.0) < 1e-13
  check abs(modes[1].chi + 1.0) < 1e-13
  check abs(modes[2].chi - 0.6) < 1e-13
  check abs(modes[3].lam - bad) < 1e-13
  check abs(modes[3].resid - 5.0) < 1e-13
  check vnorm2(vecs[2]) == n2
  let kept = measureModes(mus, vs, sig, applyD, 0.5)
  check kept.len == 3
  let s = summarize(kept, 1e-6)
  check s.nconv == 3
  check s.nreal == 3
  check s.nplus == 2
  check s.nminus == 1
  check s.qdirac == -1.0
  check abs(s.sumChiReal - 0.6) < 1e-13
  # lambda = 1 gives residual exactly 1, including the cutoff boundary.
  let edge = measureModes(@[complex64(1.0, 0.0)], @[vecs[0]], 0.0, applyD, 1.0)
  check edge.len == 1
  check edge[0].resid == 1.0
  let n = calls
  check measureModes(newSeq[Complex64](), newSeq[F](), sig, applyD).len == 0
  expect ValueError:
    discard measureModes(newSeq[Complex64](), @[vecs[0]], sig, applyD)
  expect ValueError:
    discard measureModes(@[mu], newSeq[F](), sig, applyD)
  check calls == n

suite "spectrum mode processing":
  test "honeycomb candidate modes":
    let hl = newHcLayout([4, 4, 4, 4])
    var vs = @[newHcFermion(hl), newHcFermion(hl), newHcFermion(hl)]
    # gamma5 = diag(1, 1, -1, -1); mixed amplitudes 2 and 1 give chi = 3/5.
    for i in hl.lo.sites:
      setC(vs[0].a{i}[0][0], 3.0, 0.0)
      setC(vs[0].b{i}[0][0], 3.0, 0.0)
      setC(vs[1].a{i}[2][0], 2.0, 0.0)
      setC(vs[1].b{i}[2][0], 2.0, 0.0)
      setC(vs[2].a{i}[0][0], 2.0, 0.0)
      setC(vs[2].b{i}[0][0], 2.0, 0.0)
      setC(vs[2].a{i}[2][0], 1.0, 0.0)
      setC(vs[2].b{i}[2][0], 1.0, 0.0)
    checkModes(vs)

  test "cubic candidate modes":
    let lo = newLayout(@[4, 4, 4, 4])
    var vs = @[lo.DiracFermion(), lo.DiracFermion(), lo.DiracFermion()]
    for v in vs: v := 0
    for i in lo.sites:
      setC(vs[0]{i}[0][0], 3.0, 0.0)
      setC(vs[1]{i}[2][0], 2.0, 0.0)
      setC(vs[2]{i}[0][0], 2.0, 0.0)
      setC(vs[2]{i}[2][0], 1.0, 0.0)
    checkModes(vs)

  test "spectrum report fields":
    let modes = @[
      SpecMode(lam: complex64(0.5, 0.0), chi: -0.9, resid: 2e-5),
      SpecMode(lam: complex64(1.0, 0.25), chi: 0.0, resid: 3e-6),
      SpecMode(lam: complex64(1.0, -0.25), chi: 0.0, resid: 4e-6)]
    let st = SiStats(totIts: 96, maxUsedIts: 20, nHitMax: 1, worstDirect: 3e-5)
    let rows = specLines(modes, 1e-6, 7, 0.05, 0.75, 2, 12, st, 1.25)
    check rows.len == 5
    for i in 0..2:
      let cols = rows[i].splitWhitespace
      check cols.len == 8
      check cols[0..3] == @["EIG", "7", "0.050", $i]
      check cols[4].parseFloat == modes[i].lam.re
      check cols[5].parseFloat == modes[i].lam.im
      check cols[6].parseFloat == modes[i].chi
      check cols[7].parseFloat == modes[i].resid
    let cfg = rows[3].splitWhitespace
    check cfg.len == 14
    check cfg[0..2] == @["CFG", "7", "0.050"]
    check cfg[3].parseFloat == 0.5
    check cfg[4].parseFloat == 0.0
    check cfg[5..^1] == @["3", "1", "0", "1", "1.0", "0.750000", "12", "96", "1.25"]
    let ext = rows[4].splitWhitespace
    check ext.len == 21
    check ext[0..2] == @["CFGX", "7", "0.050"]
    let keys = ["conj", "worstresid", "imgapC", "imgapR", "sumchi", "nbad", "cgmax", "cghitmax", "sidirect"]
    let vals = [0.0, 2e-5, 0.25, 0.0, -0.9, 2.0, 20.0, 1.0, 3e-5]
    for i, key in keys:
      check ext[3 + 2*i] == key
      check ext[4 + 2*i].parseFloat == vals[i]

qexFinalize()
