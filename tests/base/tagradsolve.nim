import testutils
import qex, physics/stagSolve, hmc/agradOps

qexInit()

proc run(variable: static bool) =
  let lo = latticeFromLocalLattice([4,4,4,4], nRanks).newLayout
  var rng = lo.newRngField(MRG32k3a, 34129'u)
  let g = lo.newGauge
  let dg = lo.newGauge
  let b, db, res = lo.ColorVector()
  threads:
    g.random rng
    dg.random rng
    b.gaussian rng
    db.gaussian rng
  let gs = g.newGauge
  let bs = b.newOneOf
  threads: bs := b

  const mass = 0.4
  const eps = 1e-4
  var t = newAgTape()
  let gv = t.newAgVar(g)
  let bv = t.newAgVar(b)
  let rv = t.newAgVar(lo.ColorVector())
  let loss = t.newFloatV()
  when variable:
    let m = t.newFloatV(mass)
  else:
    let m = mass
  var s = newStag(g)
  var pf = initSolverParams()
  pf.backend = sbQex
  pf.sloppySolve = SloppyNone
  pf.r2req = 1e-22
  pf.maxits = 10000
  pf.verbosity = 0
  when variable:
    var pb = pf
    s.agradSolve(gv, rv, bv, m, pf, pb)
  else:
    s.agradSolve(gv, rv, bv, m, pf)
  loss.norm2subtract(rv, 0.0)

  suite "Staggered autograd solve, variable mass " & $variable:
    test "forward solves the full equation and restores the links":
      t.run
      s.rephase
      threads: s.D(res, rv.obj, mass)
      s.rephase
      check norm2(res-b) <= 1.01*pf.r2req*norm2(b)
      for mu in 0..<g.len: check norm2(g[mu]-gs[mu]) == 0

    test "pullbacks match finite differences":
      pf.resetStats()
      when variable: pb.resetStats()
      t.run
      loss.grad = 1.0
      t.grad
      let bg = redot(db, bv.grad)
      var gg = 0.0
      for mu in 0..<g.len: gg += redot(dg[mu], gv.grad[mu])
      when variable:
        let mg = m.grad
        check pf.calls == 1
        check pb.calls == 1
      else:
        check pf.calls == 2
      for mu in 0..<g.len: check norm2(g[mu]-gs[mu]) == 0

      var values: array[2,float]
      for i, step in [-eps, eps]:
        threads: b := bs + step*db
        t.run
        values[i] = loss.obj
      threads: b := bs
      var fd = (values[1]-values[0])/(2*eps)
      check abs(bg-fd) <= 2e-6*max(1.0, abs(fd))

      for i, step in [-eps, eps]:
        threads:
          for mu in 0..<g.len: g[mu] := gs[mu] + step*dg[mu]
        t.run
        values[i] = loss.obj
      threads:
        for mu in 0..<g.len: g[mu] := gs[mu]
      fd = (values[1]-values[0])/(2*eps)
      check abs(gg-fd) <= 2e-6*max(1.0, abs(fd))

      when variable:
        for i, step in [-eps, eps]:
          m.obj = mass + step
          t.run
          values[i] = loss.obj
        m.obj = mass
        fd = (values[1]-values[0])/(2*eps)
        check abs(mg-fd) <= 2e-6*max(1.0, abs(fd))

run(false)
run(true)
qexFinalize()
