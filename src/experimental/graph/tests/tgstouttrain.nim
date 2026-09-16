#RUNCMD env OMP_NUM_THREADS=1 $RUN1
import base/globals
setVLENmax(4)

import math, unittest
import helpers
when defined(graphPlanMemory):
  import std/tables
import qex except epsilon
import base/alignedMem
import algorithms/numdiff
import ../[core, scalar, gauge, plan, functional]
import ../hmcgauge/[ftstout, trajectory, integrator, config, training]

proc runStoutTrainingTests*(localLat: seq[int], beta: float) =
  qexInit()
  defer: qexFinalize()
  letParam:
    expectRanks = nRanks
  check nRanks == expectRanks
  letParam:
    mode = "planned"
  let planned = case mode
    of "direct": false
    of "planned": true
    else: raiseValueError("mode must be direct or planned")
  let rawGcThreshold = getRawMemGcThreshold()
  defer: setRawMemGcThreshold(rawGcThreshold)
  # Retained outputs and plan arenas stay live throughout the checks.
  # Forced collection cannot reclaim them while the fixture runs.
  setRawMemGcThreshold(int.high)
  when defined(graphPlanMemory):
    letParam:
      memoryOnly = false
    if memoryOnly and not planned:
      raiseValueError("memoryOnly requires -mode:planned")
  echo "Training execution: ", mode
  letParam:
    lat = latticeFromLocalLattice(localLat,nRanks)
  let lo = lat.newLayout
  let g = lo.newGauge
  let p = lo.newGauge
  const nc = g[0][0].nrows
  threads:
    for mu in 0..<g.len:
      g[mu] := 1.0
      p[mu] := 0.0
      for e in p[mu]:
        let v = 0.35*math.sin(float(3*e+2*mu+1))
        p[mu][e][0,0].im := v
        when nc > 1:
          p[mu][e][1,1].im := -v
  let
    rt = initGraphRuntime()
    rho = scalar.toGvalue(rt,0.025)
    coeff = actWilson(scalar.toGvalue(rt,beta))
    sa = stoutAction(coeff,rho,1)
    cfg = RunConfig(dt:0.1,gsteps:1,trajs:1,trajsTrain:1,lrmax:0.0001,lrmin:0.0001,
      integratorCoeffs:IntegratorCoeffs(kind:ik2MN,lambda:0.0))
  cfg.validateRunConfig
  proc memory(label: string) =
    if planned:
      when defined(graphPlanMemory):
        let
          occupied = getOccupiedMem()
          total = getTotalMem()
          raw = getRawMemUsed()
        GC_fullCollect()
        let
          collected = getOccupiedMem()
          collectedTotal = getTotalMem()
          collectedRaw = getRawMemUsed()
        echo "graph-plan-memory stage=", label,
          " assignedIds=", rt.nextStableNodeId,
          " runtimeApplyEntries=", rt.functional.applyCacheByNode.len,
          " runtimeRunStats=", rt.runStatsByNode.len,
          " occupiedBeforeGc=", occupied, " occupiedAfterGc=", collected,
          " heapBeforeGc=", total, " heapAfterGc=", collectedTotal,
          " rawUsedBeforeGc=", raw, " rawUsedAfterGc=", collectedRaw,
          " rawAllocated=", getRawMemAllocated()
      else:
        echo "  training memory ", label, ": assigned IDs ", rt.nextStableNodeId,
          ", occupied ", getOccupiedMem(), ", total ", getTotalMem(),
          ", rawUsed ", getRawMemUsed(), ", rawAllocated ", getRawMemAllocated()
  memory("before trajectory graph")
  let graph = buildTrajectoryGraph(rt,g,p,sa.action,cfg,buildTraining=true,
    parameters = [LearnedParameter(name:"flowRho",node:rho)])
  memory("after trajectory graph")
  defer:
    rt.resetGradCache(false)
    rt.resetApplyCache
    rt.resetLdjCache
  memory("before mixed derivatives")
  let
    dt = graph.learnedParameters[0].node
    dr = graph.learnedParameters[2].gradientExpr
    mixed1Expr = Gscalar(grad(dr,dt))
    mixed2Expr = Gscalar(grad(graph.learnedParameters[0].gradientExpr,rho))
  memory("after mixed derivatives")
  var roots = @[Gvalue(graph.deltaHamiltonian),Gvalue(graph.lossExpr)]
  for lp in graph.learnedParameters:
    roots.add Gvalue(lp.gradientExpr)
  roots.add Gvalue(mixed1Expr)
  roots.add Gvalue(mixed2Expr)
  memory("source-complete")
  let joint = if planned: plan(roots) else: nil
  when not defined(graphPlanMemory):
    memory("after joint plan")
  var firstEval = true

  proc value(x: Gscalar): float =
    if joint == nil:
      return x.eval.sval
    for i, root in roots:
      if root.nodeKey == x.nodeKey:
        if firstEval:
          when defined(graphPlanMemory):
            joint.reportMemory("before-first-evaluation")
          else:
            memory("before first joint evaluation")
        discard joint.eval
        if firstEval:
          when defined(graphPlanMemory):
            joint.reportMemory("first-evaluation-complete")
          else:
            memory("after first joint evaluation")
          echo "  training first plan arenaBytes ", joint.stats.arenaBytes,
            ", peakLiveBytes ", joint.stats.peakLiveBytes,
            ", workspaceBytes ", joint.stats.workspaceBytes,
            ", workspaces ", joint.stats.workspaces,
            ", buffers ", joint.stats.buffers, ", forwards ", joint.stats.forwards,
            ", reuses ", joint.stats.reuses
          firstEval = false
        return Gscalar(joint[i]).sval
    raiseValueError("training test requested an unplanned expression")

  when defined(graphPlanMemory):
    if memoryOnly:
      let dh = value(graph.deltaHamiltonian)
      echo "  training memory-only deltaHamiltonian ", dh,
        ", loss ", Gscalar(joint[1]).sval
      check dh > 1e-7
      check joint.stats.forwards > 0
      joint.clear
      joint.reportMemory("clear-complete")
      return

  var fdPlan: GraphPlan
  var fdx: Gscalar

  proc clearFd() =
    if fdPlan != nil:
      echo "  training FD plan ", (if fdx == graph.lossExpr: "loss" else: "rho gradient"),
        " arenaBytes ", fdPlan.stats.arenaBytes,
        ", workspaceBytes ", fdPlan.stats.workspaceBytes,
        ", forwards ", fdPlan.stats.forwards, ", reuses ", fdPlan.stats.reuses
      fdPlan.clear
      fdPlan = nil
      fdx = nil

  proc sampleValue(x: Gscalar): float =
    if joint == nil:
      return x.eval.sval
    if fdx != x:
      clearFd()
      fdx = x
      fdPlan = plan(x)
    discard fdPlan.eval
    Gscalar(fdPlan[0]).sval

  defer: clearFd()

  proc checkUnchangedPlan() =
    if joint != nil:
      discard joint.eval
      let before = joint.stats
      discard joint.eval
      echo "  training plan arenaBytes ", joint.stats.arenaBytes,
        ", workspaceBytes ", joint.stats.workspaceBytes,
        ", workspaces ", joint.stats.workspaces,
        ", forwards ", joint.stats.forwards, ", reuses ", joint.stats.reuses,
        ", unchanged forwards ", joint.stats.forwards-before.forwards,
        ", unchanged reuses ", joint.stats.reuses-before.reuses
      check joint.stats.arenaBytes > 0
      check joint.stats.forwards > 0
      check joint.stats.reuses > 0
      check joint.stats.forwards == before.forwards
      check joint.stats.reuses == before.reuses
      check joint.stats.arenaBytes == before.arenaBytes
      check joint.stats.workspaceBytes == before.workspaceBytes
      check joint.stats.workspaces == before.workspaces

  proc checkDerivative(name: string, f, derivative, param: Gscalar, step: float) =
    let at = param.sval
    let an = value(derivative)
    proc val(x: float): float =
      param.update x
      sampleValue(f)
    var fd, err: float
    ndiff(fd,err,val,at,step,ordMax=3)
    param.update at
    echo "  ", name, " derivative ", an, ", residual ", abs(an-fd), ", FD estimate ", err
    check abs(an-fd) < 2e-6*max(1.0,abs(fd))

  suite "stout trajectory training":
    test "training exposes flow rho beside integrator parameters":
      check graph.learnedParameters.len == 3
      check graph.learnedParameters[0].name == "dt"
      check graph.learnedParameters[1].name == "lambda"
      check graph.learnedParameters[2].name == "flowRho"
      check graph.learnedParameters[2].node.nodeKey == rho.nodeKey

    test "positive Hamiltonian change exercises transformed force derivatives":
      let dh = value(graph.deltaHamiltonian)
      echo "  training deltaHamiltonian ", dh, ", loss ", value(graph.lossExpr)
      check dh > 1e-7
      for lp in graph.learnedParameters:
        checkDerivative("trajectory " & lp.name,graph.lossExpr,lp.gradientExpr,lp.node,0.0003)
      checkUnchangedPlan()

    test "flow rho and step size retain their mixed derivative":
      checkDerivative("trajectory rho/dt",dr,mixed1Expr,dt,0.0003)
      let mixed1 = value(mixed1Expr)
      let mixed2 = value(mixed2Expr)
      echo "  mixed rho/dt ", mixed1, ", dt/rho ", mixed2
      check abs(mixed1-mixed2) < 2e-7*max(1.0,abs(mixed1))
      checkUnchangedPlan()

    test "the optimizer updates the exposed flow parameter":
      clearFd()
      let before = rho.sval
      var learner = initTrainingState(graph,0.0)
      var gradients = newSeq[float](graph.learnedParameters.len)
      for i, lp in graph.learnedParameters:
        gradients[i] = value(lp.gradientExpr)
      learner.trainStep(cfg,1,gradients)
      echo "  flow rho before ", before, ", after ", rho.sval
      check rho.sval != before

when isMainModule:
  runStoutTrainingTests(@[4,4,4,4],5.4)
