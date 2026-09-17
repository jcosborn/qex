## Explicit application comparison; compare.py supplies all inputs and checks outputs.
import base/globals
setDefaultNc(1)
setVLENmax(4)

import qex
import io/[arrays, arrayfields]
import ../../graph/[core, scalar, gauge, plan]
import ../[flow, model, graph, io]
import ../../graph/hmcgauge/[integrator, trajectory, config]
import std/[json, math]

qexInit()
let input = strParam("input", "")
let output = strParam("output", "")
let precision = strParam("precision", "double")
if input.len == 0 or output.len == 0 or precision notin ["single", "double"]:
  raise newException(ValueError, "compare requires -input:<directory> -output:<directory> -precision:single|double")
let man = readManifest(input)

proc write(name: string; data: openArray[float64]) =
  arrays.writeArray(output,man["results"][name],data,write = myRank == 0)

proc save(g: seq[DLatticeColorMatrixV]; reName, imName: string) =
  let lo = g[0].l
  let re = @[lo.RealD(),lo.RealD()]
  let im = @[lo.RealD(),lo.RealD()]
  threads:
    for d in 0..1:
      for x in re[d]:
        re[d][x] := g[d][x][0,0].re
        im[d][x] := g[d][x][0,0].im
  if reName.len > 0: saveFields(re,output,man["results"][reName])
  saveFields(im,output,man["results"][imName])

proc run[T: SomeFloat](precision: typedesc[T]) =
  let dims = arrayShape(man["arrays"]["theta"])
  let lo = newLayout(@[dims[1],dims[2]])
  let g = lo.newGauge
  let p = lo.newGauge
  loadAngles(g,input,man["arrays"]["theta"])
  let pi = @[lo.RealD(),lo.RealD()]
  loadFields(pi,input,man["arrays"]["initial_momentum"])
  threads:
    for d in 0..1:
      for x in p[d]:
        p[d][x][0,0].re := 0
        p[d][x][0,0].im := pi[d][x]
  let rt = initGraphRuntime()
  let params = toNnftModel(rt,loadNnft(input,T))
  let c = man["config"]
  let action = learnedAction(actWilson(scalar.toGvalue(rt,c["beta"].getFloat)),params)
  let cfg = RunConfig(dt:c["dt"].getFloat,gsteps:c["steps"].getInt,
    integratorCoeffs:parseIntegratorCoeffs(ik2MNp,@[c["lambda"].getFloat]))
  let traj = buildTrajectoryGraph(rt,g,p,action.action,cfg,buildTraining=false)
  let v = traj.initialState.gauge
  let u = action.flow(v)
  let ld = logDetJ(u,v)
  let evals = plan(u,ld,traj.initialState.gaugeAction,traj.mdForces[0],
    traj.finalState.gauge,traj.finalState.momentum,traj.initialState.hamiltonian,
    traj.finalState.hamiltonian,traj.deltaHamiltonian)
  defer: evals.clear
  let values = evals.eval
  save(Ggauge(values[0]).gval,"flow_re","flow_im")
  write("logdet",[Gscalar(values[1]).sval])
  write("action",[Gscalar(values[2]).sval])
  save(Ggauge(values[3]).gval,"","force")
  save(Ggauge(values[4]).gval,"proposal_re","proposal_im")
  save(Ggauge(values[5]).gval,"","momentum")
  write("H0",[Gscalar(values[6]).sval])
  write("H1",[Gscalar(values[7]).sval])
  let dh = Gscalar(values[8]).sval
  let probability = min(1.0,exp(-dh))
  write("deltaH",[dh])
  write("probability",[probability])
  let accepted = c["uniform"].getFloat < probability
  write("accepted",[if accepted: 1.0 else: 0.0])
  if accepted: traj.commitAcceptedTrajectory(Ggauge(values[4]).gaugeSnapshot)
  save(traj.initialState.gauge.gval,"committed_re","committed_im")

if precision == "single": run(float32)
else: run(float64)
qexFinalize()
