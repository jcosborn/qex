## Fused learned stage, with ordinary graph inputs and a composed derivative replica.
import qex
import ../graph/[core, scalar, gauge, multi]
import ../graph/gauge/types
import ../graph/support/op
import ../graph/hmcgauge/flow as flowActionOps
import flow, model, expr

type GnnftForward[T: SomeFloat] = ref object of Ggauge
  stage: int
  work: NnftStage[T]

# Every node's inputs are the flow input, the parameter leaves in
# NnftLayer.inputs order, then the node's own trailing operands.
template params(v: Gvalue; trailing: int): untyped =
  v.inputs.toOpenArray(1,v.inputs.len-1-trailing)

template registerValues(T: typedesc) =
  method hasStorage(x: GnnftForward[T]): bool =
    x.work != nil and (procCall Ggauge(x).hasStorage)

  method releaseWork(x: GnnftForward[T]) =
    x.work = nil

  method releaseStorage(x: GnnftForward[T]) =
    x.releaseWork
    procCall Ggauge(x).releaseStorage

  method bufferProto(x: GnnftForward[T]): Gvalue =
    # Pullbacks retain all stage activations, beyond the gauge output buffer.
    nil

  method newOneOf(x: GnnftForward[T]): Gvalue =
    GnnftForward[T](runtime:x.runtime,gval:x.gaugeNodeLike.gval,stage:x.stage).assignStableNodeId

registerValues(float32)
registerValues(float64)

proc stageForward[T: SomeFloat](v: Gvalue) =
  let z = GnnftForward[T](v)
  let w = Ggauge(v.inputs[0])
  let q = layer[T](params(v,0)).numericalParams
  if z.work == nil:
    z.work = newNnftStage(w.gval,q,z.stage)
  z.work.output = z.gval
  discard evalStage(z.work,q,w.gval)

proc internalBackward(zb, z: Gvalue; i: int; input: Gvalue): Gvalue =
  raiseUnsupportedPath("learned stage backward", "use the grouped learned outputs")

proc logdetForward[T: SomeFloat](v: Gvalue) =
  Gscalar(v).sval = GnnftForward[T](v.inputs[0]).work.logdet

proc gradientForward[T: SomeFloat](v: Gvalue) =
  # Trailing operands: the field cotangent, the logdet cotangent, the base.
  let base = GnnftForward[T](v.inputs[^1])
  let p = layer[T](params(v,3))
  let bg = Ggauge(v.inputs[^3])
  let bl = Gscalar(v.inputs[^2])
  stageVjp(base.work,p.numericalParams,bg.gval,bl.sval,Ggauge(v).gval)

proc gradientInputView(v: Gvalue; mode: InputWalkMode; visit: GnodeVisit) =
  if mode == iwmEval:
    visit v.inputs[^1]
  else:
    for i in 0..<v.inputs.len-3: visit v.inputs[i]
  visit v.inputs[^3]
  visit v.inputs[^2]

proc stageScore[T: SomeFloat](values: openArray[Gvalue]; bg: Ggauge; bl: Gscalar; stage: int): tuple[slots: seq[Gvalue], score: Gscalar] =
  # values are the flow input and the parameter leaves. Partial pullbacks
  # separate primal and seed slots, including shared parameters.
  for v in values: result.slots.add slotVar(v)
  result.slots.add [Gvalue(slotVar(bg)),Gvalue(slotVar(bl))]
  let s = result.slots
  let n = values.len
  let replica = learnedStageGraph(Ggauge(s[0]),layer[T](s.toOpenArray(1,n-1)),stage)
  result.score = redot(Ggauge(s[n]),replica.Wnew)+Gscalar(s[n+1])*replica.lj

proc stageGradient[T: SomeFloat](w: Ggauge; p: NnftLayer[T]; bg: Ggauge; bl: Gscalar; base: GnnftForward[T]): Ggauge =
  let stage = base.stage
  proc backward(zb,z: Gvalue; i: int; input: Gvalue): Gvalue =
    let rep = stageScore[T](z.inputs.toOpenArray(0,z.inputs.len-4),Ggauge(z.inputs[^3]),Gscalar(z.inputs[^2]),stage)
    let dw = grad(rep.score,Ggauge(rep.slots[0]))
    let up = requireUpstream(zb,"learned stage gradient backward",Ggauge)
    gradSeeded(dw,rep.slots[i],up)
  var args = @[Gvalue(w)]
  args.add p.inputs
  args.add [Gvalue(bg),Gvalue(bl),Gvalue(base)]
  graphNode(w.gaugeNodeLike,args,
    Gfunc(bufferMode:bmFull,forward:gradientForward[T],backward:backward,inputView:gradientInputView,
      name:"learnedStageGrad"),"learnedStageGrad")

proc parentInputView(v: Gvalue; mode: InputWalkMode; visit: GnodeVisit) =
  for i in 0..<v.inputs.len-2: visit v.inputs[i]

proc parentForward(v: Gvalue) = discard

proc outputInputView(v: Gvalue; mode: InputWalkMode; visit: GnodeVisit) =
  visit v.inputs[0]
  if mode != iwmBackward: visit v.inputs[1]

proc parentBackward[T: SomeFloat](zb,z: Gvalue; i: int; input: Gvalue): Gvalue =
  # Trailing operands: the base and its logdet value.
  let w = Ggauge(z.inputs[0])
  let p = layer[T](params(z,2))
  let base = GnnftForward[T](z.inputs[^2])
  let up = Gmulti(rootedUpstream(zb,z))
  if i == 0:
    return stageGradient(w,p,Ggauge(up[0]),Gscalar(up[1]),base)
  let rep = stageScore[T](z.inputs.toOpenArray(0,z.inputs.len-3),Ggauge(up[0]),Gscalar(up[1]),base.stage)
  grad(rep.score,rep.slots[i])

proc updateViewForward(v: Gvalue) =
  Ggauge(v).gval = Ggauge(v.inputs[1]).gval

proc updateViewBackward(zb,z: Gvalue; i: int; input: Gvalue): Gvalue =
  let parent = Gmulti(z.inputs[0])
  let up = Ggauge(rootedUpstream(zb,z))
  multiValues("learned field cotangents",up,parent.inputs[^1].zeroLike)

proc logdetViewForward(v: Gvalue) =
  Gscalar(v).sval = Gscalar(v.inputs[1]).sval

proc logdetViewBackward(zb,z: Gvalue; i: int; input: Gvalue): Gvalue =
  let parent = Gmulti(z.inputs[0])
  let up = Gscalar(rootedUpstream(zb,z))
  multiValues("learned logdet cotangents",parent.inputs[^2].zeroLike,up)

proc updateLogdet(v: Gvalue): tuple[ld,via:Gvalue] =
  let parent = Gmulti(v.inputs[0])
  let w = parent.inputs[0]
  for i in 1..<parent.inputs.len-2:
    if parent.inputs[i].reaches(w,iwmBackward):
      raiseValueError("learned stage logDetJ requires parameters independent of the flow input")
  (v.inputs[2],w)

proc learnedStage*[T: SomeFloat](w: Ggauge; p: NnftLayer[T]; stage: int): tuple[Wnew:Ggauge,lj:Gscalar] =
  w.gval.requireNnftGauge
  if stage < 0:
    raiseValueError("learned stage index must be nonnegative")
  p.requireLayer
  var args = @[Gvalue(w)]
  args.add p.inputs
  discard sharedGraphRuntime(args,"learnedStage")
  let base = graphNode(
    GnnftForward[T](runtime:w.runtime,gval:w.gaugeNodeLike.gval,stage:stage),args,
    Gfunc(bufferMode:bmFull,forward:stageForward[T],backward:internalBackward,name:"learnedStageForward"),"learnedStageForward")
  let ld = graphNode(scalarNodeLike(w),@[Gvalue(base)],
    Gfunc(bufferMode:bmFull,forward:logdetForward[T],backward:internalBackward,name:"learnedStageLogValue"),"learnedStageLogValue")
  args.add [Gvalue(base),Gvalue(ld)]
  let parent = newMultiStructureNode(@[Gvalue(base),Gvalue(ld)],args,
    Gfunc(bufferMode:bmFull,forward:parentForward,backward:parentBackward[T],inputView:parentInputView,name:"learnedStage"),"learnedStage")
  result.lj = graphNode(scalarNodeLike(w),@[Gvalue(parent),Gvalue(ld)],
    Gfunc(bufferMode:bmFull,forward:logdetViewForward,backward:logdetViewBackward,inputView:outputInputView,name:"learnedStageLogDet"),"learnedStageLogDet")
  result.Wnew = graphNode(Ggauge(runtime:w.runtime,gval:base.gval),
    @[Gvalue(parent),Gvalue(base),Gvalue(result.lj)],
    Gfunc(bufferMode:bmAlias,aliasInputs: @[1],forward:updateViewForward,backward:updateViewBackward,inputView:outputInputView,logdet:updateLogdet,name:"learnedStageUpdate"),"learnedStageUpdate")

proc learnedFlow*[T: SomeFloat](v: Ggauge; p: NnftModel[T]): Ggauge =
  result = v
  for s in 0..<p.len: result = learnedStage(result,p[s],s).Wnew

proc learnedAction*[T: SomeFloat](gc: Gactcoeff; p: NnftModel[T]): flowActionOps.FlowAction =
  flowActionOps.flowAction(gc,proc(v: Ggauge): Ggauge = learnedFlow(v,p))
