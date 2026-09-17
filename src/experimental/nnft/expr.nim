## Learned U(1) stage as ordinary field, path, and neural graph functions.
import ../graph/[core, scalar, gauge]
import ../graph/nn as neural
import ../graph/gauge/[types, basic_ops, field_ops, transport, cfield, rfield]
import flow, model

const loops = [@[1,2,-1,-2], @[1,1,2,-1,-1,-2], @[1,2,2,-1,-2,-2]]

proc openPaths(dir: int): seq[seq[int]] =
  ## Rotate the loop to its active link. L = W^a Q gives D = Q^(-a).
  for k in 0..<6:
    let path = loops[if k < 2: 0 else: dir+1]
    var pos = stapleAnchors[dir][k]
    var at = -1
    for j,p in path:
      let mu = abs(p)-1
      if p < 0: dec pos[mu]
      if mu == dir and pos == [0,0]: at = j
      if p > 0: inc pos[mu]
    doAssert at >= 0
    var tail: seq[int]
    for j in 1..<path.len: tail.add path[(at+j) mod path.len]
    if path[at] > 0:
      var rev: seq[int]
      for j in countdown(tail.high,0): rev.add -tail[j]
      result.add rev
    else:
      result.add tail

proc flowViewForward(v: Gvalue) =
  Ggauge(v).gval = Ggauge(v.inputs[0]).gval

proc flowViewBackward(zb,z: Gvalue; i: int; input: Gvalue): Gvalue =
  rootedUpstream(zb,z)

proc flowViewInputs(v: Gvalue; mode: InputWalkMode; visit: GnodeVisit) =
  visit v.inputs[0]
  if mode == iwmReachable: visit v.inputs[2]

proc flowViewLogdet(v: Gvalue): tuple[ld,via:Gvalue] =
  let w = v.inputs[1]
  for j in 3..<v.inputs.len:
    if v.inputs[j].reaches(w,iwmBackward):
      raiseValueError("learned stage logDetJ requires parameters independent of the flow input")
  (v.inputs[2],w)

let flowView = Gfunc(bufferMode:bmAlias,aliasInputs: @[0],forward:flowViewForward,backward:flowViewBackward,
  inputView:flowViewInputs,logdet:flowViewLogdet,name:"learnedStageGraph")

proc learnedStageGraph*[T: SomeFloat](w: Ggauge; p: NnftLayer[T]; stage: int): tuple[Wnew:Ggauge,lj:Gscalar] =
  w.gval.requireNnftGauge
  let masks = nnftMasks(w.gval[0].l,stage)
  let dir = stageClass(stage).dir
  var paths: seq[seq[int]] = @loops
  paths.add openPaths(dir)
  let products = lineProducts(w,paths)
  var si,co: array[3,Greal[T]]
  for k in 0..<3:
    let loop = trace(products[k])
    let im = rfield.imag(loop,T)
    let re = rfield.real(loop,T)
    let mask = toGfield(w.runtime,masks.featureMasks[k])
    si[k] = maskedCopy(im,Greal[T](im.zeroLike),mask)
    co[k] = maskedCopy(re,Greal[T](re.oneLike),mask)
  let features = concat([si[0],co[0],si[1],si[2],co[1],co[2]])
  let coef = network(features,p)
  var ds = -rfield.scale(channel(coef,coefSlots[dir][0]),products[3])
  for k in 1..<6:
    ds = ds-rfield.scale(channel(coef,coefSlots[dir][k]),products[3+k])
  let active = neural.toGvalue(w.runtime,@[masks.active])
  ds = rfield.scale(active,ds)
  let m = linkField(w,dir)*ds.adj
  let field = exp(injectLink(projTAH(m),dir,w))*w
  let diag = rfield.real(trace(m),float64)+1.0
  result.lj = sum(ln(clipMin(diag,nnftJacFloor)))
  var args = @[Gvalue(field),Gvalue(w),Gvalue(result.lj)]
  args.add p.inputs
  result.Wnew = graphNode(w.gaugeNodeLike,args,flowView,"learnedStageGraph")

proc learnedFlowGraph*[T: SomeFloat](w: Ggauge; p: NnftModel[T]): Ggauge =
  result = w
  for s in 0..<p.len:
    result = learnedStageGraph(result,p[s],s).Wnew
