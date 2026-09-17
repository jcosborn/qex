## The stage networks as NN graph functions over parameter leaves.
import ../graph/[core, nn]
import flow
import ../../nn as numeric
import std/math

type
  NnftLayer*[T: SomeFloat] = object
    ## One stage's leaves: each convolution's weights [out,in,k0,k1] and bias
    ## [out], then the output scales.
    weights*, biases*: seq[Garray[T]]
    scale*: Garray[T]
  NnftModel*[T: SomeFloat] = seq[NnftLayer[T]]

proc toNnftModel*[T: SomeFloat](rt: GraphRuntime; p: NnftParams[T]): NnftModel[T] =
  for net in p:
    net.requireNet
    var layer: NnftLayer[T]
    for q in net.layers:
      layer.weights.add toGarray(rt,q.weights,@[q.cout,q.cin] & q.kernel)
      layer.biases.add toGarray(rt,q.bias,[q.cout])
    layer.scale = toGarray(rt,net.scale,[nnftCoefs])
    result.add layer

proc update*[T: SomeFloat](model: NnftModel[T]; p: NnftParams[T]) =
  if model.len != p.len:
    raiseValueError("learned model stage count differs from the parameters")
  for s, net in p:
    net.requireNet
    if net.layers.len != model[s].weights.len:
      raiseValueError("learned model convolution count differs from the parameters")
    for l, q in net.layers:
      model[s].weights[l].update(q.weights)
      model[s].biases[l].update(q.bias)
    model[s].scale.update(net.scale)

proc inputs*[T: SomeFloat](p: NnftLayer[T]): seq[Gvalue] =
  ## Each convolution's weights then bias, then the scales.
  for l in 0..<p.weights.len:
    result.add [Gvalue(p.weights[l]),Gvalue(p.biases[l])]
  result.add Gvalue(p.scale)

proc layer*[T: SomeFloat](values: openArray[Gvalue]): NnftLayer[T] =
  ## Restore the stage's parameter slots at an erased graph boundary.
  for l in 0..<values.len div 2:
    result.weights.add Garray[T](values[2*l])
    result.biases.add Garray[T](values[2*l+1])
  result.scale = Garray[T](values[^1])

proc requireLayer*[T: SomeFloat](p: NnftLayer[T]) =
  if p.weights.len == 0 or p.biases.len != p.weights.len:
    raiseValueError("learned layer requires a bias for every convolution")
  var cin = nnftFeatures
  for l in 0..<p.weights.len:
    let w = p.weights[l].shape
    if w.len != 4 or w[1] != cin or p.biases[l].shape != @[w[0]]:
      raiseValueError("learned layer weights must be [out,in,k0,k1] chained from " & $nnftFeatures & " channels, with a bias per output channel")
    cin = w[0]
  if cin != nnftCoefs or p.scale.shape != @[nnftCoefs]:
    raiseValueError("learned layer must end in " & $nnftCoefs & " channels and scales")

proc numericalParams*[T: SomeFloat](p: NnftLayer[T]): NnftNet[T] =
  ## Copy the current array payloads; the kernels follow the weight shapes.
  for l in 0..<p.weights.len:
    let w = p.weights[l].shape
    result.layers.add convParams(w[1],w[0],w[2..^1],p.weights[l].data,p.biases[l].data)
  result.scale = p.scale.data

proc network*[T: SomeFloat](x: Greal[T]; p: NnftLayer[T]): Greal[T] =
  p.requireLayer
  var h = x
  for l in 0..<p.weights.len:
    h = bias(conv(h,p.weights[l]),p.biases[l])
    if l < p.weights.high: h = gelu(h)
  let z = scale(h,p.scale)
  divide(divide(arctan(z),T(PI)),T(3))
