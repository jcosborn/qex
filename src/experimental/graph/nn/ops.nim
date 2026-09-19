## Separate graph operations over shared QEX field storage.
import qex
import ../../../nn as numeric
import ../[core, scalar]
import ../support/op
import types
import std/[math, tables, sequtils]

type
  ConvWork[T: SomeFloat] = ref object of Gwork
    ws: numeric.ConvWorkspace[T]
  Gconv[T: SomeFloat] = ref object of Greal[T]
    ## Numerical calls rebind every halo field before use.
    params: ConvParams[T]
    key: string
    ws: numeric.ConvWorkspace[T]
  GconvWeight[T: SomeFloat] = ref object of Garray[T]
    params: ConvParams[T]
    key: string
    ws: numeric.ConvWorkspace[T]


proc requireTaps(offsets: seq[seq[int32]], lo: Layout[VLEN]) =
  ## The numerical workspace rejects these as well, but a node may be built and
  ## never evaluated; diagnose at construction.
  for off in offsets:
    for d in 0..<off.len:
      if abs(off[d]) >= lo.physGeom[d]:
        raiseValueError("convolution tap exceeds halo periodic extent")

# Explicit instances preserve erased dispatch and avoid inferring a precision
# through the numerical field alias. Both use the same operator definitions.
template realOps(T: typedesc) {.dirty.} =
  proc `+`*(x, y: Greal[T]): Greal[T]
  proc `*`*(x, y: Greal[T]): Greal[T]
  proc `/`*(x, y: Greal[T]): Greal[T]
  proc `*`*(x: Gscalar, y: Greal[T]): Greal[T]
  proc `*`*(a: T, x: Greal[T]): Greal[T]
  proc `+`*(x: Greal[T], a: T): Greal[T]
  proc `-`*(x: Greal[T]): Greal[T]
  proc divide*(x: Greal[T], d: T): Greal[T]
  proc exp*(x: Greal[T]): Greal[T]
  proc erfc*(x: Greal[T]): Greal[T]
  proc redot*(x, y: Greal[T]): Gscalar
  proc conv*(x: Greal[T], w: Garray[T], b: Garray[T] = nil): Greal[T]
  proc convTranspose*(x: Greal[T], w: Garray[T]): Greal[T]
  proc convWeightVjp*(x, dy: Greal[T], shape: openArray[int]): Garray[T]
  proc channelSum*(x: Greal[T]): Garray[T]
  proc broadcast*(x: Garray[T], proto: Greal[T]): Greal[T]
  proc scale*(x: Greal[T], s: Garray[T]): Greal[T]
  proc `+`*(x, y: Garray[T]): Garray[T]
  proc `*`*(x: Gscalar, y: Garray[T]): Garray[T]
  proc redot*(x, y: Garray[T]): Gscalar

  proc `+`*(x, y: Greal[T]): Greal[T] =
    x.requireSameFieldShape(y,"real addition")
    proc forward(v: Gvalue) =
      let z = Greal[T](v)
      let a = Greal[T](v.inputs[0])
      let b = Greal[T](v.inputs[1])
      threads:
        for c in 0..<z.fval.len: z.fval[c] := a.fval[c]+b.fval[c]
    proc backward(zb,z: Gvalue, i: int, input: Gvalue): Gvalue = rootedUpstream(zb,z)
    graphNode(x.realNodeLike,@[Gvalue(x),Gvalue(y)],Gfunc(bufferMode:bmFull,inplace: @[0,1],forward:forward,backward:backward,name:"realAdd"),"realAdd")

  proc `*`*(x, y: Greal[T]): Greal[T] =
    x.requireSameFieldShape(y,"real multiplication")
    proc forward(v: Gvalue) =
      let z = Greal[T](v)
      let a = Greal[T](v.inputs[0])
      let b = Greal[T](v.inputs[1])
      threads:
        for c in 0..<z.fval.len: z.fval[c] := a.fval[c]*b.fval[c]
    proc backward(zb,z: Gvalue, i: int, input: Gvalue): Gvalue =
      Greal[T](rootedUpstream(zb,z))*Greal[T](z.inputs[1-i])
    graphNode(x.realNodeLike,@[Gvalue(x),Gvalue(y)],Gfunc(bufferMode:bmFull,inplace: @[0,1],forward:forward,backward:backward,name:"realMul"),"realMul")

  proc `/`*(x, y: Greal[T]): Greal[T] =
    x.requireSameFieldShape(y,"real division")
    proc forward(v: Gvalue) =
      numeric.divide(Greal[T](v).fval,Greal[T](v.inputs[0]).fval,Greal[T](v.inputs[1]).fval)
    proc backward(zb,z: Gvalue, i: int, input: Gvalue): Gvalue =
      let up = Greal[T](rootedUpstream(zb,z))
      let b = Greal[T](z.inputs[1])
      if i == 0: up/b
      else: -(up*Greal[T](z.inputs[0])/(b*b))
    graphNode(x.realNodeLike,@[Gvalue(x),Gvalue(y)],Gfunc(bufferMode:bmFull,inplace: @[0,1],forward:forward,backward:backward,name:"realDiv"),"realDiv")

  proc `*`*(x: Gscalar, y: Greal[T]): Greal[T] =
    proc forward(v: Gvalue) =
      let z = Greal[T](v)
      let a = T(Gscalar(v.inputs[0]).sval)
      let b = Greal[T](v.inputs[1])
      threads:
        for c in 0..<z.fval.len: z.fval[c] := a*b.fval[c]
    proc backward(zb,z: Gvalue, i: int, input: Gvalue): Gvalue =
      let up = Greal[T](rootedUpstream(zb,z))
      if i == 0: Gvalue(redot(up,Greal[T](z.inputs[1])))
      else: Gvalue(Gscalar(z.inputs[0])*up)
    graphNode(y.realNodeLike,@[Gvalue(x),Gvalue(y)],Gfunc(bufferMode:bmFull,inplace: @[1],forward:forward,backward:backward,name:"realScale"),"realScale")

  proc `*`*(a: T, x: Greal[T]): Greal[T] = scalar.toGvalue(x.runtime,float(a))*x
  proc `*`*(x: Greal[T], a: T): Greal[T] = a*x
  proc `-`*(x: Greal[T]): Greal[T] = T(-1)*x
  proc `-`*(x,y: Greal[T]): Greal[T] = x+(-y)

  proc `+`*(x: Greal[T], a: T): Greal[T] =
    proc forward(v: Gvalue) =
      let z = Greal[T](v)
      let x = Greal[T](v.inputs[0])
      threads:
        for c in 0..<z.fval.len: z.fval[c] := x.fval[c]+a
    proc backward(zb,z: Gvalue, i: int, input: Gvalue): Gvalue = rootedUpstream(zb,z)
    graphNode(x.realNodeLike,@[Gvalue(x)],Gfunc(bufferMode:bmFull,inplace: @[0],forward:forward,backward:backward,name:"realOffset"),"realOffset")

  proc `+`*(a: T,x: Greal[T]): Greal[T] = x+a
  proc `-`*(x: Greal[T], a: T): Greal[T] = x+(-a)
  proc `-`*(a: T,x: Greal[T]): Greal[T] = (-x)+a

  proc divide*(x: Greal[T], d: T): Greal[T] =
    proc forward(v: Gvalue) = numeric.divide(Greal[T](v).fval,Greal[T](v.inputs[0]).fval,d)
    proc backward(zb,z: Gvalue, i: int, input: Gvalue): Gvalue = divide(Greal[T](rootedUpstream(zb,z)),d)
    graphNode(x.realNodeLike,@[Gvalue(x)],Gfunc(bufferMode:bmFull,inplace: @[0],forward:forward,backward:backward,name:"realDivide"),"realDivide")

  proc exp*(x: Greal[T]): Greal[T] =
    proc forward(v: Gvalue) = numeric.exp(Greal[T](v).fval,Greal[T](v.inputs[0]).fval)
    proc backward(zb,z: Gvalue, i: int, input: Gvalue): Gvalue =
      Greal[T](rootedUpstream(zb,z))*Greal[T](z)
    graphNode(x.realNodeLike,@[Gvalue(x)],Gfunc(bufferMode:bmFull,inplace: @[0],forward:forward,backward:backward,name:"realExp"),"realExp")

  proc erfc*(x: Greal[T]): Greal[T] =
    proc forward(v: Gvalue) = numeric.erfc(Greal[T](v).fval,Greal[T](v.inputs[0]).fval)
    proc backward(zb,z: Gvalue, i: int, input: Gvalue): Gvalue =
      let a = Greal[T](z.inputs[0])
      Greal[T](rootedUpstream(zb,z))*(T(-1.12837916709551257389615890312154517)*exp(-(a*a)))
    graphNode(x.realNodeLike,@[Gvalue(x)],Gfunc(bufferMode:bmFull,inplace: @[0],forward:forward,backward:backward,name:"realErfc"),"realErfc")

  proc ln*(x: Greal[T]): Greal[T] =
    proc forward(v: Gvalue) = numeric.ln(Greal[T](v).fval,Greal[T](v.inputs[0]).fval)
    proc backward(zb,z: Gvalue, i: int, input: Gvalue): Gvalue =
      Greal[T](rootedUpstream(zb,z))/Greal[T](z.inputs[0])
    graphNode(x.realNodeLike,@[Gvalue(x)],Gfunc(bufferMode:bmFull,inplace: @[0],forward:forward,backward:backward,name:"realLog"),"realLog")

  proc gelu*(x: Greal[T]): Greal[T] =
    proc forward(v: Gvalue) = numeric.gelu(Greal[T](v).fval,Greal[T](v.inputs[0]).fval)
    proc backward(zb,z: Gvalue, i: int, input: Gvalue): Gvalue =
      let a = Greal[T](z.inputs[0])
      let phi = T(0.5)*erfc(T(-numeric.sqrt_1_2)*a)
      let pdf = T(numeric.sqrt_1_2pi)*exp((T(-0.5)*a)*a)
      Greal[T](rootedUpstream(zb,z))*(phi+a*pdf)
    graphNode(x.realNodeLike,@[Gvalue(x)],Gfunc(bufferMode:bmFull,inplace: @[0],forward:forward,backward:backward,name:"gelu"),"gelu")

  proc arctan*(x: Greal[T]): Greal[T] =
    proc forward(v: Gvalue) = numeric.arctan(Greal[T](v).fval,Greal[T](v.inputs[0]).fval)
    proc backward(zb,z: Gvalue, i: int, input: Gvalue): Gvalue =
      let a = Greal[T](z.inputs[0])
      Greal[T](rootedUpstream(zb,z))/(T(1)+a*a)
    graphNode(x.realNodeLike,@[Gvalue(x)],Gfunc(bufferMode:bmFull,inplace: @[0],forward:forward,backward:backward,name:"arctan"),"arctan")

  proc clipSlope(x: Greal[T], floor: T): Greal[T] =
    proc forward(v: Gvalue) = numeric.clipSlope(Greal[T](v).fval,Greal[T](v.inputs[0]).fval,floor)
    proc backward(zb,z: Gvalue, i: int, input: Gvalue): Gvalue = z.inputs[0].zeroLike
    graphNode(x.realNodeLike,@[Gvalue(x)],Gfunc(bufferMode:bmFull,inplace: @[0],forward:forward,backward:backward,name:"clipSlope"),"clipSlope")

  proc clipMin*(x: Greal[T], floor: T): Greal[T] =
    ## JAX maximum convention: the derivative at equality is one half.
    proc forward(v: Gvalue) = numeric.clipMin(Greal[T](v).fval,Greal[T](v.inputs[0]).fval,floor)
    proc backward(zb,z: Gvalue, i: int, input: Gvalue): Gvalue =
      Greal[T](rootedUpstream(zb,z))*clipSlope(Greal[T](z.inputs[0]),floor)
    graphNode(x.realNodeLike,@[Gvalue(x)],Gfunc(bufferMode:bmFull,inplace: @[0],forward:forward,backward:backward,name:"clipMin"),"clipMin")

  proc convNode(x: Greal[T], channels: int, p: ConvParams[T]): Gconv[T] =
    requireTaps(p.offsets, x.fval[0].l)
    var fs = newSeq[RealField[T]](channels)
    for c in 0..<channels: fs[c] = x.fval[0].newShape
    var cfg = p
    cfg.weights = @[]
    cfg.bias = @[]
    Gconv[T](runtime:x.runtime,fval:fs,params:cfg,key:"conv " & $T & " " & $p.offsets).assignStableNodeId

  method newOneOf(x: Gconv[T]): Gvalue = convNode(Greal[T](x),x.fval.len,x.params)

  method releaseWork(x: Gconv[T]) = x.ws = nil

  method releaseStorage(x: Gconv[T]) =
    x.releaseWork
    procCall Greal[T](x).releaseStorage

  proc work(ws: var numeric.ConvWorkspace[T], rt: GraphRuntime, proto: RealField[T], p: ConvParams[T], tag: string): numeric.ConvWorkspace[T] =
    if rt.evalFrame != nil and rt.workFrame != nil:
      let key: GworkKey = (kind:tag,layout:cast[pointer](proto.l),fields:p.cin,order:0)
      var entry = ConvWork[T](rt.workFrame.getOrDefault(key))
      if entry == nil:
        entry = ConvWork[T](ws:numeric.convWorkspace(proto,p))
        for h in entry.ws.halo: entry.bytes += h.halo.bytes
        rt.workFrame[key] = entry
      ws = entry.ws
    elif ws == nil:
      ws = numeric.convWorkspace(proto,p)
    ws

  method newOneOf(x: GconvWeight[T]): Gvalue =
    GconvWeight[T](runtime:x.runtime,data:newSeq[T](x.data.len),shape:x.shape.mapIt(it),
      params:x.params,key:x.key).assignStableNodeId

  method releaseWork(x: GconvWeight[T]) = x.ws = nil
  method releaseStorage(x: GconvWeight[T]) = x.releaseWork

  proc convWeightVjp*(x, dy: Greal[T], shape: openArray[int]): Garray[T] =
    ## B(x,dy)_o,i,t = sum_s dy_o(s) x_i(s+offset_t).
    let lo = x.fval[0].l
    if dy.fval[0].l != lo or shape.len != lo.nDim:
      raiseValueError("convolution weight VJP field or kernel shape differs")
    var dims = @[dy.fval.len,x.fval.len]
    dims.add shape
    var n = 1
    for d in dims: n *= d
    let node = GconvWeight[T](runtime:x.runtime,data:newSeq[T](n),shape:dims)
    node.params = numeric.convParams(x.fval.len,dy.fval.len,shape,node.data)
    node.params.weights = @[]
    node.key = "conv " & $T & " " & $node.params.offsets
    requireTaps(node.params.offsets, lo)
    proc forward(v: Gvalue) =
      let z = GconvWeight[T](v)
      let a = Greal[T](v.inputs[0])
      let b = Greal[T](v.inputs[1])
      numeric.convWeightVjp(z.data,a.fval,b.fval,z.params,z.ws.work(z.runtime,a.fval[0],z.params,z.key))
    proc backward(zb,z: Gvalue, i: int, input: Gvalue): Gvalue =
      let up = Garray[T](rootedUpstream(zb,z))
      if i == 0: convTranspose(Greal[T](z.inputs[1]),up)
      else: conv(Greal[T](z.inputs[0]),up)
    graphNode(node,@[Gvalue(x),Gvalue(dy)],Gfunc(bufferMode:bmFull,forward:forward,backward:backward,name:"convWeightVjp"),"convWeightVjp")

  proc channelSum*(x: Greal[T]): Garray[T] =
    proc forward(v: Gvalue) =
      numeric.channelSum(Garray[T](v).data,Greal[T](v.inputs[0]).fval)
    proc backward(zb,z: Gvalue, i: int, input: Gvalue): Gvalue =
      broadcast(Garray[T](rootedUpstream(zb,z)),Greal[T](z.inputs[0]))
    let node = Garray[T](runtime:x.runtime,data:newSeq[T](x.fval.len),shape: @[x.fval.len])
    graphNode(node,@[Gvalue(x)],Gfunc(bufferMode:bmFull,forward:forward,backward:backward,name:"channelSum"),"channelSum")

  proc broadcast*(x: Garray[T], proto: Greal[T]): Greal[T] =
    if x.shape != @[proto.fval.len]:
      raiseValueError("broadcast shape differs from field channels")
    proc forward(v: Gvalue) =
      numeric.broadcast(Greal[T](v).fval,Garray[T](v.inputs[0]).data)
    proc backward(zb,z: Gvalue, i: int, input: Gvalue): Gvalue =
      channelSum(Greal[T](rootedUpstream(zb,z)))
    graphNode(proto.realNodeLike,@[Gvalue(x)],Gfunc(bufferMode:bmFull,forward:forward,backward:backward,name:"broadcast"),"broadcast")

  proc conv*(x: Greal[T], w: Garray[T], b: Garray[T] = nil): Greal[T] =
    if w.shape.len != x.fval[0].l.nDim+2 or w.shape[1] != x.fval.len:
      raiseValueError("graph convolution weight shape differs from input")
    if b != nil and b.shape != @[w.shape[0]]:
      raiseValueError("graph convolution bias shape differs from output channels")
    let p = numeric.convParams(w.shape[1],w.shape[0],w.shape[2..^1],w.data)
    proc forward(v: Gvalue) =
      let z = Gconv[T](v)
      var p = z.params
      p.weights = Garray[T](v.inputs[1]).data
      if v.inputs.len == 3: p.bias = Garray[T](v.inputs[2]).data
      numeric.conv(z.fval,Greal[T](v.inputs[0]).fval,p,z.ws.work(z.runtime,z.fval[0],p,z.key))
    proc backward(zb,z: Gvalue, i: int, input: Gvalue): Gvalue =
      let up = Greal[T](rootedUpstream(zb,z))
      let w = Garray[T](z.inputs[1])
      if i == 0: Gvalue(convTranspose(up,w))
      elif i == 1: Gvalue(convWeightVjp(Greal[T](z.inputs[0]),up,w.shape[2..^1]))
      else: Gvalue(channelSum(up))
    var args = @[Gvalue(x),Gvalue(w)]
    if b != nil: args.add Gvalue(b)
    graphNode(convNode(x,w.shape[0],p),args,Gfunc(bufferMode:bmFull,forward:forward,backward:backward,name:"conv"),"conv")

  proc convTranspose*(x: Greal[T], w: Garray[T]): Greal[T] =
    if w.shape.len != x.fval[0].l.nDim+2 or w.shape[0] != x.fval.len:
      raiseValueError("graph transpose convolution weight shape differs from input")
    let p = numeric.convParams(w.shape[1],w.shape[0],w.shape[2..^1],w.data)
    proc forward(v: Gvalue) =
      let z = Gconv[T](v)
      var p = z.params
      p.weights = Garray[T](v.inputs[1]).data
      numeric.convVjp(z.fval,Greal[T](v.inputs[0]).fval,p,z.ws.work(z.runtime,z.fval[0],p,z.key))
    proc backward(zb,z: Gvalue, i: int, input: Gvalue): Gvalue =
      let up = Greal[T](rootedUpstream(zb,z))
      let w = Garray[T](z.inputs[1])
      if i == 0: Gvalue(conv(up,w))
      else: Gvalue(convWeightVjp(up,Greal[T](z.inputs[0]),w.shape[2..^1]))
    graphNode(convNode(x,w.shape[1],p),@[Gvalue(x),Gvalue(w)],Gfunc(bufferMode:bmFull,forward:forward,backward:backward,name:"convTranspose"),"convTranspose")

  proc scale*(x: Greal[T], s: Garray[T]): Greal[T] =
    if s.shape != @[x.fval.len]: raiseValueError("graph scale shape differs from channels")
    proc forward(v: Gvalue) = numeric.scale(Greal[T](v).fval,Greal[T](v.inputs[0]).fval,Garray[T](v.inputs[1]).data)
    proc backward(zb,z: Gvalue, i: int, input: Gvalue): Gvalue =
      let up = Greal[T](rootedUpstream(zb,z))
      if i == 0: Gvalue(scale(up,Garray[T](z.inputs[1])))
      else: Gvalue(channelSum(up*Greal[T](z.inputs[0])))
    graphNode(x.realNodeLike,@[Gvalue(x),Gvalue(s)],Gfunc(bufferMode:bmFull,inplace: @[0],forward:forward,backward:backward,name:"channelScale"),"channelScale")

  proc bias*(x: Greal[T], b: Garray[T]): Greal[T] =
    if b.shape != @[x.fval.len]: raiseValueError("graph bias shape differs from channels")
    proc forward(v: Gvalue) = numeric.bias(Greal[T](v).fval,Greal[T](v.inputs[0]).fval,Garray[T](v.inputs[1]).data)
    proc backward(zb,z: Gvalue, i: int, input: Gvalue): Gvalue =
      let up = Greal[T](rootedUpstream(zb,z))
      if i == 0: Gvalue(up)
      else: Gvalue(channelSum(up))
    graphNode(x.realNodeLike,@[Gvalue(x),Gvalue(b)],Gfunc(bufferMode:bmFull,inplace: @[0],forward:forward,backward:backward,name:"bias"),"bias")

  proc select(x: Greal[T], mask: Gmask, sub = "all", complement = false): Greal[T] =
    if mask.fval.l != x.fval[0].l: raiseValueError("graph mask layout differs from input")
    discard x.fval[0].l.getSubset(sub)
    proc forward(v: Gvalue) =
      numeric.maskVjp(Greal[T](v).fval,Greal[T](v.inputs[0]).fval,sub,Gmask(v.inputs[1]).fval,complement)
    proc backward(zb,z: Gvalue, i: int, input: Gvalue): Gvalue =
      if i != 0: raiseUnsupportedPath("NN mask gradient","mask samples are discrete")
      select(Greal[T](rootedUpstream(zb,z)),Gmask(z.inputs[1]),sub,complement)
    graphNode(x.realNodeLike,@[Gvalue(x),Gvalue(mask)],Gfunc(bufferMode:bmFull,inplace: @[0],forward:forward,backward:backward,name:"select"),"select")

  proc maskedCopy*(x, fill: Greal[T], mask: Gmask, sub = "all"): Greal[T] =
    x.requireSameFieldShape(fill,"graph masked copy")
    if mask.fval.l != x.fval[0].l: raiseValueError("graph mask layout differs from input")
    discard x.fval[0].l.getSubset(sub)
    proc forward(v: Gvalue) =
      let z = Greal[T](v)
      z.fval.copyFieldStorage(Greal[T](v.inputs[1]).fval)
      numeric.maskedCopy(z.fval,Greal[T](v.inputs[0]).fval,Gmask(v.inputs[2]).fval,sub)
    proc backward(zb,z: Gvalue, i: int, input: Gvalue): Gvalue =
      if i == 2: raiseUnsupportedPath("NN mask gradient","mask samples are discrete")
      select(Greal[T](rootedUpstream(zb,z)),Gmask(z.inputs[2]),sub,i==1)
    graphNode(x.realNodeLike,@[Gvalue(x),Gvalue(fill),Gvalue(mask)],Gfunc(bufferMode:bmFull,forward:forward,backward:backward,name:"maskedCopy"),"maskedCopy")

  proc redot*(x, y: Greal[T]): Gscalar =
    x.requireSameFieldShape(y,"real dot product")
    proc forward(v: Gvalue) =
      Gscalar(v).sval = numeric.redot(Greal[T](v.inputs[0]).fval,Greal[T](v.inputs[1]).fval)
    proc backward(zb,z: Gvalue, i: int, input: Gvalue): Gvalue = bilinearBackward(zb,z,i,Greal[T])
    graphNode(scalarNodeLike(x),@[Gvalue(x),Gvalue(y)],Gfunc(bufferMode:bmFull,forward:forward,backward:backward,name:"realDot"),"realDot")

  proc sum*(x: Greal[T]): Gscalar =
    redot(x,Greal[T](x.oneLike))

  proc channel*(x: Greal[T], idx: int, count = 1): Greal[T]
  proc putChannels(x: Greal[T], first: int, proto: Greal[T]): Greal[T]

  proc channel*(x: Greal[T], idx: int, count = 1): Greal[T] =
    ## Channels idx..<idx+count as a value of their own.
    if idx < 0 or count <= 0 or idx+count > x.fval.len:
      raiseValueError("channel slice out of range")
    var fs = newSeq[RealField[T]](count)
    for c in 0..<count: fs[c] = x.fval[idx+c].newShape
    proc forward(v: Gvalue) =
      let z = Greal[T](v)
      let a = Greal[T](v.inputs[0])
      threads:
        for c in 0..<count: z.fval[c] := a.fval[idx+c]
    proc backward(zb,z: Gvalue, i: int, input: Gvalue): Gvalue =
      putChannels(Greal[T](rootedUpstream(zb,z)),idx,Greal[T](z.inputs[0]))
    graphNode(Greal[T](runtime:x.runtime,fval:fs),@[Gvalue(x)],Gfunc(bufferMode:bmFull,forward:forward,backward:backward,name:"channels"),"channels")

  proc putChannels(x: Greal[T], first: int, proto: Greal[T]): Greal[T] =
    let count = x.fval.len
    proc forward(v: Gvalue) =
      let z = Greal[T](v)
      let a = Greal[T](v.inputs[0])
      threads:
        for c in 0..<count: z.fval[first+c] := a.fval[c]
    proc backward(zb,z: Gvalue, i: int, input: Gvalue): Gvalue =
      channel(Greal[T](rootedUpstream(zb,z)),first,count)
    graphNode(proto.realNodeLike,@[Gvalue(x)],Gfunc(bufferMode:bmZero,forward:forward,backward:backward,name:"putChannels"),"putChannels")

  proc concat*(xs: openArray[Greal[T]]): Greal[T] =
    if xs.len == 0: raiseValueError("channel concatenation requires inputs")
    let lo = xs[0].fval[0].l
    var fs: seq[RealField[T]]
    var args: seq[Gvalue]
    var offsets: seq[int]
    for x in xs:
      if x.fval[0].l != lo: raiseValueError("channel concatenation layouts differ")
      offsets.add fs.len
      args.add Gvalue(x)
      for f in x.fval: fs.add f.newShape
    proc forward(v: Gvalue) =
      let z = Greal[T](v)
      threads:
        for i in 0..<v.inputs.len:
          let a = Greal[T](v.inputs[i])
          for c in 0..<a.fval.len: z.fval[offsets[i]+c] := a.fval[c]
    proc backward(zb,z: Gvalue, i: int, input: Gvalue): Gvalue =
      channel(Greal[T](rootedUpstream(zb,z)),offsets[i],Greal[T](z.inputs[i]).fval.len)
    graphNode(Greal[T](runtime:xs[0].runtime,fval:fs),args,Gfunc(bufferMode:bmFull,forward:forward,backward:backward,name:"concat"),"concat")

  proc `+`*(x, y: Garray[T]): Garray[T] =
    if not x.copyCompatible(y): raiseValueError("parameter addition shapes differ")
    proc forward(v: Gvalue) =
      let z = Garray[T](v)
      let a = Garray[T](v.inputs[0])
      let b = Garray[T](v.inputs[1])
      for i in 0..<z.data.len: z.data[i] = a.data[i]+b.data[i]
    proc backward(zb,z: Gvalue, i: int, input: Gvalue): Gvalue = rootedUpstream(zb,z)
    graphNode(x.arrayNodeLike,@[Gvalue(x),Gvalue(y)],Gfunc(bufferMode:bmFull,forward:forward,backward:backward,name:"arrayAdd"),"arrayAdd")

  proc redot*(x, y: Garray[T]): Gscalar =
    if not x.copyCompatible(y): raiseValueError("parameter dot product shapes differ")
    proc forward(v: Gvalue) =
      let a = Garray[T](v.inputs[0])
      let b = Garray[T](v.inputs[1])
      var s = 0.0
      for i in 0..<a.data.len: s += float(a.data[i])*float(b.data[i])
      Gscalar(v).sval = s
    proc backward(zb,z: Gvalue, i: int, input: Gvalue): Gvalue = bilinearBackward(zb,z,i,Garray[T])
    graphNode(scalarNodeLike(x),@[Gvalue(x),Gvalue(y)],Gfunc(bufferMode:bmFull,forward:forward,backward:backward,name:"arrayDot"),"arrayDot")

  proc `*`*(x: Gscalar, y: Garray[T]): Garray[T] =
    proc forward(v: Gvalue) =
      let z = Garray[T](v)
      let a = T(Gscalar(v.inputs[0]).sval)
      let b = Garray[T](v.inputs[1])
      for i in 0..<z.data.len: z.data[i] = a*b.data[i]
    proc backward(zb,z: Gvalue, i: int, input: Gvalue): Gvalue =
      let up = Garray[T](rootedUpstream(zb,z))
      if i == 0: Gvalue(redot(up,Garray[T](z.inputs[1])))
      else: Gvalue(Gscalar(z.inputs[0])*up)
    graphNode(y.arrayNodeLike,@[Gvalue(x),Gvalue(y)],Gfunc(bufferMode:bmFull,forward:forward,backward:backward,name:"arrayScale"),"arrayScale")

  method addLike*(prototype: Greal[T], x,y: Gvalue): Gvalue =
    if not (x of Greal[T]) or not (y of Greal[T]): raiseValueError("real addition precision differs")
    Greal[T](x)+Greal[T](y)
  method scaleLike*(x: Greal[T], y: Gvalue): Gvalue =
    if y of Gscalar: return Gscalar(y)*x
    if y of Greal[T]: return Greal[T](y)*x
    raiseValueError("real gradient scaling requires a scalar or real value")
  method addLike*(prototype: Garray[T], x,y: Gvalue): Gvalue =
    if not (x of Garray[T]) or not (y of Garray[T]): raiseValueError("parameter addition precision differs")
    Garray[T](x)+Garray[T](y)
  method scaleLike*(x: Garray[T], y: Gvalue): Gvalue =
    if y of Gscalar: return Gscalar(y)*x
    raiseValueError("parameter gradient scaling requires a scalar")

realOps(float32)
realOps(float64)
