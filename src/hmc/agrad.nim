## Auto Grad implementation
import strformat

type
  ## a tape stores sequences of operations, each sequencs is a separate track
  AgTape* = ref object of RootObj
    track: int  ## default track to use
    trackops: seq[seq[AgOpBase]]  ## sequence of tracks of operations

  AgVarBase* = ref object of RootObj
    doGrad*: bool
    ctx*: AgTape

  AgVar*[T] = ref object of AgVarBase
    obj*: T
    grad*: T

  AgOpFnBase* = proc(op: AgOpBase) {.nimcall.}
  AgOpFn*[I,O] = proc(op: AgOp[I,O]) {.nimcall.}

  AgOpBase* = ref object of RootObj
    fwd: AgOpFnBase
    bck: AgOpFnBase
    #nIn: int
    #nOut: int
    #vars: seq[AgVarBase]

  AgOp*[I,O] = ref object of AgOpBase
    inputs: I
    outputs: O

# AgTape routines

proc newAgTape*: AgTape =
  result.new
  result.track = 0
  result.trackops.newSeq(1)
  result.trackops[0].newSeq(0)

proc `$`*(t: AgTape): string =
  result = &"AgTape[ track: {t.track} / {t.trackops.len}\n"
  for i in 0..<t.trackops.len:
    result &= &"  {i}: {t.trackops[i].len}\n"
  result &= "]"

proc numTracks*(t: AgTape): int =
  result = t.trackops.len

proc getTrack*(t: AgTape): int =
  result = t.track

proc setTrack*(t: var AgTape, i: int) =
  doAssert(i>=0)
  doAssert(i<t.trackops.len)
  t.track = i

proc addTrack*(t: var AgTape) =
  t.trackops.add @[]
  t.track = t.trackops.len - 1

template ops(t: AgTape): auto = t.trackops[t.track]

proc erase*(t: var AgTape) =
  t.ops.setLen 0

proc add*(t: var AgTape, o: AgOpBase) =
  t.ops.add o

proc run*(t: AgTape, verb = 0) =
  for i in 0..<t.ops.len:
    #if verb > 0: echo &"AgTape.run {i}: {t.ops[i]}"
    t.ops[i].fwd(t.ops[i])

proc grad*(t: AgTape) =
  let n = t.ops.len - 1
  # t.ops[n].vars[^1]
  for i in countdown(n, 0):
    t.ops[i].bck(t.ops[i])

# AgVar routines

proc newAgVar*[T](c: AgTape): AgVar[T] =
  result.new
  result.doGrad = true
  result.ctx = c

# AgOp routines

#proc newAgOp*[T,U](ip: T; op: U; ni,no: int; sv: seq[AgVarBase];
#                   fd,bk: untyped): AgOp[T,U] =
proc newAgOpImpl*[I,O](ip: I; op: O; fd,bk: AgOpFn[I,O]): AgOp[I,O] =
  result.new
  result.fwd = cast[AgOpFnBase](fd)
  result.bck = cast[AgOpFnBase](bk)
  #result.nIn = ni
  #result.nOut = no
  #result.vars = sv
  result.inputs = ip
  result.outputs = op
template newAgOp*[I,O](ip: I; op: O; fd,bk: untyped): AgOp[I,O] =
  newAgOpImpl(ip, op, fd[I,O], bk[I,O])

# predefined AgOp

template maybeObj(x: auto): auto =
  when x is AgVar: x.obj else: x

proc addfwd[I,O](op: AgOp[I,O]) {.nimcall.} =
  mixin add, zero
  add(op.outputs.obj, op.inputs[0].maybeObj, op.inputs[1].maybeObj)
  zero op.outputs.grad

proc addbck[I,O](op: AgOp[I,O]) {.nimcall.} =
  mixin `+=`
  when op.inputs[0] is AgVar:
    if op.inputs[0].doGrad:
      op.inputs[0].grad += op.outputs.grad
  when op.inputs[1] is AgVar:
    if op.inputs[1].doGrad:
      op.inputs[1].grad += op.outputs.grad

proc add(c: var AgTape, r: AgVar, x: auto, y: auto) =
  var op = newAgOp((x,y), r, addfwd, addbck)
  c.add op
template add(r: AgVar, x: auto, y: auto) =
  r.ctx.add(r, x, y)


proc mulfwd[I,O](op: AgOp[I,O]) {.nimcall.} =
  mixin mul, zero
  mul(op.outputs.obj, op.inputs[0].maybeObj, op.inputs[1].maybeObj)
  zero op.outputs.grad
proc mulbck[I,O](op: AgOp[I,O]) {.nimcall.} =
  mixin mul
  when op.inputs[0] is AgVar:
    if op.inputs[0].doGrad:
      muladd(op.inputs[0].grad, op.outputs.grad, op.inputs[1].maybeObj)
  when op.inputs[1] is AgVar:
    if op.inputs[1].doGrad:
      muladd(op.inputs[1].grad, op.inputs[0].maybeObj, op.outputs.grad)
proc mul(c: var AgTape, r: AgVar, x: auto, y: auto) =
  var op = newAgOp((x,y), r, mulfwd, mulbck)
  c.add op
template mul(r: AgVar, x: auto, y: auto) =
  r.ctx.mul(r, x, y)



when isMainModule:
  template assign(x: var SomeNumber; y: SomeNumber) =
    x = y
  template add(x: var SomeNumber; y,z: SomeNumber) =
    x = y + z
  template mul(x: var SomeNumber; y,z: SomeNumber) =
    x = y * z
  template muladd(x: var SomeNumber; y,z: SomeNumber) =
    x += y * z
  template zero(x: var Somenumber) =
    x = 0
  type
    FV = AgVar[float]

  template newFV(c: AgTape, x: float): untyped =
    var t = FV.new()
    #t.wantGrad = true
    t.doGrad = true
    t.ctx = c
    t.obj = x
    #t.grad = 0.0
    #zero t.grad
    t

  var c = newAgTape()
  var x = c.newFV(1.0)
  var y = c.newFV(2.0)
  var z = c.newFV(0.0)
  echo "x: ", x.obj
  echo "y: ", y.obj

  #echo z.obj
  add(z, x, y)
  c.run
  echo "z = x + y: ", z.obj
  x.grad = 0.0
  y.grad = 0.0
  z.grad = 2.0
  c.grad
  echo "z.grad: ", z.grad
  echo "x.grad: ", x.grad
  echo "y.grad: ", y.grad

  c.addTrack
  add(z, x, x)
  c.run
  x.grad = 0.0
  z.grad = 2.0
  echo "z = x + x: ", z.obj
  c.grad
  echo "z.grad: ", z.grad
  echo "x.grad: ", x.grad

  c.setTrack(0)
  c.run
  echo "z = x + y: ", z.obj
  x.grad = 0.0
  y.grad = 0.0
  z.grad = 2.0
  c.grad
  echo "z.grad: ", z.grad
  echo "x.grad: ", x.grad
  echo "y.grad: ", y.grad

  c.addTrack
  add(z, x, 1.0)
  c.run
  x.grad = 0.0
  z.grad = 5.0
  echo "z = x + 1.0: ", z.obj
  c.grad
  echo "z.grad: ", z.grad
  echo "x.grad: ", x.grad

  c.addTrack
  mul(z, x, y)
  c.run
  echo "z = x * y: ", z.obj
  x.grad = 0.0
  y.grad = 0.0
  z.grad = 3.0
  c.grad
  echo "z.grad: ", z.grad
  echo "x.grad: ", x.grad
  echo "y.grad: ", y.grad

  c.addTrack
  mul(z, x, x)
  c.run
  echo "z = x * x: ", z.obj
  x.grad = 0.0
  z.grad = 3.0
  c.grad
  echo "z.grad: ", z.grad
  echo "x.grad: ", x.grad

  c.addTrack
  mul(z, x, 3.0)
  c.run
  echo "z = x * 3.0: ", z.obj
  x.grad = 0.0
  z.grad = 4.0
  c.grad
  echo "z.grad: ", z.grad
  echo "x.grad: ", x.grad
