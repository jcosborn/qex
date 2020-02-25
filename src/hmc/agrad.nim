type
  AgTape = ref object of RootObj
    ops: seq[AgOpBase]

  AgVarBase = ref object of RootObj
    wantGrad: bool
    doGrad: bool
    ctx: AgTape

  AgVar[T] = ref object of AgVarBase
    obj: T
    grad: T

  AgOpBaseArg = proc(op: AgOpBase) {.nimcall.}

  AgOpBase = ref object of RootObj
    fwd: AgOpBaseArg
    bck: AgOpBaseArg
    nIn: int
    nOut: int
    vars: seq[AgVarBase]

  AgOp[I,O] = ref object of AgOpBase
    inputs: I
    outputs: O

proc run*(t: AgTape) =
  for i in 0..<t.ops.len:
    t.ops[i].fwd(t.ops[i])

proc grad*(t: AgTape) =
  let n = t.ops.len - 1
  # t.ops[n].vars[^1]
  for i in countdown(n, 0):
    t.ops[i].bck(t.ops[i])

template newAgOp[T,U](ip: T; op: U; ni,no: int; sv: seq[AgVarBase];
                      fd,bk: untyped): untyped =
  var o = AgOp[T,U].new()
  o.fwd = cast[AgOpBaseArg](fd[T,U])
  o.bck = cast[AgOpBaseArg](bk[T,U])
  o.nIn = ni
  o.nOut = no
  o.vars = sv
  o.inputs = ip
  o.outputs = op
  o

proc addfwd[I,O](op: AgOp[I,O]) {.nimcall.} =
  mixin add
  add(op.outputs.obj, op.inputs[0].obj, op.inputs[1].obj)

proc addbck[I,O](op: AgOp[I,O]) {.nimcall.} =
  mixin assign
  if op.inputs[0].doGrad:
    assign(op.inputs[0].grad, op.outputs.grad)
  if op.inputs[1].doGrad:
    assign(op.inputs[1].grad, op.outputs.grad)

proc add[R,T,U](c: AgTape, r: AgVar[R], x: AgVar[T], y: AgVar[U]) =
  var op = newAgOp((x,y), r, 2, 1, @[x.AgVarBase,y,r], addfwd, addbck)
  c.ops.add op
  op.fwd(op)
template add[R,T,U](r: AgVar[R], x: AgVar[T], y: AgVar[U]) =
  r.ctx.add(r, x, y)


proc mulfwd[I,O](op: AgOp[I,O]) {.nimcall.} =
  mixin mul
  mul(op.outputs.obj, op.inputs[0].obj, op.inputs[1].obj)

proc mulbck[I,O](op: AgOp[I,O]) {.nimcall.} =
  mixin mul
  if op.inputs[0].doGrad:  # mulna
    mul(op.inputs[0].grad, op.outputs.grad, op.inputs[1].obj)
  if op.inputs[1].doGrad:  # mulan
    mul(op.inputs[1].grad, op.inputs[0].obj, op.outputs.grad)

proc mul[R,T,U](c: AgTape, r: AgVar[R], x: AgVar[T], y: AgVar[U]) =
  var op = newAgOp((x,y), r, 2, 1, @[x.AgVarBase,y,r], mulfwd, mulbck)
  c.ops.add op
  op.fwd(op)
template mul[R,T,U](r: AgVar[R], x: AgVar[T], y: AgVar[U]) =
  r.ctx.mul(r, x, y)



when isMainModule:
  template assign(x: var SomeNumber; y: SomeNumber) =
    x = y
  template add(x: var SomeNumber; y,z: SomeNumber) =
    x = y + z
  template mul(x: var SomeNumber; y,z: SomeNumber) =
    x = y * z
  type
    FV = AgVar[float]

  template newFV(c: AgTape, x: float): untyped =
    var t = FV.new()
    t.wantGrad = true
    t.doGrad = true
    t.ctx = c
    t.obj = x
    t.grad = 0.0
    t

  var c = AgTape.new()
  var c2 = AgTape.new()

  var x = c.newFV(1.0)
  var y = c.newFV(2.0)
  var z = c.newFV(4.0)

  echo z.obj
  add(z, x, y)
  echo z.obj

  z.grad = 1.0
  c.grad()
  echo "x.grad: ", x.grad
  echo "y.grad: ", y.grad

  c2.mul(z, x, y)
  z.grad = 2.0
  c2.grad()
  echo z.obj
  echo "x.grad: ", x.grad
  echo "y.grad: ", y.grad
