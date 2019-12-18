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
  for i in countdown(n, 1):
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
  op.outputs.obj = op.inputs[0].obj + op.inputs[1].obj

proc addbck[I,O](op: AgOp[I,O]) {.nimcall.} =
  if op.inputs[0].doGrad:
    op.inputs[0].grad += op.outputs.grad
  if op.inputs[1].doGrad:
    op.inputs[1].grad += op.outputs.grad

proc add[R,T,U](r: AgVar[R], x: AgVar[T], y: AgVar[U]) =
  assert(x.ctx == y.ctx)
  let ctx = x.ctx
  var op = newAgOp((x,y), r, 2, 1, @[x.AgVarBase,y,r], addfwd, addbck)
  ctx.ops.add op
  op.fwd(op)


type
  FV = AgVar[float]

template newFV(c: AgTape, x: float): untyped =
  var t = FV.new()
  t.wantGrad = true
  t.ctx = c
  t.obj = x
  t

var c = AgTape.new()

var x = c.newFV(1.0)
var y = c.newFV(2.0)
var z = c.newFV(4.0)

echo z.obj
add(z, x, y)
echo z.obj

