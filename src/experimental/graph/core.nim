#[

- Graph traversals are not thread safe
- backward functions for scalar output may receive nil gradient

TODO

- function/lambda

]#

from strutils import join, toHex, strip

type
  Gfunc* {.acyclic.} = ref object
    ## Represent an functional operation: [input] -> output,
    forward: proc(z: Gvalue)
    backward: proc(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue  ## create new graph for backprop
    runCount: int
    name: string
  Gtag = enum
    gtVisited, gtGrad, gtFixedGrad
  Gtags = set[Gtag]
  Gvalue* {.acyclic.} = ref object of RootObj
    ## A Value knows its dependencies, which allows backpropagation.
    tag: Gtags
    inputs*: seq[Gvalue]
    gfunc*: Gfunc
    locals*: seq[Gvalue]  ## for sharing values between forward and among backward functions
    epoch: int

type
  GraphError* = object of Defect
  GraphValueError* = object of GraphError

template raiseError*(msg: string) =
  raise newException(GraphError, msg)

template raiseValueError*(msg: string) =
  raise newException(GraphValueError, msg)

template raiseErrorBaseMethod*(msg: string) =
  raiseError(
    "Base method invoked: " & msg &
    "\nMake sure to pass `--multimethods:on` and check there is a custom method for each derived type.")

var graphDebug* = false

proc newGfunc*(
    forward: proc(z: Gvalue) = nil,
    backward: proc(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue = nil,
    name: string): Gfunc =
  Gfunc(
    forward: forward,
    backward: backward,
    name: name)

proc runCount*(f: Gfunc): int = f.runCount

proc `$`*(x: Gfunc): string

method `$`*(x: Gvalue): string {.base.} =
  let f = x.gfunc
  result = "Gvalue(" & $x.epoch & " " & $x.tag & ")"
  if f != nil:
    result &= " " & $f

proc `$`*(x: Gfunc): string = x.name & "<" & $x.runCount & ">"

proc nodeRepr*(x: Gvalue): string =
  let f = x.gfunc
  result = $x & " (" & $x.epoch & " " & $x.tag & ")" & "@0X" & strip(toHex(cast[int](x)), trailing = false, chars = {'0'})
  if f != nil:
    result &= " " & $f & "@0X" & strip(toHex(cast[int](f)), trailing = false, chars = {'0'})

method newOneOf*(x: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("newOneOf(" & $x & ")")  ## Be sure to zero init fields
method valCopy*(z: Gvalue, x: Gvalue) {.base.} = raiseErrorBaseMethod("valCopy(" & $z & "," & $x & ")")

method isZero*(x: Gvalue): bool {.base.} = raiseErrorBaseMethod("isZero(" & $x & ")")
method update*(x: Gvalue, y: int) {.base.} = raiseErrorBaseMethod("update(" & $x & "," & $y & ")")
method update*(x: Gvalue, y: float) {.base.} = raiseErrorBaseMethod("update(" & $x & "," & $y & ")")

proc assignGvalue(z: Gvalue, x: Gvalue) =
  z.tag = x.tag
  z.inputs = x.inputs
  z.gfunc = x.gfunc
  z.epoch = x.epoch
  z.valCopy x

proc copyGvalue(x: Gvalue): Gvalue =
  result = newOneOf x
  result.assignGvalue x

let identPlaceholderGFunc = newGfunc(name = "identPlaceholder")
proc identPlaceholder(x: Gvalue): Gvalue =
  result = x.copyGvalue
  result.tag = {}
  result.inputs = @[x]
  result.gfunc = identPlaceholderGFunc
  result.epoch = 0

proc tagClearVisited(x: Gvalue) =
  ## only works after recursive proc used gtVisited for the graph traversal.
  if gtVisited in x.tag:
    x.tag.excl gtVisited
    for i in x.inputs:
      i.tagClearVisited

proc tagClear(x: Gvalue, t: Gtag) =
  proc c(v: Gvalue) =
    if gtVisited in v.tag:
      return
    v.tag.incl gtVisited
    v.tag.excl t
    for i in v.inputs:
      i.c
  x.c
  x.tagClearVisited

proc treeRepr*(v: Gvalue): string =
  var shared = newseq[Gvalue]()
  proc s(x: Gvalue) =
    if gtVisited in x.tag:
      if shared.find(x) < 0:
        shared.add x
    else:
      x.tag.incl gtVisited
      for i in x.inputs:
        i.s
  proc r(x: Gvalue): seq[string] =
    let si = shared.find x
    result = @[x.nodeRepr]
    if gtVisited in x.tag:
      result[0] &= " #" & $si
    else:
      if si >= 0:
        result[0] &= " #" & $si & "#"
      x.tag.incl gtVisited
      for i in x.inputs:
        for ir in i.r:
          result.add("  " & ir)
  v.s
  v.tagClearVisited
  result = v.r.join "\n"
  v.tagClearVisited

method `-`*(x: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("`-`(" & $x & ")")
method `+`*(x: Gvalue, y: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("`+`(" & $x & ", " & $y & ")")
method `*`*(x: Gvalue, y: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("`*`(" & $x & ", " & $y & ")")
method `-`*(x: Gvalue, y: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("`-`(" & $x & ", " & $y & ")")
method `/`*(x: Gvalue, y: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("`/`(" & $x & ", " & $y & ")")
method exp*(x: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("exp(" & $x & ")")

proc cond*(c: Gvalue, x: Gvalue, y: Gvalue): Gvalue

proc condb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  case i
  of 0:
    let r = z.inputs[0].newOneOf
    r.update 0
    return r
  of 1:
    if zb == nil:
      # the output must be a scalar, otherwise crash later
      let r1 = z.inputs[1].newOneOf
      let r0 = z.inputs[1].newOneOf
      r1.update 1
      r0.update 0
      return cond(z.inputs[0], r1, r0)
    else:
      return cond(z.inputs[0], zb, zb.newOneOf)
  of 2:
    if zb == nil:
      # the output must be a scalar, otherwise crash later
      let r0 = z.inputs[2].newOneOf
      let r1 = z.inputs[2].newOneOf
      r0.update 0
      r1.update 1
      return cond(z.inputs[0], r0, r1)
    else:
      return cond(z.inputs[0], zb.newOneOf, zb)
  else:
    raiseValueError("i must be 0 or 1, got: " & $i)

proc condf(v: Gvalue) =
  if v.inputs[0].isZero:
    v.valCopy v.inputs[2]
  else:
    v.valCopy v.inputs[1]

let gcond = newGfunc(forward = condf, backward = condb, name = "cond")

proc cond*(c: Gvalue, x: Gvalue, y: Gvalue): Gvalue =
  ## Assume the result is the same type as y, otherwise it'll throw exception later in forward valCopy.
  result = y.newOneOf
  result.inputs = @[c, x, y]
  result.gfunc = gcond

proc updated*(x: Gvalue) =
  var epoch {.global.} = 0
  inc epoch
  x.epoch = epoch

proc evaluated*(x: Gvalue) =
  ## signal up-to-date value, given inputs, useful for update value outside of eval, as in update to locals in forward
  var maxep = 0
  for i in x.inputs:
    if maxep < i.epoch:
      maxep = i.epoch
  x.epoch = maxep

proc eval*(v: Gvalue): Gvalue {.discardable.} =
  proc r(x: Gvalue) =
    if gtVisited in x.tag:
      return
    x.tag.incl gtVisited
    var maxep = 0
    if x.gfunc == gcond:
      x.inputs[0].r
      if maxep < x.inputs[0].epoch:
        maxep = x.inputs[0].epoch
      if x.inputs[0].isZero:
        x.inputs[2].r
        if maxep < x.inputs[2].epoch:
          maxep = x.inputs[2].epoch
      else:
        x.inputs[1].r
        if maxep < x.inputs[1].epoch:
          maxep = x.inputs[1].epoch
    else:
      for i in x.inputs:
        i.r
        if maxep < i.epoch:
          maxep = i.epoch
    if x.epoch < maxep:
      let f = x.gfunc
      if graphDebug:
        var s = "[graph/core] eval: " & x.nodeRepr
        for c in x.inputs:
          s &= "\n  " & c.nodeRepr
        echo s
      if f.forward != nil:
        x.epoch = maxep
        f.runCount.inc
        f.forward x
      else:
        raiseError("inputs.len: " & $x.inputs.len & ", but no forward function defined for:\n" & x.nodeRepr)
  v.r
  v.tagClearVisited
  v

type
  Grad = object
    input: Gvalue
    grad: Gvalue
  Grads = object
    output: Gvalue
    grads: seq[Grad]

var gradientList = newseq[Grads]()

proc dumpGradientList* =
  echo "# Gradient List:"
  for gs in gradientList:
    echo "## output: ",gs.output.nodeRepr
    for g in gs.grads:
      echo "### w.r.t.: ",g.input.nodeRepr
      echo g.grad.treeRepr

proc recordGrad(input: Gvalue, output: Gvalue, gradient: Gvalue) =
  for k in 0..<gradientList.len:
    if output == gradientList[k].output:
      for j in 0..<gradientList[k].grads.len:
        if input == gradientList[k].grads[j].input:
          var m = "Gradient exists for output:\n" & output.nodeRepr & "\nw.r.t. input:\n" & input.nodeRepr
          m &= "\nCurrent:\n" & gradientList[k].grads[j].grad.nodeRepr & "\nNew:\n" & gradient.nodeRepr
          raiseError m
      gradientList[k].grads.add Grad(input: input, grad: gradient)
      return
  gradientList.add Grads(output: output, grads: @[Grad(input: input, grad: gradient)])

proc findGrad*(input: Gvalue, output: Gvalue): Gvalue =
  ## Find the gradient of output with respect to input, may return nil
  var o = -1
  for k in 0..<gradientList.len:
    if output == gradientList[k].output:
      o = k
      break
  if o >= 0:
    for k in 0..<gradientList[o].grads.len:
      if input == gradientList[o].grads[k].input:
        return gradientList[o].grads[k].grad
  return nil

proc grad*(dep: Gvalue, x: Gvalue): Gvalue =
  proc t(v: Gvalue) =
    if gtVisited in v.tag:
      return
    v.tag.incl gtVisited
    var need = false
    for i in v.inputs:
      i.t
      need = need or gtGrad in i.tag
    if need:
      v.tag.incl gtGrad
  proc g(v: Gvalue) =
    let vgr = v.findGrad dep
    for i in 0..<v.inputs.len:
      let input = v.inputs[i]
      if gtGrad in input.tag:
        let f = v.gfunc
        if f.backward == nil:
          raiseError(v.nodeRepr & ":" & $i & ":" & input.nodeRepr & ": backward undefined")
        let gr = input.findGrad dep
        if gtVisited in input.tag:
          if gtFixedGrad notin input.tag:
            # We are in the process of building up the gradient, sum them up
            # Previous visit has its child grad reference gr, now copy and assign, no need to record
            if gr.gfunc == identPlaceholderGFunc:
              # we get the input of the placeholder so we don't leak identPlaceholder out of gradientList
              gr.assignGvalue(gr.inputs[0] + f.backward(vgr, v, i, dep))
            else:
              gr.assignGvalue(gr.copyGvalue + f.backward(vgr, v, i, dep))
          # else do nothing
        else:
          # first time for this child
          input.tag.incl gtVisited
          if gr == nil:
            # shared nodes will get revisited and assigned by a new node
            # use an identPlaceholder for all grad to avoid overwriting existing nodes returned from backward
            input.recordGrad(dep, f.backward(vgr, v, i, dep).identPlaceholder)
          else:
            input.tag.incl gtFixedGrad
          input.g
  x.tag.incl gtGrad
  dep.t
  dep.tagClearVisited
  dep.tag.incl gtVisited
  dep.g
  dep.tagClearVisited
  dep.tagClear gtGrad
  dep.tagClear gtFixedGrad
  # now remove the identPlaceholder
  for k in 0..<gradientList.len:
    if dep == gradientList[k].output:
      for j in 0..<gradientList[k].grads.len:
        if gradientList[k].grads[j].grad.gfunc == identPlaceholderGFunc:
          gradientList[k].grads[j].grad.assignGvalue(gradientList[k].grads[j].grad.inputs[0])
  x.findGrad dep
