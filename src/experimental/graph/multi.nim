import core

type
  Gmulti* {.final.} = ref object of Gvalue
    ## Fused multi-output carrier.
    slots: seq[Gvalue]
    shapeOnly: bool
  GmultiSelect = ref object of Gfunc
    index: int

proc requireMultiArity(dstLen: int,
                       srcLen: int,
                       label: string) =
  if dstLen != srcLen:
    raiseValueError(
      label & " arity mismatch: " &
      $dstLen & " vs " & $srcLen)

proc copySlotValues(dst: var seq[Gvalue],
                    src: openArray[Gvalue]) =
  # Both callers guarantee matching lengths: the forward path builds dst and src
  # from one seq, and cond-driven valCopy is shape-checked by copyCompatible.
  for i in 0..<dst.len:
    dst[i].valCopy(src[i])

# Construct multi-output carriers through this helper.
proc newMultiOutputNode*(slotProtos: openArray[Gvalue],
                         inputs: openArray[Gvalue],
                         gfuncValue: Gfunc,
                         label: string): Gmulti =
  if slotProtos.len == 0:
    raiseValueError(label & " requires at least one slot")
  let slotGrt = sharedGraphRuntime(slotProtos, label)
  var slotStorage = newseq[Gvalue](slotProtos.len)
  for i in 0..<slotProtos.len:
    slotStorage[i] = slotProtos[i].newOneOf
  # Multi carriers take their runtime from output slot prototypes, not inputs.
  result = Gmulti(runtime: slotGrt)
  result.slots = slotStorage
  result = graphNode(result, inputs, gfuncValue, label)

proc newMultiStructureNode*(slotProtos, inputs: openArray[Gvalue], gfuncValue: Gfunc, label: string): Gmulti =
  ## Construct a structural carrier whose slots are shape prototypes only.
  ## Its forward hook must not write them; consumers expose values through views.
  if slotProtos.len == 0:
    raiseValueError(label & " requires at least one slot")
  result = Gmulti(runtime: sharedGraphRuntime(slotProtos, label), shapeOnly: true)
  result.slots = @slotProtos
  result = graphNode(result, inputs, gfuncValue, label)

proc `[]`*(x: Gmulti, i: int): Gvalue
proc multiValues*(label: string, values: varargs[Gvalue]): Gmulti

# Concrete slot storage stays separate from symbolic `x[k]` selection.
proc storedSlot*(x: Gmulti, k: int): Gvalue =
  ## Last evaluated slot; does not build a selection node.
  if x.shapeOnly:
    raiseValueError("structural multi carrier has no stored slot values")
  x.slots[k]

proc mapSlots(x: Gmulti,
              label: string,
              f: proc(slot: Gvalue, k: int): Gvalue): Gmulti =
  ## Build a fresh multi carrier whose slot k is f(stored slot k, k).
  var values = newseq[Gvalue](x.slots.len)
  for k in 0..<values.len:
    values[k] = f(x.storedSlot(k), k)
  multiValues(label, values)

method newOneOf*(x: Gmulti): Gvalue =
  if x.shapeOnly:
    Gmulti(runtime: x.runtime, slots: x.slots, shapeOnly: true).assignStableNodeId
  else:
    newMultiOutputNode(x.slots, @[], nil, "multi prototype")

method valCopy*(z: Gmulti, x: Gvalue) =
  ## Copy slot values only; copyCompatible has checked their shapes.
  if z.shapeOnly or Gmulti(x).shapeOnly:
    raiseValueError("structural multi carrier cannot be copied")
  z.slots.copySlotValues(Gmulti(x).slots)

method copyCompatible*(prototype: Gmulti, value: Gvalue): bool =
  if not (value of Gmulti):
    return false
  let multiValue = Gmulti(value)
  if prototype.slots.len != multiValue.slots.len:
    return false
  for i in 0..<prototype.slots.len:
    if not prototype.slots[i].copyCompatible(multiValue.slots[i]):
      return false
  true

method `$`*(x: Gmulti): string =
  $x.slots

proc multiValuesForward(v: Gvalue) =
  ## Restore input aliases after generic graph cloning creates fresh slot prototypes.
  Gmulti(v).slots = v.inputs

proc multiValuesBackward(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let upstream = Gmulti(rootedUpstream(zb, z))
  upstream[i]

let multiValuesFunc = Gfunc(
  forward: multiValuesForward,
  backward: multiValuesBackward,
  name: "multiValues")

proc multiValues*(label: string, values: varargs[Gvalue]): Gmulti =
  ## Zero-copy bundle: slots alias inputs and the forward is a no-op.
  ## Use the bundle as a read source; destinations own their storage.
  if values.len == 0:
    raiseValueError(label & " requires at least one slot")
  result = Gmulti(runtime: sharedGraphRuntime(values, label))
  result.slots = @values
  result = graphNode(result, values, multiValuesFunc, label)

proc multiSelectForward(v: Gvalue) =
  let f = GmultiSelect(v.gfunc)
  let base = Gmulti(v.inputs[0])
  v.valCopy base.storedSlot(f.index)

proc multiSelectBackward(zb: Gvalue,
                         z: Gvalue,
                         i: int,
                         input: Gvalue): Gvalue =
  let f = GmultiSelect(z.gfunc)
  let base = Gmulti(z.inputs[0])
  let slotGrad = rootedUpstream(zb, z)
  base.mapSlots(f.name & " backward", proc(slot: Gvalue, j: int): Gvalue =
    if j == f.index: slotGrad else: slot.zeroLike)

proc newMultiSelectFunc(index: int): Gfunc =
  let label = "multiSelect[" & $index & "]"
  GmultiSelect(
    index: index,
    forward: multiSelectForward,
    backward: multiSelectBackward,
    name: label)

proc `[]`*(x: Gmulti, i: int): Gvalue =
  ## Return an aliased bundle input; otherwise build a slot selection.
  if x.shapeOnly:
    raiseValueError("structural multi carrier slots cannot be selected")
  if x.gfunc == multiValuesFunc:
    return x.inputs[i]
  let proto = x.storedSlot(i)
  graphNode(proto.newOneOf, @[Gvalue(x)], newMultiSelectFunc(i), "multiSelect")

method addLike*(prototype: Gmulti, x: Gvalue, y: Gvalue): Gvalue =
  let left = Gmulti(x)
  let right = Gmulti(y)
  requireMultiArity(left.slots.len, right.slots.len, "multi add")
  left.mapSlots("multi add", proc(slot: Gvalue, k: int): Gvalue =
    let
      a = left[k]
      b = right[k]
    if a.isStaticZeroLeaf:
      b
    elif b.isStaticZeroLeaf:
      a
    else:
      slot.addLike(a, b))

method oneLike*(x: Gmulti): Gvalue =
  x.mapSlots("multi one", proc(slot: Gvalue, k: int): Gvalue = slot.oneLike)

method zeroLike*(x: Gmulti): Gvalue =
  x.mapSlots("multi zero", proc(slot: Gvalue, k: int): Gvalue = slot.zeroLike)
