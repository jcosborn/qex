import core

type
  Gmulti* {.final.} = ref object of Gvalue
    ## Fused multi-output carrier.
    slots: seq[Gvalue]
    shapeOnly: bool
  GmultiSelect = ref object of Gfunc
    index: int

proc isAliasBundle(x: Gmulti): bool

method isStructuralValue*(x: Gmulti): bool = x.shapeOnly

method hasStorage*(x: Gmulti): bool =
  if x.shapeOnly:
    return true
  for slot in x.slots:
    if not slot.hasStorage:
      return false
  true

method ensureStorage*(x: Gmulti) =
  if x.shapeOnly or x.isAliasBundle:
    return
  for slot in x.slots:
    if not slot.hasStorage:
      slot.ensureStorage

method releaseStorage*(x: Gmulti) =
  if not x.shapeOnly:
    for slot in x.slots:
      if not slot.isStructuralValue:
        slot.releaseStorage

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
    if src[i].isStructuralValue:
      dst[i] = src[i]
    else:
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
    slotStorage[i] = slotProtos[i].valueLike
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

proc len*(x: Gmulti): int = x.slots.len

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
    var slots = newSeq[Gvalue](x.slots.len)
    for i, slot in x.slots:
      slots[i] = if slot.isStructuralValue: slot else: slot.newOneOf
    graphNode(Gmulti(runtime: x.runtime, slots: slots),
      newSeq[Gvalue](), nil, "multi prototype")

method valueLike*(x: Gmulti): Gvalue =
  if x.shapeOnly:
    return x.newOneOf
  var slots = newSeq[Gvalue](x.slots.len)
  for i, slot in x.slots:
    slots[i] = if slot.isStructuralValue: slot else: slot.valueLike
  Gmulti(runtime: x.runtime, slots: slots).assignStableNodeId

method bufferBytes*(x: Gmulti): int =
  if not x.shapeOnly:
    for slot in x.slots:
      result += slot.bufferBytes

method bufferProto*(x: Gmulti): Gvalue =
  if x.shapeOnly:
    return nil
  var slots = newSeq[Gvalue](x.slots.len)
  var hasBuffer = false
  for i, slot in x.slots:
    let buf = slot.bufferProto
    if buf != nil:
      slots[i] = buf
      hasBuffer = true
    elif slot.bufferBytes > 0:
      return nil
    else:
      slots[i] = if slot.isStructuralValue: slot else: slot.newOneOf
  if hasBuffer:
    result = Gmulti(runtime: x.runtime, slots: slots).assignStableNodeId

method bufferCompatible*(x: Gmulti, y: Gvalue): bool =
  if not (y of Gmulti):
    return false
  let src = Gmulti(y)
  if x.shapeOnly or src.shapeOnly or x.slots.len != src.slots.len:
    return false
  for i, slot in x.slots:
    let other = src.slots[i]
    if slot.bufferBytes == 0:
      if other.bufferBytes != 0:
        return false
    elif not slot.bufferCompatible(other) or not other.bufferCompatible(slot):
      return false
  true

method bindBuffer*(x: Gmulti, y: Gvalue) =
  let src = Gmulti(y)
  for i, slot in x.slots:
    if slot.bufferBytes > 0:
      slot.bindBuffer(src.slots[i])

method clearBuffer*(x: Gmulti) =
  if not x.shapeOnly:
    for slot in x.slots:
      if slot.bufferBytes > 0:
        slot.clearBuffer

method valCopy*(z: Gmulti, x: Gvalue) =
  ## Copy slot values only; copyCompatible has checked their shapes.
  if z.shapeOnly or Gmulti(x).shapeOnly:
    raiseValueError("structural multi carrier cannot be copied")
  z.slots.copySlotValues(Gmulti(x).slots)

method valAlias*(z: Gmulti, x: Gvalue) =
  let src = Gmulti(x)
  if z.shapeOnly and src.shapeOnly:
    return
  if z.shapeOnly or src.shapeOnly:
    raiseValueError("structural multi carrier cannot hold slot values")
  for i in 0..<z.slots.len:
    if src.slots[i].isStructuralValue:
      z.slots[i] = src.slots[i]
    else:
      z.slots[i].valAlias(src.slots[i])

method slotForward*(z: Gmulti, x: Gvalue) =
  if not z.shapeOnly:
    if z.runtime.evalFrame != nil:
      z.valAlias(x)
    else:
      z.valCopy(x)

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
  ## Structural slots retain their declared input; numerical wrappers own descriptors.
  let z = Gmulti(v)
  for i in 0..<z.slots.len:
    if v.inputs[i].isStructuralValue:
      z.slots[i] = v.inputs[i]
    else:
      z.slots[i].valAlias(v.inputs[i])

proc multiValuesBackward(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let upstream = Gmulti(rootedUpstream(zb, z))
  upstream[i]

let multiValuesFunc = Gfunc(
  forward: multiValuesForward,
  backward: multiValuesBackward,
  bufferMode: bmAlias,
  name: "multiValues")

proc isAliasBundle(x: Gmulti): bool = x.gfunc == multiValuesFunc

proc multiValues*(label: string, values: varargs[Gvalue]): Gmulti =
  ## Numerical slots own wrappers; structural slots retain their input values.
  ## Use the bundle as a read source; destinations own their storage.
  if values.len == 0:
    raiseValueError(label & " requires at least one slot")
  result = Gmulti(runtime: sharedGraphRuntime(values, label),
    slots: newSeq[Gvalue](values.len))
  for i, value in values:
    result.slots[i] = if value.isStructuralValue: value else: value.newOneOf
  result = graphNode(result, values, multiValuesFunc, label)

proc multiSelectForward(v: Gvalue) =
  let f = GmultiSelect(v.gfunc)
  let base = Gmulti(v.inputs[0])
  if v.runtime.evalFrame != nil:
    v.valAlias base.storedSlot(f.index)
  else:
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
    bufferMode: bmAlias,
    name: label)

proc `[]`*(x: Gmulti, i: int): Gvalue =
  ## Return an aliased bundle input; otherwise build a slot selection.
  if x.shapeOnly:
    raiseValueError("structural multi carrier slots cannot be selected")
  if x.gfunc == multiValuesFunc:
    return x.inputs[i]
  let proto = x.storedSlot(i)
  graphNode(proto.valueLike, @[Gvalue(x)], newMultiSelectFunc(i), "multiSelect")

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
