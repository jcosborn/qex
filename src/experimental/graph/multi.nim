import core

type
  Gmulti* {.final.} = ref object of Gvalue
    slots: seq[Gvalue]  ## Forward storage for each output slot.
                        ## Read storage with `slotValue`; build graph expressions with `[]`.

proc requireConcreteSlotIndex(index: int,
                              size: int,
                              label: string,
                              node: Gvalue = nil): int =
  if index < 0 or index >= size:
    raiseValueError(
      label & " index out of range: " & $index &
      " for size " & $size &
      node.nodeContext)
  index

proc requireMultiSlotProto(slotProto: Gvalue, label: string) =
  if slotProto == nil:
    raiseValueError(label & " slot prototype cannot be nil")

proc allocSlots(slotProtos: openArray[Gvalue],
                label: string): seq[Gvalue] =
  if slotProtos.len == 0:
    raiseValueError(label & " requires at least one slot")
  result = newseq[Gvalue](slotProtos.len)
  for i in 0..<slotProtos.len:
    let slotProto = slotProtos[i]
    slotProto.requireMultiSlotProto(label)
    result[i] = slotProto.newOneOf

proc requireMultiArity(dstLen: int,
                       srcLen: int,
                       label: string) =
  if dstLen != srcLen:
    raiseValueError(
      label & " arity mismatch: " &
      $dstLen & " vs " & $srcLen)

proc requireMultiValue*(v: Gvalue,
                        label: string): Gmulti =
  if v == nil:
    raiseValueError(label & " requires non-nil multi value")
  if not (v of Gmulti):
    raiseValueError(label & " expects multi value")
  Gmulti(v)

proc requireMultiUpstream*(zb: Gvalue,
                           label: string): Gmulti =
  if zb == nil:
    raiseValueError(label & " requires non-nil upstream gradient")
  zb.requireMultiValue(label & " upstream gradient")

proc copySlotValues(dst: var seq[Gvalue],
                    src: openArray[Gvalue],
                    label: string) =
  requireMultiArity(dst.len, src.len, label)
  for i in 0..<dst.len:
    dst[i].valCopy(src[i])

# Construct multi-output carriers through this helper.
proc newMultiOutputNode*(slotProtos: openArray[Gvalue],
                         inputs: openArray[Gvalue],
                         gfuncValue: Gfunc,
                         label: string): Gmulti =
  let checkedInputs = checkedInputValues(inputs, label)
  let slotStorage = allocSlots(slotProtos, label)
  # Multi carriers take their runtime from output slot prototypes, not inputs.
  let slotGrt = slotProtos.sharedGraphRuntime
  result = Gmulti().attachRuntime(slotGrt)
  result.defineGraphNode(checkedInputs, gfuncValue)
  result.slots = slotStorage
  let inputGrt = checkedInputs.sharedGraphRuntime
  # `nil` only means there were no inputs; values themselves have runtimes.
  if inputGrt != nil and inputGrt != result.runtime:
    raiseValueError(label & " mixes multiple graph runtimes")

proc multiCarrierFromExprs(values: openArray[Gvalue],
                           label = "multi values"): Gmulti
proc `[]`*(x: Gmulti, i: int): Gvalue

proc refreshSlotStorage(x: Gmulti, label: string) =
  x.slots.copySlotValues(x.inputs, label)

proc requireSlotIndex(x: Gmulti,
                      index: int,
                      label: string,
                      node: Gvalue = nil): int =
  result = requireConcreteSlotIndex(index, x.slots.len, label, node)

proc slotValue*(x: Gmulti, k: int): Gvalue

proc requireMultiSelectBase(v: Gvalue,
                            label: string): Gmulti =
  v.requireInputCountExactly(1, label)
  v.requireNodeInput(0, label, "base").requireMultiValue(label & " base")

# Concrete slot storage stays separate from symbolic `x[k]` selection.
proc slotValue*(x: Gmulti, k: int): Gvalue =
  ## Read last evaluated slot storage for forward hooks and tests.
  ## This does not build a graph node and is only current after evaluating the
  ## carrier or a consumer.
  x.slots[x.requireSlotIndex(k, "multi slot")]

proc scatterSlotGrad(x: Gmulti,
                     index: int,
                     slotGrad: Gvalue,
                     label: string,
                     node: Gvalue): Gvalue =
  let slotIndex = x.requireSlotIndex(index, label, node)
  var values = newseq[Gvalue](x.slots.len)
  for j in 0..<values.len:
    if j == slotIndex:
      values[j] = slotGrad
    else:
      values[j] = x.slotValue(j).zeroLike
  multiCarrierFromExprs(values, label)

method newOneOf*(x: Gmulti): Gvalue =
  newMultiOutputNode(x.slots, @[], nil, "multi prototype")

proc valCopy*(z: Gmulti, x: Gmulti) =
  ## Copy concrete slot storage only; graph structure is not cloned.
  if z.runtime != x.runtime:
    raiseValueError("multi copy mixes multiple graph runtimes")
  z.slots.copySlotValues(x.slots, "multi copy")

method valCopy*(z: Gmulti, x: Gvalue) =
  z.valCopy(x.requireMultiValue("multi copy"))

proc copyCompatible*(prototype: Gmulti, value: Gmulti): bool =
  if prototype.slots.len != value.slots.len:
    return false
  for i in 0..<prototype.slots.len:
    if not prototype.slots[i].copyCompatible(value.slots[i]):
      return false
  true

method copyCompatible*(prototype: Gmulti, value: Gvalue): bool =
  if value == nil or not (value of Gmulti):
    return false
  prototype.copyCompatible(Gmulti(value))

method `$`*(x: Gmulti): string =
  $x.slots

proc multiValuesForward(v: Gvalue) =
  let z = v.requireMultiValue("values forward")
  ## Copy input values into slot storage.
  z.refreshSlotStorage("values node storage/input")

proc multiValuesBackward(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  discard dep
  if i < 0 or i >= z.inputs.len:
    raiseValueError("values backward input index out of range: " & $i)
  let upstream = requireMultiUpstream(zb, "values backward")
  upstream[i]

let multiValuesFunc = newGfunc(
  forward = multiValuesForward,
  backward = multiValuesBackward,
  name = "multiValues")

proc multiCarrierFromExprs(values: openArray[Gvalue],
                           label = "multi values"): Gmulti =
  newMultiOutputNode(values, values, multiValuesFunc, label)

const multiSelectSignatureNamespace = 0x6d756c7469530000'u64 # "multiS"

proc multiSelectSignatureKey(index: int): uint64 =
  multiSelectSignatureNamespace xor uint64(index)

proc newMultiSelectFunc(index: int): Gfunc =
  let label = "multiSelect[" & $index & "]"

  proc multiSelectForward(v: Gvalue) =
    let base = v.requireMultiSelectBase(label)
    ## The slot index is fixed when the selection node is built.
    v.valCopy base.slotValue(index)

  proc multiSelectBackward(zb: Gvalue,
                           z: Gvalue,
                           i: int,
                           dep: Gvalue): Gvalue =
    discard dep
    if i != 0:
      raiseValueError(label & " input index must be 0, got: " & $i)
    let base = z.requireMultiSelectBase(label)
    let slotGrad =
      if zb == nil:
        z.oneLike
      else:
        zb
    base.scatterSlotGrad(index, slotGrad, label & " backward", z)

  newGfunc(
    forward = multiSelectForward,
    backward = multiSelectBackward,
    signatureKey = multiSelectSignatureKey(index),
    name = label)

proc `[]`*(x: Gmulti, i: int): Gvalue =
  ## Build a symbolic selection for a statically chosen slot.
  let k = x.requireSlotIndex(i, "[]", x)
  graphNode(x.slotValue(k).newOneOf, @[Gvalue(x)], newMultiSelectFunc(k), "multiSelect")

proc `+`*(x: Gmulti, y: Gmulti): Gmulti =
  requireMultiArity(x.slots.len, y.slots.len, "multi add")
  var values = newseq[Gvalue](x.slots.len)
  for k in 0..<values.len:
    values[k] = x.slotValue(k).addLike(x[k], y[k])
  multiCarrierFromExprs(values, "multi add")

method addLike*(prototype: Gmulti, x: Gvalue, y: Gvalue): Gvalue =
  let left = x.requireMultiValue("multi gradient add left")
  let right = y.requireMultiValue("multi gradient add right")
  requireMultiArity(left.slots.len, right.slots.len, "multi gradient add")
  requireMultiArity(prototype.slots.len, left.slots.len, "multi gradient add prototype")
  var values = newseq[Gvalue](prototype.slots.len)
  for k in 0..<values.len:
    values[k] = prototype.slotValue(k).addLike(left[k], right[k])
  multiCarrierFromExprs(values, "multi gradient add")

method oneLike*(x: Gmulti): Gvalue =
  var values = newseq[Gvalue](x.slots.len)
  for k in 0..<values.len:
    values[k] = x.slotValue(k).oneLike
  multiCarrierFromExprs(values, "multi one")

method zeroLike*(x: Gmulti): Gvalue =
  var values = newseq[Gvalue](x.slots.len)
  for k in 0..<values.len:
    values[k] = x.slotValue(k).zeroLike
  multiCarrierFromExprs(values, "multi zero")
