import core
import core/base

type
  Gmulti* {.final.} = ref object of Gvalue
    ## Multi-output carrier used by fused graph ops. Forward hooks read
    ## evaluated slot storage with `storedSlot`; backward builders use `[]` to
    ## build symbolic slot selections, allowing one fused backward graph to
    ## return all input-slot gradients together.
    slots: seq[Gvalue]

proc requireMultiArity(dstLen: int,
                       srcLen: int,
                       label: string) =
  if dstLen != srcLen:
    raiseValueError(
      label & " arity mismatch: " &
      $dstLen & " vs " & $srcLen)

proc requireMultiValue*(v: Gvalue,
                        label: string): Gmulti =
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
  if slotProtos.len == 0:
    raiseValueError(label & " requires at least one slot")
  var checkedSlots = newseq[Gvalue](slotProtos.len)
  for i in 0..<slotProtos.len:
    checkedSlots[i] =
      slotProtos[i].requireGraphValue(label & " slot prototype " & $i)
  let checkedInputs = checkedInputValues(inputs, label)
  if checkedInputs.len > 0 and gfuncValue == nil:
    raiseValueError(label & " with inputs requires a graph function")
  # Multi carriers take their runtime from output slot prototypes, not inputs.
  let slotGrt = checkedSlots.sharedGraphRuntime
  let inputGrt = checkedInputs.sharedGraphRuntime
  # `nil` only means there were no inputs; values themselves have runtimes.
  if inputGrt != nil and inputGrt != slotGrt:
    raiseValueError(label & " mixes multiple graph runtimes")
  var slotStorage = newseq[Gvalue](checkedSlots.len)
  for i in 0..<checkedSlots.len:
    slotStorage[i] = checkedSlots[i].newOneOf
  result = Gmulti().attachRuntime(slotGrt)
  result.slots = slotStorage
  result.inputs = checkedInputs
  result.gfunc = gfuncValue

proc `[]`*(x: Gmulti, i: int): Gvalue
proc multiValueNode(values: openArray[Gvalue],
                    label: string): Gmulti

proc requireSlotIndex(x: Gmulti,
                      index: int,
                      label: string,
                      node: Gvalue = nil): int =
  if index < 0 or index >= x.slots.len:
    raiseValueError(
      label & " index out of range: " & $index &
      " for size " & $x.slots.len &
      node.nodeContext)
  index

proc storedSlot*(x: Gmulti, k: int): Gvalue

proc requireMultiSelectBase(v: Gvalue,
                            label: string): Gmulti =
  v.requireInputCountExactly(1, label)
  v.requireNodeInput(0, label, "base").requireMultiValue(label & " base")

# Concrete slot storage stays separate from symbolic `x[k]` selection.
proc storedSlot*(x: Gmulti, k: int): Gvalue =
  ## Read last evaluated slot storage for forward hooks and tests.
  ## This does not build a graph node and is only current after evaluating the
  ## carrier or a consumer.
  x.slots[x.requireSlotIndex(k, "multi slot")]

proc requireMultiSlotCount(x: Gmulti,
                           expected: int,
                           label: string) =
  let actual = x.slots.len
  if actual != expected:
    raiseValueError(
      label & " expects " & $expected &
      " packed slots, got " & $actual)

proc tupleSlotCount[T: tuple](): int =
  var slots: T
  for _, _ in fieldPairs(slots):
    inc result

proc requireSlotAs[T: Gvalue](slot: Gvalue,
                              _: typedesc[T],
                              label: string): T =
  if slot == nil:
    raiseValueError(label & " expects " & $T & ", got nil")
  if not (slot of T):
    raiseValueError(label & " expects " & $T & ", got:\n" & slot.nodeRepr)
  T(slot)

proc storedSlots*[T: tuple](pack: Gmulti,
                            label: string): T =
  ## Read concrete slot storage into a typed tuple.
  pack.requireMultiSlotCount(tupleSlotCount[T](), label)
  var index = 0
  for name, field in fieldPairs(result):
    field = pack.storedSlot(index).requireSlotAs(
      typeof(field),
      label & "." & name)
    inc index

proc symbolicSlots*[T: tuple](pack: Gmulti,
                              label: string): T =
  ## Build symbolic slot selections into a typed tuple.
  pack.requireMultiSlotCount(tupleSlotCount[T](), label)
  var index = 0
  for name, field in fieldPairs(result):
    field = pack[index].requireSlotAs(
      typeof(field),
      label & "." & name)
    inc index

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
      values[j] = x.storedSlot(j).zeroLike
  multiValueNode(values, label)

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
  if not (value of Gmulti):
    return false
  prototype.copyCompatible(Gmulti(value))

method `$`*(x: Gmulti): string =
  $x.slots

proc multiValuesForward(v: Gvalue) =
  let z = v.requireMultiValue("values forward")
  ## Copy input values into slot storage.
  z.slots.copySlotValues(z.inputs, "values node storage/input")

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

proc multiValueNode(values: openArray[Gvalue],
                    label: string): Gmulti =
  ## Build a symbolic multi carrier from already-constructed graph values.
  newMultiOutputNode(values, values, multiValuesFunc, label)

proc multiValues*(label: string,
                  values: varargs[Gvalue]): Gmulti =
  ## Build a symbolic multi carrier from already-constructed graph values.
  multiValueNode(values, label)

const multiSelectSignatureNamespace = 0x6d756c7469530000'u64 # "multiS"

proc multiSelectSignatureKey(index: int): uint64 =
  multiSelectSignatureNamespace xor uint64(index)

proc newMultiSelectFunc(index: int): Gfunc =
  ## Slot index is static operator metadata. It must enter the signature key so
  ## selections of different slots cannot share cached gradient structure.
  let label = "multiSelect[" & $index & "]"

  proc multiSelectForward(v: Gvalue) =
    let base = v.requireMultiSelectBase(label)
    ## The slot index is fixed when the selection node is built.
    v.valCopy base.storedSlot(index)

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
  graphNode(x.storedSlot(k).newOneOf, @[Gvalue(x)], newMultiSelectFunc(k), "multiSelect")

proc `+`*(x: Gmulti, y: Gmulti): Gmulti =
  requireMultiArity(x.slots.len, y.slots.len, "multi add")
  var values = newseq[Gvalue](x.slots.len)
  for k in 0..<values.len:
    values[k] = x.storedSlot(k).addLike(x[k], y[k])
  multiValueNode(values, "multi add")

method addLike*(prototype: Gmulti, x: Gvalue, y: Gvalue): Gvalue =
  discard prototype
  let left = x.requireMultiValue("multi gradient add left")
  let right = y.requireMultiValue("multi gradient add right")
  left + right

method oneLike*(x: Gmulti): Gvalue =
  var values = newseq[Gvalue](x.slots.len)
  for k in 0..<values.len:
    values[k] = x.storedSlot(k).oneLike
  multiValueNode(values, "multi one")

method zeroLike*(x: Gmulti): Gvalue =
  var values = newseq[Gvalue](x.slots.len)
  for k in 0..<values.len:
    values[k] = x.storedSlot(k).zeroLike
  multiValueNode(values, "multi zero")
