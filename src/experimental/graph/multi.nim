import core

type
  Gmulti* {.final.} = ref object of Gvalue
    ## Multi-output carrier used by fused graph ops. Forward hooks read
    ## evaluated slot storage with `storedSlot`; backward builders use `[]` to
    ## build symbolic slot selections, allowing one fused backward graph to
    ## return all input-slot gradients together.
    slots: seq[Gvalue]
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

proc `[]`*(x: Gmulti, i: int): Gvalue
proc multiValues*(label: string, values: varargs[Gvalue]): Gmulti

# Concrete slot storage stays separate from symbolic `x[k]` selection.
proc storedSlot*(x: Gmulti, k: int): Gvalue =
  ## Read last evaluated slot storage for forward hooks and tests.
  ## This does not build a graph node and is only current after evaluating the
  ## carrier or a consumer.
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
  newMultiOutputNode(x.slots, @[], nil, "multi prototype")

method valCopy*(z: Gmulti, x: Gvalue) =
  ## Copy concrete slot storage only; graph structure is not cloned. The only
  ## runtime path here is `cond` selecting between carriers, which `copyCompatible`
  ## already required to have matching slot shape.
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
  let z = Gmulti(v)
  ## Copy input values into slot storage.
  z.slots.copySlotValues(z.inputs)

proc multiValuesBackward(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
  let upstream = Gmulti(rootedUpstream(zb, z))
  upstream[i]

let multiValuesFunc = Gfunc(
  forward: multiValuesForward,
  backward: multiValuesBackward,
  name: "multiValues")

proc multiValues*(label: string,
                  values: varargs[Gvalue]): Gmulti =
  ## Build a symbolic multi carrier from already-constructed graph values.
  newMultiOutputNode(values, values, multiValuesFunc, label)

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
  ## Build a symbolic selection for a statically chosen slot.
  let proto = x.storedSlot(i)
  graphNode(proto.newOneOf, @[Gvalue(x)], newMultiSelectFunc(i), "multiSelect")

method addLike*(prototype: Gmulti, x: Gvalue, y: Gvalue): Gvalue =
  let left = Gmulti(x)
  let right = Gmulti(y)
  requireMultiArity(left.slots.len, right.slots.len, "multi add")
  left.mapSlots("multi add", proc(slot: Gvalue, k: int): Gvalue =
    slot.addLike(left[k], right[k]))

method oneLike*(x: Gmulti): Gvalue =
  x.mapSlots("multi one", proc(slot: Gvalue, k: int): Gvalue = slot.oneLike)

method zeroLike*(x: Gmulti): Gvalue =
  x.mapSlots("multi zero", proc(slot: Gvalue, k: int): Gvalue = slot.zeroLike)
