import ../../[core, scalar]
import types

proc callableWrapperNode*(retPrototype: Gvalue): Gwrapper =
  let proto = retPrototype.requireGraphValue("callable wrapper result prototype")
  Gwrapper(kind: wkCallable, retProto: proto).attachRuntime(proto.runtime)

proc local*(retPrototype: Gvalue): Gvalue =
  let proto = retPrototype.requireGraphValue("local result prototype")
  result = Gwrapper(kind: wkLocal, retProto: proto).attachRuntime(proto.runtime)
  result.updated

proc localScalar*(grt: GraphRuntime): Gscalar =
  result = scalarNodeIn(grt)
  result.updated

proc localInt*(grt: GraphRuntime): Gint =
  result = intNodeIn(grt)
  result.updated

proc wrapperProtosCompatible(prototype: Gwrapper, value: Gwrapper): bool =
  if prototype.kind != value.kind:
    return false
  if prototype.retProto == nil or value.retProto == nil:
    return false
  prototype.retProto.copyCompatible(value.retProto)

proc requireCompatibleWrapperBinding(z: Gwrapper,
                                     x: Gvalue,
                                     label: string) =
  let expectedProto = z.retProto
  if expectedProto == nil:
    raiseValueError(label & " requires a result prototype:\n" & z.nodeRepr)
  let actualProto = x.bindingResultProto
  if actualProto == nil:
    raiseValueError(label & " expects callable binding, got:\n" & x.nodeRepr)
  if not expectedProto.copyCompatible(actualProto):
    raiseValueError(
      label & " binding has incompatible result prototype" &
      "\nexpected: " & expectedProto.nodeRepr &
      "\nactual: " & actualProto.nodeRepr)

method newOneOf*(x: Gwrapper): Gvalue =
  Gwrapper(kind: x.kind, retProto: x.retProto).attachRuntime(x.runtime)

proc valCopy*(z: Gwrapper, x: Gwrapper) =
  if not z.wrapperProtosCompatible(x):
    raiseValueError(
      "wrapper copy requires compatible result prototypes" &
      "\ndestination: " & z.nodeRepr &
      "\nsource: " & x.nodeRepr)
  if z.runtime != x.runtime:
    raiseValueError("wrapper copy mixes multiple graph runtimes")
  z.retProto = x.retProto
  z.bound = x.bound
  z.updated

method valCopy*(z: Gwrapper, x: Gvalue) =
  z.requireCompatibleWrapperBinding(x, "callable wrapper")
  if z.runtime != x.runtime:
    raiseValueError("callable wrapper binding mixes multiple graph runtimes")
  z.bound = x
  z.updated

method copyCompatible*(prototype: Gwrapper, value: Gvalue): bool =
  if not (value of Gwrapper):
    return false
  prototype.wrapperProtosCompatible(Gwrapper(value))

method `$`*(x: Gwrapper): string =
  let label =
    case x.kind
    of wkLocal:
      "local"
    of wkCallable:
      "callable"
  if x.bound == nil:
    return label
  label & "(" & $x.bound & ")"

method newOneOf*(x: Glambda): Gvalue =
  Glambda(param: x.param).attachRuntime(x.runtime)

method `$`*(x: Glambda): string =
  result = "lambda(" & $x.param & " -> " & $x.body & ")"
  if x.env.len > 0:
    result &= "[env:" & $x.env.len & "]"

method walkHiddenDeps*(v: Glambda,
                       mode: InputWalkMode,
                       visit: GnodeVisit) =
  if mode notin {iwmGradSignature, iwmDepend} or not v.isResolvedClosure:
    return
  for binding in v.env:
    visit binding.value
