import ../../[core, scalar]
import types

proc wrapper*(v: Gvalue): Gwrapper =
  if v == nil or not (v of Gwrapper):
    return nil
  Gwrapper(v)

proc isWrapperKind*(w: Gwrapper, kind: WrapperKind): bool {.inline.} =
  w != nil and w.kind == kind

proc localWrapper*(v: Gvalue): Gwrapper =
  let w = v.wrapper
  if not w.isWrapperKind(wkLocal):
    return nil
  w

proc callableWrapper*(v: Gvalue): Gwrapper =
  let w = v.wrapper
  if not w.isWrapperKind(wkCallable):
    return nil
  w

proc isLocalWrapper*(v: Gvalue): bool =
  v.localWrapper != nil

proc isCallableWrapper*(v: Gvalue): bool =
  v.callableWrapper != nil

proc callableKey*(fn: Glambda): NodeId =
  if fn == nil:
    return 0
  fn.stableNodeId

proc lambdaValue*(v: Gvalue): Glambda =
  if v == nil or not (v of Glambda):
    return nil
  Glambda(v)

proc resolvedClosure*(v: Gvalue): Glambda =
  result = v.lambdaValue
  if not result.isResolvedClosure:
    result = nil

proc resolvedLambdaValue*(v: Gvalue): Glambda =
  v.resolvedClosure

proc wrapperResultProto*(w: Gwrapper): Gvalue =
  if w == nil:
    return nil
  w.retProto

proc initWrapper(kind: WrapperKind,
                 grt: GraphRuntime,
                 retPrototype: Gvalue = nil): Gwrapper =
  Gwrapper(
    kind: kind,
    retProto: retPrototype).attachRuntime(grt)

proc callableWrapperNode*(retPrototype: Gvalue = nil): Gwrapper =
  if retPrototype == nil:
    raiseValueError("callable wrapper requires result prototype")
  initWrapper(
    wkCallable,
    retPrototype.runtime,
    retPrototype)

proc local*(grt: GraphRuntime): Gvalue =
  result = initWrapper(wkLocal, grt)
  result.updated

proc local*(retPrototype: Gvalue): Gvalue =
  result = initWrapper(
    wkLocal,
    retPrototype.runtime,
    retPrototype)
  result.updated

proc localValue*(prototype: Gvalue): Gvalue =
  prototype.newOneOf

proc localScalar*(grt: GraphRuntime): Gscalar =
  result = scalarNodeIn(grt)
  result.updated

proc localInt*(grt: GraphRuntime): Gint =
  result = intNodeIn(grt)
  result.updated

proc copyWrapperFields(dst: Gwrapper, src: Gwrapper) =
  dst.retProto = src.retProto
  dst.bound = src.bound

proc bindingResultProto*(v: Gvalue): Gvalue =
  if v == nil:
    return nil
  let fn = v.lambdaValue
  if fn.isResolvedClosure:
    return fn.body
  let w = v.wrapper
  if w == nil:
    return nil
  w.wrapperResultProto

proc requireCompatibleWrapperBinding(z: Gwrapper,
                                     x: Gvalue) =
  if x == nil:
    raiseValueError("callable wrapper binding cannot be nil")
  let expectedProto = z.wrapperResultProto
  if expectedProto == nil:
    return
  let actualProto = x.bindingResultProto
  if actualProto == nil:
    raiseValueError("callable wrapper expects callable binding, got:\n" & x.nodeRepr)
  if not expectedProto.copyCompatible(actualProto):
    raiseValueError(
      "callable wrapper binding has incompatible result prototype" &
      "\nexpected: " & expectedProto.nodeRepr &
      "\nactual: " & actualProto.nodeRepr)

method newOneOf*(x: Gwrapper): Gvalue =
  result = initWrapper(x.kind, x.runtime, x.retProto)

proc valCopy*(z: Gwrapper, x: Gwrapper) =
  if not z.copyCompatible(x):
    raiseValueError(
      "wrapper copy requires compatible result prototypes" &
      "\ndestination: " & z.nodeRepr &
      "\nsource: " & x.nodeRepr)
  if z.runtime != x.runtime:
    raiseValueError("wrapper copy mixes multiple graph runtimes")
  z.copyWrapperFields(x)
  if z.kind == wkLocal:
    z.updated

proc copyCompatible*(prototype: Gwrapper, value: Gwrapper): bool =
  if prototype.kind != value.kind:
    return false
  if prototype.retProto == nil or value.retProto == nil:
    return prototype.retProto == nil and value.retProto == nil
  prototype.retProto.copyCompatible(value.retProto)

method copyCompatible*(prototype: Gwrapper, value: Gvalue): bool =
  if value == nil or not (value of Gwrapper):
    return false
  prototype.copyCompatible(Gwrapper(value))

method valCopy*(z: Gwrapper, x: Gvalue) =
  z.requireCompatibleWrapperBinding(x)
  if z.runtime != x.runtime:
    raiseValueError("callable wrapper binding mixes multiple graph runtimes")
  z.bound = x
  if z.kind == wkLocal:
    z.updated

proc valCopy*(z: Gwrapper,
              x: float) =
  z.valCopy(toGvalue(z.runtime, x))

proc valCopy*(z: Gwrapper,
              x: int) =
  z.valCopy(toGvalue(z.runtime, x))

proc wrapperKindLabel(kind: WrapperKind): string =
  case kind
  of wkLocal:
    "local"
  of wkCallable:
    "callable"

proc callableShellRepr*(label: string, bound: Gvalue): string =
  if bound == nil:
    return label
  label & "(" & $bound & ")"

method `$`*(x: Gwrapper): string =
  callableShellRepr(wrapperKindLabel(x.kind), x.bound)

method newOneOf*(x: Glambda): Gvalue =
  result = Glambda(
    param: x.param).attachRuntime(x.runtime)

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
    if binding.value != nil:
      visit binding.value
