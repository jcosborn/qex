import ../../core
import ../../core/base

type
  LambdaBinding* = object
    param*: Gvalue
    value*: Gvalue
  WrapperKind* = enum
    wkLocal, wkCallable
  Gwrapper* {.final.} = ref object of Gvalue
    kind*: WrapperKind
    retProto*: Gvalue
    bound*: Gvalue
  Glambda* {.final.} = ref object of Gvalue
    param*: Gvalue
    body*: Gvalue
    env*: seq[LambdaBinding]
  CallableResolveMode* = enum
    crmShallow, crmReduced
  Bindings* = NodeTable[Gvalue]

proc isResolvedClosure*(fn: Glambda): bool {.inline.} =
  fn != nil and fn.body != nil

proc bindingResultProto*(v: Gvalue): Gvalue =
  if v of Glambda:
    let fn = Glambda(v)
    if fn.isResolvedClosure:
      return fn.body
    return nil
  if v of Gwrapper:
    return Gwrapper(v).retProto
  nil
