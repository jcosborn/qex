import ../../core

type
  LambdaBinding* = object
    param*: Gvalue
    value*: Gvalue
  WrapperKind* = enum
    wkLocal, wkCallable
  Gwrapper* {.final.} = ref object of Gvalue
    ## `wkLocal` is a directly rebound placeholder.
    ## `wkCallable` caches a function-valued binding while its producer is fresh.
    kind*: WrapperKind
    retProto*: Gvalue
    bound*: Gvalue
  Glambda* {.final.} = ref object of Gvalue
    param*: Gvalue
    body*: Gvalue
    env*: seq[LambdaBinding]
  CallableResolveMode* = enum
    crmShallow, crmReduced
  CallableInspect* = object
    directValue*: Gvalue
    directFn*: Glambda
    reducedValue*: Gvalue
    reducedFn*: Glambda
    hasReduced*: bool
  Bindings* = NodeTable[Gvalue]

proc isResolvedClosure*(fn: Glambda): bool {.inline.} =
  fn != nil and fn.body != nil

proc initBindings*(): Bindings =
  initNodeTable[Gvalue]()

proc bindNode*(binding: var Bindings, key, value: Gvalue) {.inline.} =
  binding.putNode(key, value)

proc deleteBinding*(binding: var Bindings, key: Gvalue) {.inline.} =
  if key != nil and binding.hasNode(key):
    binding.delNode(key)

proc isDeferredApplyValue*(v: Gvalue): bool =
  v.hasReduceValue

proc reduceDeferredApplyValue*(v: Gvalue): Gvalue =
  v.runReduceValue
