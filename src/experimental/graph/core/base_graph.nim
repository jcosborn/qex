import base_types, validate

template graphNode*(nodeValue,
                    inputsValue,
                    gfuncValue: untyped,
                    label = "graph node"): untyped =
  block:
    let inputs {.gensym.} = inputsValue
    let checkedInputs {.gensym.} = checkedInputValues(inputs, label)
    let gfunc {.gensym.} = gfuncValue
    if checkedInputs.len > 0 and gfunc == nil:
      raiseValueError(label & " with inputs requires a graph function")
    var node {.gensym.} = nodeValue
    discard node.requireGraphValue(label & " result")
    let inputGrt {.gensym.} = checkedInputs.sharedGraphRuntime
    # `nil` only means there were no inputs; values themselves have runtimes.
    if inputGrt != nil and node.runtime != inputGrt:
      raiseValueError(label & " mixes multiple graph runtimes")
    node.inputs = checkedInputs
    node.gfunc = gfunc
    node
