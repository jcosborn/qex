import base_types, validate

template graphNode*(nodeValue,
                    inputsValue,
                    gfuncValue: untyped,
                    label = "graph node"): untyped =
  block:
    var node {.gensym.} = nodeValue
    let inputs {.gensym.} = inputsValue
    node.inputs = checkedInputValues(inputs, label)
    let inputGrt {.gensym.} = node.inputs.sharedGraphRuntime
    # `nil` only means there were no inputs; values themselves have runtimes.
    if inputGrt != nil and node.runtime != inputGrt:
      raiseValueError(label & " mixes multiple graph runtimes")
    node.gfunc = gfuncValue
    node

template defineUnaryGraphOp*(gfuncName,
                             methodSym,
                             InputType,
                             inputSym,
                             resultNode,
                             forwardName,
                             backwardName: untyped,
                             gfuncLabel: static[string]) =
  let gfuncName = newGfunc(
    forward = forwardName,
    backward = backwardName,
    name = gfuncLabel)

  proc methodSym*(inputSym: InputType): typeof(resultNode) =
    graphNode(resultNode, @[Gvalue(inputSym)], gfuncName, gfuncLabel)

template defineBinaryGraphOp*(gfuncName,
                              methodSym,
                              LeftType,
                              RightType,
                              leftSym,
                              rightSym,
                              resultNode,
                              forwardName,
                              backwardName: untyped,
                              gfuncLabel: static[string]) =
  let gfuncName = newGfunc(
    forward = forwardName,
    backward = backwardName,
    name = gfuncLabel)

  proc methodSym*(leftSym: LeftType, rightSym: RightType): typeof(resultNode) =
    graphNode(resultNode, @[Gvalue(leftSym), Gvalue(rightSym)], gfuncName, gfuncLabel)
