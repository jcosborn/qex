import base, traverse

proc appendRawInputSignatureTokens(node: Gvalue,
                                   tokens: var GradSignature) =
  for inputIndex in 0..<node.inputs.len:
    let input = node.inputs[inputIndex]
    if input == nil:
      raiseError(
        "grad signature node has nil input at index " & $inputIndex &
        ":\n" & node.nodeRepr)
    tokens.add GradSigToken(kind: gstInput, key: input.stableNodeId)

proc appendNodeSignatureTokens(node: Gvalue,
                               tokens: var GradSignature) =
  tokens.add GradSigToken(kind: gstNode, key: node.stableNodeId)
  let graphFunc = node.gfunc
  if graphFunc != nil:
    tokens.add GradSigToken(kind: gstFunc, key: signatureKey(graphFunc))
  # Raw input identities remain part of the signature even when the
  # grad-signature walk hides internal structure. This keeps root identity
  # changes cache-visible without forcing recursive traversal into lazy branches
  # or deferred apply-partial targets.
  node.appendRawInputSignatureTokens(tokens)
  node.appendGfuncSignature(tokens)
  node.appendSignatureTokens(tokens)

proc buildGradSignature*(dep: Gvalue): GradSignature =
  var seen = initNodeSet()
  var stack: seq[tuple[node: Gvalue, expanded: bool]] = @[(dep, false)]
  var signatureOrder: seq[Gvalue] = @[]

  while stack.len > 0:
    let frame = stack[^1]
    stack.setLen(stack.len - 1)

    if frame.node == nil:
      raiseError("grad signature traversal encountered nil node")

    if frame.expanded:
      signatureOrder.add frame.node
      continue

    if not seen.markSeenNode(frame.node):
      continue

    stack.add((frame.node, true))
    let signatureDeps = frame.node.collectNodeInputs(iwmGradSignature)
    for i in countdown(signatureDeps.high, 0):
      stack.add((signatureDeps[i], false))

  for node in signatureOrder:
    node.appendNodeSignatureTokens(result)
