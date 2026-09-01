## XML parsing DSL for QEX
##
## XML parsing DSL that makes it easy to specify what parameters would be read from
## input XML file; includes default values for specified parameters.
## 
## It would be nice to have a corresponding DSL for JSON file using Nim's native
## JSON parsing tools. 
## 
## Author: Curtis Taylor Peterson

import std/[os]
import std/[macros]
import std/[strutils]
import std/[xmlparser, xmltree]

#[ XML printout ]#

type XmlParamObj = object
  name: string
  value: string

var xmlParamInfo = newSeq[XmlParamObj](0)

proc addXmlParam(name, value: string) =
  for p in xmlParamInfo.mitems:
    if p.name == name:
      p.value = value
      return
  xmlParamInfo.add XmlParamObj(name: name, value: value)

proc echoXmlX: string =
  result = "XML parameters:\n"
  for p in xmlParamInfo: result &= "  " & p.name & ": " & p.value & "\n"
  result.removeSuffix '\n'

template echoXml* = 
  mixin echo
  echo echoXmlX()

#[ XML parsing ]#

proc findPath(root: XmlNode; path: openArray[string]): XmlNode =
  result = root
  for p in path:
    result = result.child(p)
    if result.isNil: return

proc fromXml(s: string; _: string): string = s
proc fromXml(s: string; _: int): int = parseInt(s)
proc fromXml(s: string; _: float): float = parseFloat(s)
proc fromXml(s: string; _: bool): bool = parseBool(s)
proc fromXml[T: enum](s: string; _: T): T = parseEnum[T](s)
proc fromXml[T](s: string; _: seq[T]): seq[T] =
  for x in s.replace(",", " ").splitWhitespace():
    result.add fromXml(x, default(T))

proc xmlParam[T](root: XmlNode; path: openArray[string]; defaultValue: T): T = 
  let node = findPath(root, path)

  result = defaultValue
  if not node.isNil:
    result = fromXml(node.innerText.strip, defaultValue)

  addXmlParam(path.join("/"), $result)

proc xmlAttr[T](
  root: XmlNode;
  path: openArray[string];
  name: string;
  defaultValue: T
): T = 
  let node = findPath(root, path)

  result = defaultValue
  if not node.isNil:
    let value = node.attr(name)
    if value.len > 0: result = fromXml(value, defaultValue)

  let prefix = path.join("/")
  let attrPath = if prefix.len > 0: prefix & "/@" & name else: "@" & name
  addXmlParam(attrPath, $result)

#[ DSL AST parsing ]#

proc pathNode(path: seq[string]): NimNode {.compileTime.} = 
  result = newNimNode(nnkBracket)
  for name in path: result.add newLit(name)

proc isAttrCall(node: NimNode): bool {.compileTime.} =
  return (
    node.kind in CallNodes and
    node.len > 0 and
    node[0].kind in {nnkIdent, nnkSym} and
    node[0].eqIdent("attr")
  )

proc emitXml(
  root: NimNode;
  path: seq[string];
  body: NimNode;
  output: NimNode
) {.compileTime.} =
  # walk the AST of the letXml block and emit code to read XML parameters
  for node in body:
    if node.kind == nnkAsgn: # assignment branch
      let name = node[0]
      let value = node[1]

      if name.kind notin {nnkIdent, nnkSym}: 
        error "letXml expects `name = default`: ", node
    
      if isAttrCall(value): # XML attribute branch
        if value.len != 3: error "attr expects `attr(name, default)`: ", value

        let p = pathNode(path)
        let attrName = value[1]
        let defaultValue = value[2]

        output.add quote do:
          let `name` = xmlAttr(`root`, `p`, `attrName`, `defaultValue`)
      else: # XML parameter branch
        var childPath = path
        childPath.add name.strVal

        let p = pathNode(childPath)
        output.add quote do:
          let `name` = xmlParam(`root`, `p`, `value`)
    elif node.kind in CallNodes: # nested block branch
      let section = node[0]
      if node.len != 2 or node[^1].kind != nnkStmtList:
        error "malformed letXml subtree: ", node
      if section.kind notin {nnkIdent, nnkSym}: 
        error "letXml subtree name must be an identifier: ", section

      var childPath = path
      childPath.add section.strVal

      emitXml(root, childPath, node[^1], output)
    elif node.kind == nnkCommentStmt: discard
    else: error "unsupported syntax in letXml: ", node

macro letXml*(xml: typed; body: untyped): untyped = 
  ## Bind XML values to local variables
  ## 
  ## # Example XML (input.xml):
  ## # <?xml version="1.0"?>
  ## # <qex>
  ## #   <action>
  ## #     <beta>7.5</beta>
  ## #     <mass>0.005</mass>
  ## #   </action>
  ## #   <hmc>
  ## #     <tau>1.0</tau>
  ## #     <steps>10</steps>
  ## #   </hmc>
  ## # </qex>
  ## 
  ## # Example usage:
  ## let xml = "input.xml"
  ## letXml xml: 
  ##   action:
  ##     beta = 7.5
  ##     mass = 0.005
  ##   hmc:
  ##     tau = 1.0
  ##     steps = 10

  let root = genSym(nskLet, "xml")
  result = newStmtList()

  result.add quote do: 
    let `root` = loadXml(`xml`)

  emitXml(root, @[], body, result)

when isMainModule:
  let testFile = getTempDir() / "qex_xml_test.xml"

  writeFile(
    testFile, 
"""
<qex>
  <lattice>
    <dims>16 16 16 32</dims>
  </lattice>
  <action>
    <gauge>
      <beta>7.5</beta>
    </gauge>
    <ferm>
      <mass>0.005</mass>
    </ferm>
  </action>
  <hmc name="OMF2">
    <tau>1.0</tau>
    <steps>10</steps>
  </hmc>
</qex>
"""
  )

  letXml testFile:
    lattice:
      dims = @[8, 8, 8, 8]
    action:
      gauge:
        beta = 6.0
      ferm:
        mass = 0.001
    hmc:
      integrator = attr("name", "leapfrog")
      tau = 0.5
      steps = 5
      missing = 42

  doAssert dims == @[16, 16, 16, 32]
  doAssert beta == 7.5
  doAssert mass == 0.005
  doAssert integrator == "OMF2"
  doAssert tau == 1.0
  doAssert steps == 10
  doAssert missing == 42

  echoXml()

  echo "xmlparams tests passed"

  removeFile(testFile)
