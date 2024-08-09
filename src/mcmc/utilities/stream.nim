type
  MCStream* = object
    name*: string
    output*: string

proc newMCStream*(name: string; start: bool = false): MCStream = 
  result = MCStream(name: name)
  case start:
    of true: echo "<begin: " & result.name & ">"
    of false: result.output = "<begin: " & result.name & ">\n"

proc add*(self: var MCStream; text: string) =
  self.output = self.output & text & "\n"

proc finishStream*(self: var MCStream) =
  self.output = self.output & "<end: " & self.name & ">"
  echo self.output