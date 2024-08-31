type
  MCStream* = object
    name*: string
    output*: string

proc newMCStream*(name: string; start: bool = false): MCStream = 
  result = MCStream(name: name)
  discard #echo "<begin: " & result.name & ">"

proc add*(self: var MCStream; text: string) =
  discard #echo self.output & text

proc finishStream*(self: var MCStream) =
  discard #echo self.output & "<end: " & self.name & ">"