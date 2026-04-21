import strutils

let endmark = "<!--SKIP"
let replacements = [
  ("](src/", "](https://github.com/jcosborn/qex/blob/devel/src/"),
  ("](build/", "](https://github.com/jcosborn/qex/blob/devel/build/"),
  ("](bootstrap-travis", "](https://github.com/jcosborn/qex/blob/devel/bootstrap-travis"),
  ("](qex.nimble", "](https://github.com/jcosborn/qex/blob/devel/qex.nimble"),
  #("](INSTALL.md#", "](install.html#"),
  #("](BUILD.md#", "](build.html#"),
  ("INSTALL.md", "install.html"),
  ("BUILD.md", "build.html")
]

proc translateMd*(filename: string): string =
  for line in filename.lines:
    if line.startsWith(endmark):
      return result
    var newline = line.multiReplace(replacements)
    if newline[0..1]=="##":  # add anchor
      let text = 1 + newline.find(' ')
      let anchorName = newline[text..^1].toLower.replace(" ", "-")
      newline = "<a name = \"" & anchorName & "\"></a>\n" & newline
    result &= newline & '\n'
