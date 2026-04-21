import nimib, nimibook
import strutils
import bookutils
nbInit(theme = useNimibook)

let dir = currentSourcePath.parentDir()
let md = dir / "../../BUILD.md"
nbText: translateMd md
nbSave
