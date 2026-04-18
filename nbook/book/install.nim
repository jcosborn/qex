import nimib, nimibook
import strutils
import bookutils
nbInit(theme = useNimibook)

let dir = currentSourcePath.parentDir()
let md = dir / "../../INSTALL.md"
nbText: translateMd md
nbSave
