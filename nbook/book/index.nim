import nimib, nimibook
import strutils
import bookutils
nbInit(theme = useNimibook)

let dir = currentSourcePath.parentDir()
let readme = dir / "../../README.md"
nbText: translateMd readme
nbSave
