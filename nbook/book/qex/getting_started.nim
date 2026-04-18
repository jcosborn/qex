import nimib, nimibook

nbInit(theme = useNimibook)
#nb.pathToRoot = "../.."  # adjust if needed

nbText: """
## QEX Getting Started

This section demonstrates how to set up QEX in Nim.
"""

# #[
nbCode:
  import qex
  qexInit()
  let lat = [4, 4, 4, 4]
  let lo = newLayout(lat)
  var v = lo.ColorVector()
  v := 1
  if myRank == 0:
    echo "Vector element at site 0:", v[0][0]
  qexFinalize()
# ]#

nbSave
