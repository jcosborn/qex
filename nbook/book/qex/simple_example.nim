import nimib, nimibook

nbInit(theme = useNimibook)

nbText: """
## QEX Simple Example
"""

nbCode:
  import qex
  qexInit()
  let lat = [4, 4, 4, 4]     # lattice size, used 'let' since we won't change it
  let lo = newLayout(lat)    # layout object, distributes lattice across ranks
  var v = lo.ColorVector()   # default lattice color vector field (typically Nc=3, double precision)
  threads:                   # start thread block
    v := 1                   # set all components of vector to 1
  echo "norm2 of vector:", v.norm2
  echo "3 * volume: ", 3 * lo.physVol
  qexFinalize()

nbSave
