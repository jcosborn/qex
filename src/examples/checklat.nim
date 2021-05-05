import qex
import os

qexInit()
threads: echo "thread ",threadNum," / ",numThreads

let nargs = paramCount()
if nargs != 1:
  echo "Error: Requires one file argument."
  qexAbort()
let lat = paramStr(1).getFileLattice
echo "Lattice size: ",lat

var (l,g,r) = setupLattice(lat)
let p = g.plaq
echo p

qexFinalize()
