#[
  Read complex fields from files and compare them site by site

  ./fieldCmp FILE0 [ FILE1 ... ]
]#

import qex
import os, times

qexInit()

let nfile = paramCount()
if nfile == 0:
  echo "Error: Requires at least one file argument."
  qexAbort()
let lat = paramStr(1).getFileLattice
echo "Lattice size: ",lat
for i in 2..nfile:
  let l = paramStr(i).getFileLattice
  if l != lat:
    echo "File ",paramStr(i)," has lattice size: ",l
    echo "Different from file ",paramStr(0)," that has lattice size: ",lat
    qexAbort()

let nt = lat[^1]
var
  lo = lat.newLayout
  fields = newseq[type(lo.Complex)](nfile)
  traces = newseq[seq[float]](nfile)
  del = lo.Complex
  norms = newseq[float](nfile)
  diff = newseq[float](nfile)
  spatv = 1
for i in 0..<lat.len-1: spatv *= lat[i]

template loadField(field, file:untyped) =
  echo "Loading field from file: ", file
  var reader = lo.newReader file
  reader.read field
  reader.close
  echo "File metadata: ", reader.fileMetadata
  echo "Record metadata: ", reader.recordMetadata

for i in 0..<nfile:
  fields[i] = lo.Complex  # This is for Complex fields only.
  fields[i].loadField paramStr(i+1)

  # per-slice trace
  traces[i].newseq(nt)
  for s in fields[i].sites:
    var t:float
    t := fields[i]{s}.re
    traces[i][lo.coords[3][s]] += t
  traces[i].ranksum

for i in 0..<nfile:
  # norms
  norms[i] = fields[i].norm2
  echo "Field ",i," norm2: ",norms[i]

  if i>0:
    # Compare this and the previous one
    del := fields[i] - fields[i-1]
    diff[i] = del.norm2
    echo "Difference between field ",i-1," and field ",i
    echo "	Norm2 of the field difference: ",diff[i]," relative: ",diff[i]/min(norms[i], norms[i-1])
    var x = 0.0
    for t in 0..<nt:
      x = max(x, abs(traces[i][t]-traces[i-1][t]) / min(traces[i][t], traces[i-1][t]))
    echo "	Max relative difference in time slice traces: ",x

qexFinalize()
