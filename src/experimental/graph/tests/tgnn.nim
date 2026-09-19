## Standalone NN graph value and operation tests.
import base/globals
setVLENmax(4) # two-dimensional lattices
import qex
import ../../../nn as numeric
import ../nn as gnn
import ../[core, scalar, functional, plan]
import base/alignedMem
import std/[math, sequtils, unittest]

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))
include gnn/[contracts, storage, precision, lifetime, higher]

qexInit()
let rt = initGraphRuntime()
precision(rt)
contracts[float32](rt)
contracts[float64](rt)
storage[float32](rt)
storage[float64](rt)
lifetime[float32]()
lifetime[float64]()
higher[float32]()
higher[float64]()
qexFinalize()
