## Standalone numerical NN module tests.
import base/globals
setVLENmax(4) # two-dimensional lattices
import qex, nn
import std/unittest
import ./[pointwise, convolution, masks]

qexInit()
testPointwise[float32]()
testPointwise[float64]()
testConvolution[float32]()
testConvolution[float64]()
testMasks[float32]()
testMasks[float64]()
qexFinalize()
