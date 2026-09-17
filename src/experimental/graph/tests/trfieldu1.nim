import base/globals
setDefaultNc(1)
import qex
import std/unittest
import trfield

qexInit()
runBridges(float32)
runBridges(float64)
qexFinalize()
