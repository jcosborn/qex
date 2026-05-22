import math, strutils, unittest

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))

import helpers
import ../[core, scalar, functional]
import ../core/base

let grt = initGraphRuntime()

include functional/gradients
include functional/lambda
include functional/reuse
