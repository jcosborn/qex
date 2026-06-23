import math, strutils, unittest

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))

import helpers
import ../[core, scalar]
import ../core/base
import ../functional

let grt = initGraphRuntime()

include functional/gradients
include functional/lambda
include functional/reuse
