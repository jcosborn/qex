import math, strutils, unittest

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))

import helpers
import ../[core, scalar, multi]
from ../core/grad_cache import findGrad

let grt = initGraphRuntime()

include core/scalar_basic
include core/scalar_d2
include core/cond
include core/multi
