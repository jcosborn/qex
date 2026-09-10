## Higher-derivative towers for Nc=1 in 2D. The 1x1 exp kernel takes the
## scalar-exp path (not the poly-squaring scheme), so the exact analytic
## expDeriv branch needs its own coverage.
## Runs the shared path comparisons and derivative towers like tgtower does
## for Nc=3; run alongside tggauge.
## Run with OMP_NUM_THREADS=1 (see tgtower.nim).

#RUNCMD env OMP_NUM_THREADS=1 $RUN1

import base/globals
setDefaultNc(1)

import tgtower

# No smearing: the 2D U(1) towers pass on raw random configurations.
runTowerTests(@[8,8], 13579111u64, 1, 0)
