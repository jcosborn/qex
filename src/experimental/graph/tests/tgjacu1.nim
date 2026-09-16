#RUNCMD env OMP_NUM_THREADS=1 $RUN1
import base/globals
setDefaultNc(1)
setVLENmax(4)
import helpers
import tgjac
runJacTests(@[4,4])
