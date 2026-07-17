## Run the shared stout FTHMC tests in 2D with Nc=1.

import base/globals
setDefaultNc(1)
setVLENmax(4)

import tftstout

runFtStoutTests(@[8, 8], 6.5, 0.1)
