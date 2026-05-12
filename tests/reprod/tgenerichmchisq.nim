{.define: HisqSmearing.}
const defaultLat = @[8, 8, 8, 32]

import qex, hmc/hmcAction, macros

proc reldiff(a,b:float):float = 0.5*abs(a-b)/(abs(a)+abs(b))
macro name(x:untyped):auto = newLit(x.repr)

var err = 0

template test(x,y,e:float) =
  let r = reldiff(x,y)
  var s = "OK    "
  if r > e:
    inc err
    s = "FAILED"
  echo s, " ", name(x), ": ", x, " ", y, " ", r, " ", e

when defined(FUELCompat):
  #let t = [-123189.78139623265,-123193.36720113509,1.0,0.06697195768356323] # 8^4
  discard # need to add 8^3x32
else:
  #let t = [-122744.18900875491,-122747.76198003897,1.0,0.06697195768356323] # 8^4
  let t = [-492089.7593182555,-492104.19773799664,1.0,0.06697195768356323] # 8^3x32

proc genericHmcTest(hmc: HmcAction) =
  let eps = 1e-12
  test(hmc.hOld, t[0], eps)
  test(hmc.hNew, t[1], eps)
  test(hmc.pAccept, t[2], eps)
  test(hmc.rnd, t[3], eps)
  doAssert(err == 0)

include examples/generichmc
