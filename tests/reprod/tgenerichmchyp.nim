{.define: HypSmearing.}
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

proc genericHmcTest(hmc: HmcAction) =
  let eps = 1e-12
  test(hmc.hOld, -122757.42513634665, eps)
  test(hmc.hNew, -122758.33371907799, eps)
  test(hmc.pAccept, 1.0, eps)
  test(hmc.rnd, 0.06697195768356323, eps)
  doAssert(err == 0)

include examples/generichmc
