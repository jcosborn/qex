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
  let eps0 = 1e-14
  # pAccept = exp(hOld - hNew): its relative error is the absolute error of
  # hNew - hOld, which the eps0 checks bound by 4 eps0 (|hOld| + |hNew|); in
  # reldiff units (a quarter of the relative difference) eps0 (|hOld| + |hNew|).
  let eps1 = eps0*(abs(hmc.hOld)+abs(hmc.hNew))
  test(hmc.hOld, -122757.42513634665, eps0)
  test(hmc.hNew, -122755.99837730126, eps0)
  test(hmc.pAccept, 0.240085769784275, eps1)
  test(hmc.rnd, 0.06697195768356323, eps0)
  rankSum err
  if myRank==0:
    doAssert(err == 0)

import examples/generichmc

genericHmcTest(hmc)

finalize()
