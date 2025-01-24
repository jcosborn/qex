import qex, gauge/wflow
import testutils

proc EQ(g:auto,loop:int):auto =
  let
    f = g.fmunu loop
    (es,et) = f.densityE
    q = f.topoQ
  return [es,et,q]

qexInit()

suite "Test Wilson flow and topological charge":
  echo "rank ",myRank," / ",nRanks
  threads: echo "thread ",threadNum," / ",numThreads
  CT = 1e-11

  let
    seed = 17u64^13
    lat = @[8,8,8,8]
    lo = lat.newLayout
  var
    g = lo.newGauge
    r = lo.newRNGField(MRG32k3a, seed)
  g.warm(0.4, r)
  if g[0].numberType is float32: CT = 3e-3

  test "topo":
    const rl = [
      [0.0,0.0,0.0],
      [0.9278998428166274, 0.9259837153220379, -0.1798124862963527],
      [0.0,0.0,0.0],
      [2.099099182199596, 2.096760447628166, -0.5037001984505194],
      [3.627286981133991, 3.619461449116142, -0.6989980450588869],
      [2.769773748283728, 2.765062802182978, -0.6040403379495964]]
    for loop in [1,3,4,5]:
      check(g.EQ(loop) ~ rl[loop])

  test "wflow fine topo":
    const rl = [
      [0.0,0.0,0.0],
      [0.6597045206103821, 0.6563289799384344, 0.03799796090979355],
      [0.0,0.0,0.0],
      [1.458814621776772, 1.453027521663594, 0.01066005399628158],
      [2.109560159995052, 2.104682631200975, 0.08902327151785866],
      [1.749675334680512, 1.744341716556223, 0.05315490339585562]]
    g.gaugeFlow(20, 0.005):
      discard
    for loop in [1,3,4,5]:
      check(g.EQ(loop) ~ rl[loop])

  test "wflow coarse topo":
    const rl = [
      [0.0,0.0,0.0],
      [0.3330415059918272, 0.3290982076857988, 0.00251688258620615],
      [0.0,0.0,0.0],
      [0.7073852695739449, 0.6993729505801789, -0.008560219145048894],
      [0.9144248903211708, 0.9058947912540163, 0.02577408143748513],
      [0.8013051494825723, 0.7930601360032297, 0.01073616156391243]]
    g.gaugeFlow(1, 0.1):
      discard
    for loop in [1,3,4,5]:
      check(g.EQ(loop) ~ rl[loop])

qexFinalize()
