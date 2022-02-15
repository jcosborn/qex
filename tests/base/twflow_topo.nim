import qex, gauge/wflow
import testutils

proc EQ(g:auto,loop:int):auto =
  let
    f = g.fmunu loop
    (es,et) = f.densityE
    q = f.topoQ
  return [es,et,q]

suite "Test Wilson flow and topological charge":
  qexInit()
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

  test "topo":
    const rl = [
      [0.0,0.0,0.0],
      [1.053685180326327, 1.057184990709729, 0.3484160134779718],
      [0.0,0.0,0.0],
      [2.391841002230686, 2.402675980439514, 0.8365938024236883],
      [4.727855394069386, 4.765561700399072, 1.249916993983875],
      [3.406323101589064, 3.429186597981905, 1.045437791093939]]
    for loop in [1,3,4,5]:
      check(g.EQ(loop) ~ rl[loop])

  test "wflow fine topo":
    const rl = [
      [0.0,0.0,0.0],
      [1.242963724972141, 1.248154190817798, 0.5874265544292369],
      [0.0,0.0,0.0],
      [2.805927859637404, 2.820786889388478, 1.393005673033181],
      [4.791267459839846, 4.827373059918503, 1.965703793198889],
      [3.681453306961074, 3.706279375107263, 1.673630922959521]]
    g.gaugeFlow(20, 0.005):
      discard
    for loop in [1,3,4,5]:
      check(g.EQ(loop) ~ rl[loop])

  test "wflow coarse topo":
    const rl = [
      [0.0,0.0,0.0],
      [1.292616079945679, 1.30148754005825, 0.5424518406094369],
      [0.0,0.0,0.0],
      [2.882075796510504, 2.905993460018066, 1.151042571865504],
      [4.329359048051327, 4.375099479223682, 1.312160447530413],
      [3.526486784704142, 3.560799656228135, 1.220486360255304]]
    g.gaugeFlow(1, 0.1):
      discard
    for loop in [1,3,4,5]:
      check(g.EQ(loop) ~ rl[loop])

  qexFinalize()
