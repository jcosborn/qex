import qex
import testutils

qexInit()

suite "Test MRG32k3a":
  echo "rank ", myRank, "/", nRanks
  threads: echo "thread ",threadNum," / ",numThreads
  let seed = 17u64^13

  test "subsequence":
    const res = [
      [0.3268000301845387, 0.1909631348029552, 0.3976696014207036],
      [0.2491408676889959, 0.8109031896264907, 0.4171423534316965]]
    var s: MRG32k3a
    for k in 0..1:
      s.seed(seed, k)
      for i in 0..2:
        check(s.uniform ~ res[k][i])

  test "colorspinor":
    var lo = newLayout([8,8,8,16])
    var r = lo.newRNGField(MRG32k3a, seed)
    var d = lo.DiracFermion
    if d.numberType is float32: CT = 1e-6
    d.uniform r
    check(d.norm2 ~ 65517.83893610391)

qexFinalize()
