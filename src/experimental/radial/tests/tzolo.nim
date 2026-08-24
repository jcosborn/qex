#RUNCMD $RUN1

## Zolotarev rational approximation to 1/sqrt(x) -- targets T1.3a, T1.3b.

import std/[math, strformat, unittest]
import ../ops/zolotarev

const
  win = [(1.0, 50.0), (1.0, 200.0)]   ## (smin, smax); the talk's realistic condition numbers
  odd = [3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31]

func err(r: Rat, x: float): float =
  ## e(x) = 1 - sqrt(x) R(x); the minimax quantity, |e| <= maxRelErr on [smin^2, smax^2]
  1.0 - sqrt(x)*ratValue(r, x)

func prodValue(r: Rat, x: float): float =
  ## cst * prod_i (x + zero_i)/(x + pole_i), the form the partial fractions came from
  result = r.cst
  for i in 0..<r.npole: result *= (x + r.zero[i])/(x + r.pole[i])

func logGrid(r: Rat, n: int): seq[float] =
  result = newSeq[float](n)
  let
    lo = 2.0*ln(r.smin)
    dl = 2.0*(ln(r.smax) - ln(r.smin))/float(n-1)
  for i in 0..<n: result[i] = exp(lo + dl*float(i))

func peaks(r: Rat, n: int): seq[float] =
  ## Signed extremum of each constant-sign run of e(x) over n log-spaced samples.
  ## Interior peaks are refined by the vertex of the parabola through the three
  ## samples around the discrete maximum: the grid alone undershoots a peak by
  ## ~(pi h/w)^2/2 with w the lobe width, which is 3e-7 relative at order 31.
  var e = newSeq[float](n)
  for i, x in logGrid(r, n): e[i] = err(r, x)
  var i = 0
  while i < n:
    var
      j = i
      b = i
    while j < n and (e[j] >= 0.0) == (e[i] >= 0.0):
      if abs(e[j]) > abs(e[b]): b = j
      inc j
    var v = e[b]
    if b > 0 and b < n-1:
      let d = e[b-1] - 2.0*e[b] + e[b+1]
      if d != 0.0: v = e[b] - 0.125*(e[b+1] - e[b-1])*(e[b+1] - e[b-1])/d
    result.add v
    i = j

suite "zolotarev":

  test "ellipticK":
    check abs(ellipticK(0.0)/(0.5*PI) - 1.0) < 1e-15
    check abs(ellipticK(1.0/sqrt(2.0))/1.8540746773013719 - 1.0) < 1e-15
    check abs(ellipticK(0.9)/2.2805491384227703 - 1.0) < 1e-15
    check abs(ellipticK(-0.9)/2.2805491384227703 - 1.0) < 1e-15

  test "jacobiSn":
    for u in [0.0, 0.1, 0.7, 1.5, 3.0, 7.0]:
      check abs(jacobiSn(u, 0.0) - sin(u)) < 1e-15      # sn(u,0) = sin u
      check abs(jacobiSn(u, 1.0) - tanh(u)) < 1e-15     # sn(u,1) = tanh u
      check abs(jacobiSn(-u, 0.3) + jacobiSn(u, 0.3)) < 1e-15
    for k in [0.1, 0.5, 0.9, 0.99, 0.9998, 0.9999999]:
      check abs(jacobiSn(ellipticK(k), k) - 1.0) < 1e-14        # quarter period
      check abs(jacobiSn(2.0*ellipticK(k), k)) < 1e-13          # sn(2K) = 0
      check abs(jacobiSn(3.0*ellipticK(k), k) + 1.0) < 1e-13    # sn(3K) = -1

  test "newRat rejects bad input":
    expect ValueError: discard newRat(0.0, 50.0, 11)
    expect ValueError: discard newRat(-1.0, 50.0, 11)
    expect ValueError: discard newRat(50.0, 50.0, 11)
    expect ValueError: discard newRat(50.0, 1.0, 11)
    expect ValueError: discard newRat(1.0, 50.0, 10)
    expect ValueError: discard newRat(1.0, 50.0, 12)
    expect ValueError: discard newRat(1.0, 50.0, 1)
    expect ValueError: discard newRat(1.0, 50.0, 2)
    expect ValueError: discard newRat(1.0, 50.0, -3)
    expect ValueError: discard newRat(1.0, 50.0, 11, 1)   # would report maxRelErr = 0

  test "npole, interlacing, positive residues":     # T1.3a
    for (a, b) in win:
      for n in odd:
        let r = newRat(a, b, n)
        check r.npole == (n-1) div 2
        check r.pole.len == r.npole
        check r.zero.len == r.npole
        check r.res.len == r.npole
        check r.cst > 0.0
        check r.pole[0] > 0.0
        for j in 0..<r.npole:
          check r.pole[j] < r.zero[j]
          if j+1 < r.npole: check r.zero[j] < r.pole[j+1]
          check r.res[j] > 0.0

  test "partial fractions == product form":
    for (a, b) in win:
      for n in odd:
        let r = newRat(a, b, n)
        for x in r.logGrid(2001):
          let p = r.prodValue x
          check abs(r.ratValue(x)/p - 1.0) < 1e-12

  test "endpoint equioscillation e(smin^2) = -e(smax^2)":
    for (a, b) in win:
      for n in odd:
        let r = newRat(a, b, n)
        check abs(r.err(a*a) + r.err(b*b)) < 1e-12

  test "equioscillation scan":                      # T1.3b
    # e(x) is evaluated in double and the residues carry ~1e-14 relative error from the
    # log-space partial fractions, so |e| sits on an absolute noise floor.  Measure that
    # floor instead of guessing it: at order 41 on the 1/50 window the true Zolotarev
    # error is ~1e-19, so whatever maxRelErr reports there is pure roundoff.
    let efloor = 10.0*newRat(1.0, 50.0, 41).maxRelErr
    echo &"  evaluation noise floor 10x{efloor/10.0:.3e}"
    for (a, b) in win:
      for n in odd:
        # nsample large enough that maxRelErr itself resolves the peaks to <1e-8 relative
        let
          r = newRat(a, b, n, 400001)
          p = r.peaks 40001
        check p.len == n+1                          # exactly n+1 alternating extrema
        for i in 1..<p.len:
          check (p[i] >= 0.0) != (p[i-1] >= 0.0)
        for v in p:
          check abs(abs(v) - r.maxRelErr) < 2e-7*r.maxRelErr + efloor

  test "error decreases with order":
    for (a, b) in win:
      var prev = 1.0
      for n in odd:
        let e = newRat(a, b, n).maxRelErr
        check e < prev
        prev = e

  test "sigma^2 rescaling":
    # spec(X^dag X) in [smin^2, smax^2]: poles and zeros carry smax^2, residues smax,
    # cst 1/smax, so e(x) is invariant under sigma -> s sigma.  s is a power of two and
    # smin, smax are exact, so k = smin/smax and every scaled parameter come out to the
    # last bit; maxRelErr does not, it is a grid-sampled max and the grids do not line up.
    const s = 8.0
    for n in odd:
      let
        r = newRat(0.13, 6.5, n)
        q = newRat(s*0.13, s*6.5, n)
      check abs(q.cst/(r.cst/s) - 1.0) < 1e-15
      for j in 0..<r.npole:
        check abs(q.pole[j]/(r.pole[j]*s*s) - 1.0) < 1e-15
        check abs(q.zero[j]/(r.zero[j]*s*s) - 1.0) < 1e-15
        check abs(q.res[j]/(r.res[j]*s) - 1.0) < 1e-15
      for x in r.logGrid(1001):
        check abs(r.err(x) - q.err(x*s*s)) < 1e-15
    # and the window really is [smin^2, smax^2], not [smin, smax]
    let r = newRat(0.13, 6.5, 11)
    for x in [r.smin*r.smin, r.smin*r.smax, r.smax*r.smax]:
      check abs(1.0 - sqrt(x)*r.ratValue(x)) <= r.maxRelErr + 1e-15

  test "hash":
    let r = newRat(1.0, 50.0, 31)
    check newRat(1.0, 50.0, 31).hash == r.hash          # deterministic
    check newRat(1.0, 50.0, 31, 4001).hash == r.hash    # nsample is not part of the fingerprint
    check newRat(1.0, 50.1, 31).hash != r.hash          # window
    check newRat(1.02, 50.0, 31).hash != r.hash
    check newRat(1.0, 50.0, 29).hash != r.hash          # order
    check newRat(0.1, 5.0, 31).hash != r.hash           # same k, different sigma scale

  test "talk orders 31 (action) and 11 (force)":
    for b in [50.0, 100.0, 200.0, 500.0, 1000.0]:
      for n in [31, 11]:
        let r = newRat(1.0, b, n)
        echo &"  smin/smax 1/{b:<6.0f} order {n:2d} npole {r.npole:2d} " &
             &"maxRelErr {r.maxRelErr:.4e} maxAbsErr {r.maxAbsErr:.4e} " &
             &"pole [{r.pole[0]:.4e}, {r.pole[^1]:.4e}] hash {r.hash:#x}"
        check r.npole == (n-1) div 2
    check newRat(1.0, 50.0, 31).npole == 15    # slide 6: "n = 31 (15 poles)"
    check newRat(1.0, 50.0, 11).npole == 5     # slide 6 says "6 poles": 5 poles + the constant
