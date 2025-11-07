# 2N Low-Storage RK integrators, preserving Lie group structure, and commutator free

import math

#[
delta-form update (used in Luscher, 2010)
advance:
  delta <- alpha_i * delta + h * beta_i * f(y)
  y <- compose(y, delta)
In all cases, alpha_0 is 0.0. No need to zero delta before. And may avoid 1 load based on alpha being 0.

Naming convention: RK2NXY for (X, Y) 2N RK scheme with X stages and Y-th order accurate
]#

# Encapsulated 2N operator capturing state and callbacks
type
  RK2NOp*[Y, T] = object of RootObj
    y*: Y                 ## current state, updated in-place
    delta*: T             ## current delta accumulator
    advance*: proc(y: var Y, delta: var T, alpha: float, beta: float)
  RK2NAdaptiveOp*[Y, T] = object of RK2NOp[Y, T]
    y0*: Y                ## backup buffer for copy-on-reject
    assign*: proc(dst: var Y, src: Y)  ## for backup and restore
    errDelta*: proc(delta: T): float   ## error measure based on delta

type RK2NCoeffs*[S: static[int]] = object
  alpha*: array[S, float]
  beta*: array[S, float]

# Stats returned by adaptive integrator
type RKAdaptiveStats* = object
  steps*: int
  accepts*: int
  rejects*: int
  minH*: float
  maxH*: float

proc rk2nStep*[S: static[int], Y, T](op: var RK2NOp[Y, T], h: float, coeffs: RK2NCoeffs[S]) =
  ## One step using an encapsulated 2N operator (state captured in op)
  for i in 0..<S:
    op.advance(op.y, op.delta, coeffs.alpha[i], h * coeffs.beta[i])

proc rk2n*[S: static[int], Y, T](
    op: var RK2NOp[Y, T],
    coeffs: RK2NCoeffs[S],
    steps: int,
    h: float,
    measureCb: proc(t: float) {.closure.} = nil
  ) =
  ## Run 2N RK for a fixed number of steps with constant step size h.
  ## Avoids rounding the final step to match an exact target time.
  var curT = 0.0
  var i = 0
  while i < steps:
    rk2nStep(op, h, coeffs)
    curT += h
    if not measureCb.isNil:
      measureCb(curT)
    inc i

proc rk2nAdaptive*[S: static[int], Y, T](
    op: var RK2NAdaptiveOp[Y, T],
    coeffs: RK2NCoeffs[S],
    t: float,
    h0: float,
    deltaTol: float,
    safety: float = 0.95,
    maxSteps: int = 100000,
    controllerExp: float = 1.0/3.0,
    measureCb: proc(t: float) {.closure.} = nil
  ): RKAdaptiveStats =
  ## Run adaptive 2N RK using explicit coefficients and copy-on-reject pattern.
  ## Recommend an embedded low order RK 2N scheme, such as RK53_4_2N, where the last stage delta encodes the difference between the embedded and full stage.
  ## controllerExp should be the inverse of the local error order for the embedded scheme, 1/3 for RK53_4_2N.
  var stats: RKAdaptiveStats
  let dir = if t >= 0.0: 1.0 else: -1.0
  var h = abs(h0) * dir
  var curT = 0.0
  while (curT - t) * dir < 0.0:
    if stats.steps >= maxSteps: break
    inc stats.steps
    var hStep = h
    if (curT + hStep - t) * dir > 0.0:
      hStep = t - curT
    op.assign(op.y0, op.y)
    rk2nStep(op, hStep, coeffs)
    let dErr = op.errDelta(op.delta)
    let newHmag = safety * pow(deltaTol / max(dErr, 1e-32), controllerExp) * abs(hStep)
    if dErr <= deltaTol or hStep == (t - curT):
      curT += hStep
      inc stats.accepts
      let hmag = abs(hStep)
      if stats.accepts == 1:
        stats.minH = hmag; stats.maxH = hmag
      else:
        if hmag < stats.minH: stats.minH = hmag
        if hmag > stats.maxH: stats.maxH = hmag
      if not measureCb.isNil:
        measureCb(curT)
    else:
      op.assign(op.y, op.y0)
      inc stats.rejects
    h = dir * newHmag
  stats

# Canonical 2N RK coefficient sets from the literature.
# A,B coefficients are Williamson's notation (Williamson, 1980)
# and used in A-form of Bazavov (Bazavov, 2025).
# delta <- delta + h * B_i * f(y)
# y <- compose(y, A_i * delta)

type RK2NABCoeffs[S: static[int]] = object
  A*: array[S, float]
  B*: array[S, float]

func toRK2NCoeffs[S: static[int]](ab: RK2NABCoeffs[S]): RK2NCoeffs[S] =
  ## Compile-time conversion from (A,B) to (alpha', beta) for 2N delta-form update
  var alpha: array[S, float]
  alpha[0] = 0.0
  for i in 1..<S:
    alpha[i] = (ab.B[i] * ab.A[i]) / ab.B[i-1]
  RK2NCoeffs[S](alpha: alpha, beta: ab.B)

# Third order, Table I, no. 6 and 7, (Williamson, 1980)
# Table B.3 in Bazavov & Chun (2021, 2101.05320)
const RK3W6_AB = RK2NABCoeffs[3](
  A: [0.0, -17.0/32.0, -32.0/27.0],
  B: [1.0/4.0, 8.0/9.0, 3.0/4.0]
)
const RK3W7_AB = RK2NABCoeffs[3](
  A: [0.0, -5.0/9.0, -153.0/128.0],
  B: [1.0/3.0, 15.0/16.0, 8.0/15.0]
)
const RK3W6_2N* = toRK2NCoeffs(RK3W6_AB)
const RK3W7_2N* = toRK2NCoeffs(RK3W7_AB)

# (4,3) four sets from 2506.07359v1 Tables 5–6
const RK43_1_AB = RK2NABCoeffs[4](
  A: [
    0.0,
    -1.0/2.0,
    -13.0/9.0,
    -846.0/625.0
  ],
  B: [
    1.0/4.0,
    2.0/3.0,
    39.0/50.0,
    25.0/78.0
  ]
)
const RK43_2_AB = RK2NABCoeffs[4](
  A: [
    0.0,
    -7.0/15.0,
    -6.0/5.0,
    -145.0/81.0
  ],
  B: [
    1.0/5.0,
    3.0/4.0,
    20.0/27.0,
    3.0/8.0
  ]
)
const RK43_3_AB = RK2NABCoeffs[4](
  A: [
    0.0,
    -29.0/45.0,
    -9.0/5.0,
    -35.0/27.0
  ],
  B: [
    2.0/15.0,
    3.0/4.0,
    10.0/9.0,
    3.0/8.0
  ]
)
const RK43_4_AB = RK2NABCoeffs[4](
  A: [
    0.0,
    -99.0/112.0,
    -16.0/7.0,
    -427.0/648.0
  ],
  B: [
    13.0/28.0,
    12.0/13.0,
    91.0/216.0,
    3.0/13.0
  ]
)

const RK43_1_2N* = toRK2NCoeffs(RK43_1_AB)
const RK43_2_2N* = toRK2NCoeffs(RK43_2_AB)
const RK43_3_2N* = toRK2NCoeffs(RK43_3_AB)
const RK43_4_2N* = toRK2NCoeffs(RK43_4_AB)

# Fourth order AB tables (2101.05320 Table B.4)
# Originally, Carpenter & Kennedy, 1994
const RK4CK_AB = RK2NABCoeffs[5](
  A: [
    0.0,
    -567301805773.0/1357537059087.0,
    -2404267990393.0/2016746695238.0,
    -3550918686646.0/2091501179385.0,
    -1275806237668.0/842570457699.0
  ],
  B: [
    1432997174477.0/9575080441755.0,
    5161836677717.0/13612068292357.0,
    1720146321549.0/2090206949498.0,
    3134564353537.0/4481467310338.0,
    2277821191437.0/14882151754819.0
  ]
)
# Originally, Berland, Bogey, Bailly, 2005
const RK4BBB_AB = RK2NABCoeffs[6](
  A: [
    0.0,
    -0.737101392796,
    -1.634740794341,
    -0.744739003780,
    -1.469897351522,
    -2.813971388035
  ],
  B: [
    0.032918605146,
    0.823256998200,
    0.381530948900,
    0.200092213184,
    1.718581042715,
    0.27
  ]
)
const RK4CK_2N* = toRK2NCoeffs(RK4CK_AB)
const RK4BBB_2N* = toRK2NCoeffs(RK4BBB_AB)

# (5,3) families — four sets from 2506.07359v1 Tables 7–10
const RK53_1_AB = RK2NABCoeffs[5](
  A: [
    0.0,
    -17.0/32.0,
    -9856.0/5625.0,
    -1127375.0/329171.0,
    -4913.0/8800.0
  ],
  B: [
    1.0/4.0,
    136.0/225.0,
    1100.0/1139.0,
    289.0/880.0,
    10.0/47.0
  ]
)
const RK53_2_AB = RK2NABCoeffs[5](
  A: [
    0.0,
    -9.0/16.0,
    -62032.0/41503.0,
    5929.0/9234.0,
    -45.0/98.0
  ],
  B: [
    1.0/4.0,
    36.0/49.0,
    847.0/3078.0,
    3.0/14.0,
    7.0/43.0
  ]
)
const RK53_3_AB = RK2NABCoeffs[5](
  A: [
    0.0,
    -5.0/9.0,
    -14.0/9.0,
    -36.0/25.0,
    -261.0/625.0  # originally -8/25 from Table 9 of the paper 2506.07359v1
  ],
  B: [
    2.0/9.0,
    5.0/8.0,
    18.0/25.0,
    8.0/25.0,
    25.0/192.0
  ]
)
## The first 4 stages of RK53_4_AB is 2nd
const RK53_4_AB = RK2NABCoeffs[5](
  A: [
    0.0,
    -5.0/8.0,
    -4.0/3.0,
    -3.0/4.0,
    -8.0/5.0
  ],
  B: [
    1.0/4.0,
    2.0/3.0,
    1.0/2.0,
    2.0/5.0,
    1.0/9.0
  ]
)

const RK53_1_2N* = toRK2NCoeffs(RK53_1_AB)
const RK53_2_2N* = toRK2NCoeffs(RK53_2_AB)
const RK53_3_2N* = toRK2NCoeffs(RK53_3_AB)
const RK53_4_2N* = toRK2NCoeffs(RK53_4_AB)

# (6,4) from RK2N64 (high-precision decimals), 2506.07359v1, Appendix D
const RK64_AB = RK2NABCoeffs[6](
  A: [
    0.0,
    -7.371013927959100015085736294563710861301655e-01,
    -1.634740794340906961222612899974121227203739e+00,
    -7.447390037800703313971792823734483498376512e-01,
    -1.469897351521944371244484234187043583134644e+00,
    -2.813971388035238894872690695659944758090490e+00
  ],
  B: [
    3.291860514560574016139360757085052620500596e-02,
    8.232569981988439778822317874254015260794315e-01,
    3.815309489002858170631520216481864120871775e-01,
    2.000922131840258454393248810001898523823106e-01,
    1.718581042714403494253985915871400632540402e+00,
    2.700000000000000000000000000000000000000000e-01
  ]
)
const RK64_2N* = toRK2NCoeffs(RK64_AB)

when isMainModule:
  import maths
  import maths/matrixFunctions
  import maths/matlog
  import maths/groupOps
  import maths/matrixConcept
  import maths/complexNumbers
  import strutils, sequtils

  proc fmt(x: float): string = formatFloat(x, ffScientific, 6)

  echo "SU(3) flow tests"

  # SU(3) helpers and test ODE: d/dt Y = - P(H Y) Y
  type Cmplx[T] = ComplexType[T]
  type Mat3C[T] = MatrixArray[3,3,Cmplx[T]]

  # identity in SU(3)
  proc eye3[T](): Mat3C[T] =
    var I: Mat3C[T]
    I := 0
    I[0,0].re = 1.0
    I[1,1].re = 1.0
    I[2,2].re = 1.0
    I

  # Build SU(3) test and ops
  proc buildSu3Test(): tuple[H: Mat3C[float], Y0: Mat3C[float]] =
    var H: Mat3C[float]
    H[0,0].re= 0.7;  H[0,0].im= 0.1
    H[0,1].re = -0.3;  H[0,1].im = 0.2
    H[0,2].re = 0.5;  H[0,2].im = -0.4
    H[1,0].re = 0.2;  H[1,0].im = -0.6
    H[1,1].re = -0.8;  H[1,1].im = 0.3
    H[1,2].re= 0.1;  H[1,2].im= 0.9
    H[2,0].re = -0.4;  H[2,0].im = 0.7
    H[2,1].re = 0.6;  H[2,1].im = -0.2
    H[2,2].re= 0.3;  H[2,2].im= 0.5
    # Build su(3) algebra element via generator coefficients
    var v: VectorArray[8, float]
    v[0] = 0.30
    v[1] = -0.10
    v[2] = 0.05
    v[3] = -0.20
    v[4] = 0.15
    v[5] = 0.40
    v[6] = -0.35
    v[7] = 0.10
    let A = suFromVec(v)
    let Y0 = exp(A)
    (H, Y0)

  let (H, Y0) = buildSu3Test()

  proc su3Ops(y: var Mat3C[float], d: var Mat3C[float], a: float, b: float) =
    var Fy: Mat3C[float]
    projectTAH(Fy, H * y)
    d := a*d - b*Fy
    y := exp(d) * y

  proc su3OpsScaled(y: var Mat3C[float], d: var Mat3C[float], a: float, b: float) =
    const beta = 1.2
    var Fy: Mat3C[float]
    let Hy = beta * H * y
    projectTAH(Fy, Hy)
    let s = exp(Hy.trace.re / 3.0)
    d := a*d - b*(s * Fy)
    y := exp(d) * y

  var aopUnscaled = RK2NAdaptiveOp[Mat3C[float], Mat3C[float]](
    y: Y0,
    delta: default(Mat3C[float]),
    y0: default(Mat3C[float]),
    advance: su3Ops,
    assign: (proc(dst: var Mat3C[float], src: Mat3C[float]) = dst = src),
    errDelta: (proc(d: Mat3C[float]): float = sqrt(norm2(d)))
  )

  var aopScaled = RK2NAdaptiveOp[Mat3C[float], Mat3C[float]](
    y: Y0,
    delta: default(Mat3C[float]),
    y0: default(Mat3C[float]),
    advance: su3OpsScaled,
    assign: (proc(dst: var Mat3C[float], src: Mat3C[float]) = dst = src),
    errDelta: (proc(d: Mat3C[float]): float = sqrt(norm2(d)))
  )

  proc integrateRef(op: var RK2NAdaptiveOp[Mat3C[float], Mat3C[float]], tEnd: float, hRef: float): Mat3C[float] =
    op.y = Y0
    let steps = max(1, int(ceil(tEnd / hRef)))
    let heff = tEnd / float(steps)
    rk2n(op, RK64_2N, steps, heff)
    op.y

  let
    tEnd = 3.0
    Yref = aopUnscaled.integrateRef(tEnd, tEnd/65536.0)
    tEndScaled = 1.0
    YrefScaled = aopScaled.integrateRef(tEndScaled, tEndScaled/65536.0)

  proc fitPower(hs: openArray[float], es: openArray[float]): tuple[p, C: float] =
    ## Linear fit in log-space: ln(err) ≈ ln(C) + p ln(h)
    ##   a n + p Σ x_i = Σ y_i
    ##   a Σ x_i + p Σ x_i^2 = Σ x_i y_i
    ## p = (n Σ x_i y_i − (Σ x_i)(Σ y_i)) / (n Σ x_i^2 − (Σ x_i)^2)
    ## a = (Σ y_i − p Σ x_i)/n
    var sx, sy, sxx, sxy: float
    let n = hs.len.float
    var i = 0
    while i < hs.len:
      let x = ln(hs[i])
      let y = ln(es[i])
      sx += x; sy += y; sxx += x*x; sxy += x*y
      inc i
    let denom = n * sxx - sx * sx
    let p = if abs(denom) > 1e-300: (n * sxy - sx * sy) / denom else: 0.0
    let b = (sy - p * sx) / n
    (p, exp(b))

  proc logDist[T](Uref, U: Mat3C[T]): float =
    ## Right-invariant SU(3) distance: ||log(Uref† U)||_F,
    ## with a robust near-identity approximation.
    let Z = Uref.adj * U
    # Near-identity first-order: Z ≈ I + S, S ∈ su(3), and (Z - Z†)/2 ≈ S
    var S: Mat3C[T]
    projectTAH(S, Z)
    let dFirst = sqrt(max(0.0, (S.adj * S).trace.re))
    # If rotation is sufficiently small, use first-order; else use logm
    let I = eye3[T]()
    let D = Z - I
    let dLin = sqrt(max(0.0, (D.adj * D).trace.re))
    if dLin <= 1e-2:
      dFirst
    else:
      let L = logm(Z)
      sqrt(max(0.0, (L.adj * L).trace.re))

  proc summarize(name: string, coeffs: auto) =
    let hs = @[0.5, 0.25, 0.125, 0.0625]
    var errs = newSeq[float](hs.len)
    var i = 0
    while i < hs.len:
      aopUnscaled.y = Y0
      let steps = max(1, int(ceil(tEnd / hs[i])))
      let heff = tEnd / float(steps)
      rk2n(aopUnscaled, coeffs, steps, heff)
      errs[i] = logDist(Yref, aopUnscaled.y)
      inc i
    let fit = fitPower(hs, errs)
    let ferrs = errs.mapIt(fmt(it))
    echo "SU(3) flow ", name, ": p≈", fit.p, "  C≈", fit.C, "  errs: ", ferrs

  echo "\nUnscaled: dY/dt = - P(H Y) Y"
  summarize("RK3W6", RK3W6_2N)
  summarize("RK3W7", RK3W7_2N)
  summarize("RK43_1", RK43_1_2N)
  summarize("RK43_2", RK43_2_2N)
  summarize("RK43_3", RK43_3_2N)
  summarize("RK43_4", RK43_4_2N)
  summarize("RK4CK", RK4CK_2N)
  summarize("RK4BBB", RK4BBB_2N)
  summarize("RK53_1", RK53_1_2N)
  summarize("RK53_2", RK53_2_2N)
  summarize("RK53_3", RK53_3_2N)
  summarize("RK53_4", RK53_4_2N)
  summarize("RK64", RK64_2N)

  proc summarizeScaled(name: string, coeffs: auto) =
    let hs = @[0.2, 0.1, 0.05, 0.025]
    var errs = newSeq[float](hs.len)
    var i = 0
    while i < hs.len:
      aopScaled.y = Y0
      let steps = max(1, int(ceil(tEndScaled / hs[i])))
      let heff = tEndScaled / float(steps)
      rk2n(aopScaled, coeffs, steps, heff)
      errs[i] = logDist(YrefScaled, aopScaled.y)
      inc i
    let fit = fitPower(hs, errs)
    let ferrs = errs.mapIt(fmt(it))
    echo "SU(3) flow ", name, ": p≈", fit.p, "  C≈", fit.C, "  errs: ", ferrs

  echo "\nScaled: dY/dt = - exp(ReTr(H Y)) P(H Y) Y"
  summarizeScaled("RK3W6", RK3W6_2N)
  summarizeScaled("RK3W7", RK3W7_2N)
  summarizeScaled("RK43_1", RK43_1_2N)
  summarizeScaled("RK43_2", RK43_2_2N)
  summarizeScaled("RK43_3", RK43_3_2N)
  summarizeScaled("RK43_4", RK43_4_2N)
  summarizeScaled("RK4CK", RK4CK_2N)
  summarizeScaled("RK4BBB", RK4BBB_2N)
  summarizeScaled("RK53_1", RK53_1_2N)
  summarizeScaled("RK53_2", RK53_2_2N)
  summarizeScaled("RK53_3", RK53_3_2N)
  summarizeScaled("RK53_4", RK53_4_2N)
  summarizeScaled("RK64", RK64_2N)

  # Adaptive 2N harness using prebuilt op
  proc runAdaptive(label: string, coeffs: auto, op: var RK2NAdaptiveOp[Mat3C[float], Mat3C[float]], Y0, Yref: Mat3C[float], tEnd, h0, tol: float) =
    op.y = Y0
    let res = rk2nAdaptive(op, coeffs, tEnd, h0, tol, maxSteps=500000)
    let err = logDist(Yref, op.y)
    echo "Adaptive (", label, "): steps=", res.steps, " acc=", res.accepts, " rej=", res.rejects, "  err=", fmt(err), "  h∈[", res.minH, ", ", res.maxH, "]"

  echo "\nAdaptive RK(5,3)^4 comparisons"
  runAdaptive("unscaled", RK53_4_2N, aopUnscaled, Y0, Yref, tEnd, 0.4, 1e-4)
  runAdaptive("unscaled", RK53_4_2N, aopUnscaled, Y0, Yref, tEnd, 0.4, 1e-10)
  runAdaptive("unscaled", RK53_4_2N, aopUnscaled, Y0, Yref, tEnd, 0.4, 1e-11)
  runAdaptive("unscaled", RK53_4_2N, aopUnscaled, Y0, Yref, tEnd, 0.4, 1e-12)
  runAdaptive("unscaled", RK53_4_2N, aopUnscaled, Y0, Yref, tEnd, 0.4, 1e-13)
  runAdaptive("scaled", RK53_4_2N, aopScaled, Y0, YrefScaled, tEndScaled, 0.2, 1e-4)
  runAdaptive("scaled", RK53_4_2N, aopScaled, Y0, YrefScaled, tEndScaled, 0.2, 1e-10)
  runAdaptive("scaled", RK53_4_2N, aopScaled, Y0, YrefScaled, tEndScaled, 0.2, 1e-11)
  runAdaptive("scaled", RK53_4_2N, aopScaled, Y0, YrefScaled, tEndScaled, 0.2, 1e-12)
  runAdaptive("scaled", RK53_4_2N, aopScaled, Y0, YrefScaled, tEndScaled, 0.2, 1e-13)

  # Fixed-step effort helpers

  # Compare error at similar total stage budgets
  proc runFixedAtStages(label: string, coeffs: auto, mname: string, budgetStages: int) =
    const sPer = coeffs.beta.len
    let steps = max(1, int(round(budgetStages.float / sPer.float)))
    let h = (if label == "scaled": tEndScaled else: tEnd) / steps.float
    var op = if label == "scaled": aopScaled else: aopUnscaled
    op.y = Y0
    rk2n(op, coeffs, steps, h)
    let err = logDist((if label == "scaled": YrefScaled else: Yref), op.y)
    echo "Budgeted (", label, ", ", mname, "): stages=", steps * sPer, " h=", h, " err=", fmt(err)

  proc compareAtBudget(label: string, h0, tol: float) =
    var op = if label == "scaled": aopScaled else: aopUnscaled
    let tEnd = if label == "scaled": tEndScaled else: tEnd
    op.y = Y0
    let res = rk2nAdaptive(op, RK53_4_2N, tEnd, h0, tol, maxSteps=500000)
    let aderr = logDist((if label == "scaled": YrefScaled else: Yref), op.y)
    let ads = res.accepts * RK53_4_2N.beta.len
    echo "\nBudgeted comparison (", label, "): adaptive stages=", ads, " err=", fmt(aderr), "  h∈[", res.minH, ", ", res.maxH, "]"
    runFixedAtStages(label, RK3W6_2N,  "RK3W6",  ads)
    runFixedAtStages(label, RK3W7_2N,  "RK3W7",  ads)
    runFixedAtStages(label, RK43_1_2N, "RK43_1", ads)
    #runFixedAtStages(label, RK43_2_2N, "RK43_2", ads)
    #runFixedAtStages(label, RK43_3_2N, "RK43_3", ads)
    #runFixedAtStages(label, RK43_4_2N, "RK43_4", ads)
    runFixedAtStages(label, RK4CK_2N,  "RK4CK",  ads)
    runFixedAtStages(label, RK4BBB_2N, "RK4BBB", ads)
    runFixedAtStages(label, RK53_1_2N, "RK53_1", ads)
    #runFixedAtStages(label, RK53_2_2N, "RK53_2", ads)
    #runFixedAtStages(label, RK53_3_2N, "RK53_3", ads)
    #runFixedAtStages(label, RK53_4_2N, "RK53_4", ads)
    runFixedAtStages(label, RK64_2N,   "RK64",   ads)

  compareAtBudget("unscaled", 0.4, 1e-4)
  compareAtBudget("unscaled", 0.4, 1e-12)
  compareAtBudget("scaled",   0.2, 1e-4)
  compareAtBudget("scaled",   0.2, 1e-12)

  # Forward/backward reversibility tests (scaled ODE):
  proc forwardBackwardAdaptiveScaled(mname: string, coeffs: auto, tEnd: float, h0: float, tol: float) =
    aopScaled.y = Y0
    let resF = rk2nAdaptive(aopScaled, coeffs, tEnd, h0, tol, maxSteps=500000)
    let resB = rk2nAdaptive(aopScaled, coeffs, -tEnd, h0, tol, maxSteps=500000)
    let err = logDist(Y0, aopScaled.y)
    echo "Fwd/Back Adaptive (scaled, ", mname, "): fwd steps=", resF.steps, " acc=", resF.accepts, " rej=", resF.rejects,
         " h∈[", resF.minH, ", ", resF.maxH, "]  bwd steps=", resB.steps, " acc=", resB.accepts, " rej=", resB.rejects,
         " h∈[", resB.minH, ", ", resB.maxH, "]  err=", fmt(err)

  # Scaling summaries over all RK coefficients
  proc summarizeFwdBackFixedScaled(name: string, coeffs: auto) =
    let hs = @[0.2, 0.1, 0.05, 0.025]
    let tEnd = 2.0
    var errs = newSeq[float](hs.len)
    var i = 0
    while i < hs.len:
      let steps = max(1, int(ceil(tEnd / hs[i])))
      let heff = tEnd / steps.float
      aopScaled.y = Y0
      rk2n(aopScaled, coeffs, steps, heff)
      rk2n(aopScaled, coeffs, steps, -heff)
      errs[i] = logDist(Y0, aopScaled.y)
      inc i
    let fit = fitPower(hs, errs)
    let ferrs = errs.mapIt(fmt(it))
    echo "Fwd/Back Fixed (scaled, ", name, "): p≈", fit.p, "  C≈", fit.C, "  errs: ", ferrs

  proc summarizeFwdBackAdaptiveScaled(name: string, coeffs: auto) =
    let tols = @[1e-4, 1e-6, 1e-8, 1e-10]
    let tEnd = 2.0
    let h0 = 0.2
    var errs = newSeq[float](tols.len)
    var hsEff = newSeq[float](tols.len)
    var i = 0
    while i < tols.len:
      aopScaled.y = Y0
      let resF = rk2nAdaptive(aopScaled, coeffs, tEnd, h0, tols[i], maxSteps=500000)
      let resB = rk2nAdaptive(aopScaled, coeffs, -tEnd, h0, tols[i], maxSteps=500000)
      errs[i] = logDist(Y0, aopScaled.y)
      let accSteps = max(1, resF.accepts + resB.accepts)
      hsEff[i] = 2.0 * tEnd / accSteps.float
      inc i
    let fit = fitPower(hsEff, errs)
    let ferrs = errs.mapIt(fmt(it))
    echo "Fwd/Back Adaptive (scaled, ", name, "): p≈", fit.p, "  C≈", fit.C, "  errs: ", ferrs, "  h_eff: ", hsEff

  echo "\nForward/Backward reversibility (scaled)"
  forwardBackwardAdaptiveScaled("RK53_4", RK53_4_2N, 2.0, 0.2, 1e-5)
  forwardBackwardAdaptiveScaled("RK53_4", RK53_4_2N, 2.0, 0.2, 1e-10)

  # Scaling summaries over all available RK coefficients (fixed-step)
  summarizeFwdBackFixedScaled("RK3W6", RK3W6_2N)
  summarizeFwdBackFixedScaled("RK3W7", RK3W7_2N)
  summarizeFwdBackFixedScaled("RK43_1", RK43_1_2N)
  #summarizeFwdBackFixedScaled("RK43_2", RK43_2_2N)
  #summarizeFwdBackFixedScaled("RK43_3", RK43_3_2N)
  #summarizeFwdBackFixedScaled("RK43_4", RK43_4_2N)
  summarizeFwdBackFixedScaled("RK4CK", RK4CK_2N)
  summarizeFwdBackFixedScaled("RK4BBB", RK4BBB_2N)
  summarizeFwdBackFixedScaled("RK53_1", RK53_1_2N)
  #summarizeFwdBackFixedScaled("RK53_2", RK53_2_2N)
  #summarizeFwdBackFixedScaled("RK53_3", RK53_3_2N)
  #summarizeFwdBackFixedScaled("RK53_4", RK53_4_2N)
  summarizeFwdBackFixedScaled("RK64", RK64_2N)

  summarizeFwdBackAdaptiveScaled("RK53_4", RK53_4_2N)

  # Budgeted comparison for forward/backward (use adaptive RK53_4 to set budget)
  proc runFwdBackAtStages(label: string, coeffs: auto, tEnd: float, mname: string, budgetStages: int) =
    const sPer = coeffs.beta.len
    let steps = max(1, int(round(budgetStages.float / sPer.float)))
    let h = tEnd / steps.float
    var op = (if label == "scaled": aopScaled else: aopUnscaled)
    op.y = Y0
    rk2n(op, coeffs, steps, h)
    rk2n(op, coeffs, steps, -h)
    let err = logDist(Y0, op.y)
    echo "Budgeted Fwd/Back (", label, ", ", mname, "): stages=", steps * sPer, " h=", h, " err=", fmt(err)

  proc compareFwdBackAtBudget(label: string, tEnd, h0, tol: float) =
    var op = (if label == "scaled": aopScaled else: aopUnscaled)
    op.y = Y0
    let resF = rk2nAdaptive(op, RK53_4_2N, tEnd, h0, tol, maxSteps=500000)
    let resB = rk2nAdaptive(op, RK53_4_2N, -tEnd, h0, tol, maxSteps=500000)
    let aderr = logDist(Y0, op.y)
    let adsF = resF.accepts * RK53_4_2N.beta.len
    let adsB = resB.accepts * RK53_4_2N.beta.len
    echo "\nBudgeted Fwd/Back (", label, "): adaptive stages=", adsF, "+", adsB, " err=", fmt(aderr),
         "  fwd h∈[", resF.minH, ", ", resF.maxH, "]  bwd h∈[", resB.minH, ", ", resB.maxH, "]"
    let ads = (adsF+adsB) div 2
    runFwdBackAtStages(label, RK3W6_2N,  tEnd, "RK3W6",  ads)
    runFwdBackAtStages(label, RK3W7_2N,  tEnd, "RK3W7",  ads)
    runFwdBackAtStages(label, RK43_1_2N, tEnd, "RK43_1", ads)
    #runFwdBackAtStages(label, RK43_2_2N, tEnd, "RK43_2", ads)
    #runFwdBackAtStages(label, RK43_3_2N, tEnd, "RK43_3", ads)
    #runFwdBackAtStages(label, RK43_4_2N, tEnd, "RK43_4", ads)
    runFwdBackAtStages(label, RK4CK_2N,  tEnd, "RK4CK",  ads)
    runFwdBackAtStages(label, RK4BBB_2N, tEnd, "RK4BBB", ads)
    runFwdBackAtStages(label, RK53_1_2N, tEnd, "RK53_1", ads)
    #runFwdBackAtStages(label, RK53_2_2N, tEnd, "RK53_2", ads)
    #runFwdBackAtStages(label, RK53_3_2N, tEnd, "RK53_3", ads)
    #runFwdBackAtStages(label, RK53_4_2N, tEnd, "RK53_4", ads)
    runFwdBackAtStages(label, RK64_2N,   tEnd, "RK64",   ads)

  compareFwdBackAtBudget("scaled", 2.0, 0.2, 1e-4)
  compareFwdBackAtBudget("scaled", 2.0, 0.2, 1e-6)
  compareFwdBackAtBudget("scaled", 2.0, 0.2, 1e-7)
  compareFwdBackAtBudget("scaled", 2.0, 0.2, 1e-9)
  compareFwdBackAtBudget("scaled", 2.0, 0.2, 1e-10)

  # Constant-generator step factor check (commuting case):
  # For commuting f, product of exponentials reduces to exp(F * h),
  # with F = sum_i w_i * beta_i, where w_s=1 and w_i = 1 + alpha_{i+1} * w_{i+1}.
  proc commutingStepFactor[S: static[int]](c: RK2NCoeffs[S]): float =
    var w = newSeq[float](S)
    w[S-1] = 1.0
    var i = S-2
    while i >= 0:
      w[i] = 1.0 + c.alpha[i+1] * w[i+1]
      dec i
    var fac = 0.0
    i = 0
    while i < S:
      fac += w[i] * c.beta[i]
      inc i
    fac

  echo "\nCommuting-step factors (should be 1.0 for exact commuting flows):"
  echo "  RK3W6:  ", commutingStepFactor(RK3W6_2N)
  echo "  RK3W7:  ", commutingStepFactor(RK3W7_2N)
  echo "  RK43_1: ", commutingStepFactor(RK43_1_2N)
  echo "  RK43_2: ", commutingStepFactor(RK43_2_2N)
  echo "  RK43_3: ", commutingStepFactor(RK43_3_2N)
  echo "  RK43_4: ", commutingStepFactor(RK43_4_2N)
  echo "  RK4Ck:  ", commutingStepFactor(RK4CK_2N)
  echo "  RK4BBB: ", commutingStepFactor(RK4BBB_2N)
  echo "  RK53_1: ", commutingStepFactor(RK53_1_2N)
  echo "  RK53_2: ", commutingStepFactor(RK53_2_2N)
  echo "  RK53_3: ", commutingStepFactor(RK53_3_2N)
  echo "  RK53_4: ", commutingStepFactor(RK53_4_2N)
  echo "  RK64:   ", commutingStepFactor(RK64_2N)

  # Order-condition diagnostics for 2N (α,β) -> classical RK (a,b,c)
  proc classicalFrom2N[S: static[int]](c2n: RK2NCoeffs[S]): tuple[a: array[S, array[S, float]], b: array[S, float], c: array[S, float]] =
    ## Convert delta-form 2N coefficients (alpha,beta) to classical RK (a,b,c)
    ## based on the recursion:
    ##   delta_i = alpha_i delta_{i-1} + h beta_i k_i,  y^i = y^{i-1} + delta_i
    ## Let delta_i = h sum_{j<=i} d_{i,j} k_j, then
    ##   d_{1,1} = beta_1; d_{i,i} = beta_i; d_{i,j} = alpha_i * d_{i-1,j} for j<i.
    ## Let y^i = y^0 + h sum_{j<=i} aY_{i,j} k_j, with aY_{0,*}=0 and
    ##   aY_{i,j} = aY_{i-1,j} + d_{i,j}.
    ## Then classical RK has a_{i,j} = aY_{i-1,j} for i>j, b_j = aY_{S,j}, c_i = sum_{j<i} aY_{i-1,j}.
    var d: array[S, array[S, float]]
    var aY: array[S+1, array[S, float]] # aY[i] holds coefficients of y^i
    # Initialize d and aY
    for j in 0..S-1:
      d[0][j] = 0.0
      aY[0][j] = 0.0
    # Stage 1
    d[0][0] = c2n.beta[0]
    for j in 0..0:
      aY[1][j] = aY[0][j] + d[0][j]
    # Stages 2..S
    var i = 1
    while i < S:
      # diagonal term
      d[i][i] = c2n.beta[i]
      # subdiagonal terms
      var j = 0
      while j < i:
        d[i][j] = c2n.alpha[i] * d[i-1][j]
        inc j
      # accumulate y^i -> y^{i+1}
      j = 0
      while j <= i:
        aY[i+1][j] = aY[i][j] + d[i][j]
        inc j
      inc i
    # Build classical a,b,c
    var a: array[S, array[S, float]]
    var b: array[S, float]
    var c: array[S, float]
    i = 0
    while i < S:
      # a_{i+1, j+1} = aY_i,j for j < i
      var j = 0
      while j < i:
        a[i][j] = aY[i][j]
        inc j
      # b_j = aY_S,j; c_{i+1} = sum_{j<i} aY_i,j
      b[i] = aY[S][i]
      var s = 0.0
      j = 0
      while j < i:
        s += aY[i][j]
        inc j
      c[i] = s
      inc i
    (a, b, c)

  proc dot(x, y: openArray[float]): float =
    var s = 0.0
    var i = 0
    while i < x.len:
      s += x[i] * y[i]
      inc i
    s

  proc checkOrder(name: string, coeffs: RK2NCoeffs) =
    let (a,b,c) = classicalFrom2N(coeffs)
    var e = newSeq[float](b.len); for i in 0..<e.len: e[i] = 1.0
    var c2 = newSeq[float](c.len); for i in 0..<c2.len: c2[i] = c[i]*c[i]
    var c3 = newSeq[float](c.len); for i in 0..<c3.len: c3[i] = c2[i]*c[i]
    # A c, A c^2, (A c)∘c, A^2 c
    var Ac = newSeq[float](b.len)
    var Ac2 = newSeq[float](b.len)
    var A2c = newSeq[float](b.len)
    for i in 0..<b.len:
      var s1 = 0.0; var s2 = 0.0; var s3 = 0.0
      for j in 0..i-1:
        s1 += a[i][j] * c[j]
        s2 += a[i][j] * c2[j]
        # For A^2 c: sum_k a_{ik} (Ac)_k
        var t = 0.0
        for k in 0..j-1:
          t += a[j][k] * c[k]
        s3 += a[i][j] * t
      Ac[i] = s1
      Ac2[i] = s2
      A2c[i] = s3
    let r1 = dot(b, e) - 1.0
    let r2 = dot(b, c) - 0.5
    let r3 = dot(b, c2) - (1.0/3.0)
    let r4 = dot(b, c3) - 0.25
    let r5 = dot(b, Ac) - (1.0/6.0)
    # Standard order-4 constants (Butcher):
    # b^T c^3 = 1/4; b^T A c = 1/6; b^T A c^2 = 1/12; b^T (A c)∘c = 1/8; b^T A^2 c = 1/24
    let r6 = dot(b, Ac2) - (1.0/12.0)
    # (Ac)∘c
    var AcHad = newSeq[float](Ac.len); for i in 0..<Ac.len: AcHad[i] = Ac[i]*c[i]
    let r7 = dot(b, AcHad) - (1.0/8.0)
    let r8 = dot(b, A2c) - (1.0/24.0)
    echo name, " residuals, 3rd: ", r1, ", ", r2, ", ", r3, ", ", r5, "; 4th: ", r4, ", ", r6, ", ", r7, ", ", r8

  echo "\nOrder-condition diagnostics:"
  checkOrder("RK3W6", RK3W6_2N)
  checkOrder("RK3W7", RK3W7_2N)
  checkOrder("RK43_1", RK43_1_2N)
  checkOrder("RK43_2", RK43_2_2N)
  checkOrder("RK43_3", RK43_3_2N)
  checkOrder("RK43_4", RK43_4_2N)
  checkOrder("RK4Ck", RK4CK_2N)
  checkOrder("RK4BBB", RK4BBB_2N)
  checkOrder("RK53_1", RK53_1_2N)
  checkOrder("RK53_2", RK53_2_2N)
  checkOrder("RK53_3", RK53_3_2N)
  checkOrder("RK53_4", RK53_4_2N)
  checkOrder("RK64", RK64_2N)
