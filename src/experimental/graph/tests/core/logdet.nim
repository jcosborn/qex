# logDetJ: factorization-chain log-Jacobian builder.
#
# Toy step y = x + a*x^2/2: dy/dx = 1 + a*x, so ld = ln(1 + a*x) with via = x.
# For a scalar chain the defining identity is grad(u, v) == exp(logDetJ(u, v)).

var ldjHookCalls = 0

proc quadStep(x: Gscalar, a: Gscalar): Gscalar =
  proc fwd(v: Gvalue) =
    let x = Gscalar(v.inputs[0])
    let a = Gscalar(v.inputs[1])
    Gscalar(v).sval = x.sval + 0.5*a.sval*x.sval*x.sval
  proc bwd(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
    let x = Gscalar(z.inputs[0])
    let a = Gscalar(z.inputs[1])
    let up = Gscalar(rootedUpstream(zb, z))
    if i == 0:
      Gvalue(up*(1.0 + a*x))
    else:
      Gvalue(up*(0.5*(x*x)))
  proc ldj(z: Gvalue): tuple[ld, via: Gvalue] =
    inc ldjHookCalls
    let x = Gscalar(z.inputs[0])
    let a = Gscalar(z.inputs[1])
    (Gvalue(ln(1.0 + a*x)), z.inputs[0])
  graphNode(scalarNodeLike(x), @[Gvalue(x), Gvalue(a)],
    Gfunc(forward: fwd, backward: bwd, logdet: ldj, name: "quadStep"),
    "quadStep")

proc quadValue(x, a: float): float = x + 0.5*a*x*x

proc zeroLdVia0(z: Gvalue): tuple[ld, via: Gvalue] =
  ## Unit-determinant declaration: static zero ld, flow through input 0.
  (scalarNodeLike(z).zeroLike, z.inputs[0])

proc ldInput1Via0(z: Gvalue): tuple[ld, via: Gvalue] =
  ## Error-path helper: use input 1 as the erased contribution.
  (z.inputs[1], z.inputs[0])

proc probeStep(inputs: openArray[Gvalue], ldj: GlogdetHook): Gscalar =
  ## Structure-only step for engine error paths; never evaluated.
  graphNode(scalarNodeLike(inputs[0]), inputs,
    Gfunc(logdet: ldj, name: "probeStep"), "probeStep")

suite "logDetJ":
  setup:
    let vv {.used.} = 0.3
    let as3 {.used.} = [0.2, -0.4, 0.35]
    let v {.used.} = grt.toGvalue(vv)

  test "chain sum, gradients, and exp identity":
    var u = v
    var aks: seq[Gscalar]
    for ak in as3:
      aks.add grt.toGvalue(ak)
      u = quadStep(u, aks[^1])
    let ldj = logDetJ(u, v)
    var w = vv
    var lds = 0.0
    var dw = 1.0        # dw_{k-1}/dv
    var dldv = 0.0
    var dwa = 0.0       # dw_{k-1}/da_1
    var dlda = 0.0
    for k, ak in as3:
      let
        d = 1.0 + ak*w
        da = if k == 0: 1.0 else: 0.0   # d a_k/d a_1
      lds += ln(d)
      dldv += ak*dw/d
      dlda += (da*w + ak*dwa)/d
      dwa = d*dwa + da*0.5*w*w
      dw = d*dw
      w = quadValue(w, ak)
    ldj :~ lds
    u :~ w
    grad(u, v) :~ exp(ldj)
    grad(u, v) :~ dw
    # The mirror sums associate differently and the terms cancel, so compare
    # these two absolutely instead of by ULP.
    check abs(grad(ldj, v).eval.sval - dldv) < 1e-14
    check abs(grad(ldj, aks[0]).eval.sval - dlda) < 1e-14

  test "empty chain and intermediate base":
    var u = v
    var mid: Gscalar
    for k, ak in as3:
      u = quadStep(u, grt.toGvalue(ak))
      if k == 0: mid = u
    logDetJ(v, v) :~ 0.0
    check logDetJ(v, v).isStaticZeroLeaf
    var w = quadValue(vv, as3[0])
    var lds = 0.0
    for ak in as3[1..^1]:
      lds += ln(1.0 + ak*w)
      w = quadValue(w, ak)
    logDetJ(u, mid) :~ lds

  test "chain sums and hook results are cached":
    let calls0 = ldjHookCalls
    var u = v
    var mid: Gscalar
    for k, ak in as3:
      u = quadStep(u, grt.toGvalue(ak))
      if k == 0: mid = u
    let l1 = logDetJ(u, v)
    check ldjHookCalls == calls0 + 3
    let l2 = logDetJ(u, v)
    check l1 == l2
    check ldjHookCalls == calls0 + 3
    # A different base builds a different sum from the same hook results.
    discard logDetJ(u, mid)
    check ldjHookCalls == calls0 + 3
    # Extending the chain reuses the cached inner sums.
    let u2 = quadStep(u, grt.toGvalue(0.1))
    discard logDetJ(u2, v)
    check ldjHookCalls == calls0 + 4
    # Resetting the cache rebuilds hook results and sums.
    grt.resetLdjCache
    let l3 = logDetJ(u, v)
    check ldjHookCalls == calls0 + 7
    check l3 != l1
    l3 :~ l1

  test "chains compose across former endpoints":
    # V = g(V0), u = flow(V), u1 = f(u): a wider logDetJ walks through nodes
    # previously used as endpoints, reusing every hook result.
    let v0 = grt.toGvalue(0.31)
    var w = v0
    for ak in [0.21, 0.17]: w = quadStep(w, grt.toGvalue(ak))
    let vmid = w                       # V = g(V0)
    for ak in [0.13, -0.29]: w = quadStep(w, grt.toGvalue(ak))
    let u = w                          # u = flow(V)
    let u1 = quadStep(w, grt.toGvalue(0.11))   # u1 = f(u)
    let
      lg = logDetJ(vmid, v0)
      lm = logDetJ(u, vmid)
      lf = logDetJ(u1, u)
      calls0 = ldjHookCalls
      full = logDetJ(u1, v0)
    check ldjHookCalls == calls0       # hook results are base-independent
    full :~ lg + lm + lf
    grad(u1, v0) :~ exp(full)

  test "leaf updates flow through without rebuilding":
    let a = grt.toGvalue(as3[0])
    let u = quadStep(v, a)
    let ldj = logDetJ(u, v)
    ldj :~ ln(1.0 + as3[0]*vv)
    a.update 0.05
    ldj :~ ln(1.0 + 0.05*vv)
    v.update 0.7
    ldj :~ ln(1.0 + 0.05*0.7)
    a.update as3[0]
    v.update vv

  test "volume-preserving steps contribute exactly zero":
    proc shearStep(x: Gscalar, c: Gscalar): Gscalar =
      # y = x + c: dy/dx = 1.
      proc fwd(z: Gvalue) =
        Gscalar(z).sval = Gscalar(z.inputs[0]).sval + Gscalar(z.inputs[1]).sval
      graphNode(scalarNodeLike(x), @[Gvalue(x), Gvalue(c)],
        Gfunc(forward: fwd, logdet: zeroLdVia0, name: "shearStep"), "shearStep")
    # An all-symplectic chain is a static zero, not a computed 0.0.
    let s2 = shearStep(shearStep(v, grt.toGvalue(0.4)), grt.toGvalue(-0.2))
    let z2 = logDetJ(s2, v)
    check z2.isStaticZeroLeaf
    z2 :~ 0.0
    # Mixed chain: zero steps are elided, the rest see the flowed field.
    let u = quadStep(shearStep(v, grt.toGvalue(0.4)), grt.toGvalue(as3[0]))
    logDetJ(u, v) :~ ln(1.0 + as3[0]*(vv + 0.4))

  test "leapfrog on a phase-space carrier is a static zero chain":
    # Harmonic H = p^2/2 + k q^2/2, KDK leapfrog. Each half-step is a shear
    # of the (q, p) carrier — unit determinant, declared as a static zero —
    # so ln det of the whole trajectory map is structurally zero.
    const
      dt = 0.3
      k = 0.7
      nstep = 4
    proc shear(s: Gmulti, cq, cp: float): Gmulti =
      # (q, p) -> (q + cq*p, p + cp*q); cq*cp == 0 keeps it triangular.
      proc fwd(v: Gvalue) =
        let x = Gmulti(v.inputs[0])
        let q = Gscalar(x.storedSlot(0)).sval
        let p = Gscalar(x.storedSlot(1)).sval
        let z = Gmulti(v)
        Gscalar(z.storedSlot(0)).sval = q + cq*p
        Gscalar(z.storedSlot(1)).sval = p + cp*q
      newMultiOutputNode(
        [s.storedSlot(0), s.storedSlot(1)], [Gvalue(s)],
        Gfunc(forward: fwd, logdet: zeroLdVia0, name: "shearQP"), "shearQP")
    let
      q0 = grt.toGvalue(0.9)
      p0 = grt.toGvalue(-0.4)
      s0 = multiValues("leapfrog start", q0, p0)
    var s = s0
    for step in 0..<nstep:
      s = shear(s, 0.0, -0.5*dt*k)   # half kick
      s = shear(s, dt, 0.0)          # drift
      s = shear(s, 0.0, -0.5*dt*k)   # half kick
    let ld = logDetJ(s, s0)
    check ld.isStaticZeroLeaf
    ld :~ 0.0
    discard s.eval
    # The q-marginal of the same integrator is not a factorization chain:
    # the momentum context of the second step reaches the base around via.
    let
      phalf = p0 - (0.5*dt*k)*q0
      q1 = probeStep([Gvalue(q0), Gvalue(phalf)], zeroLdVia0)
      p3half = phalf - (dt*k)*q1
      q2 = probeStep([Gvalue(q1), Gvalue(p3half)], zeroLdVia0)
    expect GraphValueError:
      discard logDetJ(q2, q0)

  test "logDetJ distributes over cond":
    let
      a1 = grt.toGvalue(as3[0])
      a2 = grt.toGvalue(as3[1])
      a3 = grt.toGvalue(as3[2])
      k = grt.toGvalue(1)
      chainA = quadStep(quadStep(v, a1), a2)
      chainB = quadStep(v, a3)
      f = cond(k, chainA, chainB)
      ldj = logDetJ(f, v)
      w1 = quadValue(vv, as3[0])
      ldA = ln(1.0 + as3[0]*vv) + ln(1.0 + as3[1]*w1)
      ldB = ln(1.0 + as3[2]*vv)
    ldj :~ ldA
    k.update 0
    ldj :~ ldB
    k.update 1
    ldj :~ ldA
    # A branch that is the base contributes zero.
    let g = cond(k, chainB, v)
    let ldg = logDetJ(g, v)
    ldg :~ ldB
    k.update 0
    ldg :~ 0.0
    k.update 1

  test "missing factorization fails loudly":
    let u = quadStep(v, grt.toGvalue(as3[0]))
    let s = u + v   # scalar add declares no factorization
    expect GraphError:
      discard logDetJ(s, v)
    # A chain that dead-ends on another leaf never reaches the base.
    let othr = grt.toGvalue(0.9)
    expect GraphError:
      discard logDetJ(quadStep(othr, grt.toGvalue(0.1)), v)

  test "logDetJ rejects dependency cycles":
    var left: Gscalar
    var right: Gscalar
    proc leftInputView(node: Gvalue, mode: InputWalkMode, visit: GnodeVisit) =
      discard node
      discard mode
      visit right
    proc rightInputView(node: Gvalue, mode: InputWalkMode, visit: GnodeVisit) =
      discard node
      discard mode
      visit left
    proc leftLd(node: Gvalue): tuple[ld, via: Gvalue] =
      (scalarNodeLike(node).zeroLike, Gvalue(right))
    proc rightLd(node: Gvalue): tuple[ld, via: Gvalue] =
      (scalarNodeLike(node).zeroLike, Gvalue(left))

    left = graphNode(
      scalarNodeLike(v),
      newSeq[Gvalue](),
      Gfunc(inputView: leftInputView, logdet: leftLd, name: "cycle left"),
      "cycle left")
    right = graphNode(
      scalarNodeLike(v),
      newSeq[Gvalue](),
      Gfunc(inputView: rightInputView, logdet: rightLd, name: "cycle right"),
      "cycle right")

    expect GraphError:
      discard logDetJ(left, v)

  test "base dependence bypassing via is rejected":
    let u = probeStep([Gvalue(quadStep(v, grt.toGvalue(as3[0]))), Gvalue(v)],
      zeroLdVia0)
    expect GraphValueError:
      discard logDetJ(u, v)

  test "non-square and non-scalar declarations are rejected":
    let i = grt.toGvalue(3)
    expect GraphValueError:   # via has int shape: the step is not square
      discard logDetJ(probeStep([Gvalue(v), Gvalue(i)],
        proc(z: Gvalue): tuple[ld, via: Gvalue] = (z.inputs[0], z.inputs[1])), v)
    expect GraphValueError:   # contributions must restore to Gscalar
      discard logDetJ(probeStep([Gvalue(v), Gvalue(i)],
        proc(z: Gvalue): tuple[ld, via: Gvalue] = (z.inputs[1], z.inputs[0])), v)

  test "multi-step non-scalar declarations are rejected before composition":
    let
      scalarLd = grt.toGvalue(0.2)
      intLd = grt.toGvalue(1)
      malformedInner = probeStep(
        [Gvalue(v), Gvalue(intLd)], ldInput1Via0)
      scalarOuter = probeStep(
        [Gvalue(malformedInner), Gvalue(scalarLd)], ldInput1Via0)
    # The scalar outer term must not call Gscalar.addLike on the int term.
    expect GraphValueError:
      discard logDetJ(scalarOuter, v)

    let
      validInner = probeStep(
        [Gvalue(v), Gvalue(scalarLd)], ldInput1Via0)
      intZero = grt.toGvalue(0).zeroLike
      malformedZeroOuter = probeStep(
        [Gvalue(validInner), intZero], ldInput1Via0)
    # Validate before static-zero elision can hide the malformed outer term.
    expect GraphValueError:
      discard logDetJ(malformedZeroOuter, v)

  test "via must be a backward dependency":
    let stray = grt.toGvalue(0.5)
    let u = probeStep([Gvalue(grt.toGvalue(1.1))],
      proc(z: Gvalue): tuple[ld, via: Gvalue] =
        (Gvalue(toGvalue(z.runtime, 0.0)), Gvalue(stray)))
    expect GraphValueError:
      discard logDetJ(u, v)
