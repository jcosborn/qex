## Subset stout graph ops. For S = (parity,dir),
##
##   F_x = projectTAH(W_x ds_x†)
##   W'_x = exp(alpha F_x) W_x,  x in S;  W'_x = W_x otherwise
##   M_x = alpha W_x ds_x†
##   log J_S = sum(x in S) expProjMulLogJac(M_x).
##
## The subset Jacobian is link-local. U(1) and SU(3) dispatch by matrix size.

import ../[core, scalar, multi]
import ../support/op
import layout, gauge, physics/qcdTypes
import shared, basic_ops
from action/ops import Gactcoeff, gaugeActionDeriv
import maths/groupOps, maths/matrixFunctions

proc gaugeGradSlot(x: Gmulti, i: int): Ggauge =
  ## View slot i without copying; x remains an evaluation dependency.
  let slot = Ggauge(x.storedSlot(i))
  let view = Ggauge(runtime: slot.runtime, gval: slot.gval)
  proc fwd(v: Gvalue) =
    Ggauge(v).gval = Ggauge(Gmulti(v.inputs[0]).storedSlot(i)).gval
  proc bwd(zb: Gvalue, z: Gvalue, j: int, input: Gvalue): Gvalue =
    raiseUnsupportedPath("stout gradient slot view backward", "higher stout derivatives")
  graphNode(view, @[Gvalue(x)], Gfunc(forward: fwd, backward: bwd, name: "stoutGradView"), "stoutGradView")

# --- stoutUpdate: fused subset map -------------------------------------------

type GstoutUpdate = ref object of Ggauge
  expa: DLatticeColorMatrixV

method newOneOf(x: GstoutUpdate): Gvalue =
  let g = x.gval.newOneOf
  g.zeroGaugeStorage
  GstoutUpdate(
    runtime: x.runtime,
    gval: g,
    expa: x.expa.newOneOf).assignStableNodeId

proc stoutUpdateImpl(W, ds: Ggauge, alpha: Gscalar, parity, dir: int): GstoutUpdate =
  ## W' = exp(alpha*projectTAH(W ds†))*W on (parity,dir).
  W.requireSameGaugeShape(ds, "stoutUpdate")
  let
    sub = W.gval.paritySubset(parity)
    other = W.gval.paritySubset(1 - parity)
    args = multiValues("stoutUpdate args", W, ds, alpha)

  proc forward(v: Gvalue) =
    tic("stoutUpdate kernel")
    let args = Gmulti(v.inputs[0])
    let W = Ggauge(args.storedSlot(0))
    let ds = Ggauge(args.storedSlot(1))
    let alpha = Gscalar(args.storedSlot(2))
    let z = GstoutUpdate(v)
    threads:
      for mu in 0..<z.gval.len:
        if mu != dir:
          z.gval[mu] := W.gval[mu]
      for x in other:
        z.gval[dir][x] := W.gval[dir][x]
      for x in sub:
        var A, F, E {.noinit.}: evalType(W.gval[dir][x])
        F.projectTAH(W.gval[dir][x] * ds.gval[dir][x].adj)
        A := alpha.sval * F
        E[] := expAH(A[])
        z.expa[x] := E
        z.gval[dir][x] := E * W.gval[dir][x]
    toc("stoutUpdate kernel end")

  proc gradPair(W, ds: Ggauge, alpha: Gscalar, upstream: Ggauge, update: GstoutUpdate): Gmulti =
    # The update dependency refreshes expa before this pullback.
    let terms = upstream.gaugeAddTerms
    let nterms = terms.len
    var values = @[Gvalue(W), Gvalue(ds), Gvalue(alpha)]
    for term in terms:
      values.add Gvalue(term)
    values.add Gvalue(update)
    let gradArgs = multiValues("stoutUpdateGrad args", values)
    proc kf(v: Gvalue) =
      tic("stoutUpdateGrad kernel")
      let args = Gmulti(v.inputs[0])
      let W = Ggauge(args.storedSlot(0))
      let ds = Ggauge(args.storedSlot(1))
      let alpha = Gscalar(args.storedSlot(2))
      let update = GstoutUpdate(args.storedSlot(3 + nterms))
      let z = Gmulti(v)
      let dW = Ggauge(z.storedSlot(0))
      let dds = Ggauge(z.storedSlot(1))
      template upstreamSum(mu, x: untyped): untyped =
        block:
          var r {.noinit.}: evalType(W.gval[mu][x])
          r := Ggauge(args.storedSlot(3)).gval[mu][x]
          for j in 1..<nterms:
            r += Ggauge(args.storedSlot(3 + j)).gval[mu][x]
          r
      threads:
        for mu in 0..<dW.gval.len:
          if mu != dir:
            for x in dW.gval[mu]:
              dW.gval[mu][x] := upstreamSum(mu, x)
        for x in other:
          dW.gval[dir][x] := upstreamSum(dir, x)
        for x in sub:
          let
            w = W.gval[dir][x]
            d = ds.gval[dir][x]
            b = upstreamSum(dir, x)
          var B, P {.noinit.}: evalType(w)
          B := update.expa[x].adj * b
          expProjectTAHPullback(P[], (alpha.sval * (w * d.adj))[], (B * w.adj)[])
          let H = alpha.sval * P
          dW.gval[dir][x] := B + H * d
          dds.gval[dir][x] := H.adj * w
      toc("stoutUpdateGrad kernel end")
    proc kb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
      raiseUnsupportedPath("stoutUpdate gradient backward", "higher stout update derivatives")
    result = newMultiOutputNode(@[Gvalue(W), Gvalue(ds)], @[Gvalue(gradArgs)], Gfunc(forward: kf, backward: kb, name: "stoutUpdateGrad"), "stoutUpdateGrad")
    Ggauge(result.storedSlot(1)).zeroGaugeStorage

  proc alphaGrad(W, ds: Ggauge, alpha: Gscalar, upstream: Ggauge): Gscalar =
    let terms = upstream.gaugeAddTerms
    let nterms = terms.len
    var values = @[Gvalue(W), Gvalue(ds), Gvalue(alpha)]
    for term in terms:
      values.add Gvalue(term)
    let gradArgs = multiValues("stoutUpdateAlphaGrad args", values)
    proc kf(v: Gvalue) =
      tic("stoutUpdateAlphaGrad kernel")
      let args = Gmulti(v.inputs[0])
      let W = Ggauge(args.storedSlot(0))
      let ds = Ggauge(args.storedSlot(1))
      let alpha = Gscalar(args.storedSlot(2))
      let z = Gscalar(v)
      var da = 0.0
      threads:
        var s = 0.0
        for x in sub:
          let
            w = W.gval[dir][x]
            d = ds.gval[dir][x]
          var b {.noinit.}: evalType(w)
          b := Ggauge(args.storedSlot(3)).gval[dir][x]
          for j in 1..<nterms:
            b += Ggauge(args.storedSlot(3 + j)).gval[dir][x]
          var F, D {.noinit.}: evalType(w)
          F.projectTAH(w * d.adj)
          D := expDeriv(alpha.sval * F, b * w.adj)
          let r = redot(D, F)
          s += simdSum(r)
        s.threadRankSum
        threadSingle: da = s
      z.sval = da
      toc("stoutUpdateAlphaGrad kernel end")
    proc kb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
      raiseUnsupportedPath("stoutUpdate alpha gradient backward", "higher stout update derivatives")
    graphNode(scalarNodeLike(alpha), @[Gvalue(gradArgs)], Gfunc(forward: kf, backward: kb, name: "stoutUpdateAlphaGrad"), "stoutUpdateAlphaGrad")

  proc backward(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
    let args = Gmulti(z.inputs[0])
    let W = Ggauge(args[0])
    let ds = Ggauge(args[1])
    let alpha = Gscalar(args[2])
    let upstream = Ggauge(rootedUpstream(zb, z))
    let pair = gradPair(W, ds, alpha, upstream, GstoutUpdate(z))
    Gvalue(multiValues("stoutUpdate input gradients", gaugeGradSlot(pair, 0), gaugeGradSlot(pair, 1), alphaGrad(W, ds, alpha, upstream)))

  let g = W.gval.newOneOf
  g.zeroGaugeStorage
  graphNode(
    GstoutUpdate(
      runtime: W.runtime,
      gval: g,
      expa: W.gval[dir].newOneOf),
    @[Gvalue(args)],
    Gfunc(forward: forward, backward: backward, name: "stoutUpdate"),
    "stoutUpdate")

proc stoutUpdate*(W, ds: Ggauge, alpha: Gscalar, parity, dir: int): Ggauge =
  ## W' = exp(alpha*projectTAH(W ds†))*W on (parity,dir).
  stoutUpdateImpl(W, ds, alpha, parity, dir)

# --- stoutLogDetJ: per-substep log-Jacobian -----------------------------------

proc stoutLogDetJ*(W, ds: Ggauge, alpha: Gscalar, parity, dir: int): Gscalar =
  ## sum(x in S) log det J(M_x), M_x = alpha W_x ds_x†.
  ## ds must be the derivative used by the matching subset update.
  W.requireSameGaugeShape(ds, "stoutLogDetJ")
  let args = multiValues("stoutLogDetJ args", W, ds, alpha)
  let sub = W.gval.paritySubset(parity)

  proc forward(v: Gvalue) =
    tic("stoutLogDetJ kernel")
    let args = Gmulti(v.inputs[0])
    let W = Ggauge(args.storedSlot(0))
    let ds = Ggauge(args.storedSlot(1))
    let alpha = Gscalar(args.storedSlot(2))
    let z = Gscalar(v)
    var res = 0.0
    threads:
      var s = 0.0
      for x in sub:
        var M {.noinit.}: evalType(W.gval[dir][x])
        M := alpha.sval * (W.gval[dir][x] * ds.gval[dir][x].adj)
        s += simdSum(expProjMulLogJac(M[]))
      s.threadRankSum
      threadSingle: res = s
    z.sval = res
    toc("stoutLogDetJ kernel end")

  proc gradPair(W, ds: Ggauge, alpha, upstream: Gscalar): Gmulti =
    let gradArgs = multiValues("stoutLogDetJgrad args", W, ds, alpha, upstream)
    proc kf(v: Gvalue) =
      tic("stoutLogDetJgrad kernel")
      let args = Gmulti(v.inputs[0])
      let W = Ggauge(args.storedSlot(0))
      let ds = Ggauge(args.storedSlot(1))
      let alpha = Gscalar(args.storedSlot(2))
      let upstream = Gscalar(args.storedSlot(3))
      let z = Gmulti(v)
      let dW = Ggauge(z.storedSlot(0))
      let dds = Ggauge(z.storedSlot(1))
      threads:
        for x in sub:
          let
            w = W.gval[dir][x]
            d = ds.gval[dir][x]
          var q, M, G {.noinit.}: evalType(w)
          q := w * d.adj
          M := alpha.sval * q
          expProjMulLogJacGrad(G[], M[])
          let H = (upstream.sval * alpha.sval) * G
          dW.gval[dir][x] := H * d
          dds.gval[dir][x] := H.adj * w
      toc("stoutLogDetJgrad kernel end")
    proc kb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
      raiseUnsupportedPath("stoutLogDetJ grad kernel backward", "stout log-Jacobian gradient is not differentiated")
    let kfn = Gfunc(forward: kf, backward: kb, name: "stoutLogDetJgrad")
    result = newMultiOutputNode(@[Gvalue(W), Gvalue(ds)], @[Gvalue(gradArgs)], kfn, "stoutLogDetJgrad")
    Ggauge(result.storedSlot(0)).zeroGaugeStorage
    Ggauge(result.storedSlot(1)).zeroGaugeStorage

  proc alphaGrad(W, ds: Ggauge, alpha, upstream: Gscalar): Gscalar =
    let gradArgs = multiValues("stoutLogDetJalphaGrad args", W, ds, alpha, upstream)
    proc kf(v: Gvalue) =
      tic("stoutLogDetJalphaGrad kernel")
      let args = Gmulti(v.inputs[0])
      let W = Ggauge(args.storedSlot(0))
      let ds = Ggauge(args.storedSlot(1))
      let alpha = Gscalar(args.storedSlot(2))
      let upstream = Gscalar(args.storedSlot(3))
      let z = Gscalar(v)
      var da = 0.0
      threads:
        var s = 0.0
        for x in sub:
          let
            w = W.gval[dir][x]
            d = ds.gval[dir][x]
          var q, M, G {.noinit.}: evalType(w)
          q := w * d.adj
          M := alpha.sval * q
          expProjMulLogJacGrad(G[], M[])
          let r = redot(G, q)
          s += upstream.sval * simdSum(r)
        s.threadRankSum
        threadSingle: da = s
      z.sval = da
      toc("stoutLogDetJalphaGrad kernel end")
    proc kb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
      raiseUnsupportedPath("stoutLogDetJ alpha gradient backward", "higher stout log-Jacobian derivatives")
    graphNode(scalarNodeLike(alpha), @[Gvalue(gradArgs)], Gfunc(forward: kf, backward: kb, name: "stoutLogDetJalphaGrad"), "stoutLogDetJalphaGrad")

  proc backward(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
    let args = Gmulti(z.inputs[0])
    let W = Ggauge(args[0])
    let ds = Ggauge(args[1])
    let alpha = Gscalar(args[2])
    let upstream = Gscalar(rootedUpstream(zb, z))
    let pair = gradPair(W, ds, alpha, upstream)
    Gvalue(multiValues("stoutLogDetJ input gradients", gaugeGradSlot(pair, 0), gaugeGradSlot(pair, 1), alphaGrad(W, ds, alpha, upstream)))

  let f = Gfunc(forward: forward, backward: backward, name: "stoutLogDetJ")
  graphNode(scalarNodeLike(W), @[Gvalue(args)], f, "stoutLogDetJ")

# --- grouped stout update and log-Jacobian backward --------------------------

proc stoutStepOutputView(v: Gvalue, mode: InputWalkMode, visit: GnodeVisit) =
  case mode
  of iwmBackward:
    visit v.inputs[0]
  of iwmEval, iwmReachable:
    visit v.inputs[0]
    visit v.inputs[1]

proc stoutStepParentInputView(v: Gvalue, mode: InputWalkMode, visit: GnodeVisit) =
  visit v.inputs[0]

proc stoutActionStepParentInputView(v: Gvalue, mode: InputWalkMode, visit: GnodeVisit) =
  visit v.inputs[0]
  visit v.inputs[1]
  visit v.inputs[2]

type GstoutStepPullback = ref object of Ggauge
  hdir: DLatticeColorMatrixV

method newOneOf(x: GstoutStepPullback): Gvalue =
  let g = x.gval.newOneOf
  g.zeroGaugeStorage
  GstoutStepPullback(
    runtime: x.runtime,
    gval: g,
    hdir: x.hdir.newOneOf).assignStableNodeId

type GstoutStepForward = ref object of GstoutUpdate
  m: DLatticeColorMatrixV

method newOneOf(x: GstoutStepForward): Gvalue =
  let g = x.gval.newOneOf
  g.zeroGaugeStorage
  GstoutStepForward(
    runtime: x.runtime,
    gval: g,
    expa: x.expa.newOneOf,
    m: x.m.newOneOf).assignStableNodeId

proc stoutStepForwardImpl(W, ds: Ggauge, alpha: Gscalar, parity, dir: int): GstoutStepForward =
  W.requireSameGaugeShape(ds, "stoutStepForward")
  let
    sub = W.gval.paritySubset(parity)
    other = W.gval.paritySubset(1 - parity)
    args = multiValues("stoutStepForward args", W, ds, alpha)

  proc forward(v: Gvalue) =
    tic("stoutStepForward kernel")
    let
      args = Gmulti(v.inputs[0])
      W = Ggauge(args.storedSlot(0))
      ds = Ggauge(args.storedSlot(1))
      alpha = Gscalar(args.storedSlot(2))
      z = GstoutStepForward(v)
    threads:
      for mu in 0..<z.gval.len:
        if mu != dir:
          z.gval[mu] := W.gval[mu]
      for x in other:
        z.gval[dir][x] := W.gval[dir][x]
      for x in sub:
        let
          w = W.gval[dir][x]
          d = ds.gval[dir][x]
        var q, M, F, E {.noinit.}: evalType(w)
        q := w * d.adj
        M := alpha.sval * q
        F.projectTAH M
        E[] := expAH(F[])
        z.expa[x] := E
        z.m[x] := M
        z.gval[dir][x] := E * w
    toc("stoutStepForward kernel end")

  proc backward(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
    raiseUnsupportedPath("stoutStepForward backward", "use the grouped stout outputs")

  let g = W.gval.newOneOf
  g.zeroGaugeStorage
  graphNode(
    GstoutStepForward(
      runtime: W.runtime,
      gval: g,
      expa: W.gval[dir].newOneOf,
      m: W.gval[dir].newOneOf),
    @[Gvalue(args)],
    Gfunc(forward: forward, backward: backward, name: "stoutStepForward"),
    "stoutStepForward")

proc stoutStepLogDetValue(x: GstoutStepForward, parity, dir: int): Gscalar =
  let sub = x.gval.paritySubset(parity)
  proc forward(v: Gvalue) =
    tic("stoutStepLogDetJ kernel")
    let
      x = GstoutStepForward(v.inputs[0])
      z = Gscalar(v)
    var res = 0.0
    threads:
      var s = 0.0
      for i in sub:
        s += simdSum(expProjMulLogJac(x.m[i][]))
      s.threadRankSum
      threadSingle: res = s
    z.sval = res
    toc("stoutStepLogDetJ kernel end")
  proc backward(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
    raiseUnsupportedPath("stoutStepLogDetValue backward", "use the grouped stout outputs")
  graphNode(scalarNodeLike(x), @[Gvalue(x)], Gfunc(forward: forward, backward: backward, name: "stoutStepLogDetValue"), "stoutStepLogDetValue")

proc stoutUpdateLogDetJImpl(W, ds: Ggauge, alpha: Gscalar, c: Gactcoeff, parity, dir: int, fuseHess: static bool): tuple[Wnew: Ggauge, lj: Gscalar] =
  ## Share the update/log-Jacobian forward and pullback.
  W.requireSameGaugeShape(ds, "stoutUpdateLogDetJ")
  let
    sub = W.gval.paritySubset(parity)
    other = W.gval.paritySubset(1 - parity)
    updateBase = stoutStepForwardImpl(W, ds, alpha, parity, dir)
    logdetBase = stoutStepLogDetValue(updateBase, parity, dir)
  when not fuseHess:
    let args = multiValues("stoutStep args", W, ds, alpha)

  proc gradPair(W, ds: Ggauge, alpha: Gscalar, c: Gactcoeff, updateBase: GstoutStepForward, upUpdate: Ggauge, upLog: Gscalar): Gvalue =
    let
      hasUpdate = not upUpdate.isStaticZeroLeaf
      hasLog = not upLog.isStaticZeroLeaf
      terms = if hasUpdate: upUpdate.gaugeAddTerms else: @[]
      nterms = terms.len
    var values = @[Gvalue(W), Gvalue(ds), Gvalue(alpha)]
    when fuseHess:
      values.add Gvalue(c)
      const termIndex = 4
    else:
      const termIndex = 3
    for term in terms:
      values.add Gvalue(term)
    let logIndex = values.len
    if hasLog:
      values.add Gvalue(upLog)
    let updateIndex = values.len
    if hasUpdate or hasLog:
      values.add Gvalue(updateBase)
    let gradArgs = multiValues(
      when fuseHess: "stoutStepGradHess args"
      else: "stoutStepGrad args", values)

    proc kf(v: Gvalue) =
      when fuseHess:
        tic("stoutStepGradHess direction")
      else:
        tic("stoutStepGrad kernel")
      let
        args = Gmulti(v.inputs[0])
        W = Ggauge(args.storedSlot(0))
        ds = Ggauge(args.storedSlot(1))
        alpha = Gscalar(args.storedSlot(2))
        upLog = if hasLog: Gscalar(args.storedSlot(logIndex)) else: nil
        update = if hasUpdate or hasLog:
          GstoutStepForward(args.storedSlot(updateIndex))
        else:
          nil
      var dW: Ggauge
      when fuseHess:
        let pull = GstoutStepPullback(v)
        dW = Ggauge(pull)
      else:
        let
          z = Gmulti(v)
          dds = Ggauge(z.storedSlot(1))
        dW = Ggauge(z.storedSlot(0))
      template upstreamSum(mu, x: untyped): untyped =
        block:
          var r {.noinit.}: evalType(W.gval[mu][x])
          r := Ggauge(args.storedSlot(termIndex)).gval[mu][x]
          for j in 1..<nterms:
            r += Ggauge(args.storedSlot(termIndex + j)).gval[mu][x]
          r
      threads:
        when fuseHess:
          if not hasUpdate:
            for mu in 0..<dW.gval.len:
              dW.gval[mu] := 0.0
            # Full-field and subset loops distribute sites differently.
            threadBarrier()
        if hasUpdate:
          when fuseHess:
            if nterms != 1:
              for mu in 0..<dW.gval.len:
                if mu != dir:
                  for x in dW.gval[mu]:
                    dW.gval[mu][x] := upstreamSum(mu, x)
              for x in other:
                dW.gval[dir][x] := upstreamSum(dir, x)
          else:
            for mu in 0..<dW.gval.len:
              if mu != dir:
                for x in dW.gval[mu]:
                  dW.gval[mu][x] := upstreamSum(mu, x)
            for x in other:
              dW.gval[dir][x] := upstreamSum(dir, x)
        for x in sub:
          let
            w = W.gval[dir][x]
            d = ds.gval[dir][x]
          var rw, rd {.noinit.}: evalType(w)
          rw := 0.0
          rd := 0.0
          if hasUpdate:
            let base = upstreamSum(dir, x)
            var B, P {.noinit.}: evalType(w)
            B := update.expa[x].adj * base
            if hasLog:
              var G {.noinit.}: evalType(w)
              expProjMulLogJacGrad(G[], P[], update.m[x][], (B * w.adj)[])
              let H = alpha.sval * P + (upLog.sval * alpha.sval) * G
              rw := B + H * d
              rd := H.adj * w
            else:
              expProjectTAHPullback(P[], update.m[x][], (B * w.adj)[])
              let H = alpha.sval * P
              rw := B + H * d
              rd := H.adj * w
          elif hasLog:
            var G {.noinit.}: evalType(w)
            expProjMulLogJacGrad(G[], update.m[x][])
            let H = (upLog.sval * alpha.sval) * G
            rw := H * d
            rd := H.adj * w
          dW.gval[dir][x] := rw
          when fuseHess:
            pull.hdir[x] := rd
          else:
            dds.gval[dir][x] := rd
      when fuseHess:
        toc("stoutStepGradHess direction end")
        let coeff = Gactcoeff(args.storedSlot(3))
        if hasUpdate and nterms == 1:
          let base = Ggauge(args.storedSlot(termIndex))
          coeff.cval.gaugeDerivDeriv2SubsetAddBase(W.gval, pull.hdir, base.gval, pull.gval, parity, dir)
        else:
          coeff.cval.gaugeDerivDeriv2SubsetAdd(W.gval, pull.hdir, pull.gval, parity, dir)
      else:
        toc("stoutStepGrad kernel end")

    proc kb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
      raiseUnsupportedPath("stoutStep gradient backward", "higher stout derivatives")

    when fuseHess:
      let z = W.gaugeNodeLike
      result = graphNode(
        GstoutStepPullback(
          runtime: z.runtime,
          gval: z.gval,
          hdir: W.gval[dir].newOneOf),
        @[Gvalue(gradArgs)],
        Gfunc(forward: kf, backward: kb, name: "stoutStepGradHess"),
        "stoutStepGradHess")
    else:
      let z = newMultiOutputNode(@[Gvalue(W), Gvalue(ds)], @[Gvalue(gradArgs)], Gfunc(forward: kf, backward: kb, name: "stoutStepGrad"), "stoutStepGrad")
      Ggauge(z.storedSlot(0)).zeroGaugeStorage
      Ggauge(z.storedSlot(1)).zeroGaugeStorage
      result = Gvalue(z)

  proc alphaGrad(W, ds: Ggauge, alpha: Gscalar, updateBase: GstoutStepForward, upUpdate: Ggauge, upLog: Gscalar): Gscalar =
    let
      hasUpdate = not upUpdate.isStaticZeroLeaf
      hasLog = not upLog.isStaticZeroLeaf
      terms = if hasUpdate: upUpdate.gaugeAddTerms else: @[]
      nterms = terms.len
    if not hasUpdate and not hasLog:
      return Gscalar(alpha.zeroLike)
    var values = @[Gvalue(W), Gvalue(ds), Gvalue(alpha)]
    for term in terms:
      values.add Gvalue(term)
    let logIndex = values.len
    if hasLog:
      values.add Gvalue(upLog)
    let updateIndex = values.len
    if hasLog:
      values.add Gvalue(updateBase)
    let gradArgs = multiValues("stoutStepAlphaGrad args", values)

    proc kf(v: Gvalue) =
      tic("stoutStepAlphaGrad kernel")
      let
        args = Gmulti(v.inputs[0])
        W = Ggauge(args.storedSlot(0))
        ds = Ggauge(args.storedSlot(1))
        alpha = Gscalar(args.storedSlot(2))
        upLog = if hasLog: Gscalar(args.storedSlot(logIndex)) else: nil
        update = if hasLog:
          GstoutStepForward(args.storedSlot(updateIndex))
        else:
          nil
        z = Gscalar(v)
      var da = 0.0
      threads:
        var s = 0.0
        for x in sub:
          let
            w = W.gval[dir][x]
            d = ds.gval[dir][x]
          var q {.noinit.}: evalType(w)
          q := w * d.adj
          if hasUpdate:
            var b {.noinit.}: evalType(w)
            b := Ggauge(args.storedSlot(3)).gval[dir][x]
            for j in 1..<nterms:
              b += Ggauge(args.storedSlot(3 + j)).gval[dir][x]
            var F, D {.noinit.}: evalType(w)
            F.projectTAH q
            D := expDeriv(alpha.sval * F, b * w.adj)
            s += simdSum(redot(D, F))
          if hasLog:
            var G {.noinit.}: evalType(w)
            expProjMulLogJacGrad(G[], update.m[x][])
            s += upLog.sval * simdSum(redot(G, q))
        s.threadRankSum
        threadSingle: da = s
      z.sval = da
      toc("stoutStepAlphaGrad kernel end")

    proc kb(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
      raiseUnsupportedPath("stoutStep alpha gradient backward", "higher stout derivatives")

    graphNode(scalarNodeLike(alpha), @[Gvalue(gradArgs)], Gfunc(forward: kf, backward: kb, name: "stoutStepAlphaGrad"), "stoutStepAlphaGrad")

  proc parentForward(v: Gvalue) = discard

  proc parentBackward(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
    when fuseHess:
      let
        W = Ggauge(z.inputs[0])
        c = Gactcoeff(z.inputs[1])
        alpha = Gscalar(z.inputs[2])
        ds = Ggauge(z.inputs[3])
        updateBase = GstoutStepForward(z.inputs[4])
        upstream = Gmulti(rootedUpstream(zb, z))
        upUpdate = Ggauge(upstream[0])
        upLog = Gscalar(upstream[1])
      case i
      of 0:
        gradPair(W, ds, alpha, c, updateBase, upUpdate, upLog)
      of 1:
        raiseUnsupportedPath("stoutStep backward", "derivative with respect to gauge-action coefficients")
      of 2:
        Gvalue(alphaGrad(W, ds, alpha, updateBase, upUpdate, upLog))
      else:
        raiseValueError("stoutStep parent has three backward inputs")
    else:
      if i != 0:
        raiseValueError("stoutStep parent has one packed input")
      let
        args = Gmulti(z.inputs[0])
        W = Ggauge(args[0])
        ds = Ggauge(args[1])
        alpha = Gscalar(args[2])
        updateBase = GstoutStepForward(z.inputs[1])
        upstream = Gmulti(rootedUpstream(zb, z))
        upUpdate = Ggauge(upstream[0])
        upLog = Gscalar(upstream[1])
        pair = Gmulti(gradPair(W, ds, alpha, nil, updateBase, upUpdate, upLog))
        dW = gaugeGradSlot(pair, 0)
      Gvalue(multiValues("stoutStep input gradients", dW, gaugeGradSlot(pair, 1), alphaGrad(W, ds, alpha, updateBase, upUpdate, upLog)))

  var parentInputs: seq[Gvalue]
  var parentFunc: Gfunc
  when fuseHess:
    parentInputs = @[Gvalue(W), Gvalue(c), Gvalue(alpha), Gvalue(ds), Gvalue(updateBase), Gvalue(logdetBase)]
    parentFunc = Gfunc(
      forward: parentForward,
      backward: parentBackward,
      inputView: stoutActionStepParentInputView,
      name: "stoutStep")
    const
      updateBaseInput = 4
      logdetBaseInput = 5
  else:
    parentInputs = @[Gvalue(args), Gvalue(updateBase), Gvalue(logdetBase)]
    parentFunc = Gfunc(
      forward: parentForward,
      backward: parentBackward,
      inputView: stoutStepParentInputView,
      name: "stoutStep")
    const
      updateBaseInput = 1
      logdetBaseInput = 2
  let parent = newMultiStructureNode(@[Gvalue(updateBase), Gvalue(logdetBase)], parentInputs, parentFunc, "stoutStep")

  proc updateForward(v: Gvalue) =
    Ggauge(v).gval = Ggauge(v.inputs[1]).gval

  proc updateBackward(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
    if i != 0:
      raiseValueError("stoutStep update view has one backward input")
    let upstream = Ggauge(rootedUpstream(zb, z))
    let parent = Gmulti(z.inputs[0])
    Gvalue(multiValues("stoutStep update cotangents", upstream, parent.inputs[logdetBaseInput].zeroLike))

  let updateView = Ggauge(runtime: updateBase.runtime, gval: updateBase.gval)
  result.Wnew = graphNode(
    updateView, @[Gvalue(parent), Gvalue(updateBase)],
    Gfunc(forward: updateForward, backward: updateBackward, inputView: stoutStepOutputView, name: "stoutStepUpdate"),
    "stoutStepUpdate")

  proc logdetForward(v: Gvalue) =
    Gscalar(v).sval = Gscalar(v.inputs[1]).sval

  proc logdetBackward(zb: Gvalue, z: Gvalue, i: int, input: Gvalue): Gvalue =
    if i != 0:
      raiseValueError("stoutStep logdet view has one backward input")
    let upstream = Gscalar(rootedUpstream(zb, z))
    let parent = Gmulti(z.inputs[0])
    Gvalue(multiValues("stoutStep logdet cotangents", parent.inputs[updateBaseInput].zeroLike, upstream))

  result.lj = graphNode(
    scalarNodeLike(logdetBase), @[Gvalue(parent), Gvalue(logdetBase)],
    Gfunc(forward: logdetForward, backward: logdetBackward, inputView: stoutStepOutputView, name: "stoutStepLogDet"),
    "stoutStepLogDet")

proc stoutUpdateLogDetJ*(W, ds: Ggauge, alpha: Gscalar, parity, dir: int): tuple[Wnew: Ggauge, lj: Gscalar] =
  ## Return one subset update and its log-Jacobian.
  stoutUpdateLogDetJImpl(W, ds, alpha, nil, parity, dir, false)

proc stoutUpdateLogDetJ*(W: Ggauge, c: Gactcoeff, alpha: Gscalar, parity, dir: int): tuple[Wnew: Ggauge, lj: Gscalar] =
  ## As above; form ds internally and fuse its Hessian pullback.
  let ds = gaugeActionDeriv(c, W, parity, dir)
  stoutUpdateLogDetJImpl(W, ds, alpha, c, parity, dir, true)
