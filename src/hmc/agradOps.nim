import hmc/agrad, base/qexInternal, field, gauge, maths/matrixFunctions
export agrad, qexInternal, field, gauge
import physics/stagD

type
  FloatV* = AgVar[float]
  #GaugeF* = seq[Field]
  #GaugeFV* = AgVar[GaugeF]
  GaugeF*[V:static int,T] = seq[Field[V,T]]
  GaugeFV*[V:static int,T] = AgVar[GaugeF[V,T]]
  FieldV*[V:static int,T] = AgVar[Field[V,T]]

template newFloatV*(c: AgTape, x = 0.0): auto =
  var t = FloatV.new()
  t.doGrad = true
  t.ctx = c
  t.obj = x
  t
template newAgVar*(c: AgTape, x: auto): auto =
  var t = AgVar[typeof(x)].new()
  t.doGrad = true
  t.ctx = c
  t.obj = x
  t.grad = x.newOneOf
  t
template newGaugeFV*(c: AgTape, x: GaugeF): auto =
  var t = AgVar[typeof(x)].new()
  t.doGrad = true
  t.ctx = c
  t.obj = x
  t.grad = x.newOneOf
  t

template maybeObj(x: auto): auto =
  when x is AgVar: x.obj else: x
template `maybeObj=`(x,y: auto) =
  when x is AgVar: x.obj = y else: x = y
template zero(r: float) =  r = 0.0
proc zero[F:Field](r: F) =
  threads:
    r := 0
proc zero[G:GaugeF](r: G) =
  threads:
    for mu in 0..<r.len:
      r[mu] := 0
proc peqmul[F:Field](r: F, x: float, y: F) =
  threads:
    for s in r:
      r[s] += x * y[s]
proc peqmul[G:GaugeF](r: G, x: float, y: G) =
  threads:
    for mu in 0..<r.len:
      for s in r[mu]:
        r[mu][s] += x * y[mu][s]

proc assigngrad(c: float, r: var float) =
  r += c
proc assignfwd[I,O](op: AgOp[I,O]) {.nimcall.} =
  #mixin peq
  op.outputs.obj := op.inputs.maybeObj
  when op.inputs is AgVar:
    zero op.inputs.grad
proc assignbck[I,O](op: AgOp[I,O]) {.nimcall.} =
  #mixin peq
  when op.inputs is AgVar:
    if op.inputs.doGrad:
      assigngrad(op.outputs.grad, op.inputs.grad)
proc assign(c: var AgTape, r: AgVar, x: auto) =
  var op = newAgOp(x, r, assignfwd, assignbck)
  c.add op
template assign*(r: AgVar, x: auto) =
  r.ctx.assign(r, x)
template `:=`*(r: AgVar, x: auto) =
  r.ctx.assign(r, x)

proc addgrad1(c: float, r: var float, y: float) =
  r += c
proc addgrad2(c: float, x: float, r: var float) =
  r += c

proc addfwd[I,O](op: AgOp[I,O]) {.nimcall.} =
  mixin add
  add(op.outputs.obj, op.inputs[0].maybeObj, op.inputs[1].maybeObj)
  when op.inputs[0] is AgVar:
    zero op.inputs[0].grad
  when op.inputs[1] is AgVar:
    zero op.inputs[1].grad
proc addbck[I,O](op: AgOp[I,O]) {.nimcall.} =
  mixin add
  when op.inputs[0] is AgVar:
    if op.inputs[0].doGrad:
      addgrad1(op.outputs.grad, op.inputs[0].grad, op.inputs[1].maybeObj)
  when op.inputs[1] is AgVar:
    if op.inputs[1].doGrad:
      addgrad2(op.outputs.grad, op.inputs[0].maybeObj, op.inputs[1].grad)
proc add(c: var AgTape, r: AgVar, x: auto, y: auto) =
  var op = newAgOp((x,y), r, addfwd, addbck)
  c.add op
template add*(r: AgVar, x: auto, y: auto) =
  r.ctx.add(r, x, y)

proc subgrad1(c: float, r: var float, y: float) =
  r += c
proc subgrad2(c: float, x: float, r: var float) =
  r -= c

proc subfwd[I,O](op: AgOp[I,O]) {.nimcall.} =
  mixin sub
  sub(op.outputs.obj, op.inputs[0].maybeObj, op.inputs[1].maybeObj)
  when op.inputs[0] is AgVar:
    zero op.inputs[0].grad
  when op.inputs[1] is AgVar:
    zero op.inputs[1].grad
proc subbck[I,O](op: AgOp[I,O]) {.nimcall.} =
  mixin sub
  when op.inputs[0] is AgVar:
    if op.inputs[0].doGrad:
      subgrad1(op.outputs.grad, op.inputs[0].grad, op.inputs[1].maybeObj)
  when op.inputs[1] is AgVar:
    if op.inputs[1].doGrad:
      subgrad2(op.outputs.grad, op.inputs[0].maybeObj, op.inputs[1].grad)
proc sub(c: var AgTape, r: AgVar, x: auto, y: auto) =
  var op = newAgOp((x,y), r, subfwd, subbck)
  c.add op
template sub*(r: AgVar, x: auto, y: auto) =
  r.ctx.sub(r, x, y)

proc norm2subfwd[I,O](op: AgOp[I,O]) {.nimcall.} =
  mixin norm2subtract, zero
  op.outputs.obj = norm2subtract(op.inputs[0].obj, op.inputs[1])
  zero op.inputs[0].grad
proc norm2subbck[I,O](op: AgOp[I,O]) {.nimcall.} =
  mixin peqmul
  let x = 2.0*op.outputs.grad
  if op.inputs[0].doGrad:
    peqmul(op.inputs[0].grad, x, op.inputs[0].obj)
proc norm2subtract(c: var AgTape, r: AgVar, x: auto, y: float) =
  var op = newAgOp((x,y), r, norm2subfwd, norm2subbck)
  c.add op
template norm2subtract*(r: AgVar, x: auto, y: float) =
  r.ctx.norm2subtract(r, x, y)

proc mulgrad1(c: float, r: var float, y: float) =
  r += c * y
proc mulgrad2(c: float, x: float, r: var float) =
  r += x * c

proc mul[F:Field](r: F, x: float, y: F) =
  threads:
    r := x * y
proc mulgrad1[F:Field](c: F, r: var float, y: F) =
  var rr = 0.0
  threads:
    var l: typeof(redot(y[0][], c[0][]))
    for s in c:
      l += redot(y[s][], c[s][])
    var m = simdReduce l
    threadRankSum m
    threadSingle:
      rr = m
  r += rr
proc mulgrad2[F:Field](c: F, x: float, r: F) =
  threads:
    for s in c:
      r[s][] += x * c[s][]

proc mul[G:GaugeF](r: G, x: float, y: G) =
  threads:
    for mu in 0..<r.len:
      r[mu] := x * y[mu]
proc mulgrad1[G:GaugeF](c: G, r: var float, y: G) =
  var rr = 0.0
  threads:
    var l: typeof(redot(y[0][0][], c[0][0][]))
    for mu in 0..<c.len:
      for s in c[mu]:
        l += redot(y[mu][s][], c[mu][s][])
    var m = simdReduce l
    threadRankSum m
    threadSingle:
      rr = m
  r += rr
proc mulgrad2[G:GaugeF](c: G, x: float, r: G) =
  threads:
    for mu in 0..<c.len:
      for s in c[mu]:
        r[mu][s][] += x * c[mu][s][]

proc mul[G:GaugeF](r: G, x: G, y: G) =
  threads:
    for mu in 0..<r.len:
      r[mu] := x[mu] * y[mu]
proc mulgrad1[G:GaugeF](c: G, r: G, y: G) =
  threads:
    for mu in 0..<c.len:
      for s in c[mu]:
        r[mu][s][] += c[mu][s][] * y[mu][s][].adj
proc mulgrad2[G:GaugeF](c: G, x: G, r: G) =
  threads:
    for mu in 0..<c.len:
      for s in c[mu]:
        r[mu][s][] += x[mu][s][].adj * c[mu][s][]

proc mulfwd[I,O](op: AgOp[I,O]) {.nimcall.} =
  mixin mul
  mul(op.outputs.obj, op.inputs[0].maybeObj, op.inputs[1].maybeObj)
  when op.inputs[0] is AgVar:
    zero op.inputs[0].grad
  when op.inputs[1] is AgVar:
    zero op.inputs[1].grad
proc mulbck[I,O](op: AgOp[I,O]) {.nimcall.} =
  mixin mul
  when op.inputs[0] is AgVar:
    if op.inputs[0].doGrad:
      mulgrad1(op.outputs.grad, op.inputs[0].grad, op.inputs[1].maybeObj)
  when op.inputs[1] is AgVar:
    if op.inputs[1].doGrad:
      mulgrad2(op.outputs.grad, op.inputs[0].maybeObj, op.inputs[1].grad)
proc mul(c: var AgTape, r: AgVar, x: auto, y: auto) =
  var op = newAgOp((x,y), r, mulfwd, mulbck)
  c.add op
template mul*(r: AgVar, x: auto, y: auto) =
  r.ctx.mul(r, x, y)

proc divdgrad1(c: float, r: var float, y: float) =
  r += c / y
proc divdgrad2(c: float, x,y: float, r: var float) =
  r -= x * c / (y*y)

proc divdfwd[I,O](op: AgOp[I,O]) {.nimcall.} =
  mixin divd
  divd(op.outputs.obj, op.inputs[0].maybeObj, op.inputs[1].maybeObj)
  when op.inputs[0] is AgVar:
    zero op.inputs[0].grad
  when op.inputs[1] is AgVar:
    zero op.inputs[1].grad
proc divdbck[I,O](op: AgOp[I,O]) {.nimcall.} =
  mixin divd
  when op.inputs[0] is AgVar:
    if op.inputs[0].doGrad:
      divdgrad1(op.outputs.grad, op.inputs[0].grad, op.inputs[1].maybeObj)
  when op.inputs[1] is AgVar:
    if op.inputs[1].doGrad:
      divdgrad2(op.outputs.grad, op.inputs[0].maybeObj,
                op.inputs[1].maybeObj, op.inputs[1].grad)
proc divd(c: var AgTape, r: AgVar, x: auto, y: auto) =
  var op = newAgOp((x,y), r, divdfwd, divdbck)
  c.add op
template divd*(r: AgVar, x: auto, y: auto) =
  r.ctx.divd(r, x, y)

proc peqgrad(c: float, r: var float) =
  r += c
proc peqfwd[I,O](op: AgOp[I,O]) {.nimcall.} =
  #mixin peq
  op.outputs.obj += op.inputs.maybeObj
  when op.inputs is AgVar:
    zero op.inputs.grad
proc peqbck[I,O](op: AgOp[I,O]) {.nimcall.} =
  #mixin peq
  when op.inputs is AgVar:
    if op.inputs.doGrad:
      peqgrad(op.outputs.grad, op.inputs.grad)
proc peq(c: var AgTape, r: AgVar, x: auto) =
  var op = newAgOp(x, r, peqfwd, peqbck)
  c.add op
template peq*(r: AgVar, x: auto) =
  r.ctx.peq(r, x)
template `+=`*(r: AgVar, x: auto) =
  r.ctx.peq(r, x)

proc mulna[G:GaugeF](r: G, x: G, y: G) =
  threads:
    for mu in 0..<r.len:
      r[mu] := x[mu] * y[mu].adj
proc mulnagrad1[G:GaugeF](c: G, r: G, y: G) =
  threads:
    for mu in 0..<c.len:
      for s in c[mu]:
        r[mu][s][] += c[mu][s][] * y[mu][s][]
proc mulnagrad2[G:GaugeF](c: G, x: G, r: G) =
  threads:
    for mu in 0..<c.len:
      for s in c[mu]:
        r[mu][s][] += c[mu][s][].adj * x[mu][s][]

proc mulnafwd[I,O](op: AgOp[I,O]) {.nimcall.} =
  mulna(op.outputs.obj, op.inputs[0].maybeObj, op.inputs[1].maybeObj)
  when op.inputs[0] is AgVar:
    zero op.inputs[0].grad
  when op.inputs[1] is AgVar:
    zero op.inputs[1].grad
proc mulnabck[I,O](op: AgOp[I,O]) {.nimcall.} =
  when op.inputs[0] is AgVar:
    if op.inputs[0].doGrad:
      mulnagrad1(op.outputs.grad, op.inputs[0].grad, op.inputs[1].maybeObj)
  when op.inputs[1] is AgVar:
    if op.inputs[1].doGrad:
      mulnagrad2(op.outputs.grad, op.inputs[0].maybeObj, op.inputs[1].grad)
proc mulna(c: var AgTape, r: AgVar, x: auto, y: auto) =
  var op = newAgOp((x,y), r, mulnafwd, mulnabck)
  c.add op
template mulna*(r: AgVar, x: auto, y: auto) =
  r.ctx.mulna(r, x, y)

var ep: ExpParam
ep.scale = 20
ep.kind = ekPoly
ep.order = 4
proc expvfwd[I,O](op: AgOp[I,O]) {.nimcall.} =
  mixin exp
  let g = op.outputs.obj
  let a = op.inputs[0].obj
  op.inputs[0].grad = 0
  let p = op.inputs[1].obj
  let pg = op.inputs[1].grad
  threads:
    for mu in 0..<g.len:
      for s in g[mu]:
        g[mu][s][] := ep.exp(a*p[mu][s][])
        pg[mu][s] := 0
proc expvbck[I,O](op: AgOp[I,O]) {.nimcall.} =
  mixin `+=`
  if op.inputs[0].doGrad or op.inputs[1].doGrad:
    #let g = op.outputs.obj
    let gg = op.outputs.grad
    let a = op.inputs[0].obj
    var ag = op.inputs[0].grad
    let p = op.inputs[1].obj
    let pg = op.inputs[1].grad
    threads:
      var agl:typeof(redot(p[0][0][],p[0][0][]))
      for mu in 0..<gg.len:
        for s in gg[mu]:
          let t = a*p[mu][s][]
          let d = ep.expDeriv(t, gg[mu][s][])
          if op.inputs[0].doGrad:
            agl += redot(p[mu][s][], d)
          if op.inputs[1].doGrad:
            pg[mu][s][] += a * d
      var agls = agl.simdReduce
      threadRankSum agls
      threadSingle:
        ag += agls
    op.inputs[0].grad = ag
proc exp*[G:GaugeFV](c: var AgTape, g: G, a: FloatV, p: G) =
  var op = newAgOp((a,p), g, expvfwd, expvbck)
  c.add op
template exp*[G:GaugeFV](g: G, a: FloatV, p: G) =
  g.ctx.exp(g, a, p)

proc expfwd[I,O](op: AgOp[I,O]) {.nimcall.} =
  mixin exp
  let g = op.outputs.obj
  let p = op.inputs.obj
  let pg = op.inputs.grad
  threads:
    for mu in 0..<g.len:
      for s in g[mu]:
        g[mu][s][] := ep.exp(p[mu][s][])
        pg[mu][s] := 0
proc expbck[I,O](op: AgOp[I,O]) {.nimcall.} =
  mixin `+=`
  if op.inputs.doGrad:
    #let g = op.outputs.obj
    let gg = op.outputs.grad
    let p = op.inputs.obj
    let pg = op.inputs.grad
    threads:
      for mu in 0..<gg.len:
        for s in gg[mu]:
          let d = ep.expDeriv(p[mu][s][], gg[mu][s][])
          pg[mu][s][] += d
proc exp*[G:GaugeFV](c: var AgTape, g: G, p: G) =
  var op = newAgOp(p, g, expfwd, expbck)
  c.add op
template exp*[G:GaugeFV](g: G, p: G) =
  g.ctx.exp(g, p)

proc gactionfwd[I,O](op: AgOp[I,O]) {.nimcall.} =
  op.outputs.obj = actionA(op.inputs[0], op.inputs[1].obj)
  zero op.inputs[1].grad
proc gactionbck[I,O](op: AgOp[I,O]) {.nimcall.} =
  if op.inputs[1].doGrad:
    var gc = op.outputs.grad * op.inputs[0]
    gaugeDeriv2(gc, op.inputs[1].obj, op.inputs[1].grad)
proc gaction(c: var AgTape, gc: GaugeActionCoeffs, r: AgVar, x: auto) =
  var op = newAgOp((gc,x), r, gactionfwd, gactionbck)
  c.add op
template gaction*(gc: GaugeActionCoeffs, r: AgVar, x: auto) =
  r.ctx.gaction(gc, r, x)

proc gderivfwd[I,O](op: AgOp[I,O]) {.nimcall.} =
  zero op.outputs.obj
  gaugeDeriv2(op.inputs[0], op.inputs[1].obj, op.outputs.obj)
  zero op.inputs[1].grad
proc gderivbck[I,O](op: AgOp[I,O]) {.nimcall.} =
  if op.inputs[1].doGrad:
    gaugeDerivDeriv2(op.inputs[0], op.inputs[1].obj, op.outputs.grad, op.inputs[1].grad)
proc gderiv(c: var AgTape, gc: GaugeActionCoeffs, r: AgVar, x: auto) =
  var op = newAgOp((gc,x), r, gderivfwd, gderivbck)
  c.add op
template gderiv*(gc: GaugeActionCoeffs, r: AgVar, x: auto) =
  r.ctx.gderiv(gc, r, x)

proc projtahfwd[I,O](op: AgOp[I,O]) {.nimcall.} =
  let f = op.outputs.obj
  let g = op.inputs.obj
  threads:
    for mu in 0..<f.len:
      f[mu].projectTAH g[mu]
  zero op.inputs.grad
proc projtahbck[I,O](op: AgOp[I,O]) {.nimcall.} =
  if op.inputs.doGrad:
    let f = op.inputs.grad
    let g = op.outputs.grad
    threads:
      for mu in 0..<f.len:
        for s in f[mu]:
          var t = g[mu][s][]
          t.projectTAH
          f[mu][s][] += t
proc projtah(c: var AgTape, r: AgVar, x: auto) =
  var op = newAgOp(x, r, projtahfwd, projtahbck)
  c.add op
template projtah*(r: AgVar, x: auto) =
  r.ctx.projtah(r, x)

proc maskEvenFwd[I,O](op: AgOp[I,O]) {.nimcall.} =
  let r = op.outputs.obj
  let g = op.outputs.grad
  threads:
    r.odd := 0
    g := 0
proc maskEvenBck[I,O](op: AgOp[I,O]) {.nimcall.} =
  #let g = op.inputs.grad
  let g = op.outputs.grad
  threads:
    g.odd := 0
proc maskEven[F:FieldV](c: var AgTape, r: F) =
  var op = newAgOp(0, r, maskEvenFwd, maskEvenBck)
  c.add op
template maskEven*[F:FieldV](r: F) =
  r.ctx.maskEven(r)

proc xpayfwd[I,O](op: AgOp[I,O]) {.nimcall.} =
  let r = op.outputs.obj
  let x = op.inputs[0].obj
  let a = op.inputs[1].obj
  op.inputs[1].grad = 0
  let y = op.inputs[2].obj
  let xg = op.inputs[0].grad
  let yg = op.inputs[2].grad
  threads:
    for mu in 0..<r.len:
      for s in r[mu]:
        r[mu][s][] := x[mu][s][] + a * y[mu][s][]
        xg[mu][s][] := 0
        yg[mu][s][] := 0
proc xpaybck[I,O](op: AgOp[I,O]) {.nimcall.} =
  let rg = op.outputs.grad
  let a = op.inputs[1].obj
  let y = op.inputs[2].obj
  let xg = op.inputs[0].grad
  let yg = op.inputs[2].grad
  var ag = 0.0
  threads:
    var al: typeof(redot(y[0][0][], rg[0][0][]))
    for mu in 0..<rg.len:
      for s in rg[mu]:
        if op.inputs[0].doGrad:
          xg[mu][s][] += rg[mu][s][]
        if op.inputs[2].doGrad:
          yg[mu][s][] += a * rg[mu][s][]
        al += redot(y[mu][s][], rg[mu][s][])
    var am = simdReduce al
    threadRankSum am
    threadSingle:
      ag = am
  if op.inputs[1].doGrad:
    op.inputs[1].grad += ag
proc xpay[G:GaugeFV](c: var AgTape, r: G, x: G, a: FloatV, y: G) =
  var op = newAgOp((x,a,y), r, xpayfwd, xpaybck)
  c.add op
template xpay*[G:GaugeFV](r: G, x: G, a: FloatV, y: G) =
  r.ctx.xpay(r, x, a, y)

proc agradDfwd[I,O](op: AgOp[I,O]) {.nimcall.} =
  var s = op.inputs[0]
  let g = op.inputs[1]
  let x = op.inputs[2]
  let m = op.inputs[3]
  let r = op.outputs
  let g0 = s.g
  s.g = g.maybeObj
  s.rephase
  s.D(r.obj, x.maybeObj, m.maybeObj)
  s.rephase
  s.g = g0
  when g is AgVar:
    zero g.grad
  when x is AgVar:
    zero x.grad
  when m is AgVar:
    zero m.grad
proc agradDbck[I,O](op: AgOp[I,O]) {.nimcall.} =
  var s = op.inputs[0]
  let g = op.inputs[1]
  let x = op.inputs[2]
  let m = op.inputs[3]
  let r = op.outputs
  when g is AgVar:
    if g.doGrad:
      #g.grad += rephase [outer(c shift x') - outer(x shift c')]
      for mu in 0..<g.grad.len: g.grad[mu] *= 2.0
      s.rephase g.grad
      s.stagD2deriv(g.grad, r.grad, x.maybeObj)
      s.rephase g.grad
      for mu in 0..<g.grad.len: g.grad[mu] *= 0.5
  when x is AgVar:
    if x.doGrad:
      let g0 = s.g
      s.g = g.maybeObj
      s.rephase
      s.peqDdag(x.grad, r.grad, m.maybeObj)
      s.rephase
      s.g = g0
  when m is AgVar:
    if m.doGrad:
      m.grad += redot(x.maybeObj, r.grad)
proc agradD(c: var AgTape, s,g,r,x,m: auto) =
  var op = newAgOp((s,g,x,m), r, agradDfwd, agradDbck)
  c.add op
template agradD*(s: Staggered, g,r,x,m: auto) =
  r.ctx.agradD(s, g, r, x, m)

proc agradStagDerivfwd[I,O](op: AgOp[I,O]) {.nimcall.} =
  var s = op.inputs[0]
  let x = op.inputs[1]
  let f = op.outputs
  zero f.obj
  stagDeriv(s, f.obj, x.obj)
  zero x.grad
proc agradStagDerivbck[I,O](op: AgOp[I,O]) {.nimcall.} =
  var s = op.inputs[0]
  let x = op.inputs[1]
  let f = op.outputs
  if x.doGrad:
    let g0 = s.g
    x.obj.even *= 2
    x.obj.odd *= -2
    s.g = f.grad
    s.rephase
    s.peqDdag(x.grad, x.obj, 0.0)
    s.rephase
    s.g = g0
    x.obj.even *= 0.5
    x.obj.odd *= -0.5
proc agradStagDeriv(c: var AgTape, s,f,x: auto) =
  var op = newAgOp((s,x), f, agradStagDerivfwd, agradStagDerivbck)
  c.add op
template agradStagDeriv*(s: Staggered, f,x: auto) =
  f.ctx.agradStagDeriv(s, f, x)

proc agradSolvefwd[I,O](op: AgOp[I,O]) {.nimcall.} =
  var s = op.inputs[0]
  let g = op.inputs[1]
  let x = op.inputs[2]
  let m = op.inputs[3]
  let p = op.inputs[4]
  let r = op.outputs
  let g0 = s.g
  s.g = g.maybeObj
  s.rephase
  s.solve(r.obj, x.maybeObj, m.maybeObj, p[])
  s.rephase
  s.g = g0
  when g is AgVar:
    zero g.grad
  when x is AgVar:
    zero x.grad
  when m is AgVar:
    zero m.grad
proc agradSolvebck[I,O](op: AgOp[I,O]) {.nimcall.} =
  var s = op.inputs[0]
  let g = op.inputs[1]
  let x = op.inputs[2]
  let m = op.inputs[3]
  let p = op.inputs[5]
  let r = op.outputs
  var c = r.grad.newOneOf
  r.grad.odd *= -1
  let g0 = s.g
  s.g = g.maybeObj
  s.rephase
  s.solve(c, r.grad, m.maybeObj, p[])
  s.rephase
  s.g = g0
  r.grad.odd *= -1
  c.even *= -1
  when g is AgVar:
    if g.doGrad:
      for mu in 0..<g.grad.len: g.grad[mu] *= 2.0
      s.rephase g.grad
      s.stagD2deriv(g.grad, c, r.obj)
      s.rephase g.grad
      for mu in 0..<g.grad.len: g.grad[mu] *= 0.5
  when x is AgVar:
    if x.doGrad:
      x.grad -= c
  when m is AgVar:
    if m.doGrad:
      m.grad += redot(r.obj, c)
proc agradSolve(c: var AgTape, s,g,r,x,m,pf,pb: auto) =
  var op = newAgOp((s,g,x,m,pf,pb), r, agradSolvefwd, agradSolvebck)
  c.add op
template agradSolve*(s: Staggered, g,r,x,m,pf,pb: auto) =
  ## g: gauge, r: result, x: src, m: mass, p: solve params
  r.ctx.agradSolve(s, g, r, x, m, addr pf, addr pb)
template agradSolve*(s: Staggered, g,r,x,m,p: auto) =
  ## g: gauge, r: result, x: src, m: mass, p: solve params
  r.ctx.agradSolve(s, g, r, x, m, addr p, addr p)

when isMainModule:
  import qex, physics/stagSolve
  qexInit()
  defaultSetup()
  var rs = newRNGField(RngMilc6, lo, intParam("seed", 987654321).uint64)
  g.random rs
  echo "plaq: ", g.plaq

  proc testAgradD =
    let eps = 1e-5
    var t = newAgTape()
    var gv = t.newAgVar(g)
    var s = newStag(g)
    var m = 0.1
    var mv = t.newFloatV(m)
    var nv = t.newFloatV(0.0)
    var v1 = lo.ColorVector()
    var v2 = lo.ColorVector()
    var c = lo.ColorVector()
    var gc = lo.newGauge()
    var v1v = t.newAgVar(v1)
    var v2v = t.newAgVar(v2)
    v1.gaussian rs
    c.gaussian rs
    gc.gaussian rs
    s.agradD(gv, v2v, v1v, mv)
    norm2subtract(nv, v2v, 0.0)
    echo nv.obj
    t.run
    echo nv.obj
    let n0 = nv.obj
    nv.grad = 1.0
    t.grad
    let mg = mv.grad
    let vg = redot(c,v1v.grad)
    var gg = 0.0
    for mu in 0..<gc.len: gg += redot(gc[mu],gv.grad[mu])
    mv.obj += eps
    t.run
    echo "mg: ", mg
    echo "    ", (nv.obj-n0)/eps
    mv.obj = m
    v1 += eps * c
    t.run
    echo "vg: ", vg
    echo "    ", (nv.obj-n0)/eps
    v1 -= eps * c
    for mu in 0..<gc.len: g[mu] += eps * gc[mu]
    t.run
    echo "gg: ", gg
    echo "    ", (nv.obj-n0)/eps
    for mu in 0..<gc.len: g[mu] -= eps * gc[mu]

  proc testAgradSolve =
    let eps = 1e-5
    var t = newAgTape()
    var gv = t.newAgVar(g)
    var s = newStag(g)
    var m = 0.1
    var mv = t.newFloatV(m)
    var nv = t.newFloatV(0.0)
    var v1 = lo.ColorVector()
    var v2 = lo.ColorVector()
    var c = lo.ColorVector()
    var gc = lo.newGauge()
    var v1v = t.newAgVar(v1)
    var v2v = t.newAgVar(v2)
    var sp = initSolverParams()
    sp.r2req = 1e-8
    sp.maxits = 10000
    v1.gaussian rs
    c.gaussian rs
    gc.gaussian rs
    s.agradSolve(gv, v2v, v1v, mv, sp)
    norm2subtract(nv, v2v, 0.0)
    echo nv.obj
    t.run
    echo nv.obj
    let n0 = nv.obj
    nv.grad = 1.0
    t.grad
    let mg = mv.grad
    let vg = redot(c,v1v.grad)
    var gg = 0.0
    for mu in 0..<gc.len: gg += redot(gc[mu],gv.grad[mu])
    mv.obj += eps
    t.run
    echo "mg: ", mg
    echo "    ", (nv.obj-n0)/eps
    mv.obj = m
    v1 += eps * c
    t.run
    echo "vg: ", vg
    echo "    ", (nv.obj-n0)/eps
    v1 -= eps * c
    for mu in 0..<gc.len: g[mu] += eps * gc[mu]
    t.run
    echo "gg: ", gg
    echo "    ", (nv.obj-n0)/eps
    for mu in 0..<gc.len: g[mu] -= eps * gc[mu]

  proc norm2subtract(x: seq, y: float): float =
    for i in 0..<x.len:
      result += x[i].norm2subtract(y)
  proc testAgradStagDeriv =
    let eps = 1e-8
    var t = newAgTape()
    var gv = t.newAgVar(g)
    var s = newStag(g)
    var nv = t.newFloatV(0.0)
    var v1 = lo.ColorVector()
    var v2 = lo.ColorVector()
    var c = lo.ColorVector()
    var gc = lo.newGauge()
    var v1v = t.newAgVar(v1)
    var v2v = t.newAgVar(v2)
    v1.gaussian rs
    c.gaussian rs
    #gc.gaussian rs
    s.agradStagDeriv(gv, v1v)
    norm2subtract(nv, gv, 0.0)
    echo nv.obj
    t.run
    echo nv.obj
    let n0 = nv.obj
    nv.grad = 1.0
    t.grad
    let vg = redot(c,v1v.grad)
    v1 += eps * c
    t.run
    echo "vg: ", vg
    echo "    ", (nv.obj-n0)/eps
    v1 -= eps * c

  testAgradD()
  #testAgradSolve()
  #testAgradStagDeriv()

  qexFinalize()
