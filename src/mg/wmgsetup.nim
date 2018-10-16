import base
import field
import mg/mgargs
import eigens/svdLanczos

type SvdOp*[T,L,F] = object
  op: T
  lo: L

template newVector(op: SvdOp): untyped =
  #op.lo.ColorVector()
  op.F.new(op.lo)
  #op.lo.newField(op.F)
template apply*(sop: SvdOp, r,v: typed) =
  sop.op.apply(r, v)
template applyAdj*(sop: SvdOp, r,v: typed) =
  sop.op.applyDag(r, v)
template newRightVec*(op: SvdOp): untyped =
  #newVector(op).even
  newVector(op)
template newLeftVec*(op: SvdOp): untyped =
  #newVector(op).odd
  newVector(op)

proc mgsetupSvd*(r,p: var MgTransfer, op: any, x: Field) =
  var op = SvdOp[type(op),type(x.l),type(x)](op: op, lo: x.l)
  var src = newOneOf(x)
  src := 1
  let nv = p.v[0].len
  var sv = newSeq[float](nv)
  var qv = newSeq[type(op.newRightVec)](nv)
  var qva = newSeq[type(op.newLeftVec)](nv)
  for i in 0..<nv:
    qv[i] = op.newRightVec()
    qva[i] = op.newLeftVec()
  var maxit = intParam("maxit", 1000)
  var emin = 0.0
  var emax = 10.0
  var verb = 0

  let its = svdLanczos(op, src, sv, qv, qva, 0.0, maxit, emin, emax, verb)
  echo "SVD: its"
  for i in 0..<nv:
    echo i, ": ", $sv[i]

  r.v.wmgzero
  p.v.wmgzero
  for i in 0..<nv:
    #op.apply(qva[i], qv[i])
    #p.wmgBlockNormalizeInsert(qv[i], i, x, op.op.cb)
    qv[i].wmgProject(p)
    qv[i].normalize
    p.v.wmgInsert(qv[i], i)
    #r.wmgBlockNormalizeInsert(qva[i], i, x, op.op.cb)
    #qva[i].wmgProject(r)
    #qva[i].normalize
    #r.v.wmgInsert(qva[i], i)
    r.v.wmgInsert(qv[i], i)
